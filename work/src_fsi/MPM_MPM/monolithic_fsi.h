#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include "module/bc.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

namespace interface_geometry {
class InterfaceSDF;
}

namespace mpm_mpm_monolithic_fsi {

class MPMMPMMonolithicFSI : public MaterialPoint {
  public:
    static constexpr int max_NR_it = 100;

    stabilizedmpm::StabilizedMPM fluid_;
    implicitmpm::ImplicitSolidMPM solid_;
    CrsMat fsi_sys;

    std::vector<double> nlm_lump_local, nlm_lump;
    std::vector<double> nlambda;      // nodal multiplier force density
    std::vector<char> fixed_dof;      // constrained and inactive components
    std::vector<double> column_scale; // linear increment scaling factors

    /** @brief Configure ten coupled components and the 37 stored fluid-solid-multiplier blocks. */
    MPMMPMMonolithicFSI() {
        this->fluid_.blocks = {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18};
        this->solid_.blocks = {19, 20, 21, 23, 24, 25, 27, 28, 29};
        this->fluid_.rhs_start = 0;
        this->solid_.rhs_start = 4;
        this->integrator_.ode_order = 2;
        this->fsi_sys.ndof = 10;
        this->fsi_sys.block_row = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3,
                                   4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9};
        this->fsi_sys.block_col = {0, 1, 2, 3, 7, 0, 1, 2, 3, 8, 0, 1, 2, 3, 9, 0, 1, 2, 3,
                                   4, 5, 6, 7, 4, 5, 6, 8, 4, 5, 6, 9, 0, 4, 1, 5, 2, 6};
        this->fsi_sys.FEM_flag = false;
        this->fsi_sys.use_petsc = true;
        this->fsi_sys.use_schur_fieldsplit = true;
        this->fsi_sys.amg_rebuild_freq = 1; // rebuild AMG every step
        this->fsi_sys.owner_ = this;
    }

    /**
     * @brief Read and initialize the coupled FSI case input data.
     *
     * Loads shared global parameters, then dispatches fluid and solid
     * boundary-condition/particle input to the corresponding sub-solvers.
     */
    void DataInput();

    /** @brief Integrate the wet solid-surface distance delta into lumped nodal interface-area weights. */
    void LumpedLagrangeMultiplier();

    /**
     * @brief Add interface multiplier forces to the fluid or solid momentum RHS.
     * @param af_coeff Signed generalized-alpha force weight: positive for fluid and negative for solid.
     * @param offsets Start indices of the field's three momentum components in the coupled RHS.
     */
    void AddLagrangeMultiplierToRHS(double af_coeff, const std::vector<int> &offsets);

    /**
     * @brief Assemble the fluid field and its multiplier load into the coupled system.
     * @param nvel_k Fluid nodal velocity at the current endpoint iterate.
     * @param naccel_k Fluid nodal acceleration at the current endpoint iterate.
     */
    void AssembleFluidSystem(const std::vector<double> &nvel_k, //
                             const std::vector<double> &naccel_k);

    /**
     * @brief Assemble the solid field and its multiplier load into the coupled system.
     * @param nvel_k Solid nodal velocity at the current endpoint iterate.
     * @param naccel_k Solid nodal acceleration at the current endpoint iterate.
     * @param stress_k Particle stress state for the current Newton iterate.
     */
    void AssembleSolidSystem(const std::vector<double> &nvel_k,   //
                             const std::vector<double> &naccel_k, //
                             std::vector<std::array<double, 6>> &stress_k);

    /**
     * @brief Constrain the fluid and solid nodal displacement increments to agree on the wet interface.
     * @note Uses the current time-step ndispl fields, not accumulated particle coordinates.
     */
    void AssembleInterfaceSystem();

    /** @brief Solve the coupled fields with interface displacement-increment continuity. */
    void SolveFSISystem();

    /**
     * @brief Move fluid centres on the inside of the final solid phi=0.5 surface directly to that surface.
     * @note Call after both Node2Particle updates and solid MoveParticle, before fluid MoveParticle. Reuses the solid contour
     * and finite-triangle distance routines with an updated temporary phase field. Projection is normal
     * to a triangle interior, or to its nearest edge/vertex. Queries use the existing one-grid-spacing
     * surface band and nearest-triangle sign. Velocity, acceleration, TPIC state, pressure, mass and volume are preserved.
     * This is geometric remapping, not a physical contact-force or momentum correction.
     */
    void CorrectFluidPenetration();

    /**
     * @brief Apply residual constraints and solve for physical monolithic increments.
     * @param NR_it Newton iteration controlling preconditioner setup; each correction starts from zero.
     * @return Number of Krylov iterations used by the coupled linear solve.
     * @note KSP divergence is currently logged by CrsMat but does not stop this
     * routine from rescaling and returning the computed increment.
     */
    int SolveSystem(int NR_it);

    /** @brief Add physical Newton increments to fluid and solid displacements, fluid pressure and interface multipliers. */
    void UpdateNRIncrement() override;

    /** @brief Select active components using field-specific nodal-mass thresholds and physical boundary constraints. */
    void BuildActiveDOFs();

    /** @brief Scale the matrix, RHS and physical initial increment; fixed and inactive guesses are zero. */
    void ScaleSystem();

    /**
     * @brief Register fixed and inactive monolithic components.
     * @param mat Coupled system receiving zero-increment boundary rows.
     */
    void BuildPetscBCList(CrsMat &mat) override;

    /**
     * @brief Zero fixed and inactive residual components.
     * @param rr Coupled residual in component-major order.
     */
    void BCResidualSet(std::vector<double> &rr) override;

    /**
     * @brief Insert the 37 stored scalar blocks without expanding local storage.
     * @param mat Coupled matrix receiving rank-local contributions.
     * @param ndof Number of coupled components per control point.
     */
    void AssemblePetscMat(CrsMat &mat, int ndof) override;

    /**
     * @brief Configure nested fields, multipliers and the default pressure ASM/ILU preconditioner.
     * @param mat Coupled linear system being solved.
     * @param pc PETSc preconditioner for the coupled Krylov solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

    /** @brief Release the coupled PETSc solver and pressure preconditioner resources. */
    ~MPMMPMMonolithicFSI();

  private:
    /**
     * @brief Measure the nonlinear field RHS and direct interface displacement difference on globally owned active components.
     * @return RMS fluid momentum [N], continuity [m^3/s], solid momentum [N] and unweighted interface displacement [m] residuals.
     * @note Requires initialized PETSc row ownership and an unscaled, overlap-synchronized RHS.
     */
    std::array<double, 4> ComputeNRResidualNorms();

    /**
     * @brief Check the updated nonlinear RHS and direct solid-fluid displacement-increment continuity.
     * @param initial_norm Initial RMS residuals saved on iteration zero and used as relative references for the first three
     * fields.
     * @param NR_it Current Newton iteration.
     * @param solver_it Krylov iterations in the last coupled linear solve.
     * @return Whether all field residuals meet tolerance.
     * @note Call after reassembly, before row scaling. The interface test uses solid_.ndispl-fluid_.ndispl without area weights.
     */
    bool CheckNRConvergence(std::array<double, 4> &initial_norm, int NR_it, int solver_it);
    /**
     * @brief Exchange current interface geometry across MPI ranks.
     * @param local_geometry Rank-local solid triangles or nearby fluid particle-domain bounds.
     * @return Physical geometry available to all ranks in the current configuration.
     */
    template <size_t vertex_count>
    std::vector<std::array<std::array<double, 3>, vertex_count>>
    GatherInterfaceGeometry(const std::vector<std::array<std::array<double, 3>, vertex_count>> &local_geometry) const;

    /**
     * @brief Locate the zero contour on a tetrahedron edge by linear interpolation.
     * @param vertices Physical tetrahedron vertices in the current configuration.
     * @param level Signed phase values at the vertices: negative inside the solid.
     * @param i Index of the edge's first vertex.
     * @param j Index of the edge's second vertex.
     * @return Physical position where the interpolated phase value is zero.
     */
    std::array<double, 3> ContourIntersection(const std::array<std::array<double, 3>, 4> &vertices,
                                              const std::array<double, 4> &level, int i, int j) const;

    /**
     * @brief Orient a nonzero-area contour triangle toward the exterior and append it.
     * @param triangle Physical vertices of the reconstructed solid contour.
     * @param outward Direction from the tetrahedron's inside vertices to its outside vertices.
     * @param surface Boundary triangles receiving the oriented contour.
     */
    void AppendInterfaceTriangle(std::array<std::array<double, 3>, 3> triangle, const std::array<double, 3> &outward,
                                 std::vector<std::array<std::array<double, 3>, 3>> &surface) const;

    /**
     * @brief Extract the solid zero contour inside one tetrahedron.
     * @param vertices Physical tetrahedron vertices sampled from the current background cell.
     * @param level Signed phase values: negative inside the solid, positive outside.
     * @param surface Boundary triangles receiving the outward-oriented contour.
     */
    void AppendSolidContour(const std::array<std::array<double, 3>, 4> &vertices, const std::array<double, 4> &level,
                            std::vector<std::array<std::array<double, 3>, 3>> &surface) const;

    /**
     * @brief Interpolate solid phi to mesh vertices and reconstruct its phi=0.5 boundary without subcell sampling.
     * @param solid_phi Control-point solid volume fraction in the configuration to reconstruct.
     * @return Globally available, outward-oriented triangles in physical coordinates.
     * @note Uses xyn and nc for geometry only; field variables remain on control points.
     */
    std::vector<std::array<std::array<double, 3>, 3>> BuildSolidInterface(const std::vector<double> &solid_phi) const;

    /**
     * @brief Reconstruct and exchange finite fluid domains near the solid interface.
     * @return Globally available particle-box bounds in the current physical configuration.
     */
    std::vector<std::array<std::array<double, 3>, 2>> BuildFluidDomains() const;

    /**
     * @brief Index finite fluid boxes within the existing contact activation band.
     * @param domains Current volume-equivalent fluid particle boxes in physical coordinates.
     * @return Spatial bins covering each box expanded by 1.5 background-grid spacings.
     */
    std::map<std::array<int, 3>, std::vector<std::size_t>>
    IndexFluidDomains(const std::vector<std::array<std::array<double, 3>, 2>> &domains) const;

    /**
     * @brief Evaluate wet-interface area density chi_w * delta_epsilon(d_s).
     * @param point Current volume quadrature position.
     * @param interface Indexed solid surface in the current configuration.
     * @param domains Current volume-equivalent fluid particle boxes.
     * @param bins Candidate fluid boxes in the grid-scaled contact band.
     * @return Interface area per volume [1/m], zero outside the band or on dry surfaces.
     */
    double WetInterfaceAreaDensity(const std::array<double, 3> &point, const interface_geometry::InterfaceSDF &interface,
                                   const std::vector<std::array<std::array<double, 3>, 2>> &domains,
                                   const std::map<std::array<int, 3>, std::vector<std::size_t>> &bins) const;
};

} // namespace mpm_mpm_monolithic_fsi
