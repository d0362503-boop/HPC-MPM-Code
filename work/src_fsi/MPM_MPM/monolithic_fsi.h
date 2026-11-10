#pragma once

#include <array>
#include <vector>

#include "module/bc.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

namespace mpm_mpm_monolithic_fsi {

class MPMMPMMonolithicFSI : public MaterialPoint {
  public:
    static constexpr int kStoredBlocks = 37;
    static constexpr int max_NR_it = 100;

    stabilizedmpm::StabilizedMPM fluid_;
    implicitmpm::ImplicitSolidMPM solid_;
    CrsMat fsi_sys;

    std::vector<double> nlm_lump_local, nlm_lump;
    std::vector<double> nlambda;      // nodal multiplier force density
    std::vector<char> fixed_dof;      // constrained and inactive components
    std::vector<double> column_scale; // linear increment scaling factors

    // Driver operations and physics callbacks.
    /** @brief Configure ten coupled components and the 37 stored fluid-solid-multiplier blocks. */
    MPMMPMMonolithicFSI();

    /** @brief Release the coupled PETSc solver and compact equation resources. */
    ~MPMMPMMonolithicFSI();

    /**
     * @brief Read and initialize the coupled FSI case input data.
     *
     * Loads shared global parameters, then dispatches fluid and solid
     * boundary-condition/particle input to the corresponding sub-solvers.
     */
    void DataInput();

    /** @brief Solve the coupled fields with interface endpoint-velocity continuity. */
    void SolveFSISystem();

    /**
     * @brief Apply one final-phase normal correction to penetrated fluid centres.
     * @note Call after both Node2Particle updates and solid MoveParticle, before fluid MoveParticle.
     * Corrects phi>0.5 centres with a usable gradient once toward phi=0.5, within the global mesh.
     * Flat-gradient centres are skipped. Curved interfaces retain the single-update approximation error.
     */
    void CorrectFluidPenetration();

    /** @brief Add physical Newton increments to fluid and solid displacements, fluid pressure and interface multipliers. */
    void UpdateNRIncrement() override;

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
     */
    void AssemblePetscMat(CrsMat &mat) override;

    /**
     * @brief Configure fluid and coupled solid-multiplier field splits with solid ASM/ILU(1).
     * @param mat Coupled linear system being solved.
     * @param pc PETSc preconditioner for the coupled Krylov solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

  private:
    // Coupled assembly.
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
     * @brief Constrain the fluid and solid nodal velocities to agree at the current time-step endpoint.
     * @param nvel_f Fluid nodal velocity at the current endpoint iterate.
     * @param nvel_s Solid nodal velocity at the current endpoint iterate.
     */
    void AssembleInterfaceSystem(const std::vector<double> &nvel_f, const std::vector<double> &nvel_s);

    /**
     * @brief Add interface multiplier forces to the fluid or solid momentum RHS.
     * @param af_coeff Signed generalized-alpha force weight: positive for fluid and negative for solid.
     * @param offsets Start indices of the field's three momentum components in the coupled RHS.
     */
    void AddLagrangeMultiplierToRHS(double af_coeff, const std::vector<int> &offsets);

    // Linear correction and preconditioning.
    /** @brief Select active components using field-specific nodal-mass thresholds and physical boundary constraints. */
    void BuildActiveDOFs();

    /** @brief Scale the coupled matrix and RHS; constrained corrections stay zero. */
    void ScaleSystem();

    /**
     * @brief Apply residual constraints and solve for physical monolithic increments.
     * @param NR_it Newton iteration controlling preconditioner setup; each correction starts from zero.
     * @return Number of Krylov iterations used by the coupled linear solve.
     * @note KSP divergence is currently logged by CrsMat but does not stop this
     * routine from rescaling and returning the computed increment.
     */
    int SolveSystem(int NR_it);

    // Nonlinear residuals.
    /**
     * @brief Measure the nonlinear field RHS and interface endpoint-velocity difference on globally owned active components.
     * @return RMS fluid momentum [N], continuity [m^3/s], solid momentum [N] and unweighted interface velocity [m/s] residuals.
     * @note Requires initialized PETSc row ownership and an unscaled, overlap-synchronized RHS.
     */
    std::array<double, 4> ComputeNRResidualNorms();

    /**
     * @brief Check the updated nonlinear RHS and solid-fluid endpoint-velocity continuity.
     * @param initial_norm Initial RMS residuals saved on iteration zero and used as relative references for the first three
     * fields.
     * @param NR_it Current Newton iteration.
     * @param solver_it Krylov iterations in the last coupled linear solve.
     * @return Whether all field residuals meet tolerance.
     * @note Call after reassembly, before row scaling. The interface test divides its RHS by the active nodal area weight.
     */
    bool CheckNRConvergence(std::array<double, 4> &initial_norm, int NR_it, int solver_it);

    // Interface geometry.
    /** @brief Integrate the two-phase volume-fraction gradient measure into lumped nodal interface-area weights. */
    void LumpedLagrangeMultiplier();
};

} // namespace mpm_mpm_monolithic_fsi
