#pragma once

#include <array>
#include <vector>

#include "module/bc.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

namespace mpmfemblockfsi {

class MPMFEMBlockFSI;

class FSIFluid : public stabilizedfem::StabilizedFEM {
  public:
    MPMFEMBlockFSI &fsi_; // Reference back to the FSI coordinator

    FSIFluid(MPMFEMBlockFSI &fsi) : fsi_(fsi) {}

    // ---- Boundary Condition Overrides ----
    // --- For PETSc ---
    void BuildPetscBCList(CrsMat &mat) override;

    // --- For self-defined ---
    void BCResidualSet(std::vector<double> &rr) override;

  protected:
    /**
     * @brief Apply fluid BCs and the solid interface velocity to the n+1 Newton iterate.
     * @param ts Current time state; its load factor `facl` scales the prescribed values.
     */
    void BCNRSet(const TimeState &ts) override;
};

class FSISolid : public implicitmpm::ImplicitSolidMPM {
  public:
    MPMFEMBlockFSI &fsi_; // Reference back to the FSI coordinator

    FSISolid(MPMFEMBlockFSI &fsi) : fsi_(fsi) {}

    double ComputeNRLumpedMassMat(int pid, double sfi) const noexcept override;

    /**
     * @brief Nodal external force from one solid particle: body force on the solid plus
     *        added fluid mass, and the particle traction force.
     * @param pid Solid particle index.
     * @param sfi Shape-function value of the node at the particle.
     * @param ts Current time state; its load factor `facl` scales the body force.
     * @return Force contribution of this particle to the node.
     */
    std::array<double, 3> ComputeExternalForce(int pid, double sfi, const TimeState &ts) const noexcept override;

    /**
     * @brief Add solid inertia, added-mass inertia and FSI forces to the RHS.
     * @param mat System receiving the solid momentum residual.
     * @param naccel Nodal acceleration at the inertia evaluation time level.
     * @param offsets RHS offsets for the x, y and z momentum components.
     */
    void AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel, const std::vector<int> &offsets) override;
};

class MPMFEMBlockFSI {
  public:
    BoundaryCondition fsi_intf;

    // --- Sub-solvers (constructed with reference to this coordinator) ---
    FSIFluid fluid_;
    FSISolid solid_;

    // Block-iteration control
    const int max_block_iter = 100;
    const double phi_cut = 0.25e0;
    double relax_omega = 0.1e0;
    const double tol_ref = 1.0e-3;
    const double tol_abs = 1.0e-6;

    // -----------------------------------------------------------------
    // Constructor: pass *this to sub-solvers
    // -----------------------------------------------------------------
    MPMFEMBlockFSI() : fluid_(*this), solid_(*this) {}

    std::vector<double> added_mass;
    std::vector<double> nfsi_force;

    /**
     * @brief Read and initialize the coupled FSI case input data.
     *
     * Loads shared global parameters, then dispatches fluid and solid
     * boundary-condition/particle input to the corresponding sub-solvers.
     */
    void DataInput();

    void DetectFSIInterface();

    void UpdateFSIInterfaceBC();

    /**
     * @brief Advance one time step with partitioned (block-iterative) fluid-solid coupling.
     * @param ts Current time state of this step (load factor and physical time).
     */
    void SolveFSISystem(const TimeState &ts);

    void CalFSIResidual(double &rtr_ref, double &rtr_dof, const std::vector<double> &nvel_k);

    void CalFSIForce();

    /**
     * @brief Specifically for Turek CFD problem (Cylinder flow)
     */
    void CalDragLiftCoeffForTurekCFD();

  private:
    /**
     * @brief Check fluid-solid interface velocity agreement for the current block iteration.
     *
     * Aborts the run when the iteration reaches `max_block_iter` without converging.
     * @param block_it Current block iteration index; iteration 0 sets the reference residual.
     * @param nvel_k Solid nodal velocity at n+1 from this block iteration.
     * @param r0r_ref Interface velocity mismatch norm of the first iteration; written when `block_it == 0`.
     * @return Whether the interface velocity mismatch is within the relative or absolute tolerance.
     */
    bool CheckBlockConvergence(int block_it, const std::vector<double> &nvel_k, double &r0r_ref);
};

} // namespace mpmfemblockfsi
