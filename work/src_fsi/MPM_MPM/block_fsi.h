#pragma once

#include <array>
#include <vector>

#include "module/bc.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

namespace mpmmpmblockfsi {

class MPMMPMBlockFSI;

class FSIFluid : public stabilizedmpm::StabilizedMPM {
  public:
    MPMMPMBlockFSI &fsi_; // Reference back to the FSI coordinator

    std::vector<double> lm_lumped; // Lagrange multiplier lumped matrix

    FSIFluid(MPMMPMBlockFSI &fsi) : fsi_(fsi) {}

    /**
     * @brief Assemble the fluid tangent and residual with the volume multiplier.
     * @param nvel_k Nodal fluid velocity at the current nonlinear iterate.
     * @param naccel_k Nodal fluid acceleration at the current nonlinear iterate.
     */
    void AssembleSystem(const std::vector<double> &nvel_k, const std::vector<double> &naccel_k) override;
};

class FSISolid : public implicitmpm::ImplicitSolidMPM {
  public:
    MPMMPMBlockFSI &fsi_; // Reference back to the FSI coordinator

    std::vector<double> lm_lumped; // Lagrange multiplier lumped matrix

    FSISolid(MPMMPMBlockFSI &fsi) : fsi_(fsi) {}

    /**
     * @brief Assemble the solid tangent and residual with the volume multiplier.
     * @param naccel_k Nodal solid acceleration at the current nonlinear iterate.
     * @param nvel_k Nodal solid velocity at the current nonlinear iterate.
     * @param stress_k Particle stress at the current nonlinear iterate.
     */
    void AssembleSystem(const std::vector<double> &naccel_k, const std::vector<double> &nvel_k,
                        std::vector<std::array<double, 6>> &stress_k) override;
};

class MPMMPMBlockFSI {
  public:
    BoundaryCondition fsi_intf;

    // --- Sub-solvers (constructed with reference to this coordinator) ---
    FSIFluid fluid_;
    FSISolid solid_;

    std::vector<double> nfsi_force; // nodal multiplier force density

    // Block-iteration control
    const int max_block_iter = 100;
    const double tol_ref = 1.0e-6;
    const double tol_abs = 1.0e-10; // velocity tolerance in m/s

    // -----------------------------------------------------------------
    // Constructor: pass *this to sub-solvers
    // -----------------------------------------------------------------
    MPMMPMBlockFSI() : fluid_(*this), solid_(*this) {}

    /**
     * @brief Read and initialize the coupled FSI case input data.
     *
     * Loads shared global parameters, then dispatches fluid and solid
     * boundary-condition/particle input to the corresponding sub-solvers.
     */
    void DataInput();

    /** @brief Build common interface weights and select nodes active in both fields. */
    void DetectFSIInterface();

    /** @brief Integrate the paired fluid-solid phase gradient into common nodal interface-area weights. */
    void LumpedLagrangeMultiplier();

    /** @brief Solve the endpoint velocity constraint and check endpoint velocity continuity. */
    void SolveFSISystem();

    /**
     * @brief Measure the endpoint velocity jump on the interface.
     * @param fluid_velocity Fluid endpoint velocity reconstructed from the current displacement iterate.
     * @param solid_velocity Solid endpoint velocity reconstructed from the current displacement iterate.
     * @param rtr_ref Global Euclidean norm of the velocity jump.
     * @param rtr_dof Root-mean-square velocity jump per interface node.
     */
    void CalFSIResidual(const std::vector<double> &fluid_velocity, const std::vector<double> &solid_velocity,
                        double &rtr_ref, double &rtr_dof);

    /**
     * @brief Report displacement and velocity jumps on the current FSI interface.
     * @param fluid_velocity Fluid velocity reconstructed from the converged displacement.
     * @param solid_velocity Solid velocity reconstructed from the converged displacement.
     */
    void ReportFSIContinuity(const std::vector<double> &fluid_velocity, const std::vector<double> &solid_velocity);

    /**
     * @brief Solve the endpoint velocity Schur system and update the nodal multiplier.
     * @param block_it Current fluid-solid block iteration.
     * @param fluid_velocity Fluid endpoint velocity used to form the interface residual.
     * @param solid_velocity Solid endpoint velocity used to form the interface residual.
     */
    void UpdateFSIMultiplier(int block_it, const std::vector<double> &fluid_velocity,
                             const std::vector<double> &solid_velocity);

  private:
    std::vector<PetscInt> schur_interface_ids_; // local interface index map
    VecScatter schur_input_scatter_ = nullptr;  // multiplier gather operation
    Vec schur_input_ = nullptr;                 // gathered multiplier increment
    KSP schur_solid_ksp_ = nullptr;             // factored solid tangent response
    PetscInt schur_matvec_count_ = 0;           // Schur operator application count
    PetscInt schur_fluid_iterations_ = 0;       // fluid response iteration count
    PetscInt schur_solid_iterations_ = 0;       // solid response backsolve count

    /** @brief Report control-point pressure and multiplier data near the FSI interface. */
    void ReportFSIPressureProfile();

    /**
     * @brief Apply the exact Schur operator for the area-weighted endpoint velocity constraint.
     * @param trial_multiplier Trial increment of the nodal Lagrange multiplier.
     * @param response Endpoint velocity-constraint response produced by the fluid and solid solves.
     * @return PETSc success or error code.
     */
    PetscErrorCode ApplyExactSchur(Vec trial_multiplier, Vec response);

    /**
     * @brief Factor the constrained solid tangent once for repeated Schur responses in this block iteration.
     */
    void BuildSolidResponse();

    /**
     * @brief Solve a homogeneous-BC response with the fixed field tangent and preconditioner.
     *
     * Requires a completed field solve at the current iterate and interface loads
     * supported only on nodes active in both fields. Only the load and
     * response vectors change; Newton iteration and AMG-rebuild history remain untouched.
     * @param mat Field linear system containing the trial interface load.
     * @param ndof Number of nodal field unknowns, including fluid pressure.
     * @param solver Frozen-tangent response solver, independent of physical Newton iterations.
     * @param total_iterations Accumulated KSP iterations, one per solid backsolve.
     */
    void SolveResponse(CrsMat &mat, int ndof, KSP solver, PetscInt &total_iterations);

    /**
     * @brief Assemble nodal fluid (4x4) and solid (3x3) tangent responses into a block-lumped Schur preconditioner.
     * @param local_interface_ids Local interface-node map to global interface indices.
     * @param local_dofs Locally owned multiplier degrees of freedom.
     * @param global_dofs Global multiplier degrees of freedom.
     * @param preconditioner_mat Resulting 3x3-per-node alpha-f velocity-Schur approximation.
     */
    void BuildLumpedSchurPreconditioner(const std::vector<PetscInt> &local_interface_ids, PetscInt local_dofs,
                                        PetscInt global_dofs, Mat &preconditioner_mat);
};

/**
 * @brief Entry point of the MPM--MPM block-coupling FSI driver.
 */
void MPMBlockFSI();

} // namespace mpmmpmblockfsi
