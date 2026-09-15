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
    const double tol_ref = 1.0e-5;
    const double tol_abs = 1.0e-10; // velocity tolerance in m/s

    // -----------------------------------------------------------------
    // Constructor: pass *this to sub-solvers
    // -----------------------------------------------------------------
    MPMMPMBlockFSI() : fluid_(*this), solid_(*this) {}

    /** @brief Release the fluid response Krylov workspace. */
    ~MPMMPMBlockFSI();

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
     * @brief Solve the endpoint velocity Schur system and update the nodal multiplier.
     * @param block_it Current fluid-solid block iteration.
     * @param fluid_velocity Fluid endpoint velocity used to form the interface residual.
     * @param solid_velocity Solid endpoint velocity used to form the interface residual.
     */
    void UpdateFSIMultiplier(int block_it, const std::vector<double> &fluid_velocity,
                             const std::vector<double> &solid_velocity);

  private:
    Mat schur_fluid_coupling_ = nullptr;        // fluid multiplier load map
    Mat schur_solid_coupling_ = nullptr;        // solid multiplier load map
    Vec schur_fluid_rhs_ = nullptr;              // fluid response load
    Vec schur_fluid_solution_ = nullptr;         // fluid displacement response
    Vec schur_fluid_defect_ = nullptr;           // approximate response residual
    Vec schur_fluid_correction_ = nullptr;       // approximate response correction
    Vec schur_solid_rhs_ = nullptr;              // solid response load
    Vec schur_solid_solution_ = nullptr;         // solid displacement response
    Vec schur_solid_response_ = nullptr;         // solid interface response
    KSP schur_solid_ksp_ = nullptr;             // parallel solid response solver
    KSP schur_fluid_ksp_ = nullptr;             // fluid response Krylov solver
    PetscInt schur_fluid_iterations_ = 0;       // fluid response iteration count
    PetscInt schur_solid_iterations_ = 0;       // solid response iteration count

    /**
     * @brief Forward a PETSc shell-matrix product to the exact interface response.
     * @param mat Exact Schur shell matrix carrying the FSI object as its context.
     * @param trial Trial interface multiplier increment.
     * @param response Exact endpoint velocity response.
     * @return PETSc success or error code.
     */
    static PetscErrorCode ApplyExactSchurShell(Mat mat, Vec trial, Vec response);

    /**
     * @brief Forward a PETSc shell-matrix product to the approximate interface response.
     * @param mat Approximate Schur shell matrix carrying the FSI object as its context.
     * @param trial Trial interface multiplier increment.
     * @param response Approximate endpoint velocity response.
     * @return PETSc success or error code.
     */
    static PetscErrorCode ApplyApproximateSchurShell(Mat mat, Vec trial, Vec response);

    /**
     * @brief Apply the nested approximate Schur solve used by the outer preconditioner.
     * @param pc Shell preconditioner carrying the approximate Schur solver as its context.
     * @param load Outer Schur residual to precondition.
     * @param solution Approximate multiplier correction.
     * @return PETSc success or error code.
     */
    static PetscErrorCode SolveApproximateSchur(PC pc, Vec load, Vec solution);

    /**
     * @brief Apply the exact Schur operator for the area-weighted endpoint velocity constraint.
     * @param trial_multiplier Trial increment of the nodal Lagrange multiplier.
     * @param response Endpoint velocity-constraint response produced by the fluid and solid solves.
     * @return PETSc success or error code.
     */
    PetscErrorCode ApplyExactSchur(Vec trial_multiplier, Vec response);

    /**
     * @brief Build parallel solid AMG for repeated Schur responses in this block iteration.
     */
    void BuildSolidResponse();

    /**
     * @brief Build the isolated fluid response solver and reuse it within the current time step.
     * @param block_it Current fluid-solid block iteration.
     */
    void BuildFluidResponse(int block_it);

    /**
     * @brief Apply an AMG approximation of the interface response for preconditioning.
     * @param trial Trial interface multiplier increment.
     * @param response Approximate endpoint velocity-constraint response.
     * @return PETSc success or error code.
     */
    PetscErrorCode ApplyApproximateSchur(Vec trial, Vec response);

    /**
     * @brief Solve a homogeneous-BC response with the fixed field tangent and preconditioner.
     *
     * Requires a completed field solve at the current iterate and interface loads
     * supported only on nodes active in both fields. Only the load and
     * response vectors change; Newton iteration and AMG-rebuild history remain untouched.
     * @param coupling Sparse map from multiplier to field load.
     * @param solver Frozen-tangent response solver, independent of physical Newton iterations.
     * @param trial_multiplier Trial multiplier increment on the interface.
     * @param load Field load produced by the trial multiplier.
     * @param solution Field displacement response to the trial load.
     * @param response Area-weighted interface displacement response.
     * @return KSP iterations used for this field response.
     */
    PetscInt SolveResponse(Mat coupling, KSP solver, Vec trial_multiplier, Vec load, Vec solution, Vec response);

    /**
     * @brief Add a converged multiplier increment to the local interface force field.
     * @param delta_multiplier Distributed increment of the interface multiplier.
     * @param local_interface_ids Global multiplier block assigned to each local interface node.
     */
    void AddMultiplierIncrement(Vec delta_multiplier, const std::vector<PetscInt> &local_interface_ids);

    /** @brief Release temporary field-response matrices, vectors, and the solid response solver. */
    void DestroySchurWorkspace();

    /**
     * @brief Build one field's multiplier load map and reusable response vectors.
     * @param mat Field tangent and parallel ownership map.
     * @param weights Nodal interface-area weights.
     * @param local_interface_ids Global multiplier index of each local interface node.
     * @param local_multiplier_dofs Multiplier unknowns owned by this rank.
     * @param global_multiplier_dofs Total multiplier unknowns.
     * @param coupling Sparse multiplier-to-field load map.
     * @param load Reusable field response load.
     * @param solution Reusable field displacement response.
     */
    void BuildCouplingOperator(CrsMat &mat, const std::vector<double> &weights,
                               const std::vector<PetscInt> &local_interface_ids, PetscInt local_multiplier_dofs,
                               PetscInt global_multiplier_dofs, Mat &coupling, Vec &load, Vec &solution);

    /**
     * @brief Assemble nodal fluid (4x4) and solid (3x3) tangent responses into a block-lumped Schur preconditioner.
     * @param interface_nodes Global natural node number of each interface multiplier block.
     * @param local_dofs Locally owned multiplier degrees of freedom.
     * @param global_dofs Global multiplier degrees of freedom.
     * @param preconditioner_mat Resulting 3x3-per-node endpoint velocity-Schur approximation.
     */
    void BuildLumpedSchurPreconditioner(const std::vector<PetscInt> &interface_nodes, PetscInt local_dofs,
                                        PetscInt global_dofs, Mat &preconditioner_mat);
};

/**
 * @brief Entry point of the MPM--MPM block-coupling FSI driver.
 */
void MPMBlockFSI();

} // namespace mpmmpmblockfsi
