#pragma once

#include <array>
#include <vector>

#include "module/bc.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

namespace mpm_mpm_monolithic_fsi {

class MPMMPMMonolithicFSI : public MaterialPoint {
  public:
    static constexpr int max_NR_it = 100;

    stabilizedmpm::StabilizedMPM fluid_;
    implicitmpm::ImplicitSolidMPM solid_;
    CrsMat fsi_sys;

    std::vector<double> nlm_lump_local, nlm_lump;
    std::vector<double> nlambda; // nodal multiplier force density
    std::vector<char> fixed_dof; // constrained and inactive components
    std::vector<double> column_scale; // linear increment scaling factors

    MPMMPMMonolithicFSI() {
        this->fluid_.blocks = {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18};
        this->solid_.blocks = {19, 20, 21, 23, 24, 25, 27, 28, 29};
        this->fluid_.rhs_start = 0;
        this->solid_.rhs_start = 4;
        this->ode_order = 2;
        this->fsi_sys.ndof = 10;
        this->fsi_sys.block_row = {0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,
                                          4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,8,8,9,9};
        this->fsi_sys.block_col = {0,1,2,3,7, 0,1,2,3,8, 0,1,2,3,9, 0,1,2,3,
                                          4,5,6,7, 4,5,6,8, 4,5,6,9, 0,4,1,5,2,6};
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

    /** @brief Integrate the paired fluid-solid phase gradient into common nodal interface-area weights. */
    void LumpedLagrangeMultiplier();

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

    void AssembleInterfaceSystem(const std::vector<double> &nvel_f, //
                                 const std::vector<double> &nvel_s);

    /** @brief Solve the endpoint velocity constraint and check endpoint velocity continuity. */
    void SolveFSISystem();

    /**
     * @brief Apply residual constraints and solve for physical monolithic increments.
     * @param NR_it Current Newton iteration used for preconditioner setup.
     * @return Number of Krylov iterations used by the coupled linear solve.
     */
    int SolveSystem(int NR_it);

    /**
     * @brief Check field residuals and endpoint interface velocity continuity.
     * @param nvel_f Fluid nodal velocity at the current endpoint iterate.
     * @param nvel_s Solid nodal velocity at the current endpoint iterate.
     * @param initial_norm Initial field residuals used for relative convergence.
     * @param NR_it Current Newton iteration.
     * @param linear_iterations Accumulated Krylov iterations for this time step.
     * @return Whether all field residuals and the velocity jump meet tolerance.
     */
    bool CheckNRConvergence(const std::vector<double> &nvel_f,
                            const std::vector<double> &nvel_s,
                            std::array<double, 4> &initial_norm,
                            int NR_it, int linear_iterations);

    void UpdateNRIncrement() override;

    /** @brief Select active field components and apply physical boundary constraints. */
    void BuildActiveDOFs();

    /** @brief Equilibrate field diagonals and multiplier couplings before the linear solve. */
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
     * @brief Configure nested field and multiplier Schur preconditioning.
     * @param mat Coupled linear system being solved.
     * @param pc PETSc preconditioner for the coupled Krylov solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

    /** @brief Release the fluid and solid response solvers. */
    ~MPMMPMMonolithicFSI();
};

} // namespace mpm_mpm_monolithic_fsi
