#pragma once

#include <vector>

#include "module/dataset.h"
#include "module/implicit/generalized_alpha_integrator.h"
#include "module/implicit/newton_solver.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solver/crsmat.h"

namespace stabilizedfem {

class StabilizedFEM : public MaterialPoint, public NonlinearProblem {
  public:
    static constexpr int max_NR_it = 100;

    // Liquid density and viscosity
    double rhol, rmul;
    // Gas density and viscosity
    double rhog, rmug;
    // Free-surface height
    double fs_height;
    // Per-element SUPG/PSPG and LSIC coefficients
    std::vector<double> tau_supg, tau_lsic;
    // Element density and viscosity
    std::vector<double> rhoe, rmue;

    GeneralizedAlphaIntegrator integrator_; // time integration state
    CrsMat NS_;                             // Navier-Stokes system matrix
    CrsMat PF_;                             // Phase-field system matrix

    std::vector<int> blocks; // target scalar block indices
    int rhs_start;           // first target component index
    /**
     * @brief Configure Schur field splitting for NS and BoomerAMG for other systems.
     * @param mat Fluid or phase-field linear system being configured.
     * @param pc PETSc preconditioner associated with the system solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

    /**
     * @brief Configure the four-component Navier-Stokes system on this rank's grid.
     * @param grid   Grid of this rank (background mesh, control points and partition).
     * @param params Run parameters of this run.
     */
    StabilizedFEM(const Grid &grid, const RunParams &params)
        : MaterialPoint(grid, params), NS_(grid, params), PF_(grid, params) {
        this->blocks = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        this->rhs_start = 0;
        this->integrator_.ode_order = 1;
        this->NS_.ndof = 4;
        this->NS_.block_row = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
        this->NS_.block_col = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
        this->PF_.ndof = 1;
        this->PF_.block_row = {0};
        this->PF_.block_col = {0};
        this->NS_.FEM_flag = true;
        this->PF_.FEM_flag = true;
        this->NS_.use_petsc = true;
        this->PF_.use_petsc = true;
        this->NS_.rebuild_freq = 20; // scheduled PC refresh interval
        this->PF_.rebuild_freq = 20;
        this->NS_.owner_ = this;
        this->PF_.owner_ = this;
    }

    /**
     * @brief Read `file.dat`, the fluid parameter file and this rank's grid data for a
     *        stand-alone FEM fluid run (defined by the `work/src_fluid/FEM` driver).
     * @param grid   Grid of this rank; mesh, control-point and partition data are read into it.
     * @param params Run parameters; filled from file.dat and the parameter file.
     */
    void DataInput(Grid &grid, RunParams &params);

    /**
     * @brief Read fluid boundary-condition data from input stream.
     * @param infile Input stream.
     */
    void InputBCData(std::ifstream &infile) override;

    /** @brief Allocate and zero fluid nodal fields. */
    void InitializeMeshData() {
        VectorAssign(this->grid_.cmesh.nodec * 3, this->nvel);
        VectorAssign(this->grid_.cmesh.nodec * 3, this->naccel);
        VectorAssign(this->grid_.cmesh.nodec, this->npres);
        VectorAssign(this->grid_.cmesh.nodec, this->nphi);

        return;
    };

    /**
     * @brief Get target RHS offsets for fluid momentum and pressure.
     * @return Current control-point offsets for x, y, z momentum and pressure.
     */
    std::vector<int> RHSOffsets() const {
        return {this->grid_.cmesh.nodec * this->rhs_start, this->grid_.cmesh.nodec * (this->rhs_start + 1), this->grid_.cmesh.nodec * (this->rhs_start + 2),
                this->grid_.cmesh.nodec * (this->rhs_start + 3)};
    }

    /**
     * @brief Compute the n+1 nodal acceleration from the n+1 velocity iterate `nvel_k`,
     *        using the generalized-α relation v^{n+1} = v^n + Δt[(1-γ)a^n + γa^{n+1}].
     * @return Nodal acceleration at time level n+1, sized `nodec * 3`.
     */
    std::vector<double> ComputeNodeAccelFromVel() const {

        return this->integrator_.ComputeNodeAccelFromVel(this->nvel_k, this->nvel, this->naccel, this->params_.dt);
    }

    /**
     * @brief Solve the stabilized Navier-Stokes system for the n+1 iterate.
     *
     * Leaves `nvel`, `npres` and `naccel` at time level n; call
     * CommitNodalSolution() once the time step is accepted.
     * @param ts Current time state of this step (load factor and physical time).
     */
    void SolveNS(const TimeState &ts);

    // Newton steps driven by SolveNewton().
    /**
     * @brief Initialize the n+1 iterate with the same-velocity predictor
     *        (`nvel_k = nvel`, `npres_k = npres`) and reset the residual reference.
     */
    void BeginNewtonStep() override;

    /**
     * @brief Apply BCs to the n+1 iterate, update the n+1 acceleration and the
     *        SUPG/PSPG/LSIC coefficients at n+α_f, then assemble the Picard tangent
     *        and residual (derivatives of τ are neglected).
     * @param ts Current time state; its load factor `facl` scales BCs and body forces.
     */
    void AssembleNewtonSystem(const TimeState &ts) override;

    /**
     * @brief Check the Navier-Stokes residual.
     * @param NR_it Current Newton iteration.
     * @param ksp_its Krylov iterations of the previous correction.
     * @return Whether the Newton iteration has converged.
     */
    bool CheckNewtonConvergence(int NR_it, int ksp_its) override;

    /**
     * @brief Solve the Navier-Stokes system for the velocity-pressure correction.
     * @param NR_it Current Newton iteration.
     * @param ts Current time state; its step index sets the preconditioner rebuild cadence.
     * @return Krylov iterations used.
     */
    int SolveNewtonCorrection(int NR_it, const TimeState &ts) override;

    /** @brief Add the solved velocity/pressure correction to the n+1 iterate `nvel_k`, `npres_k`. */
    void UpdateNRIncrement() override;

    /** @brief Release the Navier-Stokes linear-system storage. */
    void EndNewtonStep() override;

    /**
     * @brief Accept the converged n+1 iterate as the new time level n.
     *
     * Updates `naccel` from the generalized-α velocity relation, then copies
     * `nvel_k` to `nvel` and `npres_k` to `npres`.
     */
    void CommitNodalSolution();

    /** @brief Project particle liquid/gas indicator to control points. */
    void Particle2NodePhi();

    /**
     * @brief Compute total liquid volume from nodal volume fraction.
     * @return Total liquid volume.
     */
    double CalLiquidVol();

    /** @brief Set element density and viscosity from nodal volume fraction. */
    void SetPFDomain();

    /** @brief Read fluid restart data. */
    void RestartInput() override;

    /** @brief Write fluid restart data. */
    void RestartOutput() override;

    /** @brief Interpolate control-point fields to visualization nodes. */
    void Cp2NodeVTK();

    /**
     * @brief Write fluid mesh data to VTK HDF5 files.
     * @param iview Output view index.
     * @param ts Current time state; its physical time is written with the data.
     */
    void OutputMeshDataVTKHDF(int iview, const TimeState &ts);

    std::vector<double> nvel_k, naccel_k, npres_k; // n+1 Newton iterate

  protected:
    double r0r = 0.0e0;

    /**
     * @brief Apply velocity and pressure boundary conditions to the n+1 iterate.
     * @param ts Current time state; its load factor `facl` scales the prescribed values.
     */
    virtual void BCNRSet(const TimeState &ts) {
        this->ApplyVelocityBC(this->nvel_k, ts);
        this->pbc.BCSetVal(0, this->npres_k, ts.facl);

        return;
    };

    /**
     * @brief Register constrained DOFs with PETSc matrix.
     * @param mat PETSc matrix.
     */
    void BuildPetscBCList(CrsMat &mat) override {
        mat.AddBCComponent(this->ubc, this->grid_.cmesh.nuc);
        mat.AddBCComponent(this->vbc, this->grid_.cmesh.nvc);
        mat.AddBCComponent(this->wbc, this->grid_.cmesh.nwc);
        mat.AddBCComponent(this->pbc, this->grid_.cmesh.npc);

        return;
    };

    /**
     * @brief Zero constrained DOFs in residual vector.
     * @param rr Residual vector.
     */
    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(this->grid_.cmesh.nuc, rr);
        this->vbc.BCSetZero(this->grid_.cmesh.nvc, rr);
        this->wbc.BCSetZero(this->grid_.cmesh.nwc, rr);
        this->pbc.BCSetZero(this->grid_.cmesh.npc, rr);

        return;
    };

    //   private:
    /**
     * @brief Compute the advection velocity at n+α_f from `nvel` (n) and `nvel_k` (n+1 iterate).
     * @return Nodal advection velocity v^{n+α_f}, sized `nodec * 3`.
     */
    std::vector<double> ComputeAdvectionVel();

    /**
     * @brief Compute per-element SUPG/PSPG (`tau_supg`) and LSIC (`tau_lsic`) coefficients.
     * @param adv_vel Nodal advection velocity at n+α_f, used for element size and |v|.
     */
    void MakeNSStabCoeff(const std::vector<double> &adv_vel);

    /**
     * @brief Assemble the stabilized Navier-Stokes Picard tangent and residual (RHS = -R).
     * @param mat Target system; fluid blocks and rows follow `blocks` and `RHSOffsets()`.
     * @param adv_vel Nodal advection velocity at n+α_f used in convection and SUPG weights.
     * @param ts Current time state; its load factor `facl` scales the body force.
     */
    virtual void AssembleSystem(CrsMat &mat, const std::vector<double> &adv_vel, const TimeState &ts);
};

} // namespace stabilizedfem
