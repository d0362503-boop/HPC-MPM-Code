#pragma once

#include <cmath>
#include <vector>

#include "module/dataset.h"
#include "module/implicit/implicit_mpm.h"
#include "module/implicit/newton_solver.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solver/crsmat.h"

namespace stabilizedmpm {

class StabilizedMPM : public MaterialPoint, public NonlinearProblem {
  public:
    static constexpr int max_NR_it = 100; // last Newton iteration index

    // Inflow particle buffer
    MaterialPoint ifp;

    enum class StabCoeff { PSPG, VMS };
    // Stabilization-coefficient selector
    static constexpr StabCoeff stab_coeff = StabCoeff::VMS;

    // Per-particle PSPG tau/rho and rho*tau_LSIC
    std::vector<double> tau_pspg, tau_lsic;

    CrsMat NS_;            // Navier-Stokes system matrix
    ImplicitMPM implicit_; // implicit time integration and kinematics

    std::vector<int> blocks; // target scalar block indices
    int rhs_start;           // first target component index

    /**
     * @brief Configure the MPM Navier-Stokes system and inflow buffer on this rank's grid.
     * @param grid   Grid of this rank (background mesh, control points and partition).
     * @param params Run parameters of this run.
     */
    StabilizedMPM(const Grid &grid, const RunParams &params)
        : MaterialPoint(grid, params), ifp(grid, params), NS_(grid, params) {
        this->blocks = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        this->rhs_start = 0;
        this->do_dlb = false;
        this->implicit_.integrator_.gamma_nb = 1.0e0;
        this->implicit_.integrator_.beta_nb = 0.5e0;
        this->implicit_.integrator_.ode_order = 2;
        this->NS_.ndof = 4;
        this->NS_.block_row = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
        this->NS_.block_col = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
        this->NS_.FEM_flag = false;
        this->NS_.use_petsc = true;
        this->NS_.rebuild_freq = 1; // refresh each step start
        this->NS_.owner_ = this;
    }

    /**
     * @brief Read fluid parameters and initialize state.
     * @param grid   Grid of this rank; mesh, control-point and partition data are read into it.
     * @param params Run parameters; filled from file.dat and the parameter file.
     */
    void DataInput(Grid &grid, RunParams &params);

    /**
     * @brief Read fluid boundary-condition data from input stream.
     * @param infile Input stream.
     */
    void InputBCData(std::ifstream &infile) override;

    /** @brief Rebuild fluid control-point and inflow BCs for the current DLB region. */
    void RebuildBC() override;

    /** @brief Initialize fluid particle state. */
    void InitializePointData() override;

    /**
     * @brief Map fluid particle mass/momentum to control points (P2G).
     * @param ts Current time state of this step (load factor and physical time).
     */
    void Particle2Node(const TimeState &ts) override;

    /**
     * @brief Commit the converged nodal kinematics and pressure, then map them back to
     *        fluid particles (G2P).
     * @param ts Current time state of this step (load factor and physical time).
     */
    void Node2Particle(const TimeState &ts) override;

    /** @brief Accept the converged pressure iterate as the nodal pressure: `npres <- npres_k`. */
    void CommitNodalPressure();

    std::vector<double> npres_k; // nodal pressure at the iterate

    /**
     * @brief Solve the stabilized Navier-Stokes system.
     * @param ts Current time state of this step (load factor and physical time).
     */
    void SolveNS(const TimeState &ts);

    // Newton steps driven by SolveNewton().
    /** @brief Freeze stabilization coefficients and reset the displacement, pressure and Newton scratch state. */
    void BeginNewtonStep() override;

    /**
     * @brief Apply iterate BCs and kinematics, then assemble the Navier-Stokes tangent and residual.
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

    /** @brief Release the Navier-Stokes linear-system storage. */
    void EndNewtonStep() override;

    /**
     * @brief Configure the velocity-pressure Schur field-split preconditioner.
     * @param mat Fluid linear system and field-split selection.
     * @param pc PETSc preconditioner associated with the fluid solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

    /**
     * @brief Get target RHS offsets for fluid momentum and pressure.
     * @return Current control-point offsets for x, y, z momentum and pressure.
     */
    std::vector<int> RHSOffsets() const {
        return {this->grid_.cmesh.nodec * this->rhs_start, this->grid_.cmesh.nodec * (this->rhs_start + 1), this->grid_.cmesh.nodec * (this->rhs_start + 2),
                this->grid_.cmesh.nodec * (this->rhs_start + 3)};
    }

    /**
     * @brief Compute VMS/PSPG stabilization coefficients.
     * @param nvel_k Nodal fluid velocity at the residual evaluation time level.
     */
    void MakeNSStabCoeff(const std::vector<double> &nvel_k);

    /**
     * @brief Read fluid point data from input stream.
     * @param infile Input stream.
     */
    void InputPointData(std::ifstream &infile) override;

    /** @brief Read fluid restart data. */
    void RestartInput() override;

    /**
     * @brief Write fluid particle data to VTK HDF5 files.
     * @param iview Output view index.
     * @param ts Current time state; its physical time is written with the data.
     */
    void OutputPointDataVTKHDF(int iview, const TimeState &ts) override;

    /** @brief Write fluid restart data. */
    void RestartOutput() override;

    /** @brief Move particles across MPI rank boundaries. */
    void MoveParticle() override;

    /**
     * @brief Apply DLB and rebuild the fluid Navier-Stokes matrix structure.
     * @param grid Writable view of the grid this object observes; repartitioned in place.
     */
    void ApplyDLB(Grid &grid) override {

        MaterialPoint::ApplyDLB(grid);

        this->NS_.BuildCrsMat(16);

        return;
    };

    /** @brief Migrate all fluid particle state using the prepared MPI communication plan. */
    void MigrateParticleData() override;

    //   protected:
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

    /**
     * @brief Apply Dirichlet displacement increments for NR iteration.
     * @param ts Current time state; its load factor `facl` scales the prescribed values.
     */
    void BCNRSet(const TimeState &ts) {
        this->ApplyDisplacementBC(this->ndispl, ts);
        this->pbc.BCSetVal(0, this->npres_k, ts.facl);

        return;
    };

    //   private:
    /**
     * @brief Assemble stabilized Navier-Stokes matrix and RHS.
     * @param mat Target system, cleared by the caller before assembly.
     * @param naccel_k Nodal acceleration at the current endpoint iterate.
     * @param nvel_k Nodal velocity at the current endpoint iterate.
     * @param ts Current time state; its load factor `facl` scales the body force.
     */
    virtual void AssembleSystem(CrsMat &mat, const std::vector<double> &nvel_k, //
                                const std::vector<double> &naccel_k, const TimeState &ts);

    /** @brief Apply converged NR increment to nodal variables. */
    void UpdateNRIncrement() override;

    /** @brief Generate inflow particles for all active directions. */
    void InflowParticles() override;

    /**
     * @brief Generate inflow particles in empty boundary cells.
     * @param dir   Inflow direction index.
     * @param ifp   Inflow particle buffer.
     * @param infbc Inflow boundary condition.
     */
    void GenerateInflowParticlesEmptyMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) override;

    /**
     * @brief Generate inflow particles by cloning existing particles.
     * @param dir   Inflow direction index.
     * @param ifp   Inflow particle buffer.
     * @param infbc Inflow boundary condition.
     */
    void GenerateInflowParticlesFilledMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) override;

  protected:
    // Newton scratch state of the current step.
    std::vector<double> nvel_k_, naccel_k_; // nodal velocity and acceleration at the iterate

  private:
    double r0r_; // initial residual reference
};

} // namespace stabilizedmpm
