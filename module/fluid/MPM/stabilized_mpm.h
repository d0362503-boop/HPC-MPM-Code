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

    StabilizedMPM() {
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

    /** @brief Read fluid parameters and initialize state. */
    void DataInput();

    /**
     * @brief Read fluid boundary-condition data from input stream.
     * @param infile Input stream.
     */
    void InputBCData(std::ifstream &infile) override;

    /** @brief Rebuild fluid control-point and inflow BCs for the current DLB region. */
    void RebuildBC() override;

    /** @brief Initialize fluid particle state. */
    void InitializePointData() override;

    /** @brief Map fluid particle mass/momentum to control points (P2G). */
    void Particle2Node() override;

    /** @brief Map nodal velocity/pressure back to fluid particles (G2P). */
    void Node2Particle() override;

    /** @brief Solve the stabilized Navier-Stokes system. */
    void SolveNS();

    // Newton steps driven by SolveNewton().
    /** @brief Freeze stabilization coefficients and reset the displacement, pressure and Newton scratch state. */
    void BeginNewtonStep() override;

    /** @brief Apply iterate BCs and kinematics, then assemble the Navier-Stokes tangent and residual. */
    void AssembleNewtonSystem() override;

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
     * @return Krylov iterations used.
     */
    int SolveNewtonCorrection(int NR_it) override;

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
        return {nodec * this->rhs_start, nodec * (this->rhs_start + 1), nodec * (this->rhs_start + 2),
                nodec * (this->rhs_start + 3)};
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
     */
    void OutputPointDataVTKHDF(int iview) override;

    /** @brief Write fluid restart data. */
    void RestartOutput() override;

    /** @brief Move particles across MPI rank boundaries. */
    void MoveParticle() override;

    /** @brief Apply DLB and rebuild the fluid Navier-Stokes matrix structure. */
    void ApplyDLB() override {

        MaterialPoint::ApplyDLB();

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
        mat.AddBCComponent(this->ubc, nuc);
        mat.AddBCComponent(this->vbc, nvc);
        mat.AddBCComponent(this->wbc, nwc);
        mat.AddBCComponent(this->pbc, npc);

        return;
    };

    /**
     * @brief Zero constrained DOFs in residual vector.
     * @param rr Residual vector.
     */
    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(nuc, rr);
        this->vbc.BCSetZero(nvc, rr);
        this->wbc.BCSetZero(nwc, rr);
        this->pbc.BCSetZero(npc, rr);

        return;
    };

    /** @brief Apply Dirichlet displacement increments for NR iteration. */
    void BCNRSet() {
        this->ApplyDisplacementBC(this->ndispl);
        this->pbc.BCSetVal(0, this->npres);

        return;
    };

    //   private:
    /**
     * @brief Assemble stabilized Navier-Stokes matrix and RHS.
     * @param mat Target system, cleared by the caller before assembly.
     * @param naccel_k Nodal acceleration at the current endpoint iterate.
     * @param nvel_k Nodal velocity at the current endpoint iterate.
     */
    virtual void AssembleSystem(CrsMat &mat, const std::vector<double> &nvel_k, //
                                const std::vector<double> &naccel_k);

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

  private:
    // Newton scratch state of the current step.
    std::vector<double> nvel_k_, naccel_k_; // nodal velocity and acceleration at the iterate
    double r0r_;                            // initial residual reference
};

} // namespace stabilizedmpm
