#pragma once

#include <cmath>
#include <vector>

#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../solver/crsmat.h"

namespace stabilizedmpm {

class StabilizedMPM : public MaterialPoint {
  public:
    // Inflow particle buffer
    MaterialPoint ifp;

    enum class StabCoeff { PSPG, VMS };
    // Stabilization-coefficient selector
    static constexpr StabCoeff stab_coeff = StabCoeff::VMS;

    // SUPG/PSPG and LSIC stabilization coefficients
    std::vector<double> tau1, tau2;
    // Navier-Stokes system matrix
    CrsMat NS_;

    StabilizedMPM() {
        this->do_dlb = false;
        this->gamma_nb = 1.0e0;
        this->beta_nb = 0.5e0;
        this->ode_order = 2;
        this->NS_.ndof = 4;
        this->NS_.FEM_flag = false;
        this->NS_.use_petsc = true;
        this->NS_.use_schur_fieldsplit = true;
        this->NS_.amg_rebuild_freq = 1; // rebuild AMG every step
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

    /**
     * @brief Read fluid point data from input stream.
     * @param inflie Input stream.
     */
    void InputPointData(std::ifstream &inflie) override;

    /** @brief Read fluid restart data. */
    void RestartInput() override;

    /**
     * @brief Write fluid particle data to VTK HDF5 files.
     * @param iview Output view index.
     * @param istep Current time step.
     */
    void OutputPointDataVTKHDF(int iview, int istep) override;

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

  private:
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
    void BCNRSet() override {
        this->ubc.BCSetDt(nuc, this->ndispl);
        this->vbc.BCSetDt(nvc, this->ndispl);
        this->wbc.BCSetDt(nwc, this->ndispl);
        this->pbc.BCSetVal(0, this->npres);

        return;
    };

    /** @brief Compute VMS/PSPG stabilization coefficients.
     * @param nvel_k   Nodal velocity vector.
     */
    void MakNSStabCoeff(const std::vector<double> &nvel_k);

    /**
     * @brief Assemble stabilized Navier-Stokes matrix and RHS.
     * @param naccel_k Nodal acceleration vector.
     * @param nvel_k   Nodal velocity vector.
     */
    void AssembleNSSystem(const std::vector<double> &nvel_k, //
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
};

} // namespace stabilizedmpm
