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
    MaterialPoint ifp; // --- Inflow particle ---

    enum class StabCoeff { PSPG, VMS };
    static constexpr StabCoeff stab_coeff = StabCoeff::VMS;

    std::vector<double> tau1, tau2;
    CrsMat NS_;

    StabilizedMPM() {
        this->ode_order = 2;
        this->NS_.ndof = 4;
        this->NS_.FEM_flag = false;
        this->NS_.use_petsc = true;
        this->NS_.use_fieldsplit = true;
        this->NS_.pressure_pc_use_amg = true;
        this->NS_.amg_rebuild_freq = 1; // rebuild fluid AMG every step
        this->NS_.owner_ = this;
    }

    /**
     * @brief Read fluid-specific parameters and initialize fluid state.
     */
    void DataInput();

    /**
     * @brief Read fluid boundary-condition data from an input stream.
     * @param infile Input stream positioned at the fluid boundary-condition section.
     */
    void InputBCData(std::ifstream &infile) override;
    /**
     * @brief Initialize fluid particle state (mass, volume, pressure) after input is read.
     */
    void InitializePointData() override;

    // --- User's responsibility ---
    void BuildPetscBCList(CrsMat &mat) override {
        mat.AddBCComponent(this->ubc, nuc);
        mat.AddBCComponent(this->vbc, nvc);
        mat.AddBCComponent(this->wbc, nwc);
        mat.AddBCComponent(this->pbc, npc);

        return;
    };

    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(nuc, rr);
        this->vbc.BCSetZero(nvc, rr);
        this->wbc.BCSetZero(nwc, rr);
        this->pbc.BCSetZero(npc, rr);

        return;
    };

    void BCNRSet() override {
        this->ubc.BCSetDt(nuc, this->ndispl);
        this->vbc.BCSetDt(nvc, this->ndispl);
        this->wbc.BCSetDt(nwc, this->ndispl);
        this->pbc.BCSetVal(0, this->npres);

        return;
    };
    // ------------------------------

    /**
     * @brief Map fluid particle mass/momentum to control points (P2G).
     */
    void Particle2Node() override;

    /**
     * @brief Map updated nodal velocity/pressure back to fluid particles (G2P).
     */
    void Node2Particle() override;

    /**
     * @brief Compute VMS/PSPG stabilization coefficients `tau1` and `tau2` for MPM fluid.
     */
    void MakNSStabCoeff();

    /**
     * @brief Assemble the stabilized Navier–Stokes MPM matrix and RHS.
     * @param naccel_k Nodal acceleration vector at the intermediate time level.
     * @param nvel_k   Nodal velocity vector at the intermediate time level.
     */
    void AssembleNSSystem(const std::vector<double> &naccel_k, //
                          const std::vector<double> &nvel_k);

    /**
     * @brief Solve the stabilized Navier–Stokes MPM system and update particle velocity/pressure.
     */
    void SolveNS();

    /**
     * @brief Apply the converged Newton–Raphson increment to nodal variables.
     */
    void UpdateNRIncrement() override;

    // --- For data IO ---
    /**
     * @brief Read fluid point data from input stream.
     * @param inflie Input file stream positioned at the point-data section.
     */
    void InputPointData(std::ifstream &inflie) override;

    /**
     * @brief Read fluid restart data from per-rank `*_re.txt` files.
     */
    void RestartInput() override;

    /**
     * @brief Write fluid particle data to VTK HDF5 visualization files.
     * @param iview Output view index.
     * @param istep Current time step.
     */
    void OutputPointDataVTKHDF(int iview, int istep) override;

    /**
     * @brief Write fluid restart data to per-rank `*_re.txt` files.
     */
    void RestartOutput() override;

    // --- MPI Particle move ---
    /**
     * @brief Move fluid particles that have crossed rank boundaries and exchange them via MPI.
     */
    void Moveparticle() override;

    // --- Inflow Particles ---
    /**
     * @brief Inject inflow particles at inflow boundaries for the current step.
     */
    void InflowParticles() override;

    /**
     * @brief Generate new inflow particles in direction `dir` using the inflow particle template.
     * @param dir  Inflow direction index (0=x, 1=y, 2=z).
     * @param ifp  Inflow particle template.
     * @param ifbc Inflow boundary condition.
     */
    void GenerateInflowParticles(int dir, MaterialPoint &ifp, //
                                 const BoundaryCondition &ifbc) override;

    // void Cp2NodeVTK();
};
} // namespace stabilizedmpm
