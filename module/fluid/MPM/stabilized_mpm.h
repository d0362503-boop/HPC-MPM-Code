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
        this->NS_.amg_rebuild_freq = 1; // rebuild fluid AMG every 20 steps
        this->NS_.owner_ = this;
    }

    void DataInput();

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

    void Particle2Node() override;

    void Node2Particle() override;

    void MakNSStabCoeff();

    void AssembleNSSystem(const std::vector<double> &naccel_k, //
                          const std::vector<double> &nvel_k);

    void SolveNS();

    void UpdateNRIncrement() override;

    // --- For data IO ---
    void InputPointData(std::ifstream &inflie) override;

    void RestartInput() override;

    void OutputPointDataVTKHDF(int iview, int istep) override;

    void RestartOutput() override;

    // --- MPI Particle move ---
    void Moveparticle() override;

    // --- Inflow Particles ---
    void InflowParticles() override;

    void GenerateInflowParticles(int dir, MaterialPoint &ifp, //
                                 const BoundaryCondition &ifbc) override;

    // void Cp2NodeVTK();
};
} // namespace stabilizedmpm
