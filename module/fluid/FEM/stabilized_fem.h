#pragma once

#include <cmath>
#include <vector>

#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../solver/crsmat.h"

namespace stabilizedfem {

class StabilizedFEM : public MaterialPoint {
  public:
    double rhol, rmul, rhog, rmug;
    double fs_height;
    std::vector<double> tau1, tau2;
    std::vector<double> rhoe, rmue;
    CrsMat NS_;
    CrsMat PF_;

    StabilizedFEM() {
        this->ode_order = 1;
        this->NS_.ndof = 4;
        this->PF_.ndof = 1;
        this->NS_.FEM_flag = true;
        this->PF_.FEM_flag = true;
        this->NS_.use_petsc = true;
        this->PF_.use_petsc = true;
        this->NS_.use_fieldsplit = true;
        this->NS_.pressure_pc_use_amg = true;
        this->NS_.amg_rebuild_freq = 20; // rebuild fluid AMG every 20 steps
        this->PF_.amg_rebuild_freq = 20;
        this->NS_.owner_ = this;
        this->PF_.owner_ = this;
    }

    void InitializeMeshData() {
        VectorAssign(nodec * 3, this->nvel);
        VectorAssign(nodec * 3, this->nvel_old);
        VectorAssign(nodec * 3, this->nvel_older);
        VectorAssign(nodec * 3, this->naccel);
        VectorAssign(nodec, this->npres);
        VectorAssign(nodec, this->npres_old);
        VectorAssign(nodec, this->nphi);

        return;
    };

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

    virtual void BCSet() {
        this->ubc.BCSetVal(nuc, this->nvel);
        this->vbc.BCSetVal(nvc, this->nvel);
        this->wbc.BCSetVal(nwc, this->nvel);
        this->ubc.BCSetZero(nuc, this->naccel);
        this->vbc.BCSetZero(nvc, this->naccel);
        this->wbc.BCSetZero(nwc, this->naccel);
        this->pbc.BCSetVal(0, this->npres);

        return;
    };
    // ------------------------------

    std::vector<double> ComputeAdvectionVel();

    void MakNSStabCoeff(const std::vector<double> &adv_vel);

    void AssembleNSSystem(const std::vector<double> &adv_vel);

    void SolveNS();

    void UpdateNodalVar();

    // --- For PhaseField ---
    void Particle2NodePhi();

    double CalLiquidVol();

    void SetPFDomain();

    // --- For data IO ---
    void RestartInput() override;

    void RestartOutput() override;

    void Cp2NodeVTK();

    void OutputMeshDataVTKHDF(int iview, int istep);
};
} // namespace stabilizedfem
