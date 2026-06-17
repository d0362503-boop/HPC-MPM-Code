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

    /**
     * @brief Compute the generalized-α advection velocity from `nvel_old` and `nvel_older`.
     * @return Advection velocity vector sized to `nodec * 3`.
     */
    std::vector<double> ComputeAdvectionVel();

    /**
     * @brief Compute SUPG/PSPG stabilization coefficients `tau1` and LSIC coefficients `tau2`.
     * @param adv_vel Advection velocity vector (size `nodec * 3`).
     */
    void MakNSStabCoeff(const std::vector<double> &adv_vel);

    /**
     * @brief Assemble the stabilized Navier–Stokes FEM matrix and RHS.
     * @param adv_vel Advection velocity vector used in the convective term.
     */
    void AssembleNSSystem(const std::vector<double> &adv_vel);

    /**
     * @brief Solve the stabilized Navier–Stokes system and update nodal velocity/pressure.
     */
    void SolveNS();

    /**
     * @brief Update time-advanced nodal variables (`nvel_old`, `npres_old`, etc.).
     */
    void UpdateNodalVar();

    // --- For PhaseField ---
    /**
     * @brief Project the liquid/gas indicator from particles to control points.
     */
    void Particle2NodePhi();

    /**
     * @brief Compute the total liquid volume from nodal volume fraction.
     * @return Total liquid volume.
     */
    double CalLiquidVol();

    /**
     * @brief Set element-wise density and viscosity based on the nodal volume fraction.
     */
    void SetPFDomain();

    // --- For data IO ---
    /**
     * @brief Read fluid restart data from per-rank `*_re.txt` files.
     */
    void RestartInput() override;

    /**
     * @brief Write fluid restart data to per-rank `*_re.txt` files.
     */
    void RestartOutput() override;

    /**
     * @brief Interpolate control-point fields to visualization nodes for VTK output.
     */
    void Cp2NodeVTK();

    /**
     * @brief Write fluid mesh data to VTK HDF5 visualization files.
     * @param iview Output view index.
     * @param istep Current time step.
     */
    void OutputMeshDataVTKHDF(int iview, int istep);
};
} // namespace stabilizedfem
