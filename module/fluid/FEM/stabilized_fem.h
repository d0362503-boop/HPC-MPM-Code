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
    // Liquid density and viscosity
    double rhol, rmul;
    // Gas density and viscosity
    double rhog, rmug;
    // Free-surface height
    double fs_height;
    // SUPG/PSPG and LSIC stabilization coefficients
    std::vector<double> tau1, tau2;
    // Element density and viscosity
    std::vector<double> rhoe, rmue;
    // Navier-Stokes system matrix
    CrsMat NS_;
    // Phase-field system matrix
    CrsMat PF_;

    StabilizedFEM() {
        this->ode_order = 1;
        this->NS_.ndof = 4;
        this->PF_.ndof = 1;
        this->NS_.FEM_flag = true;
        this->PF_.FEM_flag = true;
        this->NS_.use_petsc = true;
        this->NS_.use_schur_fieldsplit = true;
        this->PF_.use_petsc = true;
        this->NS_.amg_rebuild_freq = 20; // rebuild AMG every 20 steps
        this->PF_.amg_rebuild_freq = 20;
        this->NS_.owner_ = this;
        this->PF_.owner_ = this;
    }

    /**
     * @brief Read fluid boundary-condition data from input stream.
     * @param infile Input stream.
     */
    void InputBCData(std::ifstream &infile) override;

    /** @brief Allocate and zero fluid nodal fields. */
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

    /** @brief Solve the stabilized Navier-Stokes system. */
    void SolveNS();

    /** @brief Advance nodal velocity/pressure history arrays. */
    void UpdateNodalVar();

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
     * @param istep Current time step.
     */
    void OutputMeshDataVTKHDF(int iview, int istep);

  protected:
    /** @brief Apply velocity and pressure boundary conditions. */
    virtual void BCSet() {
        this->ApplyVelocityBC(this->nvel);
        this->ApplyAccelerationBC(this->naccel);
        this->pbc.BCSetVal(0, this->npres);

        return;
    };

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

  private:
    /**
     * @brief Generalized-alpha advection velocity.
     * @return Advection velocity vector.
     */
    std::vector<double> ComputeAdvectionVel();

    /**
     * @brief Compute SUPG/PSPG/LSIC stabilization coefficients.
     * @param adv_vel Advection velocity vector.
     */
    void MakNSStabCoeff(const std::vector<double> &adv_vel);

    /**
     * @brief Assemble stabilized Navier-Stokes matrix and RHS.
     * @param adv_vel Advection velocity vector.
     */
    void AssembleNSSystem(const std::vector<double> &adv_vel);
};

} // namespace stabilizedfem
