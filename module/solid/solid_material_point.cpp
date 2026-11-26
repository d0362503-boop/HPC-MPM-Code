#include "module/solid/solid_material_point.h"

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "module/DLB/mpm_dlb.h"
#include "module/bc.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

void SolidMaterialPointBase::InitializePointData() {
    VectorAssign(this->num, this->id);
    VectorAssign(this->num, this->matid);
    VectorAssign(this->num, this->surf_point);
    VectorAssign(this->num, this->mass);
    VectorAssign(this->num, this->vol0);
    VectorAssign(this->num, this->vol);
    VectorAssign(this->num, this->det_def_grad, 1.0e0);

    VectorAssign(this->num, this->coord);
    VectorAssign(this->num, this->vel);
    VectorAssign(this->num, this->accel);
    VectorAssign(this->num, this->trac_force);
    VectorAssign(this->num, this->stress);
    VectorAssign(this->num, this->def_grad);

    if (this->Fbar_flag) {
        VectorAssign(this->num, this->delta_def_grad_bar);
        VectorAssign(this->num, this->det_def_grad_bar, 1.0e0);
        VectorAssign(this->num, this->def_grad_bar);
    }
    if (this->solswitch == MapScheme::TPIC) {
        VectorAssign(this->num, this->tpic.vel_grad);
        VectorAssign(this->num, this->tpic.accel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        VectorAssign(this->num, this->apic.vel_Bmat);
        VectorAssign(this->num, this->apic.accel_Bmat);
        VectorAssign(this->num, this->apic.inv_Dmat);
    }

    return;
}

void SolidMaterialPointBase::RebuildBC() {

    this->ubc.RebuildLocalControlPointBC(this->grid_);
    this->vbc.RebuildLocalControlPointBC(this->grid_);
    this->wbc.RebuildLocalControlPointBC(this->grid_);

    return;
}

void SolidMaterialPointBase::MoveParticle() {
    const Partition &part = this->grid_.part;

    this->DetermineParticleRank(part.regions);

    if (this->par_comm_.nmps + this->par_comm_.nrps + this->par_comm_.nrmp == 0) return;

    this->num += this->par_comm_.nrps - this->par_comm_.nmps - this->par_comm_.nrmp;

    this->MigrateParticleData();

    return;
}

void SolidMaterialPointBase::MigrateParticleData() {

    // ---- Scalar part ----
    this->par_comm_.PointVarComm(this->num, this->id);           // ---- Point Global ID ----
    this->par_comm_.PointVarComm(this->num, this->matid);        // ---- Point material ID ----
    this->par_comm_.PointVarComm(this->num, this->surf_point);   // ---- Point surface_ID ----
    this->par_comm_.PointVarComm(this->num, this->mass);         // ---- Point mass ----
    this->par_comm_.PointVarComm(this->num, this->vol0);         // ---- Point initial volume ----
    this->par_comm_.PointVarComm(this->num, this->vol);          // ---- Point volume ----
    this->par_comm_.PointVarComm(this->num, this->det_def_grad); // ---- Point determinant of deformation gradient ----

    // ---- Vector part ----
    this->par_comm_.PointVarComm(this->num, this->coord);      // ---- Point coordinate ----
    this->par_comm_.PointVarComm(this->num, this->vel);        // ---- Point velocity ----
    this->par_comm_.PointVarComm(this->num, this->accel);      // ---- Point acceleration ---
    this->par_comm_.PointVarComm(this->num, this->trac_force); // ---- Point traction force ----
    this->par_comm_.PointVarComm(this->num, this->stress);     // ---- Point Cauchy stress ----

    // ---- Matrix part ----
    this->par_comm_.PointVarComm(this->num, this->def_grad); // --- Point deformation gradient ---

    if (Fbar_flag) {
        // --- Point determinant of deformation gradient bar ---
        this->par_comm_.PointVarComm(this->num, this->det_def_grad_bar);
        // --- Point delta deformation gradient bar ---
        this->par_comm_.PointVarComm(this->num, this->delta_def_grad_bar);
        // --- Point deformation gradient bar ---
        this->par_comm_.PointVarComm(this->num, this->def_grad_bar);
    }

    if (this->solswitch == MapScheme::TPIC) {
        // ---- TPIC acceleration gradient ----
        this->par_comm_.PointVarComm(this->num, this->tpic.accel_grad);
        // ---- TPIC velocity gradient ----
        this->par_comm_.PointVarComm(this->num, this->tpic.vel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        // ---- APIC acceleration B matrix ----
        this->par_comm_.PointVarComm(this->num, this->apic.accel_Bmat);
        // ---- APIC velocity B matrix ----
        this->par_comm_.PointVarComm(this->num, this->apic.vel_Bmat);
        // ---- APIC inverse D matrix ----
        this->par_comm_.PointVarComm(this->num, this->apic.inv_Dmat);
    }

    return;
}

void SolidMaterialPointBase::DetermineRigidBC() {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(cmesh.nodec, this->nvof);
    VectorAssign(cmesh.nodec * 3, this->nrigid_normal);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            if (this->materials_.iprop[this->matid[pid]] == -1) {
                std::array<double, 3> xyp = {this->coord[pid][0], this->coord[pid][1], this->coord[pid][2]};
                MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);
                for (int ni = 0; ni < nenode; ni++) {
                    int nid = ncm[ni];
                    double sfi = sf[ni];
                    double dsfi1 = dsf[ni][0];
                    double dsfi2 = dsf[ni][1];
                    double dsfi3 = dsf[ni][2];
                    this->nvof[nid] += sfi * this->vol[pid];
                    this->nrigid_normal[nid + cmesh.nuc] += dsfi1 * this->mass[pid];
                    this->nrigid_normal[nid + cmesh.nvc] += dsfi2 * this->mass[pid];
                    this->nrigid_normal[nid + cmesh.nwc] += dsfi3 * this->mass[pid];
                }
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->nrigid_normal, {cmesh.nuc, cmesh.nvc, cmesh.nwc});

    VectorAssign(cmesh.nodec, this->rigid_bc.nbc);
    VectorAssign(cmesh.nodec, this->rigid_bc.fbc);
    for (int n = 0; n < cmesh.nodec; n++) {
        if (this->nvof[n] > 0.0e0) {
            this->rigid_bc.nbc[this->rigid_bc.ibc] = n;
            this->rigid_bc.fbc[this->rigid_bc.ibc] = 0.0e0;
            this->rigid_bc.ibc++;
        }
        const double norm = std::sqrt(std::pow(this->nrigid_normal[n + cmesh.nuc], 2) + //
                                      std::pow(this->nrigid_normal[n + cmesh.nvc], 2) + //
                                      std::pow(this->nrigid_normal[n + cmesh.nwc], 2));
        if (norm > params.mtol) {
            this->nrigid_normal[n + cmesh.nuc] /= norm;
            this->nrigid_normal[n + cmesh.nvc] /= norm;
            this->nrigid_normal[n + cmesh.nwc] /= norm;
        } else {
            this->nrigid_normal[n + cmesh.nuc] = 0.0e0;
            this->nrigid_normal[n + cmesh.nvc] = 0.0e0;
            this->nrigid_normal[n + cmesh.nwc] = 0.0e0;
        }
    }

    return;
}

void SolidMaterialPointBase::UpdateDefGrad(int pid, int nenode, double af_coeff, const std::vector<int> &ncm,
                                           const std::vector<double> &sf, const std::vector<std::array<double, 3>> &dsf,
                                           std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                                           std::vector<std::array<std::array<double, 3>, 3>> &def_grad) {

    std::array<std::array<double, 3>, 3> ddg{};
    ddg = this->ComputeDeltaDefGrad(ncm, nenode, af_coeff, dsf);

    delta_def_grad[pid] = ddg;

    // --- Deformation gradient update: F^{n+1} = dF * F^n ---
    std::array<std::array<double, 3>, 3> F{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) { F[i][j] += ddg[i][k] * this->def_grad[pid][k][j]; }
        }
    }
    def_grad[pid] = F;

    this->det_def_grad[pid] = DetMat3(F);

    return;
}

void SolidMaterialPointBase::UpdateConstitutiveModel(int pid, std::vector<std::array<double, 6>> &stress,               //
                                                     const std::vector<double> &det_def_grad_bar,                       //
                                                     const std::vector<std::array<std::array<double, 3>, 3>> &def_grad, //
                                                     const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad) {

    std::array<std::array<double, 3>, 3> F = def_grad[pid];
    std::array<std::array<double, 3>, 3> dF = delta_def_grad[pid];

    int ipp = this->materials_.iprop[this->matid[pid]];

    cm_.UpdateStress(ipp, stress[pid], F, dF, this->materials_.mat_prop[this->matid[pid]], //
                     this->det_def_grad[pid], det_def_grad_bar[pid]);

    return;
}
