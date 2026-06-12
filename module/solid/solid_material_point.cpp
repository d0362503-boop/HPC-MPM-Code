#include "solid_material_point.h"

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "../bc.h"
#include "../dataset.h"
#include "../material_point.h"
#include "../mesh.h"
#include "../mpi_data.h"
#include "../shape_function.h"

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

void SolidMaterialPointBase::Moveparticle() {
    this->DetermineParticleRank();

    if (this->par_comm_.nmps + this->par_comm_.nrps + this->par_comm_.nrmp == 0) return;

    this->num += this->par_comm_.nrps - this->par_comm_.nmps - this->par_comm_.nrmp;

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
        // ---- Point determinant of deformation gradient bar ----
        this->par_comm_.PointVarComm(this->num, this->det_def_grad_bar);
        // --- Point deformation gradient bar ---
        this->par_comm_.PointVarComm(this->num, this->def_grad_bar);
    }

    if (this->solswitch == MapScheme::TPIC) {
        // ---- TPIC acceleration gradient ----
        this->par_comm_.PointVarComm(this->num, this->tpic.accel_grad);
        // ---- TPIC velocity gredient ----
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
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nvof);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            if (iprop[this->matid[pid]] == -1) {
                std::array<double, 3> xyp = {this->coord[pid][0], this->coord[pid][1], this->coord[pid][2]};
                MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
                for (int ni = 0; ni < nenode; ni++) {
                    int nid = ncm[ni];
                    double sfi = sf[ni];
                    this->nvof[nid] += sfi * this->vol[pid];
                }
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nvof, 0);

    VectorAssign(nodec, this->rigid_bc.nbc);
    VectorAssign(nodec, this->rigid_bc.fbc);
    for (int n = 0; n < nodec; n++) {
        if (this->nvof[n] > 0.0e0) {
            this->rigid_bc.nbc[this->rigid_bc.ibc] = n;
            this->rigid_bc.fbc[this->rigid_bc.ibc] = 0.0e0;
            this->rigid_bc.ibc++;
        }
    }

    return;
}

void SolidMaterialPointBase::UpdateConstitutiveModel(
    int pid,                                                           //
    std::vector<std::array<double, 6>> &stress,                        //
    const std::vector<double> &det_def_grad_bar,                       //
    const std::vector<std::array<std::array<double, 3>, 3>> &def_grad, //
    const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad) {

    std::array<std::array<double, 3>, 3> F = def_grad[pid];
    std::array<std::array<double, 3>, 3> dF = delta_def_grad[pid];

    int ipp = iprop[this->matid[pid]];

    cm_.UpdateStress(ipp, stress[pid], F, dF, mat_prop[this->matid[pid]], //
                     this->det_def_grad[pid], det_def_grad_bar[pid]);

    return;
}
