#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "module/bc.h"
#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/map_and_interpolate.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

using namespace implicitmpm;

void ImplicitSolidMPM::Particle2Node(const TimeState &ts) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(cmesh.nodec, this->nmass);
    VectorAssign(cmesh.nodec, this->nvof);
    VectorAssign(cmesh.nodec * 3, this->nmome);
    VectorAssign(cmesh.nodec * 3, this->nforce);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->StandardVarP2G(pid, nid, sfi, this->mass, this->nmass);
                this->StandardVarP2G(pid, nid, sfi, this->vol, this->nvof);
                this->VelP2G(pid, nid, sfi);
                this->AccelP2G(pid, nid, sfi);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nmass, 0);
    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->nmome, {cmesh.nuc, cmesh.nvc, cmesh.nwc});
    NodeVarComm(this->nforce, {cmesh.nuc, cmesh.nvc, cmesh.nwc});

    VectorAssign(cmesh.nodec, this->nphi);
    for (int n = 0; n < cmesh.nodec; n++) { this->nphi[n] = std::clamp(this->nvof[n] / cmesh.nvol[n], 0.0e0, 1.0e0); }

    this->CutOffSmallNodalVar(this->nvel, this->nmome, this->nmass, {cmesh.nuc, cmesh.nvc, cmesh.nwc});
    this->ApplyVelocityBC(this->nvel, ts);

    this->CutOffSmallNodalVar(this->naccel, this->nforce, this->nmass, {cmesh.nuc, cmesh.nvc, cmesh.nwc});
    this->ApplyAccelerationBC(this->naccel);

    return;
}

void ImplicitSolidMPM::Node2Particle(const TimeState &ts) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;

    // ---- Newmark beta velocity & acceleration ----
    std::vector<double> nvel_k(cmesh.nodec * 3), naccel_k(cmesh.nodec * 3);
    this->implicit_.ComputeNodeVelAccelFromDispl(*this, nvel_k, naccel_k);

    this->implicit_.CommitNodalKinematics(*this, nvel_k, naccel_k);

    std::vector<std::array<double, 3>> accel_old;
    VectorAssign(this->num, accel_old);
    if (this->solswitch == MapScheme::FLIP) {
        for (int n = 0; n < this->num; n++) { accel_old[n] = this->accel[n]; }
    }

    std::vector<std::array<double, 3>> displ;
    VectorAssign(this->num, displ);
    this->PICFamilyAccelReset();
    this->PICFamilyVelReset();
    if (this->solswitch == MapScheme::APIC) { VectorAssign(this->num, this->apic.inv_Dmat); }

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    std::vector<std::array<std::array<double, 3>, 3>> delta_def_grad;

    VectorAssign(this->num, delta_def_grad);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(this->grid_, m, xyp, ncm, nenode, sf, dsf);

            this->UpdateDefGrad(pid, nenode, 1.0e0, ncm, sf, dsf, delta_def_grad, this->def_grad);

            this->UpdateVolume(pid, this->det_def_grad[pid]);

            this->UpdateConstitutiveModel(pid, this->stress, this->det_def_grad, this->def_grad, delta_def_grad);

            if (iprop[this->matid[pid]] != -1) {
                std::array<std::array<double, 3>, 3> ADp{};
                for (int ni = 0; ni < nenode; ni++) {
                    int nid = ncm[ni];
                    double sfi = sf[ni];
                    this->pic.VarG2P(pid, nid, sfi, displ, this->ndispl); // ---- Displacement ----
                    this->PICFamilyAccelG2P(pid, ni, nid, sfi, dsf);
                    this->PICFamilyVelG2P(pid, ni, nid, sfi, dsf);
                    if (this->solswitch == MapScheme::APIC) { this->apic.DmatG2P(pid, nid, sfi, this->coord, ADp); }
                }
                if (this->solswitch == MapScheme::APIC) { this->apic.InvDmatG2P(pid, ADp); }
            }
            pid = this->idp2p[pid];
        }
    }

    std::vector<std::array<double, 3>> disp_corr = this->DeltaCorrectionPST();

    this->implicit_.CommitParticleKinematics(*this, accel_old, displ, disp_corr);

    return;
}
