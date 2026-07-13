#include "../../bc.h"
#include "../../cal_mat.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "../constitutive_model.h"
#include "explicit_mpm_solid.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
using namespace explicitmpm;

void ExplicitSolidMPM::Particle2Node() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nmass);
    VectorAssign(nodec, this->nvof);
    VectorAssign(nodec * 3, this->nmome);
    VectorAssign(nodec * 3, this->nforce);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            std::array<double, 6> sts = this->stress[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->StandardVarP2G(pid, nid, sfi, this->mass, this->nmass);
                this->StandardVarP2G(pid, nid, sfi, this->vol, this->nvof);
                this->VelP2G(pid, nid, sfi);
                std::array<double, 3> nfint = this->ComputeInternalForce(ni, pid, dsf, sts);
                std::array<double, 3> nfext = this->ComputeExternalForce(pid, sfi);
                this->nforce[nid + nuc] += nfint[0] + nfext[0];
                this->nforce[nid + nvc] += nfint[1] + nfext[1];
                this->nforce[nid + nwc] += nfint[2] + nfext[2];
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nmass, 0);
    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->nmome, {nuc, nvc, nwc});
    NodeVarComm(this->nforce, {nuc, nvc, nwc});

    this->CutOffSmallNodalVar(this->nvel, this->nmome, this->nmass, {nuc, nvc, nwc});
    this->ApplyVelocityBC(this->nvel);

    return;
}

void ExplicitSolidMPM::SolveSolid() {

    this->CutOffSmallNodalVar(this->naccel, this->nforce, this->nmass, {nuc, nvc, nwc});
    this->ApplyAccelerationBC(this->naccel);

    return;
}

void ExplicitSolidMPM::DoMUSL() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec * 3, this->nmome);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->VelP2G(pid, nid, sfi);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nmome, {nuc, nvc, nwc});

    this->CutOffSmallNodalVar(this->nvel, this->nmome, this->nmass, {nuc, nvc, nwc});
    this->ApplyVelocityBC(this->nvel);

    return;
}

void ExplicitSolidMPM::Node2Particle() {

    this->G2PVelocity();

    this->DoMUSL();

    this->UpdateDeformationGradient();

    this->UpdateParticlePositionAndStress();

    return;
}

void ExplicitSolidMPM::G2PVelocity() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    // --- Nodal velocity update for PIC-family schemes ---
    if (this->solswitch != MapScheme::FLIP) {
        for (int n = 0; n < nodec * 3; n++) { this->nvel[n] += dt * this->naccel[n]; }
    }

    // --- Reset particle kinematics ---
    this->PICFamilyVelReset();
    if (this->solswitch == MapScheme::APIC) { VectorAssign(this->num, this->apic.inv_Dmat); }

    // --- G2P: velocity update (PIC/TPIC/APIC/FLIP) ---
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            std::array<std::array<double, 3>, 3> ADp{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                if (this->solswitch == MapScheme::FLIP) {
                    this->flip.VarG2P(pid, nid, sfi, this->vel, this->naccel);
                } else {
                    this->PICFamilyVelG2P(pid, ni, nid, sfi, dsf);
                }
                if (this->solswitch == MapScheme::APIC) { this->apic.DmatG2P(pid, nid, sfi, this->coord, ADp); }
            }
            if (this->solswitch == MapScheme::APIC) { this->apic.InvDmatG2P(pid, ADp); }
            pid = this->idp2p[pid];
        }
    }

    return;
}

void ExplicitSolidMPM::UpdateDeformationGradient() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    // --- Nodal displacement increment for deformation-gradient update ---
    VectorAssign(nodec * 3, this->ndispl);
    for (int n = 0; n < nodec * 3; n++) { this->ndispl[n] = this->nvel[n] * dt; }

    // --- Update deformation gradient ---
    std::vector<double> det_delta_def_grad(this->num, 0.0e0);
    VectorAssign(this->num, this->delta_def_grad);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->UpdateDefGrad(pid, nenode, 1.0e0, ncm, sf, dsf, this->delta_def_grad, this->def_grad);
            det_delta_def_grad[pid] = DetMat3(this->delta_def_grad[pid]);

            pid = this->idp2p[pid];
        }
    }

    if (this->Fbar_flag) { this->ComputeDefGradBar(this->delta_def_grad, det_delta_def_grad); }

    return;
}

void ExplicitSolidMPM::UpdateParticlePositionAndStress() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    // --- Particle shifting correction ---
    std::vector<std::array<double, 3>> delta_corr;
    delta_corr = this->DeltaCorrectionParticleShifting();

    // --- Particle coordinate / volume / stress update ---
    const bool has_shift = !delta_corr.empty();
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            // 1. Position update
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->flip.VarG2P(pid, nid, sfi, this->coord, this->nvel);
            }
            if (has_shift) {
                for (int i = 0; i < 3; i++) { this->coord[pid][i] += delta_corr[pid][i]; }
            }

            // 2. Volume and stress update
            if (this->Fbar_flag) {
                this->UpdateVolume(pid, this->det_def_grad_bar[pid]);
                this->UpdateConstitutiveModel(pid, this->stress, this->det_def_grad_bar, this->def_grad_bar,
                                              this->delta_def_grad_bar);
            } else {
                this->UpdateVolume(pid, this->det_def_grad[pid]);
                this->UpdateConstitutiveModel(pid, this->stress, this->det_def_grad, this->def_grad,
                                              this->delta_def_grad);
            }

            pid = this->idp2p[pid];
        }
    }

    return;
}