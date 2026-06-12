#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../bc.h"
#include "../../cal_mat.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "implicit_mpm_solid.h"

using namespace implicitmpm;

void SolidMaterialPoint::Particle2Node() {
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
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
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
    NodeVarComm(this->nmome, {nuc, nvc, nwc});
    NodeVarComm(this->nforce, {nuc, nvc, nwc});

    VectorAssign(nodec * 3, this->nvel);
    VectorAssign(nodec * 3, this->naccel);
    for (int n = 0; n < nodec; n++) {
        if (this->nmass[n] > mtol) {
            this->nvel[n + nuc] = this->nmome[n + nuc] / this->nmass[n];
            this->nvel[n + nvc] = this->nmome[n + nvc] / this->nmass[n];
            this->nvel[n + nwc] = this->nmome[n + nwc] / this->nmass[n];
            this->naccel[n + nuc] = this->nforce[n + nuc] / this->nmass[n];
            this->naccel[n + nvc] = this->nforce[n + nvc] / this->nmass[n];
            this->naccel[n + nwc] = this->nforce[n + nwc] / this->nmass[n];
        }
    }

    this->ubc.BCSetVal(nuc, this->nvel);
    this->vbc.BCSetVal(nvc, this->nvel);
    this->wbc.BCSetVal(nwc, this->nvel);
    this->ubc.BCSetZero(nuc, this->naccel);
    this->vbc.BCSetZero(nvc, this->naccel);
    this->wbc.BCSetZero(nwc, this->naccel);

    this->rigid_bc.BCSetVal(nuc, this->nvel);
    this->rigid_bc.BCSetVal(nvc, this->nvel);
    this->rigid_bc.BCSetVal(nwc, this->nvel);
    this->rigid_bc.BCSetZero(nuc, this->naccel);
    this->rigid_bc.BCSetZero(nvc, this->naccel);
    this->rigid_bc.BCSetZero(nwc, this->naccel);

    return;
}

void SolidMaterialPoint::UpdateDefGrad(int pid, int nenode, double af_coeff, const std::vector<int> &ncm,
                                       const std::vector<double> &sf, const std::vector<std::array<double, 3>> &dsf,
                                       std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                                       std::vector<std::array<std::array<double, 3>, 3>> &def_grad) {
    std::array<std::array<double, 3>, 3> ddg{};
    ddg = this->ComputeDeltaDefGrad(ncm, nenode, af_coeff, dsf);

    delta_def_grad[pid] = ddg;

    // --- Deformation gradient update ---
    std::array<std::array<double, 3>, 3> F{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) { F[i][j] += ddg[i][k] * this->def_grad[pid][k][j]; }
        }
    }
    def_grad[pid] = F;

    this->det_def_grad[pid] = DetMat3(F);
    this->vol[pid] = this->det_def_grad[pid] * this->vol0[pid];

    return;
}

void SolidMaterialPoint::Node2Particle() {
    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);
    // ---- Newmark beta velocity & acceleration ----
    this->PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k);

    this->CommitNodalKinematics(nvel_k, naccel_k);

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
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->UpdateDefGrad(pid, nenode, 1.0e0, ncm, sf, dsf, delta_def_grad, this->def_grad);

            // std::array<std::array<double, 3>, 3> F = this->def_grad[pid];
            // this->det_def_grad[pid] = DetMat3(F);
            // this->vol[pid] = this->det_def_grad[pid] * this->vol0[pid];

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

    std::vector<std::array<double, 3>> delta_corr;
    delta_corr = this->DeltaCorrectionParticleShifting();

    this->CommitParticleKinematics(accel_old, displ, delta_corr);

    return;
}
