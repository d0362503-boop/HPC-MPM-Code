#include "../../bc.h"
#include "../../cal_mat.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "../ConstitutiveModel.h"
#include "../fbar_projection.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void Particle2NodeS() {

    double fx = bb[0] * facl;
    double fy = bb[1] * facl;
    double fz = bb[2] * facl;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, sp.nmass);
    VectorAssign(nodec, sp.nvof);
    VectorAssign(nodec * 3, sp.nmome);
    VectorAssign(nodec * 3, sp.nforce);
    for (int m = 0; m < nelem; m++) {
        std::vector<double> nfint(3), nfext(3);

        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = sp.coord[pid];
            std::array<double, 6> sts = sp.stress[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                sp.StandardVarP2G(pid, nid, sfi, sp.mass, sp.nmass);
                sp.StandardVarP2G(pid, nid, sfi, sp.vol, sp.nvof);
                if (sp.solswitch == MapScheme::FLIP || sp.solswitch == MapScheme::PIC) {
                    sp.StandardVarP2G(pid, nid, sfi, sp.mass, sp.vel, sp.nmome);
                } else if (sp.solswitch == MapScheme::TPIC) {
                    sp.tpic.TPICVarP2G(pid, nid, sfi, sp.tpic.vel_grad, sp.mass, sp.vel, sp.coord, sp.nmome);
                } else if (sp.solswitch == MapScheme::APIC) {
                    sp.apic.APICVarP2G(pid, nid, sfi, sp.apic.vel_Bmat, sp.mass, sp.vel, sp.coord, sp.nmome);
                }
                nfint[0] = -sp.vol[pid] * (sts[0] * dsfi1 + sts[5] * dsfi2 + sts[4] * dsfi3);
                nfint[1] = -sp.vol[pid] * (sts[5] * dsfi1 + sts[1] * dsfi2 + sts[3] * dsfi3);
                nfint[2] = -sp.vol[pid] * (sts[4] * dsfi1 + sts[3] * dsfi2 + sts[2] * dsfi3);
                nfext[0] = sfi * (sp.mass[pid] * fx + sp.trac_force[pid][0]);
                nfext[1] = sfi * (sp.mass[pid] * fy + sp.trac_force[pid][1]);
                nfext[2] = sfi * (sp.mass[pid] * fz + sp.trac_force[pid][2]);
                sp.nforce[nid + nuc] += nfint[0] + nfext[0];
                sp.nforce[nid + nvc] += nfint[1] + nfext[1];
                sp.nforce[nid + nwc] += nfint[2] + nfext[2];
            }
            pid = sp.idp2p[pid];
        }
    }

    NodeVarComm(sp.nmass, 0);
    NodeVarComm(sp.nvof, 0);
    NodeVarComm(sp.nmome, {nuc, nvc, nwc});
    NodeVarComm(sp.nforce, {nuc, nvc, nwc});

    VectorAssign(nodec * 3, sp.nvel);
    for (int n = 0; n < nodec; n++) {
        if (sp.nmass[n] > mtol) {
            sp.nvel[n + nuc] = sp.nmome[n + nuc] / sp.nmass[n];
            sp.nvel[n + nvc] = sp.nmome[n + nvc] / sp.nmass[n];
            sp.nvel[n + nwc] = sp.nmome[n + nwc] / sp.nmass[n];
        } else {
            sp.nvel[n + nuc] = 0.0e0;
            sp.nvel[n + nvc] = 0.0e0;
            sp.nvel[n + nwc] = 0.0e0;
        }
        if (sp.nvof[n] < mtol) { sp.nvof[n] = 0.0e0; }
    }

    BCSetVal(nuc, sp.nvel, usbc);
    BCSetVal(nvc, sp.nvel, vsbc);
    BCSetVal(nwc, sp.nvel, wsbc);
    BCSetVal(nuc, sp.nvel, rigid_bc);
    BCSetVal(nvc, sp.nvel, rigid_bc);
    BCSetVal(nwc, sp.nvel, rigid_bc);

    return;
}

void MUSL() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec * 3, sp.nmome);
    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = {sp.coord[pid][0], sp.coord[pid][1], sp.coord[pid][2]};
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                if (sp.solswitch == MapScheme::FLIP || sp.solswitch == MapScheme::PIC) {
                    sp.StandardVarP2G(pid, nid, sfi, sp.mass, sp.vel, sp.nmome);
                } else if (sp.solswitch == MapScheme::TPIC) {
                    sp.tpic.TPICVarP2G(pid, nid, sfi, sp.tpic.vel_grad, sp.mass, sp.vel, sp.coord, sp.nmome);
                } else if (sp.solswitch == MapScheme::APIC) {
                    sp.apic.APICVarP2G(pid, nid, sfi, sp.apic.vel_Bmat, sp.mass, sp.vel, sp.coord, sp.nmome);
                }
            }
            pid = sp.idp2p[pid];
        }
    }

    NodeVarComm(sp.nmome, {nuc, nvc, nwc});

    VectorAssign(nodec * 3, sp.nvel);
    for (int n = 0; n < nodec; n++) {
        if (sp.nmass[n] > mtol) {
            sp.nvel[n + nuc] = sp.nmome[n + nuc] / sp.nmass[n];
            sp.nvel[n + nvc] = sp.nmome[n + nvc] / sp.nmass[n];
            sp.nvel[n + nwc] = sp.nmome[n + nwc] / sp.nmass[n];
        }
    }

    BCSetVal(nuc, sp.nvel, usbc);
    BCSetVal(nvc, sp.nvel, vsbc);
    BCSetVal(nwc, sp.nvel, wsbc);
    BCSetVal(nuc, sp.nvel, rigid_bc);
    BCSetVal(nvc, sp.nvel, rigid_bc);
    BCSetVal(nwc, sp.nvel, rigid_bc);

    return;
}

void update_point_coord_and_stress(const std::vector<double> &det_def_grad_bar,
                                   const std::vector<std::vector<std::vector<double>>> &def_grad,
                                   const std::vector<std::vector<std::vector<double>>> &delta_def_grad) {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    std::vector<std::array<double, 3>> delta_corr;
    VectorAssign(sp.num, delta_corr);
    delta_corr = sp.DeltaCorrectionParticleShifting();

    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = {sp.coord[pid][0], sp.coord[pid][1], sp.coord[pid][2]};

            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                sp.flip.FLIPVarG2P(pid, nid, sfi, sp.coord, sp.nvel);
            }
            for (int i = 0; i < 3; i++) { sp.coord[pid][i] += delta_corr[pid][i]; }

            sp.vol[pid] = sp.det_def_grad[pid] * sp.vol0[pid];

            ConstitutiveModel(pid, sp.stress, det_def_grad_bar, def_grad, delta_def_grad);

            pid = sp.idp2p[pid];
        }
    }

    return;
}

void Node2ParticleS() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    if (sp.solswitch != MapScheme::FLIP) {
        for (int n = 0; n < nodec; n++) {
            sp.nvel[n + nuc] += dt * sp.naccel[n + nuc];
            sp.nvel[n + nvc] += dt * sp.naccel[n + nvc];
            sp.nvel[n + nwc] += dt * sp.naccel[n + nwc];
        }
    }

    if (sp.solswitch != MapScheme::FLIP) { VectorAssign(sp.num, sp.vel); }
    if (sp.solswitch == MapScheme::TPIC) {
        // ----- TPIC -----
        VectorAssign(sp.num, sp.tpic.vel_grad);
    } else if (sp.solswitch == MapScheme::APIC) {
        // ----- APIC -----
        VectorAssign(sp.num, sp.apic.vel_Bmat);
        VectorAssign(sp.num, sp.apic.inv_Dmat);
    }
    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = sp.coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            // if (wp.solswitch == MapScheme::APIC) {
            std::array<std::array<double, 3>, 3> ADp{};
            //}
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                if (sp.solswitch != MapScheme::FLIP) {
                    sp.pic.PICVarG2P(pid, nid, sfi, sp.vel, sp.nvel);
                } else if (sp.solswitch == MapScheme::FLIP) {
                    sp.flip.FLIPVarG2P(pid, nid, sfi, sp.vel, sp.naccel);
                }
                if (sp.solswitch == MapScheme::TPIC) {
                    sp.tpic.TPICVarGradG2P(pid, nid, ni, dsf, sp.tpic.vel_grad, sp.nvel);
                } else if (sp.solswitch == MapScheme::APIC) {
                    sp.apic.APICVarBmatG2P(pid, nid, sfi, sp.apic.vel_Bmat, sp.nvel, sp.coord);
                    sp.apic.APICDmat(pid, nid, sfi, sp.coord, ADp);
                }
            }
            if (sp.solswitch == MapScheme::APIC) { sp.apic.APICInvDmat(pid, ADp); }
            pid = sp.idp2p[pid];
        }
    }

    MUSL();

    std::vector<double> det_delta_def_grad(sp.num, 0.0e0);
    std::vector<std::vector<std::vector<double>>> delta_def_grad;
    VectorAssign(sp.num, delta_def_grad);
    for (int m = 0; m < nelem; m++) {
        int pid = sp.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = sp.coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                // FLIPVarG2P(pid, nid, sfi, sp.coord, sp.nvel);
                delta_def_grad[0][0][pid] += dsfi1 * sp.nvel[nid + nuc] * dt;
                delta_def_grad[0][1][pid] += dsfi2 * sp.nvel[nid + nuc] * dt;
                delta_def_grad[0][2][pid] += dsfi3 * sp.nvel[nid + nuc] * dt;
                delta_def_grad[1][0][pid] += dsfi1 * sp.nvel[nid + nvc] * dt;
                delta_def_grad[1][1][pid] += dsfi2 * sp.nvel[nid + nvc] * dt;
                delta_def_grad[1][2][pid] += dsfi3 * sp.nvel[nid + nvc] * dt;
                delta_def_grad[2][0][pid] += dsfi1 * sp.nvel[nid + nwc] * dt;
                delta_def_grad[2][1][pid] += dsfi2 * sp.nvel[nid + nwc] * dt;
                delta_def_grad[2][2][pid] += dsfi3 * sp.nvel[nid + nwc] * dt;
            }
            for (int i = 0; i < 3; i++) {
                // sp.coord[i][pid] += delta_corr[i][pid];
                delta_def_grad[i][i][pid]++;
            }

            // --- Deformation gradient update ---
            std::array<std::array<double, 3>, 3> F{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) { F[i][j] += delta_def_grad[i][k][pid] * sp.def_grad[pid][k][j]; }
                }
            }

            std::array<std::array<double, 3>, 3> dF{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sp.def_grad[pid][i][j] = F[i][j];
                    dF[i][j] = delta_def_grad[i][j][pid];
                }
            }
            sp.det_def_grad[pid] = DetMat3(F);
            det_delta_def_grad[pid] = DetMat3(dF);
            // sp.vol[pid] = sp.det_def_grad[pid] * sp.vol0[pid];

            pid = sp.idp2p[pid];
        }
    }

    if (Fbar_flag == true) {
        // --- F_bar projection ---
        cal_def_grad_bar(delta_def_grad, det_delta_def_grad);
        // --- Update point coordinate and stress with Fbar ---
        update_point_coord_and_stress(sp.det_def_grad_bar, // delta_corr,
                                      sp.def_grad_bar, delta_def_grad);
    } else {
        // --- Update point coordinate and stress without Fbar ---
        update_point_coord_and_stress(sp.det_def_grad, // delta_corr,
                                      sp.def_grad, delta_def_grad);
    }

    return;
}
