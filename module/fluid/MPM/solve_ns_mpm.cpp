#include "module/bc.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include "module/solver/crsmat.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

using namespace stabilizedmpm;

void StabilizedMPM::UpdateNRIncrement() {
    for (int n = 0; n < nodec * 3; n++) { this->ndispl[n] += this->NS_.x_lhs[n]; }
    for (int n = 0; n < nodec; n++) { this->npres[n] += this->NS_.x_lhs[n + npc]; }

    return;
}

void StabilizedMPM::SolveNS() {

    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);

    const int iter_max = 1000;

    VectorAssign(nodec * 3, this->ndispl);
    VectorAssign(nodec, this->npres);
    VectorAssign(nodec * 4, this->NS_.x_lhs); // ---- Initialize LHS x value ----
    double r0r = 0.0e0;
    for (int NR_it = 0; NR_it <= iter_max; NR_it++) {

        this->BCNRSet();

        this->ComputeNodeVelAccelFromDispl(nvel_k, naccel_k); // ---- Newmark beta velocity & acceleration ----

        this->AssembleNSSystem(nvel_k, naccel_k);

        int iter = this->NS_.SolveSystem(NR_it);

        this->UpdateNRIncrement();

        if (this->NS_.CheckNRConvergence(NR_it, iter_max, iter, r0r)) { break; }
    }

    // --- Memory clear ---
    std::vector<double>().swap(this->NS_.adiag);
    std::vector<double>().swap(this->NS_.amat);
    std::vector<double>().swap(this->NS_.b_rhs);
    std::vector<double>().swap(this->NS_.x_lhs);

    return;
}

void StabilizedMPM::MakNSStabCoeff(const std::vector<double> &nvel_k) {

    double rnu = this->rmu / this->rho;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->num, this->tau1);
    VectorAssign(this->num, this->tau2);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            double uu = 0.0e0, vv = 0.0e0, ww = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                uu += sfi * nvel_k[nid + nuc];
                vv += sfi * nvel_k[nid + nvc];
                ww += sfi * nvel_k[nid + nwc];
            }
            double uvw = uu * uu + vv * vv + ww * ww;

            double udsf = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                udsf += std::abs(uu * dsfi1 + vv * dsfi2 + ww * dsfi3);
            }

            double he;
            if (udsf < mtol) {
                he = dxy[0] * std::cbrt((6.0e0 / M_PI));
            } else {
                he = 2.0e0 * std::sqrt(uvw) / udsf;
            }

            if constexpr (this->stab_coeff == StabCoeff::VMS) { // ---- VMS ----
                double t1 = this->rho * dti;
                double t2 = 2.0e0 * this->rho * std::sqrt(uvw) / he;
                double t3 = 4.0e0 * this->rmu / (he * he);
                this->tau1[pid] = 1.0e0 / (t1 + t2 + t3);
                this->tau2[pid] = (he * he) / (4.0e0 * this->tau1[pid]);
            } else if constexpr (this->stab_coeff == StabCoeff::PSPG) { // ---- PSPG ----
                double t1 = (2.0e0 * dti) * (2.0e0 * dti);
                double t2 = 4.0e0 * uvw / (he * he);
                double t3 = std::pow((4.0e0 * rnu / (he * he)), 2);
                double tau = 1.0e0 / std::sqrt(t1 + t2 + t3);
                this->tau1[pid] = tau / this->rho;
                this->tau2[pid] = tau * uvw * this->rho;
            }
            pid = this->idp2p[pid];
        }
    }

    return;
}

void StabilizedMPM::AssembleNSSystem(const std::vector<double> &nvel_k, //
                                     const std::vector<double> &naccel_k) {

    // double A = 5.0e0 / 180.0e0 * M_PI;
    // double theta = A * std::sin(5.47e0 * real_time);
    // double fx = bb[2] * std::sin(theta);
    // double fy = bb[1] * facl;
    // double fz = bb[2] * std::cos(theta);
    double fx = bb[0] * facl;
    double fy = bb[1] * facl;
    double fz = bb[2] * facl;

    const double af = this->alpha_f;
    const double af0 = 1.0e0 - this->alpha_f;
    const double am = this->alpha_m;
    const double am0 = 1.0e0 - this->alpha_m;

    std::vector<double> nvel_af(nodec * 3), naccel_am(nodec * 3);
    for (int n = 0; n < nodec * 3; n++) {
        naccel_am[n] = am0 * this->naccel[n] + am * naccel_k[n];
        nvel_af[n] = af0 * this->nvel[n] + af * nvel_k[n];
    }
    std::vector<double> npres_af(nodec);
    for (int n = 0; n < nodec; n++) { //
        npres_af[n] = af0 * this->npres_old[n] + af * this->npres[n];
    }

    this->MakNSStabCoeff(nvel_af); // ---- Stabilized coefficient ----

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->NS_.nmata, this->NS_.amat);
    VectorAssign(nodec * 4, this->NS_.b_rhs);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->ImplicitDsfCorr(ncm, nenode, dsf);

            double t1 = this->tau1[pid];
            double t2 = this->tau2[pid];
            double volp = this->vol[pid];
            double massp = this->mass[pid];

            double pres_k = 0.0e0;
            std::array<double, 3> accel_k{};
            std::array<double, 3> grad_pres_k{};
            std::array<std::array<double, 3>, 3> grad_vel_k{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                // --- Pressure ----
                pres_k += sfi * npres_af[nid];
                // --- Acceleration ---
                accel_k[0] += sfi * naccel_am[nid + nuc];
                accel_k[1] += sfi * naccel_am[nid + nvc];
                accel_k[2] += sfi * naccel_am[nid + nwc];
                // --- Gradient of pressure ---
                grad_pres_k[0] += dsfi1 * npres_af[nid];
                grad_pres_k[1] += dsfi2 * npres_af[nid];
                grad_pres_k[2] += dsfi3 * npres_af[nid];
                // --- Gradient of velocity ---
                grad_vel_k[0][0] += dsfi1 * nvel_af[nid + nuc];
                grad_vel_k[0][1] += dsfi2 * nvel_af[nid + nuc];
                grad_vel_k[0][2] += dsfi3 * nvel_af[nid + nuc];
                grad_vel_k[1][0] += dsfi1 * nvel_af[nid + nvc];
                grad_vel_k[1][1] += dsfi2 * nvel_af[nid + nvc];
                grad_vel_k[1][2] += dsfi3 * nvel_af[nid + nvc];
                grad_vel_k[2][0] += dsfi1 * nvel_af[nid + nwc];
                grad_vel_k[2][1] += dsfi2 * nvel_af[nid + nwc];
                grad_vel_k[2][2] += dsfi3 * nvel_af[nid + nwc];
            }
            // --- Deviatoric stress ---
            std::array<std::array<double, 3>, 3> stress_k{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) { stress_k[i][j] = this->rmu * (grad_vel_k[i][j] + grad_vel_k[j][i]); }
                stress_k[i][i] -= pres_k;
            }

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                // --- Mass Matrix (Galerkin part) ---
                double emd_lu = this->ComputeNRLumpedMassMat(pid, sfi);

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];

                    int ida = this->NS_.FindIndex(nid, njd, ncol);

                    // --- Mass Matrix (Stabilized part) ---
                    double egtx = massp * dsfi1 * sfj * t1;
                    double egty = massp * dsfi2 * sfj * t1;
                    double egtz = massp * dsfi3 * sfj * t1;
                    // --- Diffusion (Galerkin part) ---
                    double su = volp * (2.0e0 * dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * this->rmu;
                    double sv = volp * (dsfi1 * dsfj1 + 2.0e0 * dsfi2 * dsfj2 + dsfi3 * dsfj3) * this->rmu;
                    double sw = volp * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + 2.0e0 * dsfi3 * dsfj3) * this->rmu;
                    double suv = volp * (dsfi2 * dsfj1) * this->rmu;
                    double suw = volp * (dsfi3 * dsfj1) * this->rmu;
                    double svu = volp * (dsfi1 * dsfj2) * this->rmu;
                    double svw = volp * (dsfi3 * dsfj2) * this->rmu;
                    double swu = volp * (dsfi1 * dsfj3) * this->rmu;
                    double swv = volp * (dsfi2 * dsfj3) * this->rmu;
                    // --- Preesure (Galerkin part) ---
                    double esgx = volp * dsfi1 * sfj;
                    double esgy = volp * dsfi2 * sfj;
                    double esgz = volp * dsfi3 * sfj;
                    // --- Preesure (Stabilized part) ---
                    double elt = volp * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * t1;
                    // --- Continuity (Galerkin part) ---
                    double Cow1 = volp * sfi * dsfj1;
                    double Cow2 = volp * sfi * dsfj2;
                    double Cow3 = volp * sfi * dsfj3;
                    // --- Shock Capturing (Stabilized part) ---
                    double scu = volp * dsfi1 * dsfj1 * t2;
                    double scv = volp * dsfi2 * dsfj2 * t2;
                    double scw = volp * dsfi3 * dsfj3 * t2;
                    double scuv = volp * dsfi1 * dsfj2 * t2;
                    double scuw = volp * dsfi1 * dsfj3 * t2;
                    double scvu = volp * dsfi2 * dsfj1 * t2;
                    double scvw = volp * dsfi2 * dsfj3 * t2;
                    double scwu = volp * dsfi3 * dsfj1 * t2;
                    double scwv = volp * dsfi3 * dsfj2 * t2;

                    if (nid == njd) {
                        this->NS_.amat[ida + this->NS_.block_id[0]] += am * emd_lu;
                        this->NS_.amat[ida + this->NS_.block_id[5]] += am * emd_lu;
                        this->NS_.amat[ida + this->NS_.block_id[10]] += am * emd_lu;
                    }
                    this->NS_.amat[ida + this->NS_.block_id[0]] += this->nb_para[0] * af * (su + scu);
                    this->NS_.amat[ida + this->NS_.block_id[1]] += this->nb_para[0] * af * (suv + scuv);
                    this->NS_.amat[ida + this->NS_.block_id[2]] += this->nb_para[0] * af * (suw + scuw);
                    this->NS_.amat[ida + this->NS_.block_id[3]] -= af * esgx;
                    this->NS_.amat[ida + this->NS_.block_id[4]] += this->nb_para[0] * af * (svu + scvu);
                    this->NS_.amat[ida + this->NS_.block_id[5]] += this->nb_para[0] * af * (sv + scv);
                    this->NS_.amat[ida + this->NS_.block_id[6]] += this->nb_para[0] * af * (svw + scvw);
                    this->NS_.amat[ida + this->NS_.block_id[7]] -= af * esgy;
                    this->NS_.amat[ida + this->NS_.block_id[8]] += this->nb_para[0] * af * (swu + scwu);
                    this->NS_.amat[ida + this->NS_.block_id[9]] += this->nb_para[0] * af * (swv + scwv);
                    this->NS_.amat[ida + this->NS_.block_id[10]] += this->nb_para[0] * af * (sw + scw);
                    this->NS_.amat[ida + this->NS_.block_id[11]] -= af * esgz;
                    this->NS_.amat[ida + this->NS_.block_id[12]] += this->nb_para[0] * af * Cow1 //
                                                                    + this->nb_para[3] * am * egtx;
                    this->NS_.amat[ida + this->NS_.block_id[13]] += this->nb_para[0] * af * Cow2 //
                                                                    + this->nb_para[3] * am * egty;
                    this->NS_.amat[ida + this->NS_.block_id[14]] += this->nb_para[0] * af * Cow3 //
                                                                    + this->nb_para[3] * am * egtz;
                    this->NS_.amat[ida + this->NS_.block_id[15]] += af * elt;
                }

                std::array<double, 4> RHS_G{}, RHS_S{};
                RHS_G[0] = volp * (dsfi1 * stress_k[0][0] + dsfi2 * stress_k[0][1] + dsfi3 * stress_k[0][2]) //
                           - sfi * massp * fx; // - volp * dsfi1 * pres_k;
                RHS_G[1] = volp * (dsfi1 * stress_k[1][0] + dsfi2 * stress_k[1][1] + dsfi3 * stress_k[1][2]) //
                           - sfi * massp * fy; // - volp * dsfi2 * pres_k;
                RHS_G[2] = volp * (dsfi1 * stress_k[2][0] + dsfi2 * stress_k[2][1] + dsfi3 * stress_k[2][2]) //
                           - sfi * massp * fz; // - volp * dsfi3 * pres_k;
                RHS_G[3] = volp * sfi * TraceMat3(grad_vel_k);

                RHS_S[0] = volp * t2 * dsfi1 * TraceMat3(grad_vel_k);
                RHS_S[1] = volp * t2 * dsfi2 * TraceMat3(grad_vel_k);
                RHS_S[2] = volp * t2 * dsfi3 * TraceMat3(grad_vel_k);
                RHS_S[3] =
                    t1 * massp * (dsfi1 * (accel_k[0] - fx) + dsfi2 * (accel_k[1] - fy) + dsfi3 * (accel_k[2] - fz)) //
                    + t1 * volp * (dsfi1 * grad_pres_k[0] + dsfi2 * grad_pres_k[1] + dsfi3 * grad_pres_k[2]);

                this->NS_.b_rhs[nid + nuc] -= (RHS_G[0] + RHS_S[0]);
                this->NS_.b_rhs[nid + nvc] -= (RHS_G[1] + RHS_S[1]);
                this->NS_.b_rhs[nid + nwc] -= (RHS_G[2] + RHS_S[2]);
                this->NS_.b_rhs[nid + npc] -= (RHS_G[3] + RHS_S[3]);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->NS_.b_rhs, {nuc, nvc, nwc, npc});

    this->AddInertialForceToRHS(this->NS_, naccel_am);

    return;
}
