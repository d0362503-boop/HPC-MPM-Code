#include "../../bc.h"
#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "../../solver/crsmat.h"
#include "../../solver/solver.h"
#include "stabilized_mpm.h"
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

    constexpr int iter_max = 10000;

    this->MakNSStabCoeff(); // ---- Stabilized coefficient ----

    VectorAssign(nodec * 3, this->ndispl);
    VectorAssign(nodec, this->npres);
    VectorAssign(nodec * 4, this->NS_.x_lhs); // ---- Initialize LHS x value ----
    double r0r = 0.0e0;
    for (int NR_it = 0; NR_it <= iter_max; NR_it++) {

        this->BCNRSet();

        this->PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k); // ---- Newmark beta velocity & acceleration ----

        this->AssembleNSSystem(nvel_k, naccel_k);

        int iter = this->SolveSystem(this->NS_, NR_it);

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

void StabilizedMPM::MakNSStabCoeff() {

    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];
    double rnu = this->rmu / this->rho;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->num, this->tau1);
    VectorAssign(this->num, this->tau2);
    for (int m = 0; m < nelem; m++) {
        int ize = m / (nex * ney);
        int iye = (m - ize * (nex * ney)) / nex;
        int ixe = m - ize * (nex * ney) - iye * nex;

        std::array<double, 3> xye;
        xye[0] = xymin[0] + dxy[0] * (double(ixe) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(iye) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ize) + 0.5e0);

        int pid = this->idepf[m];
        while (pid != -1) {
            double uu = 0.0e0, vv = 0.0e0, ww = 0.0e0;
            MakSf(m, xye, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                uu += sfi * this->nvel[nid + nuc];
                vv += sfi * this->nvel[nid + nvc];
                ww += sfi * this->nvel[nid + nwc];
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
    // double theta = A * sin(5.47e0*real_time);
    // double fx = bb[2] * sin(theta);
    // double fy = bb[1] * facl;
    // double fz = bb[2] * cos(theta);
    double fx = bb[0] * facl;
    double fy = bb[1] * facl;
    double fz = bb[2] * facl;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->NS_.nmata, this->NS_.amat);
    VectorAssign(nodec * 4, this->NS_.adiag);
    VectorAssign(nodec * 4, this->NS_.b_rhs);

    for (int m = 0; m < nelem; m++) {

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->ImplicitDsfCorr(ncm, nenode, dsf);

            double t1 = this->tau1[pid];
            double t2 = this->tau2[pid];

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
                pres_k += sfi * this->npres[nid];
                // --- Acceleration ---
                accel_k[0] += sfi * naccel_k[nid + nuc];
                accel_k[1] += sfi * naccel_k[nid + nvc];
                accel_k[2] += sfi * naccel_k[nid + nwc];
                // --- Gradient of pressure ---
                grad_pres_k[0] += dsfi1 * this->npres[nid];
                grad_pres_k[1] += dsfi2 * this->npres[nid];
                grad_pres_k[2] += dsfi3 * this->npres[nid];
                // --- Gradient of velocity ---
                grad_vel_k[0][0] += dsfi1 * nvel_k[nid + nuc];
                grad_vel_k[0][1] += dsfi2 * nvel_k[nid + nuc];
                grad_vel_k[0][2] += dsfi3 * nvel_k[nid + nuc];
                grad_vel_k[1][0] += dsfi1 * nvel_k[nid + nvc];
                grad_vel_k[1][1] += dsfi2 * nvel_k[nid + nvc];
                grad_vel_k[1][2] += dsfi3 * nvel_k[nid + nvc];
                grad_vel_k[2][0] += dsfi1 * nvel_k[nid + nwc];
                grad_vel_k[2][1] += dsfi2 * nvel_k[nid + nwc];
                grad_vel_k[2][2] += dsfi3 * nvel_k[nid + nwc];
            }
            // --- Deviatoric stress ---
            std::array<std::array<double, 3>, 3> dev_stress_k{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) { dev_stress_k[i][j] = this->rmu * (grad_vel_k[i][j] + grad_vel_k[j][i]); }
            }

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                // --- Mass Matrix (Galerkin part) ---
                double emd_lu = this->ComputeNRLumpedMat(pid, sfi);

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];

                    int ida = this->NS_.FindIndex(nid, njd, ncol);

                    // --- Mass Matrix (Stabilized part) ---
                    double egtx = this->mass[pid] * dsfi1 * sfj * t1;
                    double egty = this->mass[pid] * dsfi2 * sfj * t1;
                    double egtz = this->mass[pid] * dsfi3 * sfj * t1;
                    // --- Diffusion (Galerkin part) ---
                    double su = this->vol[pid] * (2.0 * dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * this->rmu;
                    double sv = this->vol[pid] * (dsfi1 * dsfj1 + 2.0 * dsfi2 * dsfj2 + dsfi3 * dsfj3) * this->rmu;
                    double sw = this->vol[pid] * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + 2.0 * dsfi3 * dsfj3) * this->rmu;
                    double suv = this->vol[pid] * (dsfi2 * dsfj1) * this->rmu;
                    double suw = this->vol[pid] * (dsfi3 * dsfj1) * this->rmu;
                    double svu = this->vol[pid] * (dsfi1 * dsfj2) * this->rmu;
                    double svw = this->vol[pid] * (dsfi3 * dsfj2) * this->rmu;
                    double swu = this->vol[pid] * (dsfi1 * dsfj3) * this->rmu;
                    double swv = this->vol[pid] * (dsfi2 * dsfj3) * this->rmu;
                    // --- Preesure (Galerkin part) ---
                    double esgx = this->vol[pid] * dsfi1 * sfj;
                    double esgy = this->vol[pid] * dsfi2 * sfj;
                    double esgz = this->vol[pid] * dsfi3 * sfj;
                    // --- Preesure (Stabilized part) ---
                    double elt = this->vol[pid] * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * t1;
                    // --- Continuity (Galerkin part) ---
                    double Cow1 = this->vol[pid] * sfi * dsfj1;
                    double Cow2 = this->vol[pid] * sfi * dsfj2;
                    double Cow3 = this->vol[pid] * sfi * dsfj3;
                    // --- Shock Capturing (Stabilized part) ---
                    double scu = this->vol[pid] * dsfi1 * dsfj1 * t2;
                    double scv = this->vol[pid] * dsfi2 * dsfj2 * t2;
                    double scw = this->vol[pid] * dsfi3 * dsfj3 * t2;
                    double scuv = this->vol[pid] * dsfi1 * dsfj2 * t2;
                    double scuw = this->vol[pid] * dsfi1 * dsfj3 * t2;
                    double scvu = this->vol[pid] * dsfi2 * dsfj1 * t2;
                    double scvw = this->vol[pid] * dsfi2 * dsfj3 * t2;
                    double scwu = this->vol[pid] * dsfi3 * dsfj1 * t2;
                    double scwv = this->vol[pid] * dsfi3 * dsfj2 * t2;

                    if (nid == njd) {
                        this->NS_.amat[ida + this->NS_.block_id[0]] += emd_lu;
                        this->NS_.amat[ida + this->NS_.block_id[5]] += emd_lu;
                        this->NS_.amat[ida + this->NS_.block_id[10]] += emd_lu;
                    }
                    this->NS_.amat[ida + this->NS_.block_id[0]] += this->nb_para[0] * (su + scu);
                    this->NS_.amat[ida + this->NS_.block_id[1]] += this->nb_para[0] * (suv + scuv);
                    this->NS_.amat[ida + this->NS_.block_id[2]] += this->nb_para[0] * (suw + scuw);
                    this->NS_.amat[ida + this->NS_.block_id[3]] -= esgx;
                    this->NS_.amat[ida + this->NS_.block_id[4]] += this->nb_para[0] * (svu + scvu);
                    this->NS_.amat[ida + this->NS_.block_id[5]] += this->nb_para[0] * (sv + scv);
                    this->NS_.amat[ida + this->NS_.block_id[6]] += this->nb_para[0] * (svw + scvw);
                    this->NS_.amat[ida + this->NS_.block_id[7]] -= esgy;
                    this->NS_.amat[ida + this->NS_.block_id[8]] += this->nb_para[0] * (swu + scwu);
                    this->NS_.amat[ida + this->NS_.block_id[9]] += this->nb_para[0] * (swv + scwv);
                    this->NS_.amat[ida + this->NS_.block_id[10]] += this->nb_para[0] * (sw + scw);
                    this->NS_.amat[ida + this->NS_.block_id[11]] -= esgz;
                    this->NS_.amat[ida + this->NS_.block_id[12]] += this->nb_para[0] * Cow1 + this->nb_para[3] * egtx;
                    this->NS_.amat[ida + this->NS_.block_id[13]] += this->nb_para[0] * Cow2 + this->nb_para[3] * egty;
                    this->NS_.amat[ida + this->NS_.block_id[14]] += this->nb_para[0] * Cow3 + this->nb_para[3] * egtz;
                    this->NS_.amat[ida + this->NS_.block_id[15]] += elt;
                }

                std::array<double, 4> RHS_G{}, RHS_S{};
                RHS_G[0] = this->vol[pid] *
                               (dsfi1 * dev_stress_k[0][0] + dsfi2 * dev_stress_k[0][1] + dsfi3 * dev_stress_k[0][2]) //
                           - sfi * this->mass[pid] * fx - this->vol[pid] * dsfi1 * pres_k;
                RHS_G[1] = this->vol[pid] *
                               (dsfi1 * dev_stress_k[1][0] + dsfi2 * dev_stress_k[1][1] + dsfi3 * dev_stress_k[1][2]) //
                           - sfi * this->mass[pid] * fy - this->vol[pid] * dsfi2 * pres_k;
                RHS_G[2] = this->vol[pid] *
                               (dsfi1 * dev_stress_k[2][0] + dsfi2 * dev_stress_k[2][1] + dsfi3 * dev_stress_k[2][2]) //
                           - sfi * this->mass[pid] * fz - this->vol[pid] * dsfi3 * pres_k;
                RHS_G[3] = this->vol[pid] * sfi * (grad_vel_k[0][0] + grad_vel_k[1][1] + grad_vel_k[2][2]);

                RHS_S[0] = this->vol[pid] * t2 * dsfi1 * (grad_vel_k[0][0] + grad_vel_k[1][1] + grad_vel_k[2][2]);
                RHS_S[1] = this->vol[pid] * t2 * dsfi2 * (grad_vel_k[0][0] + grad_vel_k[1][1] + grad_vel_k[2][2]);
                RHS_S[2] = this->vol[pid] * t2 * dsfi3 * (grad_vel_k[0][0] + grad_vel_k[1][1] + grad_vel_k[2][2]);
                RHS_S[3] =
                    t1 * this->mass[pid] *
                        (dsfi1 * (accel_k[0] - fx) + dsfi2 * (accel_k[1] - fy) + dsfi3 * (accel_k[2] - fz)) //
                    + t1 * this->vol[pid] * (dsfi1 * grad_pres_k[0] + dsfi2 * grad_pres_k[1] + dsfi3 * grad_pres_k[2]);

                this->NS_.b_rhs[nid + nuc] -= (RHS_G[0] + RHS_S[0]);
                this->NS_.b_rhs[nid + nvc] -= (RHS_G[1] + RHS_S[1]);
                this->NS_.b_rhs[nid + nwc] -= (RHS_G[2] + RHS_S[2]);
                this->NS_.b_rhs[nid + npc] -= (RHS_G[3] + RHS_S[3]);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->NS_.b_rhs, {nuc, nvc, nwc, npc});

    this->AddInertialForceToRHS(this->NS_, naccel_k);

    // ---- Diagional Scalling Prec ----
    if (!this->NS_.use_petsc) { this->NS_.BuildDiagonalPreconditioner(this->NS_.ndof); }

    return;
}