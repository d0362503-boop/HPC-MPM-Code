#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "../../bc.h"
#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "../../solver/crsmat.h"
#include "../../solver/solver.h"
#include "stabilized_fem.h"

using namespace stabilizedfem;

std::vector<double> StabilizedFEM::ComputeAdvectionVel() {

    std::vector<double> adv_vel(nodec * 3);
    for (int n = 0; n < nodec * 3; n++) {
        // --- Generalized-alpha prediciton ---
        adv_vel[n] = (1.0e0 + this->alpha_f) * this->nvel_old[n] //
                     - this->alpha_f * this->nvel_older[n];
    }

    return adv_vel;
}

void StabilizedFEM::SolveNS() {
    this->BCSet(); // --- user's responsibility (need to create in each src)---

    std::vector<double> adv_vel = this->ComputeAdvectionVel();

    this->MakNSStabCoeff(adv_vel); // ---- Stabilized coefficient ----

    this->AssembleNSSystem(adv_vel);

    // ---- Initialize LHS x from [nvel, npres] ----
    VectorAssign(nodec * 4, this->NS_.x_lhs);
    std::copy(this->nvel.begin(), this->nvel.end(), this->NS_.x_lhs.begin());
    std::copy(this->npres.begin(), this->npres.end(), this->NS_.x_lhs.begin() + npc);

    int iter = this->SolveSystem(NS_);

    // ---- Extract [nvel, npres] from x_lhs ----
    std::copy(this->NS_.x_lhs.begin(), this->NS_.x_lhs.begin() + npc, this->nvel.begin());
    std::copy(this->NS_.x_lhs.begin() + npc, this->NS_.x_lhs.end(), this->npres.begin());

    if (myrank == 0) { std::cout << "NS_iter:" << std::setw(15) << iter << "\n"; }

    // --- Memory clear ---
    std::vector<double>().swap(this->NS_.adiag);
    std::vector<double>().swap(this->NS_.amat);
    std::vector<double>().swap(this->NS_.b_rhs);
    std::vector<double>().swap(this->NS_.x_lhs);

    return;
}

void StabilizedFEM::MakNSStabCoeff(const std::vector<double> &adv_vel) {
    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    this->tau1.resize(nelem);
    this->tau2.resize(nelem);
    for (int m = 0; m < nelem; m++) {
        double rnue = this->rmue[m] / this->rhoe[m];

        int ize = m / (nex * ney);
        int iye = (m - ize * (nex * ney)) / nex;
        int ixe = m - ize * (nex * ney) - iye * nex;

        std::array<double, 3> xye;
        xye[0] = xymin[0] + dxy[0] * (double(ixe) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(iye) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ize) + 0.5e0);

        MakSf(m, xye, idimc, xynodec, ncm, nenode, sf, dsf);

        double uu = 0.0e0, vv = 0.0e0, ww = 0.0e0;
        for (int ni = 0; ni < nenode; ni++) {
            int nid = ncm[ni];
            uu += adv_vel[nid + nuc] / double(nenode);
            vv += adv_vel[nid + nvc] / double(nenode);
            ww += adv_vel[nid + nwc] / double(nenode);
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

        double t1 = (2.0e0 * dti) * (2.0e0 * dti);
        double t2 = 4.0e0 * uvw / (he * he);
        double t3 = std::pow((4.0e0 * rnue / (he * he)), 2);
        this->tau1[m] = 1.0e0 / std::sqrt(t1 + t2 + t3);

        if (he < mtol) {
            this->tau2[m] = 0.0e0;
        } else {
            double reu = std::sqrt(uvw) * he / (rnue * 2.0e0);
            if (reu > 3.0e0) {
                this->tau2[m] = std::sqrt(uvw) * he / 2.0e0;
            } else {
                this->tau2[m] = std::sqrt(uvw) * he * reu / 6.0e0;
            }
        }
    }

    return;
}

void StabilizedFEM::AssembleNSSystem(const std::vector<double> &adv_vel) {
    double fx = bb[0] * facl;
    double fy = bb[1] * facl;
    double fz = bb[2] * facl;

    const double am = this->alpha_m / this->gamma_nb;
    const double af = this->alpha_f;
    const double af0 = 1.0e0 - af;
    const double am0 = 1.0e0 - am;

    double G_Weight = (dxy[0] * dxy[1] * dxy[2]) / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->NS_.nmata, this->NS_.amat);
    VectorAssign(nodec * 4, this->NS_.adiag);
    VectorAssign(nodec * 4, this->NS_.b_rhs);
    for (int m = 0; m < nelem; m++) {
        double t1 = this->tau1[m];
        double t2 = this->tau2[m];

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            double phi_k = 0.0e0;
            std::array<double, 3> adv_vel_k{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                // --- Interface funtion ---
                phi_k += sfi * this->nphi[nid];
                // --- Advection velocity ---
                adv_vel_k[0] += sfi * adv_vel[nid + nuc];
                adv_vel_k[1] += sfi * adv_vel[nid + nvc];
                adv_vel_k[2] += sfi * adv_vel[nid + nwc];
            }

            double rhof = phi_k * rhol + (1.0e0 - phi_k) * rhog;
            double rmuf = phi_k * rmul + (1.0e0 - phi_k) * rmug;

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                double ust_ni = (dsfi1 * adv_vel_k[0] + dsfi2 * adv_vel_k[1] + dsfi3 * adv_vel_k[2]);
                // --- Force (Galerkin+SUPG) ---
                double fg = (sfi + ust_ni * t1) * G_Weight * rhof;
                // --- Force (PSPG) ---
                double fp = (dsfi1 * fx + dsfi2 * fy + dsfi3 * fz) * G_Weight * t1;

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];
                    double ust_nj = (dsfj1 * adv_vel_k[0] + dsfj2 * adv_vel_k[1] + dsfj3 * adv_vel_k[2]);

                    int ida = this->NS_.FindIndex(nid, njd, ncol);

                    // --- Mass Matrix (Galerkin+SUPG) ---
                    double emd = G_Weight * (sfi + ust_ni * t1) * sfj * rhof;
                    // --- Mass Matrix (PSPG) ---
                    double egtu = G_Weight * dsfi1 * sfj * t1;
                    double egtv = G_Weight * dsfi2 * sfj * t1;
                    double egtw = G_Weight * dsfi3 * sfj * t1;
                    // --- Advection Matrix (Galerkin+SUPG) ---
                    double akg = G_Weight * (sfi + ust_ni * t1) * ust_nj * rhof;
                    // --- Advection Matrix (PSPG) ---
                    double akpu = G_Weight * dsfi1 * ust_nj * t1;
                    double akpv = G_Weight * dsfi2 * ust_nj * t1;
                    double akpw = G_Weight * dsfi3 * ust_nj * t1;
                    // --- Diffusion Matrix (Galerkin part) ---
                    double su = G_Weight * (2.0 * dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * rmuf;
                    double sv = G_Weight * (dsfi1 * dsfj1 + 2.0 * dsfi2 * dsfj2 + dsfi3 * dsfj3) * rmuf;
                    double sw = G_Weight * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + 2.0 * dsfi3 * dsfj3) * rmuf;
                    double suv = G_Weight * (dsfi2 * dsfj1) * rmuf;
                    double suw = G_Weight * (dsfi3 * dsfj1) * rmuf;
                    double svu = G_Weight * (dsfi1 * dsfj2) * rmuf;
                    double svw = G_Weight * (dsfi3 * dsfj2) * rmuf;
                    double swu = G_Weight * (dsfi1 * dsfj3) * rmuf;
                    double swv = G_Weight * (dsfi2 * dsfj3) * rmuf;
                    // --- Pressure Matrix (SUPG-Galerkin) ---
                    double esgu = G_Weight * (ust_ni * dsfj1 * t1 - dsfi1 * sfj);
                    double esgv = G_Weight * (ust_ni * dsfj2 * t1 - dsfi2 * sfj);
                    double esgw = G_Weight * (ust_ni * dsfj3 * t1 - dsfi3 * sfj);
                    // --- Pressure Matrix (PSPG) ---
                    double elt = G_Weight * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) / rhof * t1;
                    // --- Continuity (Galerkin part) ---
                    double Cowu = G_Weight * sfi * dsfj1;
                    double Cowv = G_Weight * sfi * dsfj2;
                    double Coww = G_Weight * sfi * dsfj3;
                    // --- Shock Capturing ---
                    double scu = G_Weight * dsfi1 * dsfj1 * rhof * t2;
                    double scv = G_Weight * dsfi2 * dsfj2 * rhof * t2;
                    double scw = G_Weight * dsfi3 * dsfj3 * rhof * t2;
                    double scuv = G_Weight * dsfi1 * dsfj2 * rhof * t2;
                    double scuw = G_Weight * dsfi1 * dsfj3 * rhof * t2;
                    double scvu = G_Weight * dsfi2 * dsfj1 * rhof * t2;
                    double scvw = G_Weight * dsfi2 * dsfj3 * rhof * t2;
                    double scwu = G_Weight * dsfi3 * dsfj1 * rhof * t2;
                    double scwv = G_Weight * dsfi3 * dsfj2 * rhof * t2;

                    this->NS_.amat[ida + this->NS_.block_id[0]] += dti * am * emd + af * (akg + su + scu);
                    this->NS_.amat[ida + this->NS_.block_id[1]] += af * (suv + scuv);
                    this->NS_.amat[ida + this->NS_.block_id[2]] += af * (suw + scuw);
                    this->NS_.amat[ida + this->NS_.block_id[3]] += esgu * af;
                    this->NS_.amat[ida + this->NS_.block_id[4]] += af * (svu + scvu);
                    this->NS_.amat[ida + this->NS_.block_id[5]] += dti * am * emd + af * (akg + sv + scv);
                    this->NS_.amat[ida + this->NS_.block_id[6]] += af * (svw + scvw);
                    this->NS_.amat[ida + this->NS_.block_id[7]] += esgv * af;
                    this->NS_.amat[ida + this->NS_.block_id[8]] += af * (swu + scwu);
                    this->NS_.amat[ida + this->NS_.block_id[9]] += af * (swv + scwv);
                    this->NS_.amat[ida + this->NS_.block_id[10]] += dti * am * emd + af * (akg + sw + scw);
                    this->NS_.amat[ida + this->NS_.block_id[11]] += esgw * af;
                    this->NS_.amat[ida + this->NS_.block_id[12]] += dti * am * egtu + af * (Cowu + akpu);
                    this->NS_.amat[ida + this->NS_.block_id[13]] += dti * am * egtv + af * (Cowv + akpv);
                    this->NS_.amat[ida + this->NS_.block_id[14]] += dti * am * egtw + af * (Coww + akpw);
                    this->NS_.amat[ida + this->NS_.block_id[15]] += elt * af;

                    const double u_old = this->nvel_old[njd + nuc];
                    const double v_old = this->nvel_old[njd + nvc];
                    const double w_old = this->nvel_old[njd + nwc];
                    const double p_old = this->npres_old[njd];

                    const double au_old = this->naccel[njd + nuc];
                    const double av_old = this->naccel[njd + nvc];
                    const double aw_old = this->naccel[njd + nwc];

                    // x-momentum
                    this->NS_.b_rhs[nid + nuc] += (dti * am * emd - af0 * (akg + su + scu)) * u_old //
                                                  - af0 * (suv + scuv) * v_old                      //
                                                  - af0 * (suw + scuw) * w_old                      //
                                                  - af0 * esgu * p_old                              //
                                                  - am0 * emd * au_old;

                    // y-momentum
                    this->NS_.b_rhs[nid + nvc] += (dti * am * emd - af0 * (akg + sv + scv)) * v_old //
                                                  - af0 * (svu + scvu) * u_old                      //
                                                  - af0 * (svw + scvw) * w_old                      //
                                                  - af0 * esgv * p_old                              //
                                                  - am0 * emd * av_old;

                    // z-momentum
                    this->NS_.b_rhs[nid + nwc] += (dti * am * emd - af0 * (akg + sw + scw)) * w_old //
                                                  - af0 * (swu + scwu) * u_old                      //
                                                  - af0 * (swv + scwv) * v_old                      //
                                                  - af0 * esgw * p_old                              //
                                                  - am0 * emd * aw_old;

                    // continuity / pressure
                    this->NS_.b_rhs[nid + npc] += (dti * am * egtu - af0 * (Cowu + akpu)) * u_old   //
                                                  + (dti * am * egtv - af0 * (Cowv + akpv)) * v_old //
                                                  + (dti * am * egtw - af0 * (Coww + akpw)) * w_old //
                                                  - af0 * elt * p_old                               //
                                                  - am0 * (egtu * au_old + egtv * av_old + egtw * aw_old);
                }
                this->NS_.b_rhs[nid + nuc] += fg * fx;
                this->NS_.b_rhs[nid + nvc] += fg * fy;
                this->NS_.b_rhs[nid + nwc] += fg * fz;
                this->NS_.b_rhs[nid + npc] += fp;
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->NS_.b_rhs, {nuc, nvc, nwc, npc});

    // ---- Diagional Scalling Prec (only for non-PETSc solver) ----
    if (!this->NS_.use_petsc) { this->NS_.BuildDiagonalPreconditioner(this->NS_.ndof); }

    return;
}

void StabilizedFEM::UpdateNodalVar() {
    std::copy(this->nvel_old.begin(), this->nvel_old.end(), this->nvel_older.begin());
    std::copy(this->nvel.begin(), this->nvel.end(), this->nvel_old.begin());
    std::copy(this->npres.begin(), this->npres.end(), this->npres_old.begin());

    return;
}
