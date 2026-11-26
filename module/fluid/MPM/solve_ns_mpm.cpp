#include "module/bc.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/implicit/newton_solver.h"
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

void StabilizedMPM::ConfigurePreconditioner(CrsMat &mat, PC pc) {

    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, mat.HasReducedLayout() ? 1 : mat.ndof);
    mat.field_splits = {{"velocity", {0, 1, 2}}, {"pressure", {3}}};
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);

    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_pc_type", "hypre");
    PetscOptionsSetValue(nullptr, "-fieldsplit_velocity_pc_hypre_type", "boomeramg");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_type", "hypre");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_type", "boomeramg");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_coarsen_type", "hmis");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_interp_type", "ext+i");
    PetscOptionsSetValue(nullptr, "-fieldsplit_pressure_pc_hypre_boomeramg_max_iter", "2");

    return;
}

void StabilizedMPM::UpdateNRIncrement() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    for (int n = 0; n < cmesh.nodec * 3; n++) { this->ndispl[n] += this->NS_.x_lhs[n]; }
    for (int n = 0; n < cmesh.nodec; n++) { this->npres[n] += this->NS_.x_lhs[n + cmesh.npc]; }

    return;
}

void StabilizedMPM::SolveNS(const TimeState &ts) {

    SolveNewton(*this, ts, this->max_NR_it);

    return;
}

void StabilizedMPM::BeginNewtonStep() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->nvel_k_.assign(cmesh.nodec * 3, 0.0e0);
    this->naccel_k_.assign(cmesh.nodec * 3, 0.0e0);

    this->MakeNSStabCoeff(this->nvel);

    VectorAssign(cmesh.nodec * 3, this->ndispl);
    VectorAssign(cmesh.nodec, this->npres);
    this->r0r_ = 0.0e0;

    return;
}

void StabilizedMPM::AssembleNewtonSystem(const TimeState &ts) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    this->BCNRSet(ts);

    this->implicit_.ComputeNodeVelAccelFromDispl(*this, this->nvel_k_,
                                                 this->naccel_k_); // ---- Newmark beta velocity & acceleration ----

    VectorAssign(this->NS_.nmata, this->NS_.amat);
    VectorAssign(cmesh.nodec * 4, this->NS_.b_rhs);
    VectorAssign(cmesh.nodec * 4, this->NS_.x_lhs); // ---- Initialize LHS x value ----
    this->AssembleSystem(this->NS_, this->nvel_k_, this->naccel_k_, ts);

    return;
}

bool StabilizedMPM::CheckNewtonConvergence(int NR_it, int ksp_its) {

    return this->NS_.CheckNRConvergence(NR_it, max_NR_it, ksp_its, this->r0r_);
}

int StabilizedMPM::SolveNewtonCorrection(int NR_it) { return this->NS_.SolveSystem(NR_it); }

void StabilizedMPM::EndNewtonStep() {

    // --- Memory clear ---
    std::vector<double>().swap(this->NS_.adiag);
    std::vector<double>().swap(this->NS_.amat);
    std::vector<double>().swap(this->NS_.b_rhs);
    std::vector<double>().swap(this->NS_.x_lhs);

    return;
}

void StabilizedMPM::MakeNSStabCoeff(const std::vector<double> &nvel_k) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    double rnu = this->rmu / this->rho;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(this->num, this->tau_pspg);
    VectorAssign(this->num, this->tau_lsic);
    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(m, xyp, cmesh.idimc, cmesh.xynodec, ncm, nenode, sf, dsf);

            double uu = 0.0e0, vv = 0.0e0, ww = 0.0e0;
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                uu += sfi * nvel_k[nid + cmesh.nuc];
                vv += sfi * nvel_k[nid + cmesh.nvc];
                ww += sfi * nvel_k[nid + cmesh.nwc];
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
            if (udsf < params.mtol) {
                he = mesh.dxy[0] * std::cbrt((6.0e0 / M_PI));
            } else {
                he = 2.0e0 * std::sqrt(uvw) / udsf;
            }

            if constexpr (this->stab_coeff == StabCoeff::VMS) { // ---- VMS ----
                double t1 = this->rho * params.dti;
                double t2 = 2.0e0 * this->rho * std::sqrt(uvw) / he;
                double t3 = 4.0e0 * this->rmu / (he * he);
                this->tau_pspg[pid] = 1.0e0 / (t1 + t2 + t3);
                this->tau_lsic[pid] = (he * he) / (4.0e0 * this->tau_pspg[pid]);
            } else if constexpr (this->stab_coeff == StabCoeff::PSPG) { // ---- PSPG ----
                double t1 = (2.0e0 * params.dti) * (2.0e0 * params.dti);
                double t2 = 4.0e0 * uvw / (he * he);
                double t3 = std::pow((4.0e0 * rnu / (he * he)), 2);
                double tau = 1.0e0 / std::sqrt(t1 + t2 + t3);
                this->tau_pspg[pid] = tau / this->rho;
                this->tau_lsic[pid] = tau * uvw * this->rho;
            }
            pid = this->idp2p[pid];
        }
    }

    return;
}

void StabilizedMPM::AssembleSystem(CrsMat &mat, const std::vector<double> &nvel_k, //
                                   const std::vector<double> &naccel_k, const TimeState &ts) {
    const Mesh &mesh = this->grid_.mesh;
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    // double A = 5.0e0 / 180.0e0 * M_PI;
    // double theta = A * std::sin(5.47e0 * real_time);
    // double fx = bb[2] * std::sin(theta);
    // double fy = bb[1] * ts.facl;
    // double fz = bb[2] * std::cos(theta);
    double fx = params.bb[0] * ts.facl;
    double fy = params.bb[1] * ts.facl;
    double fz = params.bb[2] * ts.facl;

    const std::vector<int> offsets = this->RHSOffsets();

    const double af = this->implicit_.integrator_.alpha_f;
    const double af0 = 1.0e0 - af;
    const double am = this->implicit_.integrator_.alpha_m;
    const double am0 = 1.0e0 - am;
    const double cv = this->implicit_.integrator_.dv_du;
    const double ca = this->implicit_.integrator_.da_du;

    std::vector<double> nvel_af(cmesh.nodec * 3), naccel_am(cmesh.nodec * 3);
    for (int n = 0; n < cmesh.nodec * 3; n++) {
        naccel_am[n] = am0 * this->naccel[n] + am * naccel_k[n];
        nvel_af[n] = af0 * this->nvel[n] + af * nvel_k[n];
    }

    std::vector<double> npres_af(cmesh.nodec);
    for (int n = 0; n < cmesh.nodec; n++) { //
        npres_af[n] = af0 * this->npres_old[n] + af * this->npres[n];
    }

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    for (int m = 0; m < mesh.nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(m, xyp, cmesh.idimc, cmesh.xynodec, ncm, nenode, sf, dsf);

            double t1 = this->tau_pspg[pid];
            double t2 = this->tau_lsic[pid];
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
                accel_k[0] += sfi * naccel_am[nid + cmesh.nuc];
                accel_k[1] += sfi * naccel_am[nid + cmesh.nvc];
                accel_k[2] += sfi * naccel_am[nid + cmesh.nwc];
                // --- Gradient of pressure ---
                grad_pres_k[0] += dsfi1 * npres_af[nid];
                grad_pres_k[1] += dsfi2 * npres_af[nid];
                grad_pres_k[2] += dsfi3 * npres_af[nid];
                // --- Gradient of velocity ---
                grad_vel_k[0][0] += dsfi1 * nvel_af[nid + cmesh.nuc];
                grad_vel_k[0][1] += dsfi2 * nvel_af[nid + cmesh.nuc];
                grad_vel_k[0][2] += dsfi3 * nvel_af[nid + cmesh.nuc];
                grad_vel_k[1][0] += dsfi1 * nvel_af[nid + cmesh.nvc];
                grad_vel_k[1][1] += dsfi2 * nvel_af[nid + cmesh.nvc];
                grad_vel_k[1][2] += dsfi3 * nvel_af[nid + cmesh.nvc];
                grad_vel_k[2][0] += dsfi1 * nvel_af[nid + cmesh.nwc];
                grad_vel_k[2][1] += dsfi2 * nvel_af[nid + cmesh.nwc];
                grad_vel_k[2][2] += dsfi3 * nvel_af[nid + cmesh.nwc];
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
                double emd_lu = this->implicit_.LumpedMass(*this, pid, sfi);

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];

                    int ida = mat.FindIndex(nid, njd, ncol);

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
                    // --- Pressure (Galerkin part) ---
                    double esgx = volp * dsfi1 * sfj;
                    double esgy = volp * dsfi2 * sfj;
                    double esgz = volp * dsfi3 * sfj;
                    // --- Pressure (PSPG): ∫ (τ/ρ) ∇N_i·∇N_j ---
                    double elt = volp * (dsfi1 * dsfj1 + dsfi2 * dsfj2 + dsfi3 * dsfj3) * t1;
                    // --- Continuity (Galerkin part) ---
                    double Cow1 = volp * sfi * dsfj1;
                    double Cow2 = volp * sfi * dsfj2;
                    double Cow3 = volp * sfi * dsfj3;
                    // --- LSIC (grad-div): ∫ ρ τ_LSIC ∂N_i/∂x_k ∂N_j/∂x_l ---
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
                        mat.amat[ida + mat.block_id[this->blocks[0]]] += am * emd_lu;
                        mat.amat[ida + mat.block_id[this->blocks[5]]] += am * emd_lu;
                        mat.amat[ida + mat.block_id[this->blocks[10]]] += am * emd_lu;
                    }
                    mat.amat[ida + mat.block_id[this->blocks[0]]] += cv * af * (su + scu);
                    mat.amat[ida + mat.block_id[this->blocks[1]]] += cv * af * (suv + scuv);
                    mat.amat[ida + mat.block_id[this->blocks[2]]] += cv * af * (suw + scuw);
                    mat.amat[ida + mat.block_id[this->blocks[3]]] -= af * esgx;
                    mat.amat[ida + mat.block_id[this->blocks[4]]] += cv * af * (svu + scvu);
                    mat.amat[ida + mat.block_id[this->blocks[5]]] += cv * af * (sv + scv);
                    mat.amat[ida + mat.block_id[this->blocks[6]]] += cv * af * (svw + scvw);
                    mat.amat[ida + mat.block_id[this->blocks[7]]] -= af * esgy;
                    mat.amat[ida + mat.block_id[this->blocks[8]]] += cv * af * (swu + scwu);
                    mat.amat[ida + mat.block_id[this->blocks[9]]] += cv * af * (swv + scwv);
                    mat.amat[ida + mat.block_id[this->blocks[10]]] += cv * af * (sw + scw);
                    mat.amat[ida + mat.block_id[this->blocks[11]]] -= af * esgz;
                    mat.amat[ida + mat.block_id[this->blocks[12]]] += cv * af * Cow1 + ca * am * egtx;
                    mat.amat[ida + mat.block_id[this->blocks[13]]] += cv * af * Cow2 + ca * am * egty;
                    mat.amat[ida + mat.block_id[this->blocks[14]]] += cv * af * Cow3 + ca * am * egtz;
                    mat.amat[ida + mat.block_id[this->blocks[15]]] += af * elt;
                }

                std::array<double, 4> RHS_G{}, RHS_S{};
                RHS_G[0] = volp * (dsfi1 * stress_k[0][0] + dsfi2 * stress_k[0][1] + dsfi3 * stress_k[0][2]) //
                           - sfi * massp * fx;
                RHS_G[1] = volp * (dsfi1 * stress_k[1][0] + dsfi2 * stress_k[1][1] + dsfi3 * stress_k[1][2]) //
                           - sfi * massp * fy;
                RHS_G[2] = volp * (dsfi1 * stress_k[2][0] + dsfi2 * stress_k[2][1] + dsfi3 * stress_k[2][2]) //
                           - sfi * massp * fz;
                RHS_G[3] = volp * sfi * TraceMat3(grad_vel_k);

                RHS_S[0] = volp * t2 * dsfi1 * TraceMat3(grad_vel_k);
                RHS_S[1] = volp * t2 * dsfi2 * TraceMat3(grad_vel_k);
                RHS_S[2] = volp * t2 * dsfi3 * TraceMat3(grad_vel_k);
                RHS_S[3] = t1 * massp * (dsfi1 * (accel_k[0] - fx) + dsfi2 * (accel_k[1] - fy) + dsfi3 * (accel_k[2] - fz)) //
                           + t1 * volp * (dsfi1 * grad_pres_k[0] + dsfi2 * grad_pres_k[1] + dsfi3 * grad_pres_k[2]);

                mat.b_rhs[nid + offsets[0]] -= (RHS_G[0] + RHS_S[0]);
                mat.b_rhs[nid + offsets[1]] -= (RHS_G[1] + RHS_S[1]);
                mat.b_rhs[nid + offsets[2]] -= (RHS_G[2] + RHS_S[2]);
                mat.b_rhs[nid + offsets[3]] -= (RHS_G[3] + RHS_S[3]);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(mat.b_rhs, offsets);

    this->implicit_.AddInertialForce(*this, mat, naccel_am, offsets);

    return;
}
