#include <mpi.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "module/bc.h"
#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/implicit/newton_solver.h"
#include "module/map_and_interpolate.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solid/solid_material_point.h"
#include "module/solver/crsmat.h"

using namespace implicitmpm;

void ImplicitSolidMPM::ConfigurePreconditioner(CrsMat &mat, PC pc) {

    KSPGMRESSetRestart(mat.ksp, 200);
    PCSetType(pc, PCASM);
    PCASMSetOverlap(pc, 1);
    PCSetOptionsPrefix(pc, "solid_field_");
    PetscOptionsSetValue(nullptr, "-solid_field_sub_ksp_type", "preonly");
    PetscOptionsSetValue(nullptr, "-solid_field_sub_pc_type", "ilu");
    PetscOptionsSetValue(nullptr, "-solid_field_sub_pc_factor_levels", "1");
    PetscOptionsSetValue(nullptr, "-solid_field_sub_pc_factor_shift_type", "nonzero");
    PetscOptionsSetValue(nullptr, "-solid_field_sub_pc_factor_shift_amount", "1e-3");
    PCSetFromOptions(pc);

    return;
}

void ImplicitSolidMPM::UpdateNRIncrement() {
    for (int n = 0; n < nodec * 3; n++) { this->ndispl[n] += this->SM_.x_lhs[n]; }

    return;
}

std::array<std::array<double, 3>, 3> ImplicitSolidMPM::ComputeTangentModulus(int ni, int nj,
                                                                             const std::vector<std::array<double, 3>> &dsf,
                                                                             const std::array<double, 6> &sts_af) {
    double stsmat[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) { stsmat[i][j] = sts_af[idm[i][j]]; }
    }

    double ce[3][3][3][3] = {};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    ce[i][j][k][l] = this->cm_.stif_mat[idm[i][j]][idm[k][l]];
                    if (NR_flag) ce[i][j][k][l] += stsmat[j][l] * (i == k);
                }
            }
        }
    }

    std::array<std::array<double, 3>, 3> K_mat{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) { K_mat[i][j] += dsf[ni][k] * ce[i][k][j][l] * dsf[nj][l]; }
            }
        }
    }

    return K_mat;
}

void ImplicitSolidMPM::SolveSolid() {

    // --- NR_flag true for nonlinear elasticity and false for linear elasticity ---
    this->max_NR_it = (NR_flag) ? 100 : 0;

    SolveNewton(*this, this->max_NR_it);

    return;
}

void ImplicitSolidMPM::BeginNewtonStep() {

    this->nvel_k_.assign(nodec * 3, 0.0e0);
    this->naccel_k_.assign(nodec * 3, 0.0e0);

    this->stress_k_ = this->InitializeNRStress(); // --- Reset stress for NR ---

    VectorAssign(nodec * 3, this->ndispl); // ---- Reset displacement ----
    this->r0r_ = 0.0e0;

    return;
}

void ImplicitSolidMPM::AssembleNewtonSystem() {

    this->BCNRSet();

    this->implicit_.ComputeNodeVelAccelFromDispl(*this, this->nvel_k_, this->naccel_k_); // ---- Newmark beta velocity & acceleration ----

    VectorAssign(this->SM_.nmata, this->SM_.amat);
    VectorAssign(nodec * 3, this->SM_.b_rhs);
    VectorAssign(nodec * 3, this->SM_.x_lhs); // ---- Initialize LHS x value ----
    this->AssembleSystem(this->SM_, this->naccel_k_, this->stress_k_);

    return;
}

bool ImplicitSolidMPM::CheckNewtonConvergence(int NR_it, int ksp_its) {

    if (!NR_flag) { // --- Linear not need Newton Raphson ---
        if (myrank == 0) { std::cout << "Solid converge:" << std::setw(10) << NR_it << std::setw(10) << ksp_its << "\n"; }
        return false;
    }
    // --- Nonlinear need Newton Raphson ---
    return this->SM_.CheckNRConvergence(NR_it, this->max_NR_it, ksp_its, this->r0r_);
}

int ImplicitSolidMPM::SolveNewtonCorrection(int NR_it) { return this->SM_.SolveSystem(NR_it); }

void ImplicitSolidMPM::EndNewtonStep() {

    // --- Memory clear ---
    std::vector<double>().swap(this->SM_.adiag);
    std::vector<double>().swap(this->SM_.amat);
    std::vector<double>().swap(this->SM_.b_rhs);
    std::vector<double>().swap(this->SM_.x_lhs);

    return;
}

void ImplicitSolidMPM::AssembleSystem(CrsMat &mat, const std::vector<double> &naccel_k, //
                                      std::vector<std::array<double, 6>> &stress_k) {

    const std::vector<int> offsets = this->RHSOffsets();

    const double af = this->implicit_.integrator_.alpha_f;
    const double af0 = 1.0e0 - af;
    const double am = this->implicit_.integrator_.alpha_m;
    const double am0 = 1.0e0 - am;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    std::vector<std::array<std::array<double, 3>, 3>> delta_def_grad, def_grad_NR;

    VectorAssign(this->num, delta_def_grad);
    VectorAssign(this->num, def_grad_NR);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakeSF(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->UpdateDefGrad(pid, nenode, af, ncm, sf, dsf, delta_def_grad, def_grad_NR);

            this->UpdateVolume(pid, this->det_def_grad[pid]);

            this->UpdateConstitutiveModel(pid, stress_k, this->det_def_grad, def_grad_NR, delta_def_grad);

            this->ImplicitDsfCorr(ncm, nenode, dsf);

            // --- Transfer stress to (n+α_f) time interlevel ---
            std::array<double, 6> sts_af = stress_k[pid];
            double volp = this->vol[pid];

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];

                // --- Mass Matrix (lumpd) ---
                double emd_lu = this->ComputeNRLumpedMassMat(pid, sfi);

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];

                    int ida = mat.FindIndex(nid, njd, ncol);

                    std::array<std::array<double, 3>, 3> K_mat;
                    K_mat = this->ComputeTangentModulus(ni, nj, dsf, sts_af);

                    if (nid == njd) {
                        mat.amat[ida + mat.block_id[this->blocks[0]]] += am * emd_lu;
                        mat.amat[ida + mat.block_id[this->blocks[4]]] += am * emd_lu;
                        mat.amat[ida + mat.block_id[this->blocks[8]]] += am * emd_lu;
                    }
                    mat.amat[ida + mat.block_id[this->blocks[0]]] += af * K_mat[0][0] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[1]]] += af * K_mat[0][1] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[2]]] += af * K_mat[0][2] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[3]]] += af * K_mat[1][0] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[4]]] += af * K_mat[1][1] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[5]]] += af * K_mat[1][2] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[6]]] += af * K_mat[2][0] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[7]]] += af * K_mat[2][1] * volp;
                    mat.amat[ida + mat.block_id[this->blocks[8]]] += af * K_mat[2][2] * volp;
                }

                // --- For Generalized-α (if α_f = 1, back to Newmark-β) ---
                std::array<double, 3> nfint, nfext;
                nfint = this->ComputeInternalForce(ni, pid, dsf, sts_af);
                nfext = this->ComputeExternalForce(pid, sfi);

                mat.b_rhs[nid + offsets[0]] += nfint[0] + nfext[0];
                mat.b_rhs[nid + offsets[1]] += nfint[1] + nfext[1];
                mat.b_rhs[nid + offsets[2]] += nfint[2] + nfext[2];
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(mat.b_rhs, offsets);

    std::vector<double> naccel_am(nodec * 3);
    // --- Transfer acceleration to (n+α_m) time interlevel, if α_m = 1, back to Newmark-β) ---
    for (int n = 0; n < nodec * 3; n++) { naccel_am[n] = am * naccel_k[n] + am0 * this->naccel[n]; }
    // --------------------------------------

    this->AddInertialForceToRHS(mat, naccel_am, offsets);

    return;
}

void ImplicitSolidMPM::ImplicitDsfCorr(const std::vector<int> &nc, int nenode, //
                                    std::vector<std::array<double, 3>> &dsf) const noexcept {

    std::array<std::array<double, 3>, 3> def_grad{};
    def_grad = this->ComputeDeltaDefGrad(nc, nenode, this->implicit_.integrator_.alpha_f, dsf);

    std::array<std::array<double, 3>, 3> def_grad_inv;
    def_grad_inv = InvMat3(def_grad);

    std::array<double, 3> dsf_temp{};

    for (int ni = 0; ni < nenode; ni++) {
        for (int j = 0; j < 3; j++) {
            dsf_temp[j] = 0.0e0;
            for (int i = 0; i < 3; i++) { dsf_temp[j] += dsf[ni][i] * def_grad_inv[i][j]; }
        }
        dsf[ni] = dsf_temp;
    }

    return;
}
