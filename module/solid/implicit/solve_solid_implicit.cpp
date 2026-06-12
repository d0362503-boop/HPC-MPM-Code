#include <mpi.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "../../bc.h"
#include "../../cal_mat.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "../../solver/crsmat.h"
#include "../../solver/solver.h"
#include "../solid_material_point.h"
#include "implicit_mpm_solid.h"

using namespace implicitmpm;

void SolidMaterialPoint::UpdateNRIncrement() {
    for (int n = 0; n < nodec * 3; n++) { this->ndispl[n] += this->SM_.x_lhs[n]; }

    return;
}

std::vector<std::array<double, 6>> SolidMaterialPoint::InitializeNRStress() { return this->stress; }

auto SolidMaterialPoint::ComputeTangentModulus(int pid, int ni, int nj, const std::vector<std::array<double, 3>> &dsf,
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

void SolidMaterialPoint::SolveSolid() {
    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);

    // --- NR_flag true for nonlinear elasticity and false for linear elasticity ---
    int iter_max = (NR_flag) ? 1000 : 0;

    std::vector<std::array<double, 6>> stress_k = this->InitializeNRStress(); // --- Reset stress for NR ---

    VectorAssign(nodec * 3, this->ndispl);    // ---- Reset displacement ----
    VectorAssign(nodec * 3, this->SM_.x_lhs); // ---- Initialize LHS x value ----
    double r0r = 0.0e0;
    for (int NR_it = 0; NR_it <= iter_max; NR_it++) {
        this->BCNRSet();

        this->PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k); // ---- Newmark beta velocity & acceleration ----

        this->AssembleSystem(naccel_k, nvel_k, stress_k);

        int iter = this->SolveSystem(this->SM_, NR_it);

        this->UpdateNRIncrement();

        if (!NR_flag) { // --- Linear not need Newton Raphson ---
            if (myrank == 0) {
                std::cout << "Solid_converge:" << std::setw(10) << NR_it << std::setw(10) << iter << "\n";
            }
        } else { // --- Nonlinear need Newton Raphson ---
            if (this->SM_.CheckNRConvergence(NR_it, iter_max, iter, r0r)) { break; }
        }
    }

    // --- Memory clear ---
    std::vector<double>().swap(this->SM_.adiag);
    std::vector<double>().swap(this->SM_.amat);
    std::vector<double>().swap(this->SM_.b_rhs);
    std::vector<double>().swap(this->SM_.x_lhs);

    return;
}

void SolidMaterialPoint::AssembleSystem(const std::vector<double> &naccel_k, const std::vector<double> &nvel_k,
                                        std::vector<std::array<double, 6>> &stress_k) {
    const double af = this->alpha_f;
    const double af0 = 1.0e0 - this->alpha_f;
    const double am = this->alpha_m;
    const double am0 = 1.0e0 - this->alpha_m;

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    std::vector<std::array<std::array<double, 3>, 3>> delta_def_grad, def_grad_NR;

    VectorAssign(this->num, delta_def_grad);
    VectorAssign(this->num, def_grad_NR);
    VectorAssign(this->SM_.nmata, this->SM_.amat);
    VectorAssign(nodec * 3, this->SM_.adiag);
    VectorAssign(nodec * 3, this->SM_.b_rhs);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];

            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            this->UpdateDefGrad(pid, nenode, this->alpha_f, ncm, sf, dsf, delta_def_grad, def_grad_NR);

            this->UpdateConstitutiveModel(pid, stress_k, this->det_def_grad, def_grad_NR, delta_def_grad);

            this->ImplicitDsfCorr(ncm, nenode, dsf);

            // --- Transfer stress to (n+α_f) time interlevel ---
            std::array<double, 6> sts_af = stress_k[pid];

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];

                // --- Mass Matrix (lumpd) ---
                double emd_lu = this->ComputeNRLumpedMat(pid, sfi);

                int ncol = 0;
                for (int nj = 0; nj < nenode; nj++) {
                    int njd = ncm[nj];
                    double sfj = sf[nj];
                    double dsfj1 = dsf[nj][0];
                    double dsfj2 = dsf[nj][1];
                    double dsfj3 = dsf[nj][2];

                    int ida = this->SM_.FindIndex(nid, njd, ncol);

                    std::array<std::array<double, 3>, 3> K_mat;
                    K_mat = this->ComputeTangentModulus(pid, ni, nj, dsf, sts_af);

                    if (nid == njd) {
                        this->SM_.amat[ida + this->SM_.block_id[0]] += am * emd_lu;
                        this->SM_.amat[ida + this->SM_.block_id[4]] += am * emd_lu;
                        this->SM_.amat[ida + this->SM_.block_id[8]] += am * emd_lu;
                    }
                    this->SM_.amat[ida + this->SM_.block_id[0]] += af * K_mat[0][0] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[1]] += af * K_mat[0][1] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[2]] += af * K_mat[0][2] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[3]] += af * K_mat[1][0] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[4]] += af * K_mat[1][1] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[5]] += af * K_mat[1][2] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[6]] += af * K_mat[2][0] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[7]] += af * K_mat[2][1] * this->vol[pid];
                    this->SM_.amat[ida + this->SM_.block_id[8]] += af * K_mat[2][2] * this->vol[pid];
                }

                // --- For Generalized-α (if α_f = 1, back to Newmark-β) ---
                std::array<double, 3> nfint, nfext;
                nfint = ComputeInternalForce(ni, pid, dsf, sts_af);
                nfext = ComputeExternalForce(pid, sfi);

                this->SM_.b_rhs[nid + nuc] += nfint[0] + nfext[0];
                this->SM_.b_rhs[nid + nvc] += nfint[1] + nfext[1];
                this->SM_.b_rhs[nid + nwc] += nfint[2] + nfext[2];
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->SM_.b_rhs, {nuc, nvc, nwc});

    // --- Transfer acceleration to (n+α_m) time interlevel ---
    std::vector<double> naccel_am(nodec * 3);
    for (int n = 0; n < nodec * 3; n++) { naccel_am[n] = am * naccel_k[n] + am0 * this->naccel[n]; }
    // --------------------------------------

    this->AddInertialForceToRHS(this->SM_, naccel_am);

    // ---- Diagional Scalling Prec ----
    if (!this->SM_.use_petsc) { this->SM_.BuildDiagonalPreconditioner(this->SM_.ndof); }

    return;
}
