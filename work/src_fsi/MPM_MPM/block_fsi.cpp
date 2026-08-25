#include "work/src_fsi/MPM_MPM/block_fsi.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "module/bc.h"
#include "module/contact.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/map_and_interpolate.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/relaxation.h"
#include "module/shape_function.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solid/solid_material_point.h"
#include "module/solver/crsmat.h"

void MPMMPMBlockFSI::DetectFSIInterface() {

    VectorAssign(nodec, this->fsi_intf.nbc);
    VectorAssign(nodec, this->fsi_intf.fbc);
    VectorAssign(nodec, this->solid_.nphi);
    VectorAssign(nodec, this->fluid_.nphi);

    this->fsi_intf.ibc = 0;
    for (int n = 0; n < nodec; n++) {
        this->solid_.nphi[n] = this->solid_.nvof[n] / nvol[n];
        this->solid_.nphi[n] = std::clamp(this->solid_.nphi[n], 0.0e0, 1.0e0);

        this->fluid_.nphi[n] = this->fluid_.nvof[n] / nvol[n];
        this->fluid_.nphi[n] = std::clamp(this->fluid_.nphi[n], 0.0e0, 1.0e0);
        if (this->solid_.nphi[n] > this->phi_cut && //
            this->fluid_.nphi[n] > this->phi_cut) {
            this->fsi_intf.nbc[this->fsi_intf.ibc++] = n;
        }
    }

    return;
}

void MPMMPMBlockFSI::SolveFSISystem() {

    this->DetectFSIInterface();

    int num = this->fsi_intf.ibc;
    double rtr_fsi_ref, r0r_fsi_ref, rkr_fsi_ref, rtr_fsi_abs;
    std::vector<double> r_k_old(num * 3, 0.0e0), u_s_old(num * 3, 0.0e0);

    for (int block_it = 0; block_it <= this->max_block_iter; block_it++) {

        this->solid_.SolveSolid();

        Anderson_relaxation_M1(block_it, u_s_old, this->solid_.ndispl, r_k_old, this->fsi_intf);
        // Aitken_relaxation(block_it, this->relax_omega, this->solid_.ndispl, r_k_old, this->fsi_intf);

        this->fluid_.SolveNS();

        this->CalFSIForce();

        if (block_it == 0) {
            this->CalFSIResidual(r0r_fsi_ref, rtr_fsi_abs);
            rkr_fsi_ref = r0r_fsi_ref;
        } else {
            this->CalFSIResidual(rkr_fsi_ref, rtr_fsi_abs);
        }
        rtr_fsi_ref = (r0r_fsi_ref > 1.0e-12) ? (rkr_fsi_ref / r0r_fsi_ref) : 0.0e0;

        // --- Convergence check ---
        if (rtr_fsi_ref < this->tol_ref || r0r_fsi_ref < this->tol_ref || //
            rtr_fsi_abs < this->tol_abs) {
            if (myrank == 0) {
                std::cout << "Block_converge:" << std::setw(10) << block_it  //
                          << std::setw(15) << std::scientific << rtr_fsi_ref //
                          << std::setw(15) << std::scientific << r0r_fsi_ref //
                          << std::setw(15) << std::scientific << rtr_fsi_abs << "\n";
            }
            break;
        } else if (block_it == this->max_block_iter) {
            if (myrank == 0) {
                std::cout << "Block_it:" << std::setw(10) << block_it << "\n";
                std::cout << "Block iteration can not converge" << "\n";
            }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    return;
}

void MPMMPMBlockFSI::CalFSIResidual(double &rtr_ref, double &rtr_dof) {

    double norm = 0.0e0, intf_num = 0.0e0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        int nid = this->fsi_intf.nbc[n];

        double us1 = this->solid_.ndispl[nid + nuc];
        double us2 = this->solid_.ndispl[nid + nvc];
        double us3 = this->solid_.ndispl[nid + nwc];

        double uf1 = this->fluid_.ndispl[nid + nuc];
        double uf2 = this->fluid_.ndispl[nid + nvc];
        double uf3 = this->fluid_.ndispl[nid + nwc];

        double dr1 = us1 - uf1;
        double dr2 = us2 - uf2;
        double dr3 = us3 - uf3;

        norm += dr1 * dr1 * dbc[nid + nuc]   //
                + dr2 * dr2 * dbc[nid + nvc] //
                + dr3 * dr3 * dbc[nid + nwc];

        intf_num += dbc[nid];
    }

    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &intf_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    rtr_ref = std::sqrt(norm);
    rtr_dof = std::sqrt(norm / intf_num);

    return;
}

void MPMMPMBlockFSI::CalFSIForce() {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);
    this->fluid_.ComputeNodeVelAccelFromDispl(nvel_k, naccel_k);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    const double af = this->fluid_.alpha_f;
    const double af0 = 1.0e0 - af;

    std::vector<double> nvel_af(nodec * 3), npres_af(nodec);
    for (int n = 0; n < nodec * 3; n++) { //
        nvel_af[n] = af0 * this->fluid_.nvel[n] + af * nvel_k[n];
    }
    for (int n = 0; n < nodec; n++) { //
        npres_af[n] = af0 * this->fluid_.npres_old[n] + af * this->fluid_.npres[n];
    }
    // -----------------------------------------------------

    double g_weight = (dxy[0] * dxy[1] * dxy[2]) / (npxye[0] * npxye[1] * npxye[2]);

    VectorAssign(nodec * 3, this->nfsi_force);
    for (int m = 0; m < nelem; m++) {
        const std::array<int, 3> ijk = IndexToIJK(m, xyelem);

        std::array<double, 3> xye, xyp;
        xye[0] = xymin[0] + dxy[0] * (double(ijk[0]) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(ijk[1]) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ijk[2]) + 0.5e0);
        for (int iz = 0; iz < npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyp[0] = xye[0] + dec2p[ix][0];
                    MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

                    double pres_k = 0.0e0;
                    std::array<double, 3> grad_phi_k{};
                    std::array<std::array<double, 3>, 3> grad_vel_k{};
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        double dsfi1 = dsf[ni][0];
                        double dsfi2 = dsf[ni][1];
                        double dsfi3 = dsf[ni][2];
                        // --- Pressure ---
                        pres_k += sfi * npres_af[nid];
                        // --- Gradient of solid phi ---
                        grad_phi_k[0] += dsfi1 * this->solid_.nphi[nid];
                        grad_phi_k[1] += dsfi2 * this->solid_.nphi[nid];
                        grad_phi_k[2] += dsfi3 * this->solid_.nphi[nid];
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

                    // --- Stress ---
                    std::array<std::array<double, 3>, 3> sts_k{};
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            sts_k[i][j] = this->fluid_.rmu * (grad_vel_k[i][j] + grad_vel_k[j][i]);
                        }
                        sts_k[i][i] -= pres_k;
                    }

                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        double dsfi1 = dsf[ni][0];
                        double dsfi2 = dsf[ni][1];
                        double dsfi3 = dsf[ni][2];

                        this->nfsi_force[nid + nuc] +=
                            g_weight * sfi * //
                            (grad_phi_k[0] * sts_k[0][0] + grad_phi_k[1] * sts_k[0][1] + grad_phi_k[2] * sts_k[0][2]);
                        this->nfsi_force[nid + nvc] +=
                            g_weight * sfi * //
                            (grad_phi_k[0] * sts_k[1][0] + grad_phi_k[1] * sts_k[1][1] + grad_phi_k[2] * sts_k[1][2]);
                        this->nfsi_force[nid + nwc] +=
                            g_weight * sfi * //
                            (grad_phi_k[0] * sts_k[2][0] + grad_phi_k[1] * sts_k[2][1] + grad_phi_k[2] * sts_k[2][2]);
                    }
                }
            }
        }
    }

    NodeVarComm(this->nfsi_force, {nuc, nvc, nwc});

    for (int n = 0; n < nodec; n++) {
        if (this->solid_.nphi[n] < this->phi_cut || this->fluid_.nphi[n] < this->phi_cut) {
            this->nfsi_force[n + nuc] = 0.0e0;
            this->nfsi_force[n + nvc] = 0.0e0;
            this->nfsi_force[n + nwc] = 0.0e0;
        }
    }

    return;
}

void FSISolid::AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel) {

    implicitmpm::ImplicitSolidMPM::AddInertialForceToRHS(mat, naccel);

    for (int n = 0; n < nodec * 3; n++) {
        // --- FSI force ---
        mat.b_rhs[n] -= this->fsi_.nfsi_force[n];
    }

    return;
}
