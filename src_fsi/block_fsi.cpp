#include "block_fsi.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "../module/bc.h"
#include "../module/contact.h"
#include "../module/dataset.h"
#include "../module/fluid/FEM/stabilized_fem.h"
#include "../module/map_and_interpolate.h"
#include "../module/material_point.h"
#include "../module/mesh.h"
#include "../module/mpi_data.h"
#include "../module/relaxation.h"
#include "../module/shape_function.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"
#include "../module/solid/solid_material_point.h"
#include "../module/solver/crsmat.h"

void BlockFSI::DetectFSIInterface() {
    VectorAssign(nodec, this->fsi_intf.nbc);
    VectorAssign(nodec, this->solid_.nphi);
    this->fsi_intf.ibc = 0;
    for (int n = 0; n < nodec; n++) {
        this->solid_.nphi[n] = this->solid_.nvof[n] / nvol[n];
        if (this->solid_.nphi[n] > 0.25e0) { this->fsi_intf.nbc[this->fsi_intf.ibc++] = n; }
        if (this->solid_.nphi[n] < 0.01e0) this->solid_.nphi[n] = 0.0e0;
        if (this->solid_.nphi[n] > 0.99e0) this->solid_.nphi[n] = 1.0e0;
    }

    return;
}

void BlockFSI::UpdateFSIInterfaceBC() {
    int num = this->fsi_intf.ibc;
    if (num == 0) return;
    VectorAssign(num * 3, this->fsi_intf.fbc);
    for (int i = 0; i < num; i++) {
        int nid = this->fsi_intf.nbc[i];
        this->fsi_intf.fbc[i + 0 * num] = this->solid_.nvel[nid + nuc];
        this->fsi_intf.fbc[i + 1 * num] = this->solid_.nvel[nid + nvc];
        this->fsi_intf.fbc[i + 2 * num] = this->solid_.nvel[nid + nwc];
    }

    return;
}

void BlockFSI::SolveFSISystem() {
    this->DetectFSIInterface();

    int num = this->fsi_intf.ibc;
    double rtr_fsi_ref, r0r_fsi_ref, rkr_fsi_ref, rtr_fsi_abs;
    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3), //
        r_k_old(num * 3, 0.0e0), r_k_older(num * 3, 0.0e0),     //
        u_s_old(num * 3, 0.0e0), u_s_older(num * 3, 0.0e0);

    this->UpdateFSIInterfaceBC();

    this->fluid_.NS_.force_rebuild_next_ = true;
    this->solid_.SM_.force_rebuild_next_ = true;
    for (int block_it = 0; block_it <= this->max_block_iter; block_it++) {

        this->fluid_.SolveNS();

        this->CalFSIForce();

        this->solid_.SolveSolid();

        this->solid_.PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k);

        Anderson_relaxation_M1(block_it, u_s_old, nvel_k, r_k_old, this->fsi_intf);
        // Aitken_relaxation(block_it, this->relax_omega, nvel_k, r_k_old, this->fsi_intf);
        // Anderson_relaxation_M2(block_it, nvel_k, u_s_old, u_s_older, r_k_old, r_k_older, this->fsi_intf);

        if (block_it == 0) {
            this->CalFSIResidual(r0r_fsi_ref, rtr_fsi_abs, nvel_k);
            rkr_fsi_ref = r0r_fsi_ref;
        } else {
            this->CalFSIResidual(rkr_fsi_ref, rtr_fsi_abs, nvel_k);
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

    // this->CalDragLiftCoeff();

    this->fluid_.naccel = this->fluid_.GeneralizedAlphaNodeAccelUpdate();

    this->fluid_.UpdateNodalVar();

    return;
}

void BlockFSI::CalFSIResidual(double &rtr_ref, double &rtr_dof, const std::vector<double> &nvel_k) {
    double norm = 0.0e0, intf_num = 0.0e0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        int nid = this->fsi_intf.nbc[n];

        double us1 = nvel_k[nid + nuc];
        double us2 = nvel_k[nid + nvc];
        double us3 = nvel_k[nid + nwc];

        double uf1 = fluid_.nvel[nid + nuc];
        double uf2 = fluid_.nvel[nid + nvc];
        double uf3 = fluid_.nvel[nid + nwc];

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

void BlockFSI::CalDragLiftCoeff() {
    // --- For CD and CL coefficients ---
    double cd = 0.0e0, cl = 0.0e0;
    for (int n = 0; n < nodec; n++) {
        cd += this->nfsi_force[n + nuc] * dbc[n + nuc];
        cl += this->nfsi_force[n + nwc] * dbc[n + nwc];
    }

    cd = -cd;
    cl = -cl;
    MPI_Allreduce(MPI_IN_PLACE, &cd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &cl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    cd *= 2.0e0 / (this->fluid_.rhol * 0.1e0 * dxy[1]);
    cl *= 2.0e0 / (this->fluid_.rhol * 0.1e0 * dxy[1]);

    if (myrank == 0) {
        std::cout << "Drag and lift coefficients: " //
                  << std::setw(15) << cd            //
                  << std::setw(15) << cl << "\n";
    }

    return;
}

void BlockFSI::CalFSIForce() {
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->added_mass);
    for (int m = 0; m < nelem; m++) {
        int pid = this->solid_.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->solid_.coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                // --- water nodal mass ---
                this->added_mass[nid] += sfi * this->fluid_.rhol * this->solid_.vol[pid];
            }
            pid = this->solid_.idp2p[pid];
        }
    }

    NodeVarComm(this->added_mass, 0);

    // --- Transfet variables to (n+α_f) time interlevel ---
    const double af = this->fluid_.alpha_f;
    const double af0 = 1.0e0 - af;

    std::vector<double> nvel_af(nodec * 3), npres_af(nodec);
    for (int n = 0; n < nodec * 3; n++) {
        nvel_af[n] = af0 * this->fluid_.nvel_old[n] //
                     + af * this->fluid_.nvel[n];
    }
    for (int n = 0; n < nodec; n++) {
        npres_af[n] = af0 * this->fluid_.npres_old[n] //
                      + af * this->fluid_.npres[n];
    }
    // -----------------------------------------------------

    double G_Weight = (dxy[0] * dxy[1] * dxy[2]) / (npxye[0] * npxye[1] * npxye[2]);

    VectorAssign(nodec * 3, this->nfsi_force);
    for (int m = 0; m < nelem; m++) {
        int pid = this->fluid_.idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->fluid_.coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            double phi_k = 0.0e0, pres_k = 0.0e0;
            std::array<std::array<double, 3>, 3> grad_vel_k{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                // --- Interface funtion ---
                phi_k += sfi * this->fluid_.nphi[nid];
                // --- Pressure ---
                pres_k += sfi * npres_af[nid];
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

            double rhof = phi_k * this->fluid_.rhol + (1.0e0 - phi_k) * this->fluid_.rhog;
            double rmuf = phi_k * this->fluid_.rmul + (1.0e0 - phi_k) * this->fluid_.rmug;

            // --- Stress ---
            std::array<std::array<double, 3>, 3> sts_k{};
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) { sts_k[i][j] = rmuf * (grad_vel_k[i][j] + grad_vel_k[j][i]); }
                sts_k[i][i] -= pres_k;
            }

            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];

                this->nfsi_force[nid + nuc] += G_Weight * this->solid_.nphi[nid] *
                                               (dsfi1 * sts_k[0][0] + dsfi2 * sts_k[0][1] + dsfi3 * sts_k[0][2]);
                this->nfsi_force[nid + nvc] += G_Weight * this->solid_.nphi[nid] *
                                               (dsfi1 * sts_k[1][0] + dsfi2 * sts_k[1][1] + dsfi3 * sts_k[1][2]);
                this->nfsi_force[nid + nwc] += G_Weight * this->solid_.nphi[nid] *
                                               (dsfi1 * sts_k[2][0] + dsfi2 * sts_k[2][1] + dsfi3 * sts_k[2][2]);
            }
            pid = this->fluid_.idp2p[pid];
        }
    }

    NodeVarComm(this->nfsi_force, {nuc, nvc, nwc});

    for (int n = 0; n < nodec; n++) {
        if (this->solid_.nphi[n] < 0.25e0) {
            this->nfsi_force[n + nuc] = 0.0e0;
            this->nfsi_force[n + nvc] = 0.0e0;
            this->nfsi_force[n + nwc] = 0.0e0;
        }
    }

    return;
}

double FSISolid::ComputeNRLumpedMat(double pid, double sfi) {
    double emd = MaterialPoint::ComputeNRLumpedMat(pid, sfi);
    emd += this->nb_para[3] * sfi * this->fsi_.fluid_.rhol * this->vol[pid];

    return emd;
}

std::array<double, 3> FSISolid::ComputeExternalForce(int pid, double sfi) {
    double fx = bb[0] * facl;
    double fy = bb[1] * facl;
    double fz = bb[2] * facl;

    std::array<double, 3> nfext;
    nfext[0] = sfi * ((this->mass[pid] + this->vol[pid] * this->fsi_.fluid_.rhol) * fx + this->trac_force[pid][0]);
    nfext[1] = sfi * ((this->mass[pid] + this->vol[pid] * this->fsi_.fluid_.rhol) * fy + this->trac_force[pid][1]);
    nfext[2] = sfi * ((this->mass[pid] + this->vol[pid] * this->fsi_.fluid_.rhol) * fz + this->trac_force[pid][2]);

    return nfext;
}

void FSISolid::AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel) {
    for (int n = 0; n < nodec; n++) {
        // --- For Generalized-α (if α_m = 1, back to Newmark-β) ---
        mat.b_rhs[n + nuc] -= (this->nmass[n] + this->fsi_.added_mass[n]) * naccel[n + nuc];
        mat.b_rhs[n + nvc] -= (this->nmass[n] + this->fsi_.added_mass[n]) * naccel[n + nvc];
        mat.b_rhs[n + nwc] -= (this->nmass[n] + this->fsi_.added_mass[n]) * naccel[n + nwc];
    }

    for (int n = 0; n < nodec * 3; n++) {
        // --- FSI force ---
        mat.b_rhs[n] -= this->fsi_.nfsi_force[n];
    }

    return;
}
