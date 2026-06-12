#include "relaxation.h"

#include <mpi.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "bc.h"
#include "contact.h"
#include "dataset.h"
#include "material_point.h"
#include "mesh.h"
#include "mpi_data.h"
#include "shape_function.h"

void Anderson_relaxation_M2(int block_it, const vector<double> &nvel_k, vector<double> &u_s_old,
                            vector<double> &u_s_older, vector<double> &r_k_old, vector<double> &r_k_older,
                            BoundaryCondition &fsi_intf) {
    int num = fsi_intf.ibc;
    vector<double> u_s_new(num * 3, 0.0e0), r_k_new(num * 3, 0.0e0);
    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        u_s_new[i + 0 * num] = nvel_k[nid + nuc];
        u_s_new[i + 1 * num] = nvel_k[nid + nvc];
        u_s_new[i + 2 * num] = nvel_k[nid + nwc];

        r_k_new[i + 0 * num] = u_s_new[i + 0 * num] - fsi_intf.fbc[i + 0 * num];
        r_k_new[i + 1 * num] = u_s_new[i + 1 * num] - fsi_intf.fbc[i + 1 * num];
        r_k_new[i + 2 * num] = u_s_new[i + 2 * num] - fsi_intf.fbc[i + 2 * num];
    }

    // F0 = R_n - R_{n-1}, F1 = R_{n-1} - R_{n-2}
    vector<double> F0(num * 3), F1(num * 3);
    for (int i = 0; i < num * 3; i++) {
        F0[i] = r_k_new[i] - r_k_old[i];
        F1[i] = r_k_old[i] - r_k_older[i];
    }

    vector<double> M(4, 0.0e0); // M[0]=M00, M[1]=M01, M[2]=M10, M[3]=M11
    vector<double> rhs(2, 0.0e0);

    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        for (int d = 0; d < 3; ++d) {
            int idx = i + d * num;
            double w = 0.0e0;
            if (d == 0) w = dbc[nid + nuc];
            if (d == 1) w = dbc[nid + nvc];
            if (d == 2) w = dbc[nid + nwc];

            M[0] += F0[idx] * F0[idx] * w; // M[0][0]
            M[1] += F0[idx] * F1[idx] * w; // M[0][1]
            M[2] += F1[idx] * F0[idx] * w; // M[1][0]
            M[3] += F1[idx] * F1[idx] * w; // M[1][1]
            rhs[0] += F0[idx] * r_k_new[idx] * w;
            rhs[1] += F1[idx] * r_k_new[idx] * w;
        }
    }

    vector<double> M_all(4, 0.0e0), rhs_all(2, 0.0e0);
    MPI_Allreduce(M.data(), M_all.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(rhs.data(), rhs_all.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double det = M_all[0] * M_all[3] - M_all[1] * M_all[2];
    double lambda1 = 0.5e0, lambda2 = 0.5e0;
    if (block_it >= 2) {
        lambda1 = (M_all[3] * rhs_all[0] - M_all[1] * rhs_all[1]) / det;
        lambda2 = (M_all[0] * rhs_all[1] - M_all[2] * rhs_all[0]) / det;
    }

    for (int i = 0; i < num * 3; i++) {
        fsi_intf.fbc[i] = u_s_new[i] - lambda1 * (u_s_new[i] - u_s_old[i]) //
                          - lambda2 * (u_s_old[i] - u_s_older[i]);
    }

    for (int i = 0; i < num * 3; i++) {
        u_s_older[i] = u_s_old[i];
        u_s_old[i] = u_s_new[i];
        r_k_older[i] = r_k_old[i];
        r_k_old[i] = r_k_new[i];
    }

    return;
}

void Anderson_relaxation_M1(int block_it, vector<double> &u_s_old, const vector<double> &nvel_k,
                            vector<double> &r_k_old, BoundaryCondition &fsi_intf) {
    int num = fsi_intf.ibc;
    vector<double> u_s_new(num * 3, 0.0e0), r_k_new(num * 3, 0.0e0);
    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        u_s_new[i + 0 * num] = nvel_k[nid + nuc];
        u_s_new[i + 1 * num] = nvel_k[nid + nvc];
        u_s_new[i + 2 * num] = nvel_k[nid + nwc];

        r_k_new[i + 0 * num] = u_s_new[i + 0 * num] - fsi_intf.fbc[i + 0 * num];
        r_k_new[i + 1 * num] = u_s_new[i + 1 * num] - fsi_intf.fbc[i + 1 * num];
        r_k_new[i + 2 * num] = u_s_new[i + 2 * num] - fsi_intf.fbc[i + 2 * num];
    }

    // F0 = R_n - R_{n-1}
    vector<double> F0(num * 3);
    for (int i = 0; i < num * 3; i++) { F0[i] = r_k_new[i] - r_k_old[i]; }

    double M = 0.0e0, rhs = 0.0e0;
    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        for (int d = 0; d < 3; ++d) {
            int idx = i + d * num;
            double w = 0.0e0;
            if (d == 0) w = dbc[nid + nuc];
            if (d == 1) w = dbc[nid + nvc];
            if (d == 2) w = dbc[nid + nwc];

            M += F0[idx] * F0[idx] * w;
            rhs += F0[idx] * r_k_new[idx] * w;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &rhs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double lambda = 0.5e0;
    if (block_it >= 1) { lambda = rhs / M; }

    for (int i = 0; i < num * 3; i++) {
        fsi_intf.fbc[i] = u_s_new[i] - lambda * (u_s_new[i] - u_s_old[i]);
        u_s_old[i] = u_s_new[i];
        r_k_old[i] = r_k_new[i];
    }

    return;
}

void Aitken_relaxation(int block_it, double &omega, const vector<double> &nvel_k, vector<double> &r_k_old,
                       BoundaryCondition &fsi_intf) {
    const double omega_min = 0.05e0;
    const double omega_max = 0.8e0;

    int num = fsi_intf.ibc;
    vector<double> r_k_new(num * 3, 0.0e0);

    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        r_k_new[i + 0 * num] = nvel_k[nid + nuc] - fsi_intf.fbc[i + 0 * num];
        r_k_new[i + 1 * num] = nvel_k[nid + nvc] - fsi_intf.fbc[i + 1 * num];
        r_k_new[i + 2 * num] = nvel_k[nid + nwc] - fsi_intf.fbc[i + 2 * num];
    }

    double nume = 0.0e0, deno = 0.0e0;
    for (int i = 0; i < num; i++) {
        int nid = fsi_intf.nbc[i];
        // --- Δr_i = r^{k+1}_i - r^k_i ---
        double dr1 = r_k_new[i + 0 * num] - r_k_old[i + 0 * num];
        double dr2 = r_k_new[i + 1 * num] - r_k_old[i + 1 * num];
        double dr3 = r_k_new[i + 2 * num] - r_k_old[i + 2 * num];
        // --- ∑ r^k_i · Δr_i ---
        nume += (r_k_old[i + 0 * num] * dr1) * dbc[nid + nuc]   //
                + (r_k_old[i + 1 * num] * dr2) * dbc[nid + nvc] //
                + (r_k_old[i + 2 * num] * dr3) * dbc[nid + nwc];
        // --- ∑ |Δr_i|^2 ---
        deno += pow(dr1, 2) * dbc[nid + nuc]   //
                + pow(dr2, 2) * dbc[nid + nvc] //
                + pow(dr3, 2) * dbc[nid + nwc];
    }

    MPI_Allreduce(MPI_IN_PLACE, &nume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &deno, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    omega = 0.1e0;
    if (block_it >= 1) {
        omega *= -(nume / deno);
        // omega = min(max(omega, omega_min), omega_max);   // --- under constrain ---
    }

    for (int i = 0; i < num * 3; i++) {
        fsi_intf.fbc[i] += omega * r_k_new[i];
        r_k_old[i] = r_k_new[i];
    }

    return;
}
