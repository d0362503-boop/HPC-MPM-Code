#include "module/relaxation.h"

#include <mpi.h>

#include "module/bc.h"
#include "module/contact.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

void Anderson_relaxation_M2(int iteration, const std::vector<double> &nodal_var, std::vector<double> &var_old,
                            std::vector<double> &var_older, std::vector<double> &res_old,
                            std::vector<double> &res_older, BoundaryCondition &intf_bc) {
    int num = intf_bc.ibc;
    std::vector<double> var_new(num * 3, 0.0e0), res_new(num * 3, 0.0e0);
    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
        var_new[i + 0 * num] = nodal_var[nid + nuc];
        var_new[i + 1 * num] = nodal_var[nid + nvc];
        var_new[i + 2 * num] = nodal_var[nid + nwc];

        res_new[i + 0 * num] = var_new[i + 0 * num] - intf_bc.fbc[i + 0 * num];
        res_new[i + 1 * num] = var_new[i + 1 * num] - intf_bc.fbc[i + 1 * num];
        res_new[i + 2 * num] = var_new[i + 2 * num] - intf_bc.fbc[i + 2 * num];
    }

    // F0 = R_n - R_{n-1}, F1 = R_{n-1} - R_{n-2}
    std::vector<double> F0(num * 3), F1(num * 3);
    for (int i = 0; i < num * 3; i++) {
        F0[i] = res_new[i] - res_old[i];
        F1[i] = res_old[i] - res_older[i];
    }

    std::vector<double> M(4, 0.0e0); // M[0]=M00, M[1]=M01, M[2]=M10, M[3]=M11
    std::vector<double> rhs(2, 0.0e0);

    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
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
            rhs[0] += F0[idx] * res_new[idx] * w;
            rhs[1] += F1[idx] * res_new[idx] * w;
        }
    }

    std::vector<double> M_all(4, 0.0e0), rhs_all(2, 0.0e0);
    MPI_Allreduce(M.data(), M_all.data(), 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(rhs.data(), rhs_all.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double det = M_all[0] * M_all[3] - M_all[1] * M_all[2];
    double lambda1 = 0.5e0, lambda2 = 0.5e0;
    if (iteration >= 2) {
        lambda1 = (M_all[3] * rhs_all[0] - M_all[1] * rhs_all[1]) / det;
        lambda2 = (M_all[0] * rhs_all[1] - M_all[2] * rhs_all[0]) / det;
    }

    for (int i = 0; i < num * 3; i++) {
        intf_bc.fbc[i] = var_new[i] - lambda1 * (var_new[i] - var_old[i]) //
                         - lambda2 * (var_old[i] - var_older[i]);
    }

    for (int i = 0; i < num * 3; i++) {
        var_older[i] = var_old[i];
        var_old[i] = var_new[i];
        res_older[i] = res_old[i];
        res_old[i] = res_new[i];
    }

    return;
}

void Anderson_relaxation_M1(int iteration, std::vector<double> &var_old, const std::vector<double> &nodal_var,
                            std::vector<double> &res_old, BoundaryCondition &intf_bc) {
    int num = intf_bc.ibc;
    std::vector<double> var_new(num * 3, 0.0e0), res_new(num * 3, 0.0e0);
    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
        var_new[i + 0 * num] = nodal_var[nid + nuc];
        var_new[i + 1 * num] = nodal_var[nid + nvc];
        var_new[i + 2 * num] = nodal_var[nid + nwc];

        res_new[i + 0 * num] = var_new[i + 0 * num] - intf_bc.fbc[i + 0 * num];
        res_new[i + 1 * num] = var_new[i + 1 * num] - intf_bc.fbc[i + 1 * num];
        res_new[i + 2 * num] = var_new[i + 2 * num] - intf_bc.fbc[i + 2 * num];
    }

    // F0 = R_n - R_{n-1}
    std::vector<double> F0(num * 3);
    for (int i = 0; i < num * 3; i++) { F0[i] = res_new[i] - res_old[i]; }

    double M = 0.0e0, rhs = 0.0e0;
    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
        for (int d = 0; d < 3; ++d) {
            int idx = i + d * num;
            double w = 0.0e0;
            if (d == 0) w = dbc[nid + nuc];
            if (d == 1) w = dbc[nid + nvc];
            if (d == 2) w = dbc[nid + nwc];

            M += F0[idx] * F0[idx] * w;
            rhs += F0[idx] * res_new[idx] * w;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &M, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &rhs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double lambda = 0.5e0;
    if (iteration >= 1) { lambda = rhs / M; }

    for (int i = 0; i < num * 3; i++) {
        intf_bc.fbc[i] = var_new[i] - lambda * (var_new[i] - var_old[i]);
        var_old[i] = var_new[i];
        res_old[i] = res_new[i];
    }

    return;
}

void Aitken_relaxation(int iteration, double &relax_coef, const std::vector<double> &nodal_var,
                       std::vector<double> &res_old, BoundaryCondition &intf_bc) {
    const double relax_coef_min = 0.05e0;
    const double relax_coef_max = 0.8e0;

    int num = intf_bc.ibc;
    std::vector<double> res_new(num * 3, 0.0e0);

    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
        res_new[i + 0 * num] = nodal_var[nid + nuc] - intf_bc.fbc[i + 0 * num];
        res_new[i + 1 * num] = nodal_var[nid + nvc] - intf_bc.fbc[i + 1 * num];
        res_new[i + 2 * num] = nodal_var[nid + nwc] - intf_bc.fbc[i + 2 * num];
    }

    double nume = 0.0e0, deno = 0.0e0;
    for (int i = 0; i < num; i++) {
        int nid = intf_bc.nbc[i];
        // --- Δr_i = r^{k+1}_i - r^k_i ---
        double dr1 = res_new[i + 0 * num] - res_old[i + 0 * num];
        double dr2 = res_new[i + 1 * num] - res_old[i + 1 * num];
        double dr3 = res_new[i + 2 * num] - res_old[i + 2 * num];
        // --- ∑ r^k_i · Δr_i ---
        nume += (res_old[i + 0 * num] * dr1) * dbc[nid + nuc]   //
                + (res_old[i + 1 * num] * dr2) * dbc[nid + nvc] //
                + (res_old[i + 2 * num] * dr3) * dbc[nid + nwc];
        // --- ∑ |Δr_i|^2 ---
        deno += pow(dr1, 2) * dbc[nid + nuc]   //
                + pow(dr2, 2) * dbc[nid + nvc] //
                + pow(dr3, 2) * dbc[nid + nwc];
    }

    MPI_Allreduce(MPI_IN_PLACE, &nume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &deno, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    relax_coef = 0.1e0;
    if (iteration >= 1) {
        relax_coef *= -(nume / deno);
        relax_coef = std::clamp(relax_coef, relax_coef_min, relax_coef_max); // --- under constrain ---
    }

    for (int i = 0; i < num * 3; i++) {
        intf_bc.fbc[i] += relax_coef * res_new[i];
        res_old[i] = res_new[i];
    }

    return;
}
