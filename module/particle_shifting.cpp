#include "dataset.h"
#include "material_point.h"
#include "mpi_data.h"
#include "shape_function.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <mpi.h>
#include <vector>

std::vector<std::array<double, 3>> MaterialPoint::DeltaCorrectionParticleShifting() const {
    std::vector<double> nei(nodec, 0.0e0);
    std::vector<std::array<double, 3>> delta_corr;
    VectorAssign(this->num, delta_corr);

    double eu_norm = 0.0e0;
    for (int i = 0; i < nodec; i++) {
        nei[i] = std::max(0.0e0, -nvol[i] + this->nvof[i]);
        eu_norm += (nei[i] * nei[i]) * dbc[i];
    }

    MPI_Allreduce(MPI_IN_PLACE, &eu_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;
    std::vector<std::array<double, 3>> geup(this->num);
    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                geup[pid][0] += dsfi1 * nei[nid];
                geup[pid][1] += dsfi2 * nei[nid];
                geup[pid][2] += dsfi3 * nei[nid];
            }
            pid = this->idp2p[pid];
        }
    }

    double geup_dot = 0.0e0;
    for (int n = 0; n < this->num; n++) {
        geup[n][0] *= 2.0e0 * this->vol[n];
        geup[n][1] *= 2.0e0 * this->vol[n];
        geup[n][2] *= 2.0e0 * this->vol[n];
        for (int i = 0; i < 3; i++) { geup_dot += geup[n][i] * geup[n][i]; }
    }
    MPI_Allreduce(MPI_IN_PLACE, &geup_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double b_0 = (geup_dot < 1.0e-30) ? 0.0e0 : (eu_norm / geup_dot);

    for (int n = 0; n < this->num; n++) {
        delta_corr[n][0] = -b_0 * geup[n][0];
        delta_corr[n][1] = -b_0 * geup[n][1];
        delta_corr[n][2] = -b_0 * geup[n][2];
    }

    return delta_corr;
}

std::vector<std::array<double, 3>> MaterialPoint::SPHLikeParticleShifting() {

    std::vector<std::array<double, 3>> spring_force;
    VectorAssign(this->num, spring_force);

    // --- adjacent ghost point communication ---
    std::vector<int> mxy(3), checkmp;
    checkmp.assign(this->num, -1);
    this->par_comm_.nsp.assign(isb, 0);
    this->par_comm_.nrp.assign(isb, 0);
    for (int ip = 0; ip < this->num; ip++) {
        double radius = std::cbrt(0.75e0 * this->vol[ip] / M_PI);
        for (int i = 0; i < 3; i++) {
            if (std::abs(this->coord[ip][i] - xymin[i]) < radius) {
                mxy[i] = -1;
            } else if (std::abs(this->coord[ip][i] - xymax[i]) < radius) {
                mxy[i] = 1;
            } else {
                mxy[i] = 0;
            }
        }

        int idsb = myrank + mxy[0] + mxy[1] * nxyr[0] + mxy[2] * nxyr[0] * nxyr[1];

        for (int i = 0; i < isb; i++) {
            if (naid[i] == idsb) {
                this->par_comm_.nsp[i]++;
                checkmp[ip] = i;
                break;
            }
        }
    }

    std::vector<int> idmploc(isb);
    this->par_comm_.nmps = 0;
    for (int i = 0; i < isb; i++) {
        idmploc[i] = this->par_comm_.nmps;
        this->par_comm_.nmps += this->par_comm_.nsp[i];
    }

    VectorAssign(this->par_comm_.nmps, this->par_comm_.idmp);

    for (int ip = 0; ip < this->num; ip++) {
        int i = checkmp[ip];
        if (i != -1) { // Move
            this->par_comm_.idmp[idmploc[i]++] = ip;
        }
    }

    std::vector<MPI_Request> irqs(isb);
    std::vector<MPI_Request> irqr(isb);
    MPI_Status status;
    for (int i = 0; i < isb; i++) {
        int ncomid = naid[i];
        MPI_Isend(&this->par_comm_.nsp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqs[i]);
        MPI_Irecv(&this->par_comm_.nrp[i], 1, MPI_INT, ncomid, 1, MPI_COMM_WORLD, &irqr[i]);
    }
    for (int i = 0; i < isb; i++) {
        MPI_Wait(&irqs[i], &status);
        MPI_Wait(&irqr[i], &status);
    }

    this->par_comm_.nrps = 0;
    for (int i = 0; i < isb; i++) { this->par_comm_.nrps += this->par_comm_.nrp[i]; }

    MaterialPoint ghost_point;

    ghost_point.num = this->par_comm_.nrps;
    VectorAssign(ghost_point.num, ghost_point.coord);

    std::vector<std::array<double, 3>> bufs(this->par_comm_.nmps), bufr(this->par_comm_.nrps);
    int sip = 0;
    for (int i = 0; i < isb; i++) {
        for (int ip = 0; ip < this->par_comm_.nsp[i]; ip++) {
            int ipsip = ip + sip;
            bufs[ipsip] = this->coord[this->par_comm_.idmp[ipsip]];
        }
        sip += this->par_comm_.nsp[i];
    }

    this->par_comm_.PointVarSendrecv(bufs, bufr, 1);

    ghost_point.coord = bufr;
    // ------

    // Combine local and ghost coordinates into a single neighbor source list.
    std::vector<std::array<double, 3>> all_coords = this->coord;
    all_coords.insert(all_coords.end(), ghost_point.coord.begin(), ghost_point.coord.end());
    const int total_neighbors = this->num + ghost_point.num;

    for (int ip = 0; ip < this->num; ip++) {
        const double radius = std::cbrt(0.75e0 * this->vol[ip] / M_PI);
        const double radius_sq = radius * radius;
        const std::array<double, 3> &pi = this->coord[ip];
        std::array<double, 3> sum{0.0e0, 0.0e0, 0.0e0};

        for (int jp = 0; jp < total_neighbors; jp++) {
            if (jp == ip) { continue; }
            const std::array<double, 3> &pj = all_coords[jp];
            const double dx = pj[0] - pi[0];
            const double dy = pj[1] - pi[1];
            const double dz = pj[2] - pi[2];
            const double norm_sq = dx * dx + dy * dy + dz * dz;
            if (norm_sq > 0.0e0 && norm_sq < radius_sq) {
                const double norm = std::sqrt(norm_sq);
                const double smooth_weight = (norm < mtol) ? 0.0e0 : (1.0e0 - norm_sq / radius_sq);
                const double inv_norm = 1.0e0 / norm;
                sum[0] += dx * inv_norm * smooth_weight;
                sum[1] += dy * inv_norm * smooth_weight;
                sum[2] += dz * inv_norm * smooth_weight;
            }
        }

        const double a = 50.0e0;
        const double scale = -a * radius * dt;
        spring_force[ip][0] = scale * sum[0];
        spring_force[ip][1] = scale * sum[1];
        spring_force[ip][2] = scale * sum[2];
    }

    return spring_force;
}
