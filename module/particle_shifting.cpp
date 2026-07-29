#include "DLB/mpm_dlb.h"
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

    double b_0 = (geup_dot > 1.0e-30) ? (eu_norm / geup_dot) : 0.0e0;

    std::vector<std::array<double, 3>> disp_corr;
    VectorAssign(this->num, disp_corr);
    for (int n = 0; n < this->num; n++) {
        disp_corr[n][0] = -b_0 * geup[n][0];
        disp_corr[n][1] = -b_0 * geup[n][1];
        disp_corr[n][2] = -b_0 * geup[n][2];
    }

    return disp_corr;
}

std::vector<std::array<double, 3>> MaterialPoint::PairwiseRepulsiveParticleShifting() {

    std::vector<std::array<double, 3>> disp_corr;
    VectorAssign(this->num, disp_corr);

    // -------------------------------------------------------------------------
    // 1. Ghost particle communication (per-particle support)
    // -------------------------------------------------------------------------
    const std::vector<mpm_dlb::Region> &regions = mpm_dlb::CurrentRegions();
    std::vector<std::array<double, 3>> peer_min(isb);
    std::vector<std::array<double, 3>> peer_max(isb);
    for (int i = 0; i < isb; ++i) {
        const mpm_dlb::Region &peer = regions[naid[i]];
        for (int d = 0; d < 3; ++d) {
            peer_min[i][d] = xyminw[d] + dxy[d] * double(peer.elem_min[d]);
            peer_max[i][d] = xyminw[d] + dxy[d] * double(peer.elem_max[d] + 1);
        }
    }

    std::vector<std::vector<int>> send_ids(isb);
    for (int ip = 0; ip < this->num; ip++) {
        const double support = std::cbrt(this->vol[ip]);

        for (int i = 0; i < isb; ++i) {
            bool intersects_peer = true;
            for (int d = 0; d < 3; ++d) {
                if (this->coord[ip][d] + support <= peer_min[i][d] || this->coord[ip][d] - support >= peer_max[i][d]) {
                    intersects_peer = false;
                    break;
                }
            }
            if (intersects_peer) { send_ids[i].push_back(ip); }
        }
    }

    this->par_comm_.comm_ranks = naid;
    this->par_comm_.nsp.assign(isb, 0);
    this->par_comm_.nrp.assign(isb, 0);
    this->par_comm_.nmps = 0;
    for (int i = 0; i < isb; i++) {
        this->par_comm_.nsp[i] = static_cast<int>(send_ids[i].size());
        this->par_comm_.nmps += this->par_comm_.nsp[i];
    }

    VectorAssign(this->par_comm_.nmps, this->par_comm_.idmp);

    int idmp_pos = 0;
    for (int i = 0; i < isb; i++) {
        for (int ip : send_ids[i]) { this->par_comm_.idmp[idmp_pos++] = ip; }
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

    std::vector<std::array<double, 3>> bufs(this->par_comm_.nmps);
    int sip = 0;
    for (int i = 0; i < isb; i++) {
        for (int ip = 0; ip < this->par_comm_.nsp[i]; ip++) {
            int ipsip = ip + sip;
            bufs[ipsip] = this->coord[this->par_comm_.idmp[ipsip]];
        }
        sip += this->par_comm_.nsp[i];
    }

    std::vector<std::array<double, 3>> bufr(this->par_comm_.nrps);
    this->par_comm_.PointVarSendrecv(bufs, bufr, 1);

    const int nghost = this->par_comm_.nrps;
    const int np_total = this->num + nghost;

    // -------------------------------------------------------------------------
    // 2. Spatial hash grid using dxy as cell size, handling quasi-3D
    // -------------------------------------------------------------------------
    std::array<double, 3> box_min{};
    std::array<double, 3> box_max{};
    std::array<int, 3> ncells{};
    for (int d = 0; d < 3; d++) {
        bool has_lower_ghost = (xymin[d] > xyminw[d]);
        bool has_upper_ghost = (xymax[d] < xymaxw[d]);

        box_min[d] = xymin[d] - (has_lower_ghost ? dxy[d] : 0.0e0);
        box_max[d] = xymax[d] + (has_upper_ghost ? dxy[d] : 0.0e0);

        ncells[d] = xyelem[d];
        if (has_lower_ghost) { ncells[d]++; }
        if (has_upper_ghost) { ncells[d]++; }
        if (ncells[d] < 1) { ncells[d] = 1; }
    }
    const int ncell_total = ncells[0] * ncells[1] * ncells[2];

    std::vector<int> cell_head(ncell_total, -1);
    std::vector<int> next_particle(np_total, -1);

    auto Particle2Cell = [&](const std::array<double, 3> &x) {
        int cx = static_cast<int>((x[0] - box_min[0]) / dxy[0]);
        int cy = static_cast<int>((x[1] - box_min[1]) / dxy[1]);
        int cz = static_cast<int>((x[2] - box_min[2]) / dxy[2]);

        cx = std::max(0, std::min(cx, ncells[0] - 1));
        cy = std::max(0, std::min(cy, ncells[1] - 1));
        cz = std::max(0, std::min(cz, ncells[2] - 1));

        return cx + cy * ncells[0] + cz * ncells[0] * ncells[1];
    };

    auto InsertParticle = [&](int pid, const std::array<double, 3> &x) {
        int cell_id = Particle2Cell(x);
        next_particle[pid] = cell_head[cell_id];
        cell_head[cell_id] = pid;
    };

    for (int ip = 0; ip < this->num; ip++) { InsertParticle(ip, this->coord[ip]); }
    for (int jp = 0; jp < nghost; jp++) { InsertParticle(this->num + jp, bufr[jp]); }

    // -------------------------------------------------------------------------
    // 3. 3x3x3 cell search with per-particle support
    // -------------------------------------------------------------------------
    const double gamma_s = 50.0e0;

    for (int ip = 0; ip < this->num; ip++) {
        const double support = std::cbrt(this->vol[ip]);
        const double support_sq = support * support;
        const std::array<double, 3> &pi = this->coord[ip];

        int base_cell = Particle2Cell(pi);
        int cx = base_cell % ncells[0];
        int tmp = base_cell / ncells[0];
        int cy = tmp % ncells[1];
        int cz = tmp / ncells[1];

        std::array<double, 3> sum{0.0e0, 0.0e0, 0.0e0};

        for (int dx = -1; dx <= 1; dx++) {
            int nx = cx + dx;
            if (nx < 0 || nx >= ncells[0]) { continue; }
            for (int dy = -1; dy <= 1; dy++) {
                int ny = cy + dy;
                if (ny < 0 || ny >= ncells[1]) { continue; }
                for (int dz = -1; dz <= 1; dz++) {
                    int nz = cz + dz;
                    if (nz < 0 || nz >= ncells[2]) { continue; }

                    int cell_id = nx + ny * ncells[0] + nz * ncells[0] * ncells[1];
                    int jp = cell_head[cell_id];

                    while (jp != -1) {
                        if (jp == ip) {
                            jp = next_particle[jp];
                            continue;
                        }

                        const std::array<double, 3> &pj = (jp < this->num) ? this->coord[jp] : bufr[jp - this->num];

                        const double ddx = pj[0] - pi[0];
                        const double ddy = pj[1] - pi[1];
                        const double ddz = pj[2] - pi[2];
                        const double norm_sq = ddx * ddx + ddy * ddy + ddz * ddz;

                        if (norm_sq > mtol * mtol && norm_sq < support_sq) {
                            const double norm = std::sqrt(norm_sq);
                            const double smooth_weight = 1.0e0 - norm_sq / support_sq;
                            const double inv_norm = 1.0e0 / norm;
                            sum[0] += ddx * inv_norm * smooth_weight;
                            sum[1] += ddy * inv_norm * smooth_weight;
                            sum[2] += ddz * inv_norm * smooth_weight;
                        }

                        jp = next_particle[jp];
                    }
                }
            }
        }

        const double scale = -dt * gamma_s * support;
        disp_corr[ip][0] = scale * sum[0];
        disp_corr[ip][1] = scale * sum[1];
        disp_corr[ip][2] = scale * sum[2];
    }

    return disp_corr;
}
