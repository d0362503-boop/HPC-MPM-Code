#include "module/mpi_data.h"
#include "module/DLB/mpm_dlb.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <vector>

void MaterialPoint::DetermineParticleRank(const std::vector<mpm_dlb::Region> &regions) {

    if (int(regions.size()) != nprocs) {
        std::cerr << "particle migration error on rank " << myrank << ": invalid region count\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<int> checkmp(this->num, -1);
    this->par_comm_.ofp_id.assign(this->num, 0);
    this->par_comm_.comm_ranks = naid;
    this->par_comm_.nsp.assign(isb, 0);
    this->par_comm_.nrp.assign(isb, 0);
    this->par_comm_.nnmp = 0;
    this->par_comm_.nrmp = 0;

    for (int ip = 0; ip < this->num; ++ip) {
        std::array<int, 3> element{};
        if (!LocateGlobalElement(this->coord[ip], element)) {
            this->par_comm_.ofp_id[this->par_comm_.nrmp++] = ip;
            continue;
        }

        bool stay = true;
        for (int dir = 0; dir < 3; ++dir) {
            if (element[dir] < aelemmin[dir] || element[dir] > aelemmax[dir]) {
                stay = false;
                break;
            }
        }
        if (stay) {
            ++this->par_comm_.nnmp;
            checkmp[ip] = -2;
            continue;
        }

        int peer = -1;
        for (int i = 0; i < isb; ++i) {
            const int rank = naid[i];
            if (rank < 0 || rank >= int(regions.size())) {
                std::cerr << "particle migration error on rank " << myrank << ": invalid peer rank\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            bool inside = true;
            for (int dir = 0; dir < 3; ++dir) {
                if (element[dir] < regions[rank].elem_min[dir] || element[dir] > regions[rank].elem_max[dir]) {
                    inside = false;
                    break;
                }
            }
            if (!inside) { continue; }

            if (peer != -1) {
                std::cerr << "particle migration error on rank " << myrank << ": overlapping peer regions\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            peer = i;
        }

        if (peer == -1) {
            std::cerr << "particle migration error on rank " << myrank << ": particle crossed a non-neighbor region\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        ++this->par_comm_.nsp[peer];
        checkmp[ip] = peer;
    }

    std::vector<int> idmploc(isb);
    this->par_comm_.nmps = 0;
    for (int i = 0; i < isb; i++) {
        idmploc[i] = this->par_comm_.nmps;
        this->par_comm_.nmps += this->par_comm_.nsp[i];
    }

    VectorAssign(this->par_comm_.nmps, this->par_comm_.idmp);
    VectorAssign(this->par_comm_.nnmp, this->par_comm_.idnsp);

    int nnmp2 = 0;
    for (int ip = 0; ip < this->num; ip++) {
        int i = checkmp[ip];
        if (i == -2) { // Stay
            this->par_comm_.idnsp[nnmp2++] = ip;
        } else if (i == -1) { // Remove
            continue;
        } else { // Move
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

    return;
}

void MaterialPoint::DetermineDLBParticleRank(const std::vector<mpm_dlb::Region> &regions) {

    if (int(regions.size()) != nprocs) {
        if (myrank == 0) { std::cerr << "DLB error: region count does not match the MPI rank count\n"; }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<int> destination(this->num, -1);
    std::vector<int> send_counts(nprocs, 0);
    this->par_comm_.ofp_id.assign(this->num, 0);
    this->par_comm_.nnmp = 0;
    this->par_comm_.nrmp = 0;

    for (int ip = 0; ip < this->num; ++ip) {
        std::array<int, 3> element{};
        if (!LocateGlobalElement(this->coord[ip], element)) {
            this->par_comm_.ofp_id[this->par_comm_.nrmp++] = ip;
            continue;
        }

        int owner = -1;
        for (int rank = 0; rank < nprocs; ++rank) {
            bool inside = true;
            for (int dir = 0; dir < 3; ++dir) {
                if (element[dir] < regions[rank].elem_min[dir] || element[dir] > regions[rank].elem_max[dir]) {
                    inside = false;
                    break;
                }
            }
            if (!inside) { continue; }

            if (owner != -1) {
                if (myrank == 0) { std::cerr << "DLB error: overlapping process regions\n"; }
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            owner = rank;
        }

        if (owner == -1) {
            if (myrank == 0) { std::cerr << "DLB error: process regions do not cover the global mesh\n"; }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        destination[ip] = owner;
        if (owner == myrank) {
            ++this->par_comm_.nnmp;
        } else {
            ++send_counts[owner];
        }
    }

    std::vector<int> receive_counts(nprocs, 0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, //
                 receive_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    this->par_comm_.comm_ranks.clear();
    this->par_comm_.nsp.clear();
    this->par_comm_.nrp.clear();
    this->par_comm_.comm_ranks.reserve(nprocs - 1);
    this->par_comm_.nsp.reserve(nprocs - 1);
    this->par_comm_.nrp.reserve(nprocs - 1);
    std::vector<int> rank_to_peer(nprocs, -1);
    for (int rank = 0; rank < nprocs; ++rank) {
        if (rank == myrank || (send_counts[rank] == 0 && receive_counts[rank] == 0)) { continue; }

        rank_to_peer[rank] = int(this->par_comm_.comm_ranks.size());
        this->par_comm_.comm_ranks.push_back(rank);
        this->par_comm_.nsp.push_back(send_counts[rank]);
        this->par_comm_.nrp.push_back(receive_counts[rank]);
    }

    this->par_comm_.nmps = 0;
    this->par_comm_.nrps = 0;
    for (int peer = 0; peer < int(this->par_comm_.comm_ranks.size()); ++peer) {
        this->par_comm_.nmps += this->par_comm_.nsp[peer];
        this->par_comm_.nrps += this->par_comm_.nrp[peer];
    }

    VectorAssign(this->par_comm_.nmps, this->par_comm_.idmp);
    VectorAssign(this->par_comm_.nnmp, this->par_comm_.idnsp);

    std::vector<int> send_offset(this->par_comm_.comm_ranks.size(), 0);
    int offset = 0;
    for (int peer = 0; peer < int(this->par_comm_.comm_ranks.size()); ++peer) {
        send_offset[peer] = offset;
        offset += this->par_comm_.nsp[peer];
    }

    int local_index = 0;
    for (int ip = 0; ip < this->num; ++ip) {
        if (destination[ip] == myrank) {
            this->par_comm_.idnsp[local_index++] = ip;
        } else if (destination[ip] != -1) {
            const int peer = rank_to_peer[destination[ip]];
            this->par_comm_.idmp[send_offset[peer]++] = ip;
        }
    }

    return;
}

void MaterialPoint::RebalanceDLBParticles(const std::vector<mpm_dlb::Region> &regions) {

    this->DetermineDLBParticleRank(regions);

    if (this->par_comm_.nmps + this->par_comm_.nrps + this->par_comm_.nrmp == 0) { return; }

    this->num += this->par_comm_.nrps - this->par_comm_.nmps - this->par_comm_.nrmp;

    this->MigrateParticleData();

    return;
}

void MaterialPoint::ApplyDLB() {

    mpm_dlb::OuputDLBParticleRatio(this->num);

    this->MoveParticle();

    const auto local_samples = mpm_dlb::SelectSamples(this->coord);

    const auto regions = mpm_dlb::ComputeDLBRegions(local_samples);

    this->RebalanceDLBParticles(regions);

    mpm_dlb::UpdateDLBRegions(regions);

    BuildMesh();

    BuildControlPoint();

    ComputeNodalVol();

    this->RebuildBC();

    return;
}
