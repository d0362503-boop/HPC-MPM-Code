#include "bc.h"
#include "dataset.h"
#include "mesh.h"
#include "mpi_data.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

void BoundaryCondition::BCInput(std::ifstream &infile, bool if_fbc) {

    infile >> this->ibc;
    infile.ignore(1000, '\n');
    this->nbc.resize(this->ibc);
    if (if_fbc) this->fbc.resize(this->ibc);

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            if (if_fbc) {
                infile >> this->nbc[i] >> this->fbc[i];
            } else {
                infile >> this->nbc[i];
            }
            infile.ignore(1000, '\n');
        }
    }

    return;
}

void BoundaryCondition::BCOutput(std::ofstream &outfile, std::string bc_name, bool if_fbc) {

    outfile << std::setw(15) << this->ibc << std::setw(15) << "! ---- " << bc_name << ".ibc ----" << "\n";

    for (int i = 0; i < this->ibc; i++) {
        if (if_fbc) {
            outfile << std::setw(15) << this->nbc[i] << std::setw(15) << this->fbc[i] << "\n";
        } else {
            outfile << std::setw(15) << this->nbc[i] << "\n";
        }
    }

    return;
}

void BoundaryCondition::BCSetVal(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = this->fbc[i] * facl;
        }
    }

    return;
}

void BoundaryCondition::BCSetZero(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = 0.0e0;
        }
    }

    return;
}

void BoundaryCondition::BCSetDt(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = this->fbc[i] * dt * facl;
        }
    }

    return;
}

void BoundaryCondition::CaptureGlobalControlPointBC() {

    std::vector<int> local_global_ids(this->ibc);
    for (int m = 0; m < this->ibc; ++m) {
        const int nid = this->nbc[m];
        const int iz = nid / (xynodec[0] * xynodec[1]);
        const int iy = (nid - iz * xynodec[0] * xynodec[1]) / xynodec[0];
        const int ix = nid - iz * xynodec[0] * xynodec[1] - iy * xynodec[0];

        const int gx = ix + aelemmin[0];
        const int gy = iy + aelemmin[1];
        const int gz = iz + aelemmin[2];
        local_global_ids[m] = gx + xynodecw[0] * gy + xynodecw[0] * xynodecw[1] * gz;
    }

    this->CacheGlobalEntries(local_global_ids, true);

    return;
}

void BoundaryCondition::RebuildLocalControlPointBC() {

    this->nbc.clear();
    this->fbc.clear();
    this->nbc.reserve(this->global_nbc.size());
    this->fbc.reserve(this->global_fbc.size());

    for (int m = 0; m < int(this->global_nbc.size()); ++m) {
        const int gid = this->global_nbc[m];
        const int gz = gid / (xynodecw[0] * xynodecw[1]);
        const int gy = (gid - gz * xynodecw[0] * xynodecw[1]) / xynodecw[0];
        const int gx = gid - gz * xynodecw[0] * xynodecw[1] - gy * xynodecw[0];

        if (gx < aelemmin[0] || gx > aelemmax[0] + idimc[0] || //
            gy < aelemmin[1] || gy > aelemmax[1] + idimc[1] || //
            gz < aelemmin[2] || gz > aelemmax[2] + idimc[2]) {
            continue;
        }

        const int ix = gx - aelemmin[0];
        const int iy = gy - aelemmin[1];
        const int iz = gz - aelemmin[2];
        this->nbc.push_back(ix + xynodec[0] * iy + xynodec[0] * xynodec[1] * iz);
        this->fbc.push_back(this->global_fbc[m]);
    }
    this->ibc = int(this->nbc.size());

    return;
}

void BoundaryCondition::CaptureGlobalElementBC() {

    std::vector<int> local_global_ids(this->ibc);
    for (int m = 0; m < this->ibc; ++m) {
        const int eid = this->nbc[m];
        const int iz = eid / (xyelem[0] * xyelem[1]);
        const int iy = (eid - iz * xyelem[0] * xyelem[1]) / xyelem[0];
        const int ix = eid - iz * xyelem[0] * xyelem[1] - iy * xyelem[0];

        const int gx = ix + aelemmin[0];
        const int gy = iy + aelemmin[1];
        const int gz = iz + aelemmin[2];
        local_global_ids[m] = gx + xyelemw[0] * gy + xyelemw[0] * xyelemw[1] * gz;
    }

    this->CacheGlobalEntries(local_global_ids, false);

    return;
}

void BoundaryCondition::RebuildLocalElementBC() {

    this->nbc.clear();
    this->nbc.reserve(this->global_nbc.size());
    for (int gid : this->global_nbc) {
        const int gz = gid / (xyelemw[0] * xyelemw[1]);
        const int gy = (gid - gz * xyelemw[0] * xyelemw[1]) / xyelemw[0];
        const int gx = gid - gz * xyelemw[0] * xyelemw[1] - gy * xyelemw[0];

        if (gx < aelemmin[0] || gx > aelemmax[0] || //
            gy < aelemmin[1] || gy > aelemmax[1] || //
            gz < aelemmin[2] || gz > aelemmax[2]) {
            continue;
        }

        const int ix = gx - aelemmin[0];
        const int iy = gy - aelemmin[1];
        const int iz = gz - aelemmin[2];
        this->nbc.push_back(ix + xyelem[0] * iy + xyelem[0] * xyelem[1] * iz);
    }
    this->ibc = int(this->nbc.size());

    return;
}

void BoundaryCondition::CacheGlobalEntries(const std::vector<int> &local_global_ids, bool has_values) {

    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        this->global_nbc = local_global_ids;
        this->global_fbc = has_values ? this->fbc : std::vector<double>{};
        return;
    }

    const int local_count = int(local_global_ids.size());
    std::vector<int> counts(nprocs);
    std::vector<int> offsets(nprocs);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_count = 0;
    for (int rank = 0; rank < nprocs; ++rank) {
        offsets[rank] = total_count;
        total_count += counts[rank];
    }

    std::vector<int> all_ids(total_count);
    MPI_Allgatherv(local_global_ids.data(), local_count, MPI_INT, //
                   all_ids.data(), counts.data(), offsets.data(), MPI_INT, MPI_COMM_WORLD);

    std::vector<double> all_values;
    if (has_values) {
        if (int(this->fbc.size()) != local_count) {
            std::cerr << "BC error: prescribed-value count does not match BC-node count\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        all_values.resize(total_count);
        MPI_Allgatherv(this->fbc.data(), local_count, MPIDatatypeCheck<double>::GetType(), //
                       all_values.data(), counts.data(), offsets.data(), MPIDatatypeCheck<double>::GetType(),
                       MPI_COMM_WORLD);
    }

    this->global_nbc.clear();
    this->global_fbc.clear();
    this->global_nbc.reserve(total_count);
    if (has_values) { this->global_fbc.reserve(total_count); }
    std::unordered_map<int, int> global_index;
    global_index.reserve(total_count);

    for (int m = 0; m < total_count; ++m) {
        const auto entry = global_index.emplace(all_ids[m], int(this->global_nbc.size()));
        if (!entry.second) {
            if (has_values && std::abs(this->global_fbc[entry.first->second] - all_values[m]) > 1.0e-12) {
                std::cerr << "BC error: inconsistent prescribed values on one global boundary node\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            continue;
        }

        this->global_nbc.push_back(all_ids[m]);
        if (has_values) { this->global_fbc.push_back(all_values[m]); }
    }

    return;
}
