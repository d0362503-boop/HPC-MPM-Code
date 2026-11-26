#include "module/bc.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

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

void BoundaryCondition::BCSetVal(int nn, std::vector<double> &variable, double scale) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = this->fbc[i] * scale;
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

void BoundaryCondition::BCSetDispl(int nn, std::vector<double> &variable, double scale) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = scale * this->fbc[i];
        }
    }

    return;
}

void BoundaryCondition::BCSetFreeSlip(std::vector<double> &variable, const std::vector<double> &normal,
                                      const ControlMesh &cmesh, double tol) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            const int n1 = this->nbc[i] + cmesh.nuc;
            const int n2 = this->nbc[i] + cmesh.nvc;
            const int n3 = this->nbc[i] + cmesh.nwc;
            const double nx = normal[n1];
            const double ny = normal[n2];
            const double nz = normal[n3];
            const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > tol) {
                const double qn = variable[n1] * nx + variable[n2] * ny + variable[n3] * nz;
                variable[n1] -= qn * nx;
                variable[n2] -= qn * ny;
                variable[n3] -= qn * nz;
            } else {
                variable[n1] = 0.0e0;
                variable[n2] = 0.0e0;
                variable[n3] = 0.0e0;
            }
        }
    }

    return;
}

void BoundaryCondition::CaptureGlobalControlPointBC(const Grid &grid) {
    const Mesh &mesh = grid.mesh;
    const ControlMesh &cmesh = grid.cmesh;

    std::vector<int> local_global_ids(this->ibc);
    for (int m = 0; m < this->ibc; ++m) {
        const int nid = this->nbc[m];
        const int iz = nid / (cmesh.xynodec[0] * cmesh.xynodec[1]);
        const int iy = (nid - iz * cmesh.xynodec[0] * cmesh.xynodec[1]) / cmesh.xynodec[0];
        const int ix = nid - iz * cmesh.xynodec[0] * cmesh.xynodec[1] - iy * cmesh.xynodec[0];

        const int gx = ix + mesh.aelemmin[0];
        const int gy = iy + mesh.aelemmin[1];
        const int gz = iz + mesh.aelemmin[2];
        local_global_ids[m] = gx + cmesh.xynodecw[0] * gy + cmesh.xynodecw[0] * cmesh.xynodecw[1] * gz;
    }

    this->CacheGlobalEntries(local_global_ids, true);

    return;
}

void BoundaryCondition::RebuildLocalControlPointBC(const Grid &grid) {
    const Mesh &mesh = grid.mesh;
    const ControlMesh &cmesh = grid.cmesh;

    this->nbc.clear();
    this->fbc.clear();
    this->nbc.reserve(this->global_nbc.size());
    this->fbc.reserve(this->global_fbc.size());

    for (int m = 0; m < int(this->global_nbc.size()); ++m) {
        const int gid = this->global_nbc[m];
        const int gz = gid / (cmesh.xynodecw[0] * cmesh.xynodecw[1]);
        const int gy = (gid - gz * cmesh.xynodecw[0] * cmesh.xynodecw[1]) / cmesh.xynodecw[0];
        const int gx = gid - gz * cmesh.xynodecw[0] * cmesh.xynodecw[1] - gy * cmesh.xynodecw[0];

        if (gx < mesh.aelemmin[0] || gx > mesh.aelemmax[0] + cmesh.idimc[0] || //
            gy < mesh.aelemmin[1] || gy > mesh.aelemmax[1] + cmesh.idimc[1] || //
            gz < mesh.aelemmin[2] || gz > mesh.aelemmax[2] + cmesh.idimc[2]) {
            continue;
        }

        const int ix = gx - mesh.aelemmin[0];
        const int iy = gy - mesh.aelemmin[1];
        const int iz = gz - mesh.aelemmin[2];
        this->nbc.push_back(ix + cmesh.xynodec[0] * iy + cmesh.xynodec[0] * cmesh.xynodec[1] * iz);
        this->fbc.push_back(this->global_fbc[m]);
    }
    this->ibc = int(this->nbc.size());

    return;
}

void BoundaryCondition::CaptureGlobalElementBC(const Mesh &mesh) {

    std::vector<int> local_global_ids(this->ibc);
    for (int m = 0; m < this->ibc; ++m) {
        const int eid = this->nbc[m];
        const int iz = eid / (mesh.xyelem[0] * mesh.xyelem[1]);
        const int iy = (eid - iz * mesh.xyelem[0] * mesh.xyelem[1]) / mesh.xyelem[0];
        const int ix = eid - iz * mesh.xyelem[0] * mesh.xyelem[1] - iy * mesh.xyelem[0];

        const int gx = ix + mesh.aelemmin[0];
        const int gy = iy + mesh.aelemmin[1];
        const int gz = iz + mesh.aelemmin[2];
        local_global_ids[m] = gx + mesh.xyelemw[0] * gy + mesh.xyelemw[0] * mesh.xyelemw[1] * gz;
    }

    this->CacheGlobalEntries(local_global_ids, false);

    return;
}

void BoundaryCondition::RebuildLocalElementBC(const Mesh &mesh) {

    this->nbc.clear();
    this->nbc.reserve(this->global_nbc.size());
    for (int gid : this->global_nbc) {
        const int gz = gid / (mesh.xyelemw[0] * mesh.xyelemw[1]);
        const int gy = (gid - gz * mesh.xyelemw[0] * mesh.xyelemw[1]) / mesh.xyelemw[0];
        const int gx = gid - gz * mesh.xyelemw[0] * mesh.xyelemw[1] - gy * mesh.xyelemw[0];

        if (gx < mesh.aelemmin[0] || gx > mesh.aelemmax[0] || //
            gy < mesh.aelemmin[1] || gy > mesh.aelemmax[1] || //
            gz < mesh.aelemmin[2] || gz > mesh.aelemmax[2]) {
            continue;
        }

        const int ix = gx - mesh.aelemmin[0];
        const int iy = gy - mesh.aelemmin[1];
        const int iz = gz - mesh.aelemmin[2];
        this->nbc.push_back(ix + mesh.xyelem[0] * iy + mesh.xyelem[0] * mesh.xyelem[1] * iz);
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
                       all_values.data(), counts.data(), offsets.data(), MPIDatatypeCheck<double>::GetType(), MPI_COMM_WORLD);
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
