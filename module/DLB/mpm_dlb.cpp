#include "module/DLB/mpm_dlb.h"

#include "module/mesh.h"
#include "module/mpi_data.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <vector>

namespace mpm_dlb {

constexpr int kRootRank = 0;
constexpr double kMinSampleRate = 1.0e-4;
constexpr double kSamplesPerRank = 256.0e0;

void OuputDLBParticleRatio(int num) {

    int max_particle, min_particle;
    MPI_Allreduce(&num, &max_particle, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&num, &min_particle, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    double ratio = (min_particle > 0) ? double(max_particle) / double(min_particle) : 0.0e0;

    if (myrank == kRootRank) {
        std::cout << "Particle Max = " << std::setw(10) << max_particle     //
                  << "    Particle Min = " << std::setw(10) << min_particle //
                  << "    Particle Max / Min ratio = " << ratio << "\n";
    }

    return;
}

void AbortDLB(const char *message) {
    if (myrank == kRootRank) { std::cerr << "DLB error: " << message << "\n"; }
    MPI_Abort(MPI_COMM_WORLD, 1);
}

void CollectCurrentRegions(const Mesh &mesh, Partition &part) {

    Region local_region;
    for (int dir = 0; dir < 3; ++dir) {
        local_region.elem_min[dir] = mesh.aelemmin[dir];
        local_region.elem_max[dir] = mesh.aelemmax[dir];
    }

    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized == 0) {
        part.regions.assign(1, local_region);
        return;
    }

    int mpi_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    std::array<int, 6> local_bounds = {mesh.aelemmin[0], mesh.aelemmin[1], mesh.aelemmin[2],
                                       mesh.aelemmax[0], mesh.aelemmax[1], mesh.aelemmax[2]};
    std::vector<int> gathered_bounds(6 * mpi_size);
    MPI_Allgather(local_bounds.data(), 6, MPI_INT, gathered_bounds.data(), 6, MPI_INT, MPI_COMM_WORLD);

    part.regions.resize(mpi_size);
    for (int rank = 0; rank < mpi_size; ++rank) {
        for (int dir = 0; dir < 3; ++dir) {
            part.regions[rank].elem_min[dir] = gathered_bounds[6 * rank + dir];
            part.regions[rank].elem_max[dir] = gathered_bounds[6 * rank + 3 + dir];
        }
    }
}


// Keep samples inside one recursive child range.
std::vector<std::array<double, 3>> FilterSamples(const std::vector<std::array<double, 3>> &samples, int dir, int lo,
                                                 int hi) {
    std::vector<std::array<double, 3>> filtered;
    filtered.reserve(samples.size());

    for (const auto &sample : samples) {
        std::array<int, 3> element{};
        if (!LocateGlobalElement(g_grid.mesh, sample, element)) { AbortDLB("DLB sample lies outside the global background mesh"); }
        if (element[dir] >= lo && element[dir] <= hi) { filtered.push_back(sample); }
    }

    return filtered;
}

int ComputeSampleSkip(std::size_t local_npts) {

    unsigned long long local_count = static_cast<unsigned long long>(local_npts);
    unsigned long long global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPIDatatypeCheck<unsigned long long>::GetType(), MPI_SUM,
                  MPI_COMM_WORLD);

    if (global_count == 0) { return 1; }

    const double sample_rate = std::max(kMinSampleRate, kSamplesPerRank * double(nprocs) / double(global_count));
    const int target_samples = std::max(1, int(std::ceil(double(local_count) * sample_rate)));

    if (target_samples >= int(local_count)) { return 1; }
    return std::max(1, int(local_count) / target_samples);
}

std::vector<std::array<double, 3>> SelectSamples(const std::vector<std::array<double, 3>> &coord) {

    const int skip = ComputeSampleSkip(coord.size());
    if (skip <= 0) { AbortDLB("sample skip must be positive"); }

    std::vector<std::array<double, 3>> samples;
    samples.reserve((coord.size() + std::size_t(skip) - 1) / std::size_t(skip));

    for (std::size_t point = 0; point < coord.size(); point += std::size_t(skip)) { samples.push_back(coord[point]); }

    return samples;
}

void GatherSamplesToRoot(const std::vector<std::array<double, 3>> &local_samples,
                         std::vector<std::array<double, 3>> &global_samples) {

    const int local_count = int(local_samples.size());
    std::vector<int> sample_counts;

    if (myrank == kRootRank) { sample_counts.resize(nprocs, 0); }
    MPI_Gather(&local_count, 1, MPI_INT, //
               myrank == kRootRank ? sample_counts.data() : nullptr, 1, MPI_INT, kRootRank, MPI_COMM_WORLD);

    std::vector<double> send_buffer;
    send_buffer.reserve(local_samples.size() * 3);
    for (const auto &sample : local_samples) {
        send_buffer.push_back(sample[0]);
        send_buffer.push_back(sample[1]);
        send_buffer.push_back(sample[2]);
    }

    std::vector<int> receive_counts;
    std::vector<int> displacements;
    std::vector<double> receive_buffer;
    if (myrank == kRootRank) {
        receive_counts.resize(nprocs, 0);
        displacements.resize(nprocs, 0);

        int global_count = 0;
        for (int rank = 0; rank < nprocs; ++rank) {
            receive_counts[rank] = sample_counts[rank] * 3;
            displacements[rank] = global_count;
            global_count += receive_counts[rank];
        }

        receive_buffer.resize(global_count);
        global_samples.resize(global_count / 3);
    } else {
        global_samples.clear();
    }

    MPI_Gatherv(send_buffer.empty() ? nullptr : send_buffer.data(), int(send_buffer.size()),
                MPIDatatypeCheck<double>::GetType(), receive_buffer.empty() ? nullptr : receive_buffer.data(),
                receive_counts.empty() ? nullptr : receive_counts.data(),
                displacements.empty() ? nullptr : displacements.data(), MPIDatatypeCheck<double>::GetType(), kRootRank,
                MPI_COMM_WORLD);

    if (myrank == kRootRank) {
        for (std::size_t sample = 0; sample < global_samples.size(); ++sample) {
            global_samples[sample] = {receive_buffer[3 * sample], receive_buffer[3 * sample + 1],
                                      receive_buffer[3 * sample + 2]};
        }
    }
}

std::vector<int> CutsFromSamples(const std::vector<std::array<double, 3>> &global_samples, int dir, int lo, int hi,
                                 int parts, double prob_lo, double cell_size) {
    if (dir < 0 || dir > 2 || parts < 1 || hi < lo || cell_size <= 0.0 || hi - lo + 1 < parts) {
        AbortDLB("invalid direction, element range, cell size, or process count for DLB cut");
    }

    std::vector<int> cuts(parts + 1);
    cuts.front() = lo;
    cuts.back() = hi + 1;

    std::vector<double> coordinates;
    coordinates.reserve(global_samples.size());
    for (const auto &sample : global_samples) { coordinates.push_back(sample[dir]); }
    std::sort(coordinates.begin(), coordinates.end());

    for (int part = 1; part < parts; ++part) {
        const int min_cut = cuts[part - 1] + 1;
        const int max_cut = hi + 1 - (parts - part);
        int cut = lo + (hi + 1 - lo) * part / parts;

        if (!coordinates.empty()) {
            const int sample_index = std::clamp(int(coordinates.size()) * part / parts, 0, int(coordinates.size()) - 1);
            double cut_coordinate = coordinates[sample_index];
            if (sample_index + 1 < int(coordinates.size())) {
                cut_coordinate = 0.5e0 * (coordinates[sample_index] + coordinates[sample_index + 1]);
            }
            cut = int(std::round((cut_coordinate - prob_lo) / cell_size));
        }

        cuts[part] = std::clamp(cut, min_cut, max_cut);
    }

    return cuts;
}

std::vector<Region> ComputeDLBRegions(const Grid &grid, const std::vector<std::array<double, 3>> &local_samples) {
    if (grid.part.nxyr[0] * grid.part.nxyr[1] * grid.part.nxyr[2] != nprocs) { AbortDLB("nxyr does not match the MPI rank count"); }

    for (int dir = 0; dir < 3; ++dir) {
        if (grid.mesh.xyelemw[dir] < grid.part.nxyr[dir] || grid.mesh.dxy[dir] <= 0.0) {
            AbortDLB("the global mesh cannot provide one element for every process split");
        }
    }

    std::vector<std::array<double, 3>> global_samples;
    GatherSamplesToRoot(local_samples, global_samples);

    const int nx = grid.part.nxyr[0];
    const int ny = grid.part.nxyr[1];
    const int nz = grid.part.nxyr[2];
    std::vector<Region> regions(nprocs);

    if (myrank == kRootRank) {
        const auto z_cuts = CutsFromSamples(global_samples, 2, 0, grid.mesh.xyelemw[2] - 1, nz, grid.mesh.xyminw[2], grid.mesh.dxy[2]);

        for (int iz = 0; iz < nz; ++iz) {
            const auto z_samples = FilterSamples(global_samples, 2, z_cuts[iz], z_cuts[iz + 1] - 1);
            const auto y_cuts = CutsFromSamples(z_samples, 1, 0, grid.mesh.xyelemw[1] - 1, ny, grid.mesh.xyminw[1], grid.mesh.dxy[1]);

            for (int iy = 0; iy < ny; ++iy) {
                const auto y_samples = FilterSamples(z_samples, 1, y_cuts[iy], y_cuts[iy + 1] - 1);
                const auto x_cuts = CutsFromSamples(y_samples, 0, 0, grid.mesh.xyelemw[0] - 1, nx, grid.mesh.xyminw[0], grid.mesh.dxy[0]);

                for (int ix = 0; ix < nx; ++ix) {
                    const int rank = ix + nx * iy + nx * ny * iz;
                    regions[rank].elem_min = {x_cuts[ix], y_cuts[iy], z_cuts[iz]};
                    regions[rank].elem_max = {x_cuts[ix + 1] - 1, y_cuts[iy + 1] - 1, z_cuts[iz + 1] - 1};
                }
            }
        }
    }

    std::vector<int> packed_regions(6 * nprocs, 0);
    if (myrank == kRootRank) {
        for (int rank = 0; rank < nprocs; ++rank) {
            for (int dir = 0; dir < 3; ++dir) {
                packed_regions[6 * rank + dir] = regions[rank].elem_min[dir];
                packed_regions[6 * rank + 3 + dir] = regions[rank].elem_max[dir];
            }
        }
    }

    MPI_Bcast(packed_regions.data(), int(packed_regions.size()), MPI_INT, kRootRank, MPI_COMM_WORLD);

    if (myrank != kRootRank) {
        for (int rank = 0; rank < nprocs; ++rank) {
            for (int dir = 0; dir < 3; ++dir) {
                regions[rank].elem_min[dir] = packed_regions[6 * rank + dir];
                regions[rank].elem_max[dir] = packed_regions[6 * rank + 3 + dir];
            }
        }
    }

    return regions;
}

void UpdateDLBRegions(const std::vector<Region> &regions, Grid &grid) {
    if (int(regions.size()) != nprocs || myrank < 0 || myrank >= nprocs) {
        AbortDLB("region count or MPI rank is invalid while updating DLB regions");
    }

    for (const Region &region : regions) {
        for (int dir = 0; dir < 3; ++dir) {
            if (region.elem_min[dir] < 0 || region.elem_max[dir] < region.elem_min[dir] ||
                region.elem_max[dir] >= grid.mesh.xyelemw[dir]) {
                AbortDLB("DLB region lies outside the global background mesh");
            }
        }
    }

    grid.part.regions = regions;

    const Region &local_region = regions[myrank];
    for (int dir = 0; dir < 3; ++dir) {
        grid.mesh.aelemmin[dir] = local_region.elem_min[dir];
        grid.mesh.aelemmax[dir] = local_region.elem_max[dir];
        grid.mesh.xyelem[dir] = grid.mesh.aelemmax[dir] - grid.mesh.aelemmin[dir] + 1;
        grid.mesh.xynode[dir] = grid.mesh.xyelem[dir] + 1;
        grid.cmesh.xynodec[dir] = grid.mesh.xyelem[dir] + grid.cmesh.idimc[dir];
        grid.mesh.xymin[dir] = grid.mesh.xyminw[dir] + grid.mesh.dxy[dir] * double(grid.mesh.aelemmin[dir]);
        grid.mesh.xymax[dir] = grid.mesh.xyminw[dir] + grid.mesh.dxy[dir] * double(grid.mesh.aelemmax[dir] + 1);
    }

    grid.mesh.nelem = grid.mesh.xyelem[0] * grid.mesh.xyelem[1] * grid.mesh.xyelem[2];
    grid.mesh.node = grid.mesh.xynode[0] * grid.mesh.xynode[1] * grid.mesh.xynode[2];
    grid.cmesh.nodec = grid.cmesh.xynodec[0] * grid.cmesh.xynodec[1] * grid.cmesh.xynodec[2];
    grid.mesh.nu = 0;
    grid.mesh.nv = grid.mesh.nu + grid.mesh.node;
    grid.mesh.nw = grid.mesh.nv + grid.mesh.node;
    grid.mesh.np = grid.mesh.nw + grid.mesh.node;
    grid.cmesh.nuc = 0;
    grid.cmesh.nvc = grid.cmesh.nuc + grid.cmesh.nodec;
    grid.cmesh.nwc = grid.cmesh.nvc + grid.cmesh.nodec;
    grid.cmesh.npc = grid.cmesh.nwc + grid.cmesh.nodec;

    grid.part.naid.clear();
    grid.part.nsbc.clear();
    grid.part.nsbl.clear();
    grid.part.nsubc.clear();
    grid.part.nsubl.clear();
    grid.part.naid.reserve(nprocs - 1);
    grid.part.nsbc.reserve(nprocs - 1);
    grid.part.nsbl.reserve(nprocs - 1);

    std::vector<int> cp_shares(grid.cmesh.nodec, 1);
    std::vector<int> node_shares(grid.mesh.node, 1);

    for (int rank = 0; rank < nprocs; ++rank) {
        if (rank == myrank) { continue; }

        std::array<int, 3> cp_min{};
        std::array<int, 3> cp_max{};
        std::array<int, 3> node_min{};
        std::array<int, 3> node_max{};
        bool cp_overlap = true;
        bool node_overlap = true;
        for (int dir = 0; dir < 3; ++dir) {
            cp_min[dir] = std::max(grid.mesh.aelemmin[dir], regions[rank].elem_min[dir]);
            cp_max[dir] = std::min(grid.mesh.aelemmax[dir] + grid.cmesh.idimc[dir], regions[rank].elem_max[dir] + grid.cmesh.idimc[dir]);
            if (cp_min[dir] > cp_max[dir]) { cp_overlap = false; }

            node_min[dir] = std::max(grid.mesh.aelemmin[dir], regions[rank].elem_min[dir]);
            node_max[dir] = std::min(grid.mesh.aelemmax[dir] + 1, regions[rank].elem_max[dir] + 1);
            if (node_min[dir] > node_max[dir]) { node_overlap = false; }
        }
        if (!cp_overlap && !node_overlap) { continue; }

        int cp_count = 0;
        if (cp_overlap) {
            cp_count = (cp_max[0] - cp_min[0] + 1) * (cp_max[1] - cp_min[1] + 1) * (cp_max[2] - cp_min[2] + 1);
        }

        int node_count = 0;
        if (node_overlap) {
            node_count =
                (node_max[0] - node_min[0] + 1) * (node_max[1] - node_min[1] + 1) * (node_max[2] - node_min[2] + 1);
        }

        grid.part.naid.push_back(rank);
        grid.part.nsbc.push_back(cp_count);
        grid.part.nsbl.push_back(node_count);

        if (node_overlap) {
            for (int k = node_min[2]; k <= node_max[2]; ++k) {
                for (int j = node_min[1]; j <= node_max[1]; ++j) {
                    for (int i = node_min[0]; i <= node_max[0]; ++i) {
                        const int il = i - grid.mesh.aelemmin[0];
                        const int jl = j - grid.mesh.aelemmin[1];
                        const int kl = k - grid.mesh.aelemmin[2];
                        const int local_id = il + grid.mesh.xynode[0] * jl + grid.mesh.xynode[0] * grid.mesh.xynode[1] * kl;
                        grid.part.nsubl.push_back(local_id);
                        ++node_shares[local_id];
                    }
                }
            }
        }

        if (cp_overlap) {
            for (int k = cp_min[2]; k <= cp_max[2]; ++k) {
                for (int j = cp_min[1]; j <= cp_max[1]; ++j) {
                    for (int i = cp_min[0]; i <= cp_max[0]; ++i) {
                        const int il = i - grid.mesh.aelemmin[0];
                        const int jl = j - grid.mesh.aelemmin[1];
                        const int kl = k - grid.mesh.aelemmin[2];
                        const int local_id = il + grid.cmesh.xynodec[0] * jl + grid.cmesh.xynodec[0] * grid.cmesh.xynodec[1] * kl;
                        grid.part.nsubc.push_back(local_id);
                        ++cp_shares[local_id];
                    }
                }
            }
        }
    }

    grid.part.isb = int(grid.part.naid.size());
    grid.part.isubc = int(grid.part.nsubc.size());
    grid.part.isubl = int(grid.part.nsubl.size());
    grid.part.dbc.assign(grid.cmesh.nodec * 4, 1.0e0);
    grid.part.dbl.assign(grid.mesh.node * 4, 1.0e0);
    for (int nid = 0; nid < grid.cmesh.nodec; ++nid) {
        const double weight = 1.0e0 / double(cp_shares[nid]);
        for (int dof = 0; dof < 4; ++dof) { grid.part.dbc[nid + dof * grid.cmesh.nodec] = weight; }
    }
    for (int nid = 0; nid < grid.mesh.node; ++nid) {
        const double weight = 1.0e0 / double(node_shares[nid]);
        for (int dof = 0; dof < 4; ++dof) { grid.part.dbl[nid + dof * grid.mesh.node] = weight; }
    }
}

} // namespace mpm_dlb
