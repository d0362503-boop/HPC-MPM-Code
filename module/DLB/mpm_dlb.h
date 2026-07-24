#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace mpm_dlb {

/**
 * @brief Inclusive background-element bounds of one MPI rank.
 */
struct Region {
    std::array<int, 3> elem_min;
    std::array<int, 3> elem_max;
};

/**
 * @brief Collect and store the current MPI element regions from every rank.
 */
void CollectCurrentRegions();

/**
 * @brief Return the current MPI element regions ordered by rank.
 * @return Current MPI element regions.
 */
const std::vector<Region> &CurrentRegions();

/**
 * @brief Compute the uniform local particle sampling interval for DLB.
 * @param local_npts Number of material points currently stored on this MPI rank.
 * @return Positive index interval used to sample local material-point coordinates.
 */
int ComputeSampleSkip(std::size_t local_npts);

/**
 * @brief Uniformly sample local material-point coordinates using a fixed interval.
 * @param coord Coordinates of material points stored on this MPI rank.
 * @param skip Positive particle-index interval between consecutive samples.
 * @return Local coordinate samples that represent the particle density on this rank.
 */
std::vector<std::array<double, 3>> SelectSamples(const std::vector<std::array<double, 3>> &coord, int skip);

/**
 * @brief Gather local DLB coordinate samples to MPI rank zero.
 * @param local_samples Coordinate samples selected on this MPI rank.
 * @param global_samples Output coordinate samples on rank zero; empty on other ranks.
 */
void GatherSamplesToRoot(const std::vector<std::array<double, 3>> &local_samples,
                         std::vector<std::array<double, 3>> &global_samples);

/**
 * @brief Compute integer background-element cuts from coordinate samples in one direction.
 * @param global_samples Samples belonging to the box being subdivided.
 * @param dir Spatial direction: 0 for x, 1 for y, and 2 for z.
 * @param lo Inclusive first background-element index of the parent box.
 * @param hi Inclusive last background-element index of the parent box.
 * @param parts Number of child boxes required in the selected direction.
 * @param prob_lo Physical coordinate of global background element index zero.
 * @param cell_size Background-element size in the selected direction.
 * @return Half-open element cuts with size `parts + 1`.
 */
std::vector<int> CutsFromSamples(const std::vector<std::array<double, 3>> &global_samples, int dir, int lo, int hi,
                                 int parts, double prob_lo, double cell_size);

/**
 * @brief Compute rebalanced fixed-topology MPI regions from local particle samples.
 * @param local_samples Coordinate samples selected on this MPI rank.
 * @return One inclusive background-element region per MPI rank, ordered by rank.
 */
std::vector<Region> ComputeDLBRegions(const std::vector<std::array<double, 3>> &local_samples);

/**
 * @brief Update the MPI-region layout and rebuild its overlap communication metadata.
 * @param regions One inclusive background-element region per MPI rank, ordered by rank.
 *
 * Updates the current rank's element bounds, local mesh dimensions, and node/control-point overlap
 * tables.  The caller remains responsible for rebuilding mesh storage and physics-specific BC data.
 */
void UpdateDLBRegions(const std::vector<Region> &regions);

}  // namespace mpm_dlb
