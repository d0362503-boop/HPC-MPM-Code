#pragma once

#include "data/divide/data_partitioner.h"

#include "module/bc.h"
#include "module/material_point.h"

#include <string>

class FluidDivider : public DataPartitioner {
  public:
    /**
     * @brief Partition fluid material points into this rank's local window.
     * @param rank_id Index of the MPI rank whose partition is being built.
     *
     * Must be called before MeshPartition.
     */
    void PartitionPoints(int rank_id);

    /**
     * @brief Renumber fluid boundary conditions to this rank's local indexing.
     * @param cp_min First global control-point index of the rank window, per direction.
     * @param cp_max Last global control-point index (inclusive) of the rank window, per direction.
     *
     * Velocity/pressure BCs are control-point based; inflow BCs are element based.
     * Must be called after MeshPartition.
     */
    void PartitionBCs(const std::vector<int> &cp_min, const std::vector<int> &cp_max);

    MaterialPoint points_;            // loaded fluid particles
    MaterialPoint partition_points_;  // partitioned fluid particles

  protected:
    /** @brief Return the fluid case label. */
    std::string CaseName() const override;

    /**
     * @brief Load fluid boundary conditions.
     *
     * Reads the four control-point BC blocks (ubc/vbc/wbc/pbc) followed by the
     * three inflow BC blocks (uinfbc/vinfbc/winfbc, element IDs without values).
     */
    void LoadBoundaryData(std::ifstream &infile) override;

    /** @brief Load fluid particle records (id, material, mass, volume). */
    void LoadPointData(std::ifstream &infile) override;

    /** @brief Run point partitioning, mesh windowing, and BC renumbering for one rank. */
    void PartitionProcess(int rank_id) override;

    /** @brief Write partitioned fluid boundary conditions, inflow blocks included. */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /** @brief Write partitioned fluid particle records. */
    void WritePointData(std::ofstream &outfile) override;
};
