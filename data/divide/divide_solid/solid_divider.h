#pragma once

#include "data/divide/data_partitioner.h"

#include "module/bc.h"
#include "module/material_point.h"

#include <string>

class SolidDivider : public DataPartitioner {
  public:
    /**
     * @brief Partition solid material points into this rank's local window.
     * @param rank_id Index of the MPI rank whose partition is being built.
     *
     * Must be called before MeshPartition.
     */
    void PartitionPoints(int rank_id);

    /**
     * @brief Renumber solid boundary conditions to this rank's local control-point indexing.
     * @param cp_min First global control-point index of the rank window, per direction.
     * @param cp_max Last global control-point index (inclusive) of the rank window, per direction.
     *
     * Must be called after MeshPartition.
     */
    void PartitionBCs(const std::vector<int> &cp_min, const std::vector<int> &cp_max);

    MaterialPoint points_;            // loaded solid particles
    MaterialPoint partition_points_;  // partitioned solid particles

  protected:
    /** @brief Return the solid case label. */
    std::string CaseName() const override;

    /** @brief Load solid velocity boundary conditions (usbc/vsbc/wsbc). */
    void LoadBoundaryData(std::ifstream &infile) override;

    /** @brief Load solid particle records (id, material, surface flag, mass, volume). */
    void LoadPointData(std::ifstream &infile) override;

    /** @brief Run point partitioning, mesh windowing, and BC renumbering for one rank. */
    void PartitionProcess(int rank_id) override;

    /** @brief Write partitioned solid boundary conditions (usbc/vsbc/wsbc). */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /** @brief Write partitioned solid particle records. */
    void WritePointData(std::ofstream &outfile) override;
};
