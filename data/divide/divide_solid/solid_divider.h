#pragma once

#include "data/divide/data_partitioner.h"

#include "module/material_point.h"

#include <string>

class SolidDivider : public DataPartitioner {
  public:
    /**
     * @brief Bind the solid partitioner and its particle sets to the grid being partitioned.
     * @param grid   Grid being partitioned (global mesh, overwritten with each rank window).
     * @param params Run parameters of the case.
     */
    SolidDivider(Grid &grid, RunParams &params)
        : DataPartitioner(grid, params), points_(grid, params), partition_points_(grid, params) {}

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
    void PartitionBCs(const std::array<int, 3> &cp_min, const std::array<int, 3> &cp_max);

    MaterialPoint points_;            // loaded solid particles
    MaterialPoint partition_points_;  // partitioned solid particles

  protected:
    /** @brief Return the solid case label. */
    std::string CaseName() const override;

    /**
     * @brief Load solid velocity boundary conditions (usbc/vsbc/wsbc).
     * @param infile Input stream positioned at the global solid boundary-condition blocks.
     */
    void LoadBoundaryData(std::ifstream &infile) override;

    /**
     * @brief Load solid particle records (id, material, surface flag, mass, volume).
     * @param infile Input stream positioned at the global solid particle records.
     */
    void LoadPointData(std::ifstream &infile) override;

    /**
     * @brief Run point partitioning, mesh windowing, and BC renumbering for one rank.
     * @param rank_id MPI rank whose solid partition is being built.
     */
    void PartitionProcess(int rank_id) override;

    /**
     * @brief Write partitioned solid boundary conditions (usbc/vsbc/wsbc).
     * @param outfile Output stream for this rank's solid boundary-condition blocks.
     */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /**
     * @brief Write partitioned solid particle records.
     * @param outfile Output stream for this rank's solid particle records.
     */
    void WritePointData(std::ofstream &outfile) override;
};
