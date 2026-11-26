#pragma once

#include "data/divide/divide_solid/solid_divider.h"

#include <string>

class MPMFEMFSISolidDivider : public SolidDivider {
  public:
    using SolidDivider::SolidDivider;

    // Expose the base hook implementations for the FSI coordinator.
    using SolidDivider::LoadBoundaryData;
    using SolidDivider::LoadPointData;
    using SolidDivider::WriteBoundaryData;
    using SolidDivider::WritePointData;
};

class MPMFEMFSIFluidDivider {
  public:
    /**
     * @brief Bind the FEM-fluid BC sets to the grid being partitioned.
     * @param grid   Grid being partitioned (global mesh, overwritten with each rank window).
     * @param params Run parameters of the case.
     */
    MPMFEMFSIFluidDivider(Grid &grid, RunParams &params)
        : grid_(grid), params_(params), mesh_(grid, params), partition_mesh_(grid, params) {}

    Grid &grid_;        // grid being partitioned
    RunParams &params_; // run parameters of the case

    /** @brief Load FEM-fluid velocity/pressure BCs (uwbc/vwbc/wwbc/hpbc). */
    void LoadBoundaryData(std::ifstream &infile);

    /**
     * @brief Renumber FEM-fluid BCs to this rank's local control-point indexing.
     * @param cp_min First global control-point index of the rank window, per direction.
     * @param cp_max Last global control-point index (inclusive) of the rank window, per direction.
     *
     * Must be called after the coordinator's MeshPartition.
     */
    void PartitionBCs(const std::array<int, 3> &cp_min, const std::array<int, 3> &cp_max);

    /** @brief Write partitioned FEM-fluid BCs (uwbc/vwbc/wwbc/hpbc). */
    void WriteBoundaryData(std::ofstream &outfile);

    MaterialPoint mesh_;           // loaded FEM-fluid mesh BCs
    MaterialPoint partition_mesh_; // partitioned FEM-fluid mesh BCs
};

class MPMFEMFSIDivider : public DataPartitioner {
  public:
    /**
     * @brief Build the coupled partitioner; solid and FEM fluid share `grid` and `params`.
     * @param grid   Grid being partitioned (global mesh, overwritten with each rank window).
     * @param params Run parameters of the case.
     */
    MPMFEMFSIDivider(Grid &grid, RunParams &params)
        : DataPartitioner(grid, params), solid_(grid, params), fluid_(grid, params) {}

  protected:
    /** @brief Return the MPM-FEM FSI case label. */
    std::string CaseName() const override;

    /** @brief Load solid then FEM-fluid boundary conditions (file block order). */
    void LoadBoundaryData(std::ifstream &infile) override;

    /** @brief Load solid particle records (the FEM fluid has no particles). */
    void LoadPointData(std::ifstream &infile) override;

    /**
     * @brief Partition both physics for one rank.
     * @param rank_id Index of the MPI rank whose partition is being built.
     *
     * Solid point renumbering runs before the shared MeshPartition;
     * BC renumbering runs per physics after it.
     */
    void PartitionProcess(int rank_id) override;

    /** @brief Write solid then FEM-fluid partitioned boundary conditions. */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /** @brief Write partitioned solid particle records. */
    void WritePointData(std::ofstream &outfile) override;

  private:
    MPMFEMFSISolidDivider solid_; // solid-side partition component
    MPMFEMFSIFluidDivider fluid_; // FEM-fluid partition component
};
