#pragma once

#include "data/divide/divide_fluid/fluid_divider.h"
#include "data/divide/divide_solid/solid_divider.h"

#include <string>

class MPMMPMFSISolidDivider : public SolidDivider {
  public:
    using SolidDivider::SolidDivider;

    // Expose the base hook implementations for the FSI coordinator.
    using SolidDivider::LoadBoundaryData;
    using SolidDivider::LoadPointData;
    using SolidDivider::WriteBoundaryData;
    using SolidDivider::WritePointData;
};

class MPMMPMFSIFluidDivider : public FluidDivider {
  public:
    using FluidDivider::FluidDivider;

    // Expose the base hook implementations for the FSI coordinator.
    using FluidDivider::LoadBoundaryData;
    using FluidDivider::LoadPointData;
    using FluidDivider::WriteBoundaryData;
    using FluidDivider::WritePointData;
};

class MPMMPMFSIDivider : public DataPartitioner {
  public:
    /**
     * @brief Build the coupled partitioner; solid and fluid share `grid` and `params`.
     * @param grid   Grid being partitioned (global mesh, overwritten with each rank window).
     * @param params Run parameters of the case.
     */
    MPMMPMFSIDivider(Grid &grid, RunParams &params)
        : DataPartitioner(grid, params), solid_(grid, params), fluid_(grid, params) {}

  protected:
    /** @brief Return the MPM-MPM FSI case label. */
    std::string CaseName() const override;

    /**
     * @brief Load solid then fluid boundary conditions (file block order).
     * @param infile Input stream positioned at the first solid boundary-condition block.
     */
    void LoadBoundaryData(std::ifstream &infile) override;

    /**
     * @brief Load solid then fluid particle records (file block order).
     * @param infile Input stream positioned at the global solid particle section.
     */
    void LoadPointData(std::ifstream &infile) override;

    /**
     * @brief Partition both physics for one rank.
     * @param rank_id Index of the MPI rank whose partition is being built.
     *
     * Point renumbering runs per physics before the shared MeshPartition;
     * BC renumbering runs per physics after it.
     */
    void PartitionProcess(int rank_id) override;

    /**
     * @brief Write solid then fluid partitioned boundary conditions.
     * @param outfile Output stream for this rank's solid and fluid boundary-condition blocks.
     */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /**
     * @brief Write solid then fluid partitioned particle records.
     * @param outfile Output stream for this rank's solid and fluid particle records.
     */
    void WritePointData(std::ofstream &outfile) override;

  private:
    MPMMPMFSISolidDivider solid_; // solid-side partition component
    MPMMPMFSIFluidDivider fluid_; // fluid-side partition component
};
