#pragma once

#include "data/divide/divide_fluid/fluid_divider.h"
#include "data/divide/divide_solid/solid_divider.h"

#include <string>

class MPMMPMFSISolidDivider : public SolidDivider {
  public:
    // Expose the base hook implementations for the FSI coordinator.
    using SolidDivider::LoadBoundaryData;
    using SolidDivider::LoadPointData;
    using SolidDivider::WriteBoundaryData;
    using SolidDivider::WritePointData;
};

class MPMMPMFSIFluidDivider : public FluidDivider {
  public:
    // Expose the base hook implementations for the FSI coordinator.
    using FluidDivider::LoadBoundaryData;
    using FluidDivider::LoadPointData;
    using FluidDivider::WriteBoundaryData;
    using FluidDivider::WritePointData;
};

class MPMMPMFSIDivider : public DataPartitioner {
  protected:
    /** @brief Return the MPM-MPM FSI case label. */
    std::string CaseName() const override;

    /** @brief Load solid then fluid boundary conditions (file block order). */
    void LoadBoundaryData(std::ifstream &infile) override;

    /** @brief Load solid then fluid particle records (file block order). */
    void LoadPointData(std::ifstream &infile) override;

    /**
     * @brief Partition both physics for one rank.
     * @param rank_id Index of the MPI rank whose partition is being built.
     *
     * Point renumbering runs per physics before the shared MeshPartition;
     * BC renumbering runs per physics after it.
     */
    void PartitionProcess(int rank_id) override;

    /** @brief Write solid then fluid partitioned boundary conditions. */
    void WriteBoundaryData(std::ofstream &outfile) override;

    /** @brief Write solid then fluid partitioned particle records. */
    void WritePointData(std::ofstream &outfile) override;

  private:
    MPMMPMFSISolidDivider solid_; // solid-side partition component
    MPMMPMFSIFluidDivider fluid_; // fluid-side partition component
};
