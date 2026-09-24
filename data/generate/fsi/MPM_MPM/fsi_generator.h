#pragma once

#include "data/generate/fluid/fluid_generator.h"
#include "data/generate/solid/solid_generator.h"

#include <fstream>
#include <string>

class MPMMPMFSIFluidGenerator : public FluidGenerator {
  public:
    /** Create MPM-fluid boundary conditions for the coupled domain. */
    void CreateBCs() override;

    /** Create MPM-fluid particles in the fluid portion of the coupled domain. */
    void CreateParticles() override;

    // BC and particle records are identical to the standalone generators'.
    using FluidGenerator::WriteBCData;
    using FluidGenerator::WritePointData;
};

class MPMMPMFSISolidGenerator : public SolidGenerator {
  public:
    /** Create MPM-solid boundary conditions for the coupled domain. */
    void CreateBCs() override;

    /** Create MPM-solid particles in the solid portion of the coupled domain. */
    void CreateParticles() override;

    // BC and particle records are identical to the standalone generators'.
    using SolidGenerator::WriteBCData;
    using SolidGenerator::WritePointData;

  private:
    /** Mark the upper solid surface for FSI interface detection. */
    void MarkSurfacePoints();
};

class MPMMPMFSIGenerator : public DataGenerator {
  protected:
    /** Return the MPM-MPM FSI case label. */
    std::string CaseName() const override;

    /** Read coupled mesh, time, fluid, and solid input parameters. */
    void LoadInput() override;

    /** Create coupled fluid and solid boundary conditions. */
    void CreateBCs() override;

    /** Create coupled fluid and solid material points. */
    void CreateParticles() override;

    /**
     * Write coupled fluid and solid boundary conditions.
     *
     * @param outfile Stream receiving all coupled boundary-condition records.
     */
    void WriteBCData(std::ofstream &outfile) override;

    /**
     * Write coupled solid and fluid particle records.
     *
     * @param pointfile Stream receiving all coupled particle records.
     */
    void WritePointData(std::ofstream &pointfile) override;

    /** Write coupled mesh and particle visualization files. */
    void WriteVisualizationOutputs() override;

  private:
    MPMMPMFSIFluidGenerator fluid_; // coupled fluid data
    MPMMPMFSISolidGenerator solid_; // coupled solid data
};
