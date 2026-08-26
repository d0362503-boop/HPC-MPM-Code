#pragma once

#include "data/generate/solid/solid_generator.h"
#include "module/bc.h"

#include <fstream>
#include <string>

class MPMFEMFSIGenerator : public SolidGenerator {
  protected:
    /** Return the MPM-FEM FSI case label. */
    std::string CaseName() const override;

    /** Read coupled mesh, time, and solid input parameters. */
    void LoadInput() override;

    /** Create solid and FEM-fluid boundary conditions. */
    void CreateBCs() override;

    /** Create MPM-solid particles in the coupled domain. */
    void CreateParticles() override;

    /**
     * Write solid and FEM-fluid boundary conditions.
     *
     * @param outfile Stream receiving all coupled boundary-condition records.
     */
    void WriteBCData(std::ofstream &outfile) override;

    /**
     * Write MPM-solid particle records.
     *
     * @param pointfile Stream receiving all coupled solid-particle records.
     */
    void WritePointData(std::ofstream &pointfile) override;

    /** Write coupled mesh and solid-particle visualization files. */
    void WriteVisualizationOutputs() override;

  private:
    /** Create FEM-fluid velocity and pressure boundary conditions. */
    void CreateFluidBCs();

    /** Mark the upper solid surface for FSI interface detection. */
    void MarkSurfacePoints();

    BoundaryCondition fluid_ubc_; // FEM-fluid x-velocity BCs
    BoundaryCondition fluid_vbc_; // FEM-fluid y-velocity BCs
    BoundaryCondition fluid_wbc_; // FEM-fluid z-velocity BCs
    BoundaryCondition fluid_pbc_; // FEM-fluid pressure BCs
};
