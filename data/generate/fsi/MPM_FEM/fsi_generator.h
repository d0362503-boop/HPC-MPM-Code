#pragma once

#include "data/generate/solid/solid_generator.h"

#include <fstream>
#include <string>

class MPMFEMFSISolidGenerator : public SolidGenerator {
  public:
    using SolidGenerator::SolidGenerator;

    /** Create MPM-solid boundary conditions for the coupled domain. */
    void CreateBCs() override;

    /** Create MPM-solid particles in the solid portion of the coupled domain. */
    void CreateParticles() override;

    // BC and particle records are identical to the standalone solid generator's.
    using SolidGenerator::WriteBCData;
    using SolidGenerator::WritePointData;

  private:
    /** Mark the upper solid surface for FSI interface detection. */
    void MarkSurfacePoints();
};

class MPMFEMFSIFluidGenerator {
  public:
    /**
     * @brief Bind the FEM-fluid boundary-condition set to the grid being generated.
     * @param grid   Grid being generated (global mesh of the case).
     * @param params Run parameters read from input.txt.
     */
    MPMFEMFSIFluidGenerator(Grid &grid, RunParams &params) : grid_(grid), params_(params), points_(grid, params) {}

    Grid &grid_;           // grid being generated
    RunParams &params_;    // run parameters of the case
    MaterialPoint points_; // FEM-fluid boundary conditions

    /** Create FEM-fluid velocity and pressure boundary conditions on the mesh control points. */
    void CreateBCs();

    /**
     * Write FEM-fluid boundary conditions.
     *
     * @param outfile Stream receiving the coupled boundary-condition records.
     */
    void WriteBCData(std::ofstream &outfile);
};

class MPMFEMFSIGenerator : public DataGenerator {
  public:
    /**
     * @brief Build the coupled generator; solid and FEM fluid share `grid` and `params`.
     * @param grid   Grid being generated (global mesh of the case).
     * @param params Run parameters read from input.txt.
     */
    MPMFEMFSIGenerator(Grid &grid, RunParams &params)
        : DataGenerator(grid, params), solid_(grid, params), fluid_(grid, params) {}

  protected:
    /** Return the MPM-FEM FSI case label. */
    std::string CaseName() const override;

    /** Read coupled mesh, time, and solid input parameters. */
    void LoadInput() override;

    /** Create solid and FEM-fluid boundary conditions. */
    void CreateBCs() override;

    /** Create MPM-solid particles (the FEM fluid has no particles). */
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
    MPMFEMFSISolidGenerator solid_; // coupled solid data
    MPMFEMFSIFluidGenerator fluid_; // coupled fluid data
};
