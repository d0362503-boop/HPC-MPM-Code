#pragma once

#include "data/generate/data_generator.h"

#include "module/material_point.h"

#include <string>

class FluidGenerator : public DataGenerator, public MaterialPoint {
  protected:
    /** @brief Return the fluid case label. */
    std::string CaseName() const override;

    /** @brief Read fluid material, mesh, and time-integration parameters. */
    void LoadInput() override;

    /** @brief Create fluid control-point and inflow boundary conditions. */
    void CreateBCs() override;

    /** @brief Create fluid material points and their initial properties. */
    void CreateParticles() override;

    /**
     * @brief Write fluid control-point and inflow boundary conditions.
     * @param outfile Output stream for the generated grid-data file.
     */
    void WriteBCData(std::ofstream &outfile) override;

    /**
     * @brief Write fluid particle records to the open particle-data file.
     * @param pointfile Output stream for the fluid particle-data file.
     */
    void WritePointData(std::ofstream &pointfile) override;

    /** @brief Write fluid mesh and particle visualization files. */
    void WriteVisualizationOutputs() override;

    double fluid_point_mass_ = 0;   // mass per fluid particle
    double fluid_point_volume_ = 0; // volume per fluid particle
};
