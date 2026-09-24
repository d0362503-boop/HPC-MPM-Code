#pragma once

#include "data/generate/data_generator.h"

#include "module/material_point.h"

#include <string>

class SolidGenerator : public DataGenerator, public MaterialPoint {
  protected:
    /** @brief Return the solid case label. */
    std::string CaseName() const override;

    /** @brief Read solid material, mesh, and time-integration parameters. */
    void LoadInput() override;

    /** @brief Create solid control-point boundary conditions. */
    void CreateBCs() override;

    /** @brief Create solid material points and their initial properties. */
    void CreateParticles() override;

    /**
     * @brief Write solid control-point boundary conditions.
     * @param outfile Output stream for the generated grid-data file.
     */
    void WriteBCData(std::ofstream &outfile) override;

    /**
     * @brief Write solid particle records to the open particle-data file.
     * @param pointfile Output stream for the solid particle-data file.
     */
    void WritePointData(std::ofstream &pointfile) override;

    /** @brief Write solid mesh and particle visualization files. */
    void WriteVisualizationOutputs() override;

    double solid_point_mass_ = 0;   // mass per solid particle
    double solid_point_volume_ = 0; // volume per solid particle
};
