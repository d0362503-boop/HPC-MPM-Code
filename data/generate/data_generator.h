#pragma once

#include <fstream>
#include <string>

class DataGenerator {
  public:
    /** @brief Destroy the generator after its case data has been written. */
    virtual ~DataGenerator() = default;

    /**
     * @brief Execute the complete input-generation workflow.
     *
     * @return Zero when the case data has been generated successfully.
     */
    int Run();

  protected:
    /** @brief Return the physical case name used in generator messages. */
    virtual std::string CaseName() const = 0;

    /** @brief Read the physical parameters that define this case. */
    virtual void LoadInput() = 0;

    /**
     * @brief Build the mesh, boundary conditions, and material points.
     *
     * Common workflow shared by all cases: initialize grid geometry, build the
     * mesh and control points, then create BCs and particles. Case-specific
     * behavior is customized through the CreateBCs/CreateParticles hooks.
     */
    void BuildData();

    /** @brief Create this case's control-point and inflow boundary conditions. */
    virtual void CreateBCs();

    /** @brief Create this case's initial material points and point properties. */
    virtual void CreateParticles();

    /**
     * @brief Write this case's boundary conditions to the grid-data stream.
     * @param outfile Output stream for the generated grid-data file.
     */
    virtual void WriteBCData(std::ofstream &outfile) = 0;

    /**
     * @brief Write the common grid file and this case's particle text file.
     *
     * Coupled cases with a custom aggregate file format may override this
     * workflow directly.
     */
    virtual void WriteTextOutputs();

    /**
     * @brief Write case-specific particle records to an open output stream.
     * @param pointfile Output stream for the particle-data file.
     */
    virtual void WritePointData(std::ofstream &pointfile);

    /** @brief Write visualization files for the generated mesh and points. */
    virtual void WriteVisualizationOutputs() = 0;

    /** @brief Derive global mesh dimensions from the loaded mesh geometry. */
    void InitializeGridGeometry();

    /**
     * @brief Write the common mesh header and case boundary conditions.
     * @param filename Destination grid-data filename.
     */
    void WriteGridDataFile(const std::string &filename);
};
