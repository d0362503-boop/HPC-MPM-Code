#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "module/bc.h"
#include "module/material_point.h"

class DataPartitioner {
  public:
    virtual ~DataPartitioner() = default;

    /**
     * @brief Execute the complete partition workflow.
     * @return Zero when the partitioned data has been written successfully.
     */
    int Run();

    /**
     * @brief Filter boundary-condition entries to a rank window and convert to local indexing.
     * @param partition_bc  Output BC holding the rank-local entries.
     * @param source_bc     Input BC with globally indexed entries.
     * @param local_counts  Entries per direction in the rank window.
     * @param global_counts Entries per direction in the global mesh.
     * @param local_min     First global index of the rank window, per direction.
     * @param local_max     Last global index (inclusive) of the rank window, per direction.
     * @param compute_values Whether to also carry over the prescribed values.
     */
    static void BCRenumber(BoundaryCondition &partition_bc,
                           const BoundaryCondition &source_bc,
                           const std::vector<int> &local_counts,
                           const std::vector<int> &global_counts,
                           const std::vector<int> &local_min,
                           const std::vector<int> &local_max,
                           bool compute_values = true);

  protected:
    /** @brief Return the physical case name used in partitioner messages. */
    virtual std::string CaseName() const = 0;

    /**
     * @brief Load this case's boundary conditions from the global grid-data stream.
     * @param infile Input stream positioned at the boundary-condition blocks.
     */
    virtual void LoadBoundaryData(std::ifstream &infile) = 0;

    /**
     * @brief Load this case's material points from the global point-data stream.
     * @param infile Input stream positioned at the particle records.
     */
    virtual void LoadPointData(std::ifstream &infile) = 0;

    /**
     * @brief Partition this case's points and boundary conditions for one rank.
     * @param rank_id Index of the MPI rank whose partition is being built.
     */
    virtual void PartitionProcess(int rank_id) = 0;

    /**
     * @brief Write the partitioned boundary conditions for the current rank.
     * @param outfile Output stream for the rank's grid-data file.
     */
    virtual void WriteBoundaryData(std::ofstream &outfile) = 0;

    /**
     * @brief Write the partitioned material points for the current rank.
     * @param outfile Output stream for the rank's point-data file.
     */
    virtual void WritePointData(std::ofstream &outfile) = 0;

    /**
     * @brief Read the global mesh header (dimensions, element/node counts).
     * @param infile Input stream of the global grid-data file.
     */
    void InputMeshData(std::ifstream &infile);

    /**
     * @brief Write the per-rank mesh header and overlap information.
     * @param outfile Output stream for the rank's grid-data file.
     * @param rank_id Index of the MPI rank whose data is written.
     */
    void OutputMeshData(std::ofstream &outfile, int rank_id) const;

    /** @brief Derive per-rank window sizes from the global mesh and rank topology. */
    void PartitionInitialDataset();

    /**
     * @brief Collect the points falling inside a rank window and build the local-to-global map.
     * @param local_num Output number of points owned by the rank.
     * @param global_id Output local-to-global point index map.
     * @param point     Source material points in global indexing.
     * @param rank_id   Index of the MPI rank whose partition is being built.
     */
    void PointRenumber(int &local_num, std::vector<int> &global_id, MaterialPoint &point, int rank_id) const;

    /**
     * @brief Set up the element/control-point window and overlap metadata for one rank.
     * @param rank_id Index of the MPI rank whose partition is being built.
     */
    void MeshPartition(int rank_id);

    std::vector<int> nexyr1_{3}, nexyr2_{3};     // per-direction rank element counts (coarse/fine)
    std::vector<int> nr1_{3}, nr2_{3};           // rank counts with coarse/fine division
    std::vector<int> nxyr1_{3}, nxyr2_{3};       // per-direction rank node counts (coarse/fine)
    std::vector<double> drxy1_{3}, drxy2_{3};    // rank window sizes (coarse/fine)
    std::vector<double> bound12_{3};             // coarse/fine window boundary coordinates
    std::vector<int> inxyminc_{3}, inxymaxc_{3}; // rank control-point window (inclusive)
    std::vector<int> inxyminl_{3}, inxymaxl_{3}; // rank element window (min, max+1)

  private:
    /**
     * @brief Read the rank topology and input/output file names.
     * @param infile Input stream of the partition parameter file.
     */
    void LoadPartitionParameters(std::ifstream &infile);

    /** @brief Create the per-rank output directory if missing. */
    void CreateOutputDirectories() const;

    std::string grid_input_file_;     // global grid-data file
    std::string point_input_file_;    // global point-data file
    std::string grid_output_prefix_;  // per-rank grid-data prefix
    std::string point_output_prefix_; // per-rank point-data prefix
};
