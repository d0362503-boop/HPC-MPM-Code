#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class BoundaryCondition {
  public:
    int ibc = 0;
    std::vector<int> nbc;
    std::vector<double> fbc;
    std::vector<int> global_nbc;    // immutable global IDs
    std::vector<double> global_fbc; // immutable prescribed values

    /**
     * @brief Read boundary-condition node IDs and prescribed values from input stream.
     * @param infile  Input file stream.
     * @param if_fbc  If true, also read prescribed values into `fbc`.
     */
    void BCInput(std::ifstream &infile, bool if_fbc = true);

    /**
     * @brief Write boundary-condition node IDs and values to an output stream.
     * @param outfile   Output file stream.
     * @param bc_name   Name of the boundary condition (used in the header).
     * @param ife_fbc   If true, also write prescribed values.
     */
    void BCOutput(std::ofstream &outfile, std::string bc_name, bool ife_fbc = true);

    /**
     * @brief Apply prescribed boundary values to `variable` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Vector to be modified in-place.
     */
    void BCSetVal(int nn, std::vector<double> &variable);

    /**
     * @brief Zero out `variable` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Vector to be modified in-place.
     */
    void BCSetZero(int nn, std::vector<double> &variable);

    /**
     * @brief Apply prescribed displacement/velocity increments to `variable` at constrained DOFs.
     * @param nn       Offset into `variable` for the component being constrained.
     * @param variable Vector to be modified in-place.
     */
    void BCSetDt(int nn, std::vector<double> &variable);

    /**
     * @brief Cache this rank's control-point BCs as a global MPI-wide BC list.
     *
     * Local control-point IDs are converted to global IDs and gathered from all
     * ranks.  The resulting data remains valid after a DLB repartition.
     */
    void CaptureGlobalControlPointBC();

    /**
     * @brief Rebuild local control-point BC IDs for the current MPI region.
     */
    void RebuildLocalControlPointBC();

    /**
     * @brief Cache this rank's element BCs as a global MPI-wide BC list.
     *
     * This is used by fluid inflow BCs, whose IDs refer to background elements.
     */
    void CaptureGlobalElementBC();

    /**
     * @brief Rebuild local background-element BC IDs for the current MPI region.
     */
    void RebuildLocalElementBC();

  private:
    /**
     * @brief Gather and de-duplicate global BC entries from all MPI ranks.
     * @param local_global_ids Global IDs represented by this rank's local BC entries.
     * @param has_values Whether each BC entry has a prescribed scalar value.
     */
    void CacheGlobalEntries(const std::vector<int> &local_global_ids, bool has_values);
};
