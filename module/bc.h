#pragma once

#include <fstream>
#include <iostream>
#include <vector>

class BoundaryCondition {
  public:
    int ibc = 0;
    std::vector<int> nbc;
    std::vector<double> fbc;

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
};
