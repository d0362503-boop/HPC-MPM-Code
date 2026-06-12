#pragma once

#include <fstream>
#include <iostream>
#include <vector>

class BoundaryCondition {
  public:
    int ibc = 0;
    std::vector<int> nbc;
    std::vector<double> fbc;

    void BCInput(std::ifstream &infile, bool if_fbc = true);

    void BCOutput(std::ofstream &outfile, std::string bc_name, bool ife_fbc = true);

    void BCSetVal(int nn, std::vector<double> &variable);

    void BCSetZero(int nn, std::vector<double> &variable);

    void BCSetDt(int nn, std::vector<double> &variable); // uses mp_ if non-null
};
