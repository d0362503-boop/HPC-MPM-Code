#include "bc.h"
#include "dataset.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void BoundaryCondition::BCInput(std::ifstream &infile, bool if_fbc) {

    infile >> this->ibc;
    infile.ignore(1000, '\n');
    this->nbc.resize(this->ibc);
    if (if_fbc) this->fbc.resize(this->ibc);

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            if (if_fbc) {
                infile >> this->nbc[i] >> this->fbc[i];
            } else {
                infile >> this->nbc[i];
            }
            infile.ignore(1000, '\n');
        }
    }

    return;
}

void BoundaryCondition::BCOutput(std::ofstream &outfile, std::string bc_name, bool if_fbc) {

    outfile << std::setw(15) << this->ibc << std::setw(15) << "! ---- " << bc_name << ".ibc ----" << "\n";

    for (int i = 0; i < this->ibc; i++) {
        if (if_fbc) {
            outfile << std::setw(15) << this->nbc[i] << std::setw(15) << this->fbc[i] << "\n";
        } else {
            outfile << std::setw(15) << this->nbc[i] << "\n";
        }
    }

    return;
}

void BoundaryCondition::BCSetVal(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = this->fbc[i] * facl;
        }
    }

    return;
}

void BoundaryCondition::BCSetZero(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = 0.0e0;
        }
    }

    return;
}

void BoundaryCondition::BCSetDt(int nn, std::vector<double> &variable) {

    if (this->ibc != 0) {
        for (int i = 0; i < this->ibc; i++) {
            int n = this->nbc[i] + nn;
            variable[n] = this->fbc[i] * dt * facl;
        }
    }

    return;
}
