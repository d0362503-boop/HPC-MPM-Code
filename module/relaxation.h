#pragma once

#include "bc.h"
#include <vector>

using namespace std;

void Anderson_relaxation_M2(int block_it, const std::vector<double> &nvel_k, std::vector<double> &u_s_old,
                            std::vector<double> &u_s_older, std::vector<double> &r_k_old,
                            std::vector<double> &r_k_older, BoundaryCondition &fsi_intf);

void Anderson_relaxation_M1(int block_it, std::vector<double> &u_s_old, const std::vector<double> &nvel_k,
                            std::vector<double> &r_k_old, BoundaryCondition &fsi_intf);

void Aitken_relaxation(int block_it, double &omega, const std::vector<double> &nvel_k, std::vector<double> &r_k_old,
                       BoundaryCondition &fsi_intf);
