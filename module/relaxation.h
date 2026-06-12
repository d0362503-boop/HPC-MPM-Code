#pragma once

#include <vector>
#include "bc.h"

using namespace std;

void Anderson_relaxation_M2 (int block_it, const vector<double>& nvel_k,
                             vector<double>& u_s_old, vector<double>& u_s_older,
                             vector<double>& r_k_old, vector<double>& r_k_older,
                             BoundaryCondition& fsi_intf); 

void Anderson_relaxation_M1 (int block_it, vector<double>& u_s_old, 
                             const vector<double>& nvel_k, vector<double>& r_k_old,
                             BoundaryCondition& fsi_intf);

void Aitken_relaxation (int block_it, double& omega, 
                        const vector<double>& nvel_k, vector<double>& r_k_old,
                        BoundaryCondition& fsi_intf);
