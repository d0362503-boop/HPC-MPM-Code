#pragma once

#include "module/bc.h"
#include <vector>

using namespace std;

/**
 * @brief Apply depth-two Anderson relaxation to the coupled interface field.
 *
 * @param iteration Current block-coupling iteration index.
 * @param nodal_var Current nodal field values from the active sub-solver.
 * @param var_old Interface field from the preceding coupling iteration.
 * @param var_older Interface field from two coupling iterations earlier.
 * @param res_old Residual from the preceding coupling iteration.
 * @param res_older Residual from two coupling iterations earlier.
 * @param intf_bc Interface boundary data; its prescribed field is updated in place.
 */
void Anderson_relaxation_M2(int iteration, const std::vector<double> &nodal_var, std::vector<double> &var_old,
                            std::vector<double> &var_older, std::vector<double> &res_old,
                            std::vector<double> &res_older, BoundaryCondition &intf_bc);

/**
 * @brief Apply depth-one Anderson relaxation to the coupled interface field.
 *
 * @param iteration Current block-coupling iteration index.
 * @param var_old Interface field from the preceding coupling iteration.
 * @param nodal_var Current nodal field values from the active sub-solver.
 * @param res_old Residual from the preceding coupling iteration.
 * @param intf_bc Interface boundary data; its prescribed field is updated in place.
 */
void Anderson_relaxation_M1(int iteration, std::vector<double> &var_old, const std::vector<double> &nodal_var,
                            std::vector<double> &res_old, BoundaryCondition &intf_bc);

/**
 * @brief Apply Aitken dynamic relaxation to the coupled interface field.
 *
 * @param iteration Current block-coupling iteration index.
 * @param relax_coef Relaxation coefficient, updated from the residual history.
 * @param nodal_var Current nodal field values from the active sub-solver.
 * @param res_old Residual from the preceding coupling iteration.
 * @param intf_bc Interface boundary data; its prescribed field is updated in place.
 */
void Aitken_relaxation(int iteration, double &relax_coef, const std::vector<double> &nodal_var,
                       std::vector<double> &res_old, BoundaryCondition &intf_bc);
