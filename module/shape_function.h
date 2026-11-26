#pragma once

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include <array>
#include <vector>

/**
 * @brief Evaluate shape functions and gradients at point `xyg` for element `m`.
 * @param grid     Grid of this rank (mesh window/spacing and B-spline orders are used).
 * @param m        Rank-local element index.
 * @param xyg      Physical coordinates of the query point.
 * @param nc       Output control-point IDs supporting the element.
 * @param nenode   Output number of supporting control points.
 * @param sf       Output shape-function values.
 * @param dsf      Output shape-function gradients.
 */
void MakeSF(const Grid &grid, const int m, const std::array<double, 3> &xyg, std::vector<int> &nc, int &nenode,
            std::vector<double> &sf, std::vector<std::array<double, 3>> &dsf);
