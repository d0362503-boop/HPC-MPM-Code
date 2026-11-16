#pragma once

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include <array>
#include <vector>

/**
 * @brief Evaluate shape functions and gradients at point `xyg` for element `m`.
 * @param m        Element index.
 * @param xyg      Physical coordinates of the query point.
 * @param idimc    Shape-function order in each direction.
 * @param xynodec  Control-point count in each direction.
 * @param nc       Output node IDs of the supporting element nodes.
 * @param nenode   Output number of supporting nodes.
 * @param sf       Output shape-function values.
 * @param dsf      Output shape-function gradients.
 */
void MakeSF(const int m, const std::array<double, 3> &xyg, const std::vector<int> &idimc,
            const std::vector<int> &xynodec, std::vector<int> &nc, int &nenode, std::vector<double> &sf,
            std::vector<std::array<double, 3>> &dsf);
