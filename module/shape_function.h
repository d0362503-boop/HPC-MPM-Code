#pragma once

#include <array>
#include <vector>
#include "mesh.h"
#include "dataset.h"
#include "material_point.h"
#include "cal_mat.h"

void MakSf(const int m, const std::array<double, 3>& xyg, const std::vector<int>& idimc, const std::vector<int>& xynodec,
           std::vector<int>& nc, int& nenode, std::vector<double>& sf, std::vector<std::array<double, 3>>& dsf);
