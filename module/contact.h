#pragma once

#include "module/bc.h"
#include "module/material_point.h"
#include <vector>

void SearchClosestPoint2Node(int nid, MaterialPoint &point, double &distance);

void SearchFluidSolidIntf();
