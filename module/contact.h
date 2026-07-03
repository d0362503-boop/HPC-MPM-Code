#pragma once

#include "bc.h"
#include "material_point.h"
#include <vector>

void SearchClosestPoint2Node(int nid, MaterialPoint &point, double &distance);

void SearchFluidSolidIntf();
