#pragma once

#include <vector>
#include "material_point.h"
#include "bc.h"

void SearchClosestPoint2Node(int nid, MaterialPoint& point, double& distance);

void SearchFluidSolidIntf(); 

void RigidWallBCCheck ();
