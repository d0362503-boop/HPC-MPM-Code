#pragma once

#include <vector>
#include <cmath>
#include <mpi.h>
#include "crsmat.h"

int GPBiCGSafe(CrsMat& mat);

int GPBiCGAR(CrsMat& mat);

int GPBiCG(CrsMat& mat);

#include "../material_point.h"
