#pragma once

#include <vector>
#include <cmath>
#include <mpi.h>
#include "crsmat.h"

/**
 * @brief Solve the linear system in `mat` using the stabilized GPBiCG method.
 * @param mat Sparse matrix system containing `amat`, `b_rhs`, and `x_lhs`.
 * @return Number of iterations, or -1 on failure.
 */
int GPBiCGSafe(CrsMat& mat);

/**
 * @brief Solve the linear system in `mat` using the GPBiCG-AR method.
 * @param mat Sparse matrix system containing `amat`, `b_rhs`, and `x_lhs`.
 * @return Number of iterations, or -1 on failure.
 */
int GPBiCGAR(CrsMat& mat);

/**
 * @brief Solve the linear system in `mat` using the standard GPBiCG method.
 * @param mat Sparse matrix system containing `amat`, `b_rhs`, and `x_lhs`.
 * @return Number of iterations, or -1 on failure.
 */
int GPBiCG(CrsMat& mat);

#include "../material_point.h"
