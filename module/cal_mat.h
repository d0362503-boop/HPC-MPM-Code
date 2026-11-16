#pragma once

#include "module/dataset.h"
#include "module/mesh.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

/**
 * @brief Return the sign of a scalar value.
 * @param val Input value.
 * @return -1, 0, or 1 depending on the sign of `val`.
 */
template <typename T> T Sign(T val) { return (T(0) < val) - (val < T(0)); }

/**
 * @brief Compute the scalar product of two three-component physical vectors.
 * @param a First physical vector.
 * @param b Second physical vector.
 * @return Scalar product in the corresponding physical units.
 */
template <typename T> double DotVec3(const std::array<T, 3> &a, const std::array<T, 3> &b) {
    double product = 0.0;
    for (int d = 0; d < 3; ++d) { product += a[d] * b[d]; }

    return product;
}

/**
 * @brief Compute the Euclidean norm of a three-component physical vector.
 * @param vec Physical vector whose magnitude is required.
 * @return Vector magnitude in the same physical units.
 */
template <typename T> double NormVec3(const std::array<T, 3> &vec) {
    return std::sqrt(DotVec3(vec, vec));
}

/**
 * @brief Subtract physical positions or vectors component by component.
 * @param a Position or vector being reduced.
 * @param b Position or vector being subtracted.
 * @return Relative position or vector.
 */
template <typename T> std::array<T, 3> DifferenceVec3(const std::array<T, 3> &a, const std::array<T, 3> &b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

/**
 * @brief Compute the vector product of two three-component physical vectors.
 * @param a First physical vector or edge.
 * @param b Second physical vector or edge.
 * @return Oriented vector; its magnitude is twice the triangle area for two edges.
 */
template <typename T> std::array<T, 3> CrossVec3(const std::array<T, 3> &a, const std::array<T, 3> &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

/**
 * @brief Compute the trace of a 3x3 matrix.
 * @param mat Input matrix.
 * @return Trace of `mat`.
 */
template <typename T> double TraceMat3(const std::array<std::array<T, 3>, 3> &mat) {
    return double(mat[0][0]) + double(mat[1][1]) + double(mat[2][2]);
}

/**
 * @brief Compute the determinant of a 3x3 matrix.
 * @param mat Input matrix.
 * @return Determinant of `mat`.
 */
template <typename T> double DetMat3(const std::array<std::array<T, 3>, 3> &mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) +
           mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) +
           mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

/**
 * @brief Compute the inverse of a 3x3 matrix.
 * @param mat Input matrix.
 * @return Inverse of `mat`.
 * @note Requires a nonzero determinant.
 */
template <typename T> std::array<std::array<double, 3>, 3> InvMat3(const std::array<std::array<T, 3>, 3> &mat) {
    const double inv_det = 1.0 / DetMat3(mat);

    return {{{(mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * inv_det,
              (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * inv_det,
              (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * inv_det},
             {(mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * inv_det,
              (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * inv_det,
              (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * inv_det},
             {(mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * inv_det,
              (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * inv_det,
              (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * inv_det}}};
}

/**
 * @brief Transpose a 3x3 matrix in-place.
 * @param mat Input/output matrix.
 */
template <typename T> void TransMat3(std::array<std::array<T, 3>, 3> &mat) {
    std::swap(mat[0][1], mat[1][0]);
    std::swap(mat[0][2], mat[2][0]);
    std::swap(mat[1][2], mat[2][1]);
}
