#pragma once

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "mesh.h"
#include "dataset.h"

template <typename T>
T Sign(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
double NormVec3(const std::array<T, 3>& vec) {

    double norm_vec = 0.0e0;
    for (int n = 0; n < 3; n++) {
        norm_vec += vec[n] * vec[n];
    }
    norm_vec = std::sqrt(norm_vec);

    return norm_vec;
}

template<typename T>
double TraceMat3(const std::array<std::array<T, 3>, 3>& mat) {

    double tr_mat = 0.0e0;
    for (int i = 0; i < 3; i++) {
        tr_mat += mat[i][i]; 
    }

    return tr_mat;
}

template<typename T>
double DetMat3(const std::array<std::array<T, 3>, 3>& mat) {

    double mat_det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) +
                     mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) +
                     mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
              
    return mat_det;
}

template<typename T>
std::array<std::array<double, 3>, 3> InvMat3(const std::array<std::array<T, 3>, 3>& mat) {

    double det = DetMat3(mat);
    det = 1.0e0 / det;

    std::array<std::array<double, 3>, 3> mat_inv;
    mat_inv[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * det;
    mat_inv[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * det;
    mat_inv[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * det;
    mat_inv[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * det;
    mat_inv[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * det;
    mat_inv[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * det;
    mat_inv[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * det;
    mat_inv[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * det;
    mat_inv[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * det;                
 
    return mat_inv;
}

template<typename T>
void TransMat3(std::array<std::array<T, 3>, 3>& mat) {

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            std::swap(mat[i][j], mat[j][i]); 
        }
    }
    return;
}
