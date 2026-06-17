#pragma once

#include "dataset.h"
#include <array>
#include <vector>

// ----- Mesh & node data -----
inline int node, nodew, nelem, nelemw, nu, nv, nw, np;
inline std::vector<int> xynode(3), xyelem(3), aelemmin(3), aelemmax(3);
inline std::vector<int> xynodew(3), xyelemw(3);
inline std::vector<double> dxy(3), xymin(3), xymax(3), xyminw(3), xymaxw(3);
inline std::vector<double> rhoe, rmue;
inline std::vector<std::array<int, 8>> nc;
inline std::vector<std::array<double, 3>> xyn;
// ----------------------------

// ----- Control point data -----
inline int nodec, nodecw, nuc, nvc, nwc, npc;
inline std::vector<int> xynodec(3), xynodecw(3), idimc(3), idiml(3);
inline std::vector<double> nvol;
inline std::vector<std::vector<int>> ncc;
inline std::vector<std::array<double, 3>> xyc;
// ----------------------------

inline void InitalizeMeshAndTimeParameters() {

    dti = 1.0e0 / dt;
    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        xynodew[i] = xyelemw[i] + 1;
        xynodecw[i] = xyelemw[i] + idimc[i];
    }
    dlstep = 1.0e0 / double(nlstep);

    return;
}

/**
 * @brief Compute Gaussian quadrature points and weights for hexahedral elements.
 * @param dec2p Output array of quadrature data.
 */
void GaussianDistribution(std::array<std::array<double, 3>, 6> &dec2p);

/**
 * @brief Compute nodal/control-point volumes from element volumes.
 */
void MakNodalVol();

/**
 * @brief Build the background mesh node coordinates and element connectivity.
 */
void BuildMesh();

/**
 * @brief Build the higher-order control-point coordinates and connectivity.
 */
void BuildControlPoint();

/**
 * @brief Compute the coordinate of control point `nid` along direction `dir`.
 * @tparam dir  Spatial direction (0=x, 1=y, 2=z).
 * @param nid   Control-point index.
 * @param ijk_l Local IJK indices of the control point.
 * @param ijk_g Global IJK indices of the control point.
 */
template <int dir> void DefineCpPos(int nid, const std::vector<int> &ijk_l, const std::vector<int> &ijk_g) {

    if (idimc[dir] == 1) {
        xyc[nid][dir] = xymin[dir] + dxy[dir] * double(ijk_l[dir]);
    } else if (idimc[dir] == 2) {
        if (ijk_g[dir] == 0) {
            xyc[nid][dir] = xyminw[dir];
        } else if (ijk_g[dir] == xynodecw[dir] - 1) {
            xyc[nid][dir] = xymaxw[dir];
        } else {
            xyc[nid][dir] = xymin[dir] + dxy[dir] * (double(ijk_l[dir]) - 0.5e0);
        }
    }

    return;
}
