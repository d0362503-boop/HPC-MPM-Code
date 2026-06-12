#pragma once

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

void GaussianDistribution(std::array<std::array<double, 3>, 6> &dec2p);

void MakNodalVol();

void BuildMesh();

void BuildControlPoint();

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
