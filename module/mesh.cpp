#include "mesh.h"
#include "dataset.h"
#include "material_point.h"
#include "shape_function.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void GaussianDistribution(array<array<double, 3>, 6> &dec2p) {

    dec2p = {};

    for (int i = 0; i < 3; i++) {
        switch (npxye[i]) {
        case 1:
            dec2p[0][i] = 0.0e0;
            break;
        case 2:
            dec2p[0][i] = -0.25e0 * dxy[i];
            dec2p[1][i] = 0.25e0 * dxy[i];
            break;
        case 3:
            dec2p[0][i] = -1.0e0 / 3.0e0 * dxy[i];
            dec2p[1][i] = 0.0e0;
            dec2p[2][i] = 1.0e0 / 3.0e0 * dxy[i];
            break;
        case 4:
            dec2p[0][i] = -0.375e0 * dxy[i];
            dec2p[1][i] = -0.125e0 * dxy[i];
            dec2p[2][i] = 0.125e0 * dxy[i];
            dec2p[3][i] = 0.375e0 * dxy[i];
            break;
        case 5:
            dec2p[0][i] = -0.4e0 * dxy[i];
            dec2p[1][i] = -0.2e0 * dxy[i];
            dec2p[2][i] = 0.0e0;
            dec2p[3][i] = 0.2e0 * dxy[i];
            dec2p[4][i] = 0.4e0 * dxy[i];
            break;
        case 6:
            dec2p[0][i] = -5.0 / 12.0 * dxy[i];
            dec2p[1][i] = -3.0 / 12.0 * dxy[i]; // = -1.0/4.0
            dec2p[2][i] = -1.0 / 12.0 * dxy[i];
            dec2p[3][i] = 1.0 / 12.0 * dxy[i];
            dec2p[4][i] = 3.0 / 12.0 * dxy[i]; // = 1.0/4.0
            dec2p[5][i] = 5.0 / 12.0 * dxy[i];
            break;
        }
    }
    return;
}

void BuildMesh() {

    double xmin = xymin[0];
    double ymin = xymin[1];
    double zmin = xymin[2];
    int xelem = xyelem[0];
    int yelem = xyelem[1];
    int zelem = xyelem[2];
    double dx = dxy[0];
    double dy = dxy[1];
    double dz = dxy[2];

    int xnode = xynode[0];
    int ynode = xynode[1];
    int znode = xynode[2];

    xyn.assign(node, array<double, 3>{});
    for (int k = 0; k < znode; k++) {
        for (int j = 0; j < ynode; j++) {
            for (int i = 0; i < xnode; i++) {
                int id = i + xnode * j + xnode * ynode * k;
                xyn[id][0] = dx * i + xmin;
                xyn[id][1] = dy * j + ymin;
                xyn[id][2] = dz * k + zmin;
            }
        }
    }

    nc.assign(nelem, array<int, 8>{});
    for (int k = 0; k < zelem; k++) {
        for (int j = 0; j < yelem; j++) {
            for (int i = 0; i < xelem; i++) {
                int mid = i + xelem * j + xelem * yelem * k;
                nc[mid][0] = xnode * ynode * k + xnode * j + i;
                nc[mid][1] = xnode * ynode * k + xnode * j + i + 1;
                nc[mid][2] = xnode * ynode * k + xnode * (j + 1) + i + 1;
                nc[mid][3] = xnode * ynode * k + xnode * (j + 1) + i;
                nc[mid][4] = xnode * ynode * (k + 1) + xnode * j + i;
                nc[mid][5] = xnode * ynode * (k + 1) + xnode * j + i + 1;
                nc[mid][6] = xnode * ynode * (k + 1) + xnode * (j + 1) + i + 1;
                nc[mid][7] = xnode * ynode * (k + 1) + xnode * (j + 1) + i;
            }
        }
    }

    return;
}

// ---- Build control point coordinate ----
void BuildControlPoint() {

    int nxc = xynodec[0];
    int nyc = xynodec[1];
    int nzc = xynodec[2];
    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];

    int nenode = (idimc[0] + 1) * (idimc[1] + 1) * (idimc[2] + 1);
    ncc.assign(nelem, vector<int>(nenode));
    for (int m = 0; m < nelem; m++) {

        int ize = m / (nex * ney);
        int iye = (m - ize * (nex * ney)) / nex;
        int ixe = m - ize * (nex * ney) - iye * nex;

        int ixyn[3];
        int id = 0;
        for (int k = 0; k <= idimc[2]; k++) {
            ixyn[2] = ize + k;
            for (int j = 0; j <= idimc[1]; j++) {
                ixyn[1] = iye + j;
                for (int i = 0; i <= idimc[0]; i++) {
                    ixyn[0] = ixe + i;
                    ncc[m][id] = ixyn[0] + nxc * ixyn[1] + nxc * nyc * ixyn[2];
                    id++;
                }
            }
        }
    }

    xyc.assign(nodec, array<double, 3>{});
    for (int k = aelemmin[2]; k <= aelemmax[2] + idimc[2]; k++) {
        for (int j = aelemmin[1]; j <= aelemmax[1] + idimc[1]; j++) {
            for (int i = aelemmin[0]; i <= aelemmax[0] + idimc[0]; i++) {
                int kl = k - aelemmin[2];
                int jl = j - aelemmin[1];
                int il = i - aelemmin[0];
                int nid = il + xynodec[0] * jl + xynodec[0] * xynodec[1] * kl;
                std::vector<int> ijk_l = {il, jl, kl};
                std::vector<int> ijk_g = {i, j, k};
                // ---- X direction ----
                DefineCpPos<0>(nid, ijk_l, ijk_g);
                // ---- Y direction ----
                DefineCpPos<1>(nid, ijk_l, ijk_g);
                // ---- Z direction ----
                DefineCpPos<2>(nid, ijk_l, ijk_g);
            }
        }
    }

    return;
}

void MakNodalVol() {
    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);
    int nex = xyelem[0], ney = xyelem[1], nez = xyelem[2];

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, nvol);
    for (int m = 0; m < nelem; m++) {
        int ize = m / (nex * ney);
        int iye = (m - ize * (nex * ney)) / nex;
        int ixe = m - ize * (nex * ney) - iye * nex;

        std::array<double, 3> xye, xyp;
        xye[0] = xymin[0] + dxy[0] * (double(ixe) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(iye) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ize) + 0.5e0);
        for (int iz = 0; iz < npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyp[0] = xye[0] + dec2p[ix][0];
                    MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        nvol[nid] += sfi * volp;
                    }
                }
            }
        }
    }

    NodeVarComm(nvol, 0);

    return;
}
