#include "module/mesh.h"
#include "module/dataset.h"
#include "module/shape_function.h"
#include <algorithm>
#include <vector>

void SubcellOffsets(const std::array<int, 3> &npxye, const std::array<double, 3> &dxy,
                    std::array<std::array<double, 3>, 6> &dec2p) {

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

bool LocateGlobalElement(const Mesh &mesh, const std::array<double, 3> &xq, std::array<int, 3> &ie_g) {

    for (int d = 0; d < 3; d++) {
        if (xq[d] < mesh.xyminw[d] || xq[d] > mesh.xymaxw[d]) { return false; }

        double xi = (xq[d] - mesh.xyminw[d]) / mesh.dxy[d];
        int ie = static_cast<int>(std::floor(xi));
        ie_g[d] = std::clamp(ie, 0, mesh.xyelemw[d] - 1);
    }

    return true;
}

bool GlobalElementToLocal(const Mesh &mesh, const std::array<int, 3> &ie_g, int &m_local) {

    for (int d = 0; d < 3; d++) {
        if (ie_g[d] < mesh.aelemmin[d] || ie_g[d] > mesh.aelemmax[d]) { return false; }
    }

    int iex_l = ie_g[0] - mesh.aelemmin[0];
    int iey_l = ie_g[1] - mesh.aelemmin[1];
    int iez_l = ie_g[2] - mesh.aelemmin[2];

    m_local = iex_l + mesh.xyelem[0] * iey_l + mesh.xyelem[0] * mesh.xyelem[1] * iez_l;

    return true;
}

bool LocateLocalElement(const Mesh &mesh, const std::array<double, 3> &xq, int &m_local) {

    std::array<int, 3> ie_g{};
    if (!LocateGlobalElement(mesh, xq, ie_g)) { return false; }

    return GlobalElementToLocal(mesh, ie_g, m_local);
}

void BuildMesh(Mesh &mesh) {

    double xmin = mesh.xymin[0];
    double ymin = mesh.xymin[1];
    double zmin = mesh.xymin[2];
    int xelem = mesh.xyelem[0];
    int yelem = mesh.xyelem[1];
    int zelem = mesh.xyelem[2];
    double dx = mesh.dxy[0];
    double dy = mesh.dxy[1];
    double dz = mesh.dxy[2];

    int xnode = mesh.xynode[0];
    int ynode = mesh.xynode[1];
    int znode = mesh.xynode[2];

    std::vector<std::array<double, 3>> &xyn = mesh.xyn;
    std::vector<std::array<int, 8>> &nc = mesh.nc;
    VectorAssign(mesh.node, xyn);
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

    VectorAssign(mesh.nelem, nc);
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
void BuildControlPoint(const Mesh &mesh, ControlMesh &cm) {

    const std::array<int, 3> &idimc = cm.idimc;
    const std::array<int, 3> &aelemmin = mesh.aelemmin;
    const std::array<int, 3> &aelemmax = mesh.aelemmax;
    std::vector<std::vector<int>> &ncc = cm.ncc;
    int nxc = cm.xynodec[0];
    int nyc = cm.xynodec[1];
    int nzc = cm.xynodec[2];

    int nenode = (idimc[0] + 1) * (idimc[1] + 1) * (idimc[2] + 1);
    ncc.assign(mesh.nelem, std::vector<int>(nenode));
    for (int m = 0; m < mesh.nelem; m++) {

        const std::array<int, 3> ijk = IndexToIJK(m, mesh.xyelem);

        int ixyn[3];
        int id = 0;
        for (int k = 0; k <= idimc[2]; k++) {
            ixyn[2] = ijk[2] + k;
            for (int j = 0; j <= idimc[1]; j++) {
                ixyn[1] = ijk[1] + j;
                for (int i = 0; i <= idimc[0]; i++) {
                    ixyn[0] = ijk[0] + i;
                    ncc[m][id] = ixyn[0] + nxc * ixyn[1] + nxc * nyc * ixyn[2];
                    id++;
                }
            }
        }
    }

    cm.xyc.assign(cm.nodec, std::array<double, 3>{});
    for (int k = aelemmin[2]; k <= aelemmax[2] + idimc[2]; k++) {
        for (int j = aelemmin[1]; j <= aelemmax[1] + idimc[1]; j++) {
            for (int i = aelemmin[0]; i <= aelemmax[0] + idimc[0]; i++) {
                int kl = k - aelemmin[2];
                int jl = j - aelemmin[1];
                int il = i - aelemmin[0];
                int nid = il + cm.xynodec[0] * jl + cm.xynodec[0] * cm.xynodec[1] * kl;
                std::array<int, 3> ijk_l = {il, jl, kl};
                std::array<int, 3> ijk_g = {i, j, k};
                // ---- X direction ----
                DefineCpPos<0>(mesh, cm, nid, ijk_l, ijk_g);
                // ---- Y direction ----
                DefineCpPos<1>(mesh, cm, nid, ijk_l, ijk_g);
                // ---- Z direction ----
                DefineCpPos<2>(mesh, cm, nid, ijk_l, ijk_g);
            }
        }
    }

    return;
}

void ComputeNodalVol(const RunParams &params, Grid &grid) {

    const Mesh &mesh = grid.mesh;
    ControlMesh &cm = grid.cmesh;
    const std::array<int, 3> &npxye = params.npxye;
    const std::array<double, 3> &dxy = mesh.dxy;

    std::array<std::array<double, 3>, 6> dec2p;
    SubcellOffsets(npxye, dxy, dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(cm.nodec, cm.nvol);
    for (int m = 0; m < mesh.nelem; m++) {
        const std::array<double, 3> xye = ElementCenter(mesh, m);
        std::array<double, 3> xyp;
        for (int iz = 0; iz < npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyp[0] = xye[0] + dec2p[ix][0];
                    MakeSF(grid, m, xyp, ncm, nenode, sf, dsf);
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        cm.nvol[nid] += sfi * volp;
                    }
                }
            }
        }
    }

    NodeVarComm(grid.part, cm.nvol, 0);

    return;
}
