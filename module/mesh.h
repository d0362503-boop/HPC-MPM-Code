#pragma once

#include "module/dataset.h"
#include "module/mpi_data.h"
#include <array>
#include <vector>

/// Background (visualization) mesh: global sizes plus this rank's element window.
struct Mesh {
    int node = 0, nodew = 0, nelem = 0, nelemw = 0; // rank / global node and element counts
    int nu = 0, nv = 0, nw = 0, np = 0;               // node-field component offsets
    std::array<int, 3> xynode{}, xyelem{};            // rank nodes / elements per direction
    std::array<int, 3> aelemmin{}, aelemmax{};        // rank global element window
    std::array<int, 3> xynodew{}, xyelemw{};          // global nodes / elements per direction
    std::array<double, 3> dxy{};                      // element size
    std::array<double, 3> xymin{}, xymax{};           // rank bounding box
    std::array<double, 3> xyminw{}, xymaxw{};         // global bounding box
    std::vector<std::array<int, 8>> nc;               // element-node connectivity
    std::vector<std::array<double, 3>> xyn;           // node coordinates
};

/// B-spline control points of this rank (the computational DOF layout).
struct ControlMesh {
    int nodec = 0, nodecw = 0;                  // rank / global control points
    int nuc = 0, nvc = 0, nwc = 0, npc = 0;     // control-field component offsets
    std::array<int, 3> xynodec{}, xynodecw{};   // rank / global points per direction
    std::array<int, 3> idimc{}, idiml{};        // B-spline order per direction
    std::vector<double> nvol;                   // control-point volume
    std::vector<std::vector<int>> ncc;          // element-control-point connectivity
    std::vector<std::array<double, 3>> xyc;     // control-point coordinates
};

/// Complete computational grid of this rank: background mesh, control points and partition.
struct Grid {
    Mesh mesh;          // background mesh
    ControlMesh cmesh;  // control points
    Partition part;     // rank overlap bookkeeping
};

inline Grid g_grid; // transitional global instance (still read by CrsMat)

// ----- Transitional aliases still read by CrsMat; remove together with g_grid when CrsMat takes the grid -----
inline int &nodec = g_grid.cmesh.nodec;
inline int &nuc = g_grid.cmesh.nuc, &nvc = g_grid.cmesh.nvc, &nwc = g_grid.cmesh.nwc, &npc = g_grid.cmesh.npc;
inline std::array<int, 3> &xynodec = g_grid.cmesh.xynodec, &xynodecw = g_grid.cmesh.xynodecw;
inline std::array<int, 3> &idimc = g_grid.cmesh.idimc;
inline std::array<int, 3> &aelemmin = g_grid.mesh.aelemmin, &aelemmax = g_grid.mesh.aelemmax;
inline std::vector<double> &dbc = g_grid.part.dbc;

/// Transitional: NodeVarComm on the global partition.
template <typename T> void NodeVarComm(std::vector<T> &dat, int offset) { NodeVarComm(g_grid.part, dat, offset); }
/// Transitional: NodeVarComm on the global partition.
template <typename T> void NodeVarComm(std::vector<T> &dat, const std::vector<int> &offsets) {
    NodeVarComm(g_grid.part, dat, offsets);
}
// -----------------------------------------------------------------------------------

/**
 * @brief Initialize mesh spacings and time-integration scalars derived from the input parameters.
 *
 * Computes `dti` and `dlstep` in `params`, and `dxy`, `xynodew`, `xynodecw` in `grid`
 * from the already-read global box, element counts and B-spline orders.
 * @param params Run parameters; `dt` and `nlstep` must be set, `dti` and `dlstep` are written.
 * @param grid   Grid whose global box, element counts and `idimc` are set; spacings are written.
 */
inline void InitializeMeshAndTimeParameters(RunParams &params, Grid &grid) {

    Mesh &mesh = grid.mesh;
    ControlMesh &cm = grid.cmesh;
    params.dti = 1.0e0 / params.dt;
    for (int i = 0; i < 3; i++) {
        mesh.dxy[i] = (mesh.xymaxw[i] - mesh.xyminw[i]) / double(mesh.xyelemw[i]);
        mesh.xynodew[i] = mesh.xyelemw[i] + 1;
        cm.xynodecw[i] = mesh.xyelemw[i] + cm.idimc[i];
    }
    params.dlstep = 1.0e0 / double(params.nlstep);

    return;
}

/**
 * @brief Decompose a linear background element index into per-direction element indices.
 * @param m  Linear element index (x varies fastest, then y, then z).
 * @param ne Number of background elements along x/y/z (typically the rank `xyelem`).
 * @return Element indices along x/y/z.
 */
inline std::array<int, 3> IndexToIJK(int m, const std::array<int, 3> &ne) {

    const int ize = m / (ne[0] * ne[1]);
    const int iye = (m - ize * (ne[0] * ne[1])) / ne[0];
    const int ixe = m - ize * (ne[0] * ne[1]) - iye * ne[0];

    return {ixe, iye, ize};
}

/**
 * @brief Flatten per-direction background element indices into a linear index.
 * @param ijk Element indices along x/y/z.
 * @param ne Number of background elements along x/y/z.
 * @return Linear element index with x varying fastest, then y, then z.
 */
inline int IJKToIndex(const std::array<int, 3> &ijk, const std::array<int, 3> &ne) {
    return ijk[0] + ne[0] * (ijk[1] + ne[1] * ijk[2]);
}

/**
 * @brief Compute the center coordinate of a rank-local background element.
 * @param mesh Background mesh of this rank.
 * @param m Local element index on this rank (x varies fastest, then y, then z).
 * @return Physical coordinates of the element center.
 */
inline std::array<double, 3> ElementCenter(const Mesh &mesh, int m) {

    const std::array<int, 3> ijk = IndexToIJK(m, mesh.xyelem);
    std::array<double, 3> xye;
    for (int i = 0; i < 3; i++) { xye[i] = mesh.xymin[i] + mesh.dxy[i] * (double(ijk[i]) + 0.5e0); }

    return xye;
}

/**
 * @brief Compute sub-cell center offsets used to place points inside one background element.
 *
 * Each element is split into `npxye[i]` equal sub-cells along direction i, and
 * `dec2p[k][i]` is the offset of the k-th sub-cell center from the element center:
 * ((k + 0.5) / npxye[i] - 0.5) * dxy[i]. For example, 2 points give -dxy/4 and +dxy/4.
 * These are uniformly spaced points, not Gauss quadrature points, and no weights are computed.
 * Supports 1 to 6 points per direction; unused entries are zero.
 * Used for initial particle positions, nodal volumes, and FSI interface area weights.
 *
 * @param npxye Number of sub-cell points per element along x/y/z.
 * @param dxy   Element size along x/y/z.
 * @param dec2p Output offsets; add `dec2p[k][i]` to the element center to get the k-th point along direction i.
 */
void SubcellOffsets(const std::array<int, 3> &npxye, const std::array<double, 3> &dxy,
                    std::array<std::array<double, 3>, 6> &dec2p);

/**
 * @brief Locate the global background element indices containing point `xq`.
 * @param mesh Background mesh (global box, element size and counts are used).
 * @param xq   Query point in physical coordinates.
 * @param ie_g Output global element indices along x/y/z.
 * @return `true` if the point lies inside the global background mesh bounds.
 */
bool LocateGlobalElement(const Mesh &mesh, const std::array<double, 3> &xq, std::array<int, 3> &ie_g);

/**
 * @brief Convert global background element indices to the current rank-local element index.
 * @param mesh    Background mesh of this rank (its global element window is used).
 * @param ie_g    Global element indices along x/y/z.
 * @param m_local Output local element index on this rank.
 * @return `true` if the element belongs to the current rank partition.
 */
bool GlobalElementToLocal(const Mesh &mesh, const std::array<int, 3> &ie_g, int &m_local);

/**
 * @brief Locate the current rank-local background element containing point `xq`.
 * @param mesh    Background mesh of this rank.
 * @param xq      Query point in physical coordinates.
 * @param m_local Output local element index on this rank.
 * @return `true` if the point is inside the current rank partition.
 */
bool LocateLocalElement(const Mesh &mesh, const std::array<double, 3> &xq, int &m_local);

/**
 * @brief Compute control-point volumes from sub-cell point volumes, summed over overlap ranks.
 * @param params Run parameters (`npxye` sets the sub-cell points per element).
 * @param grid   Grid of this rank; `grid.cmesh.nvol` is written.
 */
void ComputeNodalVol(const RunParams &params, Grid &grid);

/**
 * @brief Build the background mesh node coordinates and element connectivity.
 * @param mesh Background mesh whose sizes and box are set; `xyn` and `nc` are written.
 */
void BuildMesh(Mesh &mesh);

/**
 * @brief Build the B-spline control-point coordinates and element connectivity.
 * @param mesh Background mesh of this rank (box, spacing and element window).
 * @param cm   Control points whose sizes and orders are set; `xyc` and `ncc` are written.
 */
void BuildControlPoint(const Mesh &mesh, ControlMesh &cm);

/**
 * @brief Compute the coordinate of control point `nid` along direction `dir`.
 * @tparam dir  Spatial direction (0=x, 1=y, 2=z).
 * @param mesh  Background mesh of this rank.
 * @param cm    Control points; `cm.xyc[nid][dir]` is written.
 * @param nid   Control-point index.
 * @param ijk_l Local IJK indices of the control point.
 * @param ijk_g Global IJK indices of the control point.
 */
template <int dir>
void DefineCpPos(const Mesh &mesh, ControlMesh &cm, int nid, const std::array<int, 3> &ijk_l,
                 const std::array<int, 3> &ijk_g) {

    if (cm.idimc[dir] == 1) {
        cm.xyc[nid][dir] = mesh.xymin[dir] + mesh.dxy[dir] * double(ijk_l[dir]);
    } else if (cm.idimc[dir] == 2) {
        if (ijk_g[dir] == 0) {
            cm.xyc[nid][dir] = mesh.xyminw[dir];
        } else if (ijk_g[dir] == cm.xynodecw[dir] - 1) {
            cm.xyc[nid][dir] = mesh.xymaxw[dir];
        } else {
            cm.xyc[nid][dir] = mesh.xymin[dir] + mesh.dxy[dir] * (double(ijk_l[dir]) - 0.5e0);
        }
    }

    return;
}
