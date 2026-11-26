#include "module/data_io.h"
#include "module/DLB/mpm_dlb.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include <fstream>
#include <stdexcept>
#include <vector>

std::ifstream OpenInputFile(const std::string &filename) {

    std::ifstream infile(filename);
    if (!infile.is_open()) { throw std::runtime_error("Failed to open input file: " + filename); }
    return infile;
}

std::ofstream OpenOutputFile(const std::string &filename) {

    std::ofstream outfile(filename);
    if (!outfile.is_open()) { throw std::runtime_error("Failed to open output file: " + filename); }
    outfile.flags(std::ios::right | std::ios::scientific);
    return outfile;
}

void InputParaGriddata(std::ifstream &infile, Grid &grid) {

    infile >> myrank;
    infile.ignore(1000, '\n');

    infile >> grid.mesh.nelem;
    for (int i = 0; i < 3; i++) infile >> grid.mesh.xyelem[i];
    infile.ignore(1000, '\n');

    infile >> grid.mesh.node;
    for (int i = 0; i < 3; i++) infile >> grid.mesh.xynode[i];
    infile.ignore(1000, '\n');

    infile >> grid.cmesh.nodec;
    for (int i = 0; i < 3; i++) infile >> grid.cmesh.xynodec[i];
    infile.ignore(1000, '\n');

    grid.mesh.nu = 0;
    grid.mesh.nv = grid.mesh.nu + grid.mesh.node;
    grid.mesh.nw = grid.mesh.nv + grid.mesh.node;
    grid.mesh.np = grid.mesh.nw + grid.mesh.node;

    grid.cmesh.nuc = 0;
    grid.cmesh.nvc = grid.cmesh.nuc + grid.cmesh.nodec;
    grid.cmesh.nwc = grid.cmesh.nvc + grid.cmesh.nodec;
    grid.cmesh.npc = grid.cmesh.nwc + grid.cmesh.nodec;

    for (int i = 0; i < 3; i++) infile >> grid.part.nxyr[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> grid.mesh.xymin[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> grid.mesh.xymax[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> grid.mesh.aelemmin[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> grid.mesh.aelemmax[i];
    infile.ignore(1000, '\n');

    infile >> grid.part.isb; // --- Overlapped rank number ---
    infile.ignore(1000, '\n');
    if (grid.part.isb != 0) {
        grid.part.naid.resize(grid.part.isb);
        InputVector(infile, grid.part.isb, grid.part.naid); // --- Overlapped rank ID ---

        infile >> grid.part.isubc; // --- Overlapped control point total number ---
        infile.ignore(1000, '\n');
        grid.part.nsbc.resize(grid.part.isb);
        InputVector(infile, grid.part.isb, grid.part.nsbc); // --- Overlapped control point number with each rank ---

        grid.part.nsubc.resize(grid.part.isubc);
        InputVector(infile, grid.part.isubc, grid.part.nsubc); // --- Overlapped control point ID ---

        std::vector<int> dbn(grid.cmesh.nodec);
        InputVector(infile, grid.cmesh.nodec, dbn); // --- Number of control point shares ---

        grid.part.dbc.resize(grid.cmesh.nodec * 4);
        for (int i = 0; i < grid.cmesh.nodec; ++i) {
            const double val = 1.0e0 / double(dbn[i]);
            for (int d = 0; d < 4; ++d) { grid.part.dbc[i + d * grid.cmesh.nodec] = val; }
        }

        infile >> grid.part.isubl; // --- Overlapped node total number ---
        infile.ignore(1000, '\n');
        grid.part.nsbl.resize(grid.part.isb);
        InputVector(infile, grid.part.isb, grid.part.nsbl); // --- Overlapped node number with each rank ---

        grid.part.nsubl.resize(grid.part.isubl);
        InputVector(infile, grid.part.isubl, grid.part.nsubl); // --- Overlapped node ID ---

        dbn.assign(grid.mesh.node, 0.0e0);
        InputVector(infile, grid.mesh.node, dbn); // --- Number of node shares ---

        grid.part.dbl.resize(grid.mesh.node * 4);
        for (int i = 0; i < grid.mesh.node; i++) {
            grid.part.dbl[i + grid.mesh.nu] = 1.0e0 / double(dbn[i]);
            grid.part.dbl[i + grid.mesh.nv] = 1.0e0 / double(dbn[i]);
            grid.part.dbl[i + grid.mesh.nw] = 1.0e0 / double(dbn[i]);
            grid.part.dbl[i + grid.mesh.np] = 1.0e0 / double(dbn[i]);
        }
    } else {
        grid.part.naid.resize(0);
        grid.part.nsbc.resize(0);
        grid.part.nsubc.resize(0);
        grid.part.dbc.assign(grid.cmesh.nodec * 4, 1.0e0);
        grid.part.nsbl.resize(0);
        grid.part.nsubl.resize(0);
        grid.part.dbl.assign(grid.mesh.node * 4, 1.0e0);
    }

    mpm_dlb::CollectCurrentRegions(grid.mesh, grid.part);

    return;
}
