#include "data_io.h"
#include "dataset.h"
#include "material_point.h"
#include "mesh.h"
#include "mpi_data.h"
#include "shape_function.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

ifstream OpenInputFile(const string &filename) {

    ifstream infile(filename);
    if (!infile.is_open()) { throw runtime_error("Failed to open input file: " + filename); }
    return infile;
}

ofstream OpenOutputFile(const string &filename) {

    ofstream outfile(filename);
    if (!outfile.is_open()) { throw runtime_error("Failed to open output file: " + filename); }
    outfile.flags(ios::right | ios::scientific);
    return outfile;
}

void InputParaGriddata(ifstream &infile) {

    infile >> myrank;
    infile.ignore(1000, '\n');

    infile >> nelem;
    for (int i = 0; i < 3; i++) infile >> xyelem[i];
    infile.ignore(1000, '\n');

    infile >> node;
    for (int i = 0; i < 3; i++) infile >> xynode[i];
    infile.ignore(1000, '\n');

    infile >> nodec;
    for (int i = 0; i < 3; i++) infile >> xynodec[i];
    infile.ignore(1000, '\n');

    nu = 0;
    nv = nu + node;
    nw = nv + node;
    np = nw + node;

    nuc = 0;
    nvc = nuc + nodec;
    nwc = nvc + nodec;
    npc = nwc + nodec;

    for (int i = 0; i < 3; i++) infile >> nxyr[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> xymin[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> xymax[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> aelemmin[i];
    infile.ignore(1000, '\n');

    for (int i = 0; i < 3; i++) infile >> aelemmax[i];
    infile.ignore(1000, '\n');

    infile >> isb; // --- Overlapped rank number ---
    infile.ignore(1000, '\n');
    if (isb != 0) {
        naid.resize(isb);
        InputVector(infile, isb, naid); // --- Overlapped rank ID ---

        infile >> isubc; // --- Overlapped control point total number ---
        infile.ignore(1000, '\n');
        nsbc.resize(isb);
        InputVector(infile, isb, nsbc); // --- Overlapped control point number with each rank ---

        nsubc.resize(isubc);
        InputVector(infile, isubc, nsubc); // --- Overlapped control point ID ---

        vector<int> dbn(nodec);
        InputVector(infile, nodec, dbn); // --- Number of control point shares ---

        dbc.resize(nodec * 4);
        for (int i = 0; i < nodec; ++i) {
            const double val = 1.0e0 / double(dbn[i]);
            for (int d = 0; d < 4; ++d) { dbc[i + d * nodec] = val; }
        }

        infile >> isubl; // --- Overlapped node total number ---
        infile.ignore(1000, '\n');
        nsbl.resize(isb);
        InputVector(infile, isb, nsbl); // --- Overlapped node number with each rank ---

        nsubl.resize(isubl);
        InputVector(infile, isubl, nsubl); // --- Overlapped node ID ---

        dbn.assign(node, 0.0e0);
        InputVector(infile, node, dbn); // --- Number of node shares ---

        dbl.resize(node * 4);
        for (int i = 0; i < node; i++) {
            dbl[i + nu] = 1.0e0 / double(dbn[i]);
            dbl[i + nv] = 1.0e0 / double(dbn[i]);
            dbl[i + nw] = 1.0e0 / double(dbn[i]);
            dbl[i + np] = 1.0e0 / double(dbn[i]);
        }
    } else {
        naid.resize(0);
        nsbc.resize(0);
        nsubc.resize(0);
        dbc.assign(nodec * 4, 1.0e0);
        nsbl.resize(0);
        nsubl.resize(0);
        dbl.assign(node * 4, 1.0e0);
    }

    return;
}
