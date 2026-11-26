#include "module/bc.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace stabilizedmpm;

void StabilizedMPM::DataInput(Grid &grid, RunParams &params) {

    std::ifstream infile = OpenInputFile("file.dat");
    getline(infile, params.parafile);
    getline(infile, params.gridfile);
    getline(infile, params.pointfile);
    getline(infile, params.outfile);
    infile.close();

    infile = OpenInputFile(params.parafile);
    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> this->do_dlb >> params.rstflag >> params.nlstep;
    this->solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.ista >> params.iend >> params.iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.dt >> params.mtol >> this->implicit_.integrator_.spec_rad;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> grid.mesh.xyminw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> grid.mesh.xymaxw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> grid.mesh.xyelemw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> grid.cmesh.idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->rho >> this->rmu;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.bb[i];
    infile.ignore(1000, '\n');

    infile.close();

    InitializeMeshAndTimeParameters(params, grid);

    // --- Generalized α & Newmark β parameter Init ---
    this->implicit_.integrator_.GeneralizedAlphaParaSet();
    this->implicit_.integrator_.NewmarkBetaParaSet(params.dt);
    // --------------------------------

    std::string filename = params.gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile, grid);
    // --- Input Fluid BC ---
    this->InputBCData(infile);

    infile.close();

    if (params.rstflag == 0 || params.rstflag == 2) {
        filename = params.pointfile + std::to_string(myrank) + ".txt";
        infile = OpenInputFile(filename);

        this->InputPointData(infile);

        infile.close();
    }

    return;
}
