#include "module/bc.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace implicitmpm;

void ImplicitSolidMPM::DataInput(Grid &grid, RunParams &params) {

    std::ifstream infile = OpenInputFile("file.dat");
    getline(infile, params.parafile);
    getline(infile, params.gridfile);
    getline(infile, params.pointfile);
    getline(infile, params.outfile);
    infile.close();

    infile = OpenInputFile(params.parafile);
    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> params.rstflag >> params.nlstep;
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
    infile >> this->rho >> this->materials_.nmat >> this->materials_.npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(this->materials_.nmat, this->materials_.iprop);
    VectorAssign(this->materials_.nmat, this->materials_.nprop);
    this->materials_.mat_prop.assign(this->materials_.nmat, std::vector<double>(this->materials_.npropmax));
    for (int n = 0; n < this->materials_.nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> this->materials_.iprop[n] >> this->materials_.nprop[n];
        infile.ignore(1000, '\n');
        for (int i = 0; i < this->materials_.nprop[n]; i++) { infile >> this->materials_.mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

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
    // --- solid boundary condition ---
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
