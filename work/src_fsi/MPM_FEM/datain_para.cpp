#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "module/bc.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solid/solid_material_point.h"
#include "work/src_fsi/MPM_FEM/block_fsi.h"

using namespace implicitmpm;
using namespace stabilizedfem;

using namespace mpmfemblockfsi;

void MPMFEMBlockFSI::DataInput(Grid &grid, RunParams &params) {

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
    this->solid_.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.ista >> params.iend >> params.iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.dt >> params.mtol >> this->fluid_.integrator_.spec_rad >> this->solid_.implicit_.integrator_.spec_rad;
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
    infile >> this->solid_.rho >> this->fluid_.rhol >> this->fluid_.rmul //
        >> this->fluid_.rhog >> this->fluid_.rmug >> this->fluid_.fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->solid_.materials_.nmat >> this->solid_.materials_.npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(this->solid_.materials_.nmat, this->solid_.materials_.iprop);
    VectorAssign(this->solid_.materials_.nmat, this->solid_.materials_.nprop);
    this->solid_.materials_.mat_prop.assign(this->solid_.materials_.nmat, std::vector<double>(this->solid_.materials_.npropmax));
    for (int n = 0; n < this->solid_.materials_.nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> this->solid_.materials_.iprop[n] >> this->solid_.materials_.nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < this->solid_.materials_.nprop[n]; i++) { infile >> this->solid_.materials_.mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.bb[i];
    infile.ignore(1000, '\n');
    infile.close();

    InitializeMeshAndTimeParameters(params, grid);

    // --- Generalized α & Newmark β parameter Init ---
    this->solid_.implicit_.integrator_.GeneralizedAlphaParaSet();
    this->solid_.implicit_.integrator_.NewmarkBetaParaSet(params.dt);
    this->fluid_.integrator_.GeneralizedAlphaParaSet();
    // --------------------------------

    std::string filename = params.gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile, grid);
    // --- Solid boundary condition ---
    this->solid_.InputBCData(infile);
    // --- Fluid boundary condition ---
    this->fluid_.InputBCData(infile);

    infile.close();

    if (params.rstflag == 0 || params.rstflag == 2) {
        filename = params.pointfile + std::to_string(myrank) + ".txt";
        infile = OpenInputFile(filename);

        this->solid_.InputPointData(infile);

        infile.close();

        this->fluid_.InitializeMeshData();
    }

    return;
}
