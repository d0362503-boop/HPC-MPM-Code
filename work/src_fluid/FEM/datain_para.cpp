#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include <fstream>
#include <string>

using namespace stabilizedfem;

void StabilizedFEM::DataInput(Grid &grid, RunParams &params) {

    std::ifstream infile = OpenInputFile("file.dat");
    getline(infile, params.parafile);
    getline(infile, params.gridfile);
    getline(infile, params.pointfile);
    getline(infile, params.outfile);
    infile.close();

    infile = OpenInputFile(params.parafile);
    infile.ignore(1000, '\n');
    infile >> params.rstflag >> params.nlstep;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.ista >> params.iend >> params.iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.dt >> params.mtol >> this->integrator_.spec_rad;
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
    infile >> this->rhol >> this->rmul >> this->rhog >> this->rmug >> this->fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.bb[i];
    infile.ignore(1000, '\n');

    infile.close();

    InitializeMeshAndTimeParameters(params, grid);

    // --- Generalized α parameter Init ---
    this->integrator_.GeneralizedAlphaParaSet();
    // --------------------------------

    std::string filename = params.gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile, grid);
    // --- Input Fluid BC ---
    this->InputBCData(infile);

    infile.close();

    if (params.rstflag == 0 || params.rstflag == 2) { this->InitializeMeshData(); }

    return;
}
