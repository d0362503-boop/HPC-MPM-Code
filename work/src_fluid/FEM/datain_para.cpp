#include "module/bc.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/FEM/stabilized_fem.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace stabilizedfem;

void StabilizedFEM::DataInput() {

    std::ifstream infile = OpenInputFile("file.dat");
    getline(infile, parafile);
    getline(infile, gridfile);
    getline(infile, pointfile);
    getline(infile, outfile);
    infile.close();

    infile = OpenInputFile(parafile);
    infile.ignore(1000, '\n');
    infile >> rstflag >> nlstep;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> this->integrator_.spec_rad;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xyminw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xymaxw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xyelemw[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->rhol >> this->rmul >> this->rhog >> this->rmug >> this->fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');

    infile.close();

    InitializeMeshAndTimeParameters(g_params, g_grid);

    // --- Generalized α parameter Init ---
    this->integrator_.GeneralizedAlphaParaSet();
    // --------------------------------

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile, g_grid);
    // --- Input Fluid BC ---
    this->InputBCData(infile);

    infile.close();

    if (rstflag == 0 || rstflag == 2) { this->InitializeMeshData(); }

    return;
}
