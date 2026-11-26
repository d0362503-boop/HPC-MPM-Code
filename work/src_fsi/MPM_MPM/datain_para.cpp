#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "module/bc.h"
#include "module/data_io.h"
#include "module/dataset.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/implicit/implicit_mpm_solid.h"
#include "module/solid/solid_material_point.h"
#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

using namespace implicitmpm;
using namespace stabilizedmpm;

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::DataInput() {

    std::ifstream infile = OpenInputFile("file.dat");
    getline(infile, parafile);
    getline(infile, gridfile);
    getline(infile, pointfile);
    getline(infile, outfile);
    infile.close();

    infile = OpenInputFile(parafile);
    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    infile.ignore(1000, '\n');
    this->solid_.solswitch = ParseMapScheme(solswitch_str);
    this->fluid_.solswitch = ParseMapScheme(solswitch_str);

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> this->fluid_.implicit_.integrator_.spec_rad >> this->solid_.implicit_.integrator_.spec_rad;
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
    infile >> this->solid_.rho >> this->fluid_.rho >> this->fluid_.rmu;
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
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');
    infile.close();

    InitializeMeshAndTimeParameters(g_params, g_grid);

    // --- Generalized α & Newmark β parameter Init ---
    this->solid_.implicit_.integrator_.GeneralizedAlphaParaSet();
    this->solid_.implicit_.integrator_.NewmarkBetaParaSet(dt);
    this->fluid_.implicit_.integrator_.GeneralizedAlphaParaSet();
    this->fluid_.implicit_.integrator_.NewmarkBetaParaSet(dt);
    // --------------------------------

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile, g_grid);
    // --- Solid boundary condition ---
    this->solid_.InputBCData(infile);
    // --- Fluid boundary condition ---
    this->fluid_.InputBCData(infile);

    infile.close();

    if (rstflag == 0 || rstflag == 2) {
        filename = pointfile + std::to_string(myrank) + ".txt";
        infile = OpenInputFile(filename);

        this->solid_.InputPointData(infile);

        this->fluid_.InputPointData(infile);

        infile.close();
    }

    return;
}
