#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../module/bc.h"
#include "../module/data_io.h"
#include "../module/dataset.h"
#include "../module/fluid/FEM/stabilized_fem.h"
#include "../module/material_point.h"
#include "../module/mesh.h"
#include "../module/mpi_data.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"
#include "../module/solid/solid_material_point.h"
#include "block_fsi.h"

using namespace implicitmpm;
using namespace stabilizedfem;

void BlockFSI::DataInput() {
    std::ifstream infile;
    infile.open("file.dat");
    getline(infile, parafile);
    getline(infile, gridfile);
    getline(infile, pointfile);
    getline(infile, outfile);
    infile.close();

    infile.open(parafile);
    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    this->solid_.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> this->fluid_.spec_rad >> this->solid_.spec_rad;
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
    infile >> this->solid_.rho >> this->fluid_.rhol >> this->fluid_.rmul  //
        >> this->fluid_.rhog >> this->fluid_.rmug >> this->fluid_.fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> nmat >> npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(nmat, iprop);
    VectorAssign(nmat, nprop);
    mat_prop.assign(nmat, std::vector<double>(npropmax));
    for (int n = 0; n < nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> iprop[n] >> nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < nprop[n]; i++) {
            infile >> mat_prop[n][i];
        }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');
    infile.close();

    dti = 1.0e0 / dt;
    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        xynodew[i] = xyelemw[i] + 1;
        xynodecw[i] = xyelemw[i] + idimc[i];
    }
    dlstep = 1.0e0 / double(nlstep);

    // --- Generalized α & Newmark β parameter Init ---
    this->solid_.GeneralizedAlphaParaSet();
    this->solid_.NewmarkBetaParaSet();
    this->fluid_.GeneralizedAlphaParaSet();
    // --------------------------------

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile.open(filename);

    InputParaGriddata(infile);

    // --- solid boundary condition ---
    this->solid_.ubc.BCInput(infile);
    this->solid_.vbc.BCInput(infile);
    this->solid_.wbc.BCInput(infile);
    // --- water boundary condition ---
    this->fluid_.ubc.BCInput(infile);
    this->fluid_.vbc.BCInput(infile);
    this->fluid_.wbc.BCInput(infile);
    this->fluid_.pbc.BCInput(infile);

    infile.close();

    if (rstflag == 0 || rstflag == 2) {
        filename = pointfile + std::to_string(myrank) + ".txt";
        infile.open(filename);

        this->solid_.InputPointData(infile);

        infile.close();

        this->fluid_.InitializeMeshData();
    }

    return;
}
