#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include "../../module/solid/implicit/implicit_mpm_solid.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace implicitmpm;

void ImplicitSolidMPM::DataInput() {

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
    this->solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> this->spec_rad;
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
    infile >> this->rho >> nmat >> npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(nmat, iprop);
    VectorAssign(nmat, nprop);
    mat_prop.assign(nmat, std::vector<double>(npropmax));
    for (int n = 0; n < nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> iprop[n] >> nprop[n];
        infile.ignore(1000, '\n');
        for (int i = 0; i < nprop[n]; i++) { infile >> mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');
    infile.close();

    InitializeMeshAndTimeParameters();

    // --- Generalized α & Newmark β parameter Init ---
    this->GeneralizedAlphaParaSet();
    this->NewmarkBetaParaSet();
    // --------------------------------

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile);
    // --- solid boundary condition ---
    this->InputBCData(infile);

    infile.close();

    if (rstflag == 0 || rstflag == 2) {
        filename = pointfile + std::to_string(myrank) + ".txt";
        infile = OpenInputFile(filename);

        this->InputPointData(infile);

        infile.close();
    }

    return;
}
