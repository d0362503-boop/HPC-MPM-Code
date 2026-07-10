#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/fluid/MPM/stabilized_mpm.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace stabilizedmpm;

void StabilizedMPM::DataInput() {

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
    infile >> dt >> mtol >> this->gamma_nb >> this->beta_nb;
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
    infile >> this->rho >> this->rmu;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');
    infile.close();

    InitializeMeshAndTimeParameters();

    // --- Generalized α & Newmark β parameter Init ---
    // this->GeneralizedAlphaParaSet();
    this->NewmarkBetaParaSet();
    // --------------------------------

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile = OpenInputFile(filename);

    InputParaGriddata(infile);
    // --- Input Fluid BC ---
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
