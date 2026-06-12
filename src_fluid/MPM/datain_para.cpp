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

    dti = 1.0e0 / dt;
    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        xynodew[i] = xyelemw[i] + 1;
        xynodecw[i] = xyelemw[i] + idimc[i];
    }
    dlstep = 1.0e0 / double(nlstep);

    // ----- Newmark beta parameter -----
    this->NewmarkBetaParaSet();

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile.open(filename);

    InputParaGriddata(infile);

    this->ubc.BCInput(infile);
    this->vbc.BCInput(infile);
    this->wbc.BCInput(infile);
    this->pbc.BCInput(infile);
    // --- Inflow boundary ---
    this->uinfbc.BCInput(infile, false);
    this->vinfbc.BCInput(infile, false);
    this->winfbc.BCInput(infile, false);

    infile.close();

    if (rstflag == 0 || rstflag == 2) {
        filename = pointfile + std::to_string(myrank) + ".txt";
        infile.open(filename);
        infile >> this->num;
        infile.ignore(1000, '\n');

        this->InitializePointData();

        InputVector(infile, this->num, this->coord);

        for (int ip = 0; ip < this->num; ip++) {
            infile >> this->id[ip] >> this->matid[ip] >> this->mass[ip] >> this->vol[ip];
            infile.ignore(1000, '\n');
        }
        infile.close();
    }

    return;
}
