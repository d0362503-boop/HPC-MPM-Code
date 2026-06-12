#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

void datain_para() {

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
    sp.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> sp.gamma_nb >> sp.beta_nb;
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
    infile >> sp.rho >> nmat >> npropmax;
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

    sp.NewmarkBetaParaSet();

    dti = 1.0e0 / dt;
    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymaxw[i] - xyminw[i]) / double(xyelemw[i]);
        xynodew[i] = xyelemw[i] + 1;
        xynodecw[i] = xyelemw[i] + idimc[i];
    }
    dlstep = 1.0e0 / double(nlstep);

    std::string filename = gridfile + std::to_string(myrank) + ".txt";
    infile.open(filename);

    InputParaGriddata(infile);

    BCInput(infile, usbc);
    BCInput(infile, vsbc);
    BCInput(infile, wsbc);

    infile.close();

    if (rstflag == 0 || rstflag == 2) {
        filename = pointfile + to_string(myrank) + ".txt";
        infile.open(filename);

        infile >> sp.num;
        infile.ignore(1000, '\n');

        solid_point_data_assign();

        input_twod_vector(infile, 3, sp.num, sp.coord);

        for (int ip = 0; ip < sp.num; ip++) {
            infile >> sp.id[ip] >> sp.matid[ip] >> sp.mass[ip] >> sp.vol0[ip];
            infile.ignore(1000, '\n');
        }
        for (int ip = 0; ip < sp.num; ip++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    sp.def_grad[ip][i][j] = eye_mat[i][j];
                    if (Fbar_flag == true) sp.def_grad_bar[ip][i][j] = eye_mat[i][j];
                }
            }
            sp.vol[ip] = sp.vol0[ip];
        }

        infile.close();
    }

    return;
}
