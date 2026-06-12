#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include <algorithm>
#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"



void datain_para () {

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
        infile >> dt >> mtol;
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

    string filename = gridfile + to_string(myrank) + ".txt";
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
                infile.ignore(1000,'\n');
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

    // --- Initial velocity of two rings collide ---
    //for (int i = 0; i < sp.num; i++) {
    //    if (sp.matid[i] == 0) {
    //        sp.vel[0][i] = 30.0e0;
    //    }
    //    else if (sp.matid[i] == 1) {
    //        sp.vel[0][i] = -30.0e0;
    //    }
    //}
   
    // vector<double> center_coord(3);
    // for (int i = 0; i < 3; i++) {
        // center_coord[i] = 0.5e0 * (xyminw[i] + xymaxw[i]);
    // }

    // for (int ip = 0; ip < sp.num; ip++) {
        // sp.vel[0][ip] =  100.0e0 * sin((M_PI*sp.coord[2][ip])/12.0e0) * (sp.coord[1][ip] - center_coord[1]);
        // sp.vel[1][ip] = -100.0e0 * sin((M_PI*sp.coord[2][ip])/12.0e0) * (sp.coord[0][ip] - center_coord[0]);
        // sp.vel[2][ip] =  0.0e0;
    // }

    //double max_xp2, max_xp;
    //max_xp2 = *max_element(sp.coord[0].begin(), sp.coord[0].end());
    //MPI_Allreduce(&max_xp2, &max_xp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //for (int i = 0; i < sp.num; i++) {
    //    if (abs(sp.coord[0][i] - max_xp) < mtol) {
    //        sp.trac_force[2][i] = 1.0e0;
    //    }
    //}

    return;
}
