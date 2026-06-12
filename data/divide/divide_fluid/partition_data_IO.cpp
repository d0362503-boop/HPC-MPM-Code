#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../partition_module.h"
#include "../../../module/material_point.h"
#include "../../../module/bc.h"
#include "../../../module/data_IO.h"
#include "../../../module/mesh.h"
#include "../../../module/dataset.h"

using namespace std;

void input_bc_data (ifstream& infile) {

    bc_input(infile, uwbc);
    bc_input(infile, vwbc);
    bc_input(infile, wwbc);
    bc_input(infile, hpbc);

    bc_input(infile, uifbc, false);
    bc_input(infile, vifbc, false);
    bc_input(infile, wifbc, false);

    return;
}

void output_bc_data (ofstream& outfile) {

    bc_output(outfile, uwbcr, "uwbc");
    bc_output(outfile, vwbcr, "vwbc");
    bc_output(outfile, wwbcr, "wwbc");
    bc_output(outfile, hpbcr, "hpbc");

    // --- Inflow water boundary condition ---
    bc_output(outfile, uifbcr, "uifbc", false);
    bc_output(outfile, vifbcr, "vifbc", false);
    bc_output(outfile, wifbcr, "wifbc", false);
    
    return;
}

void input_point_data (ifstream& infile) {

    infile >> wp.num; 
    infile.ignore(1000, '\n');
    oned_vector_assign(wp.num, wp.id);
    oned_vector_assign(wp.num, wp.matid);
    oned_vector_assign(wp.num, wp.mass);
    oned_vector_assign(wp.num, wp.vol0);
    twod_vector_assign(3, wp.num, wp.coord);
    if (wp.num != 0) {
        input_twod_vector(infile, 3, wp.num, wp.coord);
        for (int i = 0; i < wp.num; i++) {
            infile >> wp.id[i] >> wp.matid[i] >> wp.mass[i] >> wp.vol0[i];
            infile.ignore(1000, '\n');
        }
    }
    return;
}

void output_point_data (ofstream& outfile) {

    outfile << setw(15) << wpr.num << "\n";
    for (int i = 0; i < wpr.num; i++) {
        int ip = wpr.id[i];
        for (int j = 0; j < 3; j++) outfile << wp.coord[j][ip] << setw(15);
        outfile << "\n";
    }

    for (int i = 0; i < wpr.num; i++) {
        int ip = wpr.id[i];
        outfile << setw(15) << wp.id[ip] << setw(15) << wp.matid[ip]
                << setw(15) << wp.mass[ip] << setw(15) << wp.vol0[ip] << "\n";
    }

    return;
}