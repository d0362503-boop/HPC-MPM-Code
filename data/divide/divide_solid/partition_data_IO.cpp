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

    bc_input(infile, usbc);
    bc_input(infile, vsbc);
    bc_input(infile, wsbc);

    return;
}

void output_bc_data (ofstream& outfile) {

    bc_output(outfile, usbcr, "usbc");
    bc_output(outfile, vsbcr, "vsbc");
    bc_output(outfile, wsbcr, "wsbc");

    return;
}

void input_point_data (ifstream& infile) {

    infile >> sp.num; 
    infile.ignore(1000, '\n');
    oned_vector_assign(sp.num, sp.id);
    oned_vector_assign(sp.num, sp.matid);
    oned_vector_assign(sp.num, sp.mass);
    oned_vector_assign(sp.num, sp.vol0);
    twod_vector_assign(3, sp.num, sp.coord);
    if (sp.num != 0) {
        input_twod_vector(infile, 3, sp.num, sp.coord);
        for (int i = 0; i < sp.num; i++) {
            infile >> sp.id[i] >> sp.matid[i] >> sp.mass[i] >> sp.vol0[i];
            infile.ignore(1000, '\n');
        }
    }
    return;
}

void output_point_data (ofstream& outfile) {

    outfile << setw(15) << spr.num << "\n";
    for (int i = 0; i < spr.num; i++) {
        int ip = spr.id[i];
        for (int j = 0; j < 3; j++) outfile << sp.coord[j][ip] << setw(15);
        outfile << "\n";
    }

    for (int i = 0; i < spr.num; i++) {
        int ip = spr.id[i];
        outfile << setw(15) << sp.id[ip] << setw(15) << sp.matid[ip]
                << setw(15) << sp.mass[ip] << setw(15) << sp.vol0[ip] << "\n";
    }

    return;
}