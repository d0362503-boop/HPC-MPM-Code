#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../partition_module.h"
#include "../../../module/material_point.h"
#include "../../../module/bc.h"
#include "../../../module/data_io.h"
#include "../../../module/mesh.h"
#include "../../../module/dataset.h"
#include "../../../module/solid/implicit/implicit_mpm_solid.h"

using namespace std;

void input_bc_data (ifstream& infile) {

    bc_input(infile, usbc);
    bc_input(infile, vsbc);
    bc_input(infile, wsbc);

    bc_input(infile, uwbc);
    bc_input(infile, vwbc);
    bc_input(infile, wwbc);
    bc_input(infile, hpbc);
    
    return;
}

void output_bc_data (ofstream& outfile) {

    bc_output(outfile, usbcr, "usbc");
    bc_output(outfile, vsbcr, "vsbc");
    bc_output(outfile, wsbcr, "wsbc");

    bc_output(outfile, uwbcr, "uwbc");
    bc_output(outfile, vwbcr, "vwbc");
    bc_output(outfile, wwbcr, "wwbc");
    bc_output(outfile, hpbcr, "hpbc");

    return;
}

void input_point_data (ifstream& infile) {

    infile >> sp.num; 
    infile.ignore(1000, '\n');
    // ---- Solid point data ----
    VectorAssign(sp.num, sp.id);
    VectorAssign(sp.num, sp.matid);
    VectorAssign(sp.num, sp.surf_point);
    VectorAssign(sp.num, sp.mass);
    VectorAssign(sp.num, sp.vol0);
    VectorAssign(sp.num, sp.coord);
    if (sp.num != 0) {
        InputVector(infile, sp.num, sp.coord);
        for (int i = 0; i < sp.num; i++) {
            infile >> sp.id[i] >> sp.matid[i] >> sp.surf_point[i] >> sp.mass[i] >> sp.vol0[i];
            infile.ignore(1000, '\n');
        }
    }

    return;
}

void output_point_data (ofstream& outfile) {

    outfile << setw(15) << spr.num << "\n";

    for (int i = 0; i < spr.num; i++) {
        int ip = spr.id[i];
        for (int j = 0; j < 3; j++) outfile << sp.coord[ip][j] << setw(15);
        outfile << "\n";
    }

    for (int i = 0; i < spr.num; i++) {
        int ip = spr.id[i];
        outfile << setw(15) << sp.id[ip] << setw(15) << sp.matid[ip] << setw(15) << sp.surf_point[ip]
                << setw(15) << sp.mass[ip] << setw(15) << sp.vol0[ip] << "\n";
    }

    return;
}