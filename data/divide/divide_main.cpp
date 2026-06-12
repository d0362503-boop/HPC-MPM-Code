#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include "../../module/bc.h"
#include "../../module/mesh.h"
#include "../../module/mpi_data.h"
#include "../../module/material_point.h"
#include "../../module/data_IO.h"
#include "../../module/dataset.h"
#include "partition_module.h"

using namespace std;

// -------------- void decleartion --------------
void input_file (ifstream& infile, string& gridfile, string& pointfile,     
                                   string& gridoutfile, string& pointoutfile);
void input_bc_data (ifstream& infile);
void input_point_data (ifstream& infile);
void partition_process (int ir); 
void output_bc_data (ofstream& outfile);
void output_point_data (ofstream& outfile);
// ----------------------------------------------

int main () {
	
    string gridfile, gridoutfile, pointfile, pointoutfile;

    ifstream infile;
    infile.open ("para_input_data.txt");
        input_file(infile, gridfile, pointfile, gridoutfile, pointoutfile);
    infile.close();

    infile.open(gridfile);
        input_mesh_data(infile);
        input_bc_data(infile);     // ---- CHANGE INPUT BC DATA HERE (IF YOU NEED) ----
    infile.close();

    infile.open(pointfile);
        input_point_data(infile);    // ---- CHANGE INPUT POINT DATA HERE (IF YOU NEED) ----
    infile.close();

    cout << "--- Read Finish ---" << "\n";

    int nr = nxyr[0] * nxyr[1] * nxyr[2];
    partition_initial_dataset();

    for (int ir = 0; ir < nr; ir++) {

        partition_process(ir);     // ---- CHANGE DIVIDE PROPERTY DATA HERE (IF YOU NEED) ----

        ofstream outfile;
        outfile.flags(ios::right | ios::scientific);

        string filename = gridoutfile + to_string(ir) + ".txt";
        outfile.open(filename);
            output_mesh_data(outfile, ir);
            output_bc_data(outfile);     // ---- CHANGE OUTPUT BC DATA HERE (IF YOU NEED) ----
        outfile.close();

        filename = pointoutfile + to_string(ir) + ".txt";
        outfile.open(filename);
            output_point_data(outfile);    // ---- CHANGE OUTPUT POINT DATA HERE (IF YOU NEED) ----
        outfile.close();
    } 

    cout << "--- Partition " << nr << " Cores ----" << "\n";
    cout << "--- Finished ---" << "\n";

    return 0;
}  

void input_file (ifstream& infile, string& gridfile, string& pointfile, 
                                   string& gridoutfile, string& pointoutfile) {

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> nxyr[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    getline(infile, gridfile);

    infile.ignore(1000, '\n');
    getline(infile, pointfile);

    infile.ignore(1000, '\n');
    getline(infile, gridoutfile);

    infile.ignore(1000, '\n');
    getline(infile, pointoutfile);

    return;
}