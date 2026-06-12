#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include "../../module/Newmark_beta_data.h"
#include "../../module/mesh.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/bc.h"

using namespace std;

void create_data () {

    ifstream infile;
    infile.open("fluid/input.txt");
        infile.ignore(1000, '\n');
        std::string solswitch_str;
        infile >> solswitch_str >> rstflag >> nlstep;
        wp.solswitch = ParseMapScheme(solswitch_str);
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        infile >> ista >> iend >> iout;
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        infile >> dt >> mtol >> wp.gamma_nb >> wp.beta_nb;
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> xymin[i];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> xymax[i];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> xyelem[i];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> idimc[i];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> npxye[i];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        infile >> wp.rho >> wp.rmu;
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for(int i = 0; i < 3; i++) infile >> bb[i];
        infile.ignore(1000, '\n');
    infile.close();
 
    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymax[i] - xymin[i]) / double(xyelem[i]);
        xynode[i] = xyelem[i] + 1;
        xynodec[i] = xyelem[i] + idimc[i];
    }
    int xnode = xynode[0]; 
    int ynode = xynode[1]; 
    int znode = xynode[2];
    int xnodec = xynodec[0]; 
    int ynodec = xynodec[1]; 
    int znodec = xynodec[2];
    int xelem = xyelem[0]; 
    int yelem = xyelem[1]; 
    int zelem = xyelem[2];
    double dx = dxy[0]; 
    double dy = dxy[1]; 
    double dz = dxy[2];
    double xmin = xymin[0]; 
    double ymin = xymin[1]; 
    double zmin = xymin[2];

    nelem = xelem * yelem * zelem;
    node = xnode * ynode * znode;
    nodec = xnodec * ynodec * znodec;
    int nspe = npxye[0] * npxye[1] * npxye[2];
    double volp = dx * dy * dz / double(nspe);
    double mwp = wp.rho * volp;

    // --- Gaussian distribution ---
    vector<vector<double>> dec2p;
    gaussian_distribution(dec2p);

    // --- build mesh ---
    build_mesh();

    uwbc.nbc.resize(nodec); vwbc.nbc.resize(nodec); wwbc.nbc.resize(nodec);
    //usbc.nbc.resize(nodec); vsbc.nbc.resize(nodec); wsbc.nbc.resize(nodec);
    uwbc.fbc.resize(nodec); vwbc.fbc.resize(nodec); wwbc.fbc.resize(nodec);
    //usbc.fbc.resize(nodec); vsbc.fbc.resize(nodec); wsbc.fbc.resize(nodec);
    hpbc.nbc.resize(nodec); hpbc.fbc.resize(nodec);
    uifbc.nbc.resize(nelem); vifbc.nbc.resize(nelem); wifbc.nbc.resize(nelem);

    uwbc.ibc = 0; 
    vwbc.ibc = 0; 
    wwbc.ibc = 0; 
    hpbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec-1) {// || k == 0) {
                    uwbc.nbc[uwbc.ibc] = id;
                    uwbc.fbc[uwbc.ibc] = 0.0e0;
                    uwbc.ibc++;
                }
                if (j == 0 || j == ynodec-1) {
                    vwbc.nbc[vwbc.ibc] = id;
                    vwbc.fbc[vwbc.ibc] = 0.0e0;
                    vwbc.ibc++;
                }
                if (k == 0 || k == znodec-1) {
                    wwbc.nbc[wwbc.ibc] = id;
                    wwbc.fbc[wwbc.ibc] = 0.0e0;
                    wwbc.ibc++;
                }
            }
        }
    }

    wp.coord.resize(3, vector<double>(nelem*nspe));
    wp.matid.resize(nelem*nspe);

    uifbc.ibc = 0;
    vifbc.ibc = 0;
    wifbc.ibc = 0;
    wp.num = 0;
    for (int k = 0; k < zelem; k++) {
        for (int j = 0; j < yelem; j++) {
            for (int i = 0; i < xelem; i++) {
                double ecx = dx * (double(i) + 0.5e0) + xmin;
                double ecy = dy * (double(j) + 0.5e0) + ymin;
                double ecz = dz * (double(k) + 0.5e0) + zmin;
                // // --- inflow boundary ---
                // if (i == 0) {
                //     uifbc.nbc[uifbc.ibc++] = mid;  
                // }
                for (int kp = 0; kp < npxye[2]; kp++) {
                    double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < npxye[1]; jp++) {
                        double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < npxye[0]; ip++) {
                            double xp = ecx + dec2p[ip][0];
                            if (zp <= 1.0e0) {
                              wp.coord[wp.num][0] = xp;
                              wp.coord[wp.num][1] = yp;
                              wp.coord[wp.num][2] = zp;
                              wp.matid[wp.num] = 0;
                              wp.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "Making water point file" << "\n";
    ofstream pointfile;
    pointfile.open("fluid/wpdata.txt");

        pointfile.flags(ios::right | ios::scientific); // output format
        pointfile << setw(10) << wp.num << "\n";
        for (int i = 0; i < wp.num; i++) {
            pointfile << setw(15) << wp.coord[i][0]
                      << setw(15) << wp.coord[i][1]
                      << setw(15) << wp.coord[i][2] << "\n";
        }
        for (int i = 0; i < wp.num; i++) {
            pointfile << setw(15) << i 
                      << setw(15) << wp.matid[i] 
                      << setw(15) << mwp 
                      << setw(15) << volp << "\n";
        } 

    pointfile.close();

    return;
}

void global_bc_data_output (ofstream& outfile) {

    bc_output(outfile, uwbc, "uwbc");
    bc_output(outfile, vwbc, "vwbc");
    bc_output(outfile, wwbc, "wwbc");
    bc_output(outfile, hpbc, "hpbc");
    
    // --- Inflow water boundary condition ---
    bc_output(outfile, uifbc, "uifbc", false);
    bc_output(outfile, vifbc, "vifbc", false);
    bc_output(outfile, wifbc, "wifbc", false);

    return;
}
