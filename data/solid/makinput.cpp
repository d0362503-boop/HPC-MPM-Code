#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include "../../module/mesh.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/bc.h"

using namespace std;

void create_data () {

    ifstream infile;
    infile.open("solid/input.txt");
        infile.ignore(1000, '\n');
        std::string solswitch_str;
        infile >> solswitch_str >> NR_flag >> rstflag >> nlstep;
        sp.solswitch = ParseMapScheme(solswitch_str);
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        infile >> ista >> iend >> iout;
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        infile >> dt >> mtol >> sp.gamma_nb >> sp.beta_nb;
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
        infile >> sp.rho >> nmat >> npropmax; 
        infile.ignore(1000, '\n');

        oned_vector_assign(nmat, iprop);
        oned_vector_assign(nmat, nprop);
        twod_vector_assign(nmat, npropmax, mat_prop);
        for (int n = 0; n < nmat; n++) {
            infile.ignore(1000, '\n');
            infile >> iprop[n] >> nprop[n];
            infile.ignore(1000, '\n');
            
            infile.ignore(1000, '\n');
            for (int i = 0; i < nprop[n]; i++) {
                infile >> mat_prop[n][i];
            }
            infile.ignore(1000, '\n');
        }

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
    double msp = sp.rho * volp;

    // --- Gaussian distribution ---
    vector<vector<double>> dec2p;
    gaussian_distribution(dec2p);

    // --- build mesh ---
    build_mesh();

    //uwbc.nbc.resize(nodec); vwbc.nbc.resize(nodec); wwbc.nbc.resize(nodec);
    usbc.nbc.resize(nodec); vsbc.nbc.resize(nodec); wsbc.nbc.resize(nodec);
    //uwbc.fbc.resize(nodec); vwbc.fbc.resize(nodec); wwbc.fbc.resize(nodec);
    usbc.fbc.resize(nodec); vsbc.fbc.resize(nodec); wsbc.fbc.resize(nodec);
    //hpbc.nbc.resize(nodec); hpbc.fbc.resize(nodec);

    usbc.ibc = 0; 
    vsbc.ibc = 0; 
    wsbc.ibc = 0; 
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec-1 || k == 0 || k == znodec-1) {
                    usbc.nbc[usbc.ibc] = id;
                    usbc.fbc[usbc.ibc] = 0.0e0;
                    usbc.ibc++;
                }
                if (k == 0 || j == 0 || j == ynodec-1) {
                    vsbc.nbc[vsbc.ibc] = id;
                    vsbc.fbc[vsbc.ibc] = 0.0e0;
                    vsbc.ibc++;
                }
                if (i == 0 || i == xnodec-1 || k == 0 || k == znodec-1) {
                    wsbc.nbc[wsbc.ibc] = id;
                    wsbc.fbc[wsbc.ibc] = 0.0e0;
                    wsbc.ibc++;
                }
            }
        }
    }

    sp.coord.resize(3, vector<double>(nelem*nspe));
    sp.matid.resize(nelem*nspe);

    sp.num = 0;
    for (int k = 0; k < zelem; k++) {
        for (int j = 0; j < yelem; j++) {
            for (int i = 0; i < xelem; i++) {
                double ecx = dx * (double(i) + 0.5e0) + xmin;
                double ecy = dy * (double(j) + 0.5e0) + ymin;
                double ecz = dz * (double(k) + 0.5e0) + zmin;
                for (int kp = 0; kp < npxye[2]; kp++) {
                    double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < npxye[1]; jp++) {
                        double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < npxye[0]; ip++) {
                            double xp = ecx + dec2p[ip][0];
                            //if (pow((xp-0.05e0),2)+pow((zp-0.075e0),2) >= pow(0.03e0,2) &&
                            //    pow((xp-0.05e0),2)+pow((zp-0.075e0),2) <= pow(0.04e0,2)) {  // --- ring collide ---
                            if (zp >= 3.5e0 && zp <= 4.5e0 && xp <= 4.0e0) { 
                                sp.coord[sp.num][0] = xp;
                                sp.coord[sp.num][1] = yp;
                                sp.coord[sp.num][2] = zp;
                                sp.matid[sp.num] = 0;
                                sp.num++;
                            }
                            //else if (pow((xp-0.15e0),2)+pow((zp-0.075e0),2) >= pow(0.03e0,2) && 
                            //         pow((xp-0.15e0),2)+pow((zp-0.075e0),2) <= pow(0.04e0,2)) {
                            //    sp.coord[0][sp.num] = xp;
                            //    sp.coord[1][sp.num] = yp;
                            //    sp.coord[2][sp.num] = zp;
                            //    sp.matid[sp.num] = 1;
                            //    sp.num++;
                            //}
                        }
                    }
                }
            }
        }
    }
  
    cout << "Making solid point file" << "\n";
    ofstream pointfile;
    pointfile.open("solid/spdata.txt");
        pointfile.flags(ios::right | ios::scientific); // output format
        pointfile << setw(10) << sp.num << "\n";
        for (int i = 0; i < sp.num; i++) {
            pointfile << setw(15) << sp.coord[i][0]
                      << setw(15) << sp.coord[i][1]
                      << setw(15) << sp.coord[i][2] << "\n";
        }
        for (int i = 0; i < sp.num; i++) {
            pointfile << setw(10) << i 
                      << setw(10) << sp.matid[i] 
                      << setw(15) << msp 
                      << setw(15) << volp << "\n";
        } 
    pointfile.close();

    return;
}

void global_bc_data_output (ofstream& outfile) {

    bc_output(outfile, usbc, "usbc");
    bc_output(outfile, vsbc, "vsbc");
    bc_output(outfile, wsbc, "wsbc");

    return;
}
