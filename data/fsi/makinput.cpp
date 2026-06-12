#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/fluid/FEM/stabilized_fem.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"
#include "../../module/solid/implicit/implicit_mpm_solid.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;
using namespace implicitmpm;
using namespace stabilizedfem;

void make_surface_point();

void create_data() {
    ifstream infile;
    infile.open("fsi/input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    sp.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> wfem.spec_rad >> sp.spec_rad;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xymin[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xymax[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> xyelem[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> sp.rho >> wfem.rhol >> wfem.rmul >> wfem.rhog >> wfem.rmug >> wfem.fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> nmat >> npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(nmat, iprop);
    VectorAssign(nmat, nprop);
    mat_prop.assign(nmat, vector<double>(npropmax, 0.0));
    for (int n = 0; n < nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> iprop[n] >> nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < nprop[n]; i++) { infile >> mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');

    infile.close();

    for (int i = 0; i < 3; i++) {
        dxy[i] = (xymax[i] - xymin[i]) / double(xyelem[i]);
        aelemmin[i] = 0;
        aelemmax[i] = xyelem[i] - 1;
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
    double xmax = xymax[0];
    double ymax = xymax[1];
    double zmax = xymax[2];

    nelem = xelem * yelem * zelem;
    node = xnode * ynode * znode;
    nodec = xnodec * ynodec * znodec;
    int nspe = npxye[0] * npxye[1] * npxye[2];
    double volp = dx * dy * dz / double(nspe);
    double mwp = wp.rho * volp;
    double msp = sp.rho * volp;

    // --- Gaussian distribution ---
    array<array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    // --- build mesh ---
    BuildMesh();

    BuildControlPoint(); // --- build control point ---

    wfem.ubc.nbc.resize(nodec);
    wfem.vbc.nbc.resize(nodec);
    wfem.wbc.nbc.resize(nodec);
    sp.ubc.nbc.resize(nodec);
    sp.vbc.nbc.resize(nodec);
    sp.wbc.nbc.resize(nodec);
    wfem.ubc.fbc.resize(nodec);
    wfem.vbc.fbc.resize(nodec);
    wfem.wbc.fbc.resize(nodec);
    sp.ubc.fbc.resize(nodec);
    sp.vbc.fbc.resize(nodec);
    sp.wbc.fbc.resize(nodec);
    wfem.pbc.nbc.resize(nodec);
    wfem.pbc.fbc.resize(nodec);
    uifbc.nbc.resize(nelem);
    vifbc.nbc.resize(nelem);
    wifbc.nbc.resize(nelem);

    sp.ubc.ibc = 0;
    sp.vbc.ibc = 0;
    sp.wbc.ibc = 0;
    wfem.ubc.ibc = 0;
    wfem.vbc.ibc = 0;
    wfem.wbc.ibc = 0;
    wfem.pbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == xnodec - 1 && (k != 0 || k != znodec - 1)) {
                    // --- Turek FSI ---
                    double vel = 1.5e0 * (4.0e0 / pow(zmax, 2)) * xyc[id][2] * (zmax - xyc[id][2]);
                    // --- Beam in fluid channel ---
                    // double vel = 1.5e-2 * (-pow(xyc[id][2]*100.0e0, 2) + 2.0e0 * xyc[id][2]*100.0e0);
                    // --- Cavity hyperelastic wall ---
                    // double vel = 1.0e-2;
                    // if (xyc[id][0] >= 0.0e0 && xyc[id][0] <= 0.3e-2) {
                    //     vel = pow(sin(M_PI * xyc[id][0] / 0.6e-2), 2) * 1.0e-2;
                    // } else if (xyc[id][0] >= 1.7e-2 && xyc[id][0] <= 2.0e-2) {
                    //     vel = pow(sin(M_PI * (xyc[0][id] - 2.0e-2) / 0.6e-2), 2) * 1.0e-2;
                    // } else {
                    //     vel = 1.0e-2;
                    // }
                    wfem.ubc.nbc[wfem.ubc.ibc] = id;
                    wfem.ubc.fbc[wfem.ubc.ibc] = vel;
                    wfem.ubc.ibc++;
                    wfem.wbc.nbc[wfem.wbc.ibc] = id;
                    wfem.wbc.fbc[wfem.wbc.ibc] = 0.0e0;
                    wfem.wbc.ibc++;
                }
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    sp.ubc.nbc[sp.ubc.ibc] = id;
                    sp.ubc.fbc[sp.ubc.ibc] = 0.0e0;
                    sp.ubc.ibc++;
                    sp.wbc.nbc[sp.wbc.ibc] = id;
                    sp.wbc.fbc[sp.wbc.ibc] = 0.0e0;
                    sp.wbc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    sp.vbc.nbc[sp.vbc.ibc] = id;
                    sp.vbc.fbc[sp.vbc.ibc] = 0.0e0;
                    sp.vbc.ibc++;
                    wfem.vbc.nbc[wfem.vbc.ibc] = id;
                    wfem.vbc.fbc[wfem.vbc.ibc] = 0.0e0;
                    wfem.vbc.ibc++;
                }
                if (k == 0 || k == znodec - 1) {
                    wfem.ubc.nbc[wfem.ubc.ibc] = id;
                    wfem.ubc.fbc[wfem.ubc.ibc] = 0.0e0;
                    wfem.ubc.ibc++;
                    wfem.wbc.nbc[wfem.wbc.ibc] = id;
                    wfem.wbc.fbc[wfem.wbc.ibc] = 0.0e0;
                    wfem.wbc.ibc++;
                }
                if (i == xnodec - 1) {
                    wfem.pbc.nbc[wfem.pbc.ibc] = id;
                    wfem.pbc.fbc[wfem.pbc.ibc] = 0.0e0;
                    wfem.pbc.ibc++;
                }
            }
        }
    }

    sp.coord.resize(nelem * nspe);
    sp.matid.resize(nelem * nspe);

    sp.num = 0;
    for (int k = 0; k < zelem; k++) {
        for (int j = 0; j < yelem; j++) {
            for (int i = 0; i < xelem; i++) {
                int mid = i + xelem * j + xelem * yelem * k;
                double ecx = dx * (double(i) + 0.5e0) + xmin;
                double ecy = dy * (double(j) + 0.5e0) + ymin;
                double ecz = dz * (double(k) + 0.5e0) + zmin;
                for (int kp = 0; kp < npxye[2]; kp++) {
                    double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < npxye[1]; jp++) {
                        double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < npxye[0]; ip++) {
                            double xp = ecx + dec2p[ip][0];
                            // --- Beam ---
                            // if (xp >= 1.0e-2 && xp <= 1.04e-2 && zp <= 8.0e-3) {
                            // --- Cavity hyperelastic wall ---
                            // if (zp <= 5.0e-3) {
                            // --- Ball in Cavity ---
                            // if (pow((xp - 6.0e-3), 2) + pow((zp - 5.0e-3), 2) <= pow(2.0e-3, 2)) {
                            // --- Turek FSI ---
                            if (zp >= 0.21e0 || zp <= 0.19e0) {
                                if (pow((xp - 0.2e0), 2) + pow((zp - 0.2e0), 2) < pow(0.05e0, 2)) {
                                    sp.coord[sp.num][0] = xp;
                                    sp.coord[sp.num][1] = yp;
                                    sp.coord[sp.num][2] = zp;
                                    sp.matid[sp.num] = 0;
                                    sp.num++;
                                }
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 && xp >= 0.15e0 && xp <= 0.25e0) {
                                sp.coord[0][sp.num] = xp;
                                sp.coord[1][sp.num] = yp;
                                sp.coord[2][sp.num] = zp;
                                sp.matid[sp.num] = 0;
                                sp.num++;
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 && xp >= 0.25e0 && xp <= 0.6e0) {
                                sp.coord[0][sp.num] = xp;
                                sp.coord[1][sp.num] = yp;
                                sp.coord[2][sp.num] = zp;
                                sp.matid[sp.num] = 1;
                                sp.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    // --- surface point flag ---
    VectorAssign(sp.num, sp.surf_point);
    make_surface_point();

    cout << "Making solid point file" << "\n";
    ofstream pointfile;
    pointfile.open("fsi/pointdata.txt");

    pointfile.flags(ios::right | ios::scientific); // output format
    pointfile << setw(10) << sp.num << "\n";

    // cout << "Solid point number: " << sp.num << "\n";

    OutputVector(pointfile, sp.num, sp.coord);
    for (int i = 0; i < sp.num; i++) {
        pointfile << setw(15) << i << setw(15) << sp.matid[i] << setw(15) << sp.surf_point[i] << setw(15) << msp
                  << setw(15) << volp << "\n";
    }

    pointfile.close();

    return;
}

void make_surface_point() {
    double leftside = numeric_limits<double>::infinity();   // min x
    double rightside = -numeric_limits<double>::infinity(); // max x
    double bottside = numeric_limits<double>::infinity();   // min z
    double upside = -numeric_limits<double>::infinity();    // max z

    for (int i = 0; i < sp.num; ++i) {
        const double x = sp.coord[i][0];
        const double z = sp.coord[i][2];
        if (x < leftside) leftside = x;
        if (x > rightside) rightside = x;
        if (z < bottside) bottside = z;
        if (z > upside) upside = z;
    }

    for (int i = 0; i < sp.num; ++i) {
        const double x = sp.coord[i][0];
        const double z = sp.coord[i][2];

        // if (abs(x - leftside)  <= mtol) {sp.surf_point[i] = 1; continue;}
        // if (abs(z - bottside) <= mtol) {sp.surf_point[i] = 1; continue;}
        // if (abs(x - rightside)<= mtol) {sp.surf_point[i] = 1; continue;}
        if (abs(z - upside) <= mtol) {
            sp.surf_point[i] = 1;
            continue;
        }
    }

    return;
}

void global_bc_data_output(ofstream &outfile) {
    // --- solid boundary condition ---
    sp.ubc.BCOutput(outfile, "usbc");
    sp.vbc.BCOutput(outfile, "vsbc");
    sp.wbc.BCOutput(outfile, "wsbc");

    // --- water boundary condition ---
    wfem.ubc.BCOutput(outfile, "uwbc");
    wfem.vbc.BCOutput(outfile, "vwbc");
    wfem.wbc.BCOutput(outfile, "wwbc");
    wfem.pbc.BCOutput(outfile, "hpbc");

    return;
}
