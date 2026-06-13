#include "fsi_generator.h"

#include "../output_util.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "../../../module/data_io.h"
#include "../../../module/dataset.h"
#include "../../../module/fluid/FEM/stabilized_fem.h"
#include "../../../module/material_point.h"
#include "../../../module/mesh.h"

using namespace stabilizedfem;

namespace {

MaterialPoint g_sp;
StabilizedFEM g_wfem;
double g_msp = 0.0;
double g_volp = 0.0;

void CheckSurfacePoint() {
    double leftside = std::numeric_limits<double>::infinity();
    double rightside = -std::numeric_limits<double>::infinity();
    double bottside = std::numeric_limits<double>::infinity();
    double upside = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < g_sp.num; ++i) {
        const double x = g_sp.coord[i][0];
        const double z = g_sp.coord[i][2];
        if (x < leftside) leftside = x;
        if (x > rightside) rightside = x;
        if (z < bottside) bottside = z;
        if (z > upside) upside = z;
    }

    for (int i = 0; i < g_sp.num; ++i) {
        const double z = g_sp.coord[i][2];
        if (std::abs(z - upside) <= mtol) { g_sp.surf_point[i] = 1; }
    }
}

} // namespace

std::string FsiGenerator::CaseName() const { return "fsi"; }

void FsiGenerator::Initialize() {}

void FsiGenerator::LoadInput() {
  std::ifstream infile = this->OpenInputFile("input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    g_sp.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> g_wfem.spec_rad >> g_sp.spec_rad;
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
    infile >> g_sp.rho >> g_wfem.rhol >> g_wfem.rmul >> g_wfem.rhog >> g_wfem.rmug >> g_wfem.fs_height;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> nmat >> npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(nmat, iprop);
    VectorAssign(nmat, nprop);
    mat_prop.assign(nmat, std::vector<double>(npropmax, 0.0));
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
}

void FsiGenerator::BuildData() {

    this->InitializeGridGeometry();

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
    double zmax = xymax[2];
    int nspe = npxye[0] * npxye[1] * npxye[2];

    g_volp = dx * dy * dz / static_cast<double>(nspe);
    g_msp = g_sp.rho * g_volp;

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    BuildMesh();
    BuildControlPoint();

    g_wfem.ubc.nbc.resize(nodec);
    g_wfem.vbc.nbc.resize(nodec);
    g_wfem.wbc.nbc.resize(nodec);
    g_sp.ubc.nbc.resize(nodec);
    g_sp.vbc.nbc.resize(nodec);
    g_sp.wbc.nbc.resize(nodec);
    g_wfem.ubc.fbc.resize(nodec);
    g_wfem.vbc.fbc.resize(nodec);
    g_wfem.wbc.fbc.resize(nodec);
    g_sp.ubc.fbc.resize(nodec);
    g_sp.vbc.fbc.resize(nodec);
    g_sp.wbc.fbc.resize(nodec);
    g_wfem.pbc.nbc.resize(nodec);
    g_wfem.pbc.fbc.resize(nodec);

    g_sp.ubc.ibc = 0;
    g_sp.vbc.ibc = 0;
    g_sp.wbc.ibc = 0;
    g_wfem.ubc.ibc = 0;
    g_wfem.vbc.ibc = 0;
    g_wfem.wbc.ibc = 0;
    g_wfem.pbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == xnodec - 1 && (k != 0 || k != znodec - 1)) {
                    double vel = 1.5e0 * (4.0e0 / std::pow(zmax, 2)) * xyc[id][2] * (zmax - xyc[id][2]);
                    g_wfem.ubc.nbc[g_wfem.ubc.ibc] = id;
                    g_wfem.ubc.fbc[g_wfem.ubc.ibc] = vel;
                    g_wfem.ubc.ibc++;
                    g_wfem.wbc.nbc[g_wfem.wbc.ibc] = id;
                    g_wfem.wbc.fbc[g_wfem.wbc.ibc] = 0.0e0;
                    g_wfem.wbc.ibc++;
                }
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    g_sp.ubc.nbc[g_sp.ubc.ibc] = id;
                    g_sp.ubc.fbc[g_sp.ubc.ibc] = 0.0e0;
                    g_sp.ubc.ibc++;
                    g_sp.wbc.nbc[g_sp.wbc.ibc] = id;
                    g_sp.wbc.fbc[g_sp.wbc.ibc] = 0.0e0;
                    g_sp.wbc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    g_sp.vbc.nbc[g_sp.vbc.ibc] = id;
                    g_sp.vbc.fbc[g_sp.vbc.ibc] = 0.0e0;
                    g_sp.vbc.ibc++;
                    g_wfem.vbc.nbc[g_wfem.vbc.ibc] = id;
                    g_wfem.vbc.fbc[g_wfem.vbc.ibc] = 0.0e0;
                    g_wfem.vbc.ibc++;
                }
                if (k == 0 || k == znodec - 1) {
                    g_wfem.ubc.nbc[g_wfem.ubc.ibc] = id;
                    g_wfem.ubc.fbc[g_wfem.ubc.ibc] = 0.0e0;
                    g_wfem.ubc.ibc++;
                    g_wfem.wbc.nbc[g_wfem.wbc.ibc] = id;
                    g_wfem.wbc.fbc[g_wfem.wbc.ibc] = 0.0e0;
                    g_wfem.wbc.ibc++;
                }
                if (i == xnodec - 1) {
                    g_wfem.pbc.nbc[g_wfem.pbc.ibc] = id;
                    g_wfem.pbc.fbc[g_wfem.pbc.ibc] = 0.0e0;
                    g_wfem.pbc.ibc++;
                }
            }
        }
    }

    g_sp.coord.resize(nelem * nspe);
    g_sp.matid.resize(nelem * nspe);
    g_sp.id.resize(nelem * nspe);

    g_sp.num = 0;
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
                            if (zp >= 0.21e0 || zp <= 0.19e0) {
                                if (std::pow((xp - 0.2e0), 2) + std::pow((zp - 0.2e0), 2) < std::pow(0.05e0, 2)) {
                                    g_sp.coord[g_sp.num][0] = xp;
                                    g_sp.coord[g_sp.num][1] = yp;
                                    g_sp.coord[g_sp.num][2] = zp;
                                    g_sp.matid[g_sp.num] = 0;
                                    g_sp.id[g_sp.num] = g_sp.num;
                                    g_sp.num++;
                                }
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 && xp >= 0.15e0 && xp <= 0.25e0) {
                                g_sp.coord[g_sp.num][0] = xp;
                                g_sp.coord[g_sp.num][1] = yp;
                                g_sp.coord[g_sp.num][2] = zp;
                                g_sp.matid[g_sp.num] = 0;
                                g_sp.id[g_sp.num] = g_sp.num;
                                g_sp.num++;
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 && xp >= 0.25e0 && xp <= 0.6e0) {
                                g_sp.coord[g_sp.num][0] = xp;
                                g_sp.coord[g_sp.num][1] = yp;
                                g_sp.coord[g_sp.num][2] = zp;
                                g_sp.matid[g_sp.num] = 1;
                                g_sp.id[g_sp.num] = g_sp.num;
                                g_sp.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    g_sp.coord.resize(g_sp.num);
    g_sp.matid.resize(g_sp.num);
    g_sp.id.resize(g_sp.num);

    VectorAssign(g_sp.num, g_sp.surf_point);
    CheckSurfacePoint();
}

void FsiGenerator::WriteBcData(std::ofstream &outfile) {
    g_sp.ubc.BCOutput(outfile, "usbc");
    g_sp.vbc.BCOutput(outfile, "vsbc");
    g_sp.wbc.BCOutput(outfile, "wsbc");

    g_wfem.ubc.BCOutput(outfile, "uwbc");
    g_wfem.vbc.BCOutput(outfile, "vwbc");
    g_wfem.wbc.BCOutput(outfile, "wwbc");
    g_wfem.pbc.BCOutput(outfile, "hpbc");
}

void FsiGenerator::WriteTextOutputs() {
  std::cout << "Making grid file" << "\n";
  this->WriteGridDataFile("griddata.txt");

  std::cout << "Making solid point file" << "\n";
  std::ofstream pointfile = this->OpenOutputFile("pointdata.txt");
    pointfile << std::setw(10) << g_sp.num << "\n";
    OutputVector(pointfile, g_sp.num, g_sp.coord);
    for (int i = 0; i < g_sp.num; i++) {
        pointfile << std::setw(15) << i << std::setw(15) << g_sp.matid[i] << std::setw(15) << g_sp.surf_point[i]
                  << std::setw(15) << g_msp << std::setw(15) << g_volp << "\n";
    }
}

void FsiGenerator::WriteVisualizationOutputs() {
#ifdef HAVE_HDF5
  std::cout << "Making VTK HDF5 files" << "\n";
  WriteVtkHdf5Mesh("grid.vtkhdf");
  WriteVtkHdf5Points("sp.vtkhdf", g_sp);
#endif
}

void FsiGenerator::Finalize() {}
