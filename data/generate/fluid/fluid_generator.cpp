#include "fluid_generator.h"

#include "../output_util.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../../module/data_io.h"
#include "../../../module/dataset.h"
#include "../../../module/material_point.h"
#include "../../../module/mesh.h"

namespace {

MaterialPoint g_wp;
double g_mwp = 0.0;
double g_volp = 0.0;

}  // namespace

std::string FluidGenerator::CaseName() const { return "fluid"; }

void FluidGenerator::Initialize() {}

void FluidGenerator::LoadInput() {
  std::ifstream infile = this->OpenInputFile("input.txt");

  infile.ignore(1000, '\n');
  std::string solswitch_str;
  infile >> solswitch_str >> rstflag >> nlstep;
  g_wp.solswitch = ParseMapScheme(solswitch_str);
  infile.ignore(1000, '\n');

  infile.ignore(1000, '\n');
  infile >> ista >> iend >> iout;
  infile.ignore(1000, '\n');

  infile.ignore(1000, '\n');
  infile >> dt >> mtol >> g_wp.gamma_nb >> g_wp.beta_nb;
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
  infile >> g_wp.rho >> g_wp.rmu;
  infile.ignore(1000, '\n');

  infile.ignore(1000, '\n');
  for (int i = 0; i < 3; i++) infile >> bb[i];
  infile.ignore(1000, '\n');
}

void FluidGenerator::BuildData() {
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
  int nspe = npxye[0] * npxye[1] * npxye[2];

  g_volp = dx * dy * dz / static_cast<double>(nspe);
  g_mwp = g_wp.rho * g_volp;

  std::array<std::array<double, 3>, 6> dec2p;
  GaussianDistribution(dec2p);

  BuildMesh();

  g_wp.ubc.nbc.resize(nodec);
  g_wp.vbc.nbc.resize(nodec);
  g_wp.wbc.nbc.resize(nodec);
  g_wp.pbc.nbc.resize(nodec);
  g_wp.ubc.fbc.resize(nodec);
  g_wp.vbc.fbc.resize(nodec);
  g_wp.wbc.fbc.resize(nodec);
  g_wp.pbc.fbc.resize(nodec);

  g_wp.ubc.ibc = 0;
  g_wp.vbc.ibc = 0;
  g_wp.wbc.ibc = 0;
  g_wp.pbc.ibc = 0;
  for (int k = 0; k < znodec; k++) {
    for (int j = 0; j < ynodec; j++) {
      for (int i = 0; i < xnodec; i++) {
        int id = i + xnodec * j + xnodec * ynodec * k;
        if (i == 0 || i == xnodec - 1) {
          g_wp.ubc.nbc[g_wp.ubc.ibc] = id;
          g_wp.ubc.fbc[g_wp.ubc.ibc] = 0.0e0;
          g_wp.ubc.ibc++;
        }
        if (j == 0 || j == ynodec - 1) {
          g_wp.vbc.nbc[g_wp.vbc.ibc] = id;
          g_wp.vbc.fbc[g_wp.vbc.ibc] = 0.0e0;
          g_wp.vbc.ibc++;
        }
        if (k == 0 || k == znodec - 1) {
          g_wp.wbc.nbc[g_wp.wbc.ibc] = id;
          g_wp.wbc.fbc[g_wp.wbc.ibc] = 0.0e0;
          g_wp.wbc.ibc++;
        }
      }
    }
  }

  g_wp.coord.resize(nelem * nspe);
  g_wp.matid.resize(nelem * nspe);
  g_wp.id.resize(nelem * nspe);

  g_wp.num = 0;
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
              if (zp <= 1.0e0) {
                g_wp.coord[g_wp.num][0] = xp;
                g_wp.coord[g_wp.num][1] = yp;
                g_wp.coord[g_wp.num][2] = zp;
                g_wp.matid[g_wp.num] = 0;
                g_wp.id[g_wp.num] = g_wp.num;
                g_wp.num++;
              }
            }
          }
        }
      }
    }
  }

  g_wp.coord.resize(g_wp.num);
  g_wp.matid.resize(g_wp.num);
  g_wp.id.resize(g_wp.num);
}

void FluidGenerator::WriteBcData(std::ofstream& outfile) {
  g_wp.ubc.BCOutput(outfile, "uwbc");
  g_wp.vbc.BCOutput(outfile, "vwbc");
  g_wp.wbc.BCOutput(outfile, "wwbc");
  g_wp.pbc.BCOutput(outfile, "hpbc");
}

void FluidGenerator::WriteTextOutputs() {
  std::cout << "Making grid file" << "\n";
  this->WriteGridDataFile("griddata.txt");

  std::cout << "Making water point file" << "\n";
  std::ofstream pointfile = this->OpenOutputFile("wpdata.txt");
  pointfile << std::setw(10) << g_wp.num << "\n";
  OutputVector(pointfile, g_wp.num, g_wp.coord);
  for (int i = 0; i < g_wp.num; i++) {
    pointfile << std::setw(15) << i << std::setw(15) << g_wp.matid[i]
              << std::setw(15) << g_mwp << std::setw(15) << g_volp << "\n";
  }
}

void FluidGenerator::WriteVisualizationOutputs() {
#ifdef HAVE_HDF5
  std::cout << "Making VTK HDF5 files" << "\n";
  WriteVtkHdf5Mesh("grid.vtkhdf");
  WriteVtkHdf5Points("wp.vtkhdf", g_wp);
#endif
}

void FluidGenerator::Finalize() {
}
