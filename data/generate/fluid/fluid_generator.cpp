#include "data/generate/fluid/fluid_generator.h"

#include "data/generate/output_util.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"

std::string FluidGenerator::CaseName() const { return "fluid"; }

void FluidGenerator::LoadInput() {

    std::ifstream infile = OpenInputFile("input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> this->spec_rad;
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
    infile >> this->rho >> this->rmu;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
    infile.ignore(1000, '\n');

    return;
}

void FluidGenerator::CreateBCs() {

    const int xnodec = xynodec[0];
    const int ynodec = xynodec[1];
    const int znodec = xynodec[2];

    this->ubc.nbc.resize(nodec);
    this->vbc.nbc.resize(nodec);
    this->wbc.nbc.resize(nodec);
    this->pbc.nbc.resize(nodec);
    this->ubc.fbc.resize(nodec);
    this->vbc.fbc.resize(nodec);
    this->wbc.fbc.resize(nodec);
    this->pbc.fbc.resize(nodec);

    this->uinfbc.nbc.resize(nelem);
    this->vinfbc.nbc.resize(nelem);
    this->winfbc.nbc.resize(nelem);

    this->ubc.ibc = 0;
    this->vbc.ibc = 0;
    this->wbc.ibc = 0;
    this->pbc.ibc = 0;
    this->uinfbc.ibc = 0;
    this->vinfbc.ibc = 0;
    this->winfbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1) {
                    this->ubc.nbc[this->ubc.ibc] = id;
                    this->ubc.fbc[this->ubc.ibc] = 0.0e0;
                    this->ubc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    this->vbc.nbc[this->vbc.ibc] = id;
                    this->vbc.fbc[this->vbc.ibc] = 0.0e0;
                    this->vbc.ibc++;
                }
                if (k == 0) { // || k == znodec - 1) {
                    this->wbc.nbc[this->wbc.ibc] = id;
                    this->wbc.fbc[this->wbc.ibc] = 0.0e0;
                    this->wbc.ibc++;
                }
                // if (i == 0 && k == 0) {
                //     this->pbc.nbc[this->pbc.ibc] = id;
                //     this->pbc.fbc[this->pbc.ibc] = 0.0e0;
                //     this->pbc.ibc++;
                // }
            }
        }
    }

    return;
}

void FluidGenerator::CreateParticles() {

    const int xelem = xyelem[0];
    const int yelem = xyelem[1];
    const int zelem = xyelem[2];
    const double dx = dxy[0];
    const double dy = dxy[1];
    const double dz = dxy[2];
    const double xmin = xymin[0];
    const double ymin = xymin[1];
    const double zmin = xymin[2];
    const int nspe = npxye[0] * npxye[1] * npxye[2];

    this->fluid_point_volume_ = dx * dy * dz / static_cast<double>(nspe);
    this->fluid_point_mass_ = this->rho * this->fluid_point_volume_;

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    this->coord.resize(nelem * nspe);
    this->matid.resize(nelem * nspe);
    this->id.resize(nelem * nspe);

    this->num = 0;
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
                            if (xp <= 1.2e0 && zp <= 0.6e0) {
                                this->coord[this->num][0] = xp;
                                this->coord[this->num][1] = yp;
                                this->coord[this->num][2] = zp;
                                this->matid[this->num] = 0;
                                this->id[this->num] = this->num;
                                this->num++;
                            }
                        }
                    }
                }
            }
        }
    }

    this->coord.resize(this->num);
    this->matid.resize(this->num);
    this->id.resize(this->num);

    return;
}

void FluidGenerator::WriteBCData(std::ofstream &outfile) {

    this->ubc.BCOutput(outfile, "ubc");
    this->vbc.BCOutput(outfile, "vbc");
    this->wbc.BCOutput(outfile, "wbc");
    this->pbc.BCOutput(outfile, "pbc");

    // Inflow boundaries (node IDs only)
    this->uinfbc.BCOutput(outfile, "uinfbc", false);
    this->vinfbc.BCOutput(outfile, "vinfbc", false);
    this->winfbc.BCOutput(outfile, "winfbc", false);

    return;
}

void FluidGenerator::WritePointData(std::ofstream &pointfile) {

    pointfile << std::setw(10) << this->num << "\n";
    OutputVector(pointfile, this->num, this->coord);
    for (int i = 0; i < this->num; i++) {
        pointfile << std::setw(15) << i << std::setw(15) << this->matid[i] //
                  << std::setw(15) << this->fluid_point_mass_ << std::setw(15) << this->fluid_point_volume_ << "\n";
    }

    return;
}

void FluidGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf");
    WriteVTKHDFPoints("wp.vtkhdf", *this);
#endif

    return;
}
