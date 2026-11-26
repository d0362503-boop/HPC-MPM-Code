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
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;
    RunParams &params = this->params_;

    std::ifstream infile = OpenInputFile("input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> params.rstflag >> params.nlstep;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> params.ista >> params.iend >> params.iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    double spec_rad; // read past; not used by the generator
    infile >> params.dt >> params.mtol >> spec_rad;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> mesh.xymin[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> mesh.xymax[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> mesh.xyelem[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> cmesh.idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->points_.rho >> this->points_.rmu;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.bb[i];
    infile.ignore(1000, '\n');

    return;
}

void FluidGenerator::CreateBCs() {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;

    const int xnodec = cmesh.xynodec[0];
    const int ynodec = cmesh.xynodec[1];
    const int znodec = cmesh.xynodec[2];

    this->points_.ubc.nbc.resize(cmesh.nodec);
    this->points_.vbc.nbc.resize(cmesh.nodec);
    this->points_.wbc.nbc.resize(cmesh.nodec);
    this->points_.pbc.nbc.resize(cmesh.nodec);
    this->points_.ubc.fbc.resize(cmesh.nodec);
    this->points_.vbc.fbc.resize(cmesh.nodec);
    this->points_.wbc.fbc.resize(cmesh.nodec);
    this->points_.pbc.fbc.resize(cmesh.nodec);

    this->points_.uinfbc.nbc.resize(mesh.nelem);
    this->points_.vinfbc.nbc.resize(mesh.nelem);
    this->points_.winfbc.nbc.resize(mesh.nelem);

    this->points_.ubc.ibc = 0;
    this->points_.vbc.ibc = 0;
    this->points_.wbc.ibc = 0;
    this->points_.pbc.ibc = 0;
    this->points_.uinfbc.ibc = 0;
    this->points_.vinfbc.ibc = 0;
    this->points_.winfbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1) {
                    this->points_.ubc.nbc[this->points_.ubc.ibc] = id;
                    this->points_.ubc.fbc[this->points_.ubc.ibc] = 0.0e0;
                    this->points_.ubc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    this->points_.vbc.nbc[this->points_.vbc.ibc] = id;
                    this->points_.vbc.fbc[this->points_.vbc.ibc] = 0.0e0;
                    this->points_.vbc.ibc++;
                }
                if (k == 0) { // || k == znodec - 1) {
                    this->points_.wbc.nbc[this->points_.wbc.ibc] = id;
                    this->points_.wbc.fbc[this->points_.wbc.ibc] = 0.0e0;
                    this->points_.wbc.ibc++;
                }
                // if (i == 0 && k == 0) {
                //     this->points_.pbc.nbc[this->points_.pbc.ibc] = id;
                //     this->points_.pbc.fbc[this->points_.pbc.ibc] = 0.0e0;
                //     this->points_.pbc.ibc++;
                // }
            }
        }
    }

    return;
}

void FluidGenerator::CreateParticles() {
    Mesh &mesh = this->grid_.mesh;
    RunParams &params = this->params_;

    const int xelem = mesh.xyelem[0];
    const int yelem = mesh.xyelem[1];
    const int zelem = mesh.xyelem[2];
    const double dx = mesh.dxy[0];
    const double dy = mesh.dxy[1];
    const double dz = mesh.dxy[2];
    const double xmin = mesh.xymin[0];
    const double ymin = mesh.xymin[1];
    const double zmin = mesh.xymin[2];
    const int nspe = params.npxye[0] * params.npxye[1] * params.npxye[2];

    this->fluid_point_volume_ = dx * dy * dz / static_cast<double>(nspe);
    this->fluid_point_mass_ = this->points_.rho * this->fluid_point_volume_;

    std::array<std::array<double, 3>, 6> dec2p;
    SubcellOffsets(params.npxye, mesh.dxy, dec2p);

    this->points_.coord.resize(mesh.nelem * nspe);
    this->points_.matid.resize(mesh.nelem * nspe);
    this->points_.id.resize(mesh.nelem * nspe);

    this->points_.num = 0;
    for (int k = 0; k < zelem; k++) {
        for (int j = 0; j < yelem; j++) {
            for (int i = 0; i < xelem; i++) {
                double ecx = dx * (double(i) + 0.5e0) + xmin;
                double ecy = dy * (double(j) + 0.5e0) + ymin;
                double ecz = dz * (double(k) + 0.5e0) + zmin;
                for (int kp = 0; kp < params.npxye[2]; kp++) {
                    double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < params.npxye[1]; jp++) {
                        double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < params.npxye[0]; ip++) {
                            double xp = ecx + dec2p[ip][0];
                            if (xp <= 1.2e0 && zp <= 0.6e0) {
                                this->points_.coord[this->points_.num][0] = xp;
                                this->points_.coord[this->points_.num][1] = yp;
                                this->points_.coord[this->points_.num][2] = zp;
                                this->points_.matid[this->points_.num] = 0;
                                this->points_.id[this->points_.num] = this->points_.num;
                                this->points_.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    this->points_.coord.resize(this->points_.num);
    this->points_.matid.resize(this->points_.num);
    this->points_.id.resize(this->points_.num);

    return;
}

void FluidGenerator::WriteBCData(std::ofstream &outfile) {

    this->points_.ubc.BCOutput(outfile, "ubc");
    this->points_.vbc.BCOutput(outfile, "vbc");
    this->points_.wbc.BCOutput(outfile, "wbc");
    this->points_.pbc.BCOutput(outfile, "pbc");

    // Inflow boundaries (node IDs only)
    this->points_.uinfbc.BCOutput(outfile, "uinfbc", false);
    this->points_.vinfbc.BCOutput(outfile, "vinfbc", false);
    this->points_.winfbc.BCOutput(outfile, "winfbc", false);

    return;
}

void FluidGenerator::WritePointData(std::ofstream &pointfile) {

    pointfile << std::setw(10) << this->points_.num << "\n";
    OutputVector(pointfile, this->points_.num, this->points_.coord);
    for (int i = 0; i < this->points_.num; i++) {
        pointfile << std::setw(15) << i << std::setw(15) << this->points_.matid[i] //
                  << std::setw(15) << this->fluid_point_mass_ << std::setw(15) << this->fluid_point_volume_ << "\n";
    }

    return;
}

void FluidGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf", this->grid_.mesh);
    WriteVTKHDFPoints("wp.vtkhdf", this->points_);
#endif

    return;
}
