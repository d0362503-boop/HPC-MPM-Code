#include "data/generate/solid/solid_generator.h"

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

std::string SolidGenerator::CaseName() const { return "solid"; }

void SolidGenerator::LoadInput() {
    Mesh &mesh = this->grid_.mesh;
    ControlMesh &cmesh = this->grid_.cmesh;
    RunParams &params = this->params_;

    std::ifstream infile = OpenInputFile("input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    bool nr_flag; // read past; not used by the generator
    infile >> solswitch_str >> nr_flag >> params.rstflag >> params.nlstep;
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
    infile >> this->points_.rho >> this->materials_.nmat >> this->materials_.npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(this->materials_.nmat, this->materials_.iprop);
    VectorAssign(this->materials_.nmat, this->materials_.nprop);
    this->materials_.mat_prop.assign(this->materials_.nmat, std::vector<double>(this->materials_.npropmax, 0.0));
    for (int n = 0; n < this->materials_.nmat; n++) {
        infile.ignore(1000, '\n');
        infile >> this->materials_.iprop[n] >> this->materials_.nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < this->materials_.nprop[n]; i++) { infile >> this->materials_.mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> params.bb[i];
    infile.ignore(1000, '\n');

    return;
}

void SolidGenerator::CreateBCs() {
    ControlMesh &cmesh = this->grid_.cmesh;

    const int xnodec = cmesh.xynodec[0];
    const int ynodec = cmesh.xynodec[1];
    const int znodec = cmesh.xynodec[2];

    this->points_.ubc.nbc.resize(cmesh.nodec);
    this->points_.vbc.nbc.resize(cmesh.nodec);
    this->points_.wbc.nbc.resize(cmesh.nodec);
    this->points_.ubc.fbc.resize(cmesh.nodec);
    this->points_.vbc.fbc.resize(cmesh.nodec);
    this->points_.wbc.fbc.resize(cmesh.nodec);

    this->points_.ubc.ibc = 0;
    this->points_.vbc.ibc = 0;
    this->points_.wbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    this->points_.ubc.nbc[this->points_.ubc.ibc] = id;
                    this->points_.ubc.fbc[this->points_.ubc.ibc] = 0.0e0;
                    this->points_.ubc.ibc++;
                }
                if (k == 0 || j == 0 || j == ynodec - 1) {
                    this->points_.vbc.nbc[this->points_.vbc.ibc] = id;
                    this->points_.vbc.fbc[this->points_.vbc.ibc] = 0.0e0;
                    this->points_.vbc.ibc++;
                }
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    this->points_.wbc.nbc[this->points_.wbc.ibc] = id;
                    this->points_.wbc.fbc[this->points_.wbc.ibc] = 0.0e0;
                    this->points_.wbc.ibc++;
                }
            }
        }
    }

    return;
}

void SolidGenerator::CreateParticles() {
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

    this->solid_point_volume_ = dx * dy * dz / static_cast<double>(nspe);
    this->solid_point_mass_ = this->points_.rho * this->solid_point_volume_;

    std::array<std::array<double, 3>, 6> dec2p;
    SubcellOffsets(params.npxye, mesh.dxy, dec2p);

    this->points_.coord.resize(mesh.nelem * nspe);
    this->points_.matid.resize(mesh.nelem * nspe);
    this->points_.id.resize(mesh.nelem * nspe);
    this->points_.surf_point.resize(mesh.nelem * nspe);

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
                            if (zp >= 3.5e0 && zp <= 4.5e0 && xp <= 4.0e0) {
                                this->points_.coord[this->points_.num][0] = xp;
                                this->points_.coord[this->points_.num][1] = yp;
                                this->points_.coord[this->points_.num][2] = zp;
                                this->points_.matid[this->points_.num] = 0;
                                this->points_.id[this->points_.num] = this->points_.num;
                                this->points_.surf_point[this->points_.num] = 0;
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
    this->points_.surf_point.resize(this->points_.num);

    return;
}

void SolidGenerator::WriteBCData(std::ofstream &outfile) {

    this->points_.ubc.BCOutput(outfile, "usbc");
    this->points_.vbc.BCOutput(outfile, "vsbc");
    this->points_.wbc.BCOutput(outfile, "wsbc");

    return;
}

void SolidGenerator::WritePointData(std::ofstream &pointfile) {

    pointfile << std::setw(10) << this->points_.num << "\n";
    OutputVector(pointfile, this->points_.num, this->points_.coord);
    for (int i = 0; i < this->points_.num; i++) {
        pointfile << std::setw(15) << i << std::setw(15) << this->points_.matid[i] << std::setw(15) << this->points_.surf_point[i]
                  << std::setw(15) << this->solid_point_mass_ << std::setw(15) << this->solid_point_volume_ << "\n";
    }

    return;
}

void SolidGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf", this->grid_.mesh);
    WriteVTKHDFPoints("sp.vtkhdf", this->points_);
#endif

    return;
}
