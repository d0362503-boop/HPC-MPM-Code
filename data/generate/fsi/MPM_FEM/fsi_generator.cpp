#include "data/generate/fsi/MPM_FEM/fsi_generator.h"

#include "data/generate/output_util.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"

void MPMFEMFSISolidGenerator::CreateBCs() {
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
    for (int k = 0; k < znodec; ++k) {
        for (int j = 0; j < ynodec; ++j) {
            for (int i = 0; i < xnodec; ++i) {
                const int id = i + xnodec * j + xnodec * ynodec * k;
                if (j == 0 || j == ynodec - 1) {
                    this->points_.vbc.nbc[this->points_.vbc.ibc] = id;
                    this->points_.vbc.fbc[this->points_.vbc.ibc] = 0.0e0;
                    this->points_.vbc.ibc++;
                }
            }
        }
    }

    return;
}

void MPMFEMFSISolidGenerator::CreateParticles() {
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
    for (int k = 0; k < zelem; ++k) {
        for (int j = 0; j < yelem; ++j) {
            for (int i = 0; i < xelem; ++i) {
                const double ecx = dx * (static_cast<double>(i) + 0.5e0) + xmin;
                const double ecy = dy * (static_cast<double>(j) + 0.5e0) + ymin;
                const double ecz = dz * (static_cast<double>(k) + 0.5e0) + zmin;
                for (int kp = 0; kp < params.npxye[2]; ++kp) {
                    const double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < params.npxye[1]; ++jp) {
                        const double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < params.npxye[0]; ++ip) {
                            const double xp = ecx + dec2p[ip][0];
                            if (std::pow(xp - 1.0e-2, 2) + std::pow(zp - 4.0e-2, 2) < std::pow(1.25e-3, 2)) {
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
    this->points_.surf_point.resize(this->points_.num);
    this->MarkSurfacePoints();

    return;
}

void MPMFEMFSISolidGenerator::MarkSurfacePoints() {
    RunParams &params = this->params_;

    if (this->points_.num == 0) { return; }

    double upper_surface = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < this->points_.num; ++i) { upper_surface = std::max(upper_surface, this->points_.coord[i][2]); }
    for (int i = 0; i < this->points_.num; ++i) {
        if (std::abs(this->points_.coord[i][2] - upper_surface) <= params.mtol) { this->points_.surf_point[i] = 1; }
    }

    return;
}

void MPMFEMFSIFluidGenerator::CreateBCs() {
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

    this->points_.ubc.ibc = 0;
    this->points_.vbc.ibc = 0;
    this->points_.wbc.ibc = 0;
    this->points_.pbc.ibc = 0;
    for (int k = 0; k < znodec; ++k) {
        for (int j = 0; j < ynodec; ++j) {
            for (int i = 0; i < xnodec; ++i) {
                const int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1) {
                    this->points_.ubc.nbc[this->points_.ubc.ibc] = id;
                    this->points_.ubc.fbc[this->points_.ubc.ibc] = 0.0e0;
                    this->points_.ubc.ibc++;
                    this->points_.wbc.nbc[this->points_.wbc.ibc] = id;
                    this->points_.wbc.fbc[this->points_.wbc.ibc] = 0.0e0;
                    this->points_.wbc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    this->points_.vbc.nbc[this->points_.vbc.ibc] = id;
                    this->points_.vbc.fbc[this->points_.vbc.ibc] = 0.0e0;
                    this->points_.vbc.ibc++;
                }
                if (k == 0 || k == znodec - 1) {
                    this->points_.ubc.nbc[this->points_.ubc.ibc] = id;
                    this->points_.ubc.fbc[this->points_.ubc.ibc] = 0.0e0;
                    this->points_.ubc.ibc++;
                    this->points_.wbc.nbc[this->points_.wbc.ibc] = id;
                    this->points_.wbc.fbc[this->points_.wbc.ibc] = 0.0e0;
                    this->points_.wbc.ibc++;
                }
                if (i == 0 && k == znodec - 1) {
                    this->points_.pbc.nbc[this->points_.pbc.ibc] = id;
                    this->points_.pbc.fbc[this->points_.pbc.ibc] = 0.0e0;
                    this->points_.pbc.ibc++;
                }
            }
        }
    }

    return;
}

void MPMFEMFSIFluidGenerator::WriteBCData(std::ofstream &outfile) {

    this->points_.ubc.BCOutput(outfile, "uwbc");
    this->points_.vbc.BCOutput(outfile, "vwbc");
    this->points_.wbc.BCOutput(outfile, "wwbc");
    this->points_.pbc.BCOutput(outfile, "hpbc");

    return;
}

std::string MPMFEMFSIGenerator::CaseName() const { return "MPM-FEM FSI"; }

void MPMFEMFSIGenerator::LoadInput() {
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
    for (int i = 0; i < 3; ++i) infile >> mesh.xymin[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> mesh.xymax[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> mesh.xyelem[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> cmesh.idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> params.npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->solid_.points_.rho;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->solid_.materials_.nmat >> this->solid_.materials_.npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(this->solid_.materials_.nmat, this->solid_.materials_.iprop);
    VectorAssign(this->solid_.materials_.nmat, this->solid_.materials_.nprop);
    this->solid_.materials_.mat_prop.assign(this->solid_.materials_.nmat, std::vector<double>(this->solid_.materials_.npropmax, 0.0));
    for (int n = 0; n < this->solid_.materials_.nmat; ++n) {
        infile.ignore(1000, '\n');
        infile >> this->solid_.materials_.iprop[n] >> this->solid_.materials_.nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < this->solid_.materials_.nprop[n]; ++i) { infile >> this->solid_.materials_.mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> params.bb[i];
    infile.ignore(1000, '\n');

    return;
}

void MPMFEMFSIGenerator::CreateBCs() {

    this->solid_.CreateBCs();
    this->fluid_.CreateBCs();

    return;
}

void MPMFEMFSIGenerator::CreateParticles() {

    this->solid_.CreateParticles();

    return;
}

void MPMFEMFSIGenerator::WriteBCData(std::ofstream &outfile) {

    this->solid_.WriteBCData(outfile);
    this->fluid_.WriteBCData(outfile);

    return;
}

void MPMFEMFSIGenerator::WritePointData(std::ofstream &pointfile) {

    this->solid_.WritePointData(pointfile);

    return;
}

void MPMFEMFSIGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf", this->grid_.mesh);
    WriteVTKHDFPoints("sp.vtkhdf", this->solid_.points_);
#endif

    return;
}
