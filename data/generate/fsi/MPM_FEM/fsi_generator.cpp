#include "data/generate/fsi/MPM_FEM/fsi_generator.h"

#include "data/generate/output_util.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/mesh.h"

std::string MPMFEMFSIGenerator::CaseName() const { return "MPM-FEM FSI"; }

void MPMFEMFSIGenerator::LoadInput() {

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
    for (int i = 0; i < 3; ++i) infile >> xymin[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> xymax[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> xyelem[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> idimc[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> npxye[i];
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> this->rho;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> nmat >> npropmax;
    infile.ignore(1000, '\n');

    VectorAssign(nmat, iprop);
    VectorAssign(nmat, nprop);
    mat_prop.assign(nmat, std::vector<double>(npropmax, 0.0));
    for (int n = 0; n < nmat; ++n) {
        infile.ignore(1000, '\n');
        infile >> iprop[n] >> nprop[n];
        infile.ignore(1000, '\n');

        infile.ignore(1000, '\n');
        for (int i = 0; i < nprop[n]; ++i) { infile >> mat_prop[n][i]; }
        infile.ignore(1000, '\n');
    }

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; ++i) infile >> bb[i];
    infile.ignore(1000, '\n');

    return;
}

void MPMFEMFSIGenerator::CreateBCs() {

    const int xnodec = xynodec[0];
    const int ynodec = xynodec[1];
    const int znodec = xynodec[2];

    this->ubc.nbc.resize(nodec);
    this->vbc.nbc.resize(nodec);
    this->wbc.nbc.resize(nodec);
    this->ubc.fbc.resize(nodec);
    this->vbc.fbc.resize(nodec);
    this->wbc.fbc.resize(nodec);

    this->ubc.ibc = 0;
    this->vbc.ibc = 0;
    this->wbc.ibc = 0;
    for (int k = 0; k < znodec; ++k) {
        for (int j = 0; j < ynodec; ++j) {
            for (int i = 0; i < xnodec; ++i) {
                const int id = i + xnodec * j + xnodec * ynodec * k;
                if (j == 0 || j == ynodec - 1) {
                    this->vbc.nbc[this->vbc.ibc] = id;
                    this->vbc.fbc[this->vbc.ibc] = 0.0e0;
                    this->vbc.ibc++;
                }
            }
        }
    }

    this->CreateFluidBCs();

    return;
}

void MPMFEMFSIGenerator::CreateParticles() {

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

    this->solid_point_volume_ = dx * dy * dz / static_cast<double>(nspe);
    this->solid_point_mass_ = this->rho * this->solid_point_volume_;

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    this->coord.resize(nelem * nspe);
    this->matid.resize(nelem * nspe);
    this->id.resize(nelem * nspe);
    this->surf_point.resize(nelem * nspe);

    this->num = 0;
    for (int k = 0; k < zelem; ++k) {
        for (int j = 0; j < yelem; ++j) {
            for (int i = 0; i < xelem; ++i) {
                const double ecx = dx * (static_cast<double>(i) + 0.5e0) + xmin;
                const double ecy = dy * (static_cast<double>(j) + 0.5e0) + ymin;
                const double ecz = dz * (static_cast<double>(k) + 0.5e0) + zmin;
                for (int kp = 0; kp < npxye[2]; ++kp) {
                    const double zp = ecz + dec2p[kp][2];
                    for (int jp = 0; jp < npxye[1]; ++jp) {
                        const double yp = ecy + dec2p[jp][1];
                        for (int ip = 0; ip < npxye[0]; ++ip) {
                            const double xp = ecx + dec2p[ip][0];
                            if (std::pow(xp - 1.0e-2, 2) + std::pow(zp - 4.0e-2, 2) < std::pow(1.25e-3, 2)) {
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
    this->surf_point.resize(this->num);
    this->MarkSurfacePoints();

    return;
}

void MPMFEMFSIGenerator::WriteBCData(std::ofstream &outfile) {

    this->ubc.BCOutput(outfile, "usbc");
    this->vbc.BCOutput(outfile, "vsbc");
    this->wbc.BCOutput(outfile, "wsbc");

    this->fluid_ubc_.BCOutput(outfile, "uwbc");
    this->fluid_vbc_.BCOutput(outfile, "vwbc");
    this->fluid_wbc_.BCOutput(outfile, "wwbc");
    this->fluid_pbc_.BCOutput(outfile, "hpbc");

    return;
}

void MPMFEMFSIGenerator::WritePointData(std::ofstream &pointfile) {

    pointfile << std::setw(10) << this->num << "\n";
    OutputVector(pointfile, this->num, this->coord);
    for (int i = 0; i < this->num; ++i) {
        pointfile << std::setw(15) << i << std::setw(15) << this->matid[i] << std::setw(15) << this->surf_point[i]
                  << std::setw(15) << this->solid_point_mass_ << std::setw(15) << this->solid_point_volume_ << "\n";
    }

    return;
}

void MPMFEMFSIGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf");
    WriteVTKHDFPoints("sp.vtkhdf", *this);
#endif

    return;
}

void MPMFEMFSIGenerator::CreateFluidBCs() {

    const int xnodec = xynodec[0];
    const int ynodec = xynodec[1];
    const int znodec = xynodec[2];

    this->fluid_ubc_.nbc.resize(nodec);
    this->fluid_vbc_.nbc.resize(nodec);
    this->fluid_wbc_.nbc.resize(nodec);
    this->fluid_pbc_.nbc.resize(nodec);
    this->fluid_ubc_.fbc.resize(nodec);
    this->fluid_vbc_.fbc.resize(nodec);
    this->fluid_wbc_.fbc.resize(nodec);
    this->fluid_pbc_.fbc.resize(nodec);

    this->fluid_ubc_.ibc = 0;
    this->fluid_vbc_.ibc = 0;
    this->fluid_wbc_.ibc = 0;
    this->fluid_pbc_.ibc = 0;
    for (int k = 0; k < znodec; ++k) {
        for (int j = 0; j < ynodec; ++j) {
            for (int i = 0; i < xnodec; ++i) {
                const int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1) {
                    this->fluid_ubc_.nbc[this->fluid_ubc_.ibc] = id;
                    this->fluid_ubc_.fbc[this->fluid_ubc_.ibc] = 0.0e0;
                    this->fluid_ubc_.ibc++;
                    this->fluid_wbc_.nbc[this->fluid_wbc_.ibc] = id;
                    this->fluid_wbc_.fbc[this->fluid_wbc_.ibc] = 0.0e0;
                    this->fluid_wbc_.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    this->fluid_vbc_.nbc[this->fluid_vbc_.ibc] = id;
                    this->fluid_vbc_.fbc[this->fluid_vbc_.ibc] = 0.0e0;
                    this->fluid_vbc_.ibc++;
                }
                if (k == 0 || k == znodec - 1) {
                    this->fluid_ubc_.nbc[this->fluid_ubc_.ibc] = id;
                    this->fluid_ubc_.fbc[this->fluid_ubc_.ibc] = 0.0e0;
                    this->fluid_ubc_.ibc++;
                    this->fluid_wbc_.nbc[this->fluid_wbc_.ibc] = id;
                    this->fluid_wbc_.fbc[this->fluid_wbc_.ibc] = 0.0e0;
                    this->fluid_wbc_.ibc++;
                }
                if (i == 0 && k == znodec - 1) {
                    this->fluid_pbc_.nbc[this->fluid_pbc_.ibc] = id;
                    this->fluid_pbc_.fbc[this->fluid_pbc_.ibc] = 0.0e0;
                    this->fluid_pbc_.ibc++;
                }
            }
        }
    }

    return;
}

void MPMFEMFSIGenerator::MarkSurfacePoints() {

    if (this->num == 0) { return; }

    double upper_surface = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < this->num; ++i) { upper_surface = std::max(upper_surface, this->coord[i][2]); }
    for (int i = 0; i < this->num; ++i) {
        if (std::abs(this->coord[i][2] - upper_surface) <= mtol) { this->surf_point[i] = 1; }
    }

    return;
}
