#include "data/generate/solid/solid_generator.h"

#include "data/generate/output_util.h"

#include <array>
#include <cmath>
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

    std::ifstream infile = OpenInputFile("input.txt");

    infile.ignore(1000, '\n');
    std::string solswitch_str;
    infile >> solswitch_str >> this->NR_flag >> rstflag >> nlstep;
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
    infile >> this->rho >> nmat >> npropmax;
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

    return;
}

void SolidGenerator::CreateBCs() {

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
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    this->ubc.nbc[this->ubc.ibc] = id;
                    this->ubc.fbc[this->ubc.ibc] = 0.0e0;
                    this->ubc.ibc++;
                }
                if (k == 0 || j == 0 || j == ynodec - 1) {
                    this->vbc.nbc[this->vbc.ibc] = id;
                    this->vbc.fbc[this->vbc.ibc] = 0.0e0;
                    this->vbc.ibc++;
                }
                if (i == 0 || i == xnodec - 1 || k == 0 || k == znodec - 1) {
                    this->wbc.nbc[this->wbc.ibc] = id;
                    this->wbc.fbc[this->wbc.ibc] = 0.0e0;
                    this->wbc.ibc++;
                }
            }
        }
    }

    return;
}

void SolidGenerator::CreateParticles() {

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
                            if (zp >= 3.5e0 && zp <= 4.5e0 && xp <= 4.0e0) {
                                this->coord[this->num][0] = xp;
                                this->coord[this->num][1] = yp;
                                this->coord[this->num][2] = zp;
                                this->matid[this->num] = 0;
                                this->id[this->num] = this->num;
                                this->surf_point[this->num] = 0;
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

    return;
}

void SolidGenerator::WriteBCData(std::ofstream &outfile) {

    this->ubc.BCOutput(outfile, "usbc");
    this->vbc.BCOutput(outfile, "vsbc");
    this->wbc.BCOutput(outfile, "wsbc");

    return;
}

void SolidGenerator::WritePointData(std::ofstream &pointfile) {

    pointfile << std::setw(10) << this->num << "\n";
    OutputVector(pointfile, this->num, this->coord);
    for (int i = 0; i < this->num; i++) {
        pointfile << std::setw(15) << i << std::setw(15) << this->matid[i] << std::setw(15) << this->surf_point[i]
                  << std::setw(15) << this->solid_point_mass_ << std::setw(15) << this->solid_point_volume_ << "\n";
    }

    return;
}

void SolidGenerator::WriteVisualizationOutputs() {

#ifdef HAVE_HDF5
    std::cout << "Making VTK HDF5 files" << "\n";
    WriteVTKHDFMesh("grid.vtkhdf");
    WriteVTKHDFPoints("sp.vtkhdf", *this);
#endif

    return;
}
