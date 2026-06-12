#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"

#ifdef HAVE_HDF5
#include <mpi.h>
#include "../../module/vtk_hdf5.h"
#endif

using namespace std;

MaterialPoint wp;

void WriteGlobalMeshHeader(ofstream& outfile);
void WriteGlobalBcData(ofstream& outfile);

#ifdef HAVE_HDF5
void WriteVtkHdf5Mesh(const string& filename);
void WriteVtkHdf5Points(const string& filename, const MaterialPoint& point);
#endif

int main() {
#ifdef HAVE_HDF5
    MPI_Init(nullptr, nullptr);
#endif

    cout << " ---- Start making fluid data ----" << "\n";

    ifstream infile;
    infile.open("data/fluid/input.txt");

    infile.ignore(1000, '\n');
    string solswitch_str;
    infile >> solswitch_str >> rstflag >> nlstep;
    wp.solswitch = ParseMapScheme(solswitch_str);
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> ista >> iend >> iout;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    infile >> dt >> mtol >> wp.gamma_nb >> wp.beta_nb;
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
    infile >> wp.rho >> wp.rmu;
    infile.ignore(1000, '\n');

    infile.ignore(1000, '\n');
    for (int i = 0; i < 3; i++) infile >> bb[i];
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
    double mwp = wp.rho * volp;

    array<array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    BuildMesh();

    wp.ubc.nbc.resize(nodec);
    wp.vbc.nbc.resize(nodec);
    wp.wbc.nbc.resize(nodec);
    wp.pbc.nbc.resize(nodec);
    wp.ubc.fbc.resize(nodec);
    wp.vbc.fbc.resize(nodec);
    wp.wbc.fbc.resize(nodec);
    wp.pbc.fbc.resize(nodec);

    wp.ubc.ibc = 0;
    wp.vbc.ibc = 0;
    wp.wbc.ibc = 0;
    wp.pbc.ibc = 0;
    for (int k = 0; k < znodec; k++) {
        for (int j = 0; j < ynodec; j++) {
            for (int i = 0; i < xnodec; i++) {
                int id = i + xnodec * j + xnodec * ynodec * k;
                if (i == 0 || i == xnodec - 1) {
                    wp.ubc.nbc[wp.ubc.ibc] = id;
                    wp.ubc.fbc[wp.ubc.ibc] = 0.0e0;
                    wp.ubc.ibc++;
                }
                if (j == 0 || j == ynodec - 1) {
                    wp.vbc.nbc[wp.vbc.ibc] = id;
                    wp.vbc.fbc[wp.vbc.ibc] = 0.0e0;
                    wp.vbc.ibc++;
                }
                if (k == 0 || k == znodec - 1) {
                    wp.wbc.nbc[wp.wbc.ibc] = id;
                    wp.wbc.fbc[wp.wbc.ibc] = 0.0e0;
                    wp.wbc.ibc++;
                }
            }
        }
    }

    wp.coord.resize(nelem * nspe);
    wp.matid.resize(nelem * nspe);
    wp.id.resize(nelem * nspe);

    wp.num = 0;
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
                                wp.coord[wp.num][0] = xp;
                                wp.coord[wp.num][1] = yp;
                                wp.coord[wp.num][2] = zp;
                                wp.matid[wp.num] = 0;
                                wp.id[wp.num] = wp.num;
                                wp.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    wp.coord.resize(wp.num);
    wp.matid.resize(wp.num);
    wp.id.resize(wp.num);

    cout << "Making grid file" << "\n";
    ofstream gridfile;
    gridfile.open("data/fluid/griddata.txt");
    gridfile.flags(ios::right | ios::scientific);
    WriteGlobalMeshHeader(gridfile);
    WriteGlobalBcData(gridfile);
    gridfile.close();

    cout << "Making water point file" << "\n";
    ofstream pointfile;
    pointfile.open("data/fluid/wpdata.txt");
    pointfile.flags(ios::right | ios::scientific);
    pointfile << setw(10) << wp.num << "\n";
    OutputVector(pointfile, wp.num, wp.coord);
    for (int i = 0; i < wp.num; i++) {
        pointfile << setw(15) << i << setw(15) << wp.matid[i] << setw(15) << mwp
                  << setw(15) << volp << "\n";
    }
    pointfile.close();

#ifdef HAVE_HDF5
    cout << "Making VTK HDF5 files" << "\n";
    WriteVtkHdf5Mesh("data/fluid/grid.vtkhdf");
    WriteVtkHdf5Points("data/fluid/wp.vtkhdf", wp);
#endif

    cout << " ---- Finish making fluid data ----" << "\n";

#ifdef HAVE_HDF5
    MPI_Finalize();
#endif

    return 0;
}

void WriteGlobalMeshHeader(ofstream& outfile) {
    for (int i = 0; i < 3; ++i) outfile << setw(10) << idimc[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << setw(15) << xymin[i];
    outfile << "\n";

    for (int i = 0; i < 3; ++i) outfile << setw(15) << xymax[i];
    outfile << "\n";

    outfile << setw(10) << nelem;
    for (int i = 0; i < 3; ++i) outfile << setw(10) << xyelem[i];
    outfile << "\n";

    outfile << setw(10) << node;
    for (int i = 0; i < 3; ++i) outfile << setw(10) << xynode[i];
    outfile << "\n";

    outfile << setw(10) << nodec;
    for (int i = 0; i < 3; ++i) outfile << setw(10) << xynodec[i];
    outfile << "\n";

    return;
}

void WriteGlobalBcData(ofstream& outfile) {
    wp.ubc.BCOutput(outfile, "uwbc");
    wp.vbc.BCOutput(outfile, "vwbc");
    wp.wbc.BCOutput(outfile, "wwbc");
    wp.pbc.BCOutput(outfile, "hpbc");

    return;
}

#ifdef HAVE_HDF5
void WriteVtkHdf5Mesh(const string& filename) {
    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    vtkhdf::WriteHexMeshTopology(writer, xyn, nc);
    writer.SetTime(0.0e0);

    return;
}

void WriteVtkHdf5Points(const string& filename, const MaterialPoint& point) {
    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    auto info = vtkhdf::WriteParticleTopology(writer, point.coord);
    writer.SetTime(0.0e0);

    writer.CreatePointDataGroup();
    writer.WritePointScalar("ID", info.total_npts, info.local_npts, info.global_offset, point.id);
    writer.WritePointScalar("MatID", info.total_npts, info.local_npts, info.global_offset, point.matid);

    return;
}
#endif
