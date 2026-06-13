#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "../../module/bc.h"
#include "../../module/data_io.h"
#include "../../module/dataset.h"
#include "../../module/fluid/FEM/stabilized_fem.h"
#include "../../module/material_point.h"
#include "../../module/mesh.h"

#ifdef HAVE_HDF5
#include <mpi.h>
#include "../../module/vtk_hdf5.h"
#endif

using namespace std;
using namespace stabilizedfem;

MaterialPoint sp;
StabilizedFEM wfem;

void CreateData();
void MakeSurfacePoint();
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

    cout << " ---- Start making fsi data ----" << "\n";

    CreateData();

    cout << "Making grid file" << "\n";
    ofstream gridfile;
    gridfile.open("fsi/griddata.txt");
    gridfile.flags(ios::right | ios::scientific);
    WriteGlobalMeshHeader(gridfile);
    WriteGlobalBcData(gridfile);
    gridfile.close();

#ifdef HAVE_HDF5
    cout << "Making VTK HDF5 files" << "\n";
    WriteVtkHdf5Mesh("fsi/grid.vtkhdf");
    WriteVtkHdf5Points("fsi/sp.vtkhdf", sp);
#endif

    cout << " ---- Finish making fsi data ----" << "\n";

#ifdef HAVE_HDF5
    MPI_Finalize();
#endif

    return 0;
}

void CreateData() {
    ifstream infile;
    infile.open("fsi/input.txt");

    infile.ignore(1000, '\n');
    string solswitch_str;
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
    infile >> sp.rho >> wfem.rhol >> wfem.rmul >> wfem.rhog >> wfem.rmug >>
        wfem.fs_height;
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
    double msp = sp.rho * volp;

    array<array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    BuildMesh();
    BuildControlPoint();

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
                    double vel = 1.5e0 * (4.0e0 / pow(zmax, 2)) * xyc[id][2] *
                                 (zmax - xyc[id][2]);
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
    sp.id.resize(nelem * nspe);

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
                            if (zp >= 0.21e0 || zp <= 0.19e0) {
                                if (pow((xp - 0.2e0), 2) + pow((zp - 0.2e0), 2) <
                                    pow(0.05e0, 2)) {
                                    sp.coord[sp.num][0] = xp;
                                    sp.coord[sp.num][1] = yp;
                                    sp.coord[sp.num][2] = zp;
                                    sp.matid[sp.num] = 0;
                                    sp.num++;
                                }
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 &&
                                       xp >= 0.15e0 && xp <= 0.25e0) {
                                sp.coord[sp.num][0] = xp;
                                sp.coord[sp.num][1] = yp;
                                sp.coord[sp.num][2] = zp;
                                sp.matid[sp.num] = 0;
                                sp.id[sp.num] = sp.num;
                                sp.num++;
                            } else if (zp >= 0.19e0 && zp <= 0.21e0 &&
                                       xp >= 0.25e0 && xp <= 0.6e0) {
                                sp.coord[sp.num][0] = xp;
                                sp.coord[sp.num][1] = yp;
                                sp.coord[sp.num][2] = zp;
                                sp.matid[sp.num] = 1;
                                sp.id[sp.num] = sp.num;
                                sp.num++;
                            }
                        }
                    }
                }
            }
        }
    }

    sp.coord.resize(sp.num);
    sp.matid.resize(sp.num);

    VectorAssign(sp.num, sp.surf_point);
    MakeSurfacePoint();

    cout << "Making solid point file" << "\n";
    ofstream pointfile;
    pointfile.open("fsi/pointdata.txt");
    pointfile.flags(ios::right | ios::scientific);
    pointfile << setw(10) << sp.num << "\n";
    OutputVector(pointfile, sp.num, sp.coord);
    for (int i = 0; i < sp.num; i++) {
        pointfile << setw(15) << i << setw(15) << sp.matid[i] << setw(15)
                  << sp.surf_point[i] << setw(15) << msp << setw(15) << volp
                  << "\n";
    }
    pointfile.close();

    return;
}

void MakeSurfacePoint() {
    double leftside = numeric_limits<double>::infinity();
    double rightside = -numeric_limits<double>::infinity();
    double bottside = numeric_limits<double>::infinity();
    double upside = -numeric_limits<double>::infinity();

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
        if (abs(z - upside) <= mtol) {
            sp.surf_point[i] = 1;
        }
    }

    return;
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
    sp.ubc.BCOutput(outfile, "usbc");
    sp.vbc.BCOutput(outfile, "vsbc");
    sp.wbc.BCOutput(outfile, "wsbc");

    wfem.ubc.BCOutput(outfile, "uwbc");
    wfem.vbc.BCOutput(outfile, "vwbc");
    wfem.wbc.BCOutput(outfile, "wwbc");
    wfem.pbc.BCOutput(outfile, "hpbc");

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
    writer.WritePointScalar("SurfFlag", info.total_npts, info.local_npts, info.global_offset,
                            point.surf_point);

    return;
}
#endif
