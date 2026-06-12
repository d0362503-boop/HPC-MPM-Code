#include <cmath>
#include <iomanip>
#include <iostream>
#include <array>
#include <vector>

#include "../../dataset.h"
#include "../../mesh.h"
#include "../../shape_function.h"
#include "../../vtk_hdf5.h"
#include "stabilized_fem.h"

using namespace stabilizedfem;

void StabilizedFEM::RestartInput() {
    this->InitializeMeshData();

    std::string filename = gridfile + std::to_string(myrank) + "_re.txt";

    std::ifstream reinfile;
    reinfile.open(filename);

    for (int n = 0; n < nodec; n++) {
        reinfile >> this->naccel[n + nuc] //
            >> this->naccel[n + nvc]      //
            >> this->naccel[n + nwc]      //
            >> this->nvel[n + nuc]        //
            >> this->nvel[n + nvc]        //
            >> this->nvel[n + nwc]        //
            >> this->nvel_old[n + nuc]    //
            >> this->nvel_old[n + nvc]    //
            >> this->nvel_old[n + nwc]    //
            >> this->nvel_older[n + nuc]  //
            >> this->nvel_older[n + nvc]  //
            >> this->nvel_older[n + nwc]  //
            >> this->npres[n] >> this->npres_old[n] >> this->nphi[n];
        reinfile.ignore(1000, '\n');
    }

    reinfile.close();

    return;
}

void StabilizedFEM::RestartOutput() {
    std::string filename = gridfile + std::to_string(myrank) + "_re.txt";

    std::ofstream reoutfile;
    reoutfile.flags(std::ios::right | std::ios::scientific);
    reoutfile.open(filename);

    for (int n = 0; n < nodec; n++) {
        reoutfile << std::setw(15) << this->naccel[n + nuc]     //
                  << std::setw(15) << this->naccel[n + nvc]     //
                  << std::setw(15) << this->naccel[n + nwc]     //
                  << std::setw(15) << this->nvel[n + nuc]       //
                  << std::setw(15) << this->nvel[n + nvc]       //
                  << std::setw(15) << this->nvel[n + nwc]       //
                  << std::setw(15) << this->nvel_old[n + nuc]   //
                  << std::setw(15) << this->nvel_old[n + nvc]   //
                  << std::setw(15) << this->nvel_old[n + nwc]   //
                  << std::setw(15) << this->nvel_older[n + nuc] //
                  << std::setw(15) << this->nvel_older[n + nvc] //
                  << std::setw(15) << this->nvel_older[n + nwc] //
                  << std::setw(15) << this->npres[n]            //
                  << std::setw(15) << this->npres_old[n]        //
                  << std::setw(15) << this->nphi[n] << "\n";
    }

    reoutfile.close();

    return;
}

void StabilizedFEM::OutputMeshDataVTKHDF(int iview, int istep) {
#ifdef HAVE_HDF5
    std::string filename = outfile + std::to_string(myrank) + "-" + std::to_string(iview) + "-w.vtkhdf";

    std::vector<std::array<double, 3>> velocity(node);
    std::vector<std::array<long long, 8>> conn8(nelem);

    for (int n = 0; n < node; ++n) {
        velocity[n] = {this->nvel_vtk[n + nu], this->nvel_vtk[n + nv], this->nvel_vtk[n + nw]};
    }

    for (int m = 0; m < nelem; ++m) {
        for (int n = 0; n < 8; ++n) {
            conn8[m][n] = static_cast<long long>(nc[m][n]);
        }
    }

    vtkhdf::VTKHDFWriter writer(filename, MPI_COMM_SELF);
    auto info = vtkhdf::WriteHexMeshTopology(writer, xyn, conn8);
    writer.SetTime(real_time);

    writer.CreatePointDataGroup();
    writer.WritePointVector("Velocity", info.total_npts, info.local_npts, info.point_global_offset, velocity);
    writer.WritePointScalar("Pressure", info.total_npts, info.local_npts, info.point_global_offset, this->npres_vtk);
    writer.WritePointScalar("Phi", info.total_npts, info.local_npts, info.point_global_offset, this->nphi_vtk);
#endif

    return;
}

void StabilizedFEM::Cp2NodeVTK() {
    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(node, this->nphi_vtk);
    VectorAssign(node, this->npres_vtk);
    VectorAssign(node * 3, this->nvel_vtk);
    for (int m = 0; m < nelem; m++) {
        for (int n = 0; n < 8; n++) {
            int id = nc[m][n];
            std::array<double, 3> xyp = {xyn[id][0], xyn[id][1], xyn[id][2]};
            this->nphi_vtk[id] = 0.0e0;
            this->npres_vtk[id] = 0.0e0;
            this->nvel_vtk[id + nu] = 0.0e0;
            this->nvel_vtk[id + nv] = 0.0e0;
            this->nvel_vtk[id + nw] = 0.0e0;
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                this->nphi_vtk[id] += sfi * this->nphi[nid];
                this->npres_vtk[id] += sfi * this->npres[nid];
                this->nvel_vtk[id + nu] += sfi * this->nvel[nid + nuc];
                this->nvel_vtk[id + nv] += sfi * this->nvel[nid + nvc];
                this->nvel_vtk[id + nw] += sfi * this->nvel[nid + nwc];
            }
        }
    }

    for (int n = 0; n < node; n++) {
        if (this->nphi_vtk[n] < 0.01e0) this->nphi_vtk[n] = 0.0e0;
        if (this->nphi_vtk[n] > 0.99e0) this->nphi_vtk[n] = 1.0e0;
    }

    return;
}
