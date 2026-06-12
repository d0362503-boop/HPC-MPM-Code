#include <cmath>
#include <iomanip>
#include <iostream>
#include <array>
#include <vector>

#include "../../dataset.h"
#include "../../mesh.h"
#include "../../shape_function.h"
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

void StabilizedFEM::OutputMeshData(int iview, int istep) {
    std::ofstream ofile;
    ofile.flags(std::ios::right | std::ios::scientific);
    std::string filename = outfile + std::to_string(myrank) + "-" + std::to_string(iview) + "-w.txt";

    ofile.open(filename);
    ofile << std::setw(10) << myrank << std::setw(10) << node << std::setw(10) << istep << std::setw(10) << iview
          << "\n";
    for (int i = 0; i < 3; i++) ofile << std::setw(10) << xynodew[i];
    ofile << "\n";

    for (int i = 0; i < 3; i++) ofile << std::setw(15) << xyminw[i];
    for (int i = 0; i < 3; i++) ofile << std::setw(15) << xymaxw[i];
    ofile << "\n";

    for (int i = 0; i < 3; i++) ofile << std::setw(10) << aelemmin[i];
    for (int i = 0; i < 3; i++) ofile << std::setw(10) << aelemmax[i];
    ofile << "\n";

    for (int n = 0; n < node; n++) {
        ofile << std::setw(15) << this->nvel_vtk[n + nu] << std::setw(15) << this->nvel_vtk[n + nv] << std::setw(15)
              << this->nvel_vtk[n + nw] << std::setw(15) << this->npres_vtk[n] << std::setw(15) << this->nphi_vtk[n]
              << "\n";
    }
    ofile.close();

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
