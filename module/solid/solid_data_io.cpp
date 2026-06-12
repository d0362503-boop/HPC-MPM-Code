#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "../dataset.h"
#include "../material_point.h"
#include "../mesh.h"
#include "../mpi_data.h"
#include "solid_material_point.h"

void SolidMaterialPointBase::InputPointData(std::ifstream &infile) {
    infile >> this->num;
    infile.ignore(1000, '\n');

    this->InitializePointData();

    InputVector(infile, this->num, this->coord);

    for (int ip = 0; ip < this->num; ip++) {
        infile >> this->id[ip] >> this->matid[ip] >> this->surf_point[ip] >> this->mass[ip] >> this->vol0[ip];
        infile.ignore(1000, '\n');
    }
    for (int ip = 0; ip < this->num; ip++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this->def_grad[ip][i][j] = eye_mat[i][j];
                if (this->Fbar_flag) this->def_grad_bar[ip][i][j] = eye_mat[i][j];
            }
        }
        this->vol[ip] = this->vol0[ip];
    }

    return;
}

void SolidMaterialPointBase::OutputPointData(int iview, int istep) {
    std::ofstream ofile;
    ofile.flags(std::ios::right | std::ios::scientific);
    std::string filename = outfile + std::to_string(myrank) + "-" + std::to_string(iview) + "-s.txt";
    ofile.open(filename);
    ofile << std::setw(15) << myrank << std::setw(10) << this->num //
          << std::setw(15) << istep << std::setw(10) << iview << "\n";

    for (int pid = 0; pid < this->num; pid++) {
        double mean_stress = this->ComputeMeanStress(pid);
        double VM_stress = this->ComputeVMStress(pid);
        ofile << std::setw(15) << this->id[pid]       //
              << std::setw(15) << this->coord[pid][0] //
              << std::setw(15) << this->coord[pid][1] //
              << std::setw(15) << this->coord[pid][2] //
              << std::setw(15) << this->vel[pid][0]   //
              << std::setw(15) << this->vel[pid][1]   //
              << std::setw(15) << this->vel[pid][2]   //
              << std::setw(15) << VM_stress << "\n";
    }
    ofile.close();

    return;
}

void SolidMaterialPointBase::RestartInput() {
    std::string filename = pointfile + std::to_string(myrank) + "_re.txt";

    std::ifstream reinfile;
    reinfile.open(filename);

    reinfile >> this->num;
    reinfile.ignore(1000, '\n');

    this->InitializePointData();

    if (this->num != 0) {
        InputVector(reinfile, this->num, this->coord);

        for (int ip = 0; ip < this->num; ip++) {
            reinfile >> this->id[ip] >> this->matid[ip]   //
                >> this->surf_point[ip] >> this->mass[ip] //
                >> this->vol0[ip] >> this->vol[ip];
            this->det_def_grad[ip] = this->vol[ip] / this->vol0[ip];
            reinfile.ignore(1000, '\n');
        }

        InputVector(reinfile, this->num, this->vel);
        InputVector(reinfile, this->num, this->accel);
        InputVector(reinfile, this->num, this->stress);
        InputVector(reinfile, this->num, this->def_grad);

        if (this->solswitch == MapScheme::TPIC) {
            InputVector(reinfile, this->num, this->tpic.vel_grad);
            InputVector(reinfile, this->num, this->tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            InputVector(reinfile, this->num, this->apic.vel_Bmat);
            InputVector(reinfile, this->num, this->apic.accel_Bmat);
            InputVector(reinfile, this->num, this->apic.inv_Dmat);
        }
    }

    reinfile.close();

    return;
}

void SolidMaterialPointBase::RestartOutput() {
    std::string filename = pointfile + std::to_string(myrank) + "_re.txt";

    std::ofstream reoutfile;
    reoutfile.flags(std::ios::right | std::ios::scientific);
    reoutfile.open(filename);

    reoutfile << this->num << "     --- Number of solid material point --- " << "\n";

    if (this->num != 0) {
        OutputVector(reoutfile, this->num, this->coord);

        for (int i = 0; i < this->num; i++) {
            reoutfile << this->id[i]                          //
                      << std::setw(15) << this->matid[i]      //
                      << std::setw(15) << this->surf_point[i] //
                      << std::setw(15) << this->mass[i]       //
                      << std::setw(15) << this->vol0[i]       //
                      << std::setw(15) << this->vol[i] << "\n";
        }

        OutputVector(reoutfile, this->num, this->vel);
        OutputVector(reoutfile, this->num, this->accel);
        OutputVector(reoutfile, this->num, this->stress);
        OutputVector(reoutfile, this->num, this->def_grad);

        if (this->solswitch == MapScheme::TPIC) {
            OutputVector(reoutfile, this->num, this->tpic.vel_grad);
            OutputVector(reoutfile, this->num, this->tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            OutputVector(reoutfile, this->num, this->apic.vel_Bmat);
            OutputVector(reoutfile, this->num, this->apic.accel_Bmat);
            OutputVector(reoutfile, this->num, this->apic.inv_Dmat);
        }
    }

    reoutfile.close();

    return;
}