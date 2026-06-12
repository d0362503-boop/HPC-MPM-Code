#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "../data_io.h"
#include "../dataset.h"
#include "../material_point.h"
#include "../mesh.h"
#include "../mpi_data.h"
#include "constitutive_model.h"

class SolidMaterialPointBase : public MaterialPoint {
  public:
    bool Fbar_flag = false;
    // --- Data Initialize ---
    void InputPointData(std::ifstream &inflie) override;

    void InitializePointData() override;

    void RestartInput() override;

    // --- Data I/O ---
    void OutputPointDataVTKHDF(int iview, int istep) override;

    void RestartOutput() override;

    // --- MPI Particle move ---
    void Moveparticle() override;

    // --- Stress computation ---
    double ComputeMeanStress(int pid) {
        double mean_stress = (this->stress[pid][0] + this->stress[pid][1] + this->stress[pid][2]) / 3.0e0;
        return mean_stress;
    }

    double ComputeVMStress(int pid) {
        double VM_stress = std::sqrt(
            0.5e0 * (pow((this->stress[pid][0] - this->stress[pid][1]), 2)   //
                     + pow((this->stress[pid][1] - this->stress[pid][2]), 2) //
                     + pow((this->stress[pid][2] - this->stress[pid][0]), 2)) +
            3.0e0 * (pow(this->stress[pid][3], 2) + pow(this->stress[pid][4], 2) + pow(this->stress[pid][5], 2)));
        return VM_stress;
    }
    // -------------------------------

    std::array<double, 3> ComputeInternalForce(int ni, int pid, const std::vector<std::array<double, 3>> &dsf,
                                               const std::array<double, 6> &stress) {
        double dsfi1 = dsf[ni][0];
        double dsfi2 = dsf[ni][1];
        double dsfi3 = dsf[ni][2];

        std::array<double, 3> nfint;
        nfint[0] = -this->vol[pid] * (stress[0] * dsfi1 + stress[5] * dsfi2 + stress[4] * dsfi3);
        nfint[1] = -this->vol[pid] * (stress[5] * dsfi1 + stress[1] * dsfi2 + stress[3] * dsfi3);
        nfint[2] = -this->vol[pid] * (stress[4] * dsfi1 + stress[3] * dsfi2 + stress[2] * dsfi3);

        return nfint;
    }

    virtual std::array<double, 3> ComputeExternalForce(int pid, double sfi) {
        double fx = bb[0] * facl;
        double fy = bb[1] * facl;
        double fz = bb[2] * facl;

        std::array<double, 3> nfext;
        nfext[0] = sfi * (this->mass[pid] * fx + this->trac_force[pid][0]);
        nfext[1] = sfi * (this->mass[pid] * fy + this->trac_force[pid][1]);
        nfext[2] = sfi * (this->mass[pid] * fz + this->trac_force[pid][2]);

        return nfext;
    }

    virtual void SolveSolid() = 0;

    void DetermineRigidBC();

    // --- Constitutive model ---
    ConstitutiveModel cm_;

    void UpdateConstitutiveModel(int pid, std::vector<std::array<double, 6>> &stress,
                                 const std::vector<double> &det_def_grad_bar,
                                 const std::vector<std::array<std::array<double, 3>, 3>> &def_grad,
                                 const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad);
    // -------------------------------------
};
