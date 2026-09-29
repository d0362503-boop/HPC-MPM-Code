#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "module/data_io.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/constitutive_model.h"

class SolidMaterialPointBase : public MaterialPoint {
  public:
    bool Fbar_flag = false;
    // --- Data Initialize ---

    /**
     * @brief Read solid boundary-condition data from an input stream.
     * @param infile Input stream positioned at the solid boundary-condition section.
     */
    void InputBCData(std::ifstream &infile) override;

    /** @brief Rebuild solid control-point BCs for the current DLB region. */
    void RebuildBC() override;

    /**
     * @brief Read the solid particle count, coordinate block and attribute block.
     * @param inflie Stream positioned at count, coordinates, then id/matid/surf_point/mass/vol0 rows.
     */
    void InputPointData(std::ifstream &inflie) override;

    /**
     * @brief Initialize solid particle state (mass, volume, stress, deformation gradient) after input is read.
     */
    void InitializePointData() override;

    /**
     * @brief Read solid restart data from per-rank `*_re.txt` files.
     */
    void RestartInput() override;

    // --- Data I/O ---
    /**
     * @brief Write solid particle data to VTK HDF5 visualization files.
     * @param iview Output view index.
     * @param istep Current time step.
     */
    void OutputPointDataVTKHDF(int iview, int istep) override;

    /**
     * @brief Write solid restart data to per-rank `*_re.txt` files.
     */
    void RestartOutput() override;

    // --- MPI Particle move ---
    /**
     * @brief Move particles that have crossed rank boundaries and exchange them via MPI.
     */
    void MoveParticle() override;

    /**
     * @brief Migrate all solid particle state using the prepared MPI communication plan.
     */
    void MigrateParticleData() override;

    // --- Stress computation ---
    /**
     * @brief Compute the mean normal Cauchy stress of a particle.
     * @param pid Particle whose current stress is sampled.
     * @return One third of the stress trace, with the stored stress sign convention.
     */
    double ComputeMeanStress(int pid) const noexcept {
        double mean_stress = (this->stress[pid][0] + this->stress[pid][1] + this->stress[pid][2]) / 3.0e0;
        return mean_stress;
    }

    /**
     * @brief Compute the von Mises equivalent stress of a particle.
     * @param pid Particle whose current Cauchy stress is sampled.
     * @return Nonnegative equivalent stress in the same units as the stored stress.
     */
    double ComputeVMStress(int pid) const noexcept {
        double VM_stress =
            std::sqrt(0.5e0 * (pow((this->stress[pid][0] - this->stress[pid][1]), 2)   //
                               + pow((this->stress[pid][1] - this->stress[pid][2]), 2) //
                               + pow((this->stress[pid][2] - this->stress[pid][0]), 2)) +
                      3.0e0 * (pow(this->stress[pid][3], 2) + pow(this->stress[pid][4], 2) + pow(this->stress[pid][5], 2)));
        return VM_stress;
    }
    // -------------------------------

    /**
     * @brief Compute a particle's internal-force contribution to one control point.
     * @param ni Local supporting control-point index in the gradient array.
     * @param pid Particle supplying the current integration volume.
     * @param dsf Shape gradients in the stress evaluation configuration.
     * @param stress Cauchy stress in xx, yy, zz, yz, xz, xy order.
     * @return Three momentum residual contributions from minus volume times stress times gradient.
     */
    std::array<double, 3> ComputeInternalForce(int ni, int pid, const std::vector<std::array<double, 3>> &dsf,
                                               const std::array<double, 6> &stress) const noexcept {
        double dsfi1 = dsf[ni][0];
        double dsfi2 = dsf[ni][1];
        double dsfi3 = dsf[ni][2];

        std::array<double, 3> nfint;
        nfint[0] = -this->vol[pid] * (stress[0] * dsfi1 + stress[5] * dsfi2 + stress[4] * dsfi3);
        nfint[1] = -this->vol[pid] * (stress[5] * dsfi1 + stress[1] * dsfi2 + stress[3] * dsfi3);
        nfint[2] = -this->vol[pid] * (stress[4] * dsfi1 + stress[3] * dsfi2 + stress[2] * dsfi3);

        return nfint;
    }

    /**
     * @brief Map particle body and applied forces to one control point.
     * @param pid Particle supplying mass and applied traction-force components.
     * @param sfi Shape-function weight at the receiving control point.
     * @return Weighted body force plus applied force in the three momentum components.
     */
    virtual std::array<double, 3> ComputeExternalForce(int pid, double sfi) const noexcept {
        double fx = bb[0] * facl;
        double fy = bb[1] * facl;
        double fz = bb[2] * facl;

        std::array<double, 3> nfext;
        nfext[0] = sfi * (this->mass[pid] * fx + this->trac_force[pid][0]);
        nfext[1] = sfi * (this->mass[pid] * fy + this->trac_force[pid][1]);
        nfext[2] = sfi * (this->mass[pid] * fz + this->trac_force[pid][2]);

        return nfext;
    }

    /**
     * @brief Update the deformation gradient and its incremental correction for particle `pid`.
     * @param pid             Particle index.
     * @param nenode          Number of element nodes.
     * @param af_coeff        Generalized-α coefficient for the displacement increment.
     * @param ncm             Node IDs of the element supporting the particle.
     * @param sf              Shape-function values.
     * @param dsf             Shape-function gradients.
     * @param delta_def_grad Output incremental deformation-gradient array; entry `pid` is updated.
     * @param def_grad Output total deformation-gradient array; entry `pid` is updated.
     * @note Forms F_trial = delta_F * F_n and writes the member `det_def_grad[pid]`.
     * Passing separate output arrays does not isolate this determinant from the trial state.
     */
    void UpdateDefGrad(int pid, int nenode, double af_coeff, const std::vector<int> &ncm, const std::vector<double> &sf,
                       const std::vector<std::array<double, 3>> &dsf,
                       std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                       std::vector<std::array<std::array<double, 3>, 3>> &def_grad);

    /**
     * @brief Update particle volume from the given determinant of deformation gradient.
     * @param pid Particle index.
     * @param det Determinant of the (possibly F-bar corrected) deformation gradient.
     */
    void UpdateVolume(int pid, double det) { this->vol[pid] = det * this->vol0[pid]; }

    /**
     * @brief Apply prescribed Dirichlet boundary values to a velocity/displacement vector,
     *        including rigid-body constraints.
     * @param nvel Vector to be modified in-place at constrained DOFs.
     */
    void ApplyVelocityBC(std::vector<double> &nvel) override {

        MaterialPoint::ApplyVelocityBC(nvel);

        this->rigid_bc.BCSetZero(nuc, nvel);
        this->rigid_bc.BCSetZero(nvc, nvel);
        this->rigid_bc.BCSetZero(nwc, nvel);

        return;
    }

    /**
     * @brief Zero out a vector at constrained Dirichlet DOFs (typical for acceleration/residual),
     *        including rigid-body constraints.
     * @param naccel Vector to be zeroed in-place at constrained DOFs.
     */
    void ApplyAccelerationBC(std::vector<double> &naccel) override {

        MaterialPoint::ApplyAccelerationBC(naccel);

        this->rigid_bc.BCSetZero(nuc, naccel);
        this->rigid_bc.BCSetZero(nvc, naccel);
        this->rigid_bc.BCSetZero(nwc, naccel);

        return;
    }

    /** @brief Advance the solid nodal solve using the selected explicit or implicit method. */
    virtual void SolveSolid() = 0;

    /**
     * @brief Identify particles whose material ID corresponds to a rigid body and tag them in `rigid_bc`.
     */
    void DetermineRigidBC();

    // --- Constitutive model ---
    ConstitutiveModel cm_;

    /**
     * @brief Update the stress and (if implicit) tangent stiffness for particle `pid`.
     * @param pid              Particle index.
     * @param stress           Particle stress vector (updated in-place).
     * @param det_def_grad_bar F-bar corrected determinant of deformation gradient for each particle.
     * @param def_grad         Total deformation-gradient tensor for each particle.
     * @param delta_def_grad   Incremental deformation-gradient tensor for each particle.
     */
    void UpdateConstitutiveModel(int pid, std::vector<std::array<double, 6>> &stress, const std::vector<double> &det_def_grad_bar,
                                 const std::vector<std::array<std::array<double, 3>, 3>> &def_grad,
                                 const std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad);
    // -------------------------------------
};
