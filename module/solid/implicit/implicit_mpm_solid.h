#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../solver/crsmat.h"
#include "../solid_material_point.h"

namespace implicitmpm {

class ImplicitSolidMPM : public SolidMaterialPointBase {
  public:
    CrsMat SM_;

    ImplicitSolidMPM() {
        this->ode_order = 2;
        this->cm_.implicit_flag = true;
        this->SM_.ndof = 3;
        this->SM_.use_petsc = true;
        this->SM_.FEM_flag = false;
        this->SM_.amg_rebuild_freq = 1; // keep solid AMG alive across this 10-step window unless iterations deteriorate
        this->SM_.owner_ = this;
    }

    /**
     * @brief Read and initialize standalone implicit-solid input data.
     */
    void DataInput();

    /**
     * @brief Map solid particle mass/volume/momentum/force to control points (P2G).
     */
    void Particle2Node() override;

    /**
     * @brief Map updated nodal kinematics back to solid particles (G2P).
     */
    void Node2Particle() override;

    /**
     * @brief Driver for one implicit solid time step.
     */
    void SolveSolid() override;

  private:
    // Case-specific surface traction force.
    void SetTracForce() {
        const double t0 = 1.0e0 * dxy[1] * dxy[2] / (npxye[1] * npxye[2]);
        double t = 0.0e0;
        if (real_time < 5.0e-3) { t = -t0; }

        for (int ip = 0; ip < this->num; ip++) {
            if (this->surf_point[ip] == 1) { this->trac_force[ip][0] = t; }
        }

        return;
    };

    /**
     * @brief Stress snapshot at the start of the NR loop.
     * @return Particle stress vector.
     */
    auto InitializeNRStress() { return this->stress; };

    /**
     * @brief Apply Dirichlet displacement increments for current NR iteration.
     */
    void BCNRSet() override {
        this->ubc.BCSetDt(nuc, this->ndispl);
        this->vbc.BCSetDt(nvc, this->ndispl);
        this->wbc.BCSetDt(nwc, this->ndispl);
        this->rigid_bc.BCSetDt(nuc, this->ndispl);
        this->rigid_bc.BCSetDt(nvc, this->ndispl);
        this->rigid_bc.BCSetDt(nwc, this->ndispl);

        return;
    };

    /**
     * @brief Zero constrained DOFs in residual vector.
     * @param rr Residual vector.
     */
    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(nuc, rr);
        this->vbc.BCSetZero(nvc, rr);
        this->wbc.BCSetZero(nwc, rr);
        this->rigid_bc.BCSetZero(nuc, rr);
        this->rigid_bc.BCSetZero(nvc, rr);
        this->rigid_bc.BCSetZero(nwc, rr);

        return;
    };

    /**
     * @brief Register constrained DOFs with PETSc matrix.
     * @param mat PETSc matrix.
     */
    void BuildPetscBCList(CrsMat &mat) override {
        mat.AddBCComponent(this->ubc, nuc);
        mat.AddBCComponent(this->vbc, nvc);
        mat.AddBCComponent(this->wbc, nwc);
        mat.AddBCComponent(this->rigid_bc, nuc);
        mat.AddBCComponent(this->rigid_bc, nvc);
        mat.AddBCComponent(this->rigid_bc, nwc);

        return;
    };

    /**
     * @brief Assemble tangent matrix and residual vector.
     * @param naccel_k  Nodal acceleration at intermediate time level.
     * @param nvel_k    Nodal velocity at intermediate time level.
     * @param stress_k  Particle stress state for tangent assembly.
     */
    void AssembleSystem(const std::vector<double> &naccel_k, const std::vector<double> &nvel_k,
                        std::vector<std::array<double, 6>> &stress_k);

    /**
     * @brief Tangent-modulus contribution for one particle-node pair.
     * @param pid     Particle index.
     * @param ni      First local node index.
     * @param nj      Second local node index.
     * @param dsf     Shape-function gradients.
     * @param sts_af  Particle stress at intermediate time level.
     * @return Tangent stiffness scalar contribution.
     */
    auto ComputeTangentModulus(int pid, int ni, int nj, const std::vector<std::array<double, 3>> &dsf,
                               const std::array<double, 6> &sts_af);

    /**
     * @brief Apply converged NR displacement increment to nodal displacements.
     */
    void UpdateNRIncrement() override;
};

} // namespace implicitmpm
