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
     * @brief Read and initialize the standalone implicit-solid case input data.
     *
     * Loads global control parameters, mesh/time derived quantities, boundary
     * conditions, and particle data for the implicit solid solver path.
     */
    void DataInput();

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
     * @brief Map solid particle mass/momentum/force to control points (P2G).
     */
    void Particle2Node() override;

    /**
     * @brief Map updated nodal kinematics back to solid particles (G2P).
     */
    void Node2Particle() override;

    /**
     * @brief Driver for one implicit solid time step: reset, Newton–Raphson loop, solve, update.
     */
    void SolveSolid() override;

    /**
     * @brief Update the deformation gradient and its incremental correction for particle `pid`.
     * @param pid             Particle index.
     * @param nenode          Number of element nodes.
     * @param af_coeff        Generalized-α coefficient for the displacement increment.
     * @param ncm             Node IDs of the element supporting the particle.
     * @param sf              Shape-function values.
     * @param dsf             Shape-function gradients.
     * @param delta_def_grad  Output incremental deformation-gradient tensor for each node.
     * @param def_grad        Output total deformation-gradient tensor for each node.
     */
    void UpdateDefGrad(int pid, int nenode, double af_coeff, const std::vector<int> &ncm, const std::vector<double> &sf,
                       const std::vector<std::array<double, 3>> &dsf,
                       std::vector<std::array<std::array<double, 3>, 3>> &delta_def_grad,
                       std::vector<std::array<std::array<double, 3>, 3>> &def_grad);

    /**
     * @brief Initialize stress state at the start of the Newton–Raphson loop.
     * @return Tuple or structured result holding initial stress and related data.
     */
    auto InitializeNRStress();

    void BCNRSet() override {
        this->ubc.BCSetDt(nuc, this->ndispl);
        this->vbc.BCSetDt(nvc, this->ndispl);
        this->wbc.BCSetDt(nwc, this->ndispl);
        this->rigid_bc.BCSetDt(nuc, this->ndispl);
        this->rigid_bc.BCSetDt(nvc, this->ndispl);
        this->rigid_bc.BCSetDt(nwc, this->ndispl);

        return;
    };

    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(nuc, rr);
        this->vbc.BCSetZero(nvc, rr);
        this->wbc.BCSetZero(nwc, rr);
        this->rigid_bc.BCSetZero(nuc, rr);
        this->rigid_bc.BCSetZero(nvc, rr);
        this->rigid_bc.BCSetZero(nwc, rr);

        return;
    };

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
     * @brief Assemble the implicit solid tangent matrix and residual vector.
     * @param naccel_k Nodal acceleration vector at the intermediate time level.
     * @param nvel_k   Nodal velocity vector at the intermediate time level.
     * @param stress_k Particle stress state used for the tangent assembly.
     */
    void AssembleSystem(const std::vector<double> &naccel_k, const std::vector<double> &nvel_k,
                        std::vector<std::array<double, 6>> &stress_k);

    /**
     * @brief Compute the tangent-modulus contribution of particle `pid` for nodes `ni` and `nj`.
     * @param pid    Particle index.
     * @param ni     First local node index.
     * @param nj     Second local node index.
     * @param dsf    Shape-function gradients.
     * @param sts_af Particle stress at the intermediate time level.
     * @return Tangent stiffness scalar contribution.
     */
    auto ComputeTangentModulus(int pid, int ni, int nj, const std::vector<std::array<double, 3>> &dsf,
                               const std::array<double, 6> &sts_af);

    /**
     * @brief Apply the converged Newton–Raphson displacement increment to nodal displacements.
     */
    void UpdateNRIncrement() override;
};
} // namespace implicitmpm
