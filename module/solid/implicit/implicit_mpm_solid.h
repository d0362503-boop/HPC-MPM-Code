#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "module/dataset.h"
#include "module/map_and_interpolate.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/solid_material_point.h"
#include "module/solver/crsmat.h"

namespace implicitmpm {

class ImplicitSolidMPM : public SolidMaterialPointBase {
  public:
    CrsMat SM_;
    std::vector<int> blocks; // target scalar block indices
    int rhs_start;          // first target component index

    ImplicitSolidMPM() {
        this->blocks = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        this->rhs_start = 0;
        this->NR_flag = true;
        this->do_dlb = false;
        this->gamma_nb = 0.5e0;
        this->beta_nb = 0.25e0;
        this->ode_order = 2;
        this->cm_.implicit_flag = true;
        this->SM_.ndof = 3;
        this->SM_.block_row = {0, 0, 0, 1, 1, 1, 2, 2, 2};
        this->SM_.block_col = {0, 1, 2, 0, 1, 2, 0, 1, 2};
        this->SM_.use_petsc = true;
        this->SM_.use_schur_fieldsplit = false;
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

    /**
     * @brief Configure solid BoomerAMG with Euclid smoothing.
     * @param mat Solid linear system being configured.
     * @param pc PETSc preconditioner associated with the solid solver.
     */
    void ConfigurePreconditioner(CrsMat &mat, PC pc) override;

    /**
     * @brief Get target RHS offsets for solid momentum.
     * @return Current control-point offsets for x, y and z momentum.
     */
    std::vector<int> RHSOffsets() const {
        return {nodec * this->rhs_start, nodec * (this->rhs_start + 1), nodec * (this->rhs_start + 2)};
    }

    /** @brief Apply DLB and rebuild the implicit-solid matrix structure. */
    void ApplyDLB() override {

        MaterialPoint::ApplyDLB();

        this->SM_.BuildCrsMat(9);

        return;
    };

    //   protected:
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
     * @param mat Target system, cleared by the caller before assembly.
     * @param naccel_k Nodal acceleration at the current endpoint iterate.
     * @param nvel_k Nodal velocity at the current endpoint iterate.
     * @param stress_k  Particle stress state for tangent assembly.
     */
    virtual void AssembleSystem(CrsMat &mat, const std::vector<double> &naccel_k, //
                                const std::vector<double> &nvel_k,   //
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
    std::array<std::array<double, 3>, 3> ComputeTangentModulus(int pid, int ni, int nj,
                                                               const std::vector<std::array<double, 3>> &dsf,
                                                               const std::array<double, 6> &sts_af);

    /**
     * @brief Apply converged NR displacement increment to nodal displacements.
     */
    void UpdateNRIncrement() override;
};

/**
 * @brief Driver entry point of the implicit ULMPM solid solver.
 */
void SolidImplicitULMPM();

} // namespace implicitmpm
