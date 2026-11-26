#pragma once

#include <vector>

#include "module/dataset.h"
#include "module/implicit/implicit_mpm.h"
#include "module/implicit/newton_solver.h"
#include "module/mesh.h"
#include "module/solid/solid_material_point.h"
#include "module/solver/crsmat.h"

namespace implicitmpm {

class ImplicitSolidMPM : public SolidMaterialPointBase, public NonlinearProblem {
  public:
    CrsMat SM_;
    ImplicitMPM implicit_;        // implicit time integration and kinematics
    bool NR_flag;                 // Newton iterations; false for linear elasticity
    int max_NR_it;                // last Newton iteration index, set per step from NR_flag
    std::vector<int> blocks;      // target scalar block indices
    int rhs_start;           // first target component index

    /**
     * @brief Initialize the three-component implicit solid system and its PETSc owner.
     * @param grid   Grid of this rank (background mesh, control points and partition).
     * @param params Run parameters of this run.
     */
    ImplicitSolidMPM(const Grid &grid, const RunParams &params) : SolidMaterialPointBase(grid, params), SM_(grid, params) {
        this->blocks = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        this->rhs_start = 0;
        this->NR_flag = true;
        this->do_dlb = false;
        this->implicit_.integrator_.gamma_nb = 0.5e0;
        this->implicit_.integrator_.beta_nb = 0.25e0;
        this->implicit_.integrator_.ode_order = 2;
        this->cm_.implicit_flag = true;
        this->SM_.ndof = 3;
        this->SM_.block_row = {0, 0, 0, 1, 1, 1, 2, 2, 2};
        this->SM_.block_col = {0, 1, 2, 0, 1, 2, 0, 1, 2};
        this->SM_.use_petsc = true;
        this->SM_.FEM_flag = false;
        this->SM_.rebuild_freq = 1; // refresh each step start
        this->SM_.owner_ = this;
    }

    /**
     * @brief Read and initialize standalone implicit-solid input data.
     * @param grid   Grid of this rank; mesh, control-point and partition data are read into it.
     * @param params Run parameters; filled from file.dat and the parameter file.
     */
    void DataInput(Grid &grid, RunParams &params);

    /**
     * @brief Map solid particle mass/volume/momentum/force to control points (P2G).
     * @param ts Current time state of this step (load factor and physical time).
     */
    void Particle2Node(const TimeState &ts) override;

    /**
     * @brief Map updated nodal kinematics back to solid particles (G2P).
     * @param ts Current time state of this step (load factor and physical time).
     */
    void Node2Particle(const TimeState &ts) override;

    /**
     * @brief Driver for one implicit solid time step.
     * @param ts Current time state of this step (load factor and physical time).
     */
    void SolveSolid(const TimeState &ts) override;

    // Newton steps driven by SolveNewton().
    /** @brief Reset the displacement increment and Newton scratch state. */
    void BeginNewtonStep() override;

    /**
     * @brief Apply iterate BCs and kinematics, then assemble the solid tangent and residual.
     * @param ts Current time state; its load factor `facl` scales BCs and body forces.
     */
    void AssembleNewtonSystem(const TimeState &ts) override;

    /**
     * @brief Check the solid residual; linear elasticity only reports and never converges early.
     * @param NR_it Current Newton iteration.
     * @param ksp_its Krylov iterations of the previous correction.
     * @return Whether the Newton iteration has converged.
     */
    bool CheckNewtonConvergence(int NR_it, int ksp_its) override;

    /**
     * @brief Solve the solid system for the displacement correction.
     * @param NR_it Current Newton iteration.
     * @param ts Current time state; its step index sets the preconditioner rebuild cadence.
     * @return Krylov iterations used.
     */
    int SolveNewtonCorrection(int NR_it, const TimeState &ts) override;

    /** @brief Release the solid linear-system storage. */
    void EndNewtonStep() override;

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
        return {this->grid_.cmesh.nodec * this->rhs_start, this->grid_.cmesh.nodec * (this->rhs_start + 1), this->grid_.cmesh.nodec * (this->rhs_start + 2)};
    }

    /**
     * @brief Apply DLB and rebuild the implicit-solid matrix structure.
     * @param grid Writable view of the grid this object observes; repartitioned in place.
     */
    void ApplyDLB(Grid &grid) override {

        MaterialPoint::ApplyDLB(grid);

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
     * @param ts Current time state; its load factor `facl` scales the prescribed values.
     */
    void BCNRSet(const TimeState &ts) {
        this->ApplyDisplacementBC(this->ndispl, ts);

        return;
    };

    /**
     * @brief Zero constrained DOFs in residual vector.
     * @param rr Residual vector.
     */
    void BCResidualSet(std::vector<double> &rr) override {
        this->ubc.BCSetZero(this->grid_.cmesh.nuc, rr);
        this->vbc.BCSetZero(this->grid_.cmesh.nvc, rr);
        this->wbc.BCSetZero(this->grid_.cmesh.nwc, rr);
        this->rigid_bc.BCSetZero(this->grid_.cmesh.nuc, rr);
        this->rigid_bc.BCSetZero(this->grid_.cmesh.nvc, rr);
        this->rigid_bc.BCSetZero(this->grid_.cmesh.nwc, rr);

        return;
    };

    /**
     * @brief Register constrained DOFs with PETSc matrix.
     * @param mat PETSc matrix.
     */
    void BuildPetscBCList(CrsMat &mat) override {
        mat.AddBCComponent(this->ubc, this->grid_.cmesh.nuc);
        mat.AddBCComponent(this->vbc, this->grid_.cmesh.nvc);
        mat.AddBCComponent(this->wbc, this->grid_.cmesh.nwc);
        mat.AddBCComponent(this->rigid_bc, this->grid_.cmesh.nuc);
        mat.AddBCComponent(this->rigid_bc, this->grid_.cmesh.nvc);
        mat.AddBCComponent(this->rigid_bc, this->grid_.cmesh.nwc);

        return;
    };

    /**
     * @brief Assemble tangent matrix and residual vector.
     * @param mat Target system, cleared by the caller before assembly.
     * @param naccel_k Nodal acceleration at the current endpoint iterate.
     * @param stress_k  Particle stress state for tangent assembly.
     * @param ts        Current time state (load factor for body forces).
     */
    virtual void AssembleSystem(CrsMat &mat, const std::vector<double> &naccel_k, //
                                std::vector<std::array<double, 6>> &stress_k, const TimeState &ts);

    /**
     * @brief Tangent-modulus contribution for one particle-node pair.
     * @param ni      First local node index.
     * @param nj      Second local node index.
     * @param dsf     Shape-function gradients.
     * @param sts_af  Particle stress at intermediate time level.
     * @return Three-by-three material and geometric stiffness block before volume and alpha_f weighting.
     */
    std::array<std::array<double, 3>, 3> ComputeTangentModulus(int ni, int nj, const std::vector<std::array<double, 3>> &dsf,
                                                               const std::array<double, 6> &sts_af);

    /**
     * @brief Add the current linear correction to the endpoint nodal displacement increment.
     */
    void UpdateNRIncrement() override;

    /**
     * @brief Lumped-mass contribution of one particle to the Newton tangent.
     * @param pid Particle index.
     * @param sfi Shape-function value at the control point.
     * @return Diagonal inertia term; FSI solids add their added-mass contribution.
     */
    virtual double ComputeNRLumpedMassMat(int pid, double sfi) const noexcept {

        return this->implicit_.LumpedMass(*this, pid, sfi);
    }

    /**
     * @brief Add nodal inertial forces to the implicit system RHS.
     * @param mat System receiving the inertial force residual.
     * @param naccel Nodal acceleration at the inertia evaluation time level.
     * @param offsets RHS offsets for the x, y and z momentum components.
     */
    virtual void AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel, const std::vector<int> &offsets) {

        this->implicit_.AddInertialForce(*this, mat, naccel, offsets);

        return;
    }

    /**
     * @brief Correct shape-function gradients to refer to the current (deformed) configuration.
     * @param nc     Node IDs of the element supporting the material point.
     * @param nenode Number of nodes in the element.
     * @param dsf    Shape-function gradients in the reference configuration; overwritten in-place.
     */
    void ImplicitDsfCorr(const std::vector<int> &nc, int nenode, std::vector<std::array<double, 3>> &dsf) const noexcept;

  private:
    // Newton scratch state of the current step.
    std::vector<double> nvel_k_, naccel_k_;           // nodal velocity and acceleration at the iterate
    std::vector<std::array<double, 6>> stress_k_;     // particle stress snapshot for tangent assembly
    double r0r_;                                      // initial residual reference
};

} // namespace implicitmpm
