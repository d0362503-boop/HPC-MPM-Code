#pragma once

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "bc.h"
#include "dataset.h"
#include "map_and_interpolate.h"
#include "mesh.h"
#include "mpi_data.h"
#include "solver/crsmat.h"
#include "solver/solver.h"

class MaterialPoint {
  public:
    // --- MPI part ---
    ParticleCommunication par_comm_;

    /**
     * @brief Determine which MPI rank owns each particle based on its current coordinate.
     */
    void DetermineParticleRank();

    virtual void Moveparticle() {};
    // ------------------------------

    // --- BC set ---
    // --- Velocity or displacement BC ---
    BoundaryCondition ubc, vbc, wbc;
    // --- Pressure BC ---
    BoundaryCondition pbc;
    // --- Rigid BC ---
    BoundaryCondition rigid_bc;
    // --- Inflow BC ---
    BoundaryCondition uinfbc, vinfbc, winfbc;

    /**
     * @brief Apply prescribed Dirichlet boundary values to a velocity/displacement vector.
     * @param nvel Vector to be modified in-place at constrained DOFs.
     */
    virtual void ApplyVelocityBC(std::vector<double> &nvel) {
        this->ubc.BCSetVal(nuc, nvel);
        this->vbc.BCSetVal(nvc, nvel);
        this->wbc.BCSetVal(nwc, nvel);

        return;
    }

    /**
     * @brief Zero out a vector at constrained Dirichlet DOFs (typical for acceleration/residual).
     * @param naccel Vector to be zeroed in-place at constrained DOFs.
     */
    virtual void ApplyAccelerationBC(std::vector<double> &naccel) {
        this->ubc.BCSetZero(nuc, naccel);
        this->vbc.BCSetZero(nvc, naccel);
        this->wbc.BCSetZero(nwc, naccel);

        return;
    }

    /**
     * @brief Divide a nodal vector by a nodal weight with a small-weight cutoff.
     * @param result    Output vector (zero-initialized and overwritten in-place).
     * @param numerator Nodal vector to divide.
     * @param weight    Nodal weight vector used as denominator (e.g. mass or volume).
     * @param offsets   Component offsets into the nodal vectors.
     */
    void CutOffSmallNodalVar(std::vector<double> &result, const std::vector<double> &numerator,
                             const std::vector<double> &weight, const std::vector<int> &offsets) {

        VectorAssign(numerator.size(), result);
        for (int n = 0; n < nodec; n++) {
            if (weight[n] > mtol) {
                for (int offset : offsets) { result[n + offset] = numerator[n + offset] / weight[n]; }
            }
        }

        return;
    }

    /**
     * @brief In-place nodal division by a nodal weight with a small-weight cutoff.
     * @param var     Nodal vector to divide in-place.
     * @param weight  Nodal weight vector used as denominator.
     * @param offsets Component offsets into the nodal vector.
     */
    void CutOffSmallNodalVar(std::vector<double> &var, const std::vector<double> &weight,
                             const std::vector<int> &offsets) {
        for (int n = 0; n < nodec; n++) {
            if (weight[n] > mtol) {
                for (int offset : offsets) { var[n + offset] /= weight[n]; }
            }
        }

        return;
    }

    /**
     * @brief Read boundary-condition data for the current physics object from an input stream.
     * @param infile Input stream positioned at the boundary-condition section.
     */
    virtual void InputBCData(std::ifstream &infile) {};
    // ---------------

    // --- For Implicit MPM ---
    // --- Nonlinear time integration part ---
    // --- Generalized-α part ---
    int ode_order;
    double spec_rad = 1.0e0, alpha_f = 1.0e0, alpha_m = 1.0e0;

    /**
     * @brief Set Generalized-α time-integration parameters from `spec_rad`.
     *
     * Populates `alpha_f`, `alpha_m`, and related Newmark-β parameters in `nb_para`.
     */
    void GeneralizedAlphaParaSet();

    /**
     * @brief Compute the nodal acceleration at the Generalized-α intermediate time level.
     * @return Vector of intermediate accelerations sized to `nodec * 3`.
     */
    std::vector<double> GeneralizedAlphaNodeAccelUpdate() const noexcept;

    // --- Newmark-β part ---
    double gamma_nb, beta_nb;
    std::vector<double> nb_para;

    /**
     * @brief Set Newmark-β time-integration parameters (`gamma_nb`, `beta_nb`).
     */
    void NewmarkBetaParaSet();

    /**
     * @brief Predict velocity and acceleration at the beginning of an implicit time step.
     * @param nvel_k   Nodal velocity vector to be predicted (size `nodec * 3`).
     * @param naccel_k Nodal acceleration vector to be predicted (size `nodec * 3`).
     */
    void PredictNewmarkBetaVelAndAccel(std::vector<double> &nvel_k, std::vector<double> &naccel_k) const noexcept;
    // ------------------------------

    // --- Newton-Raphson ---
    bool NR_flag = true;

    // --- Need to define in each src folder ---
    virtual void BCNRSet() {};

    virtual void BuildPetscBCList(CrsMat &mat) {};

    virtual void BCResidualSet(std::vector<double> &rr) {};
    // ------------------------------------------

    virtual double ComputeNRLumpedMassMat(int pid, double sfi) const noexcept {
        double emd = this->nb_para[3] * sfi * this->mass[pid];

        return emd;
    }

    // --- Default for implicit MPM ---
    virtual void AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel) {
        // --- By default, only inertial forces are calculated ---
        for (int n = 0; n < nodec; n++) {
            // --- For Generalized-α (if α_m = 1, back to Newmark-β) ---
            mat.b_rhs[n + nuc] -= this->nmass[n] * naccel[n + nuc];
            mat.b_rhs[n + nvc] -= this->nmass[n] * naccel[n + nvc];
            mat.b_rhs[n + nwc] -= this->nmass[n] * naccel[n + nwc];
        };

        return;
    }

    virtual void UpdateNRIncrement() {};

    int SolveSystem(CrsMat &mat, int NR_it = -1) const {
        int iter;
        if (mat.use_petsc) {
            mat.AssemblePetscMat(mat.ndof);
            iter = mat.SolveWithPetsc(mat.ndof, NR_it);
        } else {
            iter = GPBiCGSafe(mat);
        }

        return iter;
    }

    /**
     * @brief Commit the converged nodal velocity and acceleration for the time step.
     * @param nvel_k   Converged nodal velocity vector.
     * @param naccel_k Converged nodal acceleration vector.
     */
    void CommitNodalKinematics(const std::vector<double> &nvel_k, const std::vector<double> &naccel_k);

    /**
     * @brief Commit particle kinematics: update velocity via Newmark-β, update position, and apply optional
     * particle-shifting correction.
     * @param accel_old Nodal/particle acceleration from the previous step.
     * @param disp      Nodal displacement increment applied to particle coordinates.
     * @param disp_corr Optional particle-shifting correction displacement.
     */
    void CommitParticleKinematics(const std::vector<std::array<double, 3>> &accel_old,
                                  const std::vector<std::array<double, 3>> &disp,
                                  const std::vector<std::array<double, 3>> &disp_corr = {});

    /**
     * @brief Correct shape-function gradients to refer to the current (deformed) configuration.
     * @param nc     Node IDs of the element supporting the material point.
     * @param nenode Number of nodes in the element.
     * @param dsf    Shape-function gradients in the reference configuration; overwritten in-place.
     */
    void ImplicitDsfCorr(const std::vector<int> &nc, int nenode,
                         std::vector<std::array<double, 3>> &dsf) const noexcept;
    // ------------------------
    // ------------------------

    // --- Particle part ---
    int num;                              // --- Particle number ---
    MapScheme solswitch = MapScheme::PIC; // --- Particle P2G method ---
    double rho, rmu;
    std::vector<int> idepf, idp2p, numep; // --- Particle mesh linklist ---
    // ----- Material point variable -----
    std::vector<int> id, matid, surf_point; // id: point id, matid: material id, surf_point: surface point flag
    std::vector<double> pres, phi, mass, vol, vol0, det_def_grad, det_def_grad_bar;
    std::vector<std::array<double, 3>> coord, vel, accel, displ, trac_force;
    std::vector<std::array<double, 6>> stress;
    std::vector<std::vector<double>> ustatev;
    std::array<std::array<double, 6>, 6> stif_mat;
    std::vector<std::array<std::array<double, 3>, 3>> def_grad, def_grad_bar, delta_def_grad_bar;

    virtual void InputPointData(std::ifstream &inflie) {};

    virtual void InitializePointData() {};

    virtual void RestartInput() {};

    virtual void OutputPointDataVTKHDF(int iview, int istep) {};

    virtual void RestartOutput() {};

    bool infbc_isfilled = false;
    std::vector<int> inflow_row;

    /**
     * @brief Initialize per-boundary inflow layer counters.
     *
     * Resizes `inflow_row` to the total number of inflow boundary faces across
     * all directions and resets all counters to zero.
     */
    void InitializeInflowBC() {

        int num = (this->uinfbc.ibc + this->vinfbc.ibc + this->winfbc.ibc);
        this->inflow_row.assign(num, 0);

        return;
    }

    /**
     * @brief Check whether any inflow boundary cell is sufficiently filled.
     * @param infbc Inflow boundary condition describing the active boundary faces.
     * @return True if at least one marked boundary cell has a fill ratio above
     *         the tolerance (0.95), false otherwise.
     */
    bool InflowMeshisFilled(const BoundaryCondition &infbc);

    /**
     * @brief Determine the inward/outward sign for an inflow boundary cell.
     * @param dir     Inflow direction index (0=x, 1=y, 2=z).
     * @param mesh_id Local element index along direction `dir`.
     * @return +1 for the lower domain boundary, -1 for the upper domain boundary.
     */
    int GetInflowSign(int dir, int mesh_id) {

        int sign;
        int mid = mesh_id + aelemmin[dir];
        if (mid == 0) {
            sign = 1;
        } else if (mid == xyelemw[dir] - 1) {
            sign = -1;
        }

        return sign;
    };

    /**
     * @brief Dispatch inflow particle generation for a given direction.
     * @param dir   Inflow direction index (0=x, 1=y, 2=z).
     * @param ifp   Inflow particle buffer to populate.
     * @param infbc Inflow boundary condition for direction `dir`.
     *
     * Switches between the empty-mesh and filled-mesh generators based on the
     * current fill state of the boundary cells.
     */
    void GenerateInflowParticles(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc);

    /**
     * @brief Assign globally unique IDs to newly generated inflow particles.
     * @param ifp Inflow particle buffer whose `id` entries will be overwritten.
     *
     * Computes the current global maximum particle ID across all MPI ranks, then
     * uses an MPI exclusive scan to assign contiguous unique IDs to the new
     * inflow particles.
     */
    void AssignUniqueInflowIds(MaterialPoint &ifp) const;

    virtual void InflowParticles() {};

    virtual void GenerateInflowParticlesEmptyMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) {};

    virtual void GenerateInflowParticlesFilledMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) {};

    /**
     * @brief Build element-to-particle linked lists (`idepf`, `idp2p`, `numep`).
     */
    void MeshPointLinklist();

    /**
     * @brief Generate material points inside each background element using Gaussian quadrature.
     */
    void BuildGaussianPoint();

    /**
     * @brief Compute a unit normal vector for each surface particle.
     */
    void CalPointUnitNormal();

    /**
     * @brief Compute the incremental deformation gradient at a material point.
     * @param nc        Node IDs of the element supporting the material point.
     * @param nenode    Number of nodes in the element.
     * @param af_coeff  Generalized-α interpolation coefficient for the displacement increment.
     * @param dsf       Shape-function gradients in the current configuration.
     * @return Incremental deformation gradient tensor.
     */
    virtual std::array<std::array<double, 3>, 3>
    ComputeDeltaDefGrad(const std::vector<int> &nc, int nenode, double af_coeff,
                        const std::vector<std::array<double, 3>> &dsf) const noexcept;

    /**
     * @brief Compute a particle-shifting correction to regularize particle distribution.
     *
     * Based on the particle shifting scheme in Chandra et al.,
     * "Stabilized mixed material point method for incompressible fluid flow analysis,"
     * Computer Methods in Applied Mechanics and Engineering, 419, 116644, 2024.
     *
     * @return Correction displacement vector for each particle.
     */
    std::vector<std::array<double, 3>> DeltaCorrectionParticleShifting() const;

    /**
     * @brief Compute a pairwise repulsive particle-shifting correction.
     *
     * Based on the isotropic position correction in Ando et al.,
     * "Preserving Fluid Sheets with Adaptively Sampled Anisotropic Particles,"
     * IEEE Transactions on Visualization and Computer Graphics, 2012,
     * Sec. 5.2, Eq. (18).
     *
     * @return Correction displacement vector for each particle.
     */
    std::vector<std::array<double, 3>> PairwiseRepulsiveParticleShifting();
    // ----------------------------------

    // ----- Control point variable -----
    std::vector<double> nmass, nvel_old, nvel_older, nvof, nmome, nvel, //
        ndispl, npres, npres_old, nphi, nnormal,                        //
        naccel, nforce, nvel_vtk, npres_vtk, nphi_vtk;

    // --- Mapping & Interpolation scheme ---
    FLIP flip;
    PIC pic;
    APIC apic;
    TPIC tpic;

    virtual void Particle2Node() {};

    void PICFamilyVelReset() {
        if (this->solswitch != MapScheme::FLIP) { VectorAssign(this->num, this->vel); }
        if (this->solswitch == MapScheme::TPIC) {
            VectorAssign(this->num, this->tpic.vel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            VectorAssign(this->num, this->apic.vel_Bmat);
        }

        return;
    }

    /**
     * @brief Map particle velocity to the nodal momentum vector `nmome`.
     * @param pid Particle index.
     * @param nid Control-point index.
     * @param sfi Shape-function value at the control point.
     */
    void VelP2G(int pid, int nid, double sfi);

    /**
     * @brief Gather velocity from control points to a particle using FLIP/PIC/TPIC/APIC.
     * @param pid Particle index.
     * @param ni  Local node index within the element.
     * @param nid Global control-point index.
     * @param sfi Shape-function value at the control point.
     * @param dsf Shape-function gradients at the control point (used by TPIC/APIC).
     */
    void PICFamilyVelG2P(int pid, int ni, int nid, double sfi, //
                         const std::vector<std::array<double, 3>> &dsf);

    void PICFamilyAccelReset() {
        VectorAssign(this->num, this->accel);
        if (this->solswitch == MapScheme::TPIC) {
            VectorAssign(this->num, this->tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            VectorAssign(this->num, this->apic.accel_Bmat);
        }

        return;
    };

    /**
     * @brief Map particle acceleration to the nodal force vector `nforce`.
     * @param pid Particle index.
     * @param nid Control-point index.
     * @param sfi Shape-function value at the control point.
     */
    void AccelP2G(int pid, int nid, double sfi);

    /**
     * @brief Gather acceleration from control points to a particle using FLIP/PIC/TPIC/APIC.
     * @param pid Particle index.
     * @param ni  Local node index within the element.
     * @param nid Global control-point index.
     * @param sfi Shape-function value at the control point.
     * @param dsf Shape-function gradients at the control point (used by TPIC/APIC).
     */
    void PICFamilyAccelG2P(int pid, int ni, int nid, double sfi, //
                           const std::vector<std::array<double, 3>> &dsf);

    virtual void Node2Particle() {};

    /**
     * @brief Map a per-particle scalar variable to control points without mass weighting.
     * @param pid       Particle index.
     * @param nid       Global control-point index.
     * @param sfi       Shape-function value at the control point.
     * @param point_var Per-particle scalar values.
     * @param node_var  Nodal/control-point accumulator vector.
     */
    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_var,
                        std::vector<T> &node_var) const noexcept;

    /**
     * @brief Map a per-particle scalar variable to control points with mass weighting.
     * @param pid        Particle index.
     * @param nid        Global control-point index.
     * @param sfi        Shape-function value at the control point.
     * @param point_mass Per-particle mass values.
     * @param point_var  Per-particle scalar values.
     * @param node_var   Nodal/control-point accumulator vector.
     */
    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass, const std::vector<T> &point_var,
                        std::vector<T> &node_var) const noexcept;

    /**
     * @brief Map a per-particle vector variable (one component) to control points with mass weighting.
     * @param pid        Particle index.
     * @param nid        Global control-point index.
     * @param sfi        Shape-function value at the control point.
     * @param point_mass Per-particle mass values.
     * @param point_var  Per-particle vector values.
     * @param node_var   Nodal/control-point accumulator vector (size `nodec * 3`).
     */
    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                        const std::vector<std::array<T, 3>> &point_var, std::vector<T> &node_var) const noexcept;
    // -----------------------------------
    virtual ~MaterialPoint() = default;
};

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_var,
                                   std::vector<T> &node_var) const noexcept {
    node_var[nid] += sfi * point_var[pid];

    return;
}

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                                   const std::vector<T> &point_var, std::vector<T> &node_var) const noexcept {
    node_var[nid] += sfi * point_mass[pid] * point_var[pid];

    return;
}

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                                   const std::vector<std::array<T, 3>> &point_var,
                                   std::vector<T> &node_var) const noexcept {
    node_var[nid + nuc] += sfi * point_mass[pid] * point_var[pid][0];
    node_var[nid + nvc] += sfi * point_mass[pid] * point_var[pid][1];
    node_var[nid + nwc] += sfi * point_mass[pid] * point_var[pid][2];

    return;
}
