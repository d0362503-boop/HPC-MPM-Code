#pragma once

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "module/bc.h"
#include "module/dataset.h"
#include "module/map_and_interpolate.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solver/crsmat.h"

namespace mpm_dlb {
struct Region;
}

class MaterialPoint {
  public:
    // --- MPI part ---
    ParticleCommunication par_comm_;
    bool do_dlb = false; // enable dynamic load balance

    /**
     * @brief Build a normal particle-migration plan from the current MPI regions.
     * @param regions One inclusive background-element region per MPI rank, ordered by rank.
     */
    void DetermineParticleRank(const std::vector<mpm_dlb::Region> &regions);

    /**
     * @brief Build a DLB particle-migration plan for a newly rebalanced set of MPI regions.
     * @param regions One inclusive background-element region per MPI rank, ordered by rank.
     *
     * Particles outside the global background mesh are marked for removal.  Particles inside
     * the mesh may be transferred to any rank; the resulting peer list is stored in
     * `par_comm_.comm_ranks` for the subsequent point-variable communication.
     */
    void DetermineDLBParticleRank(const std::vector<mpm_dlb::Region> &regions);

    /**
     * @brief Redistribute particles to a newly rebalanced MPI partition.
     * @param regions One inclusive background-element region per MPI rank, ordered by rank.
     *
     * The function prepares the DLB particle communication plan, updates the local particle count,
     * and transfers the physics-specific particle state.
     */
    void RebalanceDLBParticles(const std::vector<mpm_dlb::Region> &regions);

    /**
     * @brief Apply one complete dynamic load-balance update to this MPM object.
     *
     * Samples particles, repartitions MPI regions, redistributes particle data,
     * and rebuilds mesh, control-point, nodal-volume, and BC data.  Time-step
     * scheduling and solver-system rebuilds remain the responsibility of the caller.
     */
    virtual void ApplyDLB();

    /**
     * @brief Migrate all physics-specific particle state using the current particle communication plan.
     *
     * `par_comm_` and the new local particle count must be prepared before this function is called.
     */
    virtual void MigrateParticleData() {}

    virtual void MoveParticle() {};
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
     * @brief Apply prescribed Dirichlet values to a nodal displacement increment.
     * @param ndispl Displacement increment modified in-place at constrained DOFs.
     */
    virtual void ApplyDisplacementBC(std::vector<double> &ndispl) {
        this->ubc.BCSetDispl(nuc, ndispl);
        this->vbc.BCSetDispl(nvc, ndispl);
        this->wbc.BCSetDispl(nwc, ndispl);

        return;
    }

    /**
     * @brief Divide a nodal vector by a nodal weight with a small-weight cutoff.
     * @param result    Output vector (zero-initialized and overwritten in-place).
     * @param numerator Nodal vector to divide.
     * @param weight    Nodal weight vector used as denominator (e.g. mass or volume).
     * @param offsets   Component offsets into the nodal vectors.
     */
    void CutOffSmallNodalVar(std::vector<double> &result, const std::vector<double> &numerator, const std::vector<double> &weight,
                             const std::vector<int> &offsets) {

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
    void CutOffSmallNodalVar(std::vector<double> &var, const std::vector<double> &weight, const std::vector<int> &offsets) {

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

    /**
     * @brief Rebuild local boundary-condition IDs after a DLB region change.
     *
     * The default is empty because FEM objects do not participate in particle
     * DLB. MPM solid and fluid objects rebuild their physics-specific BCs.
     */
    virtual void RebuildBC() {}
    // ---------------

    // --- Linear-system callbacks used by CrsMat ---
    virtual void BuildPetscBCList(CrsMat &mat) {};

    /**
     * @brief Assemble field CSR values into the distributed matrix.
     * @param mat Target linear system with local element contributions.
     */
    virtual void AssemblePetscMat(CrsMat &mat);

    /**
     * @brief Configure the default BoomerAMG preconditioner.
     * @param mat Linear system for the current physical field.
     * @param pc PETSc preconditioner associated with the system solver.
     */
    virtual void ConfigurePreconditioner(CrsMat &mat, PC pc) {
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");

        return;
    }

    virtual void BCResidualSet(std::vector<double> &rr) {};
    // ------------------------------------------

    // --- Particle part ---
    int num = 0;                          // --- Particle number ---
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

    virtual void InputPointData(std::ifstream &infile) {};

    virtual void InitializePointData() {};

    virtual void RestartInput() {};

    virtual void OutputPointDataVTKHDF(int iview) {};

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

        int sign = 0;
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
     * @brief Generate material points at the sub-cell centers of each background element (see SubcellOffsets).
     */
    void BuildSubcellPoints();

    /**
     * @brief Compute mass-weighted nodal unit normals.
     */
    void CalNodalUnitNormal();

    /**
     * @brief Interpolate nodal unit normals and volume fractions to material points.
     * @param particle_normal Output interface-normal directions; normalization is left to the caller.
     * @param surface_weight Output volume fraction interpolated at each particle position.
     */
    void CalPointUnitNormal(std::vector<std::array<double, 3>> &particle_normal, std::vector<double> &surface_weight);

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
     * @brief Remove Cartesian boundary-normal shifts within one cell of the domain faces.
     * @param particle_coord Current material-point position used to identify Cartesian boundary cells.
     * @param disp_corr Particle-shifting displacement with boundary-normal components removed in place.
     */
    void ConstrainPSTBoundaryDisplacement(const std::array<double, 3> &particle_coord, std::array<double, 3> &disp_corr) const;


    /**
     * @brief Compute a particle-shifting correction to regularize particle distribution.
     *
     * Based on the particle shifting scheme in Chandra et al.,
     * "Stabilized mixed material point method for incompressible fluid flow analysis,"
     * Computer Methods in Applied Mechanics and Engineering, 419, 116644, 2024.
     *
     * @return Correction displacement for each particle, with Cartesian boundary-normal components removed.
     */
    std::vector<std::array<double, 3>> DeltaCorrectionPST() const;

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
    std::vector<std::array<double, 3>> PairwiseRepulsivePST(const std::vector<std::array<double, 3>> &displ);
    // ----------------------------------

    // ----- Control point variable -----
    std::vector<double> nmass, nvof, nmome, nvel,               //
        ndispl, npres, npres_old, nphi, nnormal, nrigid_normal, //
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
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_var, std::vector<T> &node_var) const noexcept;

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
                                   const std::vector<std::array<T, 3>> &point_var, std::vector<T> &node_var) const noexcept {
    node_var[nid + nuc] += sfi * point_mass[pid] * point_var[pid][0];
    node_var[nid + nvc] += sfi * point_mass[pid] * point_var[pid][1];
    node_var[nid + nwc] += sfi * point_mass[pid] * point_var[pid][2];

    return;
}
