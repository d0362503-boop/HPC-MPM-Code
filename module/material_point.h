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
    // ---------------

    // --- For Implicit MPM ---
    // --- Nonlinear time integration part ---
    // --- Generalized-α part ---
    int ode_order;
    double spec_rad = 1.0e0, alpha_f = 1.0e0, alpha_m = 1.0e0;

    void GeneralizedAlphaParaSet();

    std::vector<double> GeneralizedAlphaNodeAccelUpdate();

    // --- Newmark-β part ---
    double gamma_nb, beta_nb;
    std::vector<double> nb_para;

    void NewmarkBetaParaSet();

    void PredictNewmarkBetaVelAndAccel(std::vector<double> &nvel_k, std::vector<double> &naccel_k);
    // ------------------------------

    // --- Newton-Raphson ---
    bool NR_flag = true;

    // --- Need to define in each src folder ---
    virtual void BCNRSet() {};

    virtual void BuildPetscBCList(CrsMat &mat) {};

    virtual void BCResidualSet(std::vector<double> &rr) {};
    // ------------------------------------------

    virtual double ComputeNRLumpedMat(double pid, double sfi) {
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

    int SolveSystem(CrsMat &mat, int NR_it = -1) {
        int iter;
        if (mat.use_petsc) {
            mat.AssemblePetscMat(mat.ndof);
            iter = mat.SolveWithPetsc(mat.ndof, NR_it);
        } else {
            iter = GPBiCGSafe(mat);
        }
        return iter;
    }

    void CommitNodalKinematics(const std::vector<double> &nvel_k, const std::vector<double> &naccel_k);

    void CommitParticleKinematics(const std::vector<std::array<double, 3>> &accel_old,
                                  const std::vector<std::array<double, 3>> &disp,
                                  const std::vector<std::array<double, 3>> &delta_corr = {});

    void ImplicitDsfCorr(const std::vector<int> &nc, int nenode, std::vector<std::array<double, 3>> &dsf);
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
    std::vector<std::array<std::array<double, 3>, 3>> def_grad, def_grad_bar;

    virtual void InputPointData(std::ifstream &inflie) {};

    virtual void InitializePointData() {};

    virtual void RestartInput() {};

    virtual void OutputPointDataVTKHDF(int iview, int istep) {};

    virtual void RestartOutput() {};

    virtual void InflowParticles() {};

    virtual void GenerateInflowParticles(int dir, MaterialPoint &ifp, const BoundaryCondition &ifbc) {};

    void MeshPointLinklist();

    void BuildGaussianPoint();

    void CalPointUnitNormal();

    virtual std::array<std::array<double, 3>, 3> ComputeDeltaDefGrad(const std::vector<int> &nc, int nenode,
                                                                     double af_coeff,
                                                                     const std::vector<std::array<double, 3>> &dsf);

    std::vector<std::array<double, 3>> DeltaCorrectionParticleShifting();
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

    void PICFamilyAccelG2P(int pid, int ni, int nid, double sfi, //
                           const std::vector<std::array<double, 3>> &dsf);

    virtual void Node2Particle() {};

    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_var, std::vector<T> &node_var);

    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass, const std::vector<T> &point_var,
                        std::vector<T> &node_var);

    template <typename T>
    void StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                        const std::vector<std::array<T, 3>> &point_var, std::vector<T> &node_var);
    // -----------------------------------
    virtual ~MaterialPoint() = default;
};

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_var,
                                   std::vector<T> &node_var) {
    node_var[nid] += sfi * point_var[pid];

    return;
}

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                                   const std::vector<T> &point_var, std::vector<T> &node_var) {
    node_var[nid] += sfi * point_mass[pid] * point_var[pid];

    return;
}

template <typename T>
void MaterialPoint::StandardVarP2G(int pid, int nid, double sfi, const std::vector<T> &point_mass,
                                   const std::vector<std::array<T, 3>> &point_var, std::vector<T> &node_var) {
    node_var[nid + nuc] += sfi * point_mass[pid] * point_var[pid][0];
    node_var[nid + nvc] += sfi * point_mass[pid] * point_var[pid][1];
    node_var[nid + nwc] += sfi * point_mass[pid] * point_var[pid][2];

    return;
}
