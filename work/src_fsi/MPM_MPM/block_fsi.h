#pragma once

#include <array>
#include <vector>

#include "module/bc.h"
#include "module/fluid/MPM/stabilized_mpm.h"
#include "module/solid/implicit/implicit_mpm_solid.h"

class MPMMPMBlockFSI;

class FSIFluid : public stabilizedmpm::StabilizedMPM {
  public:
    MPMMPMBlockFSI &fsi_; // Reference back to the FSI coordinator

    FSIFluid(MPMMPMBlockFSI &fsi) : fsi_(fsi) {}

    // ---- Boundary Condition Overrides ----
    void BCNRSet() override;

    // --- For PETSc ---
    void BuildPetscBCList(CrsMat &mat) override;

    // --- For self-defined ---
    void BCResidualSet(std::vector<double> &rr) override;
};

class FSISolid : public implicitmpm::ImplicitSolidMPM {
  public:
    MPMMPMBlockFSI &fsi_; // Reference back to the FSI coordinator

    FSISolid(MPMMPMBlockFSI &fsi) : fsi_(fsi) {}

    void AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel) override;
};

class MPMMPMBlockFSI {
  public:
    BoundaryCondition fsi_intf;

    // --- Sub-solvers (constructed with reference to this coordinator) ---
    FSIFluid fluid_;
    FSISolid solid_;

    // Block-iteration control
    const int max_block_iter = 100;
    static constexpr double phi_cut = 0.15e0;
    double relax_omega = 0.1e0;
    const double tol_ref = 1.0e-5;
    const double tol_abs = 1.0e-8;

    // -----------------------------------------------------------------
    // Constructor: pass *this to sub-solvers
    // -----------------------------------------------------------------
    MPMMPMBlockFSI() : fluid_(*this), solid_(*this) {}

    std::vector<double> nfsi_force;

    /**
     * @brief Read and initialize the coupled FSI case input data.
     *
     * Loads shared global parameters, then dispatches fluid and solid
     * boundary-condition/particle input to the corresponding sub-solvers.
     */
    void DataInput();

    void DetectFSIInterface();

    void PredictFSIInterfaceBC();

    void SolveFSISystem();

    void CalFSIResidual(double &rtr_ref, double &rtr_dof);

    void CalFSIForce();
};
