#pragma once

#include <vector>

#include "../module/bc.h"
#include "../module/fluid/FEM/stabilized_fem.h"
#include "../module/solid/implicit/implicit_mpm_solid.h"

class BlockFSI;

class FSIFluid : public stabilizedfem::StabilizedFEM {
  public:
    BlockFSI &fsi_; // Reference back to the FSI coordinator

    FSIFluid(BlockFSI &fsi) : fsi_(fsi) {}

    // ---- Boundary Condition Overrides ----
    // --- For PETSc ---
    void BuildPetscBCList(CrsMat &mat) override;

    // --- For self-defined ---
    void BCResidualSet(std::vector<double> &rr) override;

    void BCSet() override;
};

class FSISolid : public implicitmpm::SolidMaterialPoint {
  public:
    BlockFSI &fsi_; // Reference back to the FSI coordinator

    FSISolid(BlockFSI &fsi) : fsi_(fsi) {}

    double ComputeNRLumpedMat(double pid, double sfi) override;

    std::array<double, 3> ComputeExternalForce(int pid, double sfi) override;

    void AddInertialForceToRHS(CrsMat &mat, const std::vector<double> &naccel) override;
};

class BlockFSI {
  public:
    BoundaryCondition fsi_intf;

    // --- Sub-solvers (constructed with reference to this coordinator) ---
    FSIFluid fluid_;
    FSISolid solid_;

    // Block-iteration control
    const int max_block_iter = 100;
    double relax_omega = 0.1e0;
    const double tol_ref = 1.0e-3;
    const double tol_abs = 1.0e-6;

    // -----------------------------------------------------------------
    // Constructor: pass *this to sub-solvers
    // -----------------------------------------------------------------
    BlockFSI() : fluid_(*this), solid_(*this) {}

    std::vector<double> added_mass;
    std::vector<double> nfsi_force;

    void DataInput();

    void DetectFSIInterface();

    void UpdateFSIInterfaceBC();

    void SolveFSISystem();

    void CalFSIResidual(double &rtr_ref, double &rtr_dof, const std::vector<double> &nvel_k);

    void CalDragLiftCoeff();

    void CalFSIForce();
};