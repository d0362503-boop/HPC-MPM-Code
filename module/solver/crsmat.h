#pragma once

#include <petsc.h>

#include <vector>

#include "../bc.h"
#include "../mesh.h"
#include "../mpi_data.h"

class MaterialPoint; // forward declaration

class CrsMat {
  public:
    // Pointer to the owning physics object (fluid or solid). Used for BC
    // setup and residual callbacks. Inactive-node detection is now performed
    // by BuildActiveRowMask() based on assembled matrix content.
    MaterialPoint *owner_ = nullptr;

    // true  → FEM path (all nodes active, no identity fill needed).
    // false → MPM path (inactive nodes are skipped and later identity-filled).
    bool FEM_flag = true;

    int ndof;        // degrees of freedom per node (fluid=4, solid=3)
    int nmata, nmat; // total scalar entries in amat, and CSR entry count
    int local_node, ghost_node;

    // CSR structure (fixed after BuildCrsMat).
    //   matrow[n]      → offset of row n in matcolid / amat.
    //   matcolid[j]    → natural column id of CSR entry j.
    //   block_id[b]    → offset of scalar block b inside each CSR entry.
    std::vector<int> matcolid, matrow, block_id;

    // Dense values stored in block-major order:
    //   amat[j + block_id[row_var*ndof + col_var]]
    // j       = CSR entry index (0 .. nmat-1)
    // block_id= offset of the ndof*ndof dense block inside that entry.
    std::vector<double> amat, adiag, b_rhs, x_lhs;

    void BuildCrsMat(int num_block);

    inline int FindIndex(int nid, int njd, int &ncol) const {
        int cole = this->matrow[nid + 1] - 1;
        int cols = this->matrow[nid] + ncol;
        for (int i = cols; i <= cole; i++) {
            if (this->matcolid[i] == njd) break;
            ++ncol;
        }

        return this->matrow[nid] + ncol;
    }

    void ExtractDiagonal(int ndof);

    void ComputeDiagonalInverseSqrt(int ndof);

    void ApplyDiagonalScaling(int ndof);

    void BuildDiagonalPreconditioner(int ndof);

    void ScaleSolution(std::vector<double> &rr);

    void RestorSolution();

    std::vector<double> MatVecMult(const std::vector<double> &xx);

    void ComputePetscResidualStats(double &res_norm, double &active_dof);

    double ComputeRefResidual();

    double ComputeAbsResidual();

    void BuildActiveRowMask();

    bool CheckNRConvergence(int NR_it, int NR_it_max, int solver_it, double &r0r);

    // --- PETSc distributed objects ---
    ISLocalToGlobalMapping lgmap = nullptr;
    std::vector<PetscInt> l2g_var_map; // scalar variable-major L2G map
    Mat petsc_mat = nullptr;
    Vec petsc_b = nullptr;
    Vec petsc_x = nullptr;
    KSP ksp = nullptr;
    std::vector<int> interior_list;      // locally-owned natural node ids
    std::vector<PetscInt> petsc_bc_gids; // cached global BC row/column IDs

    VecScatter scatter_to_all = nullptr; // local scatter: petsc_x → seq_x
    Vec seq_x = nullptr;                 // sequential copy of solution on this rank

    // AMG lifecycle state (local to this solver, independent of other solvers)
    int amg_rebuild_freq;   // rebuild preconditioner every N steps
    int prev_ksp_its_ = -1; // KSP iterations of the previous solve
    bool force_rebuild_next_ = false;
    bool use_fieldsplit = false;
    bool pressure_pc_use_amg = false;
    // Reusable buffers for RHS / initial-guess insertion (size = local_node*ndof)
    std::vector<PetscInt> petsc_indices_buf;
    std::vector<PetscScalar> petsc_values_buf;

    // Reusable column buffer for matrix assembly (size = max_row_nnz * ndof)
    std::vector<PetscInt> petsc_cols_buf;

    bool use_petsc = true;

    // Block-level L2G: one global node ID per natural node.
    std::vector<PetscInt> l2g_block_map;

    // Reusable block-oriented buffers for AssemblePetscMatBlocked.
    // Sizes are fixed in BuildPetscMat() to avoid heap allocation in the hot loop.
    std::vector<PetscInt> petsc_block_cols_buf;
    std::vector<PetscScalar> petsc_block_vals_buf;
    std::vector<char> active_row_mask;

    void BuildLGMAP(int ndof);

    void BuildPetscMat(int ndof);

    void BuildKSPSolver();
    void ConfigureVelocityPressureFieldSplit(PC pc);

    void InitPetscSolver(int ndof);

    void AssemblePetscMat(int ndof);
    void AssemblePetscMatBlocked(int ndof);

    void UpdatePetscRhs(int ndof);

    int SolveWithPetsc(int ndof, int nr_it = -1);

    void BuildPetscBCList();

    void AddBCComponent(const BoundaryCondition &bc, int offset) {
        if (bc.ibc == 0) return;
        for (int i = 0; i < bc.ibc; ++i) {
            int nid = bc.nbc[i];
            this->petsc_bc_gids.push_back(this->l2g_var_map[nid + offset]);
        }

        return;
    }
};
