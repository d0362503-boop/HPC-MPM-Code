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

    /**
     * @brief Build the local CSR sparsity pattern from mesh connectivity.
     * @param num_block Number of scalar blocks per CSR entry (typically `ndof * ndof`).
     */
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

    /**
     * @brief Extract the block-diagonal entries of `amat` into `adiag`.
     * @param ndof Degrees of freedom per node.
     */
    void ExtractDiagonal(int ndof);

    /**
     * @brief Compute `adiag[i] = 1 / sqrt(adiag[i])` for diagonal scaling.
     * @param ndof Degrees of freedom per node.
     */
    void ComputeDiagonalInverseSqrt(int ndof);

    /**
     * @brief Apply left and right diagonal scaling to the CSR matrix.
     * @param ndof Degrees of freedom per node.
     */
    void ApplyDiagonalScaling(int ndof);

    /**
     * @brief Build a diagonal preconditioner from `adiag`.
     * @param ndof Degrees of freedom per node.
     */
    void BuildDiagonalPreconditioner(int ndof);

    /**
     * @brief Scale a vector by the diagonal scaling stored in `adiag`.
     * @param rr Vector to be scaled in-place.
     */
    void ScaleSolution(std::vector<double> &rr);

    /**
     * @brief Reverse the diagonal scaling applied to `x_lhs`.
     */
    void RestorSolution();

    /**
     * @brief Compute the local matrix-vector product `y = A * x` in CSR storage.
     * @param xx Input vector.
     * @return Result vector.
     */
    std::vector<double> MatVecMult(const std::vector<double> &xx);

    /**
     * @brief Compute the PETSc residual norm and active-DOF count for convergence monitoring.
     * @param res_norm  Output L2 norm of the residual.
     * @param active_dof Output number of active degrees of freedom.
     */
    void ComputePetscResidualStats(double &res_norm, double &active_dof);

    /**
     * @brief Compute a reference residual used for Newton-Raphson convergence checks.
     * @return Reference residual value.
     */
    double ComputeRefResidual();

    /**
     * @brief Compute the absolute residual used for Newton-Raphson convergence checks.
     * @return Absolute residual value.
     */
    double ComputeAbsResidual();

    /**
     * @brief Mark rows with negligible assembled entries as inactive for MPM solves.
     */
    void BuildActiveRowMask();

    /**
     * @brief Check whether the Newton–Raphson iteration has converged.
     * @param NR_it     Current Newton–Raphson iteration index.
     * @param NR_it_max Maximum allowed Newton–Raphson iterations.
     * @param solver_it Linear solver iteration count for the current NR step.
     * @param r0r       Squared residual ratio (updated in-place).
     * @return True if converged, false otherwise.
     */
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

    /**
     * @brief Build the PETSc local-to-global mapping for distributed vectors/matrices.
     * @param ndof Degrees of freedom per node.
     */
    void BuildLGMAP(int ndof);

    /**
     * @brief Create the distributed PETSc matrix from the local CSR pattern.
     * @param ndof Degrees of freedom per node.
     */
    void BuildPetscMat(int ndof);

    /**
     * @brief Create and configure the PETSc KSP solver object.
     */
    void BuildKSPSolver();

    /**
     * @brief Configure a velocity-pressure field-split preconditioner for the fluid solver.
     * @param pc PETSc preconditioner context.
     */
    void ConfigureVelocityPressureFieldSplit(PC pc);

    /**
     * @brief Initialize all PETSc objects (matrix, vectors, KSP) for this system.
     * @param ndof Degrees of freedom per node.
     */
    void InitPetscSolver(int ndof);

    /**
     * @brief Assemble the local CSR matrix `amat` into the PETSc distributed matrix.
     * @param ndof Degrees of freedom per node.
     */
    void AssemblePetscMat(int ndof);

    /**
     * @brief Assemble the local CSR matrix into the PETSc matrix using block-level mappings.
     * @param ndof Degrees of freedom per node.
     */
    void AssemblePetscMatBlocked(int ndof);

    /**
     * @brief Copy the local RHS vector `b_rhs` into the PETSc RHS vector.
     * @param ndof Degrees of freedom per node.
     */
    void UpdatePetscRhs(int ndof);

    /**
     * @brief Solve the linear system with PETSc and scatter the solution back to `x_lhs`.
     * @param ndof Degrees of freedom per node.
     * @param nr_it Current Newton–Raphson iteration (used for AMG rebuild scheduling).
     * @return Number of KSP iterations, or -1 on failure.
     */
    int SolveWithPetsc(int ndof, int nr_it = -1);

    /**
     * @brief Collect global IDs of Dirichlet boundary rows/columns into `petsc_bc_gids`.
     *
     * Delegates to the owning physics object (`owner_->BuildPetscBCList`).
     */
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
