#pragma once

#include <petsc.h>

#include <vector>

#include "module/bc.h"
#include "module/mesh.h"
#include "module/mpi_data.h"

class MaterialPoint; // forward declaration

class CrsMat {
  public:
    // Pointer to the owning physics object. Used for BC setup, residual
    // callbacks, and MPM mass-based inactive-node detection.
    MaterialPoint *owner_ = nullptr;

    // true  → FEM path (all nodes active, no identity fill needed).
    // false → MPM path (inactive nodes are skipped and later identity-filled).
    bool FEM_flag = true;

    int ndof;                              // degrees of freedom per node (fluid=4, solid=3)
    int num_block = 0;                     // stored scalar block count
    std::vector<int> block_row, block_col; // component indices per block
    int nmata, nmat;                       // total scalar entries in amat, and CSR entry count
    int local_node, ghost_node;

    // Natural-order ownership metadata and PETSc-local numbering.
    std::vector<int> owned_natural_ids;
    std::vector<int> ghost_natural_ids;
    std::vector<char> natural_is_owned;
    std::vector<int> natural_to_owned_pos;
    std::vector<int> natural_to_ghost_pos;
    std::vector<int> natural_to_petsc_local;
    std::vector<int> petsc_local_to_natural;

    // Natural-order global IDs: index with natural local node id.
    std::vector<PetscInt> natural_block_gids;
    std::vector<PetscInt> natural_var_gids;

    // PETSc-local global IDs: index with petsc-local position (owned first, ghosts after).
    std::vector<PetscInt> petsc_local_block_gids;
    std::vector<PetscInt> petsc_local_var_gids;

    ISLocalToGlobalMapping petsc_var_lgmap_ = nullptr;
    ISLocalToGlobalMapping petsc_block_lgmap_ = nullptr;

    // CSR structure (fixed after BuildCrsMat).
    //   matrow[n]      → offset of row n in matcolid / amat.
    //   matcolid[j]    → natural column id of CSR entry j.
    //   block_id[b]    → start of scalar block b in amat.
    std::vector<int> matcolid, matrow, block_id;

    // Stored blocks share one nodal CSR graph.
    // amat[j + block_id[b]] couples block_row[b] to block_col[b].
    std::vector<double> amat, adiag, b_rhs, x_lhs;

    /**
     * @brief Build the local CSR sparsity pattern from mesh connectivity.
     * @param num_block Number of scalar blocks per CSR entry (typically `ndof * ndof`).
     */
    void BuildCrsMat(int num_block);

    /**
     * @brief Release PETSc objects associated with the current matrix layout.
     *
     * This is called before rebuilding a matrix after DLB changes the local
     * mesh and ownership layout.
     */
    void ResetPetscSolver();

    /**
     * @brief Locate a control-point pair in the local CSR row.
     * @param nid Natural local control point receiving the equation contribution.
     * @param njd Natural local control point supplying the coupled unknown.
     * @param ncol In/out search offset within the row for successive ordered lookups.
     * @return CSR entry index shared by all stored component blocks.
     * @pre The requested column exists at or after the supplied search offset.
     */
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
     */
    void ApplyDiagonalScaling();

    /**
     * @brief Build a diagonal preconditioner from `adiag`.
     * Requires field diagonals; zero-diagonal multiplier rows need a separate preconditioner.
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
     * @brief Compute global L2 and RMS norms of the assembled physical RHS.
     * @param ref_tol Output overlap-weighted L2 norm of the current equilibrium residual.
     * @param abs_tol Output residual RMS over counted DOFs; not a stopping tolerance.
     * @note Zeros constrained entries of b_rhs in place. MPM counts all components
     * at nodes with nmass > mtol; FEM has no mass filter. Constrained DOFs are
     * excluded from both norms and the count through BCResidualSet(). Both
     * outputs are zero when no free active DOFs are counted.
     */
    void ComputeResidualNormSq(double &ref_tol, double &abs_tol);

    /** @brief Mark low-mass MPM rows inactive. */
    void BuildActiveRowMask();

    /**
     * @brief Check NR convergence using the freshly assembled physical RHS.
     * @param NR_it     Current Newton–Raphson iteration index.
     * @param NR_it_max Maximum allowed Newton–Raphson iterations.
     * @param solver_it Iteration count of the preceding linear solve, or zero initially.
     * @param r0r Initial RHS L2 norm, saved on iteration zero and retained as reference.
     * @return True when the current stopping thresholds pass, false to continue.
     * @note Call after assembling at the current iterate and before solving its
     * correction. This function does not reassemble the residual or compute b-A*x.
     * @warning Failure at the maximum iteration calls MPI_Abort. Absolute thresholds
     * depend on residual units and scale; the initial-norm criterion can accept iteration zero.
     */
    bool CheckNRConvergence(int NR_it, int NR_it_max, int solver_it, double &r0r);

    // --- PETSc distributed objects ---
    Mat petsc_mat = nullptr;
    Vec petsc_b = nullptr;
    Vec petsc_x = nullptr;
    KSP ksp = nullptr;
    std::vector<PetscInt> petsc_bc_gids; // cached global BC row/column IDs

    VecScatter scatter_to_all = nullptr; // local scatter: petsc_x → seq_x
    Vec seq_x = nullptr;                 // sequential copy of solution on this rank

    // AMG lifecycle state (local to this solver, independent of other solvers)
    int amg_rebuild_freq = 1; // rebuild preconditioner every N steps
    int prev_ksp_its_ = -1;   // KSP iterations of the previous solve
    bool force_rebuild_next_ = false;
    // Reusable buffers for RHS / initial-guess insertion (size = local_node*ndof)
    std::vector<PetscInt> petsc_indices_buf;
    std::vector<PetscScalar> petsc_values_buf;

    // Reusable column buffer for matrix assembly (size = max_row_nnz * ndof)
    std::vector<PetscInt> petsc_cols_buf;

    bool use_petsc = true;
    bool use_schur_fieldsplit = false;

    std::vector<char> active_row_mask;

    // Reusable block-oriented buffers for AssemblePetscMat.
    // Sizes are fixed in BuildPetscMat() to avoid heap allocation in the hot loop.
    std::vector<PetscInt> petsc_block_cols_buf;
    std::vector<PetscScalar> petsc_block_vals_buf;

    /**
     * @brief Build the PETSc local-to-global mapping for distributed vectors/matrices.
     * @param ndof Degrees of freedom per node.
     */
    void BuildLGMAP(int ndof);

    /**
     * @brief Build natural-order global ID maps from the local node-to-global mapping.
     * @param ndof Degrees of freedom per node.
     * @param node_l2g Natural local node id to global node id mapping.
     */
    void BuildNaturalGlobalMaps(int ndof, const std::vector<int> &node_l2g);

    /**
     * @brief Classify natural local nodes into owned and ghost sets using the
     *        aelemmin tie-break rule.
     * @param node_l2g Natural local node id to global node id mapping.
     */
    void BuildOwnershipMetadata(const std::vector<int> &node_l2g);

    /**
     * @brief Build PETSc-local numbering (owned first, ghosts after) and the
     *        corresponding global ID arrays.
     * @param ndof Degrees of freedom per node.
     */
    void BuildPetscLocalMaps(int ndof);

    /**
     * @brief Verify ownership metadata invariants before assembly.
     */
    void CheckOwnershipMetadata() const;

    /**
     * @brief Convert a mesh control point to its PETSc-local block index.
     * @param natural_id Natural local control-point index, including ghosts.
     * @return Owned-first block position in the local-to-global mapping.
     */
    inline PetscInt NaturalNodeToPetscLocalBlock(int natural_id) const {
        return static_cast<PetscInt>(this->natural_to_petsc_local[natural_id]);
    }

    /**
     * @brief Convert a control-point component to its PETSc-local scalar index.
     * @param natural_id Natural local control-point index, including ghosts.
     * @param var Physical field component in the system's component ordering.
     * @return Component-major index into the scalar local-to-global mapping.
     */
    inline PetscInt NaturalNodeVarToPetscLocalScalar(int natural_id, int var) const {
        return static_cast<PetscInt>(this->natural_to_petsc_local[natural_id] + var * nodec);
    }

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
     * @brief Delegate preconditioner configuration to the owning physical field.
     * Systems without an owner use BoomerAMG.
     * @param pc PETSc preconditioner associated with `ksp`.
     */
    void ConfigurePreconditioner(PC pc);

    /**
     * @brief Initialize all PETSc objects (matrix, vectors, KSP) for this system.
     * @param ndof Degrees of freedom per node.
     */
    void InitPetscSolver(int ndof);

    /**
     * @brief Assemble the local CSR matrix into the PETSc distributed matrix using block-level mappings.
     * @param ndof Degrees of freedom per node.
     */
    void AssemblePetscMat(int ndof);

    /**
     * @brief Copy the local RHS vector `b_rhs` into the PETSc RHS vector.
     * @param ndof Degrees of freedom per node.
     */
    void UpdatePetscRhs(int ndof);

    /**
     * @brief Solve the linear system with PETSc and scatter the solution back to `x_lhs`.
     * @param ndof Degrees of freedom per node.
     * @param NR_it Current Newton–Raphson iteration (used for AMG rebuild scheduling).
     * @return Number of KSP iterations, including when KSP reports divergence.
     * @warning A negative KSP convergence reason is logged on rank zero; the
     * solution is still scattered and no failure sentinel is returned.
     */
    int SolveWithPetsc(int ndof, int NR_it = -1);

    /**
     * @brief Solve this assembled linear system with PETSc or native GPBiCGAR.
     * @param NR_it Current Newton–Raphson iteration index.
     * @return Number of linear solver iterations used.
     *
     * `use_petsc` selects the PETSc or native path for the current solve.
     */
    int SolveSystem(int NR_it = -1);

    /**
     * @brief Collect global IDs of Dirichlet boundary rows/columns into `petsc_bc_gids`.
     *
     * Delegates to the owning physics object (`owner_->BuildPetscBCList`).
     */
    void BuildPetscBCList();

    /**
     * @brief Append the global scalar IDs of a boundary-condition component to `petsc_bc_gids`.
     * @param bc     Boundary condition descriptor containing natural node ids.
     * @param offset Scalar variable offset (e.g. `nuc`, `nvc`, `nwc`, `npc`).
     */
    void AddBCComponent(const BoundaryCondition &bc, int offset) {
        if (bc.ibc == 0) return;
        for (int i = 0; i < bc.ibc; ++i) {
            const int nid = bc.nbc[i];
            this->petsc_bc_gids.push_back(this->natural_var_gids[nid + offset]);
        }

        return;
    }
};
