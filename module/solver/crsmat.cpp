#include "crsmat.h"

#include <petsc.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../bc.h"
#include "../dataset.h"
#include "../material_point.h"
#include "../mesh.h"
#include "../mpi_data.h"

/**
 * @brief Build the static CSR graph and dense block storage.
 * @param num_block Number of dense blocks per CSR entry.
 */
void CrsMat::BuildCrsMat(int num_block) {
    this->nmat = 0;

    int nxp = xynodec[0];
    int nyp = xynodec[1];
    int nzp = xynodec[2];

    int idimp2 = *max_element(begin(idimc), end(idimc));

    int nmatptmp = nodec * pow((idimp2 + idimp2 + 1), 3);
    std::vector<int> matcolidp(nmatptmp, -1);
    VectorAssign(nodec + 1, this->matrow, -1);

    for (int icnzp = 0; icnzp < nzp; ++icnzp) {
        int isnzp = icnzp - idimc[2];
        if (isnzp < 0) isnzp = 0;
        int ienzp = icnzp + idimc[2];
        if (ienzp >= nzp) ienzp = nzp - 1;

        for (int icnyp = 0; icnyp < nyp; ++icnyp) {
            int isnyp = icnyp - idimc[1];
            if (isnyp < 0) isnyp = 0;
            int ienyp = icnyp + idimc[1];
            if (ienyp >= nyp) ienyp = nyp - 1;

            for (int icnxp = 0; icnxp < nxp; ++icnxp) {
                int isnxp = icnxp - idimc[0];
                if (isnxp < 0) isnxp = 0;
                int ienxp = icnxp + idimc[0];
                if (ienxp >= nxp) ienxp = nxp - 1;

                int idn = icnxp + nxp * icnyp + nxp * nyp * icnzp;
                this->matrow[idn] = this->nmat;

                for (int inz = isnzp; inz <= ienzp; ++inz) {
                    for (int iny = isnyp; iny <= ienyp; ++iny) {
                        for (int inx = isnxp; inx <= ienxp; ++inx) {
                            int idnmat = inx + nxp * iny + nxp * nyp * inz;
                            matcolidp[this->nmat] = idnmat;
                            this->nmat++;
                        }
                    }
                }
            }
        }
    }

    VectorAssign(this->nmat, this->matcolid);
    for (int i = 0; i < this->nmat; i++) { this->matcolid[i] = matcolidp[i]; }

    this->matrow[nodec] = this->nmat;

    this->nmata = this->nmat * num_block;
    VectorAssign(num_block, this->block_id);
    for (int i = 0; i < num_block; i++) { this->block_id[i] = this->nmat * i; }
    VectorAssign(this->nmata, this->amat);

    if (this->use_petsc) { this->InitPetscSolver(this->ndof); }

    return;
}

/**
 * @brief Extract the block diagonal into adiag and synchronize overlap values.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::ExtractDiagonal(int ndof) {
    for (int i = 0; i < nodec; i++) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; j++) {
            if (i == this->matcolid[j]) {
                for (int n = 0; n < ndof; n++) {
                    const int offset = nodec * n;
                    const int bid = this->block_id[n + ndof * n];
                    this->adiag[i + offset] = this->amat[j + bid];
                }
                break;
            }
        }
    }

    std::vector<int> offsets(ndof);
    for (int n = 0; n < ndof; ++n) { offsets[n] = nodec * n; }
    NodeVarComm(this->adiag, offsets);

    return;
}

/**
 * @brief Mark rows that carry meaningful physical content in the current matrix.
 */
void CrsMat::BuildActiveRowMask() {
    this->active_row_mask.assign(nodec, 1);
    if (this->FEM_flag) return;
    constexpr double kZeroRowTol = 1.0e-12;

    for (int nid = 0; nid < nodec; ++nid) {
        double row_abs_sum = 0.0e0;
        for (int row_var = 0; row_var < this->ndof; ++row_var) {
            for (int j = this->matrow[nid]; j < this->matrow[nid + 1]; ++j) {
                for (int col_var = 0; col_var < this->ndof; ++col_var) {
                    int bid = this->block_id[row_var * this->ndof + col_var];
                    row_abs_sum += std::abs(this->amat[j + bid]);
                }
            }
        }
        if (row_abs_sum < mtol) this->active_row_mask[nid] = 0;
    }
    return;
}

/**
 * @brief Convert adiag to inverse square-root scaling factors.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::ComputeDiagonalInverseSqrt(int ndof) {
    int var_size = nodec * ndof;
    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->adiag[n] = 1.0e0 / std::sqrt(std::abs(this->adiag[n]));
        } else {
            this->adiag[n] = (std::abs(this->adiag[n]) < mtol) ? 0.0e0 : (1.0e0 / std::sqrt(std::abs(this->adiag[n])));
        }
    }

    return;
}

/**
 * @brief Apply symmetric diagonal scaling to the local matrix blocks.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::ApplyDiagonalScaling(int ndof) {
    for (int i = 0; i < nodec; i++) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; j++) {
            const int col = this->matcolid[j];

            for (int m = 0; m < ndof; m++) {
                const int row_offset = nodec * m;
                const double diag_i = this->adiag[i + row_offset];

                for (int n = 0; n < ndof; n++) {
                    const int col_offset = nodec * n;
                    const int bid = this->block_id[m * ndof + n];

                    this->amat[j + bid] *= diag_i * this->adiag[col + col_offset];
                }
            }
        }
    }

    return;
}

/**
 * @brief Build the native diagonal preconditioner in-place.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::BuildDiagonalPreconditioner(int ndof) {
    this->ExtractDiagonal(ndof);
    this->ComputeDiagonalInverseSqrt(ndof);
    this->ApplyDiagonalScaling(ndof);

    return;
}

/**
 * @brief Scale the iterate and RHS for the native solver path.
 * @param rr Output buffer for the scaled right-hand side.
 */
void CrsMat::ScaleSolution(std::vector<double> &rr) {
    const int var_size = this->x_lhs.size();

    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->x_lhs[n] /= this->adiag[n];
            rr[n] = this->b_rhs[n] * this->adiag[n];
        } else {
            if (this->adiag[n] < mtol) {
                this->x_lhs[n] = 0.0e0;
                rr[n] = 0.0e0;
            } else {
                this->x_lhs[n] /= this->adiag[n];
                rr[n] = this->b_rhs[n] * this->adiag[n];
            }
        }
    }

    return;
}

/**
 * @brief Restore the physical solution after native diagonal scaling.
 */
void CrsMat::RestorSolution() {
    const int var_size = this->x_lhs.size();

    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->x_lhs[n] *= this->adiag[n];
        } else {
            if (this->adiag[n] > mtol) {
                this->x_lhs[n] *= this->adiag[n];
            } else {
                this->x_lhs[n] = 0.0e0;
            }
        }
    }

    return;
}

/**
 * @brief Build PETSc local-to-global maps and owned-node lists.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::BuildLGMAP(int ndof) {
    std::vector<int> node_l2g(nodec);

    for (int kl = 0; kl < xynodec[2]; kl++) {
        for (int jl = 0; jl < xynodec[1]; jl++) {
            for (int il = 0; il < xynodec[0]; il++) {
                int kg = kl + aelemmin[2];
                int jg = jl + aelemmin[1];
                int ig = il + aelemmin[0];
                int nid_l = il + xynodec[0] * jl + xynodec[0] * xynodec[1] * kl;
                int nid_g = ig + xynodecw[0] * jg + xynodecw[0] * xynodecw[1] * kg;
                node_l2g[nid_l] = nid_g;
            }
        }
    }

    // Gather every rank's global element anchor so that each rank can
    // determine which shared nodes it truly owns.
    std::vector<int> all_aelemmin(nprocs * 3);
    std::vector<int> all_aelemmax(nprocs * 3);
    MPI_Allgather(aelemmin.data(), 3, MPI_INT, all_aelemmin.data(), 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(aelemmax.data(), 3, MPI_INT, all_aelemmax.data(), 3, MPI_INT, MPI_COMM_WORLD);

    // Classify local nodes into owned (interior) and ghost (overlap).
    this->interior_list.clear();

    for (int n = 0; n < nodec; ++n) {
        int global_node = node_l2g[n];
        int kg = global_node / (xynodecw[0] * xynodecw[1]);
        int rem = global_node % (xynodecw[0] * xynodecw[1]);
        int jg = rem / xynodecw[0];
        int ig = rem % xynodecw[0];

        bool is_owned = true;
        for (int p = 0; p < nprocs; ++p) {
            if (p == myrank) continue;

            bool in_p = true;
            for (int d = 0; d < 3; ++d) {
                int p_min = all_aelemmin[p * 3 + d];
                int p_max = all_aelemmax[p * 3 + d] + idimc[d];
                int coord = (d == 0) ? ig : (d == 1) ? jg : kg;
                if (coord < p_min || coord > p_max) {
                    in_p = false;
                    break;
                }
            }

            if (in_p) {
                for (int d = 0; d < 3; ++d) {
                    int my_min = aelemmin[d];
                    int p_min = all_aelemmin[p * 3 + d];
                    if (p_min < my_min) {
                        is_owned = false;
                        break;
                    } else if (p_min > my_min) {
                        break;
                    }
                }
                if (!is_owned) break;
            }
        }

        if (is_owned) { this->interior_list.push_back(n); }
    }

    this->local_node = this->interior_list.size();
    this->ghost_node = nodec - this->local_node;

    // Sanity check: the sum of owned nodes across all ranks must equal the
    // global control-point count.  If not, the tie-break logic is inconsistent.
    int global_nodec = xynodecw[0] * xynodecw[1] * xynodecw[2];
    int total_local_node = 0;
    MPI_Allreduce(&this->local_node, &total_local_node, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (total_local_node != global_nodec) { MPI_Abort(MPI_COMM_WORLD, 1); }

    // Build the scalar variable-major local-to-global map.
    // Layout: [node0_var0, node1_var0, ... nodeN_var0, node0_var1, ...]
    VectorAssign(nodec * ndof, this->l2g_var_map);
    VectorAssign(nodec, this->l2g_block_map);
    for (int n = 0; n < nodec; n++) {
        int global_node = node_l2g[n];
        this->l2g_block_map[n] = global_node;
        for (int var = 0; var < ndof; var++) {
            int local_idx = n + var * nodec;
            int global_id = global_node * ndof + var;
            this->l2g_var_map[local_idx] = global_id;
        }
    }

    // Create the PETSc LGMAP object.  It is not bound to the matrix because
    // we insert values directly via global indices using MatSetValues.
    ISLocalToGlobalMappingCreate(MPI_COMM_WORLD, 1, nodec * ndof, this->l2g_var_map.data(), PETSC_COPY_VALUES,
                                 &this->lgmap);

    return;
}

/**
 * @brief Create and preallocate the distributed PETSc matrix and vectors.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::BuildPetscMat(int ndof) {
    PetscInt local_n = static_cast<PetscInt>(this->local_node) * ndof;

    MatCreate(MPI_COMM_WORLD, &this->petsc_mat);
    MatSetSizes(this->petsc_mat, local_n, local_n, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(this->petsc_mat, MATAIJ);
    MatSetBlockSize(this->petsc_mat, ndof);

    // Compute the maximum nonzeros per row across ALL nodec (owned + ghost).
    // Ghost rows are used during assembly (remote contributions), so the
    // reusable column/value buffers must be large enough for them too.
    int max_row_nnz = 0;
    for (int ii = 0; ii < nodec; ++ii) {
        int row_nnz = this->matrow[ii + 1] - this->matrow[ii];
        if (row_nnz > max_row_nnz) max_row_nnz = row_nnz;
    }

    PetscInt prealloc = static_cast<PetscInt>(max_row_nnz) * ndof;
    MatMPIAIJSetPreallocation(this->petsc_mat, prealloc, NULL, prealloc, NULL);
    MatSeqAIJSetPreallocation(this->petsc_mat, prealloc, NULL);

    // Keep the sparsity pattern across solves so that repeated assemblies
    // do not reallocate internal PETSc data structures.
    MatSetOption(this->petsc_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    // Allow dynamic allocation for new nonzeros caused by moving material
    // points (MPM) or changing active sets, preventing PETSc malloc errors.
    MatSetOption(this->petsc_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // Create distributed RHS and solution vectors matching the matrix layout.
    VecCreateMPI(MPI_COMM_WORLD, local_n, PETSC_DETERMINE, &this->petsc_b);
    VecCreateMPI(MPI_COMM_WORLD, local_n, PETSC_DETERMINE, &this->petsc_x);

    VecZeroEntries(this->petsc_b);
    VecZeroEntries(this->petsc_x);

    this->petsc_cols_buf.resize(static_cast<size_t>(max_row_nnz) * ndof);

    // Block-oriented buffers for the optional blocked assembly path.
    this->petsc_block_cols_buf.resize(static_cast<size_t>(max_row_nnz));
    this->petsc_block_vals_buf.resize(static_cast<size_t>(max_row_nnz) * ndof * ndof);

    return;
}

/**
 * @brief Configure the velocity-pressure FieldSplit preconditioner.
 * @param pc PETSc preconditioner object.
 */
void CrsMat::ConfigureVelocityPressureFieldSplit(PC pc) {
    PetscInt fields_u[] = {0, 1, 2};
    PetscInt fields_p[] = {3};

    PCSetType(pc, PCNONE);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);
    PCFieldSplitSetBlockSize(pc, 4);
    PCFieldSplitSetFields(pc, "0", 3, fields_u, fields_u);
    PCFieldSplitSetFields(pc, "1", 1, fields_p, fields_p);
}

/**
 * @brief Build the PETSc KSP solver and default preconditioner stack.
 */
void CrsMat::BuildKSPSolver() {
    KSPCreate(MPI_COMM_WORLD, &this->ksp);

    KSPSetOperators(this->ksp, this->petsc_mat, this->petsc_mat);

    KSPSetType(this->ksp, KSPFGMRES);
    // KSPSetType(this->ksp, KSPBCGS);

    PC pc;
    KSPGetPC(this->ksp, &pc);
    if (this->use_fieldsplit && this->ndof == 4) {
        this->ConfigureVelocityPressureFieldSplit(pc);
        PetscOptionsSetValue(NULL, "-pc_fieldsplit_type", "multiplicative");
        PetscOptionsSetValue(NULL, "-pc_fieldsplit_block_size", "4");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_ksp_type", "preonly");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_type", "hypre");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_type", "boomeramg");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_coarsen_type", "HMIS");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_interp_type", "ext+i");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_agg_nl", "2");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_P_max", "2");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_strong_threshold", "0.8");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_max_levels", "8");
        PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_hypre_boomeramg_nodal_coarsen", "6");
        PetscOptionsSetValue(NULL, "-fieldsplit_1_ksp_type", "preonly");
        if (this->pressure_pc_use_amg) {
            PetscOptionsSetValue(NULL, "-fieldsplit_1_pc_type", "hypre");
            PetscOptionsSetValue(NULL, "-fieldsplit_1_pc_hypre_type", "boomeramg");
        } else {
            PetscOptionsSetValue(NULL, "-fieldsplit_1_pc_type", "jacobi");
        }
    } else {
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "HMIS");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_interp_type", "ext+i");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_agg_nl", "2");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_P_max", "2");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_strong_threshold", "0.8");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_levels", "8");
        PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_nodal_coarsen", "6");
    }

    KSPSetTolerances(this->ksp, 1.0e-8, 1.0e-10, 1.0e6, 1000);

    KSPSetFromOptions(this->ksp);

    return;
}

/**
 * @brief Rebuild and deduplicate the PETSc Dirichlet index list.
 */
void CrsMat::BuildPetscBCList() {
    this->petsc_bc_gids.clear();

    if (this->owner_) { this->owner_->BuildPetscBCList(*this); }

    std::sort(this->petsc_bc_gids.begin(), this->petsc_bc_gids.end());
    this->petsc_bc_gids.erase(std::unique(this->petsc_bc_gids.begin(), this->petsc_bc_gids.end()),
                              this->petsc_bc_gids.end());
}

/**
 * @brief Initialize the PETSc matrix, solver, and local scatter objects.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::InitPetscSolver(int ndof) {
    BuildLGMAP(ndof);

    BuildPetscMat(ndof);

    BuildKSPSolver();

    IS from_is = nullptr;
    IS to_is = nullptr;
    ISCreateGeneral(PETSC_COMM_SELF, nodec * ndof, this->l2g_var_map.data(), PETSC_COPY_VALUES, &from_is);
    VecCreateSeq(PETSC_COMM_SELF, nodec * ndof, &this->seq_x);
    ISCreateStride(PETSC_COMM_SELF, nodec * ndof, 0, 1, &to_is);
    VecScatterCreate(this->petsc_x, from_is, this->seq_x, to_is, &this->scatter_to_all);
    ISDestroy(&from_is);
    ISDestroy(&to_is);

    return;
}

/**
 * @brief Assemble the PETSc matrix through the blocked assembly path.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::AssemblePetscMat(int ndof) { this->AssemblePetscMatBlocked(ndof); }

/**
 * @brief Assemble block rows and patch inactive owned rows with identity blocks.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::AssemblePetscMatBlocked(int ndof) {
    MatZeroEntries(this->petsc_mat);
    std::vector<PetscScalar> identity_block(static_cast<size_t>(ndof) * ndof, 0.0);
    for (int var = 0; var < ndof; ++var) { identity_block[var * ndof + var] = 1.0; }
    this->BuildActiveRowMask();

    for (int i = 0; i < nodec; ++i) {
        bool is_inactive = (!this->FEM_flag && this->active_row_mask[i] == 0);
        if (is_inactive) continue;

        int natural_row = i;
        int row_start = this->matrow[natural_row];
        int row_end = this->matrow[natural_row + 1];
        int ncols = row_end - row_start;

        PetscInt block_row = this->l2g_block_map[natural_row];

        size_t block_col_idx = 0;
        for (int j = row_start; j < row_end; ++j) {
            int natural_col = this->matcolid[j];
            this->petsc_block_cols_buf[block_col_idx++] = this->l2g_block_map[natural_col];
        }

        size_t val_idx = 0;
        for (int row_var = 0; row_var < ndof; ++row_var) {
            for (int j = row_start; j < row_end; ++j) {
                for (int col_var = 0; col_var < ndof; ++col_var) {
                    this->petsc_block_vals_buf[val_idx] = this->amat[j + this->block_id[row_var * ndof + col_var]];
                    ++val_idx;
                }
            }
        }

        MatSetValuesBlocked(this->petsc_mat, 1, &block_row, ncols, this->petsc_block_cols_buf.data(),
                            this->petsc_block_vals_buf.data(), ADD_VALUES);
    }

    if (!this->FEM_flag) {
        for (int ii = 0; ii < this->local_node; ++ii) {
            int natural_id = this->interior_list[ii];
            if (this->active_row_mask[natural_id] != 0) continue;

            PetscInt block_row = this->l2g_block_map[natural_id];
            MatSetValuesBlocked(this->petsc_mat, 1, &block_row, 1, &block_row, identity_block.data(), ADD_VALUES);
        }
    }

    MatAssemblyBegin(this->petsc_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(this->petsc_mat, MAT_FINAL_ASSEMBLY);

    return;
}

/**
 * @brief Copy the owned RHS entries into the PETSc right-hand side vector.
 * @param ndof Degrees of freedom per node.
 */
void CrsMat::UpdatePetscRhs(int ndof) {
    VecZeroEntries(this->petsc_b);

    size_t buf_size = static_cast<size_t>(this->local_node) * ndof;
    this->petsc_indices_buf.resize(buf_size);
    this->petsc_values_buf.resize(buf_size);

    size_t idx = 0;
    for (int i = 0; i < this->local_node; i++) {
        int natural_id = this->interior_list[i];
        for (int var = 0; var < ndof; var++) {
            int local_idx = natural_id + var * nodec;
            this->petsc_indices_buf[idx] = this->l2g_var_map[local_idx];
            this->petsc_values_buf[idx] = this->b_rhs[local_idx];
            ++idx;
        }
    }

    VecSetValues(this->petsc_b, static_cast<PetscInt>(buf_size), this->petsc_indices_buf.data(),
                 this->petsc_values_buf.data(), INSERT_VALUES);
    VecAssemblyBegin(this->petsc_b);
    VecAssemblyEnd(this->petsc_b);

    return;
}

/**
 * @brief Solve the distributed linear system with PETSc KSP.
 * @param ndof Degrees of freedom per node.
 * @param nr_it Newton iteration index, or -1 for non-Newton solves.
 * @return Number of Krylov iterations used by KSPSolve.
 */
int CrsMat::SolveWithPetsc(int ndof, int nr_it) {
    constexpr int kSolidRebuildIterTrigger = 5;

    this->UpdatePetscRhs(ndof);

    // Initialise the PETSc guess vector from the current Newton iterate.
    VecZeroEntries(this->petsc_x);

    size_t buf_size = static_cast<size_t>(this->local_node) * ndof;
    this->petsc_indices_buf.resize(buf_size);
    this->petsc_values_buf.resize(buf_size);

    size_t idx = 0;
    for (int i = 0; i < this->local_node; i++) {
        int natural_id = this->interior_list[i];
        for (int var = 0; var < ndof; var++) {
            int local_idx = natural_id + var * nodec;
            this->petsc_indices_buf[idx] = this->l2g_var_map[local_idx];
            this->petsc_values_buf[idx] = this->x_lhs[local_idx];
            ++idx;
        }
    }

    VecSetValues(this->petsc_x, static_cast<PetscInt>(buf_size), this->petsc_indices_buf.data(),
                 this->petsc_values_buf.data(), INSERT_VALUES);
    VecAssemblyBegin(this->petsc_x);
    VecAssemblyEnd(this->petsc_x);

    this->BuildPetscBCList();
    if (!this->petsc_bc_gids.empty()) {
        MatZeroRowsColumns(this->petsc_mat, static_cast<PetscInt>(this->petsc_bc_gids.size()),
                           this->petsc_bc_gids.data(), 1.0e0, this->petsc_x, this->petsc_b);
    }

    KSPSetInitialGuessNonzero(this->ksp, PETSC_TRUE);

    // ----- AMG rebuild decision -----
    PC pc;
    KSPGetPC(this->ksp, &pc);

    bool need_rebuild = false;
    if (this->amg_rebuild_freq > 0) { need_rebuild = (istep % this->amg_rebuild_freq == 0); }

    const bool force_rebuild = this->force_rebuild_next_;

    // Solid passes nr_it >= 0. In that case we only rebuild on the first
    // Newton iteration unless a previous solve explicitly requested a
    // rebuild. Fluid keeps nr_it == -1 and uses only the step-based / forced
    // policy below.
    if (nr_it > 0 && !force_rebuild) {
        need_rebuild = false;
        if (!this->FEM_flag && this->prev_ksp_its_ >= kSolidRebuildIterTrigger) { need_rebuild = true; }
    }

    need_rebuild = need_rebuild || force_rebuild;
    this->force_rebuild_next_ = false;

    // Keep fluid / solid preconditioners independent. When this solver needs
    // a rebuild, explicitly reset only its own PC first so the old AMG
    // hierarchy is released before PETSc constructs the new one. This avoids
    // the larger peak-memory spike of "old AMG + new AMG" inside one solver
    // while still allowing fluid and solid AMGs to coexist.
    if (need_rebuild) {
        PCReset(pc);
        if (this->use_fieldsplit && this->ndof == 4) { this->ConfigureVelocityPressureFieldSplit(pc); }
        KSPSetOperators(this->ksp, this->petsc_mat, this->petsc_mat);
        PCSetReusePreconditioner(pc, PETSC_FALSE);
    } else {
        PCSetReusePreconditioner(pc, PETSC_TRUE);
    }

    KSPSolve(this->ksp, this->petsc_b, this->petsc_x);

    PetscInt its;
    KSPGetIterationNumber(this->ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(this->ksp, &reason);
    if (reason < 0 && myrank == 0) { std::cout << "PETSc KSP diverged, reason: " << reason << std::endl; }

    // Gather only the nodec*ndof entries needed by this rank (owned + ghost),
    // instead of materialising the full global solution on every rank.
    VecScatterBegin(this->scatter_to_all, this->petsc_x, this->seq_x, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(this->scatter_to_all, this->petsc_x, this->seq_x, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *seq_x_arr;
    VecGetArrayRead(this->seq_x, &seq_x_arr);
    for (int n = 0; n < nodec; ++n) {
        for (int var = 0; var < ndof; ++var) { this->x_lhs[n + var * nodec] = seq_x_arr[n + var * nodec]; }
    }
    VecRestoreArrayRead(this->seq_x, &seq_x_arr);

    // If KSP iterations grew by more than 2x compared with the previous step,
    // force a rebuild on the next step.  This applies to both fluid and solid.
    if (this->prev_ksp_its_ >= 0 && static_cast<double>(its) > this->prev_ksp_its_ * 2.0) {
        this->force_rebuild_next_ = true;
    }
    this->prev_ksp_its_ = static_cast<int>(its);

    return static_cast<int>(its);
}
