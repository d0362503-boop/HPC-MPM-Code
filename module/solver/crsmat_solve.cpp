#include "module/solver/crsmat.h"

#include <petsc.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solver/native_solver.h"

int CrsMat::SolveSystem(int NR_it) {

    int iter;
    if (this->use_petsc) {
        this->AssemblePetscMat();
        iter = this->SolveWithPetsc(NR_it);
    } else {
        this->BuildDiagonalPreconditioner();
        iter = GPBiCGAR(*this);
    }

    return iter;
}

void CrsMat::BuildActiveRowMask() {

    this->active_row_mask.assign(nodec, 1);
    if (this->FEM_flag) return;

    std::vector<int> local_row_active(nodec, 0);
    for (int nid = 0; nid < nodec; ++nid) { local_row_active[nid] = this->owner_->nmass[nid] > mtol ? 1 : 0; }

    NodeVarComm(local_row_active, 0);

    for (int nid = 0; nid < nodec; ++nid) { this->active_row_mask[nid] = (local_row_active[nid] > 0) ? 1 : 0; }

    return;
}

void CrsMat::AssemblePetscMat() { this->owner_->AssemblePetscMat(*this); }

void MaterialPoint::AssemblePetscMat(CrsMat &mat) {
    MatZeroEntries(mat.petsc_mat);
    mat.BuildActiveRowMask();

    for (int i = 0; i < nodec; ++i) {
        bool is_inactive = (!mat.FEM_flag && mat.active_row_mask[i] == 0);
        if (is_inactive) continue;

        int natural_row = i;
        int row_start = mat.matrow[natural_row];
        int row_end = mat.matrow[natural_row + 1];
        int ncols = row_end - row_start;

        const PetscInt block_row = mat.NaturalNodeToPetscLocalBlock(natural_row);

        size_t block_col_idx = 0;
        for (int j = row_start; j < row_end; ++j) {
            const int natural_col = mat.matcolid[j];
            mat.petsc_block_cols_buf[block_col_idx++] = mat.NaturalNodeToPetscLocalBlock(natural_col);
        }

        size_t val_idx = 0;
        for (int row_var = 0; row_var < mat.ndof; ++row_var) {
            for (int j = row_start; j < row_end; ++j) {
                for (int col_var = 0; col_var < mat.ndof; ++col_var) {
                    mat.petsc_block_vals_buf[val_idx] = mat.amat[j + mat.block_id[row_var * mat.ndof + col_var]];
                    ++val_idx;
                }
            }
        }

        MatSetValuesBlockedLocal(mat.petsc_mat, 1, &block_row, ncols, mat.petsc_block_cols_buf.data(),
                                 mat.petsc_block_vals_buf.data(), ADD_VALUES);
    }

    MatAssemblyBegin(mat.petsc_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat.petsc_mat, MAT_FINAL_ASSEMBLY);

    return;
}

void CrsMat::BuildPetscBCList() {
    this->petsc_bc_gids.clear();

    this->owner_->BuildPetscBCList(*this);

    std::sort(this->petsc_bc_gids.begin(), this->petsc_bc_gids.end());
    this->petsc_bc_gids.erase(std::unique(this->petsc_bc_gids.begin(), this->petsc_bc_gids.end()), this->petsc_bc_gids.end());

    return;
}

void CrsMat::SetPetscVector(Vec vec, const std::vector<double> &values) {
    VecZeroEntries(vec);

    size_t buf_size = static_cast<size_t>(this->local_node) * this->ndof;
    this->petsc_indices_buf.resize(buf_size);
    this->petsc_values_buf.resize(buf_size);

    size_t idx = 0;
    for (int i = 0; i < this->local_node; ++i) {
        const int natural_id = this->owned_natural_ids[i];
        for (int var = 0; var < this->ndof; ++var) {
            this->petsc_indices_buf[idx] = this->NaturalNodeVarToPetscLocalScalar(natural_id, var);
            this->petsc_values_buf[idx] = values[natural_id + var * nodec];
            ++idx;
        }
    }

    VecSetValuesLocal(vec, static_cast<PetscInt>(buf_size), this->petsc_indices_buf.data(), this->petsc_values_buf.data(),
                      INSERT_VALUES);
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);

    return;
}

int CrsMat::SolveWithPetsc(int NR_it) {
    this->BuildPetscBCList();
    this->SetPetscVector(this->petsc_b, this->b_rhs);
    // Prescribed values for elimination: zero for MPM corrections, current solution for FEM.
    this->SetPetscVector(this->petsc_x, this->x_lhs);

    std::vector<PetscInt> fixed_gids = this->petsc_bc_gids;
    if (!this->FEM_flag) {
        for (int i = 0; i < this->local_node; ++i) {
            const int natural_id = this->owned_natural_ids[i];
            if (this->active_row_mask[natural_id] != 0) continue;
            for (int var = 0; var < this->ndof; ++var) { fixed_gids.push_back(this->natural_var_gids[natural_id + var * nodec]); }
        }
        std::sort(fixed_gids.begin(), fixed_gids.end());
        fixed_gids.erase(std::unique(fixed_gids.begin(), fixed_gids.end()), fixed_gids.end());
    }

    // Collective even when this rank has no fixed rows.
    MatZeroRowsColumns(this->petsc_mat, static_cast<PetscInt>(fixed_gids.size()), fixed_gids.empty() ? nullptr : fixed_gids.data(),
                       1.0e0, this->petsc_x, this->petsc_b);

    if (!this->FEM_flag) {
        PetscCallAbort(MPI_COMM_WORLD, this->SolveReducedSystem(this->petsc_mat, this->petsc_b, this->petsc_x, NR_it));
    } else {
        // Start each correction from zero.
        KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE);
        VecZeroEntries(this->petsc_x);

        PC pc;
        KSPGetPC(this->ksp, &pc);
        const bool need_rebuild = this->NeedsPreconditionerRebuild(NR_it);
        if (need_rebuild) {
            PetscBool is_fieldsplit;
            PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_fieldsplit);
            if (!is_fieldsplit) { PCReset(pc); }
            KSPSetOperators(this->ksp, this->petsc_mat, this->petsc_mat);
            if (!is_fieldsplit) { this->ConfigurePreconditioner(pc); }
            PCSetReusePreconditioner(pc, PETSC_FALSE);
        } else {
            KSPSetOperators(this->ksp, this->petsc_mat, this->petsc_mat);
            PCSetReusePreconditioner(pc, PETSC_TRUE);
        }
        PetscCallAbort(MPI_COMM_WORLD, KSPSolve(this->ksp, this->petsc_b, this->petsc_x));
    }

    PetscInt its;
    KSPGetIterationNumber(this->ksp, &its);

    KSPConvergedReason reason;
    KSPGetConvergedReason(this->ksp, &reason);
    if (reason < 0 && myrank == 0) { std::cout << "PETSc KSP diverged, reason: " << reason << "\n"; }

    this->RecoverSolution();

    // If KSP iterations grew by more than 2x compared with the previous step,
    // force a rebuild on the next step.  This applies to both fluid and solid.
    if (this->rebuild_on_iteration_growth && this->prev_ksp_its_ >= 0 && static_cast<double>(its) > this->prev_ksp_its_ * 2.0) {
        this->force_rebuild_next_ = true;
    }
    this->prev_ksp_its_ = static_cast<int>(its);

    return static_cast<int>(its);
}

bool CrsMat::NeedsPreconditionerRebuild(int NR_it) {
    const bool force_rebuild = this->force_rebuild_next_;
    this->force_rebuild_next_ = false;
    if (NR_it > 0 && !force_rebuild) return false;
    return force_rebuild || (this->rebuild_freq > 0 && istep % this->rebuild_freq == 0);
}

void CrsMat::ResetReducedSystem() {
    if (this->reduced_active_is_) { ISDestroy(&this->reduced_active_is_); }
    if (this->reduced_matrix_) { MatDestroy(&this->reduced_matrix_); }
    this->retained_dof_gids_.clear();
}

PetscErrorCode CrsMat::BuildReducedDofMap(Mat a, Vec b, bool &changed) {
    PetscFunctionBeginUser;
    PetscInt lo, hi;
    PetscCall(MatGetOwnershipRange(a, &lo, &hi));
    const PetscScalar *rhs;
    PetscCall(VecGetArrayRead(b, &rhs));
    std::vector<PetscInt> retained;
    for (PetscInt row = lo; row < hi; ++row) {
        PetscInt count;
        const PetscInt *columns;
        const PetscScalar *values;
        PetscCall(MatGetRow(a, row, &count, &columns, &values));
        // Remove only equations enforcing an exactly zero solution: unit diagonal, zero off-diagonals, zero RHS.
        // Zeroed rows keep their explicit zeros (MAT_KEEP_NONZERO_PATTERN), so check values, not the entry count.
        bool zero_increment = rhs[row - lo] == PetscScalar(0.0);
        for (PetscInt k = 0; k < count && zero_increment; ++k) {
            zero_increment = values[k] == (columns[k] == row ? PetscScalar(1.0) : PetscScalar(0.0));
        }
        if (!zero_increment) { retained.push_back(row); }
        PetscCall(MatRestoreRow(a, row, &count, &columns, &values));
    }
    PetscCall(VecRestoreArrayRead(b, &rhs));

    int local_changed = !this->reduced_active_is_ || retained != this->retained_dof_gids_;
    int global_changed;
    PetscCallMPI(MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
    changed = global_changed != 0;
    if (changed) {
        if (this->reduced_active_is_) { PetscCall(ISDestroy(&this->reduced_active_is_)); }
        this->retained_dof_gids_.swap(retained);
        PetscCall(ISCreateGeneral(MPI_COMM_WORLD, this->retained_dof_gids_.size(), this->retained_dof_gids_.data(),
                                  PETSC_COPY_VALUES, &this->reduced_active_is_));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CrsMat::SolveReducedSystem(Mat a, Vec b, Vec x, int NR_it) {
    PetscFunctionBeginUser;
    bool layout_changed;
    PetscCall(this->BuildReducedDofMap(a, b, layout_changed));
    const bool scheduled_rebuild = this->NeedsPreconditionerRebuild(NR_it);
    const bool rebuild = layout_changed || !this->reduced_matrix_ || scheduled_rebuild;
    if (this->reduced_matrix_) { PetscCall(MatDestroy(&this->reduced_matrix_)); }
    PetscCall(MatCreateSubMatrix(a, this->reduced_active_is_, this->reduced_active_is_, MAT_INITIAL_MATRIX,
                                 &this->reduced_matrix_));
    PetscCall(MatEliminateZeros(this->reduced_matrix_, PETSC_TRUE));
    Vec reduced_rhs, reduced_x;
    PetscCall(VecGetSubVector(b, this->reduced_active_is_, &reduced_rhs));
    PetscCall(MatCreateVecs(this->reduced_matrix_, &reduced_x, nullptr));
    PetscCall(VecZeroEntries(reduced_x));

    if (rebuild) {
        PetscReal rtol, atol, dtol;
        PetscInt max_it;
        const char *old_type, *old_prefix;
        PetscCall(KSPGetTolerances(this->ksp, &rtol, &atol, &dtol, &max_it));
        PetscCall(KSPGetType(this->ksp, &old_type));
        PetscCall(KSPGetOptionsPrefix(this->ksp, &old_prefix));
        const std::string type = old_type ? old_type : KSPFGMRES;
        const std::string prefix = old_prefix ? old_prefix : "";
        PetscCall(KSPDestroy(&this->ksp));
        PetscCall(KSPCreate(MPI_COMM_WORLD, &this->ksp));
        PetscCall(KSPSetType(this->ksp, type.c_str()));
        PetscCall(KSPSetOptionsPrefix(this->ksp, prefix.c_str()));
        PetscCall(KSPSetTolerances(this->ksp, rtol, atol, dtol, max_it));
        PetscCall(this->ConfigureReducedPreconditioner());
    } else {
        PetscCall(KSPSetOperators(this->ksp, this->reduced_matrix_, this->reduced_matrix_));
        PetscCall(KSPSetReusePreconditioner(this->ksp, PETSC_FALSE));
    }
    PetscCall(KSPSetInitialGuessNonzero(this->ksp, PETSC_FALSE));
    PetscCall(KSPSolve(this->ksp, reduced_rhs, reduced_x));
    PetscCall(VecRestoreSubVector(b, this->reduced_active_is_, &reduced_rhs));

    // Keep the full vector/scatter layout, including ghost control points.
    PetscCall(VecZeroEntries(x));
    Vec retained_x;
    PetscCall(VecGetSubVector(x, this->reduced_active_is_, &retained_x));
    PetscCall(VecCopy(reduced_x, retained_x));
    PetscCall(VecRestoreSubVector(x, this->reduced_active_is_, &retained_x));
    PetscCall(VecDestroy(&reduced_x));

    PetscFunctionReturn(PETSC_SUCCESS);
}

void CrsMat::ConfigurePreconditioner(PC pc) {

    this->field_splits.clear();
    this->owner_->ConfigurePreconditioner(*this, pc);
    if (this->HasReducedLayout()) {
        for (const auto &split : this->field_splits) {
            PetscCallAbort(MPI_COMM_WORLD, this->SetFieldSplitComponents(pc, split.name, split.components));
        }
    }

    return;
}

PetscErrorCode CrsMat::ConfigureReducedPreconditioner() {
    PetscFunctionBeginUser;
    PC pc;
    PetscInt active_count;
    PetscCall(KSPGetPC(this->ksp, &pc));
    PetscCall(MatGetSize(this->reduced_matrix_, &active_count, nullptr));
    if (active_count) {
        this->ConfigurePreconditioner(pc);
    } else {
        PetscCall(PCSetType(pc, PCNONE));
    }
    PetscCall(KSPSetOperators(this->ksp, this->reduced_matrix_, this->reduced_matrix_));
    PetscCall(KSPSetFromOptions(this->ksp));
    if (!active_count) { PetscCall(PCSetType(pc, PCNONE)); }
    PetscCall(KSPSetUp(this->ksp));
    if (active_count) { PetscCall(this->ConfigureNestedFieldSplits(this->ksp, this->field_splits)); }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CrsMat::ConfigureNestedFieldSplits(KSP solver, const std::vector<FieldSplit> &splits) {
    PetscFunctionBeginUser;
    if (std::none_of(splits.begin(), splits.end(), [](const FieldSplit &split) { return !split.children.empty(); })) {
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PC pc;
    KSP *sub_solvers;
    PetscInt count;
    PetscCall(KSPGetPC(solver, &pc));
    PetscCall(PCFieldSplitGetSubKSP(pc, &count, &sub_solvers));
    PetscCheck(count == static_cast<PetscInt>(splits.size()), MPI_COMM_WORLD, PETSC_ERR_ARG_SIZ,
               "Physical field groups must match PETSc sub-solvers");
    for (PetscInt i = 0; i < count; ++i) {
        const auto &split = splits[i];
        if (split.children.empty()) continue;
        KSP inner = sub_solvers[i];
        PC inner_pc;
        PetscBool is_ksp;
        PetscCall(KSPGetPC(inner, &inner_pc));
        PetscCall(PetscObjectTypeCompare((PetscObject)inner_pc, PCKSP, &is_ksp));
        if (is_ksp) {
            Mat pmat;
            PetscCall(PCSetUseAmat(inner_pc, PETSC_FALSE));
            PetscCall(PCKSPGetKSP(inner_pc, &inner));
            PetscCall(KSPGetOperators(sub_solvers[i], nullptr, &pmat));
            PetscCall(KSPSetOperators(inner, pmat, pmat));
            PetscCall(KSPSetFromOptions(inner));
            PetscCall(KSPGetPC(inner, &inner_pc));
        }
        for (const auto &child : split.children) {
            PetscCall(this->SetFieldSplitComponents(inner_pc, child.name, child.components, split.components));
        }
        PetscCall(KSPSetUp(inner));
        PetscCall(this->ConfigureNestedFieldSplits(inner, split.children));
    }
    PetscCall(PetscFree(sub_solvers));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CrsMat::SetFieldSplitComponents(PC pc, const char *name, const std::vector<int> &components,
                                               const std::vector<int> &parent_components) const {
    PetscFunctionBeginUser;
    const auto contains = [](const std::vector<int> &fields, int component) {
        return std::find(fields.begin(), fields.end(), component) != fields.end();
    };
    PetscInt parent_count = 0;
    std::vector<PetscInt> indices;
    for (PetscInt original : this->retained_dof_gids_) {
        const int component = original % this->ndof;
        if (!parent_components.empty() && !contains(parent_components, component)) continue;
        if (contains(components, component)) { indices.push_back(parent_count); }
        ++parent_count;
    }
    PetscInt offset = 0;
    PetscCallMPI(MPI_Exscan(&parent_count, &offset, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD));
    if (myrank == 0) { offset = 0; }
    for (PetscInt &index : indices) { index += offset; }
    IS field_is;
    PetscCall(ISCreateGeneral(MPI_COMM_WORLD, indices.size(), indices.data(), PETSC_COPY_VALUES, &field_is));
    PetscCall(PCFieldSplitSetIS(pc, name, field_is));
    PetscCall(ISDestroy(&field_is));
    PetscFunctionReturn(PETSC_SUCCESS);
}

void CrsMat::RecoverSolution() {
    // Gather only the nodec*this->ndof entries needed by this rank (owned + ghost),
    // instead of materialising the full global solution on every rank.
    VecScatterBegin(this->scatter_to_all, this->petsc_x, this->seq_x, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(this->scatter_to_all, this->petsc_x, this->seq_x, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar *seq_x_arr;
    VecGetArrayRead(this->seq_x, &seq_x_arr);
    std::copy(seq_x_arr, seq_x_arr + nodec * this->ndof, this->x_lhs.begin());
    VecRestoreArrayRead(this->seq_x, &seq_x_arr);
}

void CrsMat::ExtractDiagonal() {

    for (int i = 0; i < nodec; i++) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; j++) {
            if (i == this->matcolid[j]) {
                for (int b = 0; b < this->num_block; b++) {
                    if (this->block_row[b] != this->block_col[b]) continue;
                    const int offset = nodec * this->block_row[b];
                    this->adiag[i + offset] = this->amat[j + this->block_id[b]];
                }
                break;
            }
        }
    }

    std::vector<int> offsets(this->ndof);
    for (int n = 0; n < this->ndof; ++n) { offsets[n] = nodec * n; }
    NodeVarComm(this->adiag, offsets);

    return;
}

void CrsMat::ComputeResidualNormSq(double &ref_tol, double &abs_tol) {

    const int var_size = int(this->x_lhs.size());

    std::vector<double> free_dof(var_size, 1.0);
    this->owner_->BCResidualSet(this->b_rhs);
    this->owner_->BCResidualSet(free_dof);

    double rtr = 0.0e0, active_dof = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        if (free_dof[n] == 0.0) { continue; }
        int nid = n % nodec;
        if (!this->FEM_flag) {
            if (this->owner_->nmass[nid] > mtol) {
                double r = this->b_rhs[n];
                rtr += r * r * dbc[n];
                active_dof += dbc[n];
            }
        } else {
            double r = this->b_rhs[n];
            rtr += r * r * dbc[n];
            active_dof += dbc[n];
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &rtr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &active_dof, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (active_dof <= 0.0e0) {
        abs_tol = 0.0e0;
        ref_tol = 0.0e0;
    } else {
        abs_tol = std::sqrt(rtr / active_dof);
        ref_tol = std::sqrt(rtr);
    }

    return;
}

bool CrsMat::CheckNRConvergence(int NR_it, int NR_it_max, int solver_it, double &r0r) {

    double rkr, rtr_abs;
    this->ComputeResidualNormSq(rkr, rtr_abs);
    if (NR_it == 0) { r0r = rkr; }

    const double rtr_ref = (r0r > 1.0e-30) ? (rkr / r0r) : 0.0e0;

    if (rtr_ref < 1.0e-8 || rtr_abs < 1.0e-12) {
        if (myrank == 0) {
            std::cout << "NR_converge:" << std::setw(15) << NR_it << std::setw(15) << solver_it << std::setw(15)
                      << std::scientific << rtr_ref << std::setw(15) << r0r << std::setw(15) << rtr_abs << "\n";
        }
        return true;
    } else if (NR_it == NR_it_max) {
        if (myrank == 0) {
            std::cout << "NR_it:" << std::setw(15) << NR_it << "\n";
            std::cout << "NR can not converge" << "\n";
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    return false;
}
