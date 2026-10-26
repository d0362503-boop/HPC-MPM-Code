#include "module/solver/crsmat.h"

#include <algorithm>
#include <string>

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
        // Remove only equations enforcing an exactly zero solution.
        const bool zero_increment = count == 1 && columns[0] == row && values[0] == PetscScalar(1.0) &&
                                    rhs[row - lo] == PetscScalar(0.0);
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

PetscErrorCode CrsMat::SetFieldSplitComponents(PC pc, const char *name, const std::vector<int> &components,
                                               const std::vector<int> &parent_components) const {
    PetscFunctionBeginUser;
    const auto contains = [](const std::vector<int> &fields, int component) {
        return std::find(fields.begin(), fields.end(), component) != fields.end();
    };
    if (!this->HasReducedLayout()) {
        PetscCheck(parent_components.empty(), MPI_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                   "Nested component splits require an initialized reduced layout");
        std::vector<PetscInt> fields(components.begin(), components.end());
        PetscCall(PCFieldSplitSetFields(pc, name, fields.size(), fields.data(), fields.data()));
        PetscFunctionReturn(PETSC_SUCCESS);
    }

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

PetscErrorCode CrsMat::SolveReducedSystem(Mat a, Vec b, Vec x, int NR_it) {
    PetscFunctionBeginUser;
    bool layout_changed;
    PetscCall(this->BuildReducedDofMap(a, b, layout_changed));
    const bool scheduled_rebuild = this->NeedsPreconditionerRebuild(NR_it);
    const bool rebuild = layout_changed || !this->reduced_matrix_ || scheduled_rebuild;
    if (this->reduced_matrix_) { PetscCall(MatDestroy(&this->reduced_matrix_)); }
    PetscCall(MatCreateSubMatrix(a, this->reduced_active_is_, this->reduced_active_is_, MAT_INITIAL_MATRIX,
                                 &this->reduced_matrix_));
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
