#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

#include <algorithm>
#include <cmath>

using namespace mpm_mpm_monolithic_fsi;

int MPMMPMMonolithicFSI::SolveSystem(int NR_it) {

    this->BCResidualSet(this->fsi_sys.b_rhs);
    std::vector<double> residual = this->fsi_sys.b_rhs;

    VectorAssign(nodec * 10, this->fsi_sys.x_lhs);
    this->ScaleSystem();

    // Rebuild for the current scaled Newton matrix.
    this->fsi_sys.force_rebuild_next_ = true;
    KSPSetTolerances(this->fsi_sys.ksp, 1.0e-7, 1.0e-15, PETSC_CURRENT, PETSC_CURRENT);
    const int iter = this->fsi_sys.SolveSystem(NR_it);

    KSPConvergedReason reason;
    KSPGetConvergedReason(this->fsi_sys.ksp, &reason);
    if (reason < 0) { MPI_Abort(MPI_COMM_WORLD, 1); }

    for (int n = 0; n < nodec * 10; n++) { this->fsi_sys.x_lhs[n] *= this->column_scale[n]; }
    this->fsi_sys.b_rhs.swap(residual);

    return iter;
}

MPMMPMMonolithicFSI::~MPMMPMMonolithicFSI() { this->fsi_sys.ResetPetscSolver(); }

void MPMMPMMonolithicFSI::BuildActiveDOFs() {

    double mass_stats[4]{};
    for (int n : this->fsi_sys.owned_natural_ids) {
        if (this->fluid_.nmass[n] > mtol) {
            mass_stats[0] += this->fluid_.nmass[n];
            mass_stats[1] += 1.0;
        }
        if (this->solid_.nmass[n] > mtol) {
            mass_stats[2] += this->solid_.nmass[n];
            mass_stats[3] += 1.0;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, mass_stats, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double fluid_cut = 1.0e-4 * mass_stats[0] / mass_stats[1];
    const double solid_cut = 1.0e-4 * mass_stats[2] / mass_stats[3];
    this->fixed_dof.assign(nodec * 10, 0);
    for (int n = 0; n < nodec; n++) {
        const bool fluid_active = this->fluid_.nmass[n] > fluid_cut;
        const bool solid_active = this->solid_.nmass[n] > solid_cut;
        for (int d = 0; d < 4; d++) { this->fixed_dof[n + d * nodec] = !fluid_active; }
        for (int d = 4; d < 7; d++) { this->fixed_dof[n + d * nodec] = !solid_active; }
        for (int d = 7; d < 10; d++) {
            this->fixed_dof[n + d * nodec] = !fluid_active || !solid_active || this->nlm_lump[n] < mtol;
        }
    }

    const BoundaryCondition *bc[] = {&this->fluid_.ubc, &this->fluid_.vbc, &this->fluid_.wbc, &this->fluid_.pbc,
                                     &this->solid_.ubc, &this->solid_.vbc, &this->solid_.wbc};
    for (int d = 0; d < 7; d++) {
        for (int i = 0; i < bc[d]->ibc; i++) { this->fixed_dof[bc[d]->nbc[i] + d * nodec] = 1; }
    }

    for (int i = 0; i < this->solid_.rigid_bc.ibc; i++) {
        const int n = this->solid_.rigid_bc.nbc[i];
        for (int d = 4; d < 7; d++) { this->fixed_dof[n + d * nodec] = 1; }
    }
    for (int n = 0; n < nodec; n++) {
        for (int d = 0; d < 3; d++) {
            if (this->fixed_dof[n + d * nodec] && this->fixed_dof[n + (d + 4) * nodec]) {
                this->fixed_dof[n + (d + 7) * nodec] = 1;
            }
        }
    }
}

void MPMMPMMonolithicFSI::BuildPetscBCList(CrsMat &mat) {

    for (int n : mat.owned_natural_ids) {
        for (int d = 0; d < 10; d++) {
            if (this->fixed_dof[n + d * nodec]) { mat.petsc_bc_gids.push_back(mat.natural_var_gids[n + d * nodec]); }
        }
    }
}

void MPMMPMMonolithicFSI::ScaleSystem() {

    std::vector<double> diagonal(nodec * 7, 0.0);
    for (int n = 0; n < nodec; n++) {
        int ncol = 0;
        const int j = this->fsi_sys.FindIndex(n, n, ncol);
        for (int b = 0; b < this->fsi_sys.num_block; b++) {
            if (this->fsi_sys.block_row[b] == this->fsi_sys.block_col[b]) {
                diagonal[n + this->fsi_sys.block_row[b] * nodec] = this->fsi_sys.amat[j + this->fsi_sys.block_id[b]];
            }
        }
    }

    NodeVarComm(diagonal, {0, nodec, 2 * nodec, 3 * nodec, 4 * nodec, 5 * nodec, 6 * nodec});

    this->column_scale.assign(nodec * 10, 1.0);
    std::vector<double> row_scale(nodec * 10, 1.0);

    for (int n = 0; n < nodec; n++) {
        for (int d = 0; d < 7; d++) {
            const int i = n + d * nodec;
            if (!this->fixed_dof[i]) { row_scale[i] = this->column_scale[i] = 1.0 / std::sqrt(std::abs(diagonal[i])); }
        }
        for (int d = 0; d < 3; d++) {
            const int i = n + (d + 7) * nodec;
            if (this->fixed_dof[i]) continue;
            const double df = this->fixed_dof[n + d * nodec] ? 0.0 : 1.0 / std::abs(diagonal[n + d * nodec]);
            const double ds = this->fixed_dof[n + (d + 4) * nodec] ? 0.0 : 1.0 / std::abs(diagonal[n + (d + 4) * nodec]);
            const double af = this->fluid_.alpha_f, as = this->solid_.alpha_f;
            const double cf = this->fluid_.nb_para[0], cs = this->solid_.nb_para[0];
            this->column_scale[i] = 1.0 / (this->nlm_lump[n] * std::sqrt(af * af * df + as * as * ds));
            row_scale[i] = 1.0 / (this->nlm_lump[n] * std::sqrt(cf * cf * df + cs * cs * ds));
        }
    }

    for (int b = 0; b < this->fsi_sys.num_block; b++) {
        for (int n = 0; n < nodec; n++) {
            for (int j = this->fsi_sys.matrow[n]; j < this->fsi_sys.matrow[n + 1]; j++) {
                this->fsi_sys.amat[j + this->fsi_sys.block_id[b]] *=
                    row_scale[n + this->fsi_sys.block_row[b] * nodec] *
                    this->column_scale[this->fsi_sys.matcolid[j] + this->fsi_sys.block_col[b] * nodec];
            }
        }
    }

    for (int n = 0; n < nodec * 10; n++) { this->fsi_sys.b_rhs[n] *= row_scale[n]; }
}

void MPMMPMMonolithicFSI::BCResidualSet(std::vector<double> &rr) {

    for (int n = 0; n < nodec * 10; n++) {
        if (this->fixed_dof[n]) { rr[n] = 0.0; }
    }
}

void MPMMPMMonolithicFSI::AssemblePetscMat(CrsMat &mat, int ndof) {

    MatZeroEntries(mat.petsc_mat);
    mat.active_row_mask.assign(nodec, 1);
    std::vector<PetscInt> columns;
    std::vector<PetscScalar> values;

    for (int n = 0; n < nodec; n++) {
        for (int d = 0; d < ndof; d++) {
            if (this->fixed_dof[n + d * nodec]) continue;
            columns.clear();
            values.clear();
            for (int b = 0; b < mat.num_block; b++) {
                if (mat.block_row[b] != d) continue;
                for (int j = mat.matrow[n]; j < mat.matrow[n + 1]; j++) {
                    const double value = mat.amat[j + mat.block_id[b]];
                    if (value == 0.0) continue;
                    columns.push_back(mat.natural_var_gids[mat.matcolid[j] + mat.block_col[b] * nodec]);
                    values.push_back(value);
                }
            }
            const PetscInt row = mat.natural_var_gids[n + d * nodec];
            MatSetValues(mat.petsc_mat, 1, &row, columns.size(), columns.data(), values.data(), ADD_VALUES);
        }
    }

    MatAssemblyBegin(mat.petsc_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat.petsc_mat, MAT_FINAL_ASSEMBLY);
}

void MPMMPMMonolithicFSI::ConfigurePreconditioner(CrsMat &mat, PC pc) {

    const PetscInt fields[] = {0, 1, 2, 3, 4, 5, 6};
    const PetscInt multiplier[] = {7, 8, 9};

    KSPSetOptionsPrefix(mat.ksp, "fsi_");
    PCSetOptionsPrefix(pc, "fsi_");
    KSPGMRESSetRestart(mat.ksp, 100);

    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, 10);
    PCFieldSplitSetFields(pc, "fields", 7, fields, fields);
    PCFieldSplitSetFields(pc, "lambda", 3, multiplier, multiplier);
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);

    const char *options[][2] = {{"-fsi_ksp_gmres_cgs_refinement_type", "refine_always"},
                                {"-fsi_fieldsplit_fields_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fields_pc_type", "fieldsplit"},
                                {"-fsi_fieldsplit_fields_pc_fieldsplit_type", "additive"},
                                {"-fsi_fieldsplit_fields_pc_fieldsplit_block_size", "7"},
                                {"-fsi_fieldsplit_fields_pc_fieldsplit_0_fields", "0,1,2,3"},
                                {"-fsi_fieldsplit_fields_pc_fieldsplit_1_fields", "4,5,6"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_ksp_rtol", "0.5"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_ksp_max_it", "100"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_pc_type", "asm"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_pc_factor_zeropivot", "1e-8"},
                                {"-fsi_fieldsplit_fields_fieldsplit_0_sub_pc_factor_shift_amount", "1e-3"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_ksp_type", "preonly"},

                                {"-fsi_fieldsplit_fields_fieldsplit_1_pc_type", "asm"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fields_fieldsplit_1_sub_pc_factor_shift_amount", "1e-3"},

                                {"-fsi_fieldsplit_lambda_ksp_type", "gmres"},
                                {"-fsi_fieldsplit_lambda_ksp_rtol", "0.2"},
                                {"-fsi_fieldsplit_lambda_ksp_max_it", "5"},
                                {"-fsi_fieldsplit_lambda_ksp_gmres_restart", "5"},
                                {"-fsi_fieldsplit_lambda_pc_type", "asm"},
                                {"-fsi_fieldsplit_lambda_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_lambda_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_lambda_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_lambda_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_lambda_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_lambda_sub_pc_factor_shift_amount", "1e-3"}};

    for (const auto &option : options) {
        PetscBool supplied;
        PetscOptionsHasName(nullptr, nullptr, option[0], &supplied);
        if (!supplied) { PetscOptionsSetValue(nullptr, option[0], option[1]); }
    }
}
