#include "work/src_fsi/MPM_MPM/monolithic_fsi.h"

#include <algorithm>
#include <cmath>

using namespace mpm_mpm_monolithic_fsi;

void MPMMPMMonolithicFSI::ConfigurePreconditioner(CrsMat &mat, PC pc) {
    KSPSetOptionsPrefix(mat.ksp, "fsi_");
    PCSetOptionsPrefix(pc, "fsi_");
    KSPGMRESSetRestart(mat.ksp, 200);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, 1);
    mat.field_splits = {{"fluid", {0, 1, 2, 3}, {{"0", {0, 1, 2}}, {"1", {3}}}},
                        {"coupled", {4, 5, 6, 7, 8, 9}, {{"solid", {4, 5, 6}}, {"lambda", {7, 8, 9}}}}};
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_LOWER);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);

    // Fluid velocity-pressure preconditioner.
    const char *options[][2] = {{"-fsi_ksp_gmres_cgs_refinement_type", "refine_always"},
                                {"-fsi_fieldsplit_fluid_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_fluid_ksp_gmres_cgs_refinement_type", "refine_always"},
                                {"-fsi_fieldsplit_fluid_ksp_gmres_restart", "100"},
                                {"-fsi_fieldsplit_fluid_ksp_rtol", "0.5"},
                                {"-fsi_fieldsplit_fluid_ksp_max_it", "200"},
                                {"-fsi_fieldsplit_fluid_pc_type", "fieldsplit"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_type", "schur"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_schur_fact_type", "lower"},
                                {"-fsi_fieldsplit_fluid_pc_fieldsplit_schur_precondition", "selfp"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_pc_type", "asm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_mat_ordering_type", "rcm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_0_sub_pc_factor_shift_amount", "1e-3"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_pc_type", "asm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_zeropivot", "1e-6"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_mat_ordering_type", "rcm"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_fluid_fieldsplit_1_sub_pc_factor_shift_amount", "1e-3"},

                                // Joint solid-multiplier iteration.
                                {"-fsi_fieldsplit_coupled_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_pc_type", "ksp"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_rtol", "0.5"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_max_it", "100"},
                                {"-fsi_fieldsplit_coupled_ksp_ksp_gmres_restart", "100"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_type", "fieldsplit"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_type", "schur"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_schur_fact_type", "lower"},
                                {"-fsi_fieldsplit_coupled_ksp_pc_fieldsplit_schur_precondition", "selfp"},

                                // Solid ASM with local ILU(1).
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_type", "fgmres"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_gmres_cgs_refinement_type", "refine_always"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_gmres_restart", "30"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_rtol", "0.1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_ksp_max_it", "100"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_pc_type", "asm"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_mat_ordering_type", "rcm"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_solid_sub_pc_factor_shift_amount", "1e-3"},

                                // Multiplier Schur preconditioner.
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_pc_type", "asm"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_pc_asm_overlap", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_ksp_type", "preonly"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_type", "ilu"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_levels", "1"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_shift_type", "nonzero"},
                                {"-fsi_fieldsplit_coupled_ksp_fieldsplit_lambda_sub_pc_factor_shift_amount", "1e-3"}};

    for (const auto &option : options) {
        PetscBool supplied;
        PetscOptionsHasName(nullptr, nullptr, option[0], &supplied);
        if (!supplied) { PetscOptionsSetValue(nullptr, option[0], option[1]); }
    }
}

int MPMMPMMonolithicFSI::SolveSystem(int NR_it) {

    this->BCResidualSet(this->fsi_sys.b_rhs);
    std::vector<double> residual = this->fsi_sys.b_rhs;

    // Solve a new correction.
    VectorAssign(nodec * this->fsi_sys.ndof, this->fsi_sys.x_lhs);
    this->ScaleSystem();

    double active_dof = 0.0e0;
    for (int n = 0; n < nodec; n++) {
        if (this->nlm_lump[n] > 0.0e0) active_dof = 1.0e0;
    }

    MPI_Allreduce(MPI_IN_PLACE, &active_dof, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double ref_tol = 1.0e-8;
    double abs_tol = 1.0e-15;
    const double div_tol = NR_it == 0 ? 1.0e6 : PETSC_UNLIMITED;
    KSPSetTolerances(this->fsi_sys.ksp, ref_tol, abs_tol, div_tol, 500);
    const int iter = this->fsi_sys.SolveSystem(NR_it);

    KSPConvergedReason reason;
    KSPGetConvergedReason(this->fsi_sys.ksp, &reason);
    // if (reason < 0) { MPI_Abort(MPI_COMM_WORLD, 1); }

    for (int n = 0; n < nodec * this->fsi_sys.ndof; n++) { this->fsi_sys.x_lhs[n] *= this->column_scale[n]; }
    this->fsi_sys.b_rhs.swap(residual);

    return iter;
}

void MPMMPMMonolithicFSI::BuildActiveDOFs() {

    this->fixed_dof.assign(nodec * this->fsi_sys.ndof, 0);
    for (int n = 0; n < nodec; n++) {
        const bool fluid_active = this->fluid_.nmass[n] > mtol;
        const bool solid_active = this->solid_.nmass[n] > mtol;
        for (int d = 0; d < 4; d++) { this->fixed_dof[n + d * nodec] = !fluid_active; }
        for (int d = 4; d < 7; d++) { this->fixed_dof[n + d * nodec] = !solid_active; }
        for (int d = 7; d < this->fsi_sys.ndof; d++) {
            this->fixed_dof[n + d * nodec] = !fluid_active || !solid_active || this->nlm_lump[n] <= 0.0e0;
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
        for (int d = 0; d < this->fsi_sys.ndof; d++) {
            if (this->fixed_dof[n + d * nodec]) { mat.petsc_bc_gids.push_back(mat.natural_var_gids[n + d * nodec]); }
        }
    }
}

void MPMMPMMonolithicFSI::ScaleSystem() {

    this->fsi_sys.adiag.assign(nodec * this->fsi_sys.ndof, 0.0);
    this->fsi_sys.ExtractDiagonal();

    // Reuse nodal-mass activity mask.
    this->column_scale.assign(nodec * this->fsi_sys.ndof, 0.0);
    std::vector<double> row_scale(nodec * this->fsi_sys.ndof, 0.0);

    for (int n = 0; n < nodec; n++) {
        for (int d = 0; d < 7; d++) {
            const int i = n + d * nodec;
            if (!this->fixed_dof[i]) { row_scale[i] = this->column_scale[i] = 1.0 / std::sqrt(std::abs(this->fsi_sys.adiag[i])); }
        }
        for (int d = 0; d < 3; d++) {
            const int i = n + (d + 7) * nodec;
            if (this->fixed_dof[i]) continue;
            const double df = this->column_scale[n + d * nodec] * this->column_scale[n + d * nodec];
            const double ds = this->column_scale[n + (d + 4) * nodec] * this->column_scale[n + (d + 4) * nodec];
            const double af = this->fluid_.integrator_.alpha_f, as = this->solid_.integrator_.alpha_f;
            this->column_scale[i] = 1.0 / (this->nlm_lump[n] * std::sqrt(af * af * df + as * as * ds));
            row_scale[i] = 1.0 / (this->nlm_lump[n] * std::sqrt(df + ds));
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

    for (int n = 0; n < nodec * this->fsi_sys.ndof; n++) {
        this->fsi_sys.b_rhs[n] *= row_scale[n];
        // Keep zero correction initial guess.
        this->fsi_sys.x_lhs[n] = this->fixed_dof[n] ? 0.0 : this->fsi_sys.x_lhs[n] / this->column_scale[n];
    }
}

void MPMMPMMonolithicFSI::BCResidualSet(std::vector<double> &rr) {

    for (int n = 0; n < nodec * this->fsi_sys.ndof; n++) {
        if (this->fixed_dof[n]) { rr[n] = 0.0; }
    }
}
