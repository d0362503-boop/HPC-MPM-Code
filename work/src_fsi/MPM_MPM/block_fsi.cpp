#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/shape_function.h"

#include "work/src_fsi/MPM_MPM/block_fsi.h"

using namespace mpmmpmblockfsi;

MPMMPMBlockFSI::~MPMMPMBlockFSI() { KSPDestroy(&this->schur_fluid_ksp_); }

void MPMMPMBlockFSI::BuildFluidResponse(int block_it) {

    if (block_it > 0) {
        KSPSetOperators(this->schur_fluid_ksp_, this->fluid_.NS_.petsc_mat, this->fluid_.NS_.petsc_mat);
        KSPSetReusePreconditioner(this->schur_fluid_ksp_, PETSC_TRUE);
        return;
    }
    KSPDestroy(&this->schur_fluid_ksp_);
    KSPCreate(PETSC_COMM_WORLD, &this->schur_fluid_ksp_);
    KSPSetOperators(this->schur_fluid_ksp_, this->fluid_.NS_.petsc_mat, this->fluid_.NS_.petsc_mat);
    KSPSetType(this->schur_fluid_ksp_, KSPFGMRES);
    KSPSetOptionsPrefix(this->schur_fluid_ksp_, "fsi_response_");
    KSPSetTolerances(this->schur_fluid_ksp_, 1.0e-8, 1.0e-15, 1.0e6, 1000);
    PC pc;
    KSPGetPC(this->schur_fluid_ksp_, &pc);
    const PetscInt velocity_fields[] = {0, 1, 2};
    const PetscInt pressure_field = 3;
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetBlockSize(pc, 4);
    PCFieldSplitSetFields(pc, "velocity", 3, velocity_fields, velocity_fields);
    PCFieldSplitSetFields(pc, "pressure", 1, &pressure_field, &pressure_field);
    PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
    PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, nullptr);
    KSPSetUp(this->schur_fluid_ksp_);
    PetscInt split_count;
    KSP *response_split;
    PCFieldSplitGetSubKSP(pc, &split_count, &response_split);
    PC velocity_pc, pressure_pc;
    KSPSetType(response_split[0], KSPPREONLY);
    KSPGetPC(response_split[0], &velocity_pc);
    PCSetType(velocity_pc, PCPBJACOBI);
    KSPSetType(response_split[1], KSPPREONLY);
    KSPGetPC(response_split[1], &pressure_pc);
    PCSetType(pressure_pc, PCTELESCOPE);
    PCTelescopeSetReductionFactor(pressure_pc, nprocs);
    KSPSetUp(response_split[1]);
    KSP serial_solver;
    PCTelescopeGetKSP(pressure_pc, &serial_solver);
    if (myrank == 0) {
        Mat serial_tangent;
        KSPGetOperators(serial_solver, &serial_tangent, nullptr);
        MatEliminateZeros(serial_tangent, PETSC_TRUE);
        KSPSetDiagonalScale(serial_solver, PETSC_TRUE);
        KSPSetDiagonalScaleFix(serial_solver, PETSC_FALSE);
        PC serial_pc;
        KSPGetPC(serial_solver, &serial_pc);
        PCSetType(serial_pc, PCILU);
        PCFactorSetLevels(serial_pc, 2);
        PCFactorSetMatOrderingType(serial_pc, MATORDERINGRCM);
        PCFactorSetZeroPivot(serial_pc, 1.0e-30);
    }
    PetscFree(response_split);

    return;
}

void MPMMPMBlockFSI::BuildSolidResponse() {

    KSPCreate(PETSC_COMM_WORLD, &this->schur_solid_ksp_);
    KSPSetOperators(this->schur_solid_ksp_, this->solid_.SM_.petsc_mat, this->solid_.SM_.petsc_mat);
    KSPSetType(this->schur_solid_ksp_, KSPPREONLY);
    PC pc;
    KSPGetPC(this->schur_solid_ksp_, &pc);
    PCSetType(pc, PCTELESCOPE);
    PCTelescopeSetReductionFactor(pc, nprocs);
    KSPSetUp(this->schur_solid_ksp_);
    KSP serial_solver;
    PCTelescopeGetKSP(pc, &serial_solver);
    if (myrank == 0) {
        Mat serial_tangent;
        KSPGetOperators(serial_solver, &serial_tangent, nullptr);
        MatEliminateZeros(serial_tangent, PETSC_TRUE);
        PC serial_pc;
        KSPGetPC(serial_solver, &serial_pc);
        PCSetType(serial_pc, PCLU);
        PCFactorSetMatOrderingType(serial_pc, MATORDERINGND);
        KSPSetUp(serial_solver);
    }

    return;
}

PetscInt MPMMPMBlockFSI::SolveResponse(Mat coupling, KSP solver, Vec trial_multiplier, Vec load, Vec solution,
                                       Vec response) {

    MatMult(coupling, trial_multiplier, load);

    KSPSetInitialGuessNonzero(solver, PETSC_FALSE);
    KSPSetReusePreconditioner(solver, PETSC_TRUE);
    KSPSolve(solver, load, solution);

    PetscInt iterations;
    KSPGetIterationNumber(solver, &iterations);
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver, &reason);
    if (reason < 0 && myrank == 0) { std::cout << "PETSc KSP diverged, reason: " << reason << std::endl; }

    MatMultTranspose(coupling, solution, response);
    return iterations;
}

void MPMMPMBlockFSI::BuildCouplingOperator(CrsMat &mat, const std::vector<double> &weights,
                                           const std::vector<PetscInt> &local_interface_ids,
                                           PetscInt local_multiplier_dofs, PetscInt global_multiplier_dofs,
                                           Mat &coupling, Vec &load, Vec &solution) {

    PetscInt local_field_dofs, global_field_dofs;
    VecGetLocalSize(mat.petsc_b, &local_field_dofs);
    VecGetSize(mat.petsc_b, &global_field_dofs);
    MatCreateAIJ(PETSC_COMM_WORLD, local_field_dofs, local_multiplier_dofs, global_field_dofs, global_multiplier_dofs,
                 1, nullptr, 1, nullptr, &coupling);

    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (mat.natural_is_owned[nid] == 0) continue;

        for (int var = 0; var < 3; var++) {
            const PetscInt row = mat.natural_var_gids[nid + var * nodec];
            if (std::binary_search(mat.petsc_bc_gids.begin(), mat.petsc_bc_gids.end(), row)) continue;
            const PetscInt column = 3 * local_interface_ids[n] + var;
            MatSetValue(coupling, row, column, weights[nid], INSERT_VALUES);
        }
    }
    MatAssemblyBegin(coupling, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(coupling, MAT_FINAL_ASSEMBLY);
    VecDuplicate(mat.petsc_b, &load);
    VecDuplicate(mat.petsc_x, &solution);

    return;
}

PetscErrorCode MPMMPMBlockFSI::ApplyExactSchur(Vec trial_multiplier, Vec response) {

    this->schur_fluid_iterations_ +=
        this->SolveResponse(this->schur_fluid_coupling_, this->schur_fluid_ksp_, trial_multiplier,
                            this->schur_fluid_rhs_, this->schur_fluid_solution_, response);
    VecScale(response, this->fluid_.nb_para[0]);

    this->schur_solid_iterations_ +=
        this->SolveResponse(this->schur_solid_coupling_, this->schur_solid_ksp_, trial_multiplier,
                            this->schur_solid_rhs_, this->schur_solid_solution_, this->schur_solid_response_);
    VecAXPY(response, this->solid_.nb_para[0], this->schur_solid_response_);

    this->schur_matvec_count_++;
    return PETSC_SUCCESS;
}

void MPMMPMBlockFSI::DetectFSIInterface() {

    VectorAssign(nodec, this->fsi_intf.nbc);
    VectorAssign(nodec, this->solid_.nphi);
    VectorAssign(nodec, this->fluid_.nphi);

    for (int n = 0; n < nodec; n++) {
        this->solid_.nphi[n] = std::clamp(this->solid_.nvof[n] / nvol[n], 0.0e0, 1.0e0);
        this->fluid_.nphi[n] = std::clamp(this->fluid_.nvof[n] / nvol[n], 0.0e0, 1.0e0);
    }

    this->LumpedLagrangeMultiplier();
    this->fluid_.NS_.BuildActiveRowMask();
    this->solid_.SM_.BuildActiveRowMask();

    this->fsi_intf.ibc = 0;
    for (int n = 0; n < nodec; n++) {
        if (this->fluid_.lm_lumped[n] > mtol && this->fluid_.NS_.active_row_mask[n] != 0 &&
            this->solid_.SM_.active_row_mask[n] != 0) {
            this->fsi_intf.nbc[this->fsi_intf.ibc++] = n;
        }
    }

    VectorAssign(nodec * 3, this->nfsi_force);

    return;
}

void MPMMPMBlockFSI::LumpedLagrangeMultiplier() {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    double volp = dxy[0] * dxy[1] * dxy[2] / (npxye[0] * npxye[1] * npxye[2]);

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->fluid_.lm_lumped);
    VectorAssign(nodec, this->solid_.lm_lumped);
    for (int m = 0; m < nelem; m++) {
        const std::array<int, 3> ijk = IndexToIJK(m, xyelem);

        std::array<double, 3> xye, xyp;
        xye[0] = xymin[0] + dxy[0] * (double(ijk[0]) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(ijk[1]) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ijk[2]) + 0.5e0);
        for (int iz = 0; iz < npxye[2]; iz++) {
            xyp[2] = xye[2] + dec2p[iz][2];
            for (int iy = 0; iy < npxye[1]; iy++) {
                xyp[1] = xye[1] + dec2p[iy][1];
                for (int ix = 0; ix < npxye[0]; ix++) {
                    xyp[0] = xye[0] + dec2p[ix][0];
                    MakeSF(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

                    double phi_f = 0.0e0, phi_s = 0.0e0;
                    std::array<double, 3> grad_phi_f{}, grad_phi_s{};
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        double dsfi1 = dsf[ni][0];
                        double dsfi2 = dsf[ni][1];
                        double dsfi3 = dsf[ni][2];
                        phi_f += sfi * this->fluid_.nphi[nid];
                        phi_s += sfi * this->solid_.nphi[nid];
                        grad_phi_f[0] += dsfi1 * this->fluid_.nphi[nid];
                        grad_phi_f[1] += dsfi2 * this->fluid_.nphi[nid];
                        grad_phi_f[2] += dsfi3 * this->fluid_.nphi[nid];
                        grad_phi_s[0] += dsfi1 * this->solid_.nphi[nid];
                        grad_phi_s[1] += dsfi2 * this->solid_.nphi[nid];
                        grad_phi_s[2] += dsfi3 * this->solid_.nphi[nid];
                    }
                    std::array<double, 3> grad_phi;
                    for (int i = 0; i < 3; i++) { grad_phi[i] = phi_s * grad_phi_f[i] - phi_f * grad_phi_s[i]; }
                    double norm_grad_phi = NormVec3(grad_phi);

                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        this->fluid_.lm_lumped[nid] += sfi * volp * norm_grad_phi;
                        this->solid_.lm_lumped[nid] += sfi * volp * norm_grad_phi;
                    }
                }
            }
        }
    }

    NodeVarComm(this->fluid_.lm_lumped, 0);
    NodeVarComm(this->solid_.lm_lumped, 0);

    return;
}

void MPMMPMBlockFSI::SolveFSISystem() {

    this->DetectFSIInterface();

    std::vector<double> fluid_velocity(nodec * 3), solid_velocity(nodec * 3), acceleration(nodec * 3);
    double rtr_fsi_ref, r0r_fsi_ref, rkr_fsi_ref, rtr_fsi_abs;
    for (int block_it = 0; block_it <= this->max_block_iter; block_it++) {

        this->fluid_.SolveNS();

        this->solid_.SolveSolid();

        this->fluid_.ComputeNodeVelAccelFromDispl(fluid_velocity, acceleration);

        this->solid_.ComputeNodeVelAccelFromDispl(solid_velocity, acceleration);

        if (block_it == 0) {
            this->CalFSIResidual(fluid_velocity, solid_velocity, r0r_fsi_ref, rtr_fsi_abs);
            rkr_fsi_ref = r0r_fsi_ref;
        } else {
            this->CalFSIResidual(fluid_velocity, solid_velocity, rkr_fsi_ref, rtr_fsi_abs);
        }
        rtr_fsi_ref = (r0r_fsi_ref > 1.0e-12) ? (rkr_fsi_ref / r0r_fsi_ref) : 0.0e0;

        // --- Convergence check ---
        if (rkr_fsi_ref < this->tol_ref * r0r_fsi_ref || //
            rtr_fsi_abs < this->tol_abs) {
            if (myrank == 0) {
                std::cout << "Block_converge:" << std::setw(10) << block_it  //
                          << std::setw(15) << std::scientific << rtr_fsi_ref //
                          << std::setw(15) << std::scientific << r0r_fsi_ref //
                          << std::setw(15) << std::scientific << rtr_fsi_abs << "\n";
            }
            break;
        } else if (block_it == this->max_block_iter) {
            if (myrank == 0) {
                std::cout << "Block_it:" << std::setw(10) << block_it << "\n";
                std::cout << "Block iteration can not converge" << "\n";
            }
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        this->UpdateFSIMultiplier(block_it, fluid_velocity, solid_velocity);
    }

    // this->ReportFSIContinuity(fluid_velocity, solid_velocity);
    // this->ReportFSIPressureProfile();

    return;
}

void MPMMPMBlockFSI::CalFSIResidual(const std::vector<double> &fluid_velocity,
                                    const std::vector<double> &solid_velocity, double &rtr_ref, double &rtr_dof) {

    double norm = 0.0e0, intf_num = 0.0e0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        int nid = this->fsi_intf.nbc[n];

        double dr1 = solid_velocity[nid + nuc] - fluid_velocity[nid + nuc];
        double dr2 = solid_velocity[nid + nvc] - fluid_velocity[nid + nvc];
        double dr3 = solid_velocity[nid + nwc] - fluid_velocity[nid + nwc];

        norm += dr1 * dr1 * dbc[nid + nuc]   //
                + dr2 * dr2 * dbc[nid + nvc] //
                + dr3 * dr3 * dbc[nid + nwc];

        intf_num += dbc[nid];
    }

    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &intf_num, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    rtr_ref = std::sqrt(norm);
    rtr_dof = std::sqrt(norm / intf_num);

    return;
}

void MPMMPMBlockFSI::ReportFSIPressureProfile() {

    const int layer_count = xynodecw[2];
    std::vector<int> active_count(layer_count, 0), interface_count(layer_count, 0);
    std::vector<double> coordinate_sum(layer_count, 0.0e0), pressure_sum(layer_count, 0.0e0);
    std::vector<double> pressure_min(layer_count, std::numeric_limits<double>::max());
    std::vector<double> pressure_max(layer_count, std::numeric_limits<double>::lowest());
    std::vector<double> fluid_phi_sum(layer_count, 0.0e0), solid_phi_sum(layer_count, 0.0e0);
    std::vector<double> gf_sum(layer_count, 0.0e0), gs_sum(layer_count, 0.0e0), lambda_z_sum(layer_count, 0.0e0);
    std::vector<unsigned char> is_interface(nodec, 0);
    for (int n = 0; n < this->fsi_intf.ibc; n++) { is_interface[this->fsi_intf.nbc[n]] = 1; }

    for (int nid = 0; nid < nodec; nid++) {
        if (this->fluid_.NS_.natural_is_owned[nid] == 0 || this->fluid_.nmass[nid] <= mtol) continue;

        const PetscInt global_node = this->fluid_.NS_.natural_block_gids[nid];
        const int layer = static_cast<int>(global_node / (xynodecw[0] * xynodecw[1]));
        active_count[layer]++;
        coordinate_sum[layer] += xyc[nid][2];
        pressure_sum[layer] += this->fluid_.npres[nid];
        pressure_min[layer] = std::min(pressure_min[layer], this->fluid_.npres[nid]);
        pressure_max[layer] = std::max(pressure_max[layer], this->fluid_.npres[nid]);
        fluid_phi_sum[layer] += this->fluid_.nphi[nid];
        solid_phi_sum[layer] += this->solid_.nphi[nid];
        if (is_interface[nid] != 0) {
            interface_count[layer]++;
            gf_sum[layer] += this->fluid_.lm_lumped[nid];
            gs_sum[layer] += this->solid_.lm_lumped[nid];
            lambda_z_sum[layer] += this->nfsi_force[nid + nwc];
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, active_count.data(), layer_count, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, interface_count.data(), layer_count, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    for (std::vector<double> *values :
         {&coordinate_sum, &pressure_sum, &fluid_phi_sum, &solid_phi_sum, &gf_sum, &gs_sum, &lambda_z_sum}) {
        MPI_Allreduce(MPI_IN_PLACE, values->data(), layer_count, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }
    MPI_Allreduce(MPI_IN_PLACE, pressure_min.data(), layer_count, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, pressure_max.data(), layer_count, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

    if (myrank == 0) {
        int first_interface_layer = layer_count;
        int last_interface_layer = -1;
        for (int layer = 0; layer < layer_count; layer++) {
            if (interface_count[layer] == 0) continue;
            first_interface_layer = std::min(first_interface_layer, layer);
            last_interface_layer = layer;
        }
        const int first_layer = std::max(0, first_interface_layer - 3);
        const int last_layer = std::min(layer_count - 1, last_interface_layer + 6);
        for (int layer = first_layer; layer <= last_layer; layer++) {
            if (active_count[layer] == 0) continue;
            const double active_scale = 1.0e0 / active_count[layer];
            const double interface_scale = interface_count[layer] > 0 ? 1.0e0 / interface_count[layer] : 0.0e0;
            std::cout << "FSI_pressure_cp:" << std::setw(8) << istep << std::setw(6) << layer << std::setw(8)
                      << active_count[layer] << std::setw(8) << interface_count[layer] << std::setw(15)
                      << std::scientific << coordinate_sum[layer] * active_scale << std::setw(15)
                      << pressure_sum[layer] * active_scale << std::setw(15) << pressure_min[layer] << std::setw(15)
                      << pressure_max[layer] << std::setw(15) << fluid_phi_sum[layer] * active_scale << std::setw(15)
                      << solid_phi_sum[layer] * active_scale << std::setw(15) << gf_sum[layer] * interface_scale
                      << std::setw(15) << gs_sum[layer] * interface_scale << std::setw(15)
                      << lambda_z_sum[layer] * interface_scale << "\n";
        }
    }

    return;
}

void MPMMPMBlockFSI::ReportFSIContinuity(const std::vector<double> &fluid_velocity,
                                         const std::vector<double> &solid_velocity) {

    const int component_offsets[3] = {nuc, nvc, nwc};
    double displacement_norm_sq = 0.0e0;
    double velocity_norm_sq = 0.0e0;
    double history_velocity_norm_sq = 0.0e0;
    double history_acceleration_norm_sq = 0.0e0;
    double displacement_max = 0.0e0;
    double velocity_max = 0.0e0;
    double history_velocity_max = 0.0e0;
    double history_acceleration_max = 0.0e0;
    double velocity_identity_error_max = 0.0e0;
    int interface_count = 0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (this->fluid_.NS_.natural_is_owned[nid] == 0) continue;

        for (int offset : component_offsets) {
            const double displacement_jump = this->solid_.ndispl[nid + offset] - this->fluid_.ndispl[nid + offset];
            const double velocity_jump = solid_velocity[nid + offset] - fluid_velocity[nid + offset];
            const double history_velocity_jump = this->solid_.nvel[nid + offset] - this->fluid_.nvel[nid + offset];
            const double history_acceleration_jump =
                this->solid_.naccel[nid + offset] - this->fluid_.naccel[nid + offset];
            const double reconstructed_velocity_jump = this->solid_.nb_para[0] * this->solid_.ndispl[nid + offset] -
                                                       this->fluid_.nb_para[0] * this->fluid_.ndispl[nid + offset] -
                                                       this->solid_.nb_para[1] * this->solid_.nvel[nid + offset] +
                                                       this->fluid_.nb_para[1] * this->fluid_.nvel[nid + offset] -
                                                       this->solid_.nb_para[2] * this->solid_.naccel[nid + offset] +
                                                       this->fluid_.nb_para[2] * this->fluid_.naccel[nid + offset];
            displacement_norm_sq += displacement_jump * displacement_jump;
            velocity_norm_sq += velocity_jump * velocity_jump;
            history_velocity_norm_sq += history_velocity_jump * history_velocity_jump;
            history_acceleration_norm_sq += history_acceleration_jump * history_acceleration_jump;
            displacement_max = std::max(displacement_max, std::abs(displacement_jump));
            velocity_max = std::max(velocity_max, std::abs(velocity_jump));
            history_velocity_max = std::max(history_velocity_max, std::abs(history_velocity_jump));
            history_acceleration_max = std::max(history_acceleration_max, std::abs(history_acceleration_jump));
            velocity_identity_error_max =
                std::max(velocity_identity_error_max, std::abs(velocity_jump - reconstructed_velocity_jump));
        }
        interface_count++;
    }

    MPI_Allreduce(MPI_IN_PLACE, &displacement_norm_sq, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &velocity_norm_sq, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &history_velocity_norm_sq, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &history_acceleration_norm_sq, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &displacement_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &velocity_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &history_velocity_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &history_acceleration_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &velocity_identity_error_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &interface_count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

    if (myrank == 0 && interface_count > 0) {
        const double component_count = 3.0e0 * interface_count;
        std::cout << "FSI_continuity:" << std::setw(8) << istep << std::setw(8) << interface_count << std::setw(15)
                  << std::scientific << std::sqrt(displacement_norm_sq / component_count) << std::setw(15)
                  << displacement_max << std::setw(15) << std::sqrt(velocity_norm_sq / component_count) << std::setw(15)
                  << velocity_max << "\n";
        std::cout << "FSI_history:" << std::setw(11) << istep << std::setw(8) << interface_count << std::setw(15)
                  << std::sqrt(history_velocity_norm_sq / component_count) << std::setw(15) << history_velocity_max
                  << std::setw(15) << std::sqrt(history_acceleration_norm_sq / component_count) << std::setw(15)
                  << history_acceleration_max << std::setw(15) << velocity_identity_error_max << "\n";
    }

    return;
}

void MPMMPMBlockFSI::BuildLumpedSchurPreconditioner(const std::vector<PetscInt> &local_interface_ids,
                                                    PetscInt local_dofs, PetscInt global_dofs,
                                                    Mat &preconditioner_mat) {

    MatCreateAIJ(PETSC_COMM_WORLD, local_dofs, local_dofs, global_dofs, global_dofs, 3, nullptr, 3, nullptr,
                 &preconditioner_mat);
    MatSetBlockSize(preconditioner_mat, 3);
    MatSetOption(preconditioner_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    std::vector<PetscInt> fluid_indices, solid_indices;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (this->fluid_.NS_.natural_is_owned[nid] == 0) continue;

        for (int var = 0; var < this->fluid_.NS_.ndof; var++) {
            fluid_indices.push_back(this->fluid_.NS_.natural_var_gids[nid + var * nodec]);
        }
        for (int var = 0; var < this->solid_.SM_.ndof; var++) {
            solid_indices.push_back(this->solid_.SM_.natural_var_gids[nid + var * nodec]);
        }
    }
    // Repartition the interface blocks to their physical-node owners.
    IS fluid_is, solid_is;
    ISCreateGeneral(PETSC_COMM_WORLD, fluid_indices.size(), fluid_indices.data(), PETSC_COPY_VALUES, &fluid_is);
    ISCreateGeneral(PETSC_COMM_WORLD, solid_indices.size(), solid_indices.data(), PETSC_COPY_VALUES, &solid_is);
    Mat fluid_interface, solid_interface;
    MatCreateSubMatrix(this->fluid_.NS_.petsc_mat, fluid_is, fluid_is, MAT_INITIAL_MATRIX, &fluid_interface);
    MatCreateSubMatrix(this->solid_.SM_.petsc_mat, solid_is, solid_is, MAT_INITIAL_MATRIX, &solid_interface);
    PetscInt fluid_start, solid_start;
    MatGetOwnershipRange(fluid_interface, &fluid_start, nullptr);
    MatGetOwnershipRange(solid_interface, &solid_start, nullptr);
    int local_index = 0;
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (this->fluid_.NS_.natural_is_owned[nid] == 0) continue;

        std::array<PetscInt, 4> fluid_rows;
        std::array<PetscInt, 3> solid_rows;
        std::array<PetscScalar, 16> fluid_block;
        std::array<PetscScalar, 9> solid_block;
        for (int var = 0; var < 4; var++) { fluid_rows[var] = fluid_start + 4 * local_index + var; }
        for (int var = 0; var < 3; var++) { solid_rows[var] = solid_start + 3 * local_index + var; }
        MatGetValues(fluid_interface, 4, fluid_rows.data(), 4, fluid_rows.data(), fluid_block.data());
        MatGetValues(solid_interface, 3, solid_rows.data(), 3, solid_rows.data(), solid_block.data());

        std::array<std::array<double, 3>, 3> fluid_velocity_block, solid_velocity_block;
        std::array<bool, 3> fluid_bc, solid_bc;
        for (int i = 0; i < 3; i++) {
            fluid_bc[i] =
                std::binary_search(this->fluid_.NS_.petsc_bc_gids.begin(), this->fluid_.NS_.petsc_bc_gids.end(),
                                   this->fluid_.NS_.natural_var_gids[nid + i * nodec]);
            solid_bc[i] =
                std::binary_search(this->solid_.SM_.petsc_bc_gids.begin(), this->solid_.SM_.petsc_bc_gids.end(),
                                   this->solid_.SM_.natural_var_gids[nid + i * nodec]);
            for (int j = 0; j < 3; j++) {
                fluid_velocity_block[i][j] = PetscRealPart(
                    fluid_block[4 * i + j] - fluid_block[4 * i + 3] * fluid_block[12 + j] / fluid_block[15]);
                solid_velocity_block[i][j] = PetscRealPart(solid_block[3 * i + j]);
            }
        }
        const auto fluid_inverse = InvMat3(fluid_velocity_block);
        const auto solid_inverse = InvMat3(solid_velocity_block);
        const double gf = this->fluid_.lm_lumped[nid];
        const double gs = this->solid_.lm_lumped[nid];
        const double fluid_factor = this->fluid_.nb_para[0] * gf * gf;
        const double solid_factor = this->solid_.nb_para[0] * gs * gs;
        std::array<PetscScalar, 9> schur_block;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                schur_block[3 * i + j] = (fluid_bc[i] || fluid_bc[j] ? 0.0e0 : fluid_factor * fluid_inverse[i][j]) +
                                         (solid_bc[i] || solid_bc[j] ? 0.0e0 : solid_factor * solid_inverse[i][j]);
            }
            // Fully prescribed components have no multiplier response.
            if (fluid_bc[i] && solid_bc[i]) { schur_block[3 * i + i] = 1.0e0; }
        }
        const PetscInt row = 3 * local_interface_ids[n];
        const std::array<PetscInt, 3> rows{row, row + 1, row + 2};
        MatSetValues(preconditioner_mat, 3, rows.data(), 3, rows.data(), schur_block.data(), INSERT_VALUES);
        local_index++;
    }

    MatAssemblyBegin(preconditioner_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(preconditioner_mat, MAT_FINAL_ASSEMBLY);

    MatDestroy(&fluid_interface);
    MatDestroy(&solid_interface);
    ISDestroy(&fluid_is);
    ISDestroy(&solid_is);

    return;
}

void MPMMPMBlockFSI::UpdateFSIMultiplier(int block_it, const std::vector<double> &fluid_velocity,
                                         const std::vector<double> &solid_velocity) {

    std::vector<PetscInt> local_global_nodes(this->fsi_intf.ibc);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        local_global_nodes[n] = this->fluid_.NS_.natural_block_gids[this->fsi_intf.nbc[n]];
    }

    std::vector<int> recv_counts(nprocs), displacements(nprocs);
    const int local_count = static_cast<int>(local_global_nodes.size());
    MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, PETSC_COMM_WORLD);
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);
    const int gathered_count = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);

    std::vector<PetscInt> interface_nodes(gathered_count);
    MPI_Allgatherv(local_global_nodes.data(), local_count, MPIU_INT, interface_nodes.data(), recv_counts.data(),
                   displacements.data(), MPIU_INT, PETSC_COMM_WORLD);
    std::sort(interface_nodes.begin(), interface_nodes.end());
    interface_nodes.erase(std::unique(interface_nodes.begin(), interface_nodes.end()), interface_nodes.end());

    const PetscInt interface_count = static_cast<PetscInt>(interface_nodes.size());
    const PetscInt global_dofs = 3 * interface_count;
    std::vector<PetscInt> local_interface_ids(this->fsi_intf.ibc);
    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        local_interface_ids[n] = static_cast<PetscInt>(
            std::lower_bound(interface_nodes.begin(), interface_nodes.end(), local_global_nodes[n]) -
            interface_nodes.begin());
    }

    const PetscInt local_interface_count = interface_count / nprocs + (myrank < interface_count % nprocs ? 1 : 0);
    const PetscInt local_dofs = 3 * local_interface_count;
    Vec schur_rhs = nullptr, delta_multiplier = nullptr;
    VecCreateMPI(PETSC_COMM_WORLD, local_dofs, global_dofs, &schur_rhs);
    VecDuplicate(schur_rhs, &delta_multiplier);
    VecZeroEntries(schur_rhs);
    VecZeroEntries(delta_multiplier);

    for (int n = 0; n < this->fsi_intf.ibc; n++) {
        const int nid = this->fsi_intf.nbc[n];
        if (this->fluid_.NS_.natural_is_owned[nid] == 0) continue;

        const PetscInt row = 3 * local_interface_ids[n];
        const std::array<PetscInt, 3> rows{row, row + 1, row + 2};
        std::array<PetscScalar, 3> values;
        const double gf = this->fluid_.lm_lumped[nid];
        const double gs = this->solid_.lm_lumped[nid];
        for (int var = 0; var < 3; var++) {
            const int id = nid + var * nodec;
            values[var] = gf * fluid_velocity[id] - gs * solid_velocity[id];
        }
        VecSetValues(schur_rhs, 3, rows.data(), values.data(), INSERT_VALUES);
    }
    VecAssemblyBegin(schur_rhs);
    VecAssemblyEnd(schur_rhs);

    this->schur_matvec_count_ = 0;
    this->schur_fluid_iterations_ = 0;
    this->schur_solid_iterations_ = 0;
    this->BuildCouplingOperator(this->fluid_.NS_, this->fluid_.lm_lumped, local_interface_ids, local_dofs, global_dofs,
                                this->schur_fluid_coupling_, this->schur_fluid_rhs_, this->schur_fluid_solution_);
    this->BuildCouplingOperator(this->solid_.SM_, this->solid_.lm_lumped, local_interface_ids, local_dofs, global_dofs,
                                this->schur_solid_coupling_, this->schur_solid_rhs_, this->schur_solid_solution_);
    VecDuplicate(schur_rhs, &this->schur_solid_response_);

    Mat preconditioner_mat = nullptr;
    this->BuildLumpedSchurPreconditioner(local_interface_ids, local_dofs, global_dofs, preconditioner_mat);

    Mat schur_mat = nullptr;
    MatCreateShell(PETSC_COMM_WORLD, local_dofs, local_dofs, global_dofs, global_dofs, this, &schur_mat);
    PetscErrorCode (*apply_schur)(Mat, Vec, Vec) = [](Mat mat, Vec trial, Vec response) -> PetscErrorCode {
        void *owner;
        MatShellGetContext(mat, &owner);
        return static_cast<MPMMPMBlockFSI *>(owner)->ApplyExactSchur(trial, response);
    };
    MatShellSetOperation(schur_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(apply_schur));

    this->BuildSolidResponse();
    this->BuildFluidResponse(block_it);

    KSP schur_ksp = nullptr;
    KSPCreate(PETSC_COMM_WORLD, &schur_ksp);
    KSPSetOptionsPrefix(schur_ksp, "fsi_schur_");
    KSPSetOperators(schur_ksp, schur_mat, preconditioner_mat);
    KSPSetType(schur_ksp, KSPFGMRES);
    KSPGMRESSetRestart(schur_ksp, 50);
    KSPSetPCSide(schur_ksp, PC_RIGHT);
    KSPSetNormType(schur_ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(schur_ksp, 1.0e-8, 1.0e-20, 1.0e6, 200);
    PC schur_pc = nullptr;
    KSPGetPC(schur_ksp, &schur_pc);
    PCSetType(schur_pc, PCPBJACOBI);
    KSPSetFromOptions(schur_ksp);

    PetscReal rhs_norm = 0.0e0;
    VecNorm(schur_rhs, NORM_2, &rhs_norm);
    KSPSolve(schur_ksp, schur_rhs, delta_multiplier);

    KSPConvergedReason reason;
    PetscInt outer_iterations;
    PetscReal outer_residual, delta_norm;
    KSPGetConvergedReason(schur_ksp, &reason);
    KSPGetIterationNumber(schur_ksp, &outer_iterations);
    KSPGetResidualNorm(schur_ksp, &outer_residual);
    VecNorm(delta_multiplier, NORM_2, &delta_norm);

    if (reason > 0) {
        std::vector<PetscInt> multiplier_rows(3 * this->fsi_intf.ibc);
        for (int n = 0; n < this->fsi_intf.ibc; n++) {
            for (int var = 0; var < 3; var++) { multiplier_rows[3 * n + var] = 3 * local_interface_ids[n] + var; }
        }
        IS from_is = nullptr, to_is = nullptr;
        Vec local_multiplier = nullptr;
        VecScatter multiplier_scatter = nullptr;
        ISCreateGeneral(PETSC_COMM_SELF, multiplier_rows.size(), multiplier_rows.data(), PETSC_COPY_VALUES, &from_is);
        ISCreateStride(PETSC_COMM_SELF, multiplier_rows.size(), 0, 1, &to_is);
        VecCreateSeq(PETSC_COMM_SELF, multiplier_rows.size(), &local_multiplier);
        VecScatterCreate(delta_multiplier, from_is, local_multiplier, to_is, &multiplier_scatter);
        VecScatterBegin(multiplier_scatter, delta_multiplier, local_multiplier, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(multiplier_scatter, delta_multiplier, local_multiplier, INSERT_VALUES, SCATTER_FORWARD);
        const PetscScalar *delta_values = nullptr;
        VecGetArrayRead(local_multiplier, &delta_values);
        for (int n = 0; n < this->fsi_intf.ibc; n++) {
            const int nid = this->fsi_intf.nbc[n];
            this->nfsi_force[nid + nuc] += PetscRealPart(delta_values[3 * n]);
            this->nfsi_force[nid + nvc] += PetscRealPart(delta_values[3 * n + 1]);
            this->nfsi_force[nid + nwc] += PetscRealPart(delta_values[3 * n + 2]);
        }
        VecRestoreArrayRead(local_multiplier, &delta_values);
        VecScatterDestroy(&multiplier_scatter);
        VecDestroy(&local_multiplier);
        ISDestroy(&from_is);
        ISDestroy(&to_is);
    }

    if (myrank == 0) {
        std::cout << "Schur_FGMRES_blockPC:" << std::setw(8) << istep << std::setw(6) << block_it << std::setw(8)
                  << interface_count << std::setw(8) << outer_iterations << std::setw(8) << this->schur_matvec_count_
                  << std::setw(12) << this->schur_fluid_iterations_ << std::setw(12) << this->schur_solid_iterations_
                  << std::setw(15) << std::scientific << rhs_norm << std::setw(15) << outer_residual << std::setw(15)
                  << delta_norm << std::setw(6) << static_cast<int>(reason) << "\n";
    }

    KSPDestroy(&schur_ksp);
    KSPDestroy(&this->schur_solid_ksp_);
    MatDestroy(&preconditioner_mat);
    MatDestroy(&schur_mat);
    MatDestroy(&this->schur_fluid_coupling_);
    MatDestroy(&this->schur_solid_coupling_);
    VecDestroy(&this->schur_fluid_rhs_);
    VecDestroy(&this->schur_fluid_solution_);
    VecDestroy(&this->schur_solid_rhs_);
    VecDestroy(&this->schur_solid_solution_);
    VecDestroy(&this->schur_solid_response_);
    VecDestroy(&delta_multiplier);
    VecDestroy(&schur_rhs);

    return;
}
