#include "solver.h"
#include "../bc.h"
#include "../dataset.h"
#include "../material_point.h"
#include "../mesh.h"
#include "../mpi_data.h"
#include "crsmat.h"
#include <cmath>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Solve the scaled linear system with the GPBiCGSafe iteration.
 * @param mat Sparse matrix and solver state.
 * @return Number of Krylov iterations performed.
 */
int GPBiCGSafe(CrsMat &mat) {

    const int var_size = mat.x_lhs.size();

    std::vector<double> rr(var_size, 0.0e0);
    std::vector<double> pp(var_size, 0.0e0);
    std::vector<double> tt(var_size, 0.0e0);
    std::vector<double> Ap(var_size, 0.0e0);
    std::vector<double> Ar(var_size, 0.0e0);
    std::vector<double> Az(var_size, 0.0e0);
    std::vector<double> Au(var_size, 0.0e0);
    std::vector<double> Ax(var_size, 0.0e0);
    std::vector<double> rs(var_size, 0.0e0);
    std::vector<double> yy(var_size, 0.0e0);
    std::vector<double> uu(var_size, 0.0e0);
    std::vector<double> zz(var_size, 0.0e0);
    std::vector<double> ww(var_size, 0.0e0);
    std::vector<double> t1(var_size, 0.0e0);

    int iter;
    const int kmax = 10000;

    mat.ScaleSolution(rr);

    if (mat.owner_) { mat.owner_->BCResidualSet(rr); }

    Ax = mat.MatVecMult(mat.x_lhs);

    double btb2 = 0.0e0, rtr2 = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        double b = rr[n];
        double r = b - Ax[n];
        rr[n] = r;
        rs[n] = r;
        btb2 += b * b * dbc[n];
        rtr2 += r * r * dbc[n];
    }

    double buf2_local[2] = {rtr2, btb2};
    double buf2_global[2];
    MPI_Allreduce(buf2_local, buf2_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rtr = buf2_global[0], btb = buf2_global[1];

    double rsr = rtr;
    btb = rtr;
    double epsbtb = 1.0e-8 * btb;
    if (btb < 1.0e-18) int kmax = 0;

    double beta = 0.0e0;
    for (iter = 0; iter < kmax; iter++) {

        Ar = mat.MatVecMult(rr);

        double rAp2 = 0.0e0;
        for (int n = 0; n < var_size; ++n) {
            pp[n] = rr[n] + beta * (pp[n] - uu[n]);
            Ap[n] = Ar[n] + beta * (Ap[n] - Au[n]);
            rAp2 += rs[n] * Ap[n] * dbc[n];
        }
        double rAp;
        MPI_Allreduce(&rAp2, &rAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double rsr0 = rsr;
        double alph = rsr0 / rAp;

        double buf5_local[5] = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};
        for (int n = 0; n < var_size; n++) {
            buf5_local[0] += Az[n] * Az[n] * dbc[n];
            buf5_local[1] += Ar[n] * rr[n] * dbc[n];
            buf5_local[2] += Az[n] * rr[n] * dbc[n];
            buf5_local[3] += Ar[n] * Az[n] * dbc[n];
            buf5_local[4] += Ar[n] * Ar[n] * dbc[n];
        }
        double buf5_global[5];
        MPI_Allreduce(buf5_local, buf5_global, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double AzAz = buf5_global[0], Arrr = buf5_global[1], Azrr = buf5_global[2];
        double ArAz = buf5_global[3], ArAr = buf5_global[4];

        double zeta, eta;
        if (iter == 0) {
            zeta = Arrr / ArAr;
            eta = 0.0e0;
        } else {
            double denom = ArAr * AzAz - ArAz * ArAz;
            zeta = (AzAz * Arrr - Azrr * ArAz) / denom;
            eta = (ArAr * Azrr - ArAz * Arrr) / denom;
        }

        for (int n = 0; n < var_size; n++) { uu[n] = zeta * Ap[n] + eta * (Az[n] + beta * uu[n]); }

        Au = mat.MatVecMult(uu);

        double rsr2 = 0.0e0, rtr2 = 0.0e0;
        for (int n = 0; n < var_size; n++) {
            tt[n] = rr[n] - alph * Ap[n];
            zz[n] = zeta * rr[n] + eta * zz[n] - alph * uu[n];
            Az[n] = zeta * Ar[n] + eta * Az[n] - alph * Au[n];
            mat.x_lhs[n] += alph * pp[n] + zz[n];
            rr[n] = tt[n] - Az[n];
            rsr2 += rs[n] * rr[n] * dbc[n];
            rtr2 += rr[n] * rr[n] * dbc[n];
        }

        double buf_r[2] = {rsr2, rtr2};
        double buf_r_g[2];
        MPI_Allreduce(buf_r, buf_r_g, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        rsr = buf_r_g[0];
        rtr = buf_r_g[1];

        if (rtr < epsbtb) break;

        beta = (alph / zeta) * (rsr / rsr0);
    }

    if (iter == kmax) {
        std::cout << "Solver can not converge" << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    mat.RestorSolution();

    return iter;
}

/**
 * @brief Solve the scaled linear system with the GPBiCGAR iteration.
 * @param mat Sparse matrix and solver state.
 * @return Number of Krylov iterations performed.
 */
int GPBiCGAR(CrsMat &mat) {

    const int var_size = mat.x_lhs.size();

    std::vector<double> rr(var_size, 0.0e0);
    std::vector<double> pp(var_size, 0.0e0);
    std::vector<double> tt(var_size, 0.0e0);
    std::vector<double> Ap(var_size, 0.0e0);
    std::vector<double> Ar(var_size, 0.0e0);
    std::vector<double> Az(var_size, 0.0e0);
    std::vector<double> Au(var_size, 0.0e0);
    std::vector<double> Ax(var_size, 0.0e0);
    std::vector<double> rs(var_size, 0.0e0);
    std::vector<double> yy(var_size, 0.0e0);
    std::vector<double> uu(var_size, 0.0e0);
    std::vector<double> zz(var_size, 0.0e0);
    std::vector<double> ww(var_size, 0.0e0);
    std::vector<double> t1(var_size, 0.0e0);

    int iter;
    const int kmax = 10000;

    mat.ScaleSolution(rr);

    if (mat.owner_) { mat.owner_->BCResidualSet(rr); }

    Ax = mat.MatVecMult(mat.x_lhs);

    double btb2 = 0.0e0, rtr2 = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        double b = rr[n];
        double r = b - Ax[n];
        rr[n] = r;
        rs[n] = r;
        btb2 += b * b * dbc[n];
        rtr2 += r * r * dbc[n];
    }

    double buf2_local[2] = {rtr2, btb2};
    double buf2_global[2];
    MPI_Allreduce(buf2_local, buf2_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double rtr = buf2_global[0], btb = buf2_global[1];

    double rsr = rtr;
    btb = rtr;
    double epsbtb = 1.0e-8 * btb;
    if (btb < 1.0e-18) int kmax = 0;

    double beta = 0.0e0;
    for (iter = 0; iter < kmax; iter++) {

        Ar = mat.MatVecMult(rr);

        double rAp2 = 0.0e0;
        for (int n = 0; n < var_size; ++n) {
            pp[n] = rr[n] + beta * (pp[n] - uu[n]);
            Ap[n] = Ar[n] + beta * (Ap[n] - Au[n]);
            rAp2 += rs[n] * Ap[n] * dbc[n];
        }
        double rAp;
        MPI_Allreduce(&rAp2, &rAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double rsr0 = rsr;
        double alph = rsr0 / rAp;

        double buf5_local[5] = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};
        for (int n = 0; n < var_size; n++) {
            buf5_local[0] += Az[n] * Az[n] * dbc[n];
            buf5_local[1] += Ar[n] * rr[n] * dbc[n];
            buf5_local[2] += Az[n] * rr[n] * dbc[n];
            buf5_local[3] += Ar[n] * Az[n] * dbc[n];
            buf5_local[4] += Ar[n] * Ar[n] * dbc[n];
        }
        double buf5_global[5];
        MPI_Allreduce(buf5_local, buf5_global, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double AzAz = buf5_global[0], Arrr = buf5_global[1], Azrr = buf5_global[2];
        double ArAz = buf5_global[3], ArAr = buf5_global[4];

        double zeta, eta;
        if (iter == 0) {
            zeta = Arrr / ArAr;
            eta = 0.0e0;
        } else {
            double denom = ArAr * AzAz - ArAz * ArAz;
            zeta = (AzAz * Arrr - Azrr * ArAz) / denom;
            eta = (ArAr * Azrr - ArAz * Arrr) / denom;
        }

        for (int n = 0; n < var_size; n++) { uu[n] = zeta * Ap[n] + eta * (tt[n] - rr[n] + beta * uu[n]); }

        Au = mat.MatVecMult(uu);

        double rsr2 = 0.0e0, rtr2 = 0.0e0;
        for (int n = 0; n < var_size; n++) {
            tt[n] = rr[n] - alph * Ap[n];
            zz[n] = zeta * rr[n] + eta * zz[n] - alph * uu[n];
            Az[n] = zeta * Ar[n] + eta * Az[n] - alph * Au[n];
            mat.x_lhs[n] += alph * pp[n] + zz[n];
            rr[n] = tt[n] - Az[n];
            rsr2 += rs[n] * rr[n] * dbc[n];
            rtr2 += rr[n] * rr[n] * dbc[n];
        }

        double buf_r[2] = {rsr2, rtr2};
        double buf_r_g[2];
        MPI_Allreduce(buf_r, buf_r_g, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        rsr = buf_r_g[0];
        rtr = buf_r_g[1];

        if (rtr < epsbtb) break;

        beta = (alph / zeta) * (rsr / rsr0);
    }

    if (iter == kmax) {
        std::cout << "Solver can not converge" << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    mat.RestorSolution();

    return iter;
}

/**
 * @brief Solve the scaled linear system with the original GPBiCG iteration.
 * @param mat Sparse matrix and solver state.
 * @return Number of Krylov iterations performed.
 */
int GPBiCG(CrsMat &mat) {

    const int var_size = mat.x_lhs.size();

    std::vector<double> rr(var_size, 0.0e0);
    std::vector<double> pp(var_size, 0.0e0);
    std::vector<double> tt(var_size, 0.0e0);
    std::vector<double> Ap(var_size, 0.0e0);
    std::vector<double> At(var_size, 0.0e0);
    std::vector<double> rs(var_size, 0.0e0);
    std::vector<double> yy(var_size, 0.0e0);
    std::vector<double> uu(var_size, 0.0e0);
    std::vector<double> zz(var_size, 0.0e0);
    std::vector<double> ww(var_size, 0.0e0);
    std::vector<double> t1(var_size, 0.0e0);

    int iter;
    const int kmax = 10000;

    // ---- Initialize xx and rr ----
    mat.ScaleSolution(rr);

    if (mat.owner_) { mat.owner_->BCResidualSet(rr); }

    Ap = mat.MatVecMult(mat.x_lhs);

    double btb2 = 0.0e0, rtr2 = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        double b = rr[n];
        double r = b - Ap[n];
        rr[n] = r;
        rs[n] = r;
        pp[n] = r;
        btb2 += b * b * dbc[n];
        rtr2 += r * r * dbc[n];
    }

    double buf2_local[2] = {btb2, rtr2};
    double buf2_global[2];
    MPI_Allreduce(buf2_local, buf2_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double btb = buf2_global[0], rtr = buf2_global[1];

    double rsr = rtr;
    btb = rtr;
    double epsbtb = 1.0e-8 * btb;

    if (btb < 1.0e-18) int kmax = 0;

    double pk = 0.0e0;
    double beta = 0.0e0;

    for (iter = 0; iter < kmax; iter++) {

        Ap = mat.MatVecMult(pp);

        double rAp2 = 0.0e0;
        for (int n = 0; n < var_size; n++) { rAp2 += rs[n] * Ap[n] * dbc[n]; }
        double rAp;
        MPI_Allreduce(&rAp2, &rAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double rsr0 = rsr;
        double alph = rsr0 / rAp;

        for (int n = 0; n < var_size; n++) {
            yy[n] = tt[n] - rr[n] - alph * (ww[n] - Ap[n]);
            uu[n] = tt[n] - rr[n] + beta * uu[n];
            tt[n] = rr[n] - alph * Ap[n];
        }

        At = mat.MatVecMult(tt);

        double buf5_local[5] = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};
        for (int n = 0; n < var_size; n++) {
            buf5_local[0] += yy[n] * yy[n] * dbc[n];
            buf5_local[1] += At[n] * tt[n] * dbc[n];
            buf5_local[2] += At[n] * At[n] * dbc[n];
            buf5_local[3] += yy[n] * tt[n] * dbc[n];
            buf5_local[4] += yy[n] * At[n] * dbc[n];
        }
        double buf5_global[5];
        MPI_Allreduce(buf5_local, buf5_global, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double yty = buf5_global[0], Att = buf5_global[1], At2 = buf5_global[2];
        double ytt = buf5_global[3], ytA = buf5_global[4];

        double zeta = (yty * Att - ytt * ytA * pk) / (At2 * yty - ytA * ytA * pk);
        double eta = ((ytt * At2 - Att * ytA) / (At2 * yty - ytA * ytA)) * pk;
        pk = 1.0e0;

        double rsr2 = 0.0e0, rtr2 = 0.0e0;
        for (int n = 0; n < var_size; n++) {
            uu[n] = zeta * Ap[n] + eta * uu[n];
            zz[n] = zeta * rr[n] + eta * zz[n] - alph * uu[n];
            mat.x_lhs[n] += alph * pp[n] + zz[n];
            rr[n] = tt[n] - zeta * At[n] - eta * yy[n];
            rsr2 += rs[n] * rr[n] * dbc[n];
            rtr2 += rr[n] * rr[n] * dbc[n];
        }

        double buf_r[2] = {rsr2, rtr2};
        double buf_r_g[2];
        MPI_Allreduce(buf_r, buf_r_g, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        rsr = buf_r_g[0];
        rtr = buf_r_g[1];

        if (rtr < epsbtb) break;

        beta = (alph / zeta) * (rsr / rsr0);
        for (int n = 0; n < var_size; n++) {
            pp[n] = rr[n] + beta * (pp[n] - uu[n]);
            ww[n] = At[n] + beta * Ap[n];
        }
    }

    if (iter == kmax) {
        std::cout << "Solver can not converge" << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    mat.RestorSolution();

    return iter;
}

/**
 * @brief Apply the local block CSR matrix and overlap accumulation to a vector.
 * @param xx Input vector in node-major block layout.
 * @return Result after matrix multiplication and overlap-node communication.
 */
std::vector<double> CrsMat::MatVecMult(const std::vector<double> &xx) {

    const int var_size = (int)this->x_lhs.size();

    std::vector<double> rr(var_size, 0.0e0);

    for (int i = 0; i < nodec; ++i) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; ++j) {
            const int col = this->matcolid[j];

            for (int m = 0; m < this->ndof; ++m) {
                const int row_off = nodec * m;
                const int block_base = m * this->ndof;

                for (int n = 0; n < this->ndof; ++n) {
                    rr[i + row_off] += this->amat[j + this->block_id[block_base + n]] * xx[col + nodec * n];
                }
            }
        }
    }

    std::vector<int> offsets(this->ndof);
    for (int n = 0; n < this->ndof; ++n) { offsets[n] = nodec * n; }
    NodeVarComm(rr, offsets);

    if (this->owner_) { this->owner_->BCResidualSet(rr); }

    return rr;
}

/**
 * @brief Compute the PETSc residual norm and the active PETSc DOF count.
 * @param res_norm Euclidean norm of the masked PETSc residual.
 * @param active_dof Number of active PETSc degrees of freedom used in the norm.
 */
void CrsMat::ComputePetscResidualStats(double &res_norm, double &active_dof) {
    Vec residual = nullptr;
    Vec active_mask = nullptr;

    VecDuplicate(this->petsc_x, &residual);
    VecDuplicate(this->petsc_x, &active_mask);

    MatMult(this->petsc_mat, this->petsc_x, residual);
    VecAYPX(residual, -1.0, this->petsc_b); // residual = petsc_b - petsc_mat * petsc_x

    VecZeroEntries(active_mask);

    const PetscInt local_size = static_cast<PetscInt>(this->local_node) * this->ndof;
    std::vector<PetscInt> indices(static_cast<size_t>(local_size));
    std::vector<PetscScalar> values(static_cast<size_t>(local_size));

    size_t idx = 0;
    for (int i = 0; i < this->local_node; ++i) {
        const int natural_id = this->interior_list[i];
        const PetscScalar is_active = (this->FEM_flag || this->active_row_mask[natural_id] != 0) ? 1.0 : 0.0;
        for (int var = 0; var < this->ndof; ++var) {
            const int local_idx = natural_id + var * nodec;
            indices[idx] = this->l2g_var_map[local_idx];
            values[idx] = is_active;
            ++idx;
        }
    }

    VecSetValues(active_mask, local_size, indices.data(), values.data(), INSERT_VALUES);
    VecAssemblyBegin(active_mask);
    VecAssemblyEnd(active_mask);

    VecPointwiseMult(residual, residual, active_mask);
    VecNorm(residual, NORM_2, &res_norm);
    VecSum(active_mask, &active_dof);

    VecDestroy(&residual);
    VecDestroy(&active_mask);
}

/**
 * @brief Compute the reference residual norm used in Newton convergence checks.
 * @return Residual norm for the current solver path.
 */
double CrsMat::ComputeRefResidual() {
    const int var_size = int(this->x_lhs.size());

    if (this->use_petsc) {
        double res_norm = 0.0e0;
        double active_dof = 0.0e0;
        this->ComputePetscResidualStats(res_norm, active_dof);
        return res_norm;
    }

    if (this->owner_) { this->owner_->BCResidualSet(this->b_rhs); }

    std::vector<double> x_pre(var_size);
    for (int n = 0; n < var_size; n++) { x_pre[n] = (this->adiag[n] > mtol) ? this->x_lhs[n] / this->adiag[n] : 0.0e0; }

    std::vector<double> Ax = this->MatVecMult(x_pre);

    double rtr = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        double Ax_phys = (this->adiag[n] > mtol) ? (Ax[n] / this->adiag[n]) : 0.0e0;
        double r = this->b_rhs[n] - Ax_phys;
        rtr += r * r * dbc[n];
    }
    MPI_Allreduce(MPI_IN_PLACE, &rtr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return std::sqrt(rtr);
}

/**
 * @brief Compute the absolute residual metric used in Newton convergence checks.
 * @return Weighted absolute residual for the current solver path.
 */
double CrsMat::ComputeAbsResidual() {
    const int var_size = int(this->x_lhs.size());

    if (this->use_petsc) {
        double res_norm = 0.0e0;
        double active_dof = 0.0e0;
        this->ComputePetscResidualStats(res_norm, active_dof);
        if (active_dof <= 1.0e-12) { return 0.0e0; }
        return res_norm / std::sqrt(active_dof);
    }

    if (this->owner_) { this->owner_->BCResidualSet(this->b_rhs); }

    std::vector<double> x_pre(var_size);
    for (int n = 0; n < var_size; n++) { x_pre[n] = (this->adiag[n] > mtol) ? this->x_lhs[n] / this->adiag[n] : 0.0e0; }

    std::vector<double> Ax = this->MatVecMult(x_pre);

    double rtr = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        double Ax_phys = (this->adiag[n] > mtol) ? (Ax[n] / this->adiag[n]) : 0.0e0;
        double r = this->b_rhs[n] - Ax_phys;
        rtr += r * r * dbc[n];
    }

    double idof = 0.0e0;
    for (int n = 0; n < var_size; n++) {
        if (this->adiag[n] > mtol) { idof += dbc[n]; }
    }

    double buf_rid[2] = {rtr, idof};
    MPI_Allreduce(MPI_IN_PLACE, buf_rid, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    rtr = buf_rid[0];
    idof = buf_rid[1];

    if (idof <= 1.0e-30) { return 0.0e0; }
    rtr = std::sqrt(rtr / idof);

    return rtr;
}

/**
 * @brief Evaluate Newton convergence and print the current residual summary.
 * @param NR_it Current Newton iteration index.
 * @param NR_it_max Maximum allowed Newton iterations.
 * @param solver_it Linear solver iteration count for this Newton step.
 * @param r0r Reference residual from the first Newton iteration.
 * @return True when the Newton process should stop.
 */
bool CrsMat::CheckNRConvergence(int NR_it, int NR_it_max, int solver_it, double &r0r) {
    double rkr;
    if (NR_it == 0) {
        r0r = this->ComputeRefResidual();
        rkr = r0r;
    } else {
        rkr = this->ComputeRefResidual();
    }
    double rtr_abs = this->ComputeAbsResidual();
    double rtr_ref = (r0r > 1.0e-30) ? (rkr / r0r) : 0.0e0;

    if (rtr_ref < 1.0e-6 || r0r < 1.0e-6 || rtr_abs < 1.0e-8) {
        if (myrank == 0) {
            std::cout << "NR_converge:" << std::setw(15) << NR_it << std::setw(15) << solver_it << std::setw(15)
                      << std::scientific << rtr_ref << std::setw(15) << std::scientific << r0r << std::setw(15)
                      << std::scientific << rtr_abs << "\n";
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
