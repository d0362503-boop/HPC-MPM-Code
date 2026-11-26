#include "module/solver/native_solver.h"
#include "module/bc.h"
#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solver/crsmat.h"
#include <cmath>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

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
        btb2 += b * b * mat.grid_.part.dbc[n];
        rtr2 += r * r * mat.grid_.part.dbc[n];
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
            rAp2 += rs[n] * Ap[n] * mat.grid_.part.dbc[n];
        }
        double rAp;
        MPI_Allreduce(&rAp2, &rAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double rsr0 = rsr;
        double alph = rsr0 / rAp;

        double buf5_local[5] = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};
        for (int n = 0; n < var_size; n++) {
            buf5_local[0] += Az[n] * Az[n] * mat.grid_.part.dbc[n];
            buf5_local[1] += Ar[n] * rr[n] * mat.grid_.part.dbc[n];
            buf5_local[2] += Az[n] * rr[n] * mat.grid_.part.dbc[n];
            buf5_local[3] += Ar[n] * Az[n] * mat.grid_.part.dbc[n];
            buf5_local[4] += Ar[n] * Ar[n] * mat.grid_.part.dbc[n];
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
            rsr2 += rs[n] * rr[n] * mat.grid_.part.dbc[n];
            rtr2 += rr[n] * rr[n] * mat.grid_.part.dbc[n];
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
        btb2 += b * b * mat.grid_.part.dbc[n];
        rtr2 += r * r * mat.grid_.part.dbc[n];
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
            rAp2 += rs[n] * Ap[n] * mat.grid_.part.dbc[n];
        }
        double rAp;
        MPI_Allreduce(&rAp2, &rAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double rsr0 = rsr;
        double alph = rsr0 / rAp;

        double buf5_local[5] = {0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0};
        for (int n = 0; n < var_size; n++) {
            buf5_local[0] += Az[n] * Az[n] * mat.grid_.part.dbc[n];
            buf5_local[1] += Ar[n] * rr[n] * mat.grid_.part.dbc[n];
            buf5_local[2] += Az[n] * rr[n] * mat.grid_.part.dbc[n];
            buf5_local[3] += Ar[n] * Az[n] * mat.grid_.part.dbc[n];
            buf5_local[4] += Ar[n] * Ar[n] * mat.grid_.part.dbc[n];
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
            rsr2 += rs[n] * rr[n] * mat.grid_.part.dbc[n];
            rtr2 += rr[n] * rr[n] * mat.grid_.part.dbc[n];
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
        btb2 += b * b * mat.grid_.part.dbc[n];
        rtr2 += r * r * mat.grid_.part.dbc[n];
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
        for (int n = 0; n < var_size; n++) { rAp2 += rs[n] * Ap[n] * mat.grid_.part.dbc[n]; }
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
            buf5_local[0] += yy[n] * yy[n] * mat.grid_.part.dbc[n];
            buf5_local[1] += At[n] * tt[n] * mat.grid_.part.dbc[n];
            buf5_local[2] += At[n] * At[n] * mat.grid_.part.dbc[n];
            buf5_local[3] += yy[n] * tt[n] * mat.grid_.part.dbc[n];
            buf5_local[4] += yy[n] * At[n] * mat.grid_.part.dbc[n];
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
            rsr2 += rs[n] * rr[n] * mat.grid_.part.dbc[n];
            rtr2 += rr[n] * rr[n] * mat.grid_.part.dbc[n];
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

std::vector<double> CrsMat::MatVecMult(const std::vector<double> &xx) {
    const ControlMesh &cmesh = this->grid_.cmesh;

    std::vector<double> rr(cmesh.nodec * this->ndof, 0.0e0);

    for (int i = 0; i < cmesh.nodec; i++) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; j++) {
            const int col = this->matcolid[j];
            for (int b = 0; b < this->num_block; b++) {
                rr[i + cmesh.nodec * this->block_row[b]] += this->amat[j + this->block_id[b]] * xx[col + cmesh.nodec * this->block_col[b]];
            }
        }
    }

    std::vector<int> offsets(this->ndof);
    for (int n = 0; n < this->ndof; n++) { offsets[n] = cmesh.nodec * n; }
    NodeVarComm(this->grid_.part, rr, offsets);

    this->owner_->BCResidualSet(rr);

    return rr;
}

void CrsMat::ComputeDiagonalInverseSqrt() {
    const ControlMesh &cmesh = this->grid_.cmesh;
    const RunParams &params = this->params_;

    int var_size = cmesh.nodec * this->ndof;
    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->adiag[n] = 1.0e0 / std::sqrt(std::abs(this->adiag[n]));
        } else {
            this->adiag[n] = (std::abs(this->adiag[n]) < params.mtol) ? 0.0e0 : (1.0e0 / std::sqrt(std::abs(this->adiag[n])));
        }
    }

    return;
}

void CrsMat::ApplyDiagonalScaling() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    for (int i = 0; i < cmesh.nodec; i++) {
        for (int j = this->matrow[i]; j < this->matrow[i + 1]; j++) {
            const int col = this->matcolid[j];

            for (int b = 0; b < this->num_block; b++) {
                const int row_offset = cmesh.nodec * this->block_row[b];
                const int col_offset = cmesh.nodec * this->block_col[b];
                this->amat[j + this->block_id[b]] *= this->adiag[i + row_offset] * this->adiag[col + col_offset];
            }
        }
    }

    return;
}

void CrsMat::BuildDiagonalPreconditioner() {
    const ControlMesh &cmesh = this->grid_.cmesh;

    VectorAssign(cmesh.nodec * this->ndof, this->adiag);
    this->ExtractDiagonal();
    this->ComputeDiagonalInverseSqrt();
    this->ApplyDiagonalScaling();

    return;
}

void CrsMat::ScaleSolution(std::vector<double> &rr) {
    const RunParams &params = this->params_;

    const int var_size = this->x_lhs.size();

    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->x_lhs[n] /= this->adiag[n];
            rr[n] = this->b_rhs[n] * this->adiag[n];
        } else {
            if (this->adiag[n] < params.mtol) {
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

void CrsMat::RestorSolution() {
    const RunParams &params = this->params_;

    const int var_size = this->x_lhs.size();

    for (int n = 0; n < var_size; n++) {
        if (this->FEM_flag) {
            this->x_lhs[n] *= this->adiag[n];
        } else {
            if (this->adiag[n] > params.mtol) {
                this->x_lhs[n] *= this->adiag[n];
            } else {
                this->x_lhs[n] = 0.0e0;
            }
        }
    }

    return;
}
