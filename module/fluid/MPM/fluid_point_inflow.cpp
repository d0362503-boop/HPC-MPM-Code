#include "../../bc.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "stabilized_mpm.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

using namespace stabilizedmpm;

void StabilizedMPM::InflowParticles() {

    this->ifp.num = 0;

    int target_per_cell = npxye[0] * npxye[1] * npxye[2];
    int target_num = (this->uinfbc.ibc + this->vinfbc.ibc + this->winfbc.ibc) * target_per_cell;

    int local_has_inflow = (target_num != 0) ? 1 : 0;
    int global_has_inflow = 0;
    MPI_Allreduce(&local_has_inflow, &global_has_inflow, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_has_inflow == 0) return;

    if (target_num != 0) {
        VectorAssign(target_num, this->ifp.id);
        VectorAssign(target_num, this->ifp.matid);
        VectorAssign(target_num, this->ifp.mass);
        VectorAssign(target_num, this->ifp.vol);
        VectorAssign(target_num, this->ifp.pres);
        VectorAssign(target_num, this->ifp.coord);
        VectorAssign(target_num, this->ifp.vel);
        VectorAssign(target_num, this->ifp.accel);

        if (this->solswitch == MapScheme::TPIC) {
            VectorAssign(target_num, this->ifp.tpic.vel_grad);
            VectorAssign(target_num, this->ifp.tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            VectorAssign(target_num, this->ifp.apic.vel_Bmat);
            VectorAssign(target_num, this->ifp.apic.accel_Bmat);
            VectorAssign(target_num, this->ifp.apic.inv_Dmat);
        }

        if (this->uinfbc.ibc != 0) { this->GenerateInflowParticles(0, this->ifp, this->uinfbc); }
        if (this->vinfbc.ibc != 0) { this->GenerateInflowParticles(1, this->ifp, this->vinfbc); }
        if (this->winfbc.ibc != 0) { this->GenerateInflowParticles(2, this->ifp, this->winfbc); }
    }

    this->AssignUniqueInflowIds(this->ifp);

    return;
}

void StabilizedMPM::GenerateInflowParticlesEmptyMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) {

    std::array<std::array<double, 3>, 6> dec2p;
    GaussianDistribution(dec2p);

    int target_per_cell = npxye[0] * npxye[1] * npxye[2];
    double cell_vol = dxy[0] * dxy[1] * dxy[2];
    double particle_vol = cell_vol / static_cast<double>(target_per_cell);
    double particle_mass = this->rho * particle_vol;

    static int other_dir[3][2] = {{1, 2}, {0, 2}, {0, 1}};

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    for (int m = 0; m < infbc.ibc; m++) {
        int mid = infbc.nbc[m];
        int row = this->inflow_row[m];
        if (row >= npxye[dir]) continue;

        const std::array<int, 3> ie_indices = ElementIndexToIJK(mid, xyelem);
        int sign = this->GetInflowSign(dir, ie_indices[dir]);

        std::array<double, 3> xye;
        xye[0] = xymin[0] + dxy[0] * (double(ie_indices[0]) + 0.5e0);
        xye[1] = xymin[1] + dxy[1] * (double(ie_indices[1]) + 0.5e0);
        xye[2] = xymin[2] + dxy[2] * (double(ie_indices[2]) + 0.5e0);

        const int dir_local = (sign > 0) ? (npxye[dir] - 1 - row) : row;

        int j_loop = npxye[other_dir[dir][1]];
        int i_loop = npxye[other_dir[dir][0]];

        for (int j_local = 0; j_local < j_loop; ++j_local) {
            for (int i_local = 0; i_local < i_loop; ++i_local) {
                // --- Mass ---
                ifp.mass[ifp.num] = particle_mass;
                // --- Volume ---
                ifp.vol[ifp.num] = particle_vol;
                // --- Material ID ---
                ifp.matid[ifp.num] = 0;

                // --- Coordinate ---
                ifp.coord[ifp.num][dir] = xye[dir] + dec2p[dir_local][dir];
                ifp.coord[ifp.num][other_dir[dir][0]] = xye[other_dir[dir][0]] + dec2p[i_local][other_dir[dir][0]];
                ifp.coord[ifp.num][other_dir[dir][1]] = xye[other_dir[dir][1]] + dec2p[j_local][other_dir[dir][1]];
                // --- Pressure ---
                ifp.pres[ifp.num] = 0.0e0;
                // --- Velocity ---
                ifp.vel[ifp.num] = {};
                // --- Acceleration ---
                ifp.accel[ifp.num] = {};

                if (this->solswitch == MapScheme::TPIC) { // --- TPIC
                    ifp.tpic.vel_grad[ifp.num] = {};
                    ifp.tpic.accel_grad[ifp.num] = {};
                } else if (this->solswitch == MapScheme::APIC) { // --- APIC ---
                    ifp.apic.vel_Bmat[ifp.num] = {};
                    ifp.apic.accel_Bmat[ifp.num] = {};
                    ifp.apic.inv_Dmat[ifp.num] = {};
                }
                ifp.num++;
            }
        }
        this->inflow_row[m]++;
    }

    return;
}

void StabilizedMPM::GenerateInflowParticlesFilledMesh(int dir, MaterialPoint &ifp, const BoundaryCondition &infbc) {

    constexpr int offset_layer = 1;

    constexpr std::array<double, 2> k1{0.99e0, 0.0e0};
    constexpr std::array<double, 2> k2{0.99e0, 0.875e0};
    std::array<double, 2> offset_coeff = (offset_layer == 2) ? k2 : k1;

    static int other_dir[3][2] = {{1, 2}, {0, 2}, {0, 1}};

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];

    for (int m = 0; m < infbc.ibc; m++) {
        int ie = infbc.nbc[m];
        int pid = this->idepf[ie];

        if (pid == -1) continue;

        const std::array<int, 3> ie_indices = ElementIndexToIJK(ie, xyelem);

        std::array<double, 3> xyr;
        std::array<int, 3> iexp{};

        int sign = this->GetInflowSign(dir, ie_indices[dir]);

        while (pid != -1) {
            for (int i = 0; i < 3; i++) {
                xyr[i] = this->coord[pid][i] - xymin[i];
                iexp[i] = static_cast<int>(std::floor(xyr[i] / dxy[i]));
            }
            int ne = nex * ney * iexp[2] + nex * iexp[1] + iexp[0];
            if ((iexp[dir] - sign == ie_indices[dir]) &&                      //
                (iexp[other_dir[dir][0]] == ie_indices[other_dir[dir][0]]) && //
                (iexp[other_dir[dir][1]] == ie_indices[other_dir[dir][1]])) {
                for (int p = 0; p < offset_layer; p++) {
                    ifp.coord[ifp.num] = this->coord[pid];
                    ifp.coord[ifp.num][dir] -= double(sign) * dxy[dir] * offset_coeff[p];

                    std::array<double, 3> xyp = ifp.coord[ifp.num];
                    MakSf(ie, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

                    std::array<std::array<double, 3>, 3> ADp{};
                    for (int ni = 0; ni < nenode; ni++) {
                        int nid = ncm[ni];
                        double sfi = sf[ni];
                        double dsfi1 = dsf[ni][0];
                        double dsfi2 = dsf[ni][1];
                        double dsfi3 = dsf[ni][2];
                        ifp.pic.VarG2P(ifp.num, nid, sfi, ifp.pres, this->npres);
                        ifp.pic.VarG2P(ifp.num, nid, sfi, ifp.accel, this->naccel);
                        if (this->solswitch != MapScheme::FLIP) {
                            ifp.pic.VarG2P(ifp.num, nid, sfi, ifp.vel, this->nvel);
                        } else if (this->solswitch == MapScheme::FLIP) {
                            ifp.flip.VarG2P(ifp.num, nid, sfi, ifp.vel, this->naccel);
                        }
                        if (this->solswitch == MapScheme::TPIC) {
                            ifp.tpic.VarGradG2P(ifp.num, nid, ni, dsf, ifp.tpic.accel_grad, this->naccel);
                            ifp.tpic.VarGradG2P(ifp.num, nid, ni, dsf, ifp.tpic.vel_grad, this->nvel);
                        } else if (this->solswitch == MapScheme::APIC) {
                            ifp.apic.VarBmatG2P(ifp.num, nid, sfi, ifp.apic.accel_Bmat, this->naccel, ifp.coord);
                            ifp.apic.VarBmatG2P(ifp.num, nid, sfi, ifp.apic.vel_Bmat, this->nvel, ifp.coord);
                            ifp.apic.DmatG2P(ifp.num, nid, sfi, ifp.coord, ADp);
                        }
                    }
                    if (this->solswitch == MapScheme::APIC) { ifp.apic.InvDmatG2P(ifp.num, ADp); }
                    // ifp.id[ifp.num] = this->id[pid];
                    ifp.matid[ifp.num] = this->matid[pid];
                    ifp.mass[ifp.num] = this->mass[pid];
                    ifp.vol[ifp.num] = this->vol[pid];
                    ifp.num++;
                }
            }
            pid = this->idp2p[pid];
        }
    }

    return;
}
