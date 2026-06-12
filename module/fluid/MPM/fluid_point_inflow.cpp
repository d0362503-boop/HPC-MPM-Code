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

    int num = (this->uinfbc.ibc + this->vinfbc.ibc + this->winfbc.ibc) * 128;

    this->ifp.num = 0;

    VectorAssign(num, this->ifp.id);
    VectorAssign(num, this->ifp.matid);
    VectorAssign(num, this->ifp.mass);
    VectorAssign(num, this->ifp.vol);
    VectorAssign(num, this->ifp.pres);
    VectorAssign(num, this->ifp.coord);
    VectorAssign(num, this->ifp.vel);
    VectorAssign(num, this->ifp.accel);

    if (this->solswitch == MapScheme::TPIC) {
        VectorAssign(num, this->ifp.tpic.vel_grad);
        VectorAssign(num, this->ifp.tpic.accel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        VectorAssign(num, this->ifp.apic.vel_Bmat);
        VectorAssign(num, this->ifp.apic.accel_Bmat);
        VectorAssign(num, this->ifp.apic.inv_Dmat);
    }

    if (this->uinfbc.ibc != 0) { GenerateInflowParticles(0, this->ifp, this->uinfbc); }
    if (this->vinfbc.ibc != 0) { GenerateInflowParticles(1, this->ifp, this->vinfbc); }
    if (this->winfbc.ibc != 0) { GenerateInflowParticles(2, this->ifp, this->winfbc); }

    return;
}

void StabilizedMPM::GenerateInflowParticles(int dir, MaterialPoint &ifp, //
                                            const BoundaryCondition &ifbc) {

    constexpr int offset_layer = 1;

    constexpr std::array<double, 2> k1{0.99e0, 0.0e0};
    constexpr std::array<double, 2> k2{0.99e0, 0.875e0};
    std::array<double, 2> offset_coeff = (offset_layer == 2) ? k2 : k1;
    // ----

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    int nex = xyelem[0];
    int ney = xyelem[1];
    int nez = xyelem[2];

    for (int m = 0; m < ifbc.ibc; m++) {
        int ie = ifbc.nbc[m];
        int pid = this->idepf[ie];

        if (pid == -1) continue;

        int ize = ie / (nex * ney);
        int iye = (ie - ize * (nex * ney)) / nex;
        int ixe = ie - ize * (nex * ney) - iye * nex;

        std::vector<double> xyr(3);
        std::vector<int> iexp(3), ie_indices = {ixe, iye, ize};

        while (pid != -1) {
            for (int i = 0; i < 3; i++) {
                xyr[i] = this->coord[pid][i] - xymin[i];
                iexp[i] = floor(xyr[i] / dxy[i]);
            }
            int ne = nex * ney * iexp[2] + nex * iexp[1] + iexp[0];
            if (ne != ie && iexp[dir] > ie_indices[dir]) {
                for (int p = 0; p < offset_layer; p++) {
                    ifp.coord[ifp.num] = this->coord[pid];
                    ifp.coord[ifp.num][dir] -= dxy[dir] * offset_coeff[p];

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
                    ifp.id[ifp.num] = this->id[pid];
                    ifp.matid[ifp.num] = this->matid[pid];
                    ifp.mass[ifp.num] = this->mass[pid];
                    ifp.vol[ifp.num] = this->vol[pid];
                    ifp.num++;
                }
            }
            pid = this->idp2p[pid];
        }
    }

    if (ifp.num < 0) ifp.num = 0;

    return;
}