#include "../../bc.h"
#include "../../dataset.h"
#include "../../map_and_interpolate.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../shape_function.h"
#include "stabilized_mpm.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
using namespace stabilizedmpm;

void StabilizedMPM::Particle2Node() {

    int nenode;
    vector<int> ncm;
    vector<double> sf;
    vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nmass);
    VectorAssign(nodec, this->nvof);
    VectorAssign(nodec * 3, this->nmome);
    VectorAssign(nodec * 3, this->nforce);
    for (int m = 0; m < nelem; m++) {

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                this->StandardVarP2G(pid, nid, sfi, this->mass, this->nmass);
                this->StandardVarP2G(pid, nid, sfi, this->vol, this->nvof);
                this->VelP2G(pid, nid, sfi);
                this->AccelP2G(pid, nid, sfi);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nmass, 0);
    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->nmome, {nuc, nvc, nwc});
    NodeVarComm(this->nforce, {nuc, nvc, nwc});

    VectorAssign(nodec * 3, this->nvel);
    VectorAssign(nodec * 3, this->naccel);
    for (int n = 0; n < nodec; n++) {
        if (this->nmass[n] > mtol) {
            this->nvel[n + nuc] = this->nmome[n + nuc] / this->nmass[n];
            this->nvel[n + nvc] = this->nmome[n + nvc] / this->nmass[n];
            this->nvel[n + nwc] = this->nmome[n + nwc] / this->nmass[n];
            this->naccel[n + nuc] = this->nforce[n + nuc] / this->nmass[n];
            this->naccel[n + nvc] = this->nforce[n + nvc] / this->nmass[n];
            this->naccel[n + nwc] = this->nforce[n + nwc] / this->nmass[n];
        }
    }

    this->ubc.BCSetVal(nuc, this->nvel);
    this->vbc.BCSetVal(nvc, this->nvel);
    this->wbc.BCSetVal(nwc, this->nvel);
    this->ubc.BCSetZero(nuc, this->naccel);
    this->vbc.BCSetZero(nvc, this->naccel);
    this->wbc.BCSetZero(nwc, this->naccel);

    return;
}

void StabilizedMPM::Node2Particle() {

    vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);
    // ---- Newmark beta velocity & acceleration ----
    this->PredictNewmarkBetaVelAndAccel(nvel_k, naccel_k);

    this->CommitNodalKinematics(nvel_k, naccel_k);

    std::vector<std::array<double, 3>> accel_old;
    VectorAssign(this->num, accel_old);
    if (this->solswitch == MapScheme::FLIP) {
        for (int n = 0; n < this->num; n++) { accel_old[n] = this->accel[n]; }
    }

    VectorAssign(this->num, this->pres);
    std::vector<std::array<double, 3>> displ;
    VectorAssign(this->num, displ);
    this->PICFamilyAccelReset();
    this->PICFamilyVelReset();
    if (this->solswitch == MapScheme::APIC) { VectorAssign(this->num, this->apic.inv_Dmat); }

    int nenode;
    vector<int> ncm;
    vector<double> sf;
    vector<std::array<double, 3>> dsf;

    for (int m = 0; m < nelem; m++) {

        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            std::array<std::array<double, 3>, 3> ADp{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
                double dsfi1 = dsf[ni][0];
                double dsfi2 = dsf[ni][1];
                double dsfi3 = dsf[ni][2];
                this->pic.VarG2P(pid, nid, sfi, this->pres, this->npres); // ---- Pressure ----
                this->pic.VarG2P(pid, nid, sfi, displ, this->ndispl);     // ---- Displacement ----
                this->PICFamilyAccelG2P(pid, ni, nid, sfi, dsf);
                this->PICFamilyVelG2P(pid, ni, nid, sfi, dsf);
                if (this->solswitch == MapScheme::APIC) { this->apic.DmatG2P(pid, nid, sfi, this->coord, ADp); }
            }
            if (this->solswitch == MapScheme::APIC) { this->apic.InvDmatG2P(pid, ADp); }
            pid = this->idp2p[pid];
        }
    }

    vector<array<double, 3>> delta_corr;
    delta_corr = this->DeltaCorrectionParticleShifting();
    // particle_shifting_SPH_like(wp, wp_move, spring_force);

    this->CommitParticleKinematics(accel_old, displ, delta_corr);

    return;
}
