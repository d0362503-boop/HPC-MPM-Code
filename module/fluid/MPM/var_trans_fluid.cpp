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

using namespace stabilizedmpm;

void StabilizedMPM::Particle2Node() {

    int nenode;
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    VectorAssign(nodec, this->nmass);
    VectorAssign(nodec, this->nvof);
    VectorAssign(nodec, this->npres_old);
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
                this->StandardVarP2G(pid, nid, sfi, this->mass, this->nmass);
                this->StandardVarP2G(pid, nid, sfi, this->vol, this->nvof);
                this->StandardVarP2G(pid, nid, sfi, this->mass, this->pres, this->npres_old);
                this->VelP2G(pid, nid, sfi);
                this->AccelP2G(pid, nid, sfi);
            }
            pid = this->idp2p[pid];
        }
    }

    NodeVarComm(this->nmass, 0);
    NodeVarComm(this->nvof, 0);
    NodeVarComm(this->npres_old, 0);
    NodeVarComm(this->nmome, {nuc, nvc, nwc});
    NodeVarComm(this->nforce, {nuc, nvc, nwc});

    this->CutOffSmallNodalVar(this->npres_old, this->nmass, {0});
    this->pbc.BCSetVal(0, this->npres_old);

    this->CutOffSmallNodalVar(this->nvel, this->nmome, this->nmass, {nuc, nvc, nwc});
    this->ApplyVelocityBC(this->nvel);

    this->CutOffSmallNodalVar(this->naccel, this->nforce, this->nmass, {nuc, nvc, nwc});
    this->ApplyAccelerationBC(this->naccel);

    return;
}

void StabilizedMPM::Node2Particle() {

    std::vector<double> nvel_k(nodec * 3), naccel_k(nodec * 3);
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
    std::vector<int> ncm;
    std::vector<double> sf;
    std::vector<std::array<double, 3>> dsf;

    for (int m = 0; m < nelem; m++) {
        int pid = this->idepf[m];
        while (pid != -1) {
            std::array<double, 3> xyp = this->coord[pid];
            MakSf(m, xyp, idimc, xynodec, ncm, nenode, sf, dsf);

            std::array<std::array<double, 3>, 3> ADp{};
            for (int ni = 0; ni < nenode; ni++) {
                int nid = ncm[ni];
                double sfi = sf[ni];
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

    std::vector<std::array<double, 3>> delta_corr;
    delta_corr = this->DeltaCorrectionParticleShifting();
    // delta_corr = this->SPHLikeParticleShifting();

    this->CommitParticleKinematics(accel_old, displ, delta_corr);

    return;
}
