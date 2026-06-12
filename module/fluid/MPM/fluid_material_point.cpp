#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "stabilized_mpm.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <vector>

using namespace stabilizedmpm;

void StabilizedMPM::InitializePointData() {
    VectorAssign(this->num, this->id);
    VectorAssign(this->num, this->matid);
    VectorAssign(this->num, this->mass);
    VectorAssign(this->num, this->pres);
    VectorAssign(this->num, this->vol);

    VectorAssign(this->num, this->coord);
    VectorAssign(this->num, this->vel);
    VectorAssign(this->num, this->accel);

    if (this->solswitch == MapScheme::TPIC) {
        VectorAssign(this->num, this->tpic.vel_grad);
        VectorAssign(this->num, this->tpic.accel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        VectorAssign(this->num, this->apic.vel_Bmat);
        VectorAssign(this->num, this->apic.accel_Bmat);
        VectorAssign(this->num, this->apic.inv_Dmat);
    }

    return;
}

void StabilizedMPM::Moveparticle() {

    this->DetermineParticleRank();

    this->InflowParticles();

    if (this->par_comm_.nmps + this->par_comm_.nrps + this->par_comm_.nrmp + this->ifp.num == 0) return;

    this->num += this->par_comm_.nrps + this->ifp.num - this->par_comm_.nrmp - this->par_comm_.nmps;

    // ---- Scalar part ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->id, this->ifp.id);       // ---- Point Global ID ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->matid, this->ifp.matid); // ---- Point material ID ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->mass, this->ifp.mass);   // ---- Point mass ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->vol, this->ifp.vol);     // ---- Point volume ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->pres, this->ifp.pres);   // ---- Point pressure ----

    // ---- Vector part ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->coord, this->ifp.coord); // ---- Point coordinate ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->vel, this->ifp.vel);     // ---- Point velocity ----
    this->par_comm_.PointVarComm(this->num, this->ifp.num, this->accel,
                                 this->ifp.accel); // ---- Point acceleration ----

    // ---- Matrix part ----
    if (this->solswitch == MapScheme::TPIC) {
        // ---- TPIC acceleration gradient ----
        this->par_comm_.PointVarComm(this->num, this->ifp.num, this->tpic.accel_grad, this->ifp.tpic.accel_grad);
        // ---- TPIC velocity gredient ----
        this->par_comm_.PointVarComm(this->num, this->ifp.num, this->tpic.vel_grad, this->ifp.tpic.vel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        // ---- APIC acceleration B matrix ----
        this->par_comm_.PointVarComm(this->num, this->ifp.num, this->apic.accel_Bmat, this->ifp.apic.accel_Bmat);
        // ---- APIC velocity B matrix ----
        this->par_comm_.PointVarComm(this->num, this->ifp.num, this->apic.vel_Bmat, this->ifp.apic.vel_Bmat);
        // ---- APIC D inverse matrix ----
        this->par_comm_.PointVarComm(this->num, this->ifp.num, this->apic.inv_Dmat, this->ifp.apic.inv_Dmat);
    }

    return;
}