#include "map_and_interpolate.h"
#include "cal_mat.h"
#include "dataset.h"
#include "material_point.h"
#include "mesh.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

void MaterialPoint::VelP2G(int pid, int nid, double sfi) {

    if (this->solswitch == MapScheme::FLIP || this->solswitch == MapScheme::PIC) {
        this->StandardVarP2G(pid, nid, sfi, this->mass, this->vel, this->nmome);
    } else if (this->solswitch == MapScheme::TPIC) {
        this->tpic.VarP2G(pid, nid, sfi, this->tpic.vel_grad, this->mass, this->vel, this->coord, this->nmome);
    } else if (this->solswitch == MapScheme::APIC) {
        this->apic.VarP2G(pid, nid, sfi, this->apic.vel_Bmat, this->mass, this->vel, this->coord, this->nmome);
    }

    return;
}

void MaterialPoint::PICFamilyVelG2P(int pid, int ni, int nid, double sfi,
                                    const std::vector<std::array<double, 3>> &dsf) {

    if (this->solswitch != MapScheme::FLIP) { this->pic.VarG2P(pid, nid, sfi, this->vel, this->nvel); }
    if (this->solswitch == MapScheme::TPIC) {
        this->tpic.VarGradG2P(pid, nid, ni, dsf, this->tpic.vel_grad, this->nvel);
    } else if (this->solswitch == MapScheme::APIC) {
        this->apic.VarBmatG2P(pid, nid, sfi, this->apic.vel_Bmat, this->nvel, this->coord);
    }

    return;
}

void MaterialPoint::AccelP2G(int pid, int nid, double sfi) {

    if (this->solswitch == MapScheme::FLIP || this->solswitch == MapScheme::PIC) {
        this->StandardVarP2G(pid, nid, sfi, this->mass, this->accel, this->nforce);
    } else if (this->solswitch == MapScheme::TPIC) {
        this->tpic.VarP2G(pid, nid, sfi, this->tpic.accel_grad, this->mass, this->accel, this->coord, this->nforce);
    } else if (this->solswitch == MapScheme::APIC) {
        this->apic.VarP2G(pid, nid, sfi, this->apic.accel_Bmat, this->mass, this->accel, this->coord, this->nforce);
    }

    return;
}

void MaterialPoint::PICFamilyAccelG2P(int pid, int ni, int nid, double sfi,
                                      const std::vector<std::array<double, 3>> &dsf) {

    this->pic.VarG2P(pid, nid, sfi, this->accel, this->naccel);
    if (this->solswitch == MapScheme::TPIC) {
        this->tpic.VarGradG2P(pid, nid, ni, dsf, this->tpic.accel_grad, this->naccel);
    } else if (this->solswitch == MapScheme::APIC) {
        this->apic.VarBmatG2P(pid, nid, sfi, this->apic.accel_Bmat, this->naccel, this->coord);
    }

    return;
}
