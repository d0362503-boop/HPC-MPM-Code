#include "module/cal_mat.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "module/dataset.h"
#include "module/material_point.h"
#include "module/mesh.h"
#include "module/mpi_data.h"
#include "module/solid/solid_material_point.h"
#include "module/vtk_hdf5.h"

void SolidMaterialPointBase::InputBCData(std::ifstream &infile) {

    // --- solid boundary condition ---
    this->ubc.BCInput(infile);
    this->vbc.BCInput(infile);
    this->wbc.BCInput(infile);

    if (this->do_dlb) {
        this->ubc.CaptureGlobalControlPointBC(this->grid_);
        this->vbc.CaptureGlobalControlPointBC(this->grid_);
        this->wbc.CaptureGlobalControlPointBC(this->grid_);
    }

    return;
}

void SolidMaterialPointBase::InputPointData(std::ifstream &infile) {

    infile >> this->num;
    infile.ignore(1000, '\n');

    this->InitializePointData();

    InputVector(infile, this->num, this->coord);

    for (int ip = 0; ip < this->num; ip++) {
        infile >> this->id[ip] >> this->matid[ip] >> this->surf_point[ip] >> this->mass[ip] >> this->vol0[ip];
        infile.ignore(1000, '\n');
    }
    VectorAssign(this->num, this->def_grad, kIdentity3);
    if (this->Fbar_flag) {
        VectorAssign(this->num, this->def_grad_bar, kIdentity3);
        VectorAssign(this->num, this->delta_def_grad_bar, kIdentity3);
    }
    this->vol = this->vol0;

    return;
}

void SolidMaterialPointBase::OutputPointDataVTKHDF(int iview, const TimeState &ts) {
    const RunParams &params = this->params_;

#ifdef HAVE_HDF5
    std::string filename = params.outfile + "-" + std::to_string(iview) + "-s.vtkhdf";

    std::vector<double> vm_stress(this->num);
    for (int i = 0; i < this->num; ++i) { vm_stress[i] = this->ComputeVMStress(i); }

    vtkhdf::VTKHDFWriter writer(filename);
    auto info = vtkhdf::WriteParticleTopology(writer, this->coord);
    writer.SetTime(ts.real_time);

    writer.CreatePointDataGroup();
    writer.WritePointVector("Velocity", info.total_npts, info.local_npts, info.global_offset, this->vel);
    writer.WritePointScalar("VMStress", info.total_npts, info.local_npts, info.global_offset, vm_stress);
    writer.WritePointScalar("ID", info.total_npts, info.local_npts, info.global_offset, this->id);
    writer.WritePointScalar("MatID", info.total_npts, info.local_npts, info.global_offset, this->matid);
    writer.WritePointScalar("SurfaceTag", info.total_npts, info.local_npts, info.global_offset, this->surf_point);
#endif

    return;
}

void SolidMaterialPointBase::RestartInput() {
    const RunParams &params = this->params_;

    std::string filename = params.pointfile + std::to_string(myrank) + "_s_re.txt";

    std::ifstream reinfile;
    reinfile.open(filename);

    reinfile >> this->num;
    reinfile.ignore(1000, '\n');

    this->InitializePointData();

    if (this->num != 0) {
        InputVector(reinfile, this->num, this->coord);

        for (int ip = 0; ip < this->num; ip++) {
            reinfile >> this->id[ip] >> this->matid[ip]   //
                >> this->surf_point[ip] >> this->mass[ip] //
                >> this->vol0[ip] >> this->vol[ip];
            this->det_def_grad[ip] = this->vol[ip] / this->vol0[ip];
            reinfile.ignore(1000, '\n');
        }

        InputVector(reinfile, this->num, this->vel);
        InputVector(reinfile, this->num, this->accel);
        InputVector(reinfile, this->num, this->stress);
        InputVector(reinfile, this->num, this->def_grad);

        if (this->solswitch == MapScheme::TPIC) {
            InputVector(reinfile, this->num, this->tpic.vel_grad);
            InputVector(reinfile, this->num, this->tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            InputVector(reinfile, this->num, this->apic.vel_Bmat);
            InputVector(reinfile, this->num, this->apic.accel_Bmat);
            InputVector(reinfile, this->num, this->apic.inv_Dmat);
        }
    }

    reinfile.close();

    return;
}

void SolidMaterialPointBase::RestartOutput() {
    const RunParams &params = this->params_;

    std::string filename = params.pointfile + std::to_string(myrank) + "_s_re.txt";

    std::ofstream reoutfile;
    reoutfile.flags(std::ios::right | std::ios::scientific);
    reoutfile.precision(17);
    reoutfile.open(filename);

    reoutfile << this->num << "     --- Number of solid material point --- " << "\n";

    if (this->num != 0) {
        OutputVector(reoutfile, this->num, this->coord);

        for (int i = 0; i < this->num; i++) {
            reoutfile << this->id[i]                          //
                      << std::setw(26) << this->matid[i]      //
                      << std::setw(26) << this->surf_point[i] //
                      << std::setw(26) << this->mass[i]       //
                      << std::setw(26) << this->vol0[i]       //
                      << std::setw(26) << this->vol[i] << "\n";
        }

        OutputVector(reoutfile, this->num, this->vel);
        OutputVector(reoutfile, this->num, this->accel);
        OutputVector(reoutfile, this->num, this->stress);
        OutputVector(reoutfile, this->num, this->def_grad);

        if (this->solswitch == MapScheme::TPIC) {
            OutputVector(reoutfile, this->num, this->tpic.vel_grad);
            OutputVector(reoutfile, this->num, this->tpic.accel_grad);
        } else if (this->solswitch == MapScheme::APIC) {
            OutputVector(reoutfile, this->num, this->apic.vel_Bmat);
            OutputVector(reoutfile, this->num, this->apic.accel_Bmat);
            OutputVector(reoutfile, this->num, this->apic.inv_Dmat);
        }
    }

    reoutfile.close();

    return;
}
