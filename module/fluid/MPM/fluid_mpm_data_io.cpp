#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "../../data_io.h"
#include "../../dataset.h"
#include "../../material_point.h"
#include "../../mesh.h"
#include "../../mpi_data.h"
#include "../../vtk_hdf5.h"
#include "stabilized_mpm.h"

using namespace stabilizedmpm;

void StabilizedMPM::InputPointData(std::ifstream &infile) {
    infile >> this->num;
    infile.ignore(1000, '\n');

    this->InitializePointData();

    InputVector(infile, this->num, this->coord);

    for (int ip = 0; ip < this->num; ip++) {
        infile >> this->id[ip] >> this->matid[ip] >> this->mass[ip] >> this->vol[ip];
        infile.ignore(1000, '\n');
    }

    return;
}

void StabilizedMPM::OutputPointDataVTKHDF(int iview, int istep) {
#ifdef HAVE_HDF5
    std::string filename = outfile + "-" + std::to_string(iview) + "-w.vtkhdf";

    vtkhdf::VTKHDFWriter writer(filename);
    auto info = vtkhdf::WriteParticleTopology(writer, this->coord);
    writer.SetTime(real_time);

    writer.CreatePointDataGroup();
    writer.WritePointVector("Velocity", info.total_npts, info.local_npts, info.global_offset, this->vel);
    writer.WritePointScalar("Pressure", info.total_npts, info.local_npts, info.global_offset, this->pres);
    writer.WritePointScalar("ID", info.total_npts, info.local_npts, info.global_offset, this->id);
    writer.WritePointScalar("MatID", info.total_npts, info.local_npts, info.global_offset, this->matid);
#endif

    return;
}

void StabilizedMPM::RestartInput() {

    std::string filename = pointfile + std::to_string(myrank) + "_re.txt";

    std::ifstream reinfile;
    reinfile.open(filename);

    reinfile >> this->num;
    reinfile.ignore(1000, '\n');

    this->InputPointData(reinfile);

    if (this->num != 0) {
        InputVector(reinfile, this->num, this->coord);

        for (int ip = 0; ip < this->num; ip++) {
            reinfile >> this->id[ip] >> this->matid[ip] //
                >> this->mass[ip] >> this->vol[ip] >> this->pres[ip];
            reinfile.ignore(1000, '\n');
        }

        InputVector(reinfile, this->num, this->vel);
        InputVector(reinfile, this->num, this->accel);

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

void StabilizedMPM::RestartOutput() {

    std::string filename = pointfile + std::to_string(myrank) + "_re.txt";

    std::ofstream reoutfile;
    reoutfile.flags(std::ios::right | std::ios::scientific);
    reoutfile.open(filename);

    reoutfile << this->num << "     --- Number of water material point --- " << "\n";

    OutputVector(reoutfile, this->num, this->coord);

    for (int i = 0; i < this->num; i++) {
        reoutfile << this->id[i] << std::setw(15)    //
                  << this->matid[i] << std::setw(15) //
                  << this->mass[i] << std::setw(15)  //
                  << this->vol[i] << std::setw(15)   //
                  << this->pres[i] << "\n";
    }

    OutputVector(reoutfile, this->num, this->vel);
    OutputVector(reoutfile, this->num, this->accel);

    if (this->solswitch == MapScheme::TPIC) {
        OutputVector(reoutfile, this->num, this->tpic.vel_grad);
        OutputVector(reoutfile, this->num, this->tpic.accel_grad);
    } else if (this->solswitch == MapScheme::APIC) {
        OutputVector(reoutfile, this->num, this->apic.vel_Bmat);
        OutputVector(reoutfile, this->num, this->apic.accel_Bmat);
        OutputVector(reoutfile, this->num, this->apic.inv_Dmat);
    }

    reoutfile.close();

    return;
}
