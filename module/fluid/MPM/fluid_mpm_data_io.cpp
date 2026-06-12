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

    hsize_t local_npts = static_cast<hsize_t>(this->num);
    hsize_t global_offset = 0;
    hsize_t total_npts = local_npts;
    vtkhdf::VTKHDFWriter::ComputeGlobalInfo(local_npts, global_offset, total_npts);

    std::vector<std::array<double, 3>> points(local_npts);
    std::vector<std::array<double, 3>> velocity(local_npts);
    std::vector<double> pressure(local_npts);
    std::vector<int> pid(local_npts);
    std::vector<int> mat_id(local_npts);

    for (int i = 0; i < this->num; ++i) {
        points[i] = this->coord[i];
        velocity[i] = this->vel[i];
        pressure[i] = this->pres[i];
        pid[i] = this->id[i];
        mat_id[i] = this->matid[i];
    }

    std::vector<long long> connectivity(local_npts);
    std::vector<unsigned char> types(local_npts, 1); // VTK_VERTEX
    for (hsize_t i = 0; i < local_npts; ++i) {
        connectivity[i] = static_cast<long long>(global_offset + i);
    }

    bool is_last_rank = (myrank == nprocs - 1);
    hsize_t local_offsets_count = local_npts + (is_last_rank ? 1 : 0);
    std::vector<long long> offsets(local_offsets_count);
    for (hsize_t i = 0; i < local_npts; ++i) {
        offsets[i] = static_cast<long long>(global_offset + i);
    }
    if (is_last_rank) {
        offsets[local_npts] = static_cast<long long>(global_offset + local_npts);
    }

    vtkhdf::VTKHDFWriter writer(filename);
    writer.CreateUnstructuredGridGroup(total_npts, total_npts, total_npts);
    writer.SetTime(real_time);
    writer.WritePoints(total_npts, local_npts, global_offset, points);
    writer.WriteConnectivity(total_npts, local_npts, global_offset, connectivity);
    writer.WriteOffsets(total_npts + 1, local_offsets_count, global_offset, offsets);
    writer.WriteTypes(total_npts, local_npts, global_offset, types);

    writer.CreatePointDataGroup();
    writer.WritePointVector("Velocity", total_npts, local_npts, global_offset, velocity);
    writer.WritePointScalar("Pressure", total_npts, local_npts, global_offset, pressure);
    writer.WritePointScalarInt("ID", total_npts, local_npts, global_offset, pid);
    writer.WritePointScalarInt("MatID", total_npts, local_npts, global_offset, mat_id);
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
