#include "fluid_divider.h"

#include <iomanip>

#include "../../../module/data_io.h"
#include "../../../module/dataset.h"

std::string FluidDivider::CaseName() const { return "fluid"; }

void FluidDivider::LoadBoundaryData(std::ifstream &infile) {
    this->fluid_points_.ubc.BCInput(infile);
    this->fluid_points_.vbc.BCInput(infile);
    this->fluid_points_.wbc.BCInput(infile);
    this->fluid_points_.pbc.BCInput(infile);
}

void FluidDivider::LoadPointData(std::ifstream &infile) {
    infile >> this->fluid_points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->fluid_points_.num, this->fluid_points_.id);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.matid);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.mass);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.vol0);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.coord);

    if (this->fluid_points_.num == 0) {
        return;
    }

    InputVector(infile, this->fluid_points_.num, this->fluid_points_.coord);
    for (int i = 0; i < this->fluid_points_.num; ++i) {
        infile >> this->fluid_points_.id[i] >> this->fluid_points_.matid[i]
               >> this->fluid_points_.mass[i] >> this->fluid_points_.vol0[i];
        infile.ignore(1000, '\n');
    }
}

void FluidDivider::PartitionProcess(int rank_id) {
    this->PointRenumber(this->partition_points_.num, this->partition_points_.id,
                        this->fluid_points_, rank_id);
    this->MeshPartition(rank_id);

    this->BcRenumber(this->partition_points_.ubc, this->fluid_points_.ubc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.vbc, this->fluid_points_.vbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.wbc, this->fluid_points_.wbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.pbc, this->fluid_points_.pbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
}

void FluidDivider::WriteBoundaryData(std::ofstream &outfile) {
    this->partition_points_.ubc.BCOutput(outfile, "uwbc");
    this->partition_points_.vbc.BCOutput(outfile, "vwbc");
    this->partition_points_.wbc.BCOutput(outfile, "wwbc");
    this->partition_points_.pbc.BCOutput(outfile, "hpbc");
}

void FluidDivider::WritePointData(std::ofstream &outfile) {
    outfile << std::setw(15) << this->partition_points_.num << "\n";

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        for (int j = 0; j < 3; ++j) {
            outfile << std::setw(15) << this->fluid_points_.coord[point_id][j];
        }
        outfile << "\n";
    }

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        outfile << std::setw(15) << this->fluid_points_.id[point_id]
                << std::setw(15) << this->fluid_points_.matid[point_id]
                << std::setw(15) << this->fluid_points_.mass[point_id]
                << std::setw(15) << this->fluid_points_.vol0[point_id] << "\n";
    }
}
