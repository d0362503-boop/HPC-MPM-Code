#include "solid_divider.h"

#include <iomanip>

#include "../../../module/data_io.h"
#include "../../../module/dataset.h"

std::string SolidDivider::CaseName() const { return "solid"; }

void SolidDivider::LoadBoundaryData(std::ifstream &infile) {
    this->solid_points_.ubc.BCInput(infile);
    this->solid_points_.vbc.BCInput(infile);
    this->solid_points_.wbc.BCInput(infile);
}

void SolidDivider::LoadPointData(std::ifstream &infile) {
    infile >> this->solid_points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->solid_points_.num, this->solid_points_.id);
    VectorAssign(this->solid_points_.num, this->solid_points_.matid);
    VectorAssign(this->solid_points_.num, this->solid_points_.mass);
    VectorAssign(this->solid_points_.num, this->solid_points_.vol0);
    VectorAssign(this->solid_points_.num, this->solid_points_.coord);

    if (this->solid_points_.num == 0) {
        return;
    }

    InputVector(infile, this->solid_points_.num, this->solid_points_.coord);
    for (int i = 0; i < this->solid_points_.num; ++i) {
        infile >> this->solid_points_.id[i] >> this->solid_points_.matid[i]
               >> this->solid_points_.mass[i] >> this->solid_points_.vol0[i];
        infile.ignore(1000, '\n');
    }
}

void SolidDivider::PartitionProcess(int rank_id) {
    this->PointRenumber(this->partition_points_.num, this->partition_points_.id,
                        this->solid_points_, rank_id);
    this->MeshPartition(rank_id);

    this->BcRenumber(this->partition_points_.ubc, this->solid_points_.ubc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.vbc, this->solid_points_.vbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.wbc, this->solid_points_.wbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
}

void SolidDivider::WriteBoundaryData(std::ofstream &outfile) {
    this->partition_points_.ubc.BCOutput(outfile, "usbc");
    this->partition_points_.vbc.BCOutput(outfile, "vsbc");
    this->partition_points_.wbc.BCOutput(outfile, "wsbc");
}

void SolidDivider::WritePointData(std::ofstream &outfile) {
    outfile << std::setw(15) << this->partition_points_.num << "\n";

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        for (int j = 0; j < 3; ++j) {
            outfile << std::setw(15) << this->solid_points_.coord[point_id][j];
        }
        outfile << "\n";
    }

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        outfile << std::setw(15) << this->solid_points_.id[point_id]
                << std::setw(15) << this->solid_points_.matid[point_id]
                << std::setw(15) << this->solid_points_.mass[point_id]
                << std::setw(15) << this->solid_points_.vol0[point_id] << "\n";
    }
}
