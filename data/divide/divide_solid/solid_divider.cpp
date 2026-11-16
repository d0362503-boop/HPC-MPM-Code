#include "data/divide/divide_solid/solid_divider.h"

#include <iomanip>

#include "module/data_io.h"
#include "module/dataset.h"

std::string SolidDivider::CaseName() const { return "solid"; }

void SolidDivider::LoadBoundaryData(std::ifstream &infile) {

    this->points_.ubc.BCInput(infile);
    this->points_.vbc.BCInput(infile);
    this->points_.wbc.BCInput(infile);

    return;
}

void SolidDivider::LoadPointData(std::ifstream &infile) {

    infile >> this->points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->points_.num, this->points_.id);
    VectorAssign(this->points_.num, this->points_.matid);
    VectorAssign(this->points_.num, this->points_.surf_point);
    VectorAssign(this->points_.num, this->points_.mass);
    VectorAssign(this->points_.num, this->points_.vol0);
    VectorAssign(this->points_.num, this->points_.coord);

    if (this->points_.num == 0) { return; }

    InputVector(infile, this->points_.num, this->points_.coord);
    for (int i = 0; i < this->points_.num; ++i) {
        infile >> this->points_.id[i]      //
            >> this->points_.matid[i]      //
            >> this->points_.surf_point[i] //
            >> this->points_.mass[i]       //
            >> this->points_.vol0[i];
        infile.ignore(1000, '\n');
    }

    return;
}

void SolidDivider::PartitionProcess(int rank_id) {

    this->PartitionPoints(rank_id);
    this->MeshPartition(rank_id);
    this->PartitionBCs(this->inxyminc_, this->inxymaxc_);

    return;
}

void SolidDivider::PartitionPoints(int rank_id) {

    this->PointRenumber(this->partition_points_.num, this->partition_points_.id, this->points_, rank_id);

    return;
}

void SolidDivider::PartitionBCs(const std::vector<int> &cp_min, const std::vector<int> &cp_max) {

    this->BCRenumber(this->partition_points_.ubc, this->points_.ubc, xynodec, xynodecw, cp_min, cp_max);
    this->BCRenumber(this->partition_points_.vbc, this->points_.vbc, xynodec, xynodecw, cp_min, cp_max);
    this->BCRenumber(this->partition_points_.wbc, this->points_.wbc, xynodec, xynodecw, cp_min, cp_max);

    return;
}

void SolidDivider::WriteBoundaryData(std::ofstream &outfile) {

    this->partition_points_.ubc.BCOutput(outfile, "usbc");
    this->partition_points_.vbc.BCOutput(outfile, "vsbc");
    this->partition_points_.wbc.BCOutput(outfile, "wsbc");

    return;
}

void SolidDivider::WritePointData(std::ofstream &outfile) {

    outfile << std::setw(15) << this->partition_points_.num << "\n";

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        for (int j = 0; j < 3; ++j) { outfile << std::setw(15) << this->points_.coord[point_id][j]; }
        outfile << "\n";
    }

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        outfile << std::setw(15) << this->points_.id[point_id]         //
                << std::setw(15) << this->points_.matid[point_id]      //
                << std::setw(15) << this->points_.surf_point[point_id] //
                << std::setw(15) << this->points_.mass[point_id]       //
                << std::setw(15) << this->points_.vol0[point_id] << "\n";
    }

    return;
}
