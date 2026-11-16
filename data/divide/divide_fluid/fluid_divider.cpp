#include "data/divide/divide_fluid/fluid_divider.h"

#include <iomanip>

#include "module/data_io.h"
#include "module/dataset.h"

std::string FluidDivider::CaseName() const { return "fluid"; }

void FluidDivider::LoadBoundaryData(std::ifstream &infile) {

    this->points_.ubc.BCInput(infile);
    this->points_.vbc.BCInput(infile);
    this->points_.wbc.BCInput(infile);
    this->points_.pbc.BCInput(infile);

    // Inflow boundaries (element IDs, no values)
    this->points_.uinfbc.BCInput(infile, false);
    this->points_.vinfbc.BCInput(infile, false);
    this->points_.winfbc.BCInput(infile, false);

    return;
}

void FluidDivider::LoadPointData(std::ifstream &infile) {

    infile >> this->points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->points_.num, this->points_.id);
    VectorAssign(this->points_.num, this->points_.matid);
    VectorAssign(this->points_.num, this->points_.mass);
    VectorAssign(this->points_.num, this->points_.vol0);
    VectorAssign(this->points_.num, this->points_.coord);

    if (this->points_.num == 0) { return; }

    InputVector(infile, this->points_.num, this->points_.coord);
    for (int i = 0; i < this->points_.num; ++i) {
        infile >> this->points_.id[i] >> this->points_.matid[i] //
            >> this->points_.mass[i] >> this->points_.vol0[i];
        infile.ignore(1000, '\n');
    }

    return;
}

void FluidDivider::PartitionProcess(int rank_id) {

    this->PartitionPoints(rank_id);
    this->MeshPartition(rank_id);
    this->PartitionBCs(this->inxyminc_, this->inxymaxc_);

    return;
}

void FluidDivider::PartitionPoints(int rank_id) {

    this->PointRenumber(this->partition_points_.num, this->partition_points_.id, this->points_, rank_id);

    return;
}

void FluidDivider::PartitionBCs(const std::vector<int> &cp_min, const std::vector<int> &cp_max) {

    this->BCRenumber(this->partition_points_.ubc, this->points_.ubc, xynodec, xynodecw, cp_min, cp_max);
    this->BCRenumber(this->partition_points_.vbc, this->points_.vbc, xynodec, xynodecw, cp_min, cp_max);
    this->BCRenumber(this->partition_points_.wbc, this->points_.wbc, xynodec, xynodecw, cp_min, cp_max);
    this->BCRenumber(this->partition_points_.pbc, this->points_.pbc, xynodec, xynodecw, cp_min, cp_max);

    // Inflow BCs are indexed by background element, not control point.
    this->BCRenumber(this->partition_points_.uinfbc, this->points_.uinfbc, xyelem, xyelemw, aelemmin, aelemmax,
                     false);
    this->BCRenumber(this->partition_points_.vinfbc, this->points_.vinfbc, xyelem, xyelemw, aelemmin, aelemmax,
                     false);
    this->BCRenumber(this->partition_points_.winfbc, this->points_.winfbc, xyelem, xyelemw, aelemmin, aelemmax,
                     false);

    return;
}

void FluidDivider::WriteBoundaryData(std::ofstream &outfile) {

    this->partition_points_.ubc.BCOutput(outfile, "uwbc");
    this->partition_points_.vbc.BCOutput(outfile, "vwbc");
    this->partition_points_.wbc.BCOutput(outfile, "wwbc");
    this->partition_points_.pbc.BCOutput(outfile, "hpbc");

    // Inflow boundaries (element IDs, no values)
    this->partition_points_.uinfbc.BCOutput(outfile, "uinfbc", false);
    this->partition_points_.vinfbc.BCOutput(outfile, "vinfbc", false);
    this->partition_points_.winfbc.BCOutput(outfile, "winfbc", false);

    return;
}

void FluidDivider::WritePointData(std::ofstream &outfile) {

    outfile << std::setw(15) << this->partition_points_.num << "\n";

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        for (int j = 0; j < 3; ++j) { outfile << std::setw(15) << this->points_.coord[point_id][j]; }
        outfile << "\n";
    }

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        outfile << std::setw(15) << this->points_.id[point_id] << std::setw(15) << this->points_.matid[point_id]
                << std::setw(15) << this->points_.mass[point_id] << std::setw(15) << this->points_.vol0[point_id]
                << "\n";
    }

    return;
}
