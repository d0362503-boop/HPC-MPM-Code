#include "data/divide/divide_fluid/fluid_divider.h"

#include <iomanip>

#include "module/data_io.h"
#include "module/dataset.h"

std::string FluidDivider::CaseName() const { return "fluid"; }

void FluidDivider::LoadBoundaryData(std::ifstream &infile) {
    this->fluid_points_.ubc.BCInput(infile);
    this->fluid_points_.vbc.BCInput(infile);
    this->fluid_points_.wbc.BCInput(infile);
    this->fluid_points_.pbc.BCInput(infile);

    // Inflow boundaries (node IDs only)
    this->fluid_points_.uinfbc.BCInput(infile, false);
    this->fluid_points_.vinfbc.BCInput(infile, false);
    this->fluid_points_.winfbc.BCInput(infile, false);
}

void FluidDivider::LoadPointData(std::ifstream &infile) {
    infile >> this->fluid_points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->fluid_points_.num, this->fluid_points_.id);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.matid);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.mass);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.vol0);
    VectorAssign(this->fluid_points_.num, this->fluid_points_.coord);

    if (this->fluid_points_.num == 0) { return; }

    InputVector(infile, this->fluid_points_.num, this->fluid_points_.coord);
    for (int i = 0; i < this->fluid_points_.num; ++i) {
        infile >> this->fluid_points_.id[i] >> this->fluid_points_.matid[i] >> this->fluid_points_.mass[i] >>
            this->fluid_points_.vol0[i];
        infile.ignore(1000, '\n');
    }
}

void FluidDivider::PartitionProcess(int rank_id) {
    this->PointRenumber(this->partition_points_.num, this->partition_points_.id, this->fluid_points_, rank_id);
    this->MeshPartition(rank_id);

    this->BcRenumber(this->partition_points_.ubc, this->fluid_points_.ubc, xynodec, xynodecw, this->inxyminc_,
                     this->inxymaxc_);
    this->BcRenumber(this->partition_points_.vbc, this->fluid_points_.vbc, xynodec, xynodecw, this->inxyminc_,
                     this->inxymaxc_);
    this->BcRenumber(this->partition_points_.wbc, this->fluid_points_.wbc, xynodec, xynodecw, this->inxyminc_,
                     this->inxymaxc_);
    this->BcRenumber(this->partition_points_.pbc, this->fluid_points_.pbc, xynodec, xynodecw, this->inxyminc_,
                     this->inxymaxc_);

    // Inflow BCs are indexed by background element, not control point.
    this->BcRenumber(this->partition_points_.uinfbc, this->fluid_points_.uinfbc, xyelem, xyelemw, aelemmin,
                     aelemmax, false);
    this->BcRenumber(this->partition_points_.vinfbc, this->fluid_points_.vinfbc, xyelem, xyelemw, aelemmin,
                     aelemmax, false);
    this->BcRenumber(this->partition_points_.winfbc, this->fluid_points_.winfbc, xyelem, xyelemw, aelemmin,
                     aelemmax, false);
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
}

void FluidDivider::WritePointData(std::ofstream &outfile) {
    outfile << std::setw(15) << this->partition_points_.num << "\n";

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        for (int j = 0; j < 3; ++j) { outfile << std::setw(15) << this->fluid_points_.coord[point_id][j]; }
        outfile << "\n";
    }

    for (int i = 0; i < this->partition_points_.num; ++i) {
        const int point_id = this->partition_points_.id[i];
        outfile << std::setw(15) << this->fluid_points_.id[point_id] << std::setw(15)
                << this->fluid_points_.matid[point_id] << std::setw(15) << this->fluid_points_.mass[point_id]
                << std::setw(15) << this->fluid_points_.vol0[point_id] << "\n";
    }
}
