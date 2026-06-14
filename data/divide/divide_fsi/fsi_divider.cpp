#include "fsi_divider.h"

#include <iomanip>

#include "../../../module/data_io.h"
#include "../../../module/dataset.h"

std::string FsiDivider::CaseName() const { return "fsi"; }

void FsiDivider::LoadBoundaryData(std::ifstream &infile) {
    this->solid_points_.ubc.BCInput(infile);
    this->solid_points_.vbc.BCInput(infile);
    this->solid_points_.wbc.BCInput(infile);

    this->fluid_bc_points_.ubc.BCInput(infile);
    this->fluid_bc_points_.vbc.BCInput(infile);
    this->fluid_bc_points_.wbc.BCInput(infile);
    this->fluid_bc_points_.pbc.BCInput(infile);
}

void FsiDivider::LoadPointData(std::ifstream &infile) {
    infile >> this->solid_points_.num;
    infile.ignore(1000, '\n');

    VectorAssign(this->solid_points_.num, this->solid_points_.id);
    VectorAssign(this->solid_points_.num, this->solid_points_.matid);
    VectorAssign(this->solid_points_.num, this->solid_points_.surf_point);
    VectorAssign(this->solid_points_.num, this->solid_points_.mass);
    VectorAssign(this->solid_points_.num, this->solid_points_.vol0);
    VectorAssign(this->solid_points_.num, this->solid_points_.coord);

    if (this->solid_points_.num == 0) {
        return;
    }

    InputVector(infile, this->solid_points_.num, this->solid_points_.coord);
    for (int i = 0; i < this->solid_points_.num; ++i) {
        infile >> this->solid_points_.id[i] >> this->solid_points_.matid[i]
               >> this->solid_points_.surf_point[i]
               >> this->solid_points_.mass[i]
               >> this->solid_points_.vol0[i];
        infile.ignore(1000, '\n');
    }
}

void FsiDivider::PartitionProcess(int rank_id) {
    this->PointRenumber(this->partition_points_.num, this->partition_points_.id,
                        this->solid_points_, rank_id);
    this->MeshPartition(rank_id);

    this->BcRenumber(this->partition_points_.uinfbc, this->solid_points_.ubc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.vinfbc, this->solid_points_.vbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.winfbc, this->solid_points_.wbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);

    this->BcRenumber(this->partition_points_.ubc, this->fluid_bc_points_.ubc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.vbc, this->fluid_bc_points_.vbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.wbc, this->fluid_bc_points_.wbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
    this->BcRenumber(this->partition_points_.pbc, this->fluid_bc_points_.pbc, xynodec, xynodecw,
                     this->inxyminc_, this->inxymaxc_);
}

void FsiDivider::WriteBoundaryData(std::ofstream &outfile) {
    this->partition_points_.uinfbc.BCOutput(outfile, "usbc");
    this->partition_points_.vinfbc.BCOutput(outfile, "vsbc");
    this->partition_points_.winfbc.BCOutput(outfile, "wsbc");

    this->partition_points_.ubc.BCOutput(outfile, "uwbc");
    this->partition_points_.vbc.BCOutput(outfile, "vwbc");
    this->partition_points_.wbc.BCOutput(outfile, "wwbc");
    this->partition_points_.pbc.BCOutput(outfile, "hpbc");
}

void FsiDivider::WritePointData(std::ofstream &outfile) {
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
                << std::setw(15) << this->solid_points_.surf_point[point_id]
                << std::setw(15) << this->solid_points_.mass[point_id]
                << std::setw(15) << this->solid_points_.vol0[point_id] << "\n";
    }
}
