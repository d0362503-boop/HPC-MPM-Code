#include "data/divide/divide_fsi/MPM_FEM/fsi_divider.h"

#include "module/dataset.h"
#include "module/mesh.h"

void MPMFEMFSIFluidDivider::LoadBoundaryData(std::ifstream &infile) {

    this->mesh_.ubc.BCInput(infile);
    this->mesh_.vbc.BCInput(infile);
    this->mesh_.wbc.BCInput(infile);
    this->mesh_.pbc.BCInput(infile);

    return;
}

void MPMFEMFSIFluidDivider::PartitionBCs(const std::vector<int> &cp_min, const std::vector<int> &cp_max) {

    DataPartitioner::BCRenumber(this->partition_mesh_.ubc, this->mesh_.ubc, xynodec, xynodecw, cp_min, cp_max);
    DataPartitioner::BCRenumber(this->partition_mesh_.vbc, this->mesh_.vbc, xynodec, xynodecw, cp_min, cp_max);
    DataPartitioner::BCRenumber(this->partition_mesh_.wbc, this->mesh_.wbc, xynodec, xynodecw, cp_min, cp_max);
    DataPartitioner::BCRenumber(this->partition_mesh_.pbc, this->mesh_.pbc, xynodec, xynodecw, cp_min, cp_max);

    return;
}

void MPMFEMFSIFluidDivider::WriteBoundaryData(std::ofstream &outfile) {

    this->partition_mesh_.ubc.BCOutput(outfile, "uwbc");
    this->partition_mesh_.vbc.BCOutput(outfile, "vwbc");
    this->partition_mesh_.wbc.BCOutput(outfile, "wwbc");
    this->partition_mesh_.pbc.BCOutput(outfile, "hpbc");

    return;
}

std::string MPMFEMFSIDivider::CaseName() const { return "MPM-FEM FSI"; }

void MPMFEMFSIDivider::LoadBoundaryData(std::ifstream &infile) {

    this->solid_.LoadBoundaryData(infile);
    this->fluid_.LoadBoundaryData(infile);

    return;
}

void MPMFEMFSIDivider::LoadPointData(std::ifstream &infile) {

    this->solid_.LoadPointData(infile);

    return;
}

void MPMFEMFSIDivider::PartitionProcess(int rank_id) {

    this->PointRenumber(this->solid_.partition_points_.num, this->solid_.partition_points_.id, this->solid_.points_,
                        rank_id);

    this->MeshPartition(rank_id);

    this->solid_.PartitionBCs(this->inxyminc_, this->inxymaxc_);
    this->fluid_.PartitionBCs(this->inxyminc_, this->inxymaxc_);

    return;
}

void MPMFEMFSIDivider::WriteBoundaryData(std::ofstream &outfile) {

    this->solid_.WriteBoundaryData(outfile);
    this->fluid_.WriteBoundaryData(outfile);

    return;
}

void MPMFEMFSIDivider::WritePointData(std::ofstream &outfile) {

    this->solid_.WritePointData(outfile);

    return;
}
