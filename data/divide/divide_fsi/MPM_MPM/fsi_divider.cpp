#include "data/divide/divide_fsi/MPM_MPM/fsi_divider.h"

std::string MPMMPMFSIDivider::CaseName() const { return "MPM-MPM FSI"; }

void MPMMPMFSIDivider::LoadBoundaryData(std::ifstream &infile) {

    this->solid_.LoadBoundaryData(infile);
    this->fluid_.LoadBoundaryData(infile);

    return;
}

void MPMMPMFSIDivider::LoadPointData(std::ifstream &infile) {

    this->solid_.LoadPointData(infile);
    this->fluid_.LoadPointData(infile);

    return;
}

void MPMMPMFSIDivider::PartitionProcess(int rank_id) {

    this->PointRenumber(this->solid_.partition_points_.num, this->solid_.partition_points_.id, this->solid_.points_,
                        rank_id);
    this->PointRenumber(this->fluid_.partition_points_.num, this->fluid_.partition_points_.id, this->fluid_.points_,
                        rank_id);

    this->MeshPartition(rank_id);

    this->solid_.PartitionBCs(this->inxyminc_, this->inxymaxc_);
    this->fluid_.PartitionBCs(this->inxyminc_, this->inxymaxc_);

    return;
}

void MPMMPMFSIDivider::WriteBoundaryData(std::ofstream &outfile) {

    this->solid_.WriteBoundaryData(outfile);
    this->fluid_.WriteBoundaryData(outfile);

    return;
}

void MPMMPMFSIDivider::WritePointData(std::ofstream &outfile) {

    this->solid_.WritePointData(outfile);
    this->fluid_.WritePointData(outfile);

    return;
}
