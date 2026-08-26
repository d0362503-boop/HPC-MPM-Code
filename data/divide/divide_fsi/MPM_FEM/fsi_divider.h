#ifndef DATA_DIVIDE_FSI_MPM_FEM_FSI_DIVIDER_H_
#define DATA_DIVIDE_FSI_MPM_FEM_FSI_DIVIDER_H_

#include "data/divide/data_partitioner.h"

#include "module/bc.h"
#include "module/material_point.h"

#include <string>

class MPMFEMFSIDivider : public DataPartitioner {
  protected:
    std::string CaseName() const override;

    void LoadBoundaryData(std::ifstream &infile) override;

    void LoadPointData(std::ifstream &infile) override;

    void PartitionProcess(int rank_id) override;

    void WriteBoundaryData(std::ofstream &outfile) override;

    void WritePointData(std::ofstream &outfile) override;

  private:
    MaterialPoint solid_points_;            // loaded solid particles
    MaterialPoint fluid_mesh_;              // loaded fluid mesh BCs
    MaterialPoint solid_partition_points_;  // partitioned solid particles
    MaterialPoint fluid_partition_mesh_;    // partitioned fluid mesh BCs
};

#endif // DATA_DIVIDE_FSI_FSI_DIVIDER_H_
