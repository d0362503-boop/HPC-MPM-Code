#ifndef DATA_DIVIDE_FLUID_FLUID_DIVIDER_H_
#define DATA_DIVIDE_FLUID_FLUID_DIVIDER_H_

#include "../data_partitioner.h"

#include "../../../module/bc.h"
#include "../../../module/material_point.h"

#include <string>

class FluidDivider : public DataPartitioner {
  protected:
    std::string CaseName() const override;

    void LoadBoundaryData(std::ifstream &infile) override;

    void LoadPointData(std::ifstream &infile) override;

    void PartitionProcess(int rank_id) override;

    void WriteBoundaryData(std::ofstream &outfile) override;

    void WritePointData(std::ofstream &outfile) override;

  private:
    MaterialPoint fluid_points_;
};

#endif // DATA_DIVIDE_FLUID_FLUID_DIVIDER_H_
