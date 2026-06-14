#ifndef DATA_DIVIDE_DATA_PARTITIONER_H_
#define DATA_DIVIDE_DATA_PARTITIONER_H_

#include <fstream>
#include <string>
#include <vector>

#include "../../module/bc.h"
#include "../../module/material_point.h"

class DataPartitioner {
  public:
    virtual ~DataPartitioner() = default;

    int Run();

  protected:
    virtual std::string CaseName() const = 0;

    virtual void LoadBoundaryData(std::ifstream &infile) = 0;

    virtual void LoadPointData(std::ifstream &infile) = 0;

    virtual void PartitionProcess(int rank_id) = 0;

    virtual void WriteBoundaryData(std::ofstream &outfile) = 0;

    virtual void WritePointData(std::ofstream &outfile) = 0;

    void InputMeshData(std::ifstream &infile);

    void OutputMeshData(std::ofstream &outfile, int rank_id) const;

    void PartitionInitialDataset();

    void BcRenumber(BoundaryCondition &partition_bc,
                    const BoundaryCondition &source_bc,
                    const std::vector<int> &local_counts,
                    const std::vector<int> &global_counts,
                    const std::vector<int> &local_min,
                    const std::vector<int> &local_max,
                    bool compute_values = true) const;

    void PointRenumber(int &local_num, std::vector<int> &global_id,
                       MaterialPoint &point, int rank_id) const;

    void MeshPartition(int rank_id);

    MaterialPoint partition_points_;
    std::vector<int> nexyr1_{3}, nexyr2_{3};
    std::vector<int> nr1_{3}, nr2_{3};
    std::vector<int> nxyr1_{3}, nxyr2_{3};
    std::vector<double> drxy1_{3}, drxy2_{3}, bound12_{3};
    std::vector<int> inxyminc_{3}, inxymaxc_{3}, inxyminl_{3}, inxymaxl_{3};

  private:
    void LoadPartitionParameters(std::ifstream &infile);

    void CreateOutputDirectories() const;

    std::string grid_input_file_;
    std::string point_input_file_;
    std::string grid_output_prefix_;
    std::string point_output_prefix_;
};

#endif // DATA_DIVIDE_DATA_PARTITIONER_H_
