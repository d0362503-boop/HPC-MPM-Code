#ifndef DATA_FSI_FSI_GENERATOR_H_
#define DATA_FSI_FSI_GENERATOR_H_

#include "../data_generator.h"

#include <string>

class FsiGenerator : public DataGenerator {
  protected:
    std::string CaseName() const override;
    void LoadInput() override;
    void BuildData() override;
    void WriteBcData(std::ofstream &outfile) override;
    void WriteTextOutputs() override;
    void WriteVisualizationOutputs() override;
};

#endif // DATA_FSI_FSI_GENERATOR_H_
