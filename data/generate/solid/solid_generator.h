#ifndef DATA_SOLID_SOLID_GENERATOR_H_
#define DATA_SOLID_SOLID_GENERATOR_H_

#include "../data_generator.h"

#include <string>

class SolidGenerator : public DataGenerator {
 protected:
  std::string CaseName() const override;
  void LoadInput() override;
  void BuildData() override;
  void WriteBcData(std::ofstream& outfile) override;
  void WriteTextOutputs() override;
  void WriteVisualizationOutputs() override;
};

#endif  // DATA_SOLID_SOLID_GENERATOR_H_
