#ifndef DATA_FLUID_FLUID_GENERATOR_H_
#define DATA_FLUID_FLUID_GENERATOR_H_

#include "../data_generator.h"

#include <string>

class FluidGenerator : public DataGenerator {
  protected:
    std::string CaseName() const override;

    void LoadInput() override;

    void BuildData() override;

    void WriteBcData(std::ofstream &outfile) override;

    void WriteTextOutputs() override;

    void WriteVisualizationOutputs() override;
};

#endif // DATA_FLUID_FLUID_GENERATOR_H_
