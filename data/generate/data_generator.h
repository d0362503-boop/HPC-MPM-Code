#ifndef DATA_DATA_GENERATOR_H_
#define DATA_DATA_GENERATOR_H_

#include <fstream>
#include <string>

class DataGenerator {
  public:
    virtual ~DataGenerator() = default;

    int Run();

  protected:
    virtual std::string CaseName() const = 0;

    virtual void LoadInput() = 0;

    virtual void BuildData() = 0;

    virtual void WriteBcData(std::ofstream &outfile) = 0;

    virtual void WriteTextOutputs() = 0;

    virtual void WriteVisualizationOutputs() = 0;

    void InitializeGridGeometry();

    void WriteGridDataFile(const std::string &filename);
};

#endif // DATA_DATA_GENERATOR_H_
