#include "data/generate/fluid/fluid_generator.h"

int main() {

  Grid grid;

  RunParams params;

  FluidGenerator generator(grid, params);

  return generator.Run();
}
