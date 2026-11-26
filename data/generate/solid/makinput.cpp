#include "data/generate/solid/solid_generator.h"

int main() {

  Grid grid;

  RunParams params;

  SolidGenerator generator(grid, params);

  return generator.Run();
}
