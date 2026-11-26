#include "data/divide/divide_fluid/fluid_divider.h"

int main() {
    Grid grid;
    RunParams params;
    FluidDivider divider(grid, params);
    return divider.Run();
}
