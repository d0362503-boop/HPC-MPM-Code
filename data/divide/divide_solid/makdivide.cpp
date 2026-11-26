#include "data/divide/divide_solid/solid_divider.h"

int main() {
    Grid grid;
    RunParams params;
    SolidDivider divider(grid, params);
    return divider.Run();
}
