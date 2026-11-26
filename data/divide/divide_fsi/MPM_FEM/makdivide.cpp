#include "data/divide/divide_fsi/MPM_FEM/fsi_divider.h"

int main() {
    Grid grid;
    RunParams params;
    MPMFEMFSIDivider divider(grid, params);
    return divider.Run();
}
