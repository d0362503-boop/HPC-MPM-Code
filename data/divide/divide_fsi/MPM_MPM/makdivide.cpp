#include "data/divide/divide_fsi/MPM_MPM/fsi_divider.h"

int main() {
    Grid grid;
    RunParams params;
    MPMMPMFSIDivider divider(grid, params);
    return divider.Run();
}
