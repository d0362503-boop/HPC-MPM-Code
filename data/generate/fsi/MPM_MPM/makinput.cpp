#include "data/generate/fsi/MPM_MPM/fsi_generator.h"

int main() {

    Grid grid;

    RunParams params;

    MPMMPMFSIGenerator generator(grid, params);

    return generator.Run();
}
