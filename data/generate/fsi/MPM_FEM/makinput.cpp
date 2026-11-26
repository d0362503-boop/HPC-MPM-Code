#include "data/generate/fsi/MPM_FEM/fsi_generator.h"

int main() {

    Grid grid;

    RunParams params;

    MPMFEMFSIGenerator generator(grid, params);

    return generator.Run();
}
