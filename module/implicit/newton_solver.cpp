#include "module/implicit/newton_solver.h"

void SolveNewton(NonlinearProblem &problem, const TimeState &ts, int max_it) {

    problem.BeginNewtonStep();

    int ksp_its = 0;
    for (int NR_it = 0; NR_it <= max_it; NR_it++) {

        problem.AssembleNewtonSystem(ts);

        if (problem.CheckNewtonConvergence(NR_it, ksp_its)) { break; }

        ksp_its = problem.SolveNewtonCorrection(NR_it, ts);

        problem.UpdateNRIncrement();
    }

    problem.EndNewtonStep();

    return;
}
