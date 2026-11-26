# Implicit core

Shared building blocks of the implicit (generalized-α / Newton) solvers. The
physics-specific Newton steps stay in each solver; linear algebra stays in
[`module/solver/`](../solver/README.md).

| File | Contents | Used by |
|------|----------|---------|
| `generalized_alpha_integrator.h/.cpp` | `GeneralizedAlphaIntegrator`: generalized-α / Newmark-β coefficients and displacement-to-velocity/acceleration relations | `ImplicitMPM`, `StabilizedFEM` |
| `implicit_mpm.h/.cpp` | `ImplicitMPM` component (owned as `implicit_`): endpoint kinematics, lumped mass, inertial force, nodal and particle commits | `ImplicitSolidMPM`, `StabilizedMPM` |
| `newton_solver.h/.cpp` | `NonlinearProblem` interface and `SolveNewton()` loop | `ImplicitSolidMPM`, `StabilizedMPM`, `MPMMPMMonolithicFSI` |

`ImplicitMPM` holds only implicit state; each function takes the owning
`MaterialPoint` and reads (`const`) or writes its arrays.

## Newton driver

```cpp
problem.BeginNewtonStep();
for NR_it = 0 .. max_it
    problem.AssembleNewtonSystem(ts);
    if problem.CheckNewtonConvergence(NR_it, ksp_its) break;
    ksp_its = problem.SolveNewtonCorrection(NR_it);
    problem.UpdateNRIncrement();
problem.EndNewtonStep();
```

`ksp_its` is the Krylov count of the previous correction (0 on the first pass).
Failure handling at the last iteration belongs to each `CheckNewtonConvergence`
(the CrsMat or FSI convergence check aborts). `EndNewtonStep()` defaults to
nothing; the single-field solvers use it to release their system storage.
