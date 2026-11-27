#pragma once

#include "module/dataset.h"

/**
 * @brief Interface of a nonlinear system solved by Newton-Raphson iteration.
 *
 * Pure interface without data. Each implementation keeps its own Newton scratch state
 * (iterate kinematics, stress snapshot, reference residuals) and decides how to assemble,
 * check convergence, solve the linear correction and apply it. SolveNewton() fixes the order.
 */
class NonlinearProblem {
  public:
    /** @brief Destroy a nonlinear problem through the shared solver interface. */
    virtual ~NonlinearProblem() = default;

    /** @brief Prepare the step: reset increments and Newton scratch state before the first iteration. */
    virtual void BeginNewtonStep() = 0;

    /**
     * @brief Apply iterate BCs and kinematics, then assemble the tangent matrix and residual.
     * @param ts Current time state (load factor for BCs and body forces).
     */
    virtual void AssembleNewtonSystem(const TimeState &ts) = 0;

    /**
     * @brief Check the assembled residual.
     * @param NR_it Current Newton iteration.
     * @param ksp_its Krylov iterations of the previous linear correction (0 on the first iteration).
     * @return Whether the iteration has converged; failure handling at the last iteration is the implementation's.
     */
    virtual bool CheckNewtonConvergence(int NR_it, int ksp_its) = 0;

    /**
     * @brief Solve the linear system for the Newton correction.
     * @param NR_it Current Newton iteration.
     * @param ts Current time state (its step index drives preconditioner reuse).
     * @return Number of Krylov iterations used.
     */
    virtual int SolveNewtonCorrection(int NR_it, const TimeState &ts) = 0;

    /** @brief Add the linear correction to the nonlinear unknowns. */
    virtual void UpdateNRIncrement() = 0;

    /** @brief Finish the step, e.g. release linear-system storage. */
    virtual void EndNewtonStep() {}
};

/**
 * @brief Run Newton-Raphson iterations on a nonlinear problem.
 * @param problem Problem providing the Newton steps.
 * @param ts Current time state passed to every assembly.
 * @param max_it Last Newton iteration index; the loop performs at most max_it + 1 assemblies.
 */
void SolveNewton(NonlinearProblem &problem, const TimeState &ts, int max_it);
