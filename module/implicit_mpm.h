#pragma once

#include <array>
#include <vector>

#include "module/generalized_alpha_integrator.h"
#include "module/solver/crsmat.h"

class MaterialPoint;

/**
 * @brief Implicit (generalized-alpha / Newton) time integration state and kinematics of an MPM solver.
 *
 * Owned as `implicit_` by the implicit MPM solvers (ImplicitSolidMPM, StabilizedMPM). It holds only
 * implicit state; each function takes the owning MaterialPoint whose particle and control-point arrays
 * it reads (`const`) or writes. Explicit solvers do not carry it.
 */
class ImplicitMPM {
  public:
    GeneralizedAlphaIntegrator integrator_; // time integration state

    /**
     * @brief Compute nodal velocity and acceleration from the nodal displacement
     *        increment by inverting the Newmark-beta relations.
     * @param mp       Material-point object providing `ndispl`, `nvel` and `naccel`.
     * @param nvel_k   Output nodal velocity vector (size `nodec * 3`).
     * @param naccel_k Output nodal acceleration vector (size `nodec * 3`).
     */
    void EndpointKinematics(const MaterialPoint &mp, std::vector<double> &nvel_k, std::vector<double> &naccel_k) const noexcept;

    /**
     * @brief Apply prescribed Dirichlet values to a nodal displacement increment.
     * @param mp     Material-point object providing the velocity BCs and nodal kinematics.
     * @param ndispl Nodal displacement increment modified in-place at constrained DOFs.
     */
    void ApplyDisplacementBC(MaterialPoint &mp, std::vector<double> &ndispl);

    /**
     * @brief Lumped-mass contribution of one particle to the Newton tangent.
     * @param mp  Material-point object providing particle masses.
     * @param pid Particle index.
     * @param sfi Shape-function value at the control point.
     * @return Diagonal inertia term `da_du * sfi * mass[pid]`.
     */
    double LumpedMass(const MaterialPoint &mp, int pid, double sfi) const noexcept;

    /**
     * @brief Add nodal inertial forces to the implicit system RHS.
     * @param mp Material-point object providing nodal masses.
     * @param mat System receiving the inertial force residual.
     * @param naccel Nodal acceleration at the inertia evaluation time level.
     * @param offsets RHS offsets for the x, y and z momentum components.
     */
    void AddInertialForce(const MaterialPoint &mp, CrsMat &mat, const std::vector<double> &naccel,
                          const std::vector<int> &offsets) const;

    /**
     * @brief Commit the converged nodal velocity and acceleration for the time step.
     * @param mp       Material-point object whose `nvel` and `naccel` are overwritten.
     * @param nvel_k   Converged nodal velocity vector.
     * @param naccel_k Converged nodal acceleration vector.
     */
    void CommitNodalKinematics(MaterialPoint &mp, const std::vector<double> &nvel_k, const std::vector<double> &naccel_k) const;

    /**
     * @brief Commit particle kinematics: update velocity via the Newmark-beta acceleration blend (FLIP),
     *        update position, and apply the optional particle-shifting correction.
     * @param mp Material-point object whose particle velocities and coordinates are updated.
     * @param accel_old Particle acceleration from the previous step.
     * @param disp Particle displacement increments interpolated from control points.
     * @param disp_corr Per-particle shifting displacements; empty or zero disables shifting.
     */
    void CommitParticleKinematics(MaterialPoint &mp, const std::vector<std::array<double, 3>> &accel_old,
                                  const std::vector<std::array<double, 3>> &disp,
                                  const std::vector<std::array<double, 3>> &disp_corr = {}) const;
};
