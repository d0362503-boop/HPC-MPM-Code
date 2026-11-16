#pragma once

#include <array>
#include <vector>

/** @brief Generalized-alpha parameters and Newmark kinematics in the n+alpha convention. */
class GeneralizedAlphaIntegrator {
  public:
    int ode_order = 2;               // temporal equation order
    double spec_rad = 1.0e0;         // high-frequency spectral radius
    double alpha_f = 1.0e0;          // force evaluation weight
    double alpha_m = 1.0e0;          // inertia evaluation weight
    double gamma_nb = 0.5e0;         // Newmark velocity weight
    double beta_nb = 0.25e0;         // Newmark displacement weight
    std::array<double, 6> nb_para{}; // displacement inversion coefficients

    /**
     * @brief Compute alpha and Newmark weights from the spectral radius and equation order.
     * @pre ode_order is 1 or 2 and spec_rad is in [0, 1].
     * @note Call NewmarkBetaParaSet after this method to refresh time-step coefficients.
     */
    void GeneralizedAlphaParaSet();

    /**
     * @brief Compute displacement-to-velocity and displacement-to-acceleration coefficients.
     * @param dt Physical time-step duration; must be positive.
     */
    void NewmarkBetaParaSet(double dt);

    /**
     * @brief Recover endpoint acceleration from the Newmark velocity relation.
     * @param velocity Endpoint velocity components.
     * @param old_velocity Previous-step velocity components.
     * @param old_acceleration Previous-step acceleration components.
     * @param dt Physical time-step duration; must be positive.
     * @return Endpoint accelerations in the input component ordering.
     * @pre All input arrays have equal lengths and gamma_nb is nonzero.
     */
    std::vector<double> ComputeNodeAccelFromVel(const std::vector<double> &velocity, const std::vector<double> &old_velocity,
                                                const std::vector<double> &old_acceleration, double dt) const;

    /**
     * @brief Recover endpoint velocity and acceleration from a step displacement increment.
     * @param displacement Displacement from the previous step to the endpoint iterate.
     * @param old_velocity Previous-step velocity components.
     * @param old_acceleration Previous-step acceleration components.
     * @param velocity Output endpoint velocity components, already sized by the caller.
     * @param acceleration Output endpoint acceleration components, already sized by the caller.
     * @pre All arrays have equal lengths, outputs do not alias inputs, and nb_para
     * has been initialized for the current time step.
     */
    void ComputeNodeVelAccelFromDispl(const std::vector<double> &displacement, const std::vector<double> &old_velocity,
                                      const std::vector<double> &old_acceleration, std::vector<double> &velocity,
                                      std::vector<double> &acceleration) const noexcept;
};
