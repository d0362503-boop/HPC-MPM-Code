#pragma once

#include <vector>
#include <array>

class ConstitutiveModel {
public:
    bool implicit_flag = false;
    std::array<std::array<double, 6>, 6> stif_mat;

    /**
     * @brief Dispatch to the selected constitutive model and update particle stress.
     * @param model_type Constitutive model ID (-1 rigid, 0 linear elastic, 1 Neo-Hookean, 2 St. Venant-Kirchhoff, 3 Mooney-Rivlin).
     * @param stress     Voigt stress vector (updated in-place).
     * @param F          Total deformation gradient.
     * @param dF         Incremental deformation gradient.
     * @param mat_prop   Material parameters for the selected model.
     * @param jac        Determinant of `F`.
     * @param jac_bar    F-bar corrected determinant.
     */
    void UpdateStress(int model_type, std::array<double, 6>& stress,
                      const std::array<std::array<double, 3>, 3>& F,
                      const std::array<std::array<double, 3>, 3>& dF,
                      const std::vector<double>& mat_prop,
                      double jac, double jac_bar);

    /**
     * @brief Rigid-body stress model: set stress to zero.
     * @param stress Voigt stress vector (updated in-place).
     */
    void MatRigid(std::array<double, 6>& stress);

    /**
     * @brief Linear-elastic stress update.
     * @param stress          Voigt stress vector (updated in-place).
     * @param delta_def_grad  Incremental deformation gradient.
     * @param mat_prop        Material parameters `[E, nu]`.
     */
    void MatLinElast(std::array<double, 6>& stress,
                     const std::array<std::array<double, 3>, 3>& delta_def_grad,
                     const std::vector<double>& mat_prop);

    /**
     * @brief Neo-Hookean hyperelastic stress update.
     * @param stress   Voigt stress vector (updated in-place).
     * @param FF       Deformation gradient.
     * @param mat_prop Material parameters `[E, nu]`.
     * @param jac      Determinant of `FF`.
     * @param jac_bar  F-bar corrected determinant.
     */
    void MatNeoHookean(std::array<double, 6>& stress,
                       const std::array<std::array<double, 3>, 3>& FF,
                       const std::vector<double>& mat_prop,
                       double jac, double jac_bar);

    /**
     * @brief St. Venant-Kirchhoff hyperelastic stress update.
     * @param stress   Voigt stress vector (updated in-place).
     * @param FF       Deformation gradient.
     * @param mat_prop Material parameters `[E, nu]`.
     * @param jac      Determinant of `FF`.
     */
    void MatStVk(std::array<double, 6>& stress,
                 const std::array<std::array<double, 3>, 3>& FF,
                 const std::vector<double>& mat_prop, double jac);

    /**
     * @brief Mooney-Rivlin hyperelastic stress update.
     * @param stress   Voigt stress vector (updated in-place).
     * @param FF       Deformation gradient.
     * @param mat_prop Material parameters `[C1, C2, K]`.
     * @param jac      Determinant of `FF`.
     */
    void MatMooneyRivlin(std::array<double, 6>& stress,
                         const std::array<std::array<double, 3>, 3>& FF,
                         const std::vector<double>& mat_prop, double jac);
};
