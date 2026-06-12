#pragma once

#include <vector>
#include <array>

class ConstitutiveModel {
public:
    bool implicit_flag = false;
    std::array<std::array<double, 6>, 6> stif_mat;

    void UpdateStress(int model_type, std::array<double, 6>& stress,
                      const std::array<std::array<double, 3>, 3>& F,
                      const std::array<std::array<double, 3>, 3>& dF,
                      const std::vector<double>& mat_prop,
                      double jac, double jac_bar);

    void MatRigid(std::array<double, 6>& stress);

    void MatLinElast(std::array<double, 6>& stress,
                     const std::array<std::array<double, 3>, 3>& delta_def_grad,
                     const std::vector<double>& mat_prop);

    void MatNeoHookean(std::array<double, 6>& stress,
                       const std::array<std::array<double, 3>, 3>& FF,
                       const std::vector<double>& mat_prop,
                       double jac, double jac_bar);

    void MatStVk(std::array<double, 6>& stress,
                 const std::array<std::array<double, 3>, 3>& FF,
                 const std::vector<double>& mat_prop, double jac);

    void MatMooneyRivlin(std::array<double, 6>& stress,
                         const std::array<std::array<double, 3>, 3>& FF,
                         const std::vector<double>& mat_prop, double jac);
};
