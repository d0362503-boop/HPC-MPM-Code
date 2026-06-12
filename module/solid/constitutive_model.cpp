#include "constitutive_model.h"
#include "../cal_mat.h"
#include "../dataset.h"
#include <cmath>
#include <vector>

void ConstitutiveModel::MatRigid(std::array<double, 6> &stress) {

    for (int n = 0; n < 6; n++) { stress[n] = 0.0e0; }

    if (this->implicit_flag == true) { this->stif_mat = {}; }

    return;
}

void ConstitutiveModel::MatLinElast(std::array<double, 6> &stress,
                                    const std::array<std::array<double, 3>, 3> &delta_def_grad,
                                    const std::vector<double> &mat_prop) {

    double EE = mat_prop[0];
    double nu = mat_prop[1];

    std::array<double, 3> CC{};
    CC[0] = (EE * (1.0e0 - nu)) / ((1.0e0 - 2.0e0 * nu) * (1.0e0 + nu));
    CC[1] = (EE * nu) / ((1.0e0 - 2.0e0 * nu) * (1.0e0 + nu));
    CC[2] = EE / (2.0e0 * (1.0e0 + nu));

    std::array<double, 6> deps{};
    deps[0] = delta_def_grad[0][0] - 1.0e0;
    deps[1] = delta_def_grad[1][1] - 1.0e0;
    deps[2] = delta_def_grad[2][2] - 1.0e0;
    deps[3] = delta_def_grad[1][2] + delta_def_grad[2][1];
    deps[4] = delta_def_grad[2][0] + delta_def_grad[0][2];
    deps[5] = delta_def_grad[0][1] + delta_def_grad[1][0];

    std::array<double, 6> dsig{};
    dsig[0] = CC[0] * deps[0] + CC[1] * deps[1] + CC[1] * deps[2];
    dsig[1] = CC[1] * deps[0] + CC[0] * deps[1] + CC[1] * deps[2];
    dsig[2] = CC[1] * deps[0] + CC[1] * deps[1] + CC[0] * deps[2];
    dsig[3] = CC[2] * deps[3];
    dsig[4] = CC[2] * deps[4];
    dsig[5] = CC[2] * deps[5];

    for (int n = 0; n < 6; n++) { stress[n] += dsig[n]; }

    if (this->implicit_flag == true) {
        this->stif_mat = {};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i != j) this->stif_mat[i][j] = CC[1]; // --- lamda ---
            }
            this->stif_mat[i][i] = CC[0]; // --- lamda + 2mu ---
        }
        this->stif_mat[3][3] = CC[2]; // --- mu ---
        this->stif_mat[4][4] = CC[2]; // --- mu ---
        this->stif_mat[5][5] = CC[2]; // --- mu ---
    }

    return;
}

void ConstitutiveModel::MatNeoHookean(std::array<double, 6> &stress, const std::array<std::array<double, 3>, 3> &FF,
                                      const std::vector<double> &mat_prop, double jac, double jac_bar) {

    double yng = mat_prop[0];
    double nu = mat_prop[1];

    double lmd = nu * yng / ((1.0e0 + nu) * (1.0e0 - 2.0e0 * nu));
    double mu = yng / (2.0e0 * (1.0e0 + nu));
    double lnjac = log(jac_bar);

    std::array<std::array<double, 3>, 3> be{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { be[i][j] = FF[i][0] * FF[j][0] + FF[i][1] * FF[j][1] + FF[i][2] * FF[j][2]; }
    }

    std::array<std::array<double, 3>, 3> sig{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            sig[i][j] = (mu * (be[i][j] - eye_mat[i][j]) + lmd * lnjac * eye_mat[i][j]) / jac;
        }
    }

    stress[0] = sig[0][0];
    stress[1] = sig[1][1];
    stress[2] = sig[2][2];
    stress[3] = sig[1][2];
    stress[4] = sig[0][2];
    stress[5] = sig[0][1];

    if (this->implicit_flag == true) {
        // --- Material Stiffness ---
        std::array<double, 2> CC{};
        CC[0] = lmd / jac;                // --- lamda' ---
        CC[1] = (mu - lmd * lnjac) / jac; // --- mu' ---

        std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> ccc{};
        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) {
                        ccc[i][j][k][l] = CC[0] * (eye_mat[i][j] * eye_mat[k][l]) +
                                          CC[1] * (eye_mat[i][k] * eye_mat[j][l] + eye_mat[i][l] * eye_mat[j][k]);
                    }
                }
            }
        }

        this->stif_mat = {};
        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) { this->stif_mat[idm[i][j]][idm[k][l]] = ccc[i][j][k][l]; }
                }
            }
        }
    }

    return;
}

void ConstitutiveModel::MatStVk(std::array<double, 6> &stress, const std::array<std::array<double, 3>, 3> &FF,
                                const std::vector<double> &mat_prop, double jac) {

    double yng = mat_prop[0];
    double nu = mat_prop[1];

    double lmd = nu * yng / ((1.0e0 + nu) * (1.0e0 - 2.0e0 * nu));
    double mu = yng / (2.0e0 * (1.0e0 + nu));

    std::array<std::array<double, 3>, 3> be{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { be[i][j] = FF[i][0] * FF[j][0] + FF[i][1] * FF[j][1] + FF[i][2] * FF[j][2]; }
    }

    std::array<std::array<double, 3>, 3> be2{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { be2[i][j] = be[i][0] * be[0][j] + be[i][1] * be[1][j] + be[i][2] * be[2][j]; }
    }

    double tr_be = TraceMat3(be);
    double tr_GLS = 0.5e0 * (tr_be - 3.0e0);

    std::array<std::array<double, 3>, 3> sig{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { sig[i][j] = (mu * be2[i][j] + (lmd * tr_GLS - mu) * be[i][j]) / jac; }
    }

    stress[0] = sig[0][0];
    stress[1] = sig[1][1];
    stress[2] = sig[2][2];
    stress[3] = sig[1][2];
    stress[4] = sig[0][2];
    stress[5] = sig[0][1];

    if (this->implicit_flag == true) {
        // --- Material Stiffness ---
        std::array<double, 2> CC{};
        CC[0] = lmd / jac; // --- lamda' ---
        CC[1] = mu / jac;  // --- mu' ---

        std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> ccc{};
        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) {
                        ccc[i][j][k][l] =
                            CC[0] * (be[i][j] * be[k][l]) + CC[1] * (be[i][k] * be[j][l] + be[i][l] * be[j][k]);
                    }
                }
            }
        }

        this->stif_mat = {};
        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) { this->stif_mat[idm[i][j]][idm[k][l]] = ccc[i][j][k][l]; }
                }
            }
        }
    }

    return;
}

void ConstitutiveModel::MatMooneyRivlin(std::array<double, 6> &stress, const std::array<std::array<double, 3>, 3> &FF,
                                        const std::vector<double> &mat_prop, double jac) {

    double c1 = mat_prop[0];
    double c2 = mat_prop[1];
    double bulk = mat_prop[2];

    std::array<std::array<double, 3>, 3> FF_iso{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) { FF_iso[i][j] = 1.0e0 / std::cbrt(jac) * FF[i][j]; }
    }

    std::array<std::array<double, 3>, 3> be_iso{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            be_iso[i][j] = FF_iso[i][0] * FF_iso[j][0] + FF_iso[i][1] * FF_iso[j][1] + FF_iso[i][2] * FF_iso[j][2];
        }
    }

    std::array<std::array<double, 3>, 3> be2_iso{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            be2_iso[i][j] = be_iso[i][0] * be_iso[0][j] + be_iso[i][1] * be_iso[1][j] + be_iso[i][2] * be_iso[2][j];
        }
    }

    double first_inv = TraceMat3(be_iso);
    double second_inv = 0.5e0 * (pow(first_inv, 2) - TraceMat3(be2_iso));
    double hp = bulk * log(jac) / jac;

    std::array<std::array<double, 3>, 3> sig{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            sig[i][j] = 2.0e0 * ((c1 + c2 * first_inv) * be_iso[i][j] - c2 * be2_iso[i][j]) / jac;
        }
        sig[i][i] += hp - 2.0e0 / (3.0e0 * jac) * (c1 * first_inv + 2.0e0 * c2 * second_inv);
    }

    stress[0] = sig[0][0];
    stress[1] = sig[1][1];
    stress[2] = sig[2][2];
    stress[3] = sig[1][2];
    stress[4] = sig[0][2];
    stress[5] = sig[0][1];

    if (this->implicit_flag == true) {

        std::array<std::array<double, 3>, 3> dev_sig{};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) { dev_sig[i][j] = sig[i][j] - TraceMat3(sig) / 3.0e0 * eye_mat[i][j]; }
        }

        std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> ccc_hp{}, ccc_dev{}, ccc_w{}, ccc{};

        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) {
                        double I4_sym = 0.5e0 * (eye_mat[i][k] * eye_mat[j][l] + eye_mat[i][l] * eye_mat[j][k]);

                        double Aij = first_inv * be_iso[i][j] - be2_iso[i][j];
                        double Akl = first_inv * be_iso[k][l] - be2_iso[k][l];

                        ccc_hp[i][j][k][l] = bulk / jac * (eye_mat[i][j] * eye_mat[k][l]) - 2.0e0 * hp * I4_sym;

                        ccc_dev[i][j][k][l] =
                            -2.0e0 / 3.0e0 * (dev_sig[i][j] * eye_mat[k][l] + eye_mat[i][j] * dev_sig[k][l]);

                        double term1 = 4.0e0 / (3.0e0 * jac) * (c1 * first_inv + 2.0e0 * c2 * second_inv) *
                                       (I4_sym - (1.0e0 / 3.0e0) * eye_mat[i][j] * eye_mat[k][l]);

                        double term2 = 4.0e0 * c2 / jac *
                                       (be_iso[i][j] * be_iso[k][l] -
                                        0.5e0 * (be_iso[i][k] * be_iso[j][l] + be_iso[i][l] * be_iso[j][k]));

                        double term3 = -4.0e0 * c2 / (3.0e0 * jac) * (Aij * eye_mat[k][l] + eye_mat[i][j] * Akl);

                        double term4 = 8.0e0 * c2 * second_inv / (9.0e0 * jac) * eye_mat[i][j] * eye_mat[k][l];

                        ccc_w[i][j][k][l] = term1 + term2 + term3 + term4;

                        ccc[i][j][k][l] = ccc_hp[i][j][k][l] + ccc_dev[i][j][k][l] + ccc_w[i][j][k][l];
                    }
                }
            }
        }

        this->stif_mat = {};
        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) { this->stif_mat[idm[i][j]][idm[k][l]] = ccc[i][j][k][l]; }
                }
            }
        }
    }

    return;
}

void ConstitutiveModel::UpdateStress(int model_type, std::array<double, 6> &stress,
                                     const std::array<std::array<double, 3>, 3> &F,
                                     const std::array<std::array<double, 3>, 3> &dF,
                                     const std::vector<double> &mat_prop, double jac, double jac_bar) {

    if (model_type == -1) { // --- Rigid body ---
        this->MatRigid(stress);
    } else if (model_type == 0) { // --- Linear elasticity ---
        this->MatLinElast(stress, dF, mat_prop);
    } else if (model_type == 1) { // --- Neo Hookean ---
        this->MatNeoHookean(stress, F, mat_prop, jac, jac_bar);
    } else if (model_type == 2) { // --- St.Venant-Kirchhoff ---
        this->MatStVk(stress, F, mat_prop, jac);
    } else if (model_type == 3) { // --- Mooney-Rivlin ---
        this->MatMooneyRivlin(stress, F, mat_prop, jac);
    }

    return;
}
