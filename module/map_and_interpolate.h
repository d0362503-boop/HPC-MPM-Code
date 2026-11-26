#pragma once

#include "module/cal_mat.h"
#include "module/dataset.h"
#include "module/mesh.h"
#include <array>
#include <vector>

class FLIP {
  public:
    /**
     * @brief Bind the scheme to its owner's control points and run parameters.
     * @param cmesh  Control points of the owning material-point set (field offsets).
     * @param params Run parameters of the run (time step for FLIP increments).
     */
    FLIP(const ControlMesh &cmesh, const RunParams &params) : cmesh_(cmesh), params_(params) {}

    // --- For explicit ---
    template <typename T>
    void VarG2P(int pid, int nid, double sfi, std::vector<T> &point_var, const std::vector<T> &node_var) {

        point_var[pid] += sfi * this->params_.dt * node_var[nid];

        return;
    }

    template <typename T>
    void VarG2P(int pid, int nid, double sfi, std::vector<std::array<T, 3>> &point_var,
                const std::vector<T> &node_var) {

        point_var[pid][0] += sfi * this->params_.dt * node_var[nid + this->cmesh_.nuc];
        point_var[pid][1] += sfi * this->params_.dt * node_var[nid + this->cmesh_.nvc];
        point_var[pid][2] += sfi * this->params_.dt * node_var[nid + this->cmesh_.nwc];

        return;
    }

  private:
    const ControlMesh &cmesh_; // owner control points
    const RunParams &params_;  // run parameters
};

class PIC {
  public:
    /**
     * @brief Bind the scheme to its owner's control points.
     * @param cmesh Control points of the owning material-point set (field offsets and coordinates).
     */
    explicit PIC(const ControlMesh &cmesh) : cmesh_(cmesh) {}

    template <typename T>
    void VarG2P(int pid, int nid, double sfi, std::vector<T> &point_var, const std::vector<T> &node_var) {

        point_var[pid] += sfi * node_var[nid];

        return;
    }

    template <typename T>
    void VarG2P(int pid, int nid, double sfi, std::vector<std::array<T, 3>> &point_var,
                const std::vector<T> &node_var) {

        point_var[pid][0] += sfi * node_var[nid + this->cmesh_.nuc];
        point_var[pid][1] += sfi * node_var[nid + this->cmesh_.nvc];
        point_var[pid][2] += sfi * node_var[nid + this->cmesh_.nwc];

        return;
    }

  private:
    const ControlMesh &cmesh_; // owner control points
};

class TPIC {
  public:
    /**
     * @brief Bind the scheme to its owner's control points.
     * @param cmesh Control points of the owning material-point set (field offsets and coordinates).
     */
    explicit TPIC(const ControlMesh &cmesh) : cmesh_(cmesh) {}

    std::vector<std::array<std::array<double, 3>, 3>> vel_grad, accel_grad; // ----- For TPIC -----

    template <typename T>
    void VarP2G(int pid, int nid, double sfi, std::vector<std::array<std::array<T, 3>, 3>> &point_var_grad,
                const std::vector<T> &point_mass, const std::vector<std::array<T, 3>> &point_var,
                const std::vector<std::array<T, 3>> &point_coord, std::vector<T> &node_var) {

        std::array<double, 3> vel_grad_temp{};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                vel_grad_temp[i] += point_var_grad[pid][i][j] * (this->cmesh_.xyc[nid][j] - point_coord[pid][j]);
            }
        }

        node_var[nid + this->cmesh_.nuc] += sfi * point_mass[pid] * (point_var[pid][0] + vel_grad_temp[0]);
        node_var[nid + this->cmesh_.nvc] += sfi * point_mass[pid] * (point_var[pid][1] + vel_grad_temp[1]);
        node_var[nid + this->cmesh_.nwc] += sfi * point_mass[pid] * (point_var[pid][2] + vel_grad_temp[2]);

        return;
    }

    template <typename T>
    void VarGradG2P(int pid, int nid, int ni, const std::vector<std::array<T, 3>> &dsf,
                    std::vector<std::array<std::array<T, 3>, 3>> &point_var_grad, const std::vector<T> &node_var) {

        double dsf1 = dsf[ni][0];
        double dsf2 = dsf[ni][1];
        double dsf3 = dsf[ni][2];

        point_var_grad[pid][0][0] += dsf1 * node_var[nid + this->cmesh_.nuc];
        point_var_grad[pid][0][1] += dsf2 * node_var[nid + this->cmesh_.nuc];
        point_var_grad[pid][0][2] += dsf3 * node_var[nid + this->cmesh_.nuc];
        point_var_grad[pid][1][0] += dsf1 * node_var[nid + this->cmesh_.nvc];
        point_var_grad[pid][1][1] += dsf2 * node_var[nid + this->cmesh_.nvc];
        point_var_grad[pid][1][2] += dsf3 * node_var[nid + this->cmesh_.nvc];
        point_var_grad[pid][2][0] += dsf1 * node_var[nid + this->cmesh_.nwc];
        point_var_grad[pid][2][1] += dsf2 * node_var[nid + this->cmesh_.nwc];
        point_var_grad[pid][2][2] += dsf3 * node_var[nid + this->cmesh_.nwc];

        return;
    }

  private:
    const ControlMesh &cmesh_; // owner control points
};

class APIC {
  public:
    /**
     * @brief Bind the scheme to its owner's control points.
     * @param cmesh Control points of the owning material-point set (field offsets and coordinates).
     */
    explicit APIC(const ControlMesh &cmesh) : cmesh_(cmesh) {}

    std::vector<std::array<std::array<double, 3>, 3>> vel_Bmat, accel_Bmat, inv_Dmat; // ----- For APIC -----

    template <typename T>
    void VarP2G(int pid, int nid, double sfi, std::vector<std::array<std::array<T, 3>, 3>> &point_var_Bmat,
                const std::vector<T> &point_mass, const std::vector<std::array<T, 3>> &point_var,
                const std::vector<std::array<T, 3>> &point_coord, std::vector<T> &node_var) {

        std::array<std::array<double, 3>, 3> BD_vel_temp{};
        std::array<double, 3> vel_grad_temp{};

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    BD_vel_temp[i][j] += point_var_Bmat[pid][i][k] * this->inv_Dmat[pid][k][j];
                }
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) { vel_grad_temp[i] += BD_vel_temp[i][j] * (this->cmesh_.xyc[nid][j] - point_coord[pid][j]); }
        }

        node_var[nid + this->cmesh_.nuc] += sfi * point_mass[pid] * (point_var[pid][0] + vel_grad_temp[0]);
        node_var[nid + this->cmesh_.nvc] += sfi * point_mass[pid] * (point_var[pid][1] + vel_grad_temp[1]);
        node_var[nid + this->cmesh_.nwc] += sfi * point_mass[pid] * (point_var[pid][2] + vel_grad_temp[2]);

        return;
    }

    template <typename T>
    void VarBmatG2P(int pid, int nid, double sfi, std::vector<std::array<std::array<T, 3>, 3>> &point_var_Bmat,
                    const std::vector<T> &node_var, const std::vector<std::array<T, 3>> &point_coord) {

        point_var_Bmat[pid][0][0] += sfi * node_var[nid + this->cmesh_.nuc] * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        point_var_Bmat[pid][0][1] += sfi * node_var[nid + this->cmesh_.nuc] * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        point_var_Bmat[pid][0][2] += sfi * node_var[nid + this->cmesh_.nuc] * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);
        point_var_Bmat[pid][1][0] += sfi * node_var[nid + this->cmesh_.nvc] * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        point_var_Bmat[pid][1][1] += sfi * node_var[nid + this->cmesh_.nvc] * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        point_var_Bmat[pid][1][2] += sfi * node_var[nid + this->cmesh_.nvc] * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);
        point_var_Bmat[pid][2][0] += sfi * node_var[nid + this->cmesh_.nwc] * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        point_var_Bmat[pid][2][1] += sfi * node_var[nid + this->cmesh_.nwc] * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        point_var_Bmat[pid][2][2] += sfi * node_var[nid + this->cmesh_.nwc] * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);

        return;
    }

    template <typename T>
    void DmatG2P(int pid, int nid, double sfi, const std::vector<std::array<T, 3>> &point_coord,
                 std::array<std::array<T, 3>, 3> &ADp) {

        ADp[0][0] += sfi * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]) * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        ADp[0][1] += sfi * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]) * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        ADp[0][2] += sfi * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]) * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);
        ADp[1][0] += sfi * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]) * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        ADp[1][1] += sfi * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]) * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        ADp[1][2] += sfi * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]) * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);
        ADp[2][0] += sfi * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]) * (this->cmesh_.xyc[nid][0] - point_coord[pid][0]);
        ADp[2][1] += sfi * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]) * (this->cmesh_.xyc[nid][1] - point_coord[pid][1]);
        ADp[2][2] += sfi * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]) * (this->cmesh_.xyc[nid][2] - point_coord[pid][2]);

        return;
    }

    template <typename T> void InvDmatG2P(int pid, const std::array<std::array<T, 3>, 3> &ADp) {

        std::array<std::array<T, 3>, 3> inv_Dmat_temp{};

        inv_Dmat_temp = InvMat3(ADp);
        this->inv_Dmat[pid] = inv_Dmat_temp;

        return;
    }

  private:
    const ControlMesh &cmesh_; // owner control points
};
