#pragma once

#include "dataset.h"
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

inline int nprocs, myrank;
inline int isb, isubc, isubl;
inline std::vector<int> naid, nsbc, nsubc, nsbl, nsubl, nxyr(3);
inline std::vector<double> dbc, dbl;

template <typename T> struct MPIDatatypeCheck {
    static MPI_Datatype GetType() {
        std::cout << "MPI_datatype_check_error" << "\n";
        exit(1);
    }
};

template <> struct MPIDatatypeCheck<int> {
    static MPI_Datatype GetType() { return MPI_INT; }
};

template <> struct MPIDatatypeCheck<double> {
    static MPI_Datatype GetType() { return MPI_DOUBLE; }
};

/**
 * @brief Exchange one nodal field component on overlap control points and add received values back into `dat`.
 * @param dat Nodal data array that stores all components in contiguous blocks.
 * @param offset Component offset inside `dat`, such as `0`, `nuc`, `nvc`, or `nwc`.
 */
template <typename T> void NodeVarComm(std::vector<T> &dat, int offset) {

    std::vector<T> sdn(isubc);
    std::vector<T> rvn(isubc);
    std::vector<MPI_Request> irqs(isb);
    std::vector<MPI_Request> irqr(isb);
    MPI_Status status;

    int ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubc[j];
            sdn[j] = dat[jj + offset];
        }
    }

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        int npro = naid[i];
        MPI_Isend(&sdn[io], nsbc[i], MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqs[i]);
        MPI_Irecv(&rvn[io], nsbc[i], MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqr[i]);
    }
    for (int i = 0; i < isb; i++) {
        MPI_Wait(&irqs[i], &status);
        MPI_Wait(&irqr[i], &status);
    }

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubc[j];
            dat[jj + offset] += rvn[j];
        }
    }

    return;
}

/**
 * @brief Exchange multiple nodal field components on overlap control points and add received values back into `dat`.
 * @param dat Nodal data array that stores all components in contiguous blocks.
 * @param offsets Component offsets inside `dat` for the fields packed into one communication step.
 */
template <typename T> void NodeVarComm(std::vector<T> &dat, const std::vector<int> &offsets) {
    int ndof = static_cast<int>(offsets.size());
    if (ndof == 1) {
        NodeVarComm(dat, offsets[0]);
        return;
    }

    std::vector<T> sdn(isubc * ndof);
    std::vector<T> rvn(isubc * ndof);
    std::vector<MPI_Request> irqs(isb);
    std::vector<MPI_Request> irqr(isb);
    MPI_Status status;

    int ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubc[j];
            for (int d = 0; d < ndof; d++) { sdn[j * ndof + d] = dat[jj + offsets[d]]; }
        }
    }

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        int npro = naid[i];
        MPI_Isend(&sdn[io * ndof], nsbc[i] * ndof, MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqs[i]);
        MPI_Irecv(&rvn[io * ndof], nsbc[i] * ndof, MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqr[i]);
    }
    for (int i = 0; i < isb; i++) {
        MPI_Wait(&irqs[i], &status);
        MPI_Wait(&irqr[i], &status);
    }

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbc[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubc[j];
            for (int d = 0; d < ndof; d++) { dat[jj + offsets[d]] += rvn[j * ndof + d]; }
        }
    }

    return;
}

template <typename T> void NodeVarCommLowerOrder(std::vector<T> &dat, int nn) {

    std::vector<T> sdn(isubc);
    std::vector<T> rvn(isubc);
    std::vector<MPI_Request> irqs(isb);
    std::vector<MPI_Request> irqr(isb);
    MPI_Status status;

    int ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbl[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubl[j];
            sdn[j] = dat[jj + nn];
        }
    }

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbl[i];
        int npro = naid[i];
        MPI_Isend(&sdn[io], nsbl[i], MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqs[i]);
        MPI_Irecv(&rvn[io], nsbl[i], MPIDatatypeCheck<T>::GetType(), npro, 1, MPI_COMM_WORLD, &irqr[i]);
    }
    for (int i = 0; i < isb; i++) {
        MPI_Wait(&irqs[i], &status);
        MPI_Wait(&irqr[i], &status);
    }
    // MPI_Waitall(isb, irqs.data(), statuses.data());
    // MPI_Waitall(isb, irqr.data(), statuses.data());

    ii = 0;
    for (int i = 0; i < isb; i++) {
        int io = ii;
        ii = io + nsbl[i];
        for (int j = io; j < ii; j++) {
            int jj = nsubl[j];
            dat[jj + nn] += rvn[j];
        }
    }

    return;
}

class ParticleCommunication {
  public:
    // --- Particle migration counters (see DetermineParticleRank in mpi_data.cpp) ---
    int nnmp; // Native Material Points: particles staying on this rank (no MPI transfer)
    int nrmp; // Removed Material Points: particles leaving this rank (out-of-bounds, removed)
    int nmps; // Material Points to Send: total particles sent to neighbor ranks (sum of nsp)
    int nrps; // Material Points to Receive: total particles received from neighbor ranks (sum of nrp)
    // -----------------------------------------------------------------------------
    std::vector<int> nsp, nrp, ofp_id;
    std::vector<int> idmp, idnsp;

    template <typename T> void PointVarSendrecv(std::vector<T> &bufs, std::vector<T> &bufr, const int idm) {

        std::vector<MPI_Request> irqs(isb);
        std::vector<MPI_Request> irqr(isb);
        MPI_Status status;

        int sicounts = 0;
        int sicountr = 0;
        for (int i = 0; i < isb; i++) {
            int ncomid = naid[i];
            int icounts = this->nsp[i] * idm;
            int icountr = this->nrp[i] * idm;
            if (icounts != 0) {
                MPI_Isend(&bufs[sicounts], icounts, MPIDatatypeCheck<T>::GetType(), ncomid, 1, MPI_COMM_WORLD,
                          &irqs[i]);
            }
            if (icountr != 0) {
                MPI_Irecv(&bufr[sicountr], icountr, MPIDatatypeCheck<T>::GetType(), ncomid, 1, MPI_COMM_WORLD,
                          &irqr[i]);
            }
            sicounts += icounts;
            sicountr += icountr;
        }

        for (int i = 0; i < isb; i++) {
            if (this->nsp[i] != 0) { MPI_Wait(&irqs[i], &status); }
            if (this->nrp[i] != 0) { MPI_Wait(&irqr[i], &status); }
        }

        return;
    }

    template <typename T, std::size_t N>
    void PointVarSendrecv(std::vector<std::array<T, N>> &bufs, std::vector<std::array<T, N>> &bufr, const int idm) {

        std::vector<MPI_Request> irqs(isb);
        std::vector<MPI_Request> irqr(isb);
        MPI_Status status;

        int sicounts = 0;
        int sicountr = 0;
        for (int i = 0; i < isb; i++) {
            int ncomid = naid[i];
            int icounts = this->nsp[i] * idm;
            int icountr = this->nrp[i] * idm;
            if (icounts != 0) {
                MPI_Isend(bufs.data() + sicounts, icounts * sizeof(std::array<T, N>), MPI_BYTE, ncomid, 1,
                          MPI_COMM_WORLD, &irqs[i]);
            }
            if (icountr != 0) {
                MPI_Irecv(bufr.data() + sicountr, icountr * sizeof(std::array<T, N>), MPI_BYTE, ncomid, 1,
                          MPI_COMM_WORLD, &irqr[i]);
            }
            sicounts += icounts;
            sicountr += icountr;
        }

        for (int i = 0; i < isb; i++) {
            if (this->nsp[i] != 0) { MPI_Wait(&irqs[i], &status); }
            if (this->nrp[i] != 0) { MPI_Wait(&irqr[i], &status); }
        }

        return;
    }

    template <typename T> void PointVarComm(int point_num, std::vector<T> &point_var) {

        std::vector<T> bufs(this->nmps, T()), bufr(this->nrps, T());
        std::vector<T> tmp(this->nnmp, T());

        int sip = 0;
        for (int i = 0; i < isb; i++) {
            for (int ip = 0; ip < this->nsp[i]; ip++) {
                int ipsip = ip + sip;
                bufs[ipsip] = point_var[this->idmp[ipsip]];
            }
            sip += this->nsp[i];
        }

        this->PointVarSendrecv(bufs, bufr, 1);

        for (int i = 0; i < this->nnmp; i++) { tmp[i] = point_var[this->idnsp[i]]; }

        VectorAssign(point_num, point_var);

        for (int i = 0; i < this->nnmp; i++) { point_var[i] = tmp[i]; }

        for (int i = 0; i < this->nrps; i++) { point_var[this->nnmp + i] = bufr[i]; }

        return;
    }

    template <typename T>
    void PointVarComm(int point_num, int ifp_num, std::vector<T> &point_var, const std::vector<T> &ifp_var) {

        this->PointVarComm(point_num, point_var);

        if (ifp_num != 0) {
            int num = this->nnmp + this->nrps;
            for (int i = 0; i < ifp_num; i++) { point_var[num + i] = ifp_var[i]; }
        }

        return;
    }

    template <typename T, std::size_t N> void PointVarComm(int point_num, std::vector<std::array<T, N>> &point_var) {

        std::vector<std::array<T, N>> bufs(this->nmps), bufr(this->nrps), tmp(this->nnmp);

        int sip = 0;
        for (int i = 0; i < isb; i++) {
            for (int ip = 0; ip < this->nsp[i]; ip++) {
                int ipsip = ip + sip;
                bufs[ipsip] = point_var[this->idmp[ipsip]];
            }
            sip += this->nsp[i];
        }

        this->PointVarSendrecv(bufs, bufr, 1);

        for (int i = 0; i < this->nnmp; i++) { tmp[i] = point_var[this->idnsp[i]]; }

        VectorAssign(point_num, point_var);

        for (int i = 0; i < this->nnmp; i++) { point_var[i] = tmp[i]; }

        for (int i = 0; i < this->nrps; i++) { point_var[this->nnmp + i] = bufr[i]; }

        return;
    }

    template <typename T, std::size_t N>
    void PointVarComm(int point_num, int ifp_num, std::vector<std::array<T, N>> &point_var,
                      const std::vector<std::array<T, N>> &ifp_var) {

        this->PointVarComm(point_num, point_var);

        if (ifp_num != 0) {
            int num = this->nnmp + this->nrps;
            for (int i = 0; i < ifp_num; i++) { point_var[num + i] = ifp_var[i]; }
        }

        return;
    }

    // ---------------------------------
};
