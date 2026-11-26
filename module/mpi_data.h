#pragma once

#include "module/DLB/mpm_dlb.h"
#include "module/dataset.h"
#include <array>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

inline int nprocs, myrank;

/// Rank topology and overlap (ghost) bookkeeping of this rank.
struct Partition {
    std::array<int, 3> nxyr{};            // ranks per direction
    int isb = 0;                          // overlapping neighbour ranks
    int isubc = 0, isubl = 0;             // shared control points / nodes in total
    std::vector<int> naid;                // neighbour rank ids
    std::vector<int> nsbc, nsubc;         // shared control points per rank / ids
    std::vector<int> nsbl, nsubl;         // shared nodes per rank / ids
    std::vector<double> dbc, dbl;         // control-point / node overlap weights
    std::vector<mpm_dlb::Region> regions; // element regions of all ranks
};

/// MPI datatype for `T`; unsupported types fail at compile time (incomplete type).
template <typename T> struct MPIDatatypeCheck;

template <> struct MPIDatatypeCheck<int> {
    static MPI_Datatype GetType() { return MPI_INT; }
    static constexpr int kCount = 1; // scalars per element
};

template <> struct MPIDatatypeCheck<unsigned long long> {
    static MPI_Datatype GetType() { return MPI_UNSIGNED_LONG_LONG; }
    static constexpr int kCount = 1; // scalars per element
};

template <> struct MPIDatatypeCheck<double> {
    static MPI_Datatype GetType() { return MPI_DOUBLE; }
    static constexpr int kCount = 1; // scalars per element
};

/// `std::array<T, N>` is sent as N values of T: array<double, 3> is 3 MPI_DOUBLE, a 3x3 matrix is 9.
template <typename T, std::size_t N> struct MPIDatatypeCheck<std::array<T, N>> {
    static MPI_Datatype GetType() { return MPIDatatypeCheck<T>::GetType(); }
    static constexpr int kCount = static_cast<int>(N) * MPIDatatypeCheck<T>::kCount; // scalars per element
};

namespace internal {

/**
 * @brief Exchange contiguous per-neighbor blocks with every neighbor rank.
 * @param peers Neighbor rank IDs.
 * @param send_counts Number of entries sent to each neighbor.
 * @param recv_counts Number of entries received from each neighbor.
 * @param per_count Number of buffer elements per entry, such as components per node.
 * @param sbuf Send buffer with the blocks for each neighbor stored in `peers` order.
 * @param rbuf Receive buffer, already sized by the caller, filled in `peers` order.
 * @note Neighbors with zero count are skipped; their requests stay MPI_REQUEST_NULL for MPI_Waitall.
 */
template <typename E>
void PeerSendrecv(const std::vector<int> &peers, const std::vector<int> &send_counts, const std::vector<int> &recv_counts,
                  int per_count, const std::vector<E> &sbuf, std::vector<E> &rbuf) {

    const int npeer = static_cast<int>(peers.size());
    std::vector<MPI_Request> irq(2 * npeer, MPI_REQUEST_NULL);

    int so = 0;
    int ro = 0;
    for (int i = 0; i < npeer; i++) {
        const int sc = send_counts[i] * per_count; // elements to send
        const int rc = recv_counts[i] * per_count; // elements to receive
        if (sc != 0) {
            MPI_Isend(sbuf.data() + so, sc * MPIDatatypeCheck<E>::kCount, MPIDatatypeCheck<E>::GetType(), peers[i], 1,
                      MPI_COMM_WORLD, &irq[2 * i]);
        }
        if (rc != 0) {
            MPI_Irecv(rbuf.data() + ro, rc * MPIDatatypeCheck<E>::kCount, MPIDatatypeCheck<E>::GetType(), peers[i], 1,
                      MPI_COMM_WORLD, &irq[2 * i + 1]);
        }
        so += sc;
        ro += rc;
    }
    MPI_Waitall(2 * npeer, irq.data(), MPI_STATUSES_IGNORE);

    return;
}

/**
 * @brief Shared implementation of NodeVarComm; call NodeVarComm instead.
 * @param part Rank overlap bookkeeping (neighbour ranks and shared control points).
 * @param dat Nodal data array that stores all components in contiguous blocks.
 * @param offsets Pointer to `ndof` component offsets inside `dat`, such as `0`, `nuc`, `nvc`, or `nwc`.
 * @param ndof Number of components exchanged in one communication step.
 * @note Send/receive buffers are static and reused across calls, so this function is not reentrant or thread-safe.
 */
template <typename T>
void NodeVarCommImpl(const Partition &part, std::vector<T> &dat, const int *offsets, int ndof) {

    static std::vector<T> sdn, rvn; // reused send/recv buffers
    sdn.resize(static_cast<std::size_t>(part.isubc) * ndof);
    rvn.resize(static_cast<std::size_t>(part.isubc) * ndof);

    int io = 0;
    for (int i = 0; i < part.isb; i++) {
        for (int j = io; j < io + part.nsbc[i]; j++) {
            for (int d = 0; d < ndof; d++) { sdn[j * ndof + d] = dat[part.nsubc[j] + offsets[d]]; }
        }
        io += part.nsbc[i];
    }

    PeerSendrecv(part.naid, part.nsbc, part.nsbc, ndof, sdn, rvn);

    io = 0;
    for (int i = 0; i < part.isb; i++) {
        for (int j = io; j < io + part.nsbc[i]; j++) {
            for (int d = 0; d < ndof; d++) { dat[part.nsubc[j] + offsets[d]] += rvn[j * ndof + d]; }
        }
        io += part.nsbc[i];
    }

    return;
}

} // namespace internal

/**
 * @brief Exchange one nodal field component on overlap control points and add received values back into `dat`.
 * @param part Rank overlap bookkeeping (neighbour ranks and shared control points).
 * @param dat Nodal data array that stores all components in contiguous blocks.
 * @param offset Component offset inside `dat`, such as `0`, `nuc`, `nvc`, or `nwc`.
 */
template <typename T> void NodeVarComm(const Partition &part, std::vector<T> &dat, int offset) {

    internal::NodeVarCommImpl(part, dat, &offset, 1);

    return;
}

/**
 * @brief Exchange multiple nodal field components on overlap control points and add received values back into `dat`.
 * @param part Rank overlap bookkeeping (neighbour ranks and shared control points).
 * @param dat Nodal data array that stores all components in contiguous blocks.
 * @param offsets Component offsets inside `dat` for the fields packed into one communication step.
 */
template <typename T> void NodeVarComm(const Partition &part, std::vector<T> &dat, const std::vector<int> &offsets) {

    internal::NodeVarCommImpl(part, dat, offsets.data(), static_cast<int>(offsets.size()));

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

    std::vector<int> comm_ranks; // current particle peers
    std::vector<int> nsp, nrp, ofp_id;
    std::vector<int> idmp, idnsp;

    /**
     * @brief Send/receive a scalar or array-valued particle variable to/from neighbor ranks.
     * @param bufs Send buffer packed by neighbor.
     * @param bufr Receive buffer packed by neighbor.
     * @param idm  Number of values per particle.
     */
    template <typename E> void PointVarSendrecv(std::vector<E> &bufs, std::vector<E> &bufr, const int idm) {

        internal::PeerSendrecv(this->comm_ranks, this->nsp, this->nrp, idm, bufs, bufr);

        return;
    }

    /**
     * @brief Migrate a scalar or array-valued particle variable across rank boundaries and rebuild the local vector.
     * @param point_num  New local particle count after communication.
     * @param point_var  Particle variable vector (resized and repopulated in-place).
     */
    template <typename E> void PointVarComm(int point_num, std::vector<E> &point_var) {

        std::vector<E> bufs(this->nmps), bufr(this->nrps), tmp(this->nnmp);

        int sip = 0;
        for (int i = 0; i < int(this->comm_ranks.size()); i++) {
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

    /**
     * @brief Migrate a scalar or array-valued particle variable and append inflow-particle values.
     * @param point_num  New local particle count after communication.
     * @param ifp_num    Number of inflow particles to append.
     * @param point_var  Particle variable vector (resized and repopulated in-place).
     * @param ifp_var    Inflow-particle values to append.
     */
    template <typename E>
    void PointVarComm(int point_num, int ifp_num, std::vector<E> &point_var, const std::vector<E> &ifp_var) {

        this->PointVarComm(point_num, point_var);

        if (ifp_num != 0) {
            int num = this->nnmp + this->nrps;
            for (int i = 0; i < ifp_num; i++) { point_var[num + i] = ifp_var[i]; }
        }

        return;
    }

    // ---------------------------------
};
