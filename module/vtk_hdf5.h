#pragma once

#ifdef HAVE_HDF5

#include <hdf5.h>
#include <mpi.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
// VTK HDF5 2.0 writer class
//
// All members are public to match the requested style.  The class wraps an
// HDF5 file and exposes helpers for writing VTK UnstructuredGrid data using
// MPI-IO collective I/O.
//
// Output files (.vtkhdf) can be opened directly in ParaView 5.10+.
// For multi-rank runs, all ranks write to a single shared file.
// ----------------------------------------------------------------------------

namespace vtkhdf {

// Type mapping from C++ types to HDF5 file/memory types.
template <typename T>
struct HDF5TypeMap;

template <>
struct HDF5TypeMap<int> {
    static hid_t fileType() { return H5T_STD_I32LE; }
    static hid_t memType() { return H5T_NATIVE_INT; }
};

template <>
struct HDF5TypeMap<long long> {
    static hid_t fileType() { return H5T_STD_I64LE; }
    static hid_t memType() { return H5T_NATIVE_LLONG; }
};

template <>
struct HDF5TypeMap<double> {
    static hid_t fileType() { return H5T_IEEE_F64LE; }
    static hid_t memType() { return H5T_NATIVE_DOUBLE; }
};

template <>
struct HDF5TypeMap<unsigned char> {
    static hid_t fileType() { return H5T_STD_U8LE; }
    static hid_t memType() { return H5T_NATIVE_UCHAR; }
};

class VTKHDFWriter {
  public:
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank = 0;
    int n_ranks = 1;

    hid_t file = -1;
    hid_t vtkhdf_group = -1;
    hid_t pointdata_group = -1;
    hid_t celldata_group = -1;

    VTKHDFWriter() = default;

    explicit VTKHDFWriter(const std::string& filename, MPI_Comm comm_arg = MPI_COMM_WORLD) { Open(filename, comm_arg); }

    ~VTKHDFWriter() { Close(); }

    // Disable copy; allow move.
    VTKHDFWriter(const VTKHDFWriter&) = delete;
    VTKHDFWriter& operator=(const VTKHDFWriter&) = delete;
    VTKHDFWriter(VTKHDFWriter&& other) noexcept { *this = std::move(other); }
    VTKHDFWriter& operator=(VTKHDFWriter&& other) noexcept {
        if (this != &other) {
            Close();
            comm = other.comm;
            my_rank = other.my_rank;
            n_ranks = other.n_ranks;
            file = other.file;
            vtkhdf_group = other.vtkhdf_group;
            pointdata_group = other.pointdata_group;
            celldata_group = other.celldata_group;
            other.file = -1;
            other.vtkhdf_group = -1;
            other.pointdata_group = -1;
            other.celldata_group = -1;
        }
        return *this;
    }

    void Open(const std::string& filename, MPI_Comm comm_arg = MPI_COMM_WORLD) {
        comm = comm_arg;
        MPI_Comm_rank(comm, &my_rank);
        MPI_Comm_size(comm, &n_ranks);

        // Create parent directories if they do not exist.
        std::filesystem::path p(filename);
        if (p.has_parent_path()) {
            std::filesystem::create_directories(p.parent_path());
        }

        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
        H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
        H5Pset_all_coll_metadata_ops(fapl, false);
        H5Pset_coll_metadata_write(fapl, false);
        file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5Pclose(fapl);
        if (file < 0) throw std::runtime_error("H5Fcreate failed for " + filename);
    }

    void Close() {
        if (celldata_group >= 0) {
            H5Gclose(celldata_group);
            celldata_group = -1;
        }
        if (pointdata_group >= 0) {
            H5Gclose(pointdata_group);
            pointdata_group = -1;
        }
        if (vtkhdf_group >= 0) {
            H5Gclose(vtkhdf_group);
            vtkhdf_group = -1;
        }
        if (file >= 0) {
            H5Fclose(file);
            file = -1;
        }
    }

    // ------------------------------------------------------------------------
    // Global offset helper
    // ------------------------------------------------------------------------
    static void ComputeGlobalInfo(hsize_t local_count, hsize_t& global_offset, hsize_t& total_count,
                                  MPI_Comm comm = MPI_COMM_WORLD) {
        global_offset = 0;
        MPI_Exscan(&local_count, &global_offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        total_count = local_count;
        MPI_Allreduce(MPI_IN_PLACE, &total_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    }

    // ------------------------------------------------------------------------
    // Attribute helpers
    // ------------------------------------------------------------------------
    template <typename T>
    static void WriteAttribute(hid_t loc, const char* name, int n, const T* data) {
        hsize_t dim = static_cast<hsize_t>(n);
        hid_t attr_space = H5Screate_simple(1, &dim, nullptr);
        hid_t h5type = HDF5TypeMap<T>::memType();
        hid_t attr = H5Acreate(loc, name, h5type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        if (attr < 0) throw std::runtime_error(std::string("H5Acreate failed for attribute ") + name);
        H5Awrite(attr, h5type, data);
        H5Aclose(attr);
        H5Sclose(attr_space);
    }

    static void WriteAttributeString(hid_t loc, const char* name, const char* str) {
        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, strlen(str));
        hid_t attr_space = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate(loc, name, str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        if (attr < 0) throw std::runtime_error(std::string("H5Acreate failed for attribute ") + name);
        H5Awrite(attr, str_type, str);
        H5Aclose(attr);
        H5Sclose(attr_space);
        H5Tclose(str_type);
    }

    // ------------------------------------------------------------------------
    // Dataset creation / collective write helpers
    // ------------------------------------------------------------------------
    static hid_t CreateDataset1D(hid_t loc, const char* name, hsize_t total_size, hid_t h5type, hid_t& out_filespace) {
        out_filespace = H5Screate_simple(1, &total_size, nullptr);
        hid_t dset = H5Dcreate(loc, name, h5type, out_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dset < 0) throw std::runtime_error(std::string("H5Dcreate failed for ") + name);
        return dset;
    }

    void WriteArrayCollective(hid_t dset, hid_t filespace, const void* data, hsize_t local_count,
                              hsize_t global_offset, hid_t h5type) {
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

        if (local_count > 0) {
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &global_offset, nullptr, &local_count, nullptr);
            hid_t memspace = H5Screate_simple(1, &local_count, nullptr);
            H5Dwrite(dset, h5type, memspace, filespace, dxpl, data);
            H5Sclose(memspace);
        } else {
            H5Sselect_none(filespace);
            hsize_t zero = 0;
            hid_t memspace = H5Screate_simple(1, &zero, nullptr);
            H5Dwrite(dset, h5type, memspace, filespace, dxpl, nullptr);
            H5Sclose(memspace);
        }

        H5Pclose(dxpl);
    }

    void CreateAndWriteScalar(hid_t loc, const char* name, hsize_t total_size, hid_t filetype, hid_t memtype,
                              const void* data, hsize_t local_count, hsize_t global_offset) {
        hid_t filespace;
        hid_t dset = CreateDataset1D(loc, name, total_size, filetype, filespace);
        WriteArrayCollective(dset, filespace, data, local_count, global_offset, memtype);
        H5Sclose(filespace);
        H5Dclose(dset);
    }

    void CreateAndWriteVector(hid_t loc, const char* name, hsize_t npts, int ncomp, hid_t filetype, hid_t memtype,
                              const void* data, hsize_t local_npts, hsize_t global_offset) {
        hsize_t dims[2] = {npts, static_cast<hsize_t>(ncomp)};
        hid_t filespace = H5Screate_simple(2, dims, nullptr);
        hid_t dset = H5Dcreate(loc, name, filetype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dset < 0) throw std::runtime_error(std::string("H5Dcreate failed for ") + name);

        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

        if (local_npts > 0) {
            hsize_t start[2] = {global_offset, 0};
            hsize_t count[2] = {local_npts, static_cast<hsize_t>(ncomp)};
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

            hsize_t mem_dims[2] = {local_npts, static_cast<hsize_t>(ncomp)};
            hid_t memspace = H5Screate_simple(2, mem_dims, nullptr);
            H5Dwrite(dset, memtype, memspace, filespace, dxpl, data);
            H5Sclose(memspace);
        } else {
            H5Sselect_none(filespace);
            hsize_t zero = 0;
            hid_t memspace = H5Screate_simple(1, &zero, nullptr);
            H5Dwrite(dset, memtype, memspace, filespace, dxpl, nullptr);
            H5Sclose(memspace);
        }

        H5Pclose(dxpl);
        H5Sclose(filespace);
        H5Dclose(dset);
    }

    // ------------------------------------------------------------------------
    // UnstructuredGrid helpers
    // ------------------------------------------------------------------------
    void CreateUnstructuredGridGroup(hsize_t total_npts, hsize_t total_ncells, hsize_t total_connectivity) {
        vtkhdf_group = H5Gcreate(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (vtkhdf_group < 0) throw std::runtime_error("H5Gcreate VTKHDF failed");

        int version[2] = {2, 0};
        WriteAttribute(vtkhdf_group, "Version", 2, version);
        WriteAttributeString(vtkhdf_group, "Type", "UnstructuredGrid");

        const char* count_names[3] = {"NumberOfPoints", "NumberOfCells", "NumberOfConnectivityIds"};
        long long values[3] = {static_cast<long long>(total_npts), static_cast<long long>(total_ncells),
                               static_cast<long long>(total_connectivity)};

        for (int i = 0; i < 3; ++i) {
            hid_t filespace;
            hid_t dset = CreateDataset1D(vtkhdf_group, count_names[i], 1, H5T_STD_I64LE, filespace);
            hsize_t count = (my_rank == 0) ? 1 : 0;
            WriteArrayCollective(dset, filespace, &values[i], count, 0, H5T_NATIVE_LLONG);
            H5Sclose(filespace);
            H5Dclose(dset);
        }
    }

    void SetTime(double time) { WriteAttribute(vtkhdf_group, "Time", 1, &time); }

    void CreatePointDataGroup() {
        pointdata_group = H5Gcreate(vtkhdf_group, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    void CreateCellDataGroup() {
        celldata_group = H5Gcreate(vtkhdf_group, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    void WritePoints(hsize_t total_npts, hsize_t local_npts, hsize_t global_offset,
                     const std::vector<std::array<double, 3>>& points) {
        CreateAndWriteVector(vtkhdf_group, "Points", total_npts, 3, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, points.data(),
                             local_npts, global_offset);
    }

    void WriteConnectivity(hsize_t total_size, hsize_t local_count, hsize_t global_offset,
                           const std::vector<long long>& connectivity) {
        CreateAndWriteScalar(vtkhdf_group, "Connectivity", total_size, H5T_STD_I64LE, H5T_NATIVE_LLONG,
                             connectivity.data(), local_count, global_offset);
    }

    void WriteOffsets(hsize_t total_size, hsize_t local_count, hsize_t global_offset,
                      const std::vector<long long>& offsets) {
        CreateAndWriteScalar(vtkhdf_group, "Offsets", total_size, H5T_STD_I64LE, H5T_NATIVE_LLONG, offsets.data(),
                             local_count, global_offset);
    }

    void WriteTypes(hsize_t total_size, hsize_t local_count, hsize_t global_offset,
                    const std::vector<unsigned char>& types) {
        CreateAndWriteScalar(vtkhdf_group, "Types", total_size, H5T_STD_U8LE, H5T_NATIVE_UCHAR, types.data(),
                             local_count, global_offset);
    }

    void WriteUnstructuredGridTopology(hsize_t total_npts, hsize_t local_npts, hsize_t global_offset,
                                       const std::vector<std::array<double, 3>>& points,
                                       const std::vector<long long>& connectivity,
                                       const std::vector<long long>& offsets,
                                       const std::vector<unsigned char>& types) {
        WritePoints(total_npts, local_npts, global_offset, points);
        WriteConnectivity(total_npts, local_npts, global_offset, connectivity);

        // Offsets size depends on whether this rank holds the closing boundary value.
        bool is_last_rank = (my_rank == n_ranks - 1);
        hsize_t local_offsets_count = local_npts + (is_last_rank ? 1 : 0);
        WriteOffsets(total_npts + 1, local_offsets_count, global_offset, offsets);

        WriteTypes(total_npts, local_npts, global_offset, types);
    }

    // ------------------------------------------------------------------------
    // Point / cell data helpers
    // ------------------------------------------------------------------------
    void WritePointScalar(const char* name, hsize_t total_npts, hsize_t local_npts, hsize_t global_offset,
                          const std::vector<double>& field) {
        CreateAndWriteScalar(pointdata_group, name, total_npts, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, field.data(),
                             local_npts, global_offset);
    }

    void WritePointScalarInt(const char* name, hsize_t total_npts, hsize_t local_npts, hsize_t global_offset,
                             const std::vector<int>& field) {
        CreateAndWriteScalar(pointdata_group, name, total_npts, H5T_STD_I32LE, H5T_NATIVE_INT, field.data(), local_npts,
                             global_offset);
    }

    void WritePointVector(const char* name, hsize_t total_npts, hsize_t local_npts, hsize_t global_offset,
                          const std::vector<std::array<double, 3>>& field) {
        CreateAndWriteVector(pointdata_group, name, total_npts, 3, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, field.data(),
                             local_npts, global_offset);
    }

    void WriteCellScalar(const char* name, hsize_t total_ncells, hsize_t local_ncells, hsize_t global_offset,
                         const std::vector<double>& field) {
        CreateAndWriteScalar(celldata_group, name, total_ncells, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, field.data(),
                             local_ncells, global_offset);
    }

    void WriteCellVector(const char* name, hsize_t total_ncells, hsize_t local_ncells, hsize_t global_offset,
                         const std::vector<std::array<double, 3>>& field) {
        CreateAndWriteVector(celldata_group, name, total_ncells, 3, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, field.data(),
                             local_ncells, global_offset);
    }
};

}  // namespace vtkhdf

#endif  // HAVE_HDF5
