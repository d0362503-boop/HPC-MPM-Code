// HDF5 smoke test for build configuration verification.
// Compile and run with: make test-hdf5

#include <hdf5.h>
#include <iostream>

int main(int argc, char** argv) {
  std::cout << "HDF5 version: " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "."
            << H5_VERS_RELEASE << std::endl;

  hid_t file_id =
      H5Fcreate("test_hdf5.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "Failed to create HDF5 file" << std::endl;
    return 1;
  }

  hsize_t dims[1] = {4};
  hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
  hid_t dataset_id = H5Dcreate2(file_id, "test_dataset", H5T_NATIVE_DOUBLE,
                                dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  double data[4] = {1.0, 2.0, 3.0, 4.0};
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Fclose(file_id);

  std::cout << "HDF5 test file created successfully" << std::endl;
  return 0;
}
