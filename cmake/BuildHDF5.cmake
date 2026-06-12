# -----------------------------------------------------------------------------
# Build HDF5 from Ext/hdf5-hdf5_1.14.5.tar.gz into external/hdf5
# This script is idempotent: if external/hdf5 already exists, it does nothing.
# -----------------------------------------------------------------------------

set(HDF5_TARBALL "${SRC_DIR}/Ext/hdf5-hdf5_1.14.5.tar.gz")
set(HDF5_INSTALL_DIR "${SRC_DIR}/external/hdf5")
set(MARKER "${BIN_DIR}/hdf5_is_built")

if(EXISTS "${HDF5_INSTALL_DIR}/include/hdf5.h")
    message(STATUS "HDF5 already installed at ${HDF5_INSTALL_DIR}; skipping build.")
    file(TOUCH "${MARKER}")
    return()
endif()

if(NOT EXISTS ${HDF5_TARBALL})
    message(FATAL_ERROR "HDF5 tarball not found: ${HDF5_TARBALL}")
endif()

set(TEMP_HDF5_DIR "${BIN_DIR}/hdf5")

message(STATUS "Building HDF5 from ${HDF5_TARBALL} ...")

file(MAKE_DIRECTORY ${TEMP_HDF5_DIR})
file(REMOVE_RECURSE ${TEMP_HDF5_DIR}/src)
file(MAKE_DIRECTORY ${TEMP_HDF5_DIR}/src)

execute_process(
    COMMAND tar zxvf ${HDF5_TARBALL} -C ${TEMP_HDF5_DIR}/src --strip-components=1
    RESULT_VARIABLE tar_result)
if(NOT tar_result EQUAL 0)
    message(FATAL_ERROR "Failed to extract HDF5 tarball")
endif()

file(MAKE_DIRECTORY ${TEMP_HDF5_DIR}/build)

execute_process(
    COMMAND
        ${CMAKE_COMMAND} -S ${TEMP_HDF5_DIR}/src
        -B ${TEMP_HDF5_DIR}/build
        -DHDF5_ENABLE_SZIP_SUPPORT=OFF
        -DHDF5_BUILD_UTILS=OFF
        -DBUILD_SHARED_LIBS=OFF
        -DHDF5_BUILD_TESTING=OFF
        -DHDF5_BUILD_EXAMPLES=OFF
        -DHDF5_ENABLE_PARALLEL=${HDF5_ENABLE_PARALLEL}
        -DHDF5_BUILD_TOOLS=${HDF5_BUILD_TOOLS}
        -DHDF5_BUILD_HL_LIB=${HDF5_BUILD_HL_LIB}
        -DHDF5_BUILD_CPP_LIB=OFF
        -DH5EX_BUILD_EXAMPLES=OFF
        -DH5EX_BUILD_TESTING=OFF
        -DH5EX_BUILD_HL_LIB=OFF
        -DTEST_SHELL_SCRIPTS=OFF
        -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    RESULT_VARIABLE configure_result)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "HDF5 configure failed")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} --build ${TEMP_HDF5_DIR}/build -- -j8
    RESULT_VARIABLE build_result)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "HDF5 build failed")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} --install ${TEMP_HDF5_DIR}/build
    RESULT_VARIABLE install_result)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "HDF5 install failed")
endif()

file(TOUCH "${MARKER}")
message(STATUS "HDF5 installed at ${HDF5_INSTALL_DIR}")
