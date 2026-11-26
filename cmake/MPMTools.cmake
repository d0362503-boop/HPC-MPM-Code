# ----------------------------------------------------------------------------
# Common CMake helpers for MPM-Code
# ----------------------------------------------------------------------------

include_guard(GLOBAL)

# ----------------------------------------------------------------------------
# Add a data-generation executable (makinput_*).
#
# Usage:
#   mpm_add_data_tool(
#       TARGET <name>
#       SOURCES <source>...
#       [COMMON_TARGET <common_object_target>]
#   )
#
# The executable is linked against the common object target, which is expected
# to PUBLIC-link mpm_cxx_compat and any other dependencies it needs.
# ----------------------------------------------------------------------------
function(mpm_add_data_tool)
    set(options "")
    set(one_value_args TARGET COMMON_TARGET)
    set(multi_value_args SOURCES)
    cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(NOT ARG_TARGET)
        message(FATAL_ERROR "mpm_add_data_tool called without TARGET")
    endif()

    if(NOT ARG_COMMON_TARGET)
        set(ARG_COMMON_TARGET mpm_generate_common)
    endif()

    add_executable(${ARG_TARGET}
        ${ARG_SOURCES}
        $<TARGET_OBJECTS:${ARG_COMMON_TARGET}>
    )

    # The common object target is an OBJECT library; it can propagate include
    # directories and compile definitions, but the final executable must still
    # explicitly link mpm_modules so that the module object files are pulled in.
    target_link_libraries(${ARG_TARGET} PRIVATE ${ARG_COMMON_TARGET} mpm_modules)

    set_target_properties(${ARG_TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

# ----------------------------------------------------------------------------
# Add a data-partition executable (makdivide_*).
#
# Usage:
#   mpm_add_divide_tool(
#       TARGET <name>
#       SOURCES <source>...
#       [COMMON_TARGET <common_object_target>]
#   )
# ----------------------------------------------------------------------------
function(mpm_add_divide_tool)
    set(options "")
    set(one_value_args TARGET COMMON_TARGET)
    set(multi_value_args SOURCES)
    cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(NOT ARG_TARGET)
        message(FATAL_ERROR "mpm_add_divide_tool called without TARGET")
    endif()

    if(NOT ARG_COMMON_TARGET)
        set(ARG_COMMON_TARGET mpm_divide_common)
    endif()

    add_executable(${ARG_TARGET}
        ${ARG_SOURCES}
        $<TARGET_OBJECTS:${ARG_COMMON_TARGET}>
    )

    # The common object target is an OBJECT library; it can propagate include
    # directories and compile definitions, but the final executable must still
    # explicitly link mpm_modules so that the module object files are pulled in.
    target_link_libraries(${ARG_TARGET} PRIVATE ${ARG_COMMON_TARGET} mpm_modules)

    set_target_properties(${ARG_TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endfunction()

# ----------------------------------------------------------------------------
# Add one solver driver executable built from all .cpp files in the calling
# directory. The executable is written to the driver's own build directory,
# for example build/work/src_fsi/MPM_MPM/MPM_MPM_FSI, next to its file.dat.
#
# Usage:
#   mpm_add_driver(<target>)
# ----------------------------------------------------------------------------
function(mpm_add_driver target)
    file(GLOB driver_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    add_executable(${target} ${driver_sources})
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/module)
    target_link_libraries(${target} PRIVATE mpm_modules)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_property(GLOBAL APPEND PROPERTY MPM_DRIVERS ${target})
endfunction()

# ----------------------------------------------------------------------------
# Copy a runtime parameter file into the build tree only if it is not there
# yet, so edits made under build/ survive later re-configures. Delete the
# build copy and re-configure to get the source version again.
#
# Usage:
#   mpm_copy_input_once(<file>)
# ----------------------------------------------------------------------------
function(mpm_copy_input_once file)
    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${file})
        configure_file(${file} ${file} COPYONLY)
    endif()
endfunction()
