# -----------------------------------------------------------------------------
# Data partition tools
# -----------------------------------------------------------------------------

set(g_selected_divide_case_count 0)
if(USE_DIVIDE_FLUID)
    math(EXPR g_selected_divide_case_count "${g_selected_divide_case_count} + 1")
endif()
if(USE_DIVIDE_SOLID)
    math(EXPR g_selected_divide_case_count "${g_selected_divide_case_count} + 1")
endif()
if(USE_DIVIDE_FSI)
    math(EXPR g_selected_divide_case_count "${g_selected_divide_case_count} + 1")
endif()

if(NOT g_selected_divide_case_count EQUAL 1)
    message(FATAL_ERROR
        "Exactly one divide case must be enabled. "
        "Set exactly one of USE_DIVIDE_FLUID, USE_DIVIDE_SOLID, USE_DIVIDE_FSI to ON.")
endif()

add_library(mpm_divide_common OBJECT
    data_partitioner.cpp
)

target_link_libraries(mpm_divide_common PUBLIC mpm_modules)

if(NOT USE_DIVIDE_FLUID)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/divide/fluid)
endif()
if(NOT USE_DIVIDE_SOLID)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/divide/solid)
endif()
if(NOT USE_DIVIDE_FSI)
    file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/data/divide/fsi)
endif()

if(USE_DIVIDE_FLUID)
    add_subdirectory(divide_fluid fluid)
endif()

if(USE_DIVIDE_SOLID)
    add_subdirectory(divide_solid solid)
endif()

if(USE_DIVIDE_FSI)
    add_subdirectory(divide_fsi fsi)
endif()
