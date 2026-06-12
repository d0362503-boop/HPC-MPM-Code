# ----------------------------------------------------------------------------
# Solver source selection
#
# This directory holds the top-level solver entry points. Which ones are
# compiled into the final executable is controlled by the options in the
# root CMakeLists.txt.
# ----------------------------------------------------------------------------

if(USE_SRC_FLUID)
    add_subdirectory(src_fluid)
endif()

if(USE_SRC_SOLID)
    add_subdirectory(src_solid)
endif()

if(USE_SRC_FSI)
    add_subdirectory(src_fsi)
endif()
