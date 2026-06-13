# -----------------------------------------------------------------------------
# Data tools
# -----------------------------------------------------------------------------

if(USE_DATA_FLUID OR USE_DATA_SOLID OR USE_DATA_FSI)
    add_subdirectory(generate)
endif()
