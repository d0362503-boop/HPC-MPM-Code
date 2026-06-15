# -----------------------------------------------------------------------------
# Solid Solver Source - Main solver implementations for solid problems
# -----------------------------------------------------------------------------

# Keep all candidate subdirectories visible here.
# Uncomment exactly one case. This must stay consistent with SOLID_METHOD in
# cmake/options.cmake because module/solid still uses that option to select the
# matching solid mechanics implementation.
#
# add_subdirectory(explicit)
add_subdirectory(implicit)
