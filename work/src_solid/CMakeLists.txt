# -----------------------------------------------------------------------------
# Solid Solver Source - Main solver implementations for solid problems
# -----------------------------------------------------------------------------

# Keep all candidate subdirectories visible here.
# Uncomment exactly one case. The matching implementation is already compiled
# into mpm_modules, so this choice only selects which driver main is linked.
#
# add_subdirectory(explicit)
add_subdirectory(implicit)
