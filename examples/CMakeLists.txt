include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# Coulomb Example
add_executable(
    coulomb_example
    coulomb_example.f90
)
link_library(coulomb_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(coulomb_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})