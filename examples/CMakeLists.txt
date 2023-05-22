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

# Lu-Gre Example
add_executable(
    lugre_example
    lugre_example.f90
)
link_library(lugre_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(lugre_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(lugre_example ${diffeq_LIBRARY} ${diffeq_INCLUDE_DIR})