include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)
add_subdirectory(fortran_csv_module)

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

# Coulomb Fitting Example
add_executable(
    coulomb_fit_example
    coulomb_fit_example.f90
)
link_library(coulomb_fit_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(coulomb_fit_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(coulomb_fit_example ${fortran-csv-module_LIBRARY} ${fortran-csv-module_INCLUDE_DIR})

# Lu-Gre Fitting Example
add_executable(
    lugre_fit_example
    lugre_fit_example.f90
)
link_library(lugre_fit_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(lugre_fit_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(lugre_fit_example ${fortran-csv-module_LIBRARY} ${fortran-csv-module_INCLUDE_DIR})

# Maxwell Example
add_executable(
    maxwell_example
    maxwell_example.f90
)
link_library(maxwell_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(maxwell_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Maxwell Fitting Example
add_executable(
    maxwell_fit_example
    maxwell_fit_example.f90
)
link_library(maxwell_fit_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(maxwell_fit_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(maxwell_fit_example ${fortran-csv-module_LIBRARY} ${fortran-csv-module_INCLUDE_DIR})

# GMSM Fitting Example
add_executable(
    gmsm_fit_example
    gmsm_fit_example.f90
)
link_library(gmsm_fit_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(gmsm_fit_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(gmsm_fit_example ${fortran-csv-module_LIBRARY} ${fortran-csv-module_INCLUDE_DIR})

# Stribeck Example
add_executable(
    stribeck_example
    stribeck_example.f90
)
link_library(stribeck_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(stribeck_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
