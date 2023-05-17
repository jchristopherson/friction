include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# RK4 Example
add_executable(
    rk4_example
    rk4_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(rk4_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rk4_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Adams Method - Fixed - Example
add_executable(
    adams_fixed_example
    adams_fixed_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(adams_fixed_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(adams_fixed_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# DPRK45 Example
add_executable(
    dprk45_example
    dprk45_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(dprk45_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(dprk45_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# DPRK45 Dense Output Example
add_executable(
    dprk45_dense_example
    dprk45_dense_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(dprk45_dense_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(dprk45_dense_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# BSRK32 Example
add_executable(
    bsrk32_example
    bsrk32_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(bsrk32_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(bsrk32_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# FRF Example
add_executable(
    frf_example
    frf_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(frf_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(frf_example ${spectrum_LIBRARY} ${spectrum_INCLUDE_DIR})

# FRF Sweep Example
add_executable(
    frf_sweep_example
    frf_sweep_example.f90
)
link_library(frf_sweep_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(frf_sweep_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})