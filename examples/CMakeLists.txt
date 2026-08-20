# Get Dependencies
add_subdirectory(fortran_csv_module)
include(FetchContent)
find_package(fplot QUIET)

if (NOT fplot_FOUND)
    FetchContent_Declare(
        fplot
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/fplot
    )
    FetchContent_MakeAvailable(fplot)
    set(FPLOT_LIBRARIES fplot)
else()
    set(FPLOT_LIBRARIES fplot::fplot)
endif()

# Coulomb Example
add_executable(
    coulomb_example
    coulomb_example.f90
)
target_link_libraries(
    coulomb_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Lu-Gre Example
add_executable(
    lugre_example
    lugre_example.f90
)
target_link_libraries(
    lugre_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Maxwell Example
add_executable(
    maxwell_example
    maxwell_example.f90
)
target_link_libraries(
    maxwell_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Maxwell Fitting Example
add_executable(
    maxwell_fit_example
    maxwell_fit_example.f90
)
target_link_libraries(
    maxwell_fit_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${fortran-csv-module_LIBRARY}
)
target_include_directories(maxwell_fit_example PUBLIC ${fortran-csv-module_INCLUDE_DIR})

# Stribeck Example
add_executable(
    stribeck_example
    stribeck_example.f90
)
target_link_libraries(
    stribeck_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Modified Stribeck Example
add_executable(
    modified_stribeck_example
    modified_stribeck_example.f90
)
target_link_libraries(
    modified_stribeck_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
)

# Modified Stribeck Fitting Example
add_executable(
    modified_stribeck_fit_example
    modified_stribeck_fit_example.f90
)
target_link_libraries(
    modified_stribeck_fit_example
    ${PROJECT_NAME}
    ${FPLOT_LIBRARIES}
    ${fortran-csv-module_LIBRARY}
)
target_include_directories(modified_stribeck_fit_example PUBLIC ${fortran-csv-module_INCLUDE_DIR})

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET coulomb_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:coulomb_example> 
        $<TARGET_FILE_DIR:coulomb_example>
        COMMAND_EXPAND_LISTS
    )
endif()

configure_file(
    ${CMAKE_SOURCE_DIR}/examples/data/friction_data_1.csv
    ${CMAKE_BINARY_DIR}/examples/friction_data_1.csv
    COPYONLY
)