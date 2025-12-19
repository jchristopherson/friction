# Get Dependencies
add_subdirectory(fortran_csv_module)
include(FetchContent)
find_package(ferror QUIET)
find_package(fplot QUIET)
find_package(collections QUIET)
find_package(diffeq QUIET)
find_package(linalg QUIET)

if (NOT ferror_FOUND)
    FetchContent_Declare(
        ferror
        GIT_TAG "origin/master"
        GIT_REPOSITORY https://github.com/jchristopherson/ferror
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
else()
    set(FERROR_LIBRARIES ferror::ferror)
endif()

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

if (NOT collections_FOUND)
    message(STATUS "COLLECTIONS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        collections
        GIT_REPOSITORY "https://github.com/jchristopherson/collections"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(collections)
    set(COLLECTIONS_LIBRARIES collections)
else()
    set(COLLECTIONS_LIBRARIES collections::collections)
endif()

if (NOT diffeq_FOUND)
    message(STATUS "DIFFEQ could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        diffeq
        GIT_REPOSITORY "https://github.com/jchristopherson/diffeq"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(diffeq)
    set(DIFFEQ_LIBRARIES diffeq)
else()
    set(DIFFEQ_LIBRARIES diffeq::diffeq)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
else()
    set(LINALG_LIBRARIES linalg::linalg)
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
    ${DIFFEQ_LIBRARIES}
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
