# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FRICTION_SOURCES
    ${dir}/friction.f90
    ${dir}/friction_core.f90
    ${dir}/friction_errors.f90
    ${dir}/friction_coulomb.f90
    ${dir}/friction_lugre.f90
    ${dir}/friction_maxwell.f90
    ${dir}/friction_gmsm.f90
    ${dir}/friction_stribeck.f90
    ${dir}/friction_modified_stribeck.f90
    ${dir}/fitpack/fitpack_core.f90
    ${dir}/fitpack/fitpack.f90
)

# Dependencies
include(FetchContent)
find_package(ferror QUIET)
find_package(fstats QUIET)
find_package(collections QUIET)
find_package(diffeq QUIET)
find_package(linalg QUIET)

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

if (NOT fstats_FOUND)
    message(STATUS "FSTATS could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        fstats
        GIT_REPOSITORY "https://github.com/jchristopherson/fstats"
        GIT_TAG main
    )
    FetchContent_MakeAvailable(fstats)
    set(FSTATS_LIBRARIES fstats)
    set(BUILD_FSTATS TRUE)
else()
    set(FSTATS_LIBRARIES fstats::fstats)
    set(BUILD_FSTATS FALSE)
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
    set(BUILD_DIFFEQ TRUE)
else()
    set(DIFFEQ_LIBRARIES diffeq::diffeq)
    set(BUILD_DIFFEQ FALSE)
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
    set(BUILD_COLLECTIONS TRUE)
else()
    set(COLLECTIONS_LIBRARIES collections::collections)
    set(BUILD_COLLECTIONS FALSE)
endif()

# Build
add_library(${PROJECT_NAME} ${FRICTION_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${FERROR_LIBRARIES}
    ${FSTATS_LIBRARIES}
)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    ${LINALG_LIBRARIES}
    ${DIFFEQ_LIBRARIES}
    ${COLLECTIONS_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(FRICTION_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${FRICTION_MOD_DIR}")
    file(MAKE_DIRECTORY "${FRICTION_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${FRICTION_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${FRICTION_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

# Install
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${FRICTION_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})

if (${BUILD_FERROR} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FERROR_LIBRARIES})
endif()

if (${BUILD_LINALG} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${LINALG_LIBRARIES})
endif()

if (${BUILD_FSTATS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${FSTATS_LIBRARIES})
endif()

if (${BUILD_DIFFEQ} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${DIFFEQ_LIBRARIES})
endif()

if (${BUILD_COLLECTIONS} AND ${BUILD_SHARED_LIBS})
    install(IMPORTED_RUNTIME_ARTIFACTS ${COLLECTIONS_LIBRARIES})
endif()
