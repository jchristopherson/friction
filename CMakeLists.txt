# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    friction
    LANGUAGES Fortran
    VERSION 1.1.1
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)

# Installation Instructions
install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build FRICTION tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Examples
option(BUILD_FRICTION_EXAMPLES "Build FRICTION examples")
message(STATUS "Build FRICTION examples: ${BUILD_FRICTION_EXAMPLES}")
if (BUILD_FRICTION_EXAMPLES)
    add_subdirectory(examples)
endif()