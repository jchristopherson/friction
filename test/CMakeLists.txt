include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

add_subdirectory(fortran_test_helper)

set(friction_test_sources
    friction_test.f90
)

add_executable(friction_tests ${friction_test_sources})
link_library(friction_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(friction_tests ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME friction_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:friction_tests>
)