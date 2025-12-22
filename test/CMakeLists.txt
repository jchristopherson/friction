set(friction_test_sources
    friction_test.f90
    friction_model_tests.f90
)

# Dependencies
include(FetchContent)
find_package(fortran_test_helper QUIET)

if (NOT fortran_test_helper_FOUND)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

# Build the tests
add_executable(friction_tests ${friction_test_sources})
target_link_libraries(
    friction_tests
    friction
    ${FTH_LIBRARIES}
)
add_test(
    NAME friction_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:friction_tests>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET friction_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:friction_tests> 
        $<TARGET_FILE_DIR:friction_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
