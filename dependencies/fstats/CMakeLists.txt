# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fstats
    GIT_REPOSITORY "https://github.com/jchristopherson/fstats"
    GIT_TAG main
)

FetchContent_MakeAvailable(fstats)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:fstats>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(fstats_INCLUDE_DIR ${fstats_BINARY_DIR}/include)
set(fstats_INCLUDE_DIR ${fstats_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fstats_LIBRARY fstats)
set(fstats_LIBRARY ${fstats_LIBRARY} PARENT_SCOPE)
