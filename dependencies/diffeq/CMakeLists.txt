# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    diffeq
    GIT_REPOSITORY "https://github.com/jchristopherson/diffeq"
    GIT_TAG main
)

FetchContent_MakeAvailable(diffeq)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:diffeq>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(diffeq_INCLUDE_DIR ${diffeq_BINARY_DIR}/include)
set(diffeq_INCLUDE_DIR ${diffeq_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(diffeq_LIBRARY diffeq)
set(diffeq_LIBRARY ${diffeq_LIBRARY} PARENT_SCOPE)
