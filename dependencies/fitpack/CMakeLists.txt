include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fitpack
    GIT_REPOSITORY "https://github.com/jchristopherson/fitpack"
    GIT_TAG main
)

FetchContent_MakeAvailable(fitpack)

set(fitpack_INCLUDE_DIR ${fitpack_BINARY_DIR}/include)
set(fitpack_INCLUDE_DIR ${fitpack_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fitpack_LIBRARY fitpack)
set(fitpack_LIBRARY ${fitpack_LIBRARY} PARENT_SCOPE)
