include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fstats
    GIT_REPOSITORY "https://github.com/jchristopherson/fstats"
    GIT_TAG main
)

FetchContent_MakeAvailable(fstats)

set(fstats_INCLUDE_DIR ${fstats_BINARY_DIR}/include)
set(fstats_INCLUDE_DIR ${fstats_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fstats_LIBRARY fstats)
set(fstats_LIBRARY ${fstats_LIBRARY} PARENT_SCOPE)
