include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    diffeq
    GIT_REPOSITORY "https://github.com/jchristopherson/diffeq"
    GIT_TAG main
)

FetchContent_MakeAvailable(diffeq)

set(diffeq_INCLUDE_DIR ${diffeq_BINARY_DIR}/include)
set(diffeq_INCLUDE_DIR ${diffeq_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(diffeq_LIBRARY diffeq)
set(diffeq_LIBRARY ${diffeq_LIBRARY} PARENT_SCOPE)
