include(FetchContent)
FetchContent_Declare(
    fplot
    GIT_TAG "origin/master"
    GIT_REPOSITORY https://github.com/jchristopherson/fplot
)
FetchContent_MakeAvailable(fplot)
set(fplot_INCLUDE_DIR ${fplot_BINARY_DIR}/include)
set(fplot_LIBRARY fplot)

set(fplot_INCLUDE_DIR ${fplot_INCLUDE_DIR} PARENT_SCOPE)
set(fplot_LIBRARY ${fplot_LIBRARY} PARENT_SCOPE)