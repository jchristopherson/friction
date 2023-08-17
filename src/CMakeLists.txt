# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(FRICTION_SOURCES
    ${dir}/friction.f90
    ${dir}/friction_coulomb.f90
    ${dir}/friction_lugre.f90
    ${dir}/friction_fitting.f90
    ${dir}/friction_maxwell.f90
)
set(FRICTION_SOURCES ${FRICTION_SOURCES} PARENT_SCOPE)