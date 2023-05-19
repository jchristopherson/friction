# Get FERROR
add_subdirectory(ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)

# Get DIFFEQ
add_subdirectory(diffeq)
set(diffeq_LIBRARY ${diffeq_LIBRARY} PARENT_SCOPE)
set(diffeq_INCLUDE_DIR ${diffeq_INCLUDE_DIR} PARENT_SCOPE)