module friction_errors
    use iso_fortran_env
    implicit none

    integer(int32), parameter :: FRICTION_ARRAY_SIZE_ERROR = 100000
        !! Defines an array size error.
    integer(int32), parameter :: FRICTION_MEMORY_ERROR = 100001
        !! Defines a memory allocation error.
    integer(int32), parameter :: FRICTION_INVALID_OPERATION_ERROR = 100002
        !! Defines an error within the opration of a routine.
end module