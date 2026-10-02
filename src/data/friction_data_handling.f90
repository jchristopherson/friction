module friction_data_handling
    use iso_fortran_env
    implicit none
    private
    public :: friction_data

    type friction_data
        real(real64), allocatable, dimension(:) :: time
            !! An N-element array containing the time points at which the 
            !! data was sampled.
        real(real64), allocatable, dimension(:) :: position
            !! An N-element array containing the relative position of the
            !! two contacting bodies.
        real(real64), allocatable, dimension(:) :: velocity
            !! An N-element array containing the relative velocity of the
            !! two contacting bodies.
        real(real64), allocatable, dimension(:) :: normal_force
            !! An N-element array containing the normal force between the
            !! two contacting bodies.
        real(real64), allocatable, dimension(:) :: friction_force
            !! An N-element array containing the friction force between the
            !! two contacting bodies.
    end type
end module