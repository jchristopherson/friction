! MIT License
!
! Copyright (c) 2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
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