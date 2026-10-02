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
program example
    use iso_fortran_env
    use friction
    use fplot_core
    use diffeq
    implicit none

    ! Model Parameters
    real(real64), parameter :: mu_c = 0.4d0
    real(real64), parameter :: mu_s = 0.8d0
    real(real64), parameter :: bv = 0.15d0
    real(real64), parameter :: vs = 0.5d0
    real(real64), parameter :: ki(4) = [12.0d6, 2.5d6, 7.1d5, 7.9d5]
    real(real64), parameter :: bi(4) = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: vi(4) = [0.067d0, 0.056d0, 0.047d0, 0.83d0]
    real(real64), parameter :: stiffness = 1.0d6
    real(real64), parameter :: mass = 1.5d2

    ! Misc Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Excitation Parameters
    real(real64), parameter :: nrm = 1.0d3
    real(real64), parameter :: amp = 1.0d-1
    real(real64), parameter :: freq = 5.0d0
    real(real64), parameter :: omega = 2.0d0 * pi * freq

    ! Local Variables
    integer(int32) :: i, neqn, n
    logical :: check
    type(generalized_maxwell_slip_model) :: mdl
    type(runge_kutta_45) :: integrator
    type(ode_container) :: sys
    real(real64), allocatable, dimension(:) :: ic, F
    real(real64), allocatable, dimension(:,:) :: sol
    type(plot_2d) :: plt, fvplt

    ! Set up the model
    call mdl%initialize(4)
    mdl%coulomb_coefficient = mu_c
    mdl%static_coefficient = mu_s
    mdl%stribeck_velocity = vs
    mdl%viscous_damping = bv
    mdl%transition_sharpness = 1.0d2
    mdl%sliding_margin = 0.97d0
    mdl%reversal_sharpness = 1.0d1
    do i = 1, mdl%get_element_count()
        check = mdl%set_element_stiffness(i, ki(i))
        check = mdl%set_element_damping(i, bi(i))
        check = mdl%set_element_scaling(i, vi(i))
    end do

    ! Set up the integrator & solve
    neqn = 2 + mdl%get_state_variable_count()
    allocate(ic(neqn), source = 0.0d0)
    sys%fcn => ode_eqns
    call integrator%solve(sys, [0.0d0, 1.0d0], ic)
    sol = integrator%get_solution()

    ! Evaluate the model at each solution point
    n = size(sol, 1)
    allocate(F(n))
    do i = 1, n
        F(i) = mdl%evaluate(sol(i,1), sol(i,2), sol(i,3), nrm, sol(i,4:))
    end do

    ! Plot the solution
    call plt%initialize()
    call plt%push(sol(:,1), sol(:,2) * 1.0d3)
    call plt%set_x_axis_title("t [s]")
    call plt%set_y_axis_title("x(t) [mm]")
    call plt%draw()

    ! Plot the force-velocity solution
    call plt%clear_all()
    call plt%push(sol(:,3), F)
    call plt%set_x_axis_title("v [m/s]")
    call plt%set_y_axis_title("F [N]")
    call plt%draw()
contains
    subroutine ode_eqns(t, x, dxdt, args)
        !! The ODE's to solve.  The system is a base-excited spring-mass system
        !! with frictional damping.
        real(real64), intent(in) :: t
            !! The current simulation time.
        real(real64), intent(in), dimension(:) :: x
            !! The current state vector.
        real(real64), intent(out), dimension(:) :: dxdt
            !! The derivatives.
        class(*), intent(inout), optional :: args
            !! User arguments

        ! Local Variables
        integer(int32) :: k, j, nstates
        real(real64) :: y, dydt, Ff
        real(real64), allocatable, dimension(:) :: s, dsdt

        ! Define the state vector
        nstates = mdl%get_state_variable_count()
        allocate(s(nstates), dsdt(nstates))
        j = 2
        do k = 1, nstates
            j = j + 1
            s(k) = x(j)
        end do

        ! Define the excitation
        y = amp * sin(omega * t)
        dydt = amp * omega * cos(omega * t)

        ! Evaluate the friction model
        call mdl%state(t, x(1), x(2), nrm, s, dsdt)
        Ff = mdl%evaluate(t, x(1), x(2), nrm, s)

        ! Define the ODE's
        dxdt(1) = x(2)
        dxdt(2) = (stiffness * y - (Ff + stiffness * x(1))) / mass
        j = 2
        do k = 1, nstates
            j = j + 1
            dxdt(j) = dsdt(k)
        end do
    end subroutine
end program