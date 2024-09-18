module model
    use iso_fortran_env
    use friction
    implicit none

    ! Lu-Gre Friction Model
    type(lugre_model) :: mdl

    ! Dynamic System Parameters
    real(real64), parameter :: m = 5.0d0
    real(real64), parameter :: k = 5.0d3
    real(real64), parameter :: g = 9.81d0
    real(real64), parameter :: v = 0.1d0

contains
    subroutine equations(t, x, dxdt)
        ! Arguments
        real(real64), intent(in) :: t, x(:)
        real(real64), intent(out) :: dxdt(:)

        ! Local Variables
        real(real64) :: vr, F, N, s(1), dsdt(1)

        ! Compute the friction model
        N = m * g       ! normal force
        vr = v - x(2)   ! relative velocity
        s(1) = x(3)
        call mdl%state(t, 0.0d0, vr, N, s, dsdt)
        F = mdl%evaluate(t, 0.0d0, vr, N, s)

        ! Compute the model
        dxdt(1) = x(2)
        dxdt(2) = (F - k * x(1)) / m
        dxdt(3) = dsdt(1)
    end subroutine
end module

program example
    use iso_fortran_env
    use diffeq
    use fplot_core
    use model
    implicit none
    
    ! Local Variables
    type(runge_kutta_23) :: integrator
    type(ode_container) :: sys
    real(real64), allocatable :: sol(:,:)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd

    ! Set up the friction model
    mdl%static_coefficient = 0.598d0
    mdl%coulomb_coefficient = 0.40d0
    mdl%stribeck_velocity = 3.5d-3
    mdl%stiffness = 45.089d3
    mdl%damping = 2.6234d3
    mdl%viscous_damping = 0.0d0
    mdl%shape_parameter = 2.0d0

    ! Set up the integrator and solve the differential equations
    sys%fcn => equations
    call integrator%solve(sys, [0.0d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0])
    sol = integrator%get_solution()

    ! Plot the solution
    call plt%initialize()
    call plt%set_use_y2_axis(.true.)

    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call y2Axis%set_title("v_{drive} - v(t)")
    call lgnd%set_is_visible(.true.)

    call pd%define_data(sol(:,1), sol(:,2))
    call pd%set_name("Position")
    call pd%set_line_width(2.0)
    call plt%push(pd)

    call pd%define_data(sol(:,1), v - sol(:,3))
    call pd%set_draw_against_y2(.true.)
    call pd%set_name("Velocity")
    call plt%push(pd)
    call plt%draw()
end program