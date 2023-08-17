program example
    use iso_fortran_env
    use friction
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: k = 1.0d3
    real(real64), parameter :: mu = 0.15d0
    real(real64), parameter :: Fnrm = 1.0d2
    real(real64), parameter :: amp = 1.0d-1
    real(real64), parameter :: freq = 2.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: dt = 1.0d-3

    ! Local Variables
    integer(int32) :: i
    real(real64) :: t(npts), x(npts), v(npts), F(npts)
    type(maxwell_model) :: mdl

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define the motion profiles
    t = (/ (dt * i, i = 0, npts - 1) /)
    x = amp * cos(2.0d0 * pi * freq * t)
    v = -2.0d0 * pi * freq * amp * sin(2.0d0 * pi * freq * t)

    ! Compute the friction force
    mdl%stiffness = k
    mdl%friction_coefficient = mu
    F = (/ (mdl%evaluate(t(i), x(i), v(i), Fnrm), i = 1, npts) /)

    ! Plot the resulting friction force - velocity curve
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("v(t)")
    call yAxis%set_title("F(t)")
    call yAxis%set_autoscale(.false.)
    call yAxis%set_limits(-1.5d0 * mu * Fnrm, 1.5d0 * mu * Fnrm)

    call pd%define_data(v, F)
    call pd%set_line_width(2.0)
    call plt%push(pd)
    call plt%draw()
    call plt%clear_all()

    ! Plot the friction force - time curve
    call xAxis%set_title("t")
    call pd%define_data(t, F)
    call plt%push(pd)
    call plt%draw()
end program