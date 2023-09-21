program example
    use iso_fortran_env
    use fplot_core
    use csv_module
    use friction
    use fstats
    implicit none

    ! Local Variables
    type(csv_file) :: file
    logical :: ok
    integer(int32) :: npts, nfit, nstart, nend
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: t, x, v, nrm, frc, fmod
    type(generalized_maxwell_slip_model) :: mdl
    type(regression_statistics), allocatable, dimension(:) :: stats

    ! Plot Variables
    type(plot_2d) :: plt
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd
    type(plot_data_2d) :: pd1

    ! Read the data file
    call file%read("data/friction_data_1.csv", header_row = 1, status_ok = ok)
    if (.not.ok) then
        print *, "Could not open file."
        stop -1
    end if
    call file%get(1, t, ok)
    call file%get(2, x, ok)
    call file%get(3, v, ok)
    call file%get(4, nrm, ok)
    call file%get(5, frc, ok)

    ! Establish initial guesses at each parameter
    call mdl%initialize(2)  ! assume 2 elements - arbitrary at this point
    mdl%static_coefficient = 0.5d0
    mdl%coulomb_coefficient = 0.2d0
    mdl%stribeck_velocity = 1.0d-4
    mdl%attraction_coefficient = 1.0d0
    mdl%stiffness = 1.0d0
    mdl%viscous_damping = 0.0d0
    mdl%shape_parameter = 2.0d0
    
    ! Establish the element parameters - initial guesses
    ok = mdl%set_element_stiffness(1, 1.0d0)
    ok = mdl%set_element_damping(1, 1.0d0)
    ok = mdl%set_element_scaling(1, 0.5d0)

    ok = mdl%set_element_stiffness(2, 1.0d0)
    ok = mdl%set_element_damping(2, 1.0d0)
    ok = mdl%set_element_scaling(2, 0.5d0)

    ! The data set is very large.  For practical and efficiency reasons, reduce
    ! the amount of the original data set used for fitting.
    npts = size(t)
    dt = t(2) - t(1)
    nstart = floor(1.0d2 / dt)  ! start well into the data set to avoid start-up effects
    nend = min(nstart + floor(2.0d0 / dt), npts)
    nfit = nend - nstart + 1

    ! Attempt to fit the data
    allocate(fmod(nfit), stats(mdl%parameter_count()))
    call mdl%fit( &
        t(nstart:nend) - t(nstart), & ! time values at which data was sampled
        x(nstart:nend), &       ! relative displacement data
        v(nstart:nend), &       ! relative velocity data
        frc(nstart:nend), &     ! friction force data
        nrm(nstart:nend), &     ! normal force data
        stats = stats, &        ! retrieve fitting statistics of each parameter
        fmod = fmod &           ! retrieve fitted model results
    )
    
    ! Display the results
    call print_stats("Static Coefficient", mdl%static_coefficient, stats(1))
    call print_stats("Coulomb Coefficient", mdl%coulomb_coefficient, stats(2))
    call print_stats("Attraction Coefficient", mdl%attraction_coefficient, stats(3))
    call print_stats("Stiffness", mdl%stiffness, stats(4))
    call print_stats("Viscous Damping", mdl%viscous_damping, stats(5))
    call print_stats("Stribeck Velocity", mdl%stribeck_velocity, stats(6))
    call print_stats("Shape Parameter", mdl%shape_parameter, stats(7))
    call print_stats("Element 1 Stiffness", mdl%get_element_stiffness(1), stats(8))
    call print_stats("Element 1 Damping", mdl%get_element_damping(1), stats(9))
    call print_stats("Element 1 Scaling", mdl%get_element_scaling(1), stats(10))
    call print_stats("Element 2 Stiffness", mdl%get_element_stiffness(2), stats(11))
    call print_stats("Element 2 Damping", mdl%get_element_damping(2), stats(12))
    call print_stats("Element 2 Scaling", mdl%get_element_scaling(2), stats(13))

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t [s]")
    call yAxis%set_title("F_{f} [N]")
    
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t(nstart:nend) - t(nstart), frc(nstart:nend))
    call pd1%set_name("Raw Data")
    call plt%push(pd1)

    call pd1%define_data(t(nstart:nend) - t(nstart), fmod)
    call pd1%set_name("Model")
    call plt%push(pd1)

    call plt%draw()

contains
subroutine print_stats(paramname, paramval, stat)
    ! Arguments
    character(len = *), intent(in) :: paramname
    real(real64), intent(in) :: paramval
    type(regression_statistics), intent(in) :: stat

    ! Parameters
    character(len = *), parameter :: tab = achar(9)

    ! Print Output
    print "(A)", paramname
    print 100, tab // "Value: ", paramval
    print 100, tab // "Confidence Interval (95%): +/- ", stat%confidence_interval
    print 100, tab // "P-Value: ", stat%probability
    print 100, tab // "Standard Error: ", stat%standard_error
    print 100, tab // "T-Statistic: ", stat%t_statistic

    ! Formatting
100 format(A, G9.3)
end subroutine
end program