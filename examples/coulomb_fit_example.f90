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
    integer(int32) :: npts
    real(real64), allocatable, dimension(:) :: t, x, v, nrm, frc, fmod
    type(coulomb_model) :: mdl
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

    ! Attempt to fit the data
    npts = size(t)
    allocate(fmod(npts), stats(mdl%parameter_count()))
    mdl%friction_coefficient = 0.5d0    ! establish an initial guess
    call mdl%fit(t, x, v, frc, nrm, &
        minp = [0.0d0], &       ! the friction coefficient must be positive
        stats = stats, &        ! retrieve fitting statistics of each parameter
        fmod = fmod &           ! fitted model results
    )
    
    ! Display the results
    print 100, "Friction Coefficient: ", mdl%friction_coefficient
    print 101, "Confidence Interval (95%): +/- ", stats(1)%confidence_interval
    print 100, "P-Value: ", stats(1)%probability
    print 101, "Standard Error: ", stats(1)%standard_error
    print 101, "T-Statistic: ", stats(1)%t_statistic

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t [s]")
    call yAxis%set_title("F_{f} [N]")

    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(1.0d2, 1.1d2) ! limit range as the data set is large
    
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t, frc)
    call pd1%set_name("Raw Data")
    call plt%push(pd1)

    call pd1%define_data(t, fmod)
    call pd1%set_name("Model")
    call plt%push(pd1)

    call plt%draw()

    ! Formatting
100 format(A, F6.4)
101 format(A, G9.3)
end program