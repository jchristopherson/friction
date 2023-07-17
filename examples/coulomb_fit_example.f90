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
    real(real64), allocatable, dimension(:) :: t, x, v, nrm, frc, fmod, resid
    type(coulomb_model) :: mdl
    type(regression_statistics), allocatable, dimension(:) :: stats

    ! Plot Variables
    type(plot_2d) :: plt
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd
    type(plot_data_2d) :: pd1, pd2

    ! Read the data file
    call file%read("friction_data_1.csv", header_row = 1, status_ok = ok)
    call file%get(1, t, ok)
    call file%get(2, x, ok)
    call file%get(3, v, ok)
    call file%get(4, nrm, ok)
    call file%get(5, frc, ok)

    ! Attempt to fit the data
    npts = size(t)
    allocate(fmod(npts), resid(npts), stats(mdl%parameter_count()))
    mdl%friction_coefficient = 0.5d0    ! establish an initial guess
    call mdl%fit(t, x, v, frc, nrm, &
        usevel = .false., &     ! fit force-displacement
        stats = stats, &        ! retrieve fitting statistics of each parameter
        fmod = fmod, &          ! fitted model results
        resid = resid&          ! residual errors (fit - actual)
    )

    ! Plot the results
    call plt%initialize()
    call plt%set_use_y2_axis(.true.)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t [s]")
    call yAxis%set_title("F_{f} [N]")
    call y2Axis%set_title("F_{model} - F_{actual} [N]")
    
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t, frc)
    call pd1%set_name("Raw Data")
    call plt%push(pd1)

    call pd1%define_data(t, fmod)
    call pd1%set_name("Model")
    call plt%push(pd1)

    call pd2%set_draw_against_y2(.true.)
    call pd2%define_data(t, resid)
    call pd2%set_name("Residual (Y2)")
    call plt%push(pd2)

    call plt%draw()
end program