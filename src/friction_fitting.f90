submodule (friction) friction_fitting
    use fstats
contains
! ------------------------------------------------------------------------------
module subroutine fmdl_fit(this, t, x, v, f, n, usevel, weights, maxp, minp, &
    alpha, controls, settings, info, fmod, resid, err)
    ! Arguments
    class(friction_model), intent(inout) :: this
    real(real64), intent(in), target, dimension(:) :: t, x, v, f, n
    logical, intent(in), optional :: usevel
    real(real64), intent(in), optional, dimension(:) :: weights, maxp, minp
    real(real64), intent(in), optional :: alpha
    type(iteration_controls), intent(in), optional :: controls
    type(lm_solver_options), intent(in), optional :: settings
    type(convergence_info), intent(out), optional :: info
    real(real64), intent(out), optional, target, dimension(:) :: fmod, resid
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    logical :: uv
    integer(int32) :: npts, nparams, flag
    real(real64), allocatable, dimension(:) :: params
    real(real64), pointer, dimension(:) :: fmodptr, residptr, xptr
    real(real64), allocatable, target, dimension(:) :: fmoddef, residdef
    procedure(regression_function), pointer :: fcn
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(usevel)) then
        uv = usevel
    else
        uv = .true.
    end if
    npts = size(t)
    nparams = this%parameter_count()
    fcn => fit_fcn
    if (uv) then
        xptr(1:npts) => v(1:npts)
    else
        xptr(1:npts) => x(1:npts)
    end if

    ! Input Checking
    if (size(x) /= npts) go to 10
    if (size(v) /= npts) go to 11
    if (size(f) /= npts) go to 12
    if (size(n) /= npts) go to 13

    ! Memory Allocations
    allocate(params(nparams), stat = flag)
    if (flag /= 0) go to 30
    call this%to_array(params)

    if (present(fmod)) then
        if (size(fmod) /= npts) go to 14
        fmodptr(1:npts) => fmod(1:npts)
    else
        allocate(fmoddef(npts), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 30
        fmodptr(1:npts) => fmoddef(1:npts)
    end if

    if (present(resid)) then
        if (size(resid) /= npts) go to 15
        residptr(1:npts) => resid(1:npts)
    else
        allocate(residdef(npts), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 30
        residptr(1:npts) => residdef(1:npts)
    end if

    ! Compute the fit
    call nonlinear_least_squares(fcn, xptr, f, params, fmodptr, residptr, &
        weights = weights, maxp = maxp, minp = minp, alpha = alpha, &
        controls = controls, settings = settings, info = info, err = errmgr)
    if (errmgr%has_error_occurred()) return
    call this%from_array(params)

    ! End
    return

    ! X Array Size Error
10  continue
    call write_array_size_error("fmdl_fit", "x", npts, size(x), errmgr)
    return

    ! V Array Size Error
11  continue
    call write_array_size_error("fmdl_fit", "v", npts, size(x), errmgr)
    return

    ! F Array Size Error
12  continue
    call write_array_size_error("fmdl_fit", "f", npts, size(x), errmgr)
    return

    ! N Array Size Error
13  continue
    call write_array_size_error("fmdl_fit", "n", npts, size(x), errmgr)
    return

    ! FMod Array Size Error
14  continue
    call write_array_size_error("fmdl_fit", "fmod", npts, size(x), errmgr)
    return

    ! Resid Array Size Error
15  continue
    call write_array_size_error("fmdl_fit", "resid", npts, size(x), errmgr)
    return

    ! Memory Error
30  continue
    call write_memory_error("fmdl_fit", flag, errmgr)
    return

! ------------------------------------------------------------------------------
contains
    ! Define the function to fit
    subroutine fit_fcn(x_, p_, f_, stop_)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x_, p_
        real(real64), intent(out), dimension(:) :: f_
        logical, intent(out) :: stop_

        ! Local Variables
        integer(int32) :: i_

        ! Assign the model parameters
        call this%from_array(p_)

        ! Evaluate the friction model and compare the results
        do i_ = 1, size(x_)
            f_(i) = this%evaluate(t(i), x(i), v(i), n(i)) - f(i)
        end do

        ! No need to stop
        stop_ = .false.
    end subroutine
end subroutine

! ------------------------------------------------------------------------------
subroutine write_array_size_error(fcn, arrayname, nexpect, nactual, err)
    ! Arguments
    character(len = *), intent(in) :: fcn, arrayname
    integer(int32), intent(in) :: nexpect, nactual
    class(errors), intent(inout) :: err

    ! Local Variables
    character(len = 256) :: errmsg

    ! Process
    write(100, errmsg) "Expected " // arrayname // " to be ", nexpect, &
        " in size, but found it to be ", nactual, " in size."
    call err%report_error(fcn, trim(errmsg), FRICTION_ARRAY_SIZE_ERROR)

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine write_memory_error(fcn, flag, err)
    ! Arguments
    character(len = *), intent(in) :: fcn
    integer(int32), intent(in) :: flag
    class(errors), intent(inout) :: err

    ! Local Variables
    character(len = 256) :: errmsg

    ! Process
    write(100, errmsg) "Memory allocation error flag ", flag, " encountered."
    call err%report_error(fcn, trim(errmsg), FRICTION_MEMORY_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule