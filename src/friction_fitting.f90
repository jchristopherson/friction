submodule (friction) friction_fitting
    use fstats
    use fitpack
    use diffeq
    implicit none

! ------------------------------------------------------------------------------
    ! Variables specific to the fitting process
    real(real64), pointer, dimension(:) :: t_
    real(real64), pointer, dimension(:) :: x_
    real(real64), pointer, dimension(:) :: v_
    real(real64), pointer, dimension(:) :: f_
    real(real64), pointer, dimension(:) :: n_
    real(real64), pointer, dimension(:) :: initstate_
    type(fitpack_curve), pointer :: xinterp_
    type(fitpack_curve), pointer :: vinterp_
    type(fitpack_curve), pointer :: ninterp_
    type(ode_container), pointer :: mdl_
    class(friction_model), pointer :: fmdl_
    class(ode_integrator), pointer :: integrate_
contains
! ------------------------------------------------------------------------------
! Routine for fitting the friction model - uses module-level variables
subroutine fit_fcn(x, p, f, stop_)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, p
    real(real64), intent(out), dimension(:) :: f
    logical, intent(out) :: stop_

    ! Local Variables
    integer(int32) :: i, n, npts

    ! Initialization
    n = size(x)
    npts = n - fmdl_%get_constraint_equation_count()

    ! Assign the model parameters
    call fmdl_%from_array(p)

    ! Evaluate the friction model and compare the results
    call fmdl_%reset()
    do i = 1, npts
        f(i) = fmdl_%evaluate(t_(i), x_(i), v_(i), n_(i)) - f_(i)
    end do

    ! Evaluate constraints
    if (fmdl_%get_constraint_equation_count() > 0) then
        call fmdl_%constraint_equations(t_(:npts), x_(:npts), v_(:npts), &
            n_(:npts), f_(:npts), f(npts+1:))
    end if

    ! No need to stop
    stop_ = .false.
end subroutine

! Routine for fitting if internal variables are used by the model
subroutine internal_var_fit_fcn(x, p, f, stop_)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, p
    real(real64), intent(out), dimension(:) :: f
    logical, intent(out) :: stop_

    ! Local Variables
    integer(int32) :: i, n, npts
    real(real64), allocatable, dimension(:,:) :: dzdt

    ! Initialization
    n = size(x)
    npts = n - fmdl_%get_constraint_equation_count()

    ! Assign the model parameters
    call fmdl_%from_array(p)

    ! Integrate to determine the state variables
    call integrate_%solve(mdl_, t_, initstate_)
    dzdt = integrate_%get_solution()

    ! Evaluate the friction model and compare the results
    call fmdl_%reset()
    do i = 1, npts
        f(i) = fmdl_%evaluate(t_(i), x_(i), v_(i), n_(i), dzdt(i,2:)) - f_(i)
    end do

    ! Evaluate constraints
    if (fmdl_%get_constraint_equation_count() > 0) then
        call fmdl_%constraint_equations(t_(:npts), x_(:npts), v_(:npts), &
            n_(:npts), f_(:npts), f(npts+1:))
    end if

    ! No need to stop
    stop_ = .false.
end subroutine

! ODE Routine
subroutine internal_state_odes(t, z, dzdt)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: z
    real(real64), intent(out), dimension(:) :: dzdt

    ! Local Variables
    real(real64) :: x, v, n

    ! Interpolate to obtain the position, velocity, and normal force values
    ! corresponding to time t
    x = xinterp_%eval(t)
    v = vinterp_%eval(t)
    n = ninterp_%eval(t)

    ! Evaluate the friction model state equation
    call fmdl_%state(t, x, v, n, z, dzdt)
end subroutine

! ------------------------------------------------------------------------------
module subroutine fmdl_fit(this, t, x, v, f, n, weights, maxp, minp, &
    alpha, integrator, controls, settings, info, stats, fmod, resid, err)
    ! Arguments
    class(friction_model), intent(inout), target :: this
    real(real64), intent(in), target, dimension(:) :: t, x, v, f, n
    real(real64), intent(in), optional, dimension(:) :: weights, maxp, minp
    real(real64), intent(in), optional :: alpha
    class(ode_integrator), intent(inout), target, optional :: integrator
    type(iteration_controls), intent(in), optional :: controls
    type(lm_solver_options), intent(in), optional :: settings
    type(convergence_info), intent(out), optional :: info
    type(regression_statistics), intent(out), optional, dimension(:) :: stats
    real(real64), intent(out), optional, target, dimension(:) :: fmod, resid
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, npts, nparams, flag, np
    real(real64), allocatable, target, dimension(:) :: params, initstate, &
        tc, fc, fmc, rc
    real(real64), allocatable, dimension(:,:) :: dzdt
    real(real64), pointer, dimension(:) :: fmodptr, residptr, tptr, fptr
    real(real64), allocatable, target, dimension(:) :: fmoddef, residdef
    procedure(regression_function), pointer :: fcn
    type(fitpack_curve), target :: xinterp, vinterp, ninterp
    type(rosenbrock), target :: def_integrator
    type(ode_container), target :: mdl
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(t)
    nparams = this%parameter_count()
    np = npts + this%get_constraint_equation_count()
    if (present(integrator)) then
        integrate_ => integrator
    else
        integrate_ => def_integrator
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
        allocate(fmoddef(np), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 30
        fmodptr(1:np) => fmoddef(1:np)
    end if

    if (present(resid)) then
        if (size(resid) /= npts) go to 15
        residptr(1:npts) => resid(1:npts)
    else
        allocate(residdef(np), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 30
        residptr(1:np) => residdef(1:np)
    end if

    ! Are we using any additional constraints?
    if (this%get_constraint_equation_count() > 0) then
        allocate(tc(np), fc(np), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 30
        tptr(1:np) => tc(1:np)
        fptr(1:np) => fc(1:np)
        do i = 1, npts
            tptr(i) = t(i)
            fptr(i) = f(i)
        end do

        if (present(fmod)) then
            allocate(fmc(np), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 30
            fmodptr(1:np) => fmc(1:np)
        end if

        if (present(resid)) then
            allocate(rc(np), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 30
            residptr(1:np) => rc(1:np)
        end if
    else
        tptr(1:npts) => t
        fptr(1:npts) => f
    end if

    ! Assign pointers
    t_(1:npts) => t
    x_(1:npts) => x
    v_(1:npts) => v
    f_(1:npts) => f
    n_(1:npts) => n
    fmdl_ => this

    ! Compute the fit
    if (this%has_internal_state()) then
        fcn => internal_var_fit_fcn

        ! Define the interpolation objects & generate the fit
        flag = xinterp%new_fit(t, x)
        if (flag > 0) go to 40
        flag = vinterp%new_fit(t, v)
        if (flag > 0) go to 40
        flag = ninterp%new_fit(t, n)
        if (flag > 0) go to 40

        ! Set up the integrator
        mdl%fcn => internal_state_odes
        allocate(initstate(this%get_state_variable_count()), source = 0.0d0, &
            stat = flag)
        if (flag /= 0) go to 30

        ! Assign pointers
        mdl_ => mdl
        initstate_ => initstate
        xinterp_ => xinterp
        vinterp_ => vinterp
        ninterp_ => ninterp
    else
        fcn => fit_fcn
    end if

    call nonlinear_least_squares(fcn, tptr, fptr, params, fmodptr, residptr, &
        weights = weights, maxp = maxp, minp = minp, alpha = alpha, &
        controls = controls, settings = settings, info = info, stats = stats, &
        err = errmgr)
    if (errmgr%has_error_occurred()) return
    call this%from_array(params)

    ! Handle outputs, if constraints are employed
    if (this%get_constraint_equation_count() > 0) then
        if (present(fmod)) fmod = fmodptr(1:npts)
        if (present(resid)) resid = residptr(1:npts)
    end if

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

    ! Interpolation Error
40  continue
    call write_interpolation_error("fmdl_fit", flag, errmgr)
    return

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
    write(errmsg, 100) "Expected " // arrayname // " to be ", nexpect, &
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
    write(errmsg, 100) "Memory allocation error flag ", flag, " encountered."
    call err%report_error(fcn, trim(errmsg), FRICTION_MEMORY_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine write_interpolation_error(fcn, flag, err)
    ! Arguments
    character(len = *), intent(in) :: fcn
    integer(int32), intent(in) :: flag
    class(errors), intent(inout) :: err

    ! Local Variables
    character(len = 256) :: errmsg

    ! Process
    write(errmsg, 100) "Interpolation error flag ", flag, " encountered."
    call err%report_error(fcn, trim(errmsg), FRICTION_INVALID_OPERATION_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine fmdl_constraints(this, t, x, dxdt, nrm, f, rst)
    class(friction_model), intent(in) :: this
    real(real64), intent(in), dimension(:) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in), dimension(:) :: dxdt
    real(real64), intent(in), dimension(:) :: nrm
    real(real64), intent(in), dimension(:) :: f
    real(real64), intent(out), dimension(:) :: rst
    if (size(rst) > 0) rst = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
pure module function fmdl_get_constraint_count(this) result(rst)
    class(friction_model), intent(in) :: this
    integer(int32) :: rst
    rst = 0
end function

! ------------------------------------------------------------------------------
module subroutine fmdl_reset(this)
    class(friction_model), intent(inout) :: this
end subroutine

! ------------------------------------------------------------------------------
end submodule