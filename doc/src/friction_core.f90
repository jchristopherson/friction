module friction_core
    use iso_fortran_env
    use ferror
    use fstats
    use fitpack
    use diffeq
    use friction_errors
    implicit none
    private
    public :: friction_model
    public :: friction_evaluation
    public :: friction_logical_query
    public :: friction_state_model
    public :: friction_model_to_array
    public :: friction_model_from_array
    public :: friction_integer_query
    public :: regression_statistics
    
    type, abstract :: friction_model
        !! Defines a generic friction model.
    contains
        procedure(friction_evaluation), deferred, public :: evaluate
        procedure(friction_logical_query), deferred, public :: &
            has_internal_state
        procedure(friction_state_model), deferred, public :: state
        procedure(friction_model_to_array), deferred, public :: to_array
        procedure(friction_model_from_array), deferred, public :: from_array
        procedure(friction_integer_query), deferred, public :: parameter_count
        procedure(friction_integer_query), deferred, public :: &
            get_state_variable_count
        procedure, public :: fit => fmdl_fit
        procedure, public :: constraint_equations => fmdl_constraints
        procedure, public :: get_constraint_equation_count => &
            fmdl_get_constraint_count
        procedure, public :: reset => fmdl_reset
    end type

    interface
        function friction_evaluation(this, t, x, dxdt, nrm, svars) result(rst)
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(inout) :: this
                !! The friction_model object.
            real(real64), intent(in) :: t
                !! The current simulation time value.
            real(real64), intent(in) :: x
                !! The current value of the relative position between
                !! the contacting bodies.
            real(real64), intent(in) :: dxdt
                !! The current value of the relative velocity between
                !! the contacting bodies.
            real(real64), intent(in) :: nrm
                !! The current normal force between the contacting 
                !! bodies.
            real(real64), intent(in), optional, dimension(:) :: svars
                !! An optional array containing any internal state
                !! variables the model may rely upon.
            real(real64) :: rst
                !! The friction force.
        end function

        pure function friction_logical_query(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            import friction_model
            class(friction_model), intent(in) :: this
                !! The friction_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        subroutine friction_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            !! Evaluates the time derivatives of the internal friction state
            !! model.
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(inout) :: this
                !! The friction_model object.
            real(real64), intent(in) :: t
                !! The current simulation time value.
            real(real64), intent(in) :: x
                !! The current value of the relative position between
                !! the contacting bodies.
            real(real64), intent(in) :: dxdt
                !! The current value of the relative velocity between
                !! the contacting bodies.
            real(real64), intent(in) :: nrm
                !! The current normal force between the contacting 
                !! bodies.
            real(real64), intent(in), dimension(:) :: svars
                !! An N-element array containing any internal state
                !! variables the model may rely upon.
            real(real64), intent(out), dimension(:) :: dsdt
                !! An N-element array where the state variable 
                !! derivatives are to be written.
        end subroutine

        subroutine friction_model_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            use iso_fortran_env, only : real64
            use ferror
            import friction_model
            class(friction_model), intent(in) :: this
                !! The friction_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See @ref
                !! parameter_count to determine the size of this array.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        subroutine friction_model_from_array(this, x, err)
            !!  Converts an array into the parameters for the friction model.
            use iso_fortran_env, only : real64
            use ferror
            import friction_model
            class(friction_model), intent(inout) :: this
                !! The friction_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count to 
                !! determine the size of this array.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure function friction_integer_query(this) result(rst)
            !! Gets an integer-valued parameter from the model
            use iso_fortran_env, only : int32
            import friction_model
            class(friction_model), intent(in) :: this
                !! The friction_model object.
            integer(int32) :: rst
                !! The model parameter.
        end function
    end interface

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
subroutine fmdl_fit(this, t, x, v, f, n, weights, maxp, minp, &
    alpha, integrator, controls, settings, info, stats, fmod, resid, err)
    !! Attempts to fit a friction model to the supplied data using a 
    !! Levenberg-Marquardt solver.
    class(friction_model), intent(inout), target :: this
        !! The friction model.  On output, the model is updated with the
        !! final, fitted parameters.
    real(real64), intent(in), target, dimension(:) :: t
        !! An N-element array containing the time points at which
        !! the friction data was sampled.  This array must contain 
        !! monotonically increasing data.
    real(real64), intent(in), target, dimension(:) :: x
        !! An N-element array containing the relative position
        !! data.
    real(real64), intent(in), target, dimension(:) :: v
        !! An N-element array containing the relative velocity
        !! data.
    real(real64), intent(in), target, dimension(:) :: f
        !! An N-element array containing the friction force data.
    real(real64), intent(in), target, dimension(:) :: n
        !! An N-element array containing the normal force data.
    real(real64), intent(in), optional, dimension(:) :: weights
        !! An optional N-element array that can be used to
        !!  weight specific data points.  The default is an array of 
        !! all ones such that all points are weighted equally.
    real(real64), intent(in), optional, dimension(:) :: maxp
        !! An M-element array (M = the number of model 
        !! parameters) containing a maximum limit for each model 
        !! parameter.
    real(real64), intent(in), optional, dimension(:) :: minp
        !! An M-element array containing the minimum limit for
        !! each model parameter.
    real(real64), intent(in), optional :: alpha
        !! An optional input that defines the significance 
        !! level at which to evaluate the confidence intervals. The 
        !! default value is 0.05 such that a 95% confidence interval 
        !! is calculated.
    class(ode_integrator), intent(inout), target, optional :: integrator
        !! An optional input, used in the event the model has internal 
        !! state variables, that provides integration of the state 
        !! equations.  The defaults is a 4th order Rosenbrock method.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for 
        !! the solver.
    type(convergence_info), intent(out), optional :: info
        !! An optional output that can be used to gain 
        !! information about the iterative solution and the nature of 
        !! the convergence.
    type(regression_statistics), intent(out), optional, dimension(:) :: stats
        !! An optional output array of M-elements that can be
        !! used to retrieve statistical information regarding the fit of
        !! each of the M model parameters.
    real(real64), intent(out), optional, target, dimension(:) :: fmod
        !! An optional N-element array used to provide the fitted model 
        !! results.
    real(real64), intent(out), optional, target, dimension(:) :: resid
        !! An optional N-element array containing the fitted residuals.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to
        !! provide error handling.

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
subroutine fmdl_constraints(this, t, x, dxdt, nrm, f, rst)
    !! Overload this routine to establish constraings for the model to
    !! be enforced as part of the fitting operation.
    class(friction_model), intent(in) :: this
        !! The friction_model object.
    real(real64), intent(in), dimension(:) :: t
        !! An N-element array containing the time points at which the
        !! data to be fit was sampled.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the relative motion data.
    real(real64), intent(in), dimension(:) :: dxdt
        !! An N-element array containing the relative velocity data.
    real(real64), intent(in), dimension(:) :: nrm
        !! An N-element array containing the normal force data.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the friction force data.
    real(real64), intent(out), dimension(:) :: rst
        !! An M-element array where the results of the constraint 
        !! equations will be written.  M must be equal to the 
        !! number of constraint equations for the model.
    if (size(rst) > 0) rst = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
pure function fmdl_get_constraint_count(this) result(rst)
    !! Gets the number of constraint equations the model requires to
    !! be satisfied when fitting to data.
    class(friction_model), intent(in) :: this
        !! The friction_model object.
    integer(int32) :: rst
        !! The number of constraint equations.
    rst = 0
end function

! ------------------------------------------------------------------------------
subroutine fmdl_reset(this)
    !! Resets the friction model to it's original state.
    class(friction_model), intent(inout) :: this
        !! The friction_model object.
end subroutine

! ------------------------------------------------------------------------------
end module