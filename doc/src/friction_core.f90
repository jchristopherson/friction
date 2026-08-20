module friction_core
    use iso_fortran_env
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

        subroutine friction_model_to_array(this, x)
            !! Converts the parameters of the friction model into an array.
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(in) :: this
                !! The friction_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See @ref
                !! parameter_count to determine the size of this array.
        end subroutine

        subroutine friction_model_from_array(this, x)
            !!  Converts an array into the parameters for the friction model.
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(inout) :: this
                !! The friction_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count to 
                !! determine the size of this array.
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
    type fit_data
        real(real64), pointer, dimension(:) :: t
        real(real64), pointer, dimension(:) :: x
        real(real64), pointer, dimension(:) :: v
        real(real64), pointer, dimension(:) :: f
        real(real64), pointer, dimension(:) :: n
        real(real64), pointer, dimension(:) :: initstate
        type(fitpack_curve), pointer :: xinterp
        type(fitpack_curve), pointer :: vinterp
        type(fitpack_curve), pointer :: ninterp
        type(ode_container), pointer :: mdl
        class(friction_model), pointer :: fmdl
        class(ode_integrator), pointer :: integrate
    end type

contains
! ------------------------------------------------------------------------------
! Routine for fitting the friction model - uses module-level variables
subroutine fit_fcn(x, p, f, stop_, args)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, p
    real(real64), intent(out), dimension(:) :: f
    logical, intent(out) :: stop_
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: i, n, npts
    real(real64), pointer, dimension(:) :: t_, x_, v_, n_, f_
    class(friction_model), pointer :: fmdl_

    ! Initialization
    n = size(x)    
    if (.not.present(args)) then
        stop_ = .true.
        return
    end if
    select type (args)
    class is (fit_data)
        t_ => args%t
        x_ => args%x
        v_ => args%v
        n_ => args%n
        f_ => args%f
        fmdl_ => args%fmdl
    end select
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
subroutine internal_var_fit_fcn(x, p, f, stop_, args)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x, p
    real(real64), intent(out), dimension(:) :: f
    logical, intent(out) :: stop_
    class(*), intent(inout), optional :: args

    ! Local Variables
    integer(int32) :: i, n, npts
    real(real64), allocatable, dimension(:,:) :: dzdt
    real(real64), pointer, dimension(:) :: t_, x_, v_, n_, f_, initstate_
    class(friction_model), pointer :: fmdl_
    class(ode_integrator), pointer :: integrate_
    type(ode_container), pointer :: mdl_

    ! Initialization
    n = size(x)
    if (.not.present(args)) then
        stop_ = .true.
        return
    end if
    select type (args)
    class is (fit_data)
        t_ => args%t
        x_ => args%x
        v_ => args%v
        n_ => args%n
        f_ => args%f
        initstate_ => args%initstate
        fmdl_ => args%fmdl
        integrate_ => args%integrate
        mdl_ => args%mdl
    end select
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
subroutine internal_state_odes(t, z, dzdt, args)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: z
    real(real64), intent(out), dimension(:) :: dzdt
    class(*), intent(inout), optional :: args

    ! Local Variables
    real(real64) :: x, v, n
    type(fitpack_curve), pointer :: xinterp_, vinterp_, ninterp_
    class(friction_model), pointer :: fmdl_

    ! Initialization
    if (.not.present(args)) then
        return
    end if
    select type (args)
    class is (fit_data)
        xinterp_ => args%xinterp
        vinterp_ => args%vinterp
        ninterp_ => args%ninterp
        fmdl_ => args%fmdl
    end select

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
    alpha, integrator, controls, settings, info, stats, fmod, resid, &
    initial_state)
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
    real(real64), intent(in), optional, dimension(:) :: initial_state
        !! An optional array containing the initial conditions for the
        !! model's internal state variables.  Its size must match 
        !! @ref get_state_variable_count.  Only used if the model relies
        !! upon internal state variables.  The default is an array of
        !! all zeros.

    ! Local Variables
    integer(int32) :: i, npts, nparams, np, flag
    real(real64), allocatable, target, dimension(:) :: params, initstate, &
        tc, fc, fmc, rc, wc
    real(real64), allocatable, dimension(:,:) :: dzdt
    real(real64), pointer, dimension(:) :: fmodptr, residptr, tptr, fptr, wptr
    real(real64), allocatable, target, dimension(:) :: fmoddef, residdef
    procedure(regression_function), pointer :: fcn
    type(fitpack_curve), target :: xinterp, vinterp, ninterp
    type(rosenbrock), target :: def_integrator
    type(ode_container), target :: mdl
    type(fit_data) :: args
    
    ! Initialization
    npts = size(t)
    nparams = this%parameter_count()
    np = npts + this%get_constraint_equation_count()
    if (present(integrator)) then
        args%integrate => integrator
    else
        args%integrate => def_integrator
    end if

    ! Input Checking
    if (size(x) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
    if (size(v) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
    if (size(f) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
    if (size(n) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
    if (present(weights)) then
        if (size(weights) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
    end if
    if (present(initial_state)) then
        if (size(initial_state) /= this%get_state_variable_count()) &
            error stop FRICTION_ARRAY_SIZE_ERROR
    end if

    ! Memory Allocations
    allocate(params(nparams))
    call this%to_array(params)

    if (present(fmod)) then
        if (size(fmod) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
        fmodptr(1:npts) => fmod(1:npts)
    else
        allocate(fmoddef(np), source = 0.0d0)
        fmodptr(1:np) => fmoddef(1:np)
    end if

    if (present(resid)) then
        if (size(resid) /= npts) error stop FRICTION_ARRAY_SIZE_ERROR
        residptr(1:npts) => resid(1:npts)
    else
        allocate(residdef(np), source = 0.0d0)
        residptr(1:np) => residdef(1:np)
    end if

    ! Are we using any additional constraints?
    if (this%get_constraint_equation_count() > 0) then
        allocate(tc(np), fc(np), source = 0.0d0)
        tptr(1:np) => tc(1:np)
        fptr(1:np) => fc(1:np)
        do i = 1, npts
            tptr(i) = t(i)
            fptr(i) = f(i)
        end do

        if (present(fmod)) then
            allocate(fmc(np), source = 0.0d0)
            fmodptr(1:np) => fmc(1:np)
        end if

        if (present(resid)) then
            allocate(rc(np), source = 0.0d0)
            residptr(1:np) => rc(1:np)
        end if

        ! The constraint rows carry no data weighting of their own; pad
        ! the user-supplied (or default) weights out to np elements so
        ! the sizes agree with tptr/fptr as required by the solver.
        allocate(wc(np), source = 1.0d0)
        if (present(weights)) wc(1:npts) = weights
        wptr(1:np) => wc(1:np)
    else
        tptr(1:npts) => t
        fptr(1:npts) => f
        if (present(weights)) then
            allocate(wc(npts), source = weights)
            wptr(1:npts) => wc(1:npts)
        else
            wptr => null()
        end if
    end if

    ! Assign pointers
    args%t(1:npts) => t
    args%x(1:npts) => x
    args%v(1:npts) => v
    args%f(1:npts) => f
    args%n(1:npts) => n
    args%fmdl => this

    ! Compute the fit
    if (this%has_internal_state()) then
        fcn => internal_var_fit_fcn

        ! Define the interpolation objects & generate the fit
        flag = xinterp%new_fit(t, x)
        if (flag > 0) error stop FRICTION_INVALID_OPERATION_ERROR
        flag = vinterp%new_fit(t, v)
        if (flag > 0) error stop FRICTION_INVALID_OPERATION_ERROR
        flag = ninterp%new_fit(t, n)
        if (flag > 0) error stop FRICTION_INVALID_OPERATION_ERROR

        ! Set up the integrator
        mdl%fcn => internal_state_odes
        if (present(initial_state)) then
            initstate = initial_state
        else
            allocate(initstate(this%get_state_variable_count()), source = 0.0d0)
        end if

        ! Assign pointers
        args%mdl => mdl
        args%initstate => initstate
        args%xinterp => xinterp
        args%vinterp => vinterp
        args%ninterp => ninterp
    else
        fcn => fit_fcn
    end if

    if (associated(wptr)) then
        call nonlinear_least_squares(fcn, tptr, fptr, params, fmodptr, &
            residptr, weights = wptr, maxp = maxp, minp = minp, &
            alpha = alpha, controls = controls, settings = settings, &
            info = info, stats = stats, args = args)
    else
        call nonlinear_least_squares(fcn, tptr, fptr, params, fmodptr, &
            residptr, maxp = maxp, minp = minp, alpha = alpha, &
            controls = controls, settings = settings, info = info, &
            stats = stats, args = args)
    end if
    call this%from_array(params)

    ! Handle outputs, if constraints are employed
    if (this%get_constraint_equation_count() > 0) then
        if (present(fmod)) fmod = fmodptr(1:npts)
        if (present(resid)) resid = residptr(1:npts)
    end if
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