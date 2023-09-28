module friction
    !! Provides a collection of routines for modeling frictional behaviors
    !! of contacting bodies.
    use iso_fortran_env
    use fstats, only : convergence_info, regression_statistics, &
        iteration_controls, lm_solver_options
    use diffeq
    use ferror
    implicit none
    private
    public :: FRICTION_ARRAY_SIZE_ERROR
    public :: FRICTION_MEMORY_ERROR
    public :: FRICTION_INVALID_OPERATION_ERROR
    public :: friction_model
    public :: coulomb_model
    public :: lugre_model
    public :: maxwell_model
    public :: generalized_maxwell_slip_model
    public :: stribeck_model
    public :: modified_stribeck_model

    integer(int32), parameter :: FRICTION_ARRAY_SIZE_ERROR = 100000
        !! Defines an array size error.
    integer(int32), parameter :: FRICTION_MEMORY_ERROR = 100001
        !! Defines a memory allocation error.
    integer(int32), parameter :: FRICTION_INVALID_OPERATION_ERROR = 100002
        !! Defines an error within the opration of a routine.

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

    ! friction_fitting.f90
    interface
        module subroutine fmdl_fit(this, t, x, v, f, n, weights, maxp, &
            minp, alpha, integrator, controls, settings, info, stats, fmod, &
            resid, err)
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
                !! equations.  The defaults is a singly diagonally implicit
                !! Runge-Kutta method (4th order) that is suitable for 
                !! stiff ODE's.
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
        end subroutine

        module subroutine fmdl_constraints(this, t, x, dxdt, nrm, f, rst)
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
        end subroutine

        pure module function fmdl_get_constraint_count(this) result(rst)
            !! Gets the number of constraint equations the model requires to
            !! be satisfied when fitting to data.
            class(friction_model), intent(in) :: this
                !! The friction model object.
            integer(int32) :: rst
                !! The number of constraint equations.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(friction_model) :: coulomb_model
        !! Defines the basic Coulomb friction model.
        !!
        !! The Coulomb model is defined as follows.
        !!
        !! $$ F = sgn{ \left( v \right)} \mu_{c} N $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !! 
        !! \( v = \) Velocity
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        real(real64) :: friction_coefficient
            !! The Coulomb friction coefficient.
    contains
        procedure, public :: evaluate => cf_eval
        procedure, public :: has_internal_state => cf_has_state_vars
        procedure, public :: state => cf_state_model
        procedure, public :: to_array => cf_to_array
        procedure, public :: from_array => cf_from_array
        procedure, public :: parameter_count => cf_parameter_count
        procedure, public :: get_state_variable_count => cf_get_state_var_count
    end type

    ! friction_coulomb.f90
    interface
        module function cf_eval(this, t, x, dxdt, nrm, svars) result(rst)
            !! Evaluates the friction model given the defined parameter
            !! state.
            class(coulomb_model), intent(inout) :: this
                !! The coulomb_model object.
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

        pure module function cf_has_state_vars(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            class(coulomb_model), intent(in) :: this
                !! The coulomb_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        module subroutine cf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            !! Evaluates the time derivatives of the internal friction state
            !! model.
            class(coulomb_model), intent(inout) :: this
                !! The coulomb_model object.
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

        module subroutine cf_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(coulomb_model), intent(in) :: this
                !! The coulomb_model object.
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

        module subroutine cf_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(coulomb_model), intent(inout) :: this
                !! The coulomb_model object.
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

        pure module function cf_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(coulomb_model), intent(in) :: this
                !! The coulomb_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function

        pure module function cf_get_state_var_count(this) result(rst)
            !! Gets the number of internal state variables used by the model.
            class(coulomb_model), intent(in) :: this
                !! The coulomb_model object.
            integer(int32) :: rst
                !! The internal state variable count.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(friction_model) :: lugre_model
        !! Defines the Lu-Gre friction model.
        !!
        !! The Lu-Gre model is a bristle-type model that attempts to describe 
        !! friction using a bristle interpretation of the frictional surfaces.  
        !! The bristle-type models assume that the frictional behavior is 
        !! represented by an average deflection of elastic springs.  These 
        !! springs have their own stiffness and damping properties and act as a 
        !! typical spring-damper pair under small velocities; however, once 
        !! sufficient velocity occurs, the bristles slip resulting in 
        !! Coulomb-like sliding behavior.
        !!
        !! The Lu-Gre model is defined as follows.
        !!
        !! $$ F = \sigma_{0} z + \sigma_{1} \frac{dz}{dt} + \sigma_{2} v $$
        !! $$ \frac{dz}{dt} = v - \frac{\left| v \right| z}{g(v)} $$
        !! $$ g(v) = a_{1} + \frac{a_2}{1 + s^{\alpha}} $$
        !! $$ a_{1} = \frac{\mu_c N}{\sigma_{0}} $$
        !! $$ a_{2} = \frac{\mu_s N - \mu_c N}{\sigma_{0}} $$
        !! $$ s = \frac{\left| v \right|}{v_s} $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !! 
        !! \( v = \) Velocity
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        !!
        !! \( \mu_s = \) Static Friction Coefficient
        !!
        !! \( \sigma_{0} = \) Bristle Stiffness
        !!
        !! \( \sigma_{1} = \) Bristle Damping Coefficient
        !!
        !! \( \sigma_{2} = \) Viscous Damping Coefficient
        !! 
        !! \( \alpha = \) Stribeck Curve Shape Factor
        !!
        !! \( v_s = \) Stribeck Velocity Coefficient
        real(real64) :: static_coefficient
            !! The static friction coefficient.
        real(real64) :: coulomb_coefficient
            !! The Coulomb (dynamic) friction coefficient.
        real(real64) :: stribeck_velocity
            !! The Stribeck velocity term.
        real(real64) :: stiffness
            !! The bristle stiffness.
        real(real64) :: damping
            !! The bristle damping coefficient.
        real(real64) :: viscous_damping
            !! The viscous damping coefficient.
        real(real64) :: shape_parameter
            !! The Stribeck curve shape parameter.
    contains
        procedure, public :: evaluate => lg_eval
        procedure, public :: has_internal_state => lg_has_state_vars
        procedure, public :: state => lg_state_model
        procedure, public :: to_array => lg_to_array
        procedure, public :: from_array => lg_from_array
        procedure, public :: parameter_count => lg_parameter_count
        procedure, public :: get_state_variable_count => lg_get_state_var_count
    end type

    ! friction_lugre.f90
    interface
        module function lg_eval(this, t, x, dxdt, nrm, svars) result(rst)
            class(lugre_model), intent(inout) :: this
                !! The lugre_model object.
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

        pure module function lg_has_state_vars(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            class(lugre_model), intent(in) :: this
                !! The lugre_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        module subroutine lg_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            class(lugre_model), intent(inout) :: this
                !! The lugre_model object.
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

        module subroutine lg_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(lugre_model), intent(in) :: this
                !! The lugre_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See @ref
                !! parameter_count to determine the size of this array.  The 
                !! parameter order is as follows:
                !!
                !!  1. static_coefficient
                !!
                !!  2. coulomb_coefficient
                !!
                !!  3. stribeck_velocity
                !!
                !!  4. stiffness
                !!
                !!  5. damping
                !!
                !!  6. viscous_damping
                !!
                !!  7. shape_parameter
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        module subroutine lg_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(lugre_model), intent(inout) :: this
                !! The lugre_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count to 
                !! determine the size of this array.  The parameter order is as 
                !!  follows:
                !!
                !!  1. static_coefficient
                !!
                !!  2. coulomb_coefficient
                !!
                !!  3. stribeck_velocity
                !!
                !!  4. stiffness
                !!
                !!  5. damping
                !!
                !!  6. viscous_damping
                !!
                !!  7. shape_parameter
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function lg_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(lugre_model), intent(in) :: this
                !! The lugre_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function

        pure module function lg_get_state_var_count(this) result(rst)
            !! Gets the number of internal state variables used by the model.
            class(lugre_model), intent(in) :: this
                !! The lugre_model object.
            integer(int32) :: rst
                !! The internal state variable count.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(friction_model) :: maxwell_model
        !! Defines a single-element, Maxwell model.
        !!
        !! The signle-element, Maxwell model is defined as follows.
        !!
        !! $$ F = k \delta $$
        !! $$ \delta_{i+1} = sgn \left( x_{i+1} - x_{i} + \delta_{i} \right) \min \left( \left| x_{i+1} - x_{i} + \delta_{i} \right|, \Delta \right) $$
        !! $$ \Delta = \frac{N \mu_c}{k} $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !! 
        !! \( x = \) Position
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        !!
        !! \( k = \) Stiffness
        real(real64) :: stiffness
            !! The pre-sliding stiffness term.
        real(real64) :: friction_coefficient
            !! The Coulomb friction coefficient.

        ! Private, internal variables
        real(real64), private :: x_prev = 0.0d0
        real(real64), private :: d_prev = 0.0d0
    contains
        procedure, public :: evaluate => mx_eval
        procedure, public :: has_internal_state => mx_has_state_vars
        procedure, public :: state => mx_state_model
        procedure, public :: to_array => mx_to_array
        procedure, public :: from_array => mx_from_array
        procedure, public :: parameter_count => mx_parameter_count
        procedure, public :: get_state_variable_count => mx_get_state_var_count
    end type

    ! friction_maxwell.f90
    interface
        module function mx_eval(this, t, x, dxdt, nrm, svars) result(rst)
            !! Evaluates the friction model given the defined parameter
            !! state.
            class(maxwell_model), intent(inout) :: this
                !! The maxwell_model object.
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

        pure module function mx_has_state_vars(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            class(maxwell_model), intent(in) :: this
                !! The maxwell_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        module subroutine mx_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            !! Evaluates the time derivatives of the internal friction state
            !! model.
            class(maxwell_model), intent(inout) :: this
                !! The maxwell_model object.
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

        module subroutine mx_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(maxwell_model), intent(in) :: this
                !! The maxwell_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. stiffness
                !!
                !! 2. friction_coefficient
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        module subroutine mx_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(maxwell_model), intent(inout) :: this
                !! The maxwell_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. stiffness
                !!
                !! 2. friction_coefficient
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function mx_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(maxwell_model), intent(in) :: this
                !! The maxwell_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function

        pure module function mx_get_state_var_count(this) result(rst)
            !! Gets the number of internal state variables used by the model.
            class(maxwell_model), intent(in) :: this
                !! The maxwell_model object.
            integer(int32) :: rst
                !! The internal state variable count.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(friction_model) :: generalized_maxwell_slip_model
        !! A representation of the Generalized Maxwell Slip model.
        !!
        !! The Generalized Maxwell Slip model is defined as follows.
        !!
        !! $$ F = \sum\limits_{i=1}^{n} \left( k_i z_i + b_i \frac{dz_i}{dt} \right) + b_v v $$
        !! $$ \frac{dz_i}{dt} = \begin{cases} v & \text{if } |z_i| \le g(v) \\ sgn{ \left( v \right)} \nu_i C \left( 1 - \frac{z_i}{\nu_i g(v)} \right) & \text{otherwise} \end{cases} $$
        !! $$ g(v) = a_{1} + \frac{a_2}{1 + s^{\alpha}} $$
        !! $$ a_{1} = \frac{\mu_c N}{\sigma_{0}} $$
        !! $$ a_{2} = \frac{\mu_s N - \mu_c N}{\sigma_{0}} $$
        !! $$ s = \frac{\left| v \right|}{v_s} $$
        !! $$ \sum\limits_{i=1}^n {\nu_i} = 1 $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !!
        !! \( C = \) Attraction Coefficient
        !!
        !! \( x = \) Position
        !! 
        !! \( v = \) Velocity
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        !!
        !! \( \mu_s = \) Static Friction Coefficient
        !!
        !! \( k_i = \) i-th Element Stiffness
        !!
        !! \( b_i = \) i-th Element Damping Coefficient
        !!
        !! \( b_v = \) Viscous Damping Coefficient
        !!
        !! \( \sigma_0 = \) Frictional Stiffness
        !! 
        !! \( \alpha = \) Stribeck Curve Shape Factor
        !!
        !! \( v_s = \) Stribeck Velocity Coefficient
        !!
        !! \( \nu_i = \) i-th Element Scaling Factor
        integer(int32), private :: m_nModels = 0
            !! The number of elements in the model
        real(real64), private, allocatable, dimension(:) :: m_params
            !! An array containing the model parameters.
        real(real64) :: static_coefficient
            !! The static friction coefficient.
        real(real64) :: coulomb_coefficient
            !! The Coulomb (dynamic) friction coefficient.
        real(real64) :: stribeck_velocity
            !! The Stribeck velocity parameter.
        real(real64) :: shape_parameter
            !! The Stribeck curve shape parameter.
        real(real64) :: attraction_coefficient
            !! The attraction coefficient.
        real(real64) :: viscous_damping
            !! The viscous damping coefficient.
        real(real64) :: stiffness
            !! The frictional stiffness.
    contains
        procedure, public :: evaluate => gmsm_eval
        procedure, public :: has_internal_state => gmsm_has_state_vars
        procedure, public :: state => gmsm_state_model
        procedure, public :: to_array => gmsm_to_array
        procedure, public :: from_array => gmsm_from_array
        procedure, public :: parameter_count => gmsm_parameter_count
        procedure, public :: get_state_variable_count => gmsm_get_state_var_count
        procedure, public :: get_element_count => gmsm_get_element_count
        procedure, public :: initialize => gmsm_initialize
        procedure, public :: get_element_stiffness => gmsm_get_element_stiffness
        procedure, public :: set_element_stiffness => gmsm_set_element_stiffness
        procedure, public :: get_element_damping => gmsm_get_element_damping
        procedure, public :: set_element_damping => gmsm_set_element_damping
        procedure, public :: get_element_scaling => gmsm_get_element_scaling
        procedure, public :: set_element_scaling => gmsm_set_element_scaling
        procedure, public :: stribeck_function => gmsm_stribeck_curve
        procedure, public :: element_state => gmsm_element_state_model
        procedure, public :: get_constraint_equation_count => &
            gmsm_get_constraint_count
    end type

    ! friction_gmsm.f90
    interface
        module function gmsm_eval(this, t, x, dxdt, nrm, svars) result(rst)
            !! Evaluates the friction model given the defined parameter
            !! state.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
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

        pure module function gmsm_has_state_vars(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        module subroutine gmsm_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            !! Evaluates the time derivatives of the internal friction state
            !! model.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
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

        module subroutine gmsm_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is 
                !! as follows:
                !!
                !!  1. static_coefficient
                !!
                !!  2. coulomb_coefficient
                !!
                !!  3. attraction_coefficient
                !!
                !!  4. stiffness
                !!
                !!  5. viscous_damping
                !!
                !!  6. stribeck_velocity
                !!
                !!  7. shape_parameter
                !!
                !!  8. element stiffness
                !!  
                !!  9. element damping
                !!
                !!  10. element scaling ...
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        module subroutine gmsm_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count to 
                !! determine the size of this array. The parameter order is as 
                !! follows:
                !!
                !!  1. static_coefficient
                !!
                !!  2. coulomb_coefficient
                !!
                !!  3. attraction_coefficient
                !!
                !!  4. stiffness
                !!
                !!  5. viscous_damping
                !!
                !!  6. stribeck_velocity
                !!
                !!  7. shape_parameter
                !!
                !!  8. element stiffness
                !!  
                !!  9. element damping
                !!
                !!  10. element scaling ...
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function gmsm_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function

        pure module function gmsm_get_state_var_count(this) result(rst)
            !! Gets the number of internal state variables used by the model.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32) :: rst
                !! The internal state variable count.
        end function

        pure module function gmsm_get_element_count(this) result(rst)
            !! Gets the number of friction elements in the model.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32) :: rst
                !! The number of friction elements in the model.
        end function

        module subroutine gmsm_initialize(this, n, err)
            !! Initializes the model.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: n
                !! The number of friction elements.  This value must be at
                !! least 1.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function gmsm_get_element_stiffness(this, i) result(rst)
            !! Gets the stiffness term for the specified element.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64) :: rst
                !! The requested value.
        end function

        module function gmsm_set_element_stiffness(this, i, x) result(rst)
            !! Sets the stiffness term for the specified element.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64), intent(in) :: x
                !! The value.
            logical :: rst
                !! Returns true if the operation was successful; else, returns
                !! false if the object has not yet been initialized.
        end function

        pure module function gmsm_get_element_damping(this, i) result(rst)
            !! Gets the damping term for the specified element.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64) :: rst
                !! The requested value.
        end function

        module function gmsm_set_element_damping(this, i, x) result(rst)
            !! Sets the damping term for the specified element.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64), intent(in) :: x
                !! The value.
            logical :: rst
                !! Returns true if the operation was successful; else, returns
                !! false if the object has not yet been initialized.
        end function

        pure module function gmsm_get_element_scaling(this, i) result(rst)
            !! Gets the scaling factor for the specified element.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64) :: rst
                !! The requested value.
        end function

        module function gmsm_set_element_scaling(this, i, x) result(rst)
            !! Sets the scaling factor for the specified element.
            class(generalized_maxwell_slip_model), intent(inout) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
            real(real64), intent(in) :: x
                !! The value.
            logical :: rst
                !! Returns true if the operation was successful; else, returns
                !! false if the object has not yet been initialized.
        end function

        pure module function gmsm_stribeck_curve(this, dxdt, nrm) result(rst)
            !! Evaluates the Stribeck function for the model.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            real(real64), intent(in) :: dxdt
                !! The relative velocity between the contacting bodies.
            real(real64), intent(in) :: nrm
                !! The normal force between the contacting bodies.
            real(real64) :: rst
                !! The value of the Stribeck function.  The units are units of
                !! position.
        end function

        pure module function gmsm_element_state_model(this, i, t, x, dxdt, &
            nrm, z) result(rst)
            !! Computes the state equation for a single element.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32), intent(in) :: i
                !! The index of the element.
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
            real(real64), intent(in) :: z
                !! The current value of the state variable for the element.
            real(real64) :: rst
                !! The value of the state equation.
        end function

        module subroutine gmsm_constraints(this, t, x, dxdt, nrm, f, rst)
            !! Evaluates the constraint equation for the GMSM model.
            class(generalized_maxwell_slip_model), intent(in) :: this
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
        end subroutine

        pure module function gmsm_get_constraint_count(this) result(rst)
            !! Gets the number of constraint equations the model requires to
            !! be satisfied when fitting to data.
            class(generalized_maxwell_slip_model), intent(in) :: this
                !! The generalized_maxwell_slip_model object.
            integer(int32) :: rst
                !! The number of constraint equations.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(friction_model) :: stribeck_model
        !! This type defines a basic Stribeck-based friction model.
        !!
        !! This model is defined as follows.
        !!
        !! $$ F = \text{sgn}\left( v \right) f(v) $$
        !! $$ f(v) = \mu_c N + N (\mu_s - \mu_c) e^{-|v / v_s|^2} + b_v v $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !!
        !! \( x = \) Position
        !! 
        !! \( v = \) Velocity
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        !!
        !! \( \mu_s = \) Static Friction Coefficient
        !!
        !! \( b_v = \) Viscous Damping Coefficient
        !!
        !! \( v_s = \) Stribeck Velocity Coefficient
        real(real64) :: static_friction_coefficient
        real(real64) :: coulomb_friction_coefficient
        real(real64) :: stribeck_velocity
        real(real64) :: viscous_damping
    contains
        procedure, public :: evaluate => sf_eval
        procedure, public :: has_internal_state => sf_has_state_vars
        procedure, public :: state => sf_state_model
        procedure, public :: to_array => sf_to_array
        procedure, public :: from_array => sf_from_array
        procedure, public :: parameter_count => sf_parameter_count
        procedure, public :: get_state_variable_count => sf_get_state_var_count
    end type

    ! friction_stribeck.f90
    interface
        module function sf_eval(this, t, x, dxdt, nrm, svars) result(rst)
            !! Evaluates the friction model given the defined parameter
            !! state.
            class(stribeck_model), intent(inout) :: this
                !! The stribeck_model object.
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

        pure module function sf_has_state_vars(this) result(rst)
            !! Returns a value stating if the model relies upon internal
            !! state variables.
            class(stribeck_model), intent(in) :: this
                !! The stribeck_model object.
            logical :: rst
                !! Returns true if the model utilizes internal state variables;
                !! else, returns false.
        end function

        module subroutine sf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            !! Evaluates the time derivatives of the internal friction state
            !! model.
            class(stribeck_model), intent(inout) :: this
                !! The stribeck_model object.
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

        module subroutine sf_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(stribeck_model), intent(in) :: this
                !! The stribeck_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. static_friction_coefficient
                !!
                !! 2. coulomb_friction_coefficient
                !!
                !! 3. stribeck_velocity
                !!
                !! 4. viscous_damping
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        module subroutine sf_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(stribeck_model), intent(inout) :: this
                !! The stribeck_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. static_friction_coefficient
                !!
                !! 2. coulomb_friction_coefficient
                !!
                !! 3. stribeck_velocity
                !!
                !! 4. viscous_damping
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function sf_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(stribeck_model), intent(in) :: this
                !! The stribeck_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function

        pure module function sf_get_state_var_count(this) result(rst)
            !! Gets the number of internal state variables used by the model.
            class(stribeck_model), intent(in) :: this
                !! The stribeck_model object.
            integer(int32) :: rst
                !! The internal state variable count.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, extends(stribeck_model) :: modified_stribeck_model
        !! Defines a modification of the Stribeck model to account for 
        !! presliding displacement.  The presliding region utilizes a Maxwell
        !! type model then transitions to a traditional Stribeck model as
        !! slipping occurs.
        !!
        !! The model is defined as follows.
        !!
        !! $$ F = k \delta $$
        !! $$ delta_{i+1} = \text{sgn} \left( x_{i+1} - x_{i} + \delta_{i} \right) \min \left( x_{i+1} - x_{i} + \delta_{i}, g(v) \right) $$
        !! $$ g(v) = a_1 + a_2 e^{-|v / v_s|^2} $$
        !! $$ a_{1} = \frac{\mu_c N}{k} $$
        !! $$ a_{2} = \frac{\mu_s N - \mu_c N}{k} $$
        !!
        !! where:
        !!    
        !! \( F = \) Friction Force 
        !!    
        !! \( N = \) Normal Force
        !! 
        !! \( v = \) Velocity
        !!
        !! \( \mu_c = \) Coulomb Friction Coefficient
        !!
        !! \( \mu_s = \) Static Friction Coefficient
        !!
        !! \( k = \) Stiffness
        !!
        !! \( v_s = \) Stribeck Velocity Coefficient

        real(real64) :: stiffness
            !! The stiffness term.

        ! Private, internal variables
        real(real64), private :: x_prev = 0.0d0
        real(real64), private :: d_prev = 0.0d0
    contains
        procedure, public :: evaluate => msf_eval
        procedure, public :: to_array => msf_to_array
        procedure, public :: from_array => msf_from_array
        procedure, public :: parameter_count => msf_parameter_count
    end type

    ! friction_stribeck.f90
    interface
        module function msf_eval(this, t, x, dxdt, nrm, svars) result(rst)
            !! Evaluates the friction model given the defined parameter
            !! state.
            class(modified_stribeck_model), intent(inout) :: this
                !! The modified_stribeck_model object.
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

        module subroutine msf_to_array(this, x, err)
            !! Converts the parameters of the friction model into an array.
            class(modified_stribeck_model), intent(in) :: this
                !! The modified_stribeck_model object.
            real(real64), intent(out), dimension(:) :: x
                !! The array used to store the parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. static_friction_coefficient
                !!
                !! 2. coulomb_friction_coefficient
                !!
                !! 3. stribeck_velocity
                !!
                !! 4. viscous_damping
                !!
                !! 5. stiffness
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        module subroutine msf_from_array(this, x, err)
            !! Converts an array into the parameters for the friction model.
            class(modified_stribeck_model), intent(inout) :: this
                !! The modified_stribeck_model object.
            real(real64), intent(in), dimension(:) :: x
                !! The array of parameters.  See parameter_count 
                !! to determine the size of this array.  The parameter order is
                !! as follows.
                !!
                !! 1. static_friction_coefficient
                !!
                !! 2. coulomb_friction_coefficient
                !!
                !! 3. stribeck_velocity
                !!
                !! 4. viscous_damping
                !!
                !! 5. stiffness
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to
                !! provide error handling.
        end subroutine

        pure module function msf_parameter_count(this) result(rst)
            !! Gets the number of model parameters.
            class(modified_stribeck_model), intent(in) :: this
                !! The modified_stribeck_model object.
            integer(int32) :: rst
                !! The number of model parameters.
        end function
    end interface

! ------------------------------------------------------------------------------
end module