!> @brief Provides a collection of routines for modeling frictional behaviors
!! of contacting bodies.
module friction
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

    !> Defines an array size error.
    integer(int32), parameter :: FRICTION_ARRAY_SIZE_ERROR = 100000
    !> Defines a memory allocation error.
    integer(int32), parameter :: FRICTION_MEMORY_ERROR = 100001
    !> Defines an error within the opration of a routine.
    integer(int32), parameter :: FRICTION_INVALID_OPERATION_ERROR = 100002

    !> @brief Defines a generic friction model.
    type, abstract :: friction_model
    contains
        !> @brief Evaluates the friction model given the defined parameter
        !! state.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function evaluate( &
        !!  class(friction_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  optional real(real64) svars(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref friction_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An optional array containing any internal state
        !!  variables the model may rely upon.
        !!
        !! @return The friction force.
        procedure(friction_evaluation), deferred, public :: evaluate
        !> @brief Returns a value stating if the model relies upon internal
        !! state variables.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function has_internal_state( &
        !!  class(friction_model) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref friction_model object.
        !! @return Returns true if the model utilizes internal state variables;
        !!  else, returns false.
        procedure(friction_logical_query), deferred, public :: &
            has_internal_state
        !> @brief Evaluates the time derivatives of the internal friction state
        !! model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine state( &
        !!  class(friction_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  real(real64) svars(:), &
        !!  real(real64) dsdt(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref friction_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An N-element array containing any internal state
        !!  variables the model may rely upon.
        !! @param[out] dsdt An N-element array where the state variable 
        !!  derivatives are to be written.
        procedure(friction_state_model), deferred, public :: state
        !> @brief Converts the parameters of the friction model into an array.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine to_array( &
        !!  class(friction_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref friction_model object.
        !! @param[out] x The array used to store the parameters.  See @ref
        !!  parameter_count to determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure(friction_model_to_array), deferred, public :: to_array
        !> @brief Converts an array into the parameters for the friction model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine from_array( &
        !!  class(friction_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref friction_model object.
        !! @param[in] x The array of parameters.  See @ref parameter_count to 
        !!  determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure(friction_model_from_array), deferred, public :: from_array
        !> @brief Gets the number of model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function parameter_count(class(friction_model) this)
        !! @endcode
        !!
        !! @param[in] this The @ref friction_model object.
        !! @return The number of model parameters.
        procedure(friction_integer_query), deferred, public :: parameter_count
        !> @brief Gets the number of internal state variables used by the model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_state_variable_count( &
        !!  class(friction_model) this &
        !! )
        !! 
        !! @param[in] this The @ref friction_model object.
        !! @return The internal state variable count.
        procedure(friction_integer_query), deferred, public :: &
            get_state_variable_count
        !> @brief
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine fit( &
        !!  class(friction_model) this, &
        !!  real(real64) t(:), &
        !!  real(real64) x(:), &
        !!  real(real64) v(:), &
        !!  real(real64) f(:), &
        !!  real(real64) n(:), &
        !!  optional real(real64) weights(:), &
        !!  optional real(real64) maxp(:), &
        !!  optional real(real64) minp(:), &
        !!  optional real(real64) alpha, &
        !!  optional class(ode_integrator) integrator, &
        !!  optional type(iteration_controls) controls, &
        !!  optional type(lm_solver_options) settings, &
        !!  optional type(convergence_info) info, &
        !!  optional type(regression_statistics) stats(:), &
        !!  optional real(real64) fmod(:), &
        !!  optional real(real64) resid(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref friction_model.  On output, the model
        !!  is updated with the final, fitted parameters.
        !! @param[in] t An N-element array containing the time points at which
        !!  the friction data was sampled.  This array must contain 
        !!  monotonically increasing data.
        !! @param[in] x An N-element array containing the relative position
        !!  data.
        !! @param[in] v An N-element array containing the relative velocity
        !!  data.
        !! @param[in] f An N-element array containing the friction force data.
        !! @param[in] n An N-element array containing the normal force data.
        !! @param[in] weights An optional N-element array that can be used to
        !!  weight specific data points.  The default is an array of all ones
        !!  such that all points are weighted equally.
        !! @param[in] maxp An M-element array (M = the number of model 
        !!  parameters) containing a maximum limit for each model parameter.
        !! @param[in] minp An M-element array containing the minimum limit for
        !!  each model parameter.
        !! @param[in,out] integrator An optional input, used in the event the
        !!  model has internal state variables, that provides integration of the
        !!  state equations.  The defaults is a singly diagonally implicit
        !!  Runge-Kutta method (4th order) that is suitable for stiff ODE's.
        !! @param[in] alpha An optional input that defines the significance 
        !!  level at which to evaluate the confidence intervals. The default 
        !!  value is 0.05 such that a 95% confidence interval is calculated.
        !! @param[in] controls An optional input providing custom iteration 
        !!  controls.
        !! @param[in] settings An optional input providing custom settings for 
        !!  the solver.
        !! @param[out] info An optional output that can be used to gain 
        !!  information about the iterative solution and the nature of the 
        !!  convergence.
        !! @param[out] stats An optional output array of M-elements that can be
        !!  used to retrieve statistical information regarding the fit of each
        !!  of the M model parameters.
        !! @param[out] fmod An optional N-element array used to provide the
        !!  fitted model results.
        !! @param[out] resid An optional N-element array containing the fitted
        !!  residuals.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if any of the arrays are not
        !!      sized correctly.
        !!  - FRICTION_MEMORY_ERROR: Occurs if there is a memory allocation
        !!      error.
        !!  - FRICTION_INVALID_OPERATION_ERROR: Occurs if there is an error
        !!      interpolating within the provided data set.
        !!
        procedure, public :: fit => fmdl_fit
        !! @par Example
        !! The following example illustrates how to fit a Coulomb model to 
        !! experimental data.
        !!
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use fplot_core
        !!     use csv_module
        !!     use friction
        !!     use fstats
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(csv_file) :: file
        !!     logical :: ok
        !!     integer(int32) :: npts
        !!     real(real64), allocatable, dimension(:) :: t, x, v, nrm, frc, fmod
        !!     type(coulomb_model) :: mdl
        !!     type(regression_statistics), allocatable, dimension(:) :: stats
        !!
        !!     ! Plot Variables
        !!     type(plot_2d) :: plt
        !!     class(plot_axis), pointer :: xAxis, yAxis
        !!     class(legend), pointer :: lgnd
        !!     type(plot_data_2d) :: pd1
        !!
        !!     ! Read the data file
        !!     call file%read("data/friction_data_1.csv", header_row = 1, status_ok = ok)
        !!     if (.not.ok) then
        !!         print *, "Could not open file."
        !!         stop -1
        !!     end if
        !!     call file%get(1, t, ok)
        !!     call file%get(2, x, ok)
        !!     call file%get(3, v, ok)
        !!     call file%get(4, nrm, ok)
        !!     call file%get(5, frc, ok)
        !!
        !!     ! Attempt to fit the data
        !!     npts = size(t)
        !!     allocate(fmod(npts), stats(mdl%parameter_count()))
        !!     mdl%friction_coefficient = 0.5d0    ! establish an initial guess
        !!     call mdl%fit(t, x, v, frc, nrm, &
        !!         minp = [0.0d0], &       ! the friction coefficient must be positive
        !!         stats = stats, &        ! retrieve fitting statistics of each parameter
        !!         fmod = fmod &           ! fitted model results
        !!     )
        !!
        !!     ! Display the results
        !!     print 100, "Friction Coefficient: ", mdl%friction_coefficient
        !!     print 101, "Confidence Interval (95%): +/- ", stats(1)%confidence_interval
        !!     print 100, "P-Value: ", stats(1)%probability
        !!     print 101, "Standard Error: ", stats(1)%standard_error
        !!     print 101, "T-Statistic: ", stats(1)%t_statistic
        !!
        !!     ! Plot the results
        !!     call plt%initialize()
        !!     xAxis => plt%get_x_axis()
        !!     yAxis => plt%get_y_axis()
        !!     lgnd => plt%get_legend()
        !!
        !!     call xAxis%set_title("t [s]")
        !!     call yAxis%set_title("F_{f} [N]")
        !!
        !!     call xAxis%set_autoscale(.false.)
        !!     call xAxis%set_limits(1.0d2, 1.1d2) ! limit range as the data set is large
        !!
        !!     call lgnd%set_is_visible(.true.)
        !!
        !!     call pd1%define_data(t, frc)
        !!     call pd1%set_name("Raw Data")
        !!     call plt%push(pd1)
        !!
        !!     call pd1%define_data(t, fmod)
        !!     call pd1%set_name("Model")
        !!     call plt%push(pd1)
        !!
        !!     call plt%draw()
        !!
        !!     ! Formatting
        !! 100 format(A, F6.4)
        !! 101 format(A, G9.3)
        !! end program
        !! @endcode
        !! The program produces the following output.
        !! @code{.txt}
        !! Friction Coefficient: 0.3868
        !! Confidence Interval (95%): +/- 0.496E-03
        !! P-Value: 0.0000
        !! Standard Error: 0.253E-03
        !! T-Statistic: 0.153E+04
        !! @endcode
        !! @image html coulomb_fit_example.png
    end type

    interface
        function friction_evaluation(this, t, x, dxdt, nrm, svars) result(rst)
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), optional, dimension(:) :: svars
            real(real64) :: rst
        end function

        pure function friction_logical_query(this) result(rst)
            import friction_model
            class(friction_model), intent(in) :: this
            logical :: rst
        end function

        subroutine friction_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            use iso_fortran_env, only : real64
            import friction_model
            class(friction_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), dimension(:) :: svars
            real(real64), intent(out), dimension(:) :: dsdt
        end subroutine

        subroutine friction_model_to_array(this, x, err)
            use iso_fortran_env, only : real64
            use ferror
            import friction_model
            class(friction_model), intent(in) :: this
            real(real64), intent(out), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        subroutine friction_model_from_array(this, x, err)
            use iso_fortran_env, only : real64
            use ferror
            import friction_model
            class(friction_model), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure function friction_integer_query(this) result(rst)
            use iso_fortran_env, only : int32
            import friction_model
            class(friction_model), intent(in) :: this
            integer(int32) :: rst
        end function
    end interface

    ! friction_fitting.f90
    interface
        module subroutine fmdl_fit(this, t, x, v, f, n, weights, maxp, &
            minp, alpha, integrator, controls, settings, info, stats, fmod, &
            resid, err)
            class(friction_model), intent(inout), target :: this
            real(real64), intent(in), target, dimension(:) :: t, x, v, f, n
            real(real64), intent(in), optional, dimension(:) :: weights, maxp, &
                minp
            real(real64), intent(in), optional :: alpha
            class(ode_integrator), intent(inout), target, optional :: integrator
            type(iteration_controls), intent(in), optional :: controls
            type(lm_solver_options), intent(in), optional :: settings
            type(convergence_info), intent(out), optional :: info
            type(regression_statistics), intent(out), optional, dimension(:) :: stats
            real(real64), intent(out), optional, target, dimension(:) :: fmod, &
                resid
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the basic Coulomb friction model.
    !!
    !! @par Example
    !! The following example illustrates the evaluation of the Coulomb friction
    !! model for a system exposed to a sinusoidal velocity with a constant
    !! normal force.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use friction
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: npts = 1000
    !!     real(real64), parameter :: mu = 0.15d0
    !!     real(real64), parameter :: Fnrm = 1.0d2
    !!     real(real64), parameter :: amp = 1.0d-1
    !!     real(real64), parameter :: freq = 2.0d0
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!     real(real64), parameter :: dt = 1.0d-3
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i
    !!     real(real64) :: t(npts), x(npts), v(npts), F(npts)
    !!     type(coulomb_model) :: mdl
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!
    !!     ! Define the motion profiles
    !!     t = (/ (dt * i, i = 0, npts - 1) /)
    !!     x = amp * cos(2.0d0 * pi * freq * t)
    !!     v = -2.0d0 * pi * freq * amp * sin(2.0d0 * pi * freq * t)
    !!
    !!     ! Compute the friction force
    !!     mdl%friction_coefficient = mu
    !!     F = (/ (mdl%evaluate(t(i), x(i), v(i), Fnrm), i = 1, npts) /)
    !!
    !!     ! Plot the resulting friction force - velocity curve
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!
    !!     call xAxis%set_title("v(t)")
    !!     call yAxis%set_title("F(t)")
    !!     call yAxis%set_autoscale(.false.)
    !!     call yAxis%set_limits(-1.5d0 * mu * Fnrm, 1.5d0 * mu * Fnrm)
    !!
    !!     call pd%define_data(v, F)
    !!     call pd%set_line_width(2.0)
    !!     call plt%push(pd)
    !!     call plt%draw()
    !!     call plt%clear_all()
    !!
    !!     ! Plot the friction force - time curve
    !!     call xAxis%set_title("t")
    !!     call pd%define_data(t, F)
    !!     call plt%push(pd)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html coulomb_force_velocity.png
    !! @image html coulomb_force_time.png
    type, extends(friction_model) :: coulomb_model
        !> @brief The Coulomb friction coefficient.
        real(real64) :: friction_coefficient
    contains
        !> @brief Evaluates the friction model given the defined parameter
        !! state.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function evaluate( &
        !!  class(coulomb_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  optional real(real64) svars(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref coulomb_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An optional array containing any internal state
        !!  variables the model may rely upon.
        !!
        !! @return The friction force.
        procedure, public :: evaluate => cf_eval
        !> @brief Returns a value stating if the model relies upon internal
        !! state variables.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function has_internal_state( &
        !!  class(coulomb_model) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref coulomb_model object.
        !! @return Returns true if the model utilizes internal state variables;
        !!  else, returns false.
        procedure, public :: has_internal_state => cf_has_state_vars
        !> @brief Evaluates the time derivatives of the internal friction state
        !! model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine state( &
        !!  class(coulomb_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  real(real64) svars(:), &
        !!  real(real64) dsdt(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref coulomb_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An N-element array containing any internal state
        !!  variables the model may rely upon.
        !! @param[out] dsdt An N-element array where the state variable 
        !!  derivatives are to be written.
        procedure, public :: state => cf_state_model
        !> @brief Converts the parameters of the friction model into an array.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine to_array( &
        !!  class(coulomb_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref coulomb_model object.
        !! @param[out] x The array used to store the parameters.  See @ref
        !!  parameter_count to determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: to_array => cf_to_array
        !> @brief Converts an array into the parameters for the friction model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine from_array( &
        !!  class(coulomb_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref coulomb_model object.
        !! @param[in] x The array of parameters.  See @ref parameter_count to 
        !!  determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: from_array => cf_from_array
        !> @brief Gets the number of model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function parameter_count(class(coulomb_model) this)
        !! @endcode
        !!
        !! @param[in] this The @ref coulomb_model object.
        !! @return The number of model parameters.
        procedure, public :: parameter_count => cf_parameter_count
        !> @brief Gets the number of internal state variables used by the model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_state_variable_count( &
        !!  class(coulomb_model) this &
        !! )
        !! 
        !! @param[in] this The @ref coulomb_model object.
        !! @return The internal state variable count.
        procedure, public :: get_state_variable_count => cf_get_state_var_count
    end type

    ! friction_coulomb.f90
    interface
        module function cf_eval(this, t, x, dxdt, nrm, svars) result(rst)
            class(coulomb_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), optional, dimension(:) :: svars
            real(real64) :: rst
        end function

        pure module function cf_has_state_vars(this) result(rst)
            class(coulomb_model), intent(in) :: this
            logical :: rst
        end function

        module subroutine cf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            class(coulomb_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), dimension(:) :: svars
            real(real64), intent(out), dimension(:) :: dsdt
        end subroutine

        module subroutine cf_to_array(this, x, err)
            class(coulomb_model), intent(in) :: this
            real(real64), intent(out), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine cf_from_array(this, x, err)
            class(coulomb_model), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function cf_parameter_count(this) result(rst)
            class(coulomb_model), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function cf_get_state_var_count(this) result(rst)
            class(coulomb_model), intent(in) :: this
            integer(int32) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the Lu-Gre friction model.
    !!
    !! @par Remarks
    !! The Lu-Gre model is a heuristic-type model that attempts to describe 
    !! friction using a bristle-type interpretation of the frictional surfaces.
    !! The bristle-type models assume that the frictional behavior is 
    !! represented by an average deflection of elastic springs.  These springs
    !! have their own stiffness and damping properties and act as a typical
    !! spring-damper pair under small velocities; however, once sufficient 
    !! velocity occurs, the bristles slip resulting in Coulomb-like sliding 
    !! behavior.
    !!
    !! @par Example
    !! The following example illustrates the use of the Lu-Gre model to 
    !! define the behavior of a spring-restrained mass riding on a frictional
    !! belt moving at a constant velocity.  The normal force is the weight of
    !! the sprung mass.
    !!
    !! @par
    !! The system has the equation of motion 
    !! \f$ m \frac{d^2x}{dt^2} + k x = F_{f}(t) \f$, where \f$ F_{f}(t) \f$ is
    !! the friction force output by the Lu-Gre model.
    !!
    !! @image html schematic.png
    !!
    !! @par
    !! This example utilizes the 
    !! [DIFFEQ](https://github.com/jchristopherson/diffeq) library to solve the
    !! differential equations using a Dormand-Prince Runge-Kutta integrator.
    !! @code{.f90}
    !! module model
    !!     use iso_fortran_env
    !!     use friction
    !!     implicit none
    !!
    !!     ! Lu-Gre Friction Model
    !!     type(lugre_model) :: mdl
    !!
    !!     ! Dynamic System Parameters
    !!     real(real64), parameter :: m = 5.0d0
    !!     real(real64), parameter :: k = 5.0d3
    !!     real(real64), parameter :: g = 9.81d0
    !!     real(real64), parameter :: v = 0.1d0
    !!
    !! contains
    !!     subroutine equations(t, x, dxdt)
    !!         ! Arguments
    !!         real(real64), intent(in) :: t, x(:)
    !!         real(real64), intent(out) :: dxdt(:)
    !!
    !!         ! Local Variables
    !!         real(real64) :: vr, F, N, s(1), dsdt(1)
    !!
    !!         ! Compute the friction model
    !!         N = m * g       ! normal force
    !!         vr = v - x(2)   ! relative velocity
    !!         s(1) = x(3)
    !!         call mdl%state(t, 0.0d0, vr, N, s, dsdt)
    !!         F = mdl%evaluate(t, 0.0d0, vr, N, s)
    !!
    !!         ! Compute the model
    !!         dxdt(1) = x(2)
    !!         dxdt(2) = (F - k * x(1)) / m
    !!         dxdt(3) = dsdt(1)
    !!     end subroutine
    !! end module
    !!
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq
    !!     use fplot_core
    !!     use model
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(bsrk32_integrator) :: integrator
    !!     type(ode_container) :: sys
    !!     real(real64), allocatable :: sol(:,:)
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd
    !!     class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Set up the friction model
    !!     mdl%static_coefficient = 0.598d0
    !!     mdl%coulomb_coefficient = 0.40d0
    !!     mdl%stribeck_velocity = 3.5d-3
    !!     mdl%stiffness = 45.089d3
    !!     mdl%damping = 2.6234d3
    !!     mdl%viscous_damping = 0.0d0
    !!     mdl%shape_parameter = 2.0d0
    !!
    !!     ! Set up the integrator and solve the differential equations
    !!     sys%fcn => equations
    !!     sol = integrator%solve(sys, [0.0d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0])
    !!
    !!     ! Plot the solution
    !!     call plt%initialize()
    !!     call plt%set_use_y2_axis(.true.)
    !!
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     y2Axis => plt%get_y2_axis()
    !!     lgnd => plt%get_legend()
    !!
    !!     call xAxis%set_title("t")
    !!     call yAxis%set_title("x(t)")
    !!     call y2Axis%set_title("v_{drive} - v(t)")
    !!     call lgnd%set_is_visible(.true.)
    !!
    !!     call pd%define_data(sol(:,1), sol(:,2))
    !!     call pd%set_name("Position")
    !!     call pd%set_line_width(2.0)
    !!     call plt%push(pd)
    !!
    !!     call pd%define_data(sol(:,1), v - sol(:,3))
    !!     call pd%set_draw_against_y2(.true.)
    !!     call pd%set_name("Velocity")
    !!     call plt%push(pd)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html lugre_example.png
    type, extends(friction_model) :: lugre_model
        !> @brief The static friction coefficient.
        real(real64) :: static_coefficient
        !> @brief The Coulomb (dynamic) friction coefficient.
        real(real64) :: coulomb_coefficient
        !> @brief The Stribeck velocity.
        real(real64) :: stribeck_velocity
        !> @brief The frictional stiffness.
        real(real64) :: stiffness
        !> @brief The frictional damping coefficient.
        real(real64) :: damping
        !> @brief The viscous damping coefficient.
        real(real64) :: viscous_damping
        !> @brief The Stribeck curve shape parameter.
        real(real64) :: shape_parameter
    contains
        !> @brief Evaluates the friction model given the defined parameter
        !! state.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function evaluate( &
        !!  class(lugre_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  optional real(real64) svars(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref lugre_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An optional array containing any internal state
        !!  variables the model may rely upon.
        !!
        !! @return The friction force.
        procedure, public :: evaluate => lg_eval
        !> @brief Returns a value stating if the model relies upon internal
        !! state variables.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function has_internal_state( &
        !!  class(lugre_model) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref lugre_model object.
        !! @return Returns true if the model utilizes internal state variables;
        !!  else, returns false.
        procedure, public :: has_internal_state => lg_has_state_vars
        !> @brief Evaluates the time derivatives of the internal friction state
        !! model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine state( &
        !!  class(lugre_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  real(real64) svars(:), &
        !!  real(real64) dsdt(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref lugre_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An N-element array containing any internal state
        !!  variables the model may rely upon.
        !! @param[out] dsdt An N-element array where the state variable 
        !!  derivatives are to be written.
        procedure, public :: state => lg_state_model
        !> @brief Converts the parameters of the friction model into an array.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine to_array( &
        !!  class(lugre_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref lugre_model object.
        !! @param[out] x The array used to store the parameters.  See @ref
        !!  parameter_count to determine the size of this array.  The parameter 
        !!  order is as follows:
        !!  1. static_coefficient
        !!  2. coulomb_coefficient
        !!  3. stribeck_velocity
        !!  4. stiffness
        !!  5. damping
        !!  6. viscous_damping
        !!  7. shape_parameter
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: to_array => lg_to_array
        !> @brief Converts an array into the parameters for the friction model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine from_array( &
        !!  class(lugre_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref lugre_model object.
        !! @param[in] x The array of parameters.  See @ref parameter_count to 
        !!  determine the size of this array.  The parameter order is as 
        !!  follows:
        !!  1. static_coefficient
        !!  2. coulomb_coefficient
        !!  3. stribeck_velocity
        !!  4. stiffness
        !!  5. damping
        !!  6. viscous_damping
        !!  7. shape_parameter
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: from_array => lg_from_array
        !> @brief Gets the number of model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function parameter_count(class(lugre_model) this)
        !! @endcode
        !!
        !! @param[in] this The @ref lugre_model object.
        !! @return The number of model parameters.
        procedure, public :: parameter_count => lg_parameter_count
        !> @brief Gets the number of internal state variables used by the model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_state_variable_count( &
        !!  class(friction_model) this &
        !! )
        !! 
        !! @param[in] this The @ref friction_model object.
        !! @return The internal state variable count.
        procedure, public :: get_state_variable_count => lg_get_state_var_count
    end type

    ! friction_lugre.f90
    interface
        module function lg_eval(this, t, x, dxdt, nrm, svars) result(rst)
            class(lugre_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), optional, dimension(:) :: svars
            real(real64) :: rst
        end function

        pure module function lg_has_state_vars(this) result(rst)
            class(lugre_model), intent(in) :: this
            logical :: rst
        end function

        module subroutine lg_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            class(lugre_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), dimension(:) :: svars
            real(real64), intent(out), dimension(:) :: dsdt
        end subroutine

        module subroutine lg_to_array(this, x, err)
            class(lugre_model), intent(in) :: this
            real(real64), intent(out), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine lg_from_array(this, x, err)
            class(lugre_model), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function lg_parameter_count(this) result(rst)
            class(lugre_model), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function lg_get_state_var_count(this) result(rst)
            class(lugre_model), intent(in) :: this
            integer(int32) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a single-element, Maxwell Slip model.
    !!
    !! @par Example
    !! The following example illustrates the evaluation of the Maxwell friction
    !! model for a system exposed to a sinusoidal velocity with a constant
    !! normal force.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use friction
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: npts = 1000
    !!     real(real64), parameter :: k = 1.0d3
    !!     real(real64), parameter :: mu = 0.15d0
    !!     real(real64), parameter :: Fnrm = 1.0d2
    !!     real(real64), parameter :: amp = 1.0d-1
    !!     real(real64), parameter :: freq = 2.0d0
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!     real(real64), parameter :: dt = 1.0d-3
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i
    !!     real(real64) :: t(npts), x(npts), v(npts), F(npts)
    !!     type(maxwell_model) :: mdl
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!
    !!     ! Define the motion profiles
    !!     t = (/ (dt * i, i = 0, npts - 1) /)
    !!     x = amp * cos(2.0d0 * pi * freq * t)
    !!     v = -2.0d0 * pi * freq * amp * sin(2.0d0 * pi * freq * t)
    !!
    !!     ! Compute the friction force
    !!     mdl%stiffness = k
    !!     mdl%friction_coefficient = mu
    !!     F = (/ (mdl%evaluate(t(i), x(i), v(i), Fnrm), i = 1, npts) /)
    !!
    !!     ! Plot the resulting friction force - velocity curve
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!
    !!     call xAxis%set_title("v(t)")
    !!     call yAxis%set_title("F(t)")
    !!     call yAxis%set_autoscale(.false.)
    !!     call yAxis%set_limits(-1.5d0 * mu * Fnrm, 1.5d0 * mu * Fnrm)
    !!
    !!     call pd%define_data(v, F)
    !!     call pd%set_line_width(2.0)
    !!     call plt%push(pd)
    !!     call plt%draw()
    !!     call plt%clear_all()
    !!
    !!     ! Plot the friction force - time curve
    !!     call xAxis%set_title("t")
    !!     call pd%define_data(t, F)
    !!     call plt%push(pd)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html maxwell_force_velocity.png
    !! @image html maxwell_force_time.png
    type, extends(friction_model) :: maxwell_model
        !> @brief The presliding stiffness term.
        real(real64) :: stiffness
        !> @brief The Coulomb friction coefficient.
        real(real64) :: friction_coefficient
        ! Private, internal variables
        real(real64), private :: x_prev = 0.0d0
        real(real64), private :: d_prev = 0.0d0
    contains
        !> @brief Evaluates the friction model given the defined parameter
        !! state.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function evaluate( &
        !!  class(maxwell_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  optional real(real64) svars(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref maxwell_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An optional array containing any internal state
        !!  variables the model may rely upon.
        !!
        !! @return The friction force.
        procedure, public :: evaluate => mx_eval
        !> @brief Returns a value stating if the model relies upon internal
        !! state variables.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function has_internal_state( &
        !!  class(maxwell_model) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref maxwell_model object.
        !! @return Returns true if the model utilizes internal state variables;
        !!  else, returns false.
        procedure, public :: has_internal_state => mx_has_state_vars
        !> @brief Evaluates the time derivatives of the internal friction state
        !! model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine state( &
        !!  class(maxwell_model) this, &
        !!  real(real64) t, &
        !!  real(real64) x, &
        !!  real(real64) dxdt, &
        !!  real(real64) nrm, &
        !!  real(real64) svars(:), &
        !!  real(real64) dsdt(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref maxwell_model object.
        !! @param[in] t The current simulation time value.
        !! @param[in] x The current value of the relative position between
        !!  the contacting bodies.
        !! @param[in] dxdt The current value of the relative velocity between
        !!  the contacting bodies.
        !! @param[in] nrm The current normal force between the contacting 
        !!  bodies.
        !! @param[in] svars An N-element array containing any internal state
        !!  variables the model may rely upon.
        !! @param[out] dsdt An N-element array where the state variable 
        !!  derivatives are to be written.
        procedure, public :: state => mx_state_model
        !> @brief Converts the parameters of the friction model into an array.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine to_array( &
        !!  class(maxwell_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref maxwell_model object.
        !! @param[out] x The array used to store the parameters.  See @ref
        !!  parameter_count to determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: to_array => mx_to_array
        !> @brief Converts an array into the parameters for the friction model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine from_array( &
        !!  class(maxwell_model) this, &
        !!  real(real64) x(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref maxwell_model object.
        !! @param[in] x The array of parameters.  See @ref parameter_count to 
        !!  determine the size of this array.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - FRICTION_ARRAY_SIZE_ERROR: Occurs if @p x is not sized 
        !!      appropriately.
        procedure, public :: from_array => mx_from_array
        !> @brief Gets the number of model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function parameter_count(class(maxwell_model) this)
        !! @endcode
        !!
        !! @param[in] this The @ref maxwell_model object.
        !! @return The number of model parameters.
        procedure, public :: parameter_count => mx_parameter_count
        !> @brief Gets the number of internal state variables used by the model.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_state_variable_count( &
        !!  class(maxwell_model) this &
        !! )
        !! 
        !! @param[in] this The @ref maxwell_model object.
        !! @return The internal state variable count.
        procedure, public :: get_state_variable_count => mx_get_state_var_count
    end type

    ! friction_maxwell.f90
    interface
        module function mx_eval(this, t, x, dxdt, nrm, svars) result(rst)
            class(maxwell_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), optional, dimension(:) :: svars
            real(real64) :: rst
        end function

        pure module function mx_has_state_vars(this) result(rst)
            class(maxwell_model), intent(in) :: this
            logical :: rst
        end function

        module subroutine mx_state_model(this, t, x, dxdt, nrm, svars, dsdt)
            class(maxwell_model), intent(inout) :: this
            real(real64), intent(in) :: t, x, dxdt, nrm
            real(real64), intent(in), dimension(:) :: svars
            real(real64), intent(out), dimension(:) :: dsdt
        end subroutine

        module subroutine mx_to_array(this, x, err)
            class(maxwell_model), intent(in) :: this
            real(real64), intent(out), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine mx_from_array(this, x, err)
            class(maxwell_model), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function mx_parameter_count(this) result(rst)
            class(maxwell_model), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function mx_get_state_var_count(this) result(rst)
            class(maxwell_model), intent(in) :: this
            integer(int32) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
end module