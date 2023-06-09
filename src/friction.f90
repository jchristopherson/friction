!> @brief Provides a collection of routines for modeling frictional behaviors
!! of contacting bodies.
module friction
    use iso_fortran_env
    implicit none
    private
    public :: friction_model
    public :: coulomb_model
    public :: lugre_model

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
    end interface

! ------------------------------------------------------------------------------
end module