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