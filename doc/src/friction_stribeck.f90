module friction_stribeck
    use iso_fortran_env
    use friction_core
    use ferror
    use friction_errors
    implicit none
    private
    public :: stribeck_model

    type, extends(friction_model) :: stribeck_model
        !! This type defines a basic Stribeck-based friction model.
        !!
        !! This model is defined as follows.
        !!
        !! $$ F = \text{sgn}\left( v \right) \left( \mu_c N + N (\mu_s - \mu_c) e^{-|v / v_s|^2} \right) + b_v v $$
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

contains
! ------------------------------------------------------------------------------
function sf_eval(this, t, x, dxdt, nrm, svars) result(rst)
    !! Evaluates the friction model given the defined parameter state.
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

    real(real64) :: fv, Fc, Fs

    Fc = this%coulomb_friction_coefficient * nrm
    Fs = this%static_friction_coefficient * nrm
    fv = Fc + (Fs - Fc) * exp(-abs(dxdt / this%stribeck_velocity)**2) + &
        this%viscous_damping * dxdt
    rst = sign(1.0d0, dxdt) * fv
end function

! ------------------------------------------------------------------------------
pure function sf_has_state_vars(this) result(rst)
    !! Returns a value stating if the model relies upon internal
    !! state variables.
    class(stribeck_model), intent(in) :: this
        !! The stribeck_model object.
    logical :: rst
        !! Returns true if the model utilizes internal state variables;
        !! else, returns false.
    rst = .false.
end function

! ------------------------------------------------------------------------------
subroutine sf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    !! Evaluates the time derivatives of the internal friction state model.
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
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine sf_to_array(this, x, err)
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
    x(1) = this%static_friction_coefficient
    x(2) = this%coulomb_friction_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%viscous_damping
end subroutine

! ------------------------------------------------------------------------------
subroutine sf_from_array(this, x, err)
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
    this%static_friction_coefficient = x(1)
    this%coulomb_friction_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%viscous_damping = x(4)
end subroutine

! ------------------------------------------------------------------------------
pure function sf_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(stribeck_model), intent(in) :: this
        !! The stribeck_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function sf_get_state_var_count(this) result(rst)
    !! Gets the number of internal state variables used by the model.
    class(stribeck_model), intent(in) :: this
        !! The stribeck_model object.
    integer(int32) :: rst
        !! The internal state variable count.
    rst = 0
end function

! ------------------------------------------------------------------------------
end module