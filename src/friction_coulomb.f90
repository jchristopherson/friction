module friction_coulomb
    use iso_fortran_env
    use friction_core
    use ferror
    implicit none
    private
    public :: coulomb_model

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

contains
! ------------------------------------------------------------------------------
function cf_eval(this, t, x, dxdt, nrm, svars) result(rst)
    !! Evaluates the friction model given the defined parameter state.
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

    if (dxdt /= 0.0d0) then
        rst = this%friction_coefficient * nrm * sign(1.0d0, dxdt)
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure function cf_has_state_vars(this) result(rst)
    !! Returns a value stating if the model relies upon internal
    !! state variables.
    class(coulomb_model), intent(in) :: this
        !! The coulomb_model object.
    logical :: rst
        !! Returns true if the model utilizes internal state variables;
        !! else, returns false.
    rst = .false.
end function

! ------------------------------------------------------------------------------
subroutine cf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    !! Evaluates the time derivatives of the internal friction state model.
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
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine cf_to_array(this, x, err)
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
    x(1) = this%friction_coefficient
end subroutine

! ------------------------------------------------------------------------------
subroutine cf_from_array(this, x, err)
    !! Converts the parameters of the friction model into an array.
    class(coulomb_model), intent(inout) :: this
        !! The coulomb_model object.
    real(real64), intent(in), dimension(:) :: x
        !! The array used to store the parameters.  See
        !! parameter_count to determine the size of this array.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to
        !! provide error handling.
    this%friction_coefficient = x(1)
end subroutine

! ------------------------------------------------------------------------------
pure function cf_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(coulomb_model), intent(in) :: this
        !! The coulomb_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = 1
end function

! ------------------------------------------------------------------------------
pure function cf_get_state_var_count(this) result(rst)
    !! Gets the number of internal state variables used by the model.
    class(coulomb_model), intent(in) :: this
        !! The coulomb_model object.
    integer(int32) :: rst
        !! The internal state variable count.
    rst = 0
end function

! ------------------------------------------------------------------------------
end module