module friction_maxwell
    use iso_fortran_env
    use friction_core
    use ferror
    use friction_errors
    implicit none
    private
    public :: maxwell_model

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
        procedure, public :: reset => mx_reset
    end type

contains
! ------------------------------------------------------------------------------
function mx_eval(this, t, x, dxdt, nrm, svars) result(rst)
    !! Evaluates the friction model given the defined parameter state.
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

    real(real64) :: s, d, delta

    delta = nrm * this%friction_coefficient / this%stiffness
    s = x - this%x_prev + this%d_prev
    d = sign(1.0d0, s) * min(abs(s), delta)
    rst = this%stiffness * d
    
    this%d_prev = d
    this%x_prev = x
end function

! ------------------------------------------------------------------------------
pure function mx_has_state_vars(this) result(rst)
    !! Returns a value stating if the model relies upon internal
    !! state variables.
    class(maxwell_model), intent(in) :: this
        !! The maxwell_model object.
    logical :: rst
        !! Returns true if the model utilizes internal state variables;
        !! else, returns false.
    rst = .false.
end function

! ------------------------------------------------------------------------------
subroutine mx_state_model(this, t, x, dxdt, nrm, svars, dsdt)
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
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine mx_to_array(this, x, err)
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
    x(1) = this%stiffness
    x(2) = this%friction_coefficient
end subroutine

! ------------------------------------------------------------------------------
subroutine mx_from_array(this, x, err)
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
    this%stiffness = x(1)
    this%friction_coefficient = x(2)
end subroutine

! ------------------------------------------------------------------------------
pure function mx_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(maxwell_model), intent(in) :: this
        !! The maxwell_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = 2
end function

! ------------------------------------------------------------------------------
pure function mx_get_state_var_count(this) result(rst)
    !! Gets the number of internal state variables used by the model.
    class(maxwell_model), intent(in) :: this
        !! The maxwell_model object.
    integer(int32) :: rst
        !! The internal state variable count.
    rst = 0
end function

! ------------------------------------------------------------------------------
subroutine mx_reset(this)
    !! Resets the friction model to it's original state.
    class(maxwell_model), intent(inout) :: this
        !! The maxwell_model object.
    this%x_prev = 0.0d0
    this%d_prev = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
end module