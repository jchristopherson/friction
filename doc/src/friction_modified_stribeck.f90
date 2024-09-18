module friction_modified_stribeck
    use iso_fortran_env
    use friction_core
    use friction_stribeck
    use friction_errors
    use ferror
    implicit none
    private
    public :: modified_stribeck_model

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
        procedure, public :: reset => msf_reset
    end type

contains
! ------------------------------------------------------------------------------
function msf_eval(this, t, x, dxdt, nrm, svars) result(rst)
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

    real(real64) :: gv, a1, a2, zeta, delta

    a1 = this%coulomb_friction_coefficient * nrm / this%stiffness
    a2 = nrm * (this%static_friction_coefficient - &
        this%coulomb_friction_coefficient) / this%stiffness
    gv = a1 + a2 * exp(-abs(dxdt / this%stribeck_velocity)**2)
    zeta = x - this%x_prev + this%d_prev
    delta = sign(1.0d0, zeta) * min(abs(zeta), gv)
    rst = this%stiffness * delta + this%viscous_damping * dxdt
    
    this%d_prev = delta
    this%x_prev = x
end function

! ------------------------------------------------------------------------------
subroutine msf_to_array(this, x, err)
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
    x(1) = this%static_friction_coefficient
    x(2) = this%coulomb_friction_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%viscous_damping
    x(5) = this%stiffness
end subroutine

! ------------------------------------------------------------------------------
subroutine msf_from_array(this, x, err)
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
    this%static_friction_coefficient = x(1)
    this%coulomb_friction_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%viscous_damping = x(4)
    this%stiffness = x(5)
end subroutine

! ------------------------------------------------------------------------------
pure function msf_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(modified_stribeck_model), intent(in) :: this
        !! The modified_stribeck_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = 5
end function

! ------------------------------------------------------------------------------
subroutine msf_reset(this)
    !! Resets the friction model to it's original state.
    class(modified_stribeck_model), intent(inout) :: this
        !! The modified_stribeck_model object.
    this%x_prev = 0.0d0
    this%d_prev = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
end module