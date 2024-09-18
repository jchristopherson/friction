module friction_lugre
    use iso_fortran_env
    use friction_core
    use fstats
    use ferror
    use friction_errors
    implicit none
    private
    public :: lugre_model

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

contains
! ------------------------------------------------------------------------------
function lg_eval(this, t, x, dxdt, nrm, svars) result(rst)
    !! Evaluates the Lu-Gre friction model.
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

    ! Local Variables
    ! real(real64) :: s1, dsdt(1)
    real(real64) :: dsdt(1)

    ! Process
    call this%state(t, x, dxdt, nrm, svars, dsdt)
    ! s1 = this%damping * exp(-(dxdt / this%stribeck_velocity)**2)
    rst = this%stiffness * svars(1) + this%damping * dsdt(1) + &
        this%viscous_damping * dxdt
end function

! ------------------------------------------------------------------------------
pure function lg_has_state_vars(this) result(rst)
    !! Returns a value stating if the model relies upon internal
    !! state variables.
    class(lugre_model), intent(in) :: this
        !! The lugre_model object.
    logical :: rst
        !! Returns true if the model utilizes internal state variables;
        !! else, returns false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine lg_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    !! Evaluates the time derivatives of the internal friction state model.
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

    ! Local Variables
    real(real64) :: g, a1, a2, Fc, Fs

    ! Initialization
    Fc = nrm * this%coulomb_coefficient
    Fs = nrm * this%static_coefficient
    a1 = Fc / this%stiffness
    a2 = (Fs - Fc) / this%stiffness

    ! Compute the state variable derivative
    g = a1 + a2 / &
        (1.0d0 + (abs(dxdt) / this%stribeck_velocity)**this%shape_parameter)
    dsdt(1) = dxdt - abs(dxdt) * svars(1) / g
end subroutine

! ------------------------------------------------------------------------------
subroutine lg_to_array(this, x, err)
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
    x(1) = this%static_coefficient
    x(2) = this%coulomb_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%stiffness
    x(5) = this%damping
    x(6) = this%viscous_damping
    x(7) = this%shape_parameter
end subroutine

! ------------------------------------------------------------------------------
subroutine lg_from_array(this, x, err)
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
    this%static_coefficient = x(1)
    this%coulomb_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%stiffness = x(4)
    this%damping = x(5)
    this%viscous_damping = x(6)
    this%shape_parameter = x(7)
end subroutine

! ------------------------------------------------------------------------------
pure function lg_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(lugre_model), intent(in) :: this
        !! The lugre_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = 7
end function

! ------------------------------------------------------------------------------
pure function lg_get_state_var_count(this) result(rst)
    !! Gets the number of internal state variables used by the model.
    class(lugre_model), intent(in) :: this
        !! The lugre_model object.
    integer(int32) :: rst
        !! The internal state variable count.
    rst = 1
end function

! ------------------------------------------------------------------------------
end module
