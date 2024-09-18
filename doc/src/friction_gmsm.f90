module friction_gmsm
    use iso_fortran_env
    use friction_core
    use ferror
    use friction_errors
    use :: ieee_arithmetic, only : ieee_value, IEEE_QUIET_NAN
    implicit none
    private
    public :: generalized_maxwell_slip_model

    ! The number of model parameters per element
    integer(int32), parameter :: PER_ELEMENT_COUNT = 3

    ! The number of common model parameters
    integer(int32), parameter :: COMMON_PARAMETER_COUNT = 7

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
contains
! ------------------------------------------------------------------------------
function gmsm_eval(this, t, x, dxdt, nrm, svars) result(rst)
    !! Evaluates the friction model given the defined parameter state.
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

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: dzdt, ki, bi

    ! Process
    n = this%get_element_count()
    rst = 0.0d0
    do i = 1, n
        ki = this%get_element_stiffness(i)
        bi = this%get_element_damping(i)
        dzdt = this%element_state(i, t, x, dxdt, nrm, svars(i))
        rst = rst + ki * svars(i) + bi * dzdt
    end do
    rst = rst + this%viscous_damping * dxdt
end function

! ------------------------------------------------------------------------------
pure function gmsm_has_state_vars(this) result(rst)
    !! Returns a value stating if the model relies upon internal
    !! state variables.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    logical :: rst
        !! Returns true if the model utilizes internal state variables;
        !! else, returns false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine gmsm_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    !! Evaluates the time derivatives of the internal friction state model.
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

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = this%get_element_count()
    do i = 1, n
        dsdt(i) = this%element_state(i, t, x, dxdt, nrm, svars(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine gmsm_to_array(this, x, err)
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

    ! Process
    if (size(x) /= this%parameter_count()) return
    x(1) = this%static_coefficient
    x(2) = this%coulomb_coefficient
    x(3) = this%attraction_coefficient
    x(4) = this%stiffness
    x(5) = this%viscous_damping
    x(6) = this%stribeck_velocity
    x(7) = this%shape_parameter
    x(8:) = this%m_params
end subroutine

! ------------------------------------------------------------------------------
subroutine gmsm_from_array(this, x, err)
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

    ! Process
    if (.not.allocated(this%m_params)) return
    if (size(x) /= this%parameter_count()) return
    this%static_coefficient = x(1)
    this%coulomb_coefficient = x(2)
    this%attraction_coefficient = x(3)
    this%stiffness = x(4)
    this%viscous_damping = x(5)
    this%stribeck_velocity = x(6)
    this%shape_parameter = x(7)
    this%m_params = x(8:)
end subroutine

! ------------------------------------------------------------------------------
pure function gmsm_parameter_count(this) result(rst)
    !! Gets the number of model parameters.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32) :: rst
        !! The number of model parameters.
    rst = this%get_element_count() * PER_ELEMENT_COUNT + COMMON_PARAMETER_COUNT
end function

! ------------------------------------------------------------------------------
pure function gmsm_get_state_var_count(this) result(rst)
    !! Gets the number of internal state variables used by the model.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32) :: rst
        !! The internal state variable count.
    rst = this%get_element_count()
end function

! ------------------------------------------------------------------------------
pure function gmsm_get_element_count(this) result(rst)
    !! Gets the number of friction elements in the model.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32) :: rst
        !! The number of friction elements in the model.
    if (.not.allocated(this%m_params)) then
        rst = 0
    else
        rst = this%m_nModels
    end if
end function

! ------------------------------------------------------------------------------
subroutine gmsm_initialize(this, n, err)
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

    ! Local Variables
    integer(int32) :: m, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    m = n * PER_ELEMENT_COUNT

    ! Input Checking
    if (n < 1) then
        ! TO DO: Handle error
    end if

    ! Process
    if (.not.allocated(this%m_params)) then
        allocate(this%m_params(m), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (size(this%m_params) /= m) then
        deallocate(this%m_params)
        allocate(this%m_params(m), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    this%m_nModels = n

    ! End
    return

    ! Memory Error Handling
10  continue
    return
end subroutine

! ------------------------------------------------------------------------------
pure function gmsm_get_element_stiffness(this, i) result(rst)
    !! Gets the stiffness term for the specified element.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32), intent(in) :: i
        !! The index of the element.
    real(real64) :: rst
        !! The requested value.
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i - 2)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
function gmsm_set_element_stiffness(this, i, x) result(rst)
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

    ! Initialization
    rst = .true.

    ! Input Checking
    if (.not.allocated(this%m_params)) then
        rst = .false.
        return
    end if
    if (i < 1 .or. i > this%get_element_count()) then
        rst = .false.
        return
    end if

    ! Process
    this%m_params(PER_ELEMENT_COUNT * i - 2) = x
end function

! ------------------------------------------------------------------------------
pure function gmsm_get_element_damping(this, i) result(rst)
    !! Gets the damping term for the specified element.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32), intent(in) :: i
        !! The index of the element.
    real(real64) :: rst
        !! The requested value.
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i - 1)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
function gmsm_set_element_damping(this, i, x) result(rst)
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

    ! Initialization
    rst = .true.

    ! Input Checking
    if (.not.allocated(this%m_params)) then
        rst = .false.
        return
    end if
    if (i < 1 .or. i > this%get_element_count()) then
        rst = .false.
        return
    end if

    ! Process
    this%m_params(PER_ELEMENT_COUNT * i - 1) = x
end function

! ------------------------------------------------------------------------------
pure function gmsm_get_element_scaling(this, i) result(rst)
    !! Gets the scaling factor for the specified element.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32), intent(in) :: i
        !! The index of the element.
    real(real64) :: rst
        !! The requested value.
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
function gmsm_set_element_scaling(this, i, x) result(rst)
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

    ! Initialization
    rst = .true.

    ! Input Checking
    if (.not.allocated(this%m_params)) then
        rst = .false.
        return
    end if
    if (i < 1 .or. i > this%get_element_count()) then
        rst = .false.
        return
    end if

    ! Process
    this%m_params(PER_ELEMENT_COUNT * i) = x
end function

! ------------------------------------------------------------------------------
pure function gmsm_stribeck_curve(this, dxdt, nrm) result(rst)
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

    ! Local Variables
    real(real64) :: a1, a2, s

    ! Process
    a1 = this%coulomb_coefficient * nrm / this%stiffness
    a2 = nrm * (this%static_coefficient - this%coulomb_coefficient) / &
        this%stiffness
    s = abs(dxdt) / this%stribeck_velocity
    rst = a1 + a2 / (1.0d0 + s**this%shape_parameter)
end function

! ------------------------------------------------------------------------------
pure function gmsm_element_state_model(this, i, t, x, dxdt, nrm, z) &
    result(rst)
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

    ! Local Variables
    real(real64) :: s, vi, C

    ! Compute the Stribeck function
    s = this%stribeck_function(dxdt, nrm)

    ! Process
    if (abs(z) <= s) then
        rst = dxdt
    else
        C = this%attraction_coefficient
        vi = this%get_element_scaling(i)
        rst = sign(1.0d0, dxdt) * vi * C * (1.0d0 - z / (vi * s))
    end if
end function

! ------------------------------------------------------------------------------
subroutine gmsm_constraints(this, t, x, dxdt, nrm, f, rst)
    !! Evaluates the constraint equation for the GMSM model.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
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

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: sm

    ! Process
    if (size(rst) /= 1) return
    sm = 0.0d0
    do i = 1, this%get_element_count()
        sm = sm + this%get_element_scaling(i)
    end do
    rst(1) = sm - 1.0d0 ! the sum of the scaling terms must be equal to 1
end subroutine

! ------------------------------------------------------------------------------
pure function gmsm_get_constraint_count(this) result(rst)
    !! Gets the number of constraint equations the model requires to
    !! be satisfied when fitting to data.
    class(generalized_maxwell_slip_model), intent(in) :: this
        !! The generalized_maxwell_slip_model object.
    integer(int32) :: rst
        !! The number of constraint equations.
    rst = 1
end function

! ------------------------------------------------------------------------------
end module