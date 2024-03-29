submodule (friction) friction_gmsm
    use :: ieee_arithmetic, only : ieee_value, IEEE_QUIET_NAN

    ! The number of model parameters per element
    integer(int32), parameter :: PER_ELEMENT_COUNT = 3

    ! The number of common model parameters
    integer(int32), parameter :: COMMON_PARAMETER_COUNT = 7
contains
! ------------------------------------------------------------------------------
module function gmsm_eval(this, t, x, dxdt, nrm, svars) result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

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
pure module function gmsm_has_state_vars(this) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
module subroutine gmsm_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), dimension(:) :: svars
    real(real64), intent(out), dimension(:) :: dsdt

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = this%get_element_count()
    do i = 1, n
        dsdt(i) = this%element_state(i, t, x, dxdt, nrm, svars(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine gmsm_to_array(this, x, err)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err

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
module subroutine gmsm_from_array(this, x, err)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err

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
pure module function gmsm_parameter_count(this) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32) :: rst
    rst = this%get_element_count() * PER_ELEMENT_COUNT + COMMON_PARAMETER_COUNT
end function

! ------------------------------------------------------------------------------
pure module function gmsm_get_state_var_count(this) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32) :: rst
    rst = this%get_element_count()
end function

! ------------------------------------------------------------------------------
pure module function gmsm_get_element_count(this) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32) :: rst
    if (.not.allocated(this%m_params)) then
        rst = 0
    else
        rst = this%m_nModels
    end if
end function

! ------------------------------------------------------------------------------
module subroutine gmsm_initialize(this, n, err)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    integer(int32), intent(in) :: n
    class(errors), intent(inout), optional, target :: err

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
pure module function gmsm_get_element_stiffness(this, i) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i - 2)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
module function gmsm_set_element_stiffness(this, i, x) result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in) :: x
    logical :: rst

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
pure module function gmsm_get_element_damping(this, i) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i - 1)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
module function gmsm_set_element_damping(this, i, x) result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in) :: x
    logical :: rst

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
pure module function gmsm_get_element_scaling(this, i) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    if (this%get_element_count() >= i) then
        rst = this%m_params(PER_ELEMENT_COUNT * i)
    else
        rst = ieee_value(rst, IEEE_QUIET_NAN)
    end if
end function

! --------------------
module function gmsm_set_element_scaling(this, i, x) result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(inout) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in) :: x
    logical :: rst

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
pure module function gmsm_stribeck_curve(this, dxdt, nrm) result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(in) :: this
    real(real64), intent(in) :: dxdt, nrm
    real(real64) :: rst

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
pure module function gmsm_element_state_model(this, i, t, x, dxdt, nrm, z) &
    result(rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64), intent(in) :: t
    real(real64), intent(in) :: x
    real(real64), intent(in) :: dxdt
    real(real64), intent(in) :: nrm
    real(real64), intent(in) :: z
    real(real64) :: rst

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
module subroutine gmsm_constraints(this, t, x, dxdt, nrm, f, rst)
    ! Arguments
    class(generalized_maxwell_slip_model), intent(in) :: this
    real(real64), intent(in), dimension(:) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in), dimension(:) :: dxdt
    real(real64), intent(in), dimension(:) :: nrm
    real(real64), intent(in), dimension(:) :: f
    real(real64), intent(out), dimension(:) :: rst

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
pure module function gmsm_get_constraint_count(this) result(rst)
    class(generalized_maxwell_slip_model), intent(in) :: this
    integer(int32) :: rst
    rst = 1
end function

! ------------------------------------------------------------------------------
end submodule