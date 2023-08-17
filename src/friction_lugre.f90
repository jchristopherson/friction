submodule (friction) friction_lugre
    use fstats
    use diffeq
contains
! ------------------------------------------------------------------------------
module function lg_eval(this, t, x, dxdt, nrm, svars) result(rst)
    ! Arguments
    class(lugre_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

    ! Local Variables
    real(real64) :: s1, dsdt(1)

    ! Process
    call this%state(t, x, dxdt, nrm, svars, dsdt)
    s1 = this%damping * exp(-(dxdt / this%stribeck_velocity)**2)
    rst = this%stiffness * svars(1) + s1 * dsdt(1) + &
        this%viscous_damping * dxdt
end function

! ------------------------------------------------------------------------------
pure module function lg_has_state_vars(this) result(rst)
    class(lugre_model), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
module subroutine lg_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    ! Arguments
    class(lugre_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), dimension(:) :: svars
    real(real64), intent(out), dimension(:) :: dsdt

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
module subroutine lg_to_array(this, x, err)
    class(lugre_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    x(1) = this%static_coefficient
    x(2) = this%coulomb_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%stiffness
    x(5) = this%damping
    x(6) = this%viscous_damping
    x(7) = this%shape_parameter
end subroutine

! ------------------------------------------------------------------------------
module subroutine lg_from_array(this, x, err)
    class(lugre_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    this%static_coefficient = x(1)
    this%coulomb_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%stiffness = x(4)
    this%damping = x(5)
    this%viscous_damping = x(6)
    this%shape_parameter = x(7)
end subroutine

! ------------------------------------------------------------------------------
pure module function lg_parameter_count(this) result(rst)
    class(lugre_model), intent(in) :: this
    integer(int32) :: rst
    rst = 7
end function

! ------------------------------------------------------------------------------
pure module function lg_get_state_var_count(this) result(rst)
    class(lugre_model), intent(in) :: this
    integer(int32) :: rst
    rst = 1
end function

! ------------------------------------------------------------------------------
end submodule
