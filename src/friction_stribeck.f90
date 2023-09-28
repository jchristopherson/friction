submodule (friction) friction_stribeck
contains
! ------------------------------------------------------------------------------
module function sf_eval(this, t, x, dxdt, nrm, svars) result(rst)
    class(stribeck_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

    real(real64) :: fv, Fc, Fs

    Fc = this%coulomb_friction_coefficient * nrm
    Fs = this%static_friction_coefficient * nrm
    fv = Fc + (Fs - Fc) * exp(-abs(dxdt / this%stribeck_velocity)**2) + &
        this%viscous_damping * dxdt
    rst = sign(1.0d0, dxdt) * fv
    
    this%d_prev = delta
    this%x_prev = x
end function

! ------------------------------------------------------------------------------
pure module function sf_has_state_vars(this) result(rst)
    class(stribeck_model), intent(in) :: this
    logical :: rst
    rst = .false.
end function

! ------------------------------------------------------------------------------
module subroutine sf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    class(stribeck_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), dimension(:) :: svars
    real(real64), intent(out), dimension(:) :: dsdt
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
module subroutine sf_to_array(this, x, err)
    class(stribeck_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    x(1) = this%static_friction_coefficient
    x(2) = this%coulomb_friction_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%viscous_damping
end subroutine

! ------------------------------------------------------------------------------
module subroutine sf_from_array(this, x, err)
    class(stribeck_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    this%static_friction_coefficient = x(1)
    this%coulomb_friction_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%viscous_damping = x(4)
end subroutine

! ------------------------------------------------------------------------------
pure module function sf_parameter_count(this) result(rst)
    class(stribeck_model), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure module function sf_get_state_var_count(this) result(rst)
    class(stribeck_model), intent(in) :: this
    integer(int32) :: rst
    rst = 0
end function

! ------------------------------------------------------------------------------
end submodule