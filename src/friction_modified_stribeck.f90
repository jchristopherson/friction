submodule (friction) friction_modified_stribeck
contains
! ------------------------------------------------------------------------------
module function msf_eval(this, t, x, dxdt, nrm, svars) result(rst)
    class(modified_stribeck_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

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
module subroutine msf_to_array(this, x, err)
    class(modified_stribeck_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    x(1) = this%static_friction_coefficient
    x(2) = this%coulomb_friction_coefficient
    x(3) = this%stribeck_velocity
    x(4) = this%viscous_damping
    x(5) = this%stiffness
end subroutine

! ------------------------------------------------------------------------------
module subroutine msf_from_array(this, x, err)
    class(modified_stribeck_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    this%static_friction_coefficient = x(1)
    this%coulomb_friction_coefficient = x(2)
    this%stribeck_velocity = x(3)
    this%viscous_damping = x(4)
    this%stiffness = x(5)
end subroutine

! ------------------------------------------------------------------------------
pure module function msf_parameter_count(this) result(rst)
    class(modified_stribeck_model), intent(in) :: this
    integer(int32) :: rst
    rst = 5
end function

! ------------------------------------------------------------------------------
end submodule