submodule (friction) friction_lugre
contains
! ------------------------------------------------------------------------------
module function lg_eval(this, t, x, dxdt, nrm, svars) result(rst)
    ! Arguments
    class(lugre_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

    ! Local Variables
    real(real64) :: dsdt(1)

    ! Process
    call this%state(t, x, dxdt, nrm, svars, dsdt)
    rst = nrm * (this%stiffness * svars(1) + this%damping * dsdt(1) + &
        this%viscous_damping * dxdt)
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
    real(real64) :: g

    ! Compute the state variable derivative
    g = this%coulomb_coefficient + (this%static_coefficient - &
        this%coulomb_coefficient) * &
        exp(-abs(dxdt / this%stribeck_velocity)**this%shape_parameter)
    dsdt(1) = dxdt - this%stiffness * abs(dxdt) * svars(1) / g
end subroutine

! ------------------------------------------------------------------------------
end submodule