submodule (friction) friction_coulomb
contains
! ------------------------------------------------------------------------------
module function cf_eval(this, t, x, dxdt, nrm, svars) result(rst)
    class(coulomb_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

    if (dxdt /= 0.0d0) then
        rst = this%friction_coefficient * nrm * sign(1.0d0, dxdt)
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure module function cf_has_state_vars(this) result(rst)
    class(coulomb_model), intent(in) :: this
    logical :: rst
    rst = .false.
end function

! ------------------------------------------------------------------------------
module subroutine cf_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    class(coulomb_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), dimension(:) :: svars
    real(real64), intent(out), dimension(:) :: dsdt
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
module subroutine cf_to_array(this, x, err)
    class(coulomb_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    x(1) = this%friction_coefficient
end subroutine

! ------------------------------------------------------------------------------
module subroutine cf_from_array(this, x, err)
    class(coulomb_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    this%friction_coefficient = x(1)
end subroutine

! ------------------------------------------------------------------------------
pure module function cf_parameter_count(this) result(rst)
    class(coulomb_model), intent(in) :: this
    integer(int32) :: rst
    rst = 1
end function

! ------------------------------------------------------------------------------
pure module function cf_get_state_var_count(this) result(rst)
    class(coulomb_model), intent(in) :: this
    integer(int32) :: rst
    rst = 0
end function

! ------------------------------------------------------------------------------
end submodule