submodule (friction) friction_maxwell
contains
! ------------------------------------------------------------------------------
module function mx_eval(this, t, x, dxdt, nrm, svars) result(rst)
    class(maxwell_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), optional, dimension(:) :: svars
    real(real64) :: rst

    real(real64) :: s, d, delta

    delta = this%friction_coefficient / this%stiffness
    s = x - this%x_prev + this%d_prev
    d = sign(1.0d0, s) * min(abs(s), delta)
    rst = nrm * this%stiffness * d
    
    this%d_prev = d
    this%x_prev = x
end function

! ------------------------------------------------------------------------------
pure module function mx_has_state_vars(this) result(rst)
    class(maxwell_model), intent(in) :: this
    logical :: rst
    rst = .false.
end function

! ------------------------------------------------------------------------------
module subroutine mx_state_model(this, t, x, dxdt, nrm, svars, dsdt)
    class(maxwell_model), intent(inout) :: this
    real(real64), intent(in) :: t, x, dxdt, nrm
    real(real64), intent(in), dimension(:) :: svars
    real(real64), intent(out), dimension(:) :: dsdt
    dsdt = 0.0d0
end subroutine

! ------------------------------------------------------------------------------
module subroutine mx_to_array(this, x, err)
    class(maxwell_model), intent(in) :: this
    real(real64), intent(out), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    x(1) = this%stiffness
    x(2) = this%friction_coefficient
end subroutine

! ------------------------------------------------------------------------------
module subroutine mx_from_array(this, x, err)
    class(maxwell_model), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: x
    class(errors), intent(inout), optional, target :: err
    this%stiffness = x(1)
    this%friction_coefficient = x(2)
end subroutine

! ------------------------------------------------------------------------------
pure module function mx_parameter_count(this) result(rst)
    class(maxwell_model), intent(in) :: this
    integer(int32) :: rst
    rst = 2
end function

! ------------------------------------------------------------------------------
pure module function mx_get_state_var_count(this) result(rst)
    class(maxwell_model), intent(in) :: this
    integer(int32) :: rst
    rst = 0
end function

! ------------------------------------------------------------------------------
end submodule