submodule (friction) friction_model_routines
contains
! ------------------------------------------------------------------------------
pure module function fm_get_state_var_count(this) result(rst)
    class(friction_model), intent(in) :: this
    integer(int32) :: rst
    rst = 0
end function

! ------------------------------------------------------------------------------
end submodule