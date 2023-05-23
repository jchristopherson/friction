program test
    use friction_model_tests
    implicit none

    ! Local Variables
    logical :: check
    integer(int32) :: flag

    ! Initialization
    flag = 0

    ! Process
    check = test_coulomb()
    if (.not.check) flag = 1


    ! End
    stop flag
end program