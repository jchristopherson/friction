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

    check = test_lugre()
    if (.not.check) flag = 2

    check = test_maxwell()
    if (.not.check) flag = 3

    check = test_gmsm()
    if (.not.check) flag = 4
    
    ! End
    stop flag
end program