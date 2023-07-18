module friction_model_tests
    use iso_fortran_env
    use friction
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_coulomb() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: normal, coeff, vel, ans, f
    type(coulomb_model) :: mdl

    ! Initialization
    rst = .true.
    call random_number(normal)
    call random_number(coeff)
    call random_number(vel)
    vel = vel - 0.5d0
    mdl%friction_coefficient = coeff

    ! Compute the actual solution
    if (vel == 0.0d0) then
        ans = 0.0d0
    else
        ans = coeff * normal * sign(1.0d0, vel)
    end if

    ! Test
    f = mdl%evaluate(0.0d0, 0.0d0, vel, normal)
    if (.not.assert(f, ans)) then
        rst = .false.
        print *, "TEST FAILED: test_coulomb -1"
    end if
    if (mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_lugre -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_lugre() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: mus, mud, k, b, bv, vs, a, s, dsdt(1), v, n, f, &
        fans, dsans, g
    type(lugre_model) :: mdl

    ! Initialization
    rst = .true.
    call random_number(mus)
    call random_number(mud)
    call random_number(k)
    call random_number(b)
    call random_number(bv)
    call random_number(vs)
    call random_number(a)
    call random_number(s)
    call random_number(v)
    call random_number(n)
    v = v - 0.5d0
    mdl%static_coefficient = mus
    mdl%coulomb_coefficient = mud
    mdl%stribeck_velocity = vs
    mdl%shape_parameter = a
    mdl%stiffness = k
    mdl%damping = b
    mdl%viscous_damping = bv

    ! Compute the solution
    g = mud + (mus - mud) * exp(-abs(v / vs)**a)
    dsans = v - k * abs(v) * s / g
    fans = n * (k * s + b * dsans + bv * v)

    call mdl%state(0.0d0, 0.0d0, v, n, [s], dsdt)
    f = mdl%evaluate(0.0d0, 0.0d0, v, n, [s])

    ! Test
    if (.not.assert(dsans, dsdt(1))) then
        rst = .false.
        print *, "TEST FAILED: test_lugre -1"
    end if
    if (.not.assert(fans, f)) then
        rst = .false.
        print *, "TEST FAILED: test_lugre -2"
    end if
    if (.not.mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_lugre -3"
    end if
end function

! ------------------------------------------------------------------------------
end module