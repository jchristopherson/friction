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
        print *, "TEST FAILED: test_coulomb -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_lugre() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: mus, mud, k, b, bv, vs, a, s, dsdt(1), v, n, f, &
        fans, dsans, g, a1, a2, Fc, Fs
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
    Fc = mdl%coulomb_coefficient * n
    Fs = mdl%static_coefficient * n
    a1 = Fc / mdl%stiffness
    a2 = (Fs - Fc) / mdl%stiffness
    g = a1 + a2 / (1.0d0 + (abs(v) / mdl%stribeck_velocity)**mdl%shape_parameter)
    dsans = v - abs(v) * s / g
    fans = k * s + b * dsans + bv * v

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
function test_maxwell() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: stiff, normal, coeff, pos, ans, f, s, delta
    type(maxwell_model) :: mdl

    ! Initialization
    rst = .true.
    call random_number(stiff)
    call random_number(normal)
    call random_number(coeff)
    call random_number(pos)
    pos = pos - 0.5d0
    mdl%stiffness = 1.0d3 * stiff
    mdl%friction_coefficient = coeff

    ! Compute the actual solution
    delta = mdl%friction_coefficient / mdl%stiffness
    s = min(abs(pos), delta)
    ans = normal * mdl%stiffness * min(abs(s), delta) * sign(1.0d0, s)

    ! Test
    f = mdl%evaluate(0.0d0, pos, 0.0d0, normal)
    if (.not.assert(f, ans)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell -1"
    end if
    if (mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell -2"
    end if
end function

! ------------------------------------------------------------------------------
end module