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
    real(real64) :: stiff, normal, coeff, pos, ans, f, sdelta, delta
    type(maxwell_model) :: mdl

    ! Initialization
    rst = .true.
    call random_number(stiff)
    call random_number(normal)
    call random_number(coeff)
    call random_number(pos)
    pos = pos - 0.5d0
    mdl%stiffness = stiff
    mdl%friction_coefficient = coeff

    ! Compute the actual solution
    delta = mdl%friction_coefficient / mdl%stiffness
    sdelta = min(abs(pos), delta) * sign(1.0d0, pos)
    ans = normal * mdl%stiffness * sdelta

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
function test_gmsm() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    logical :: check
    integer(int32) :: i, j
    type(generalized_maxwell_slip_model) :: mdl
    real(real64) :: f, nrm, c, x, v, muc, mus, bv, s0, alpha, vs
    real(real64) :: dzdt(3), z(3), f_ans, g, a1, a2, s
    real(real64) :: args(9)

    ! Initialization
    rst = .true.
    call mdl%initialize(3)
    call random_number(nrm)
    call random_number(c)
    call random_number(x)
    call random_number(v)
    call random_number(muc)
    call random_number(mus)
    call random_number(bv)
    call random_number(s0)
    call random_number(alpha)
    call random_number(vs)
    call random_number(args)
    call random_number(z)

    ! Ensure sum(vi) = 1
    args(3) = 1.0d0 - (args(6) + args(9))

    ! Define the model
    j = 0
    do i = 1, size(args), 3
        j = j + 1
        check = mdl%set_element_stiffness(j, args(i))
        check = mdl%set_element_damping(j, args(i+1))
        check = mdl%set_element_scaling(j, args(i+2))
    end do
    mdl%attraction_coefficient = c
    mdl%coulomb_coefficient = muc
    mdl%shape_parameter = alpha
    mdl%static_coefficient = mus
    mdl%stiffness = s0
    mdl%stribeck_velocity = vs
    mdl%viscous_damping = bv

    ! Compute the solution
    a1 = mdl%coulomb_coefficient * nrm / mdl%stiffness
    a2 = nrm * (mdl%static_coefficient - mdl%coulomb_coefficient) / &
        mdl%stiffness
    s = abs(v) / mdl%stribeck_velocity
    g = a1 + a2 / (1.0d0 + s**mdl%shape_parameter)
    f_ans = 0.0d0
    do i = 1, 3
        if (abs(z(i)) < g) then
            dzdt(i) = v
        else
            dzdt(i) = sign(1.0d0, v) * mdl%get_element_scaling(i) * &
                mdl%attraction_coefficient * &
                (1.0d0 - z(i) / (mdl%get_element_scaling(i) * g))
        end if
        f_ans = f_ans + mdl%get_element_stiffness(i) * z(i) + &
            mdl%get_element_damping(i) * dzdt(i)
    end do
    f_ans = f_ans + mdl%viscous_damping * v

    ! Test
    f = mdl%evaluate(0.0d0, x, v, nrm, z)
    if (.not.assert(f, f_ans))  then
        rst = .false.
        print *, "TEST FAILED: test_gmsm -1"
    end if
    if (.not.mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm -2"
    end if
end function

! ------------------------------------------------------------------------------
end module