! MIT License
!
! Copyright (c) 2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
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
    delta = normal * mdl%friction_coefficient / mdl%stiffness
    sdelta = min(abs(pos), delta) * sign(1.0d0, pos)
    ans = mdl%stiffness * sdelta

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
    type(generalized_maxwell_slip_model) :: mdl_copy
    real(real64) :: f, nrm, x, v, muc, mus, bv, vs, splus
    real(real64), dimension(3) :: dzdt, z, z_scale, expected_rate
    real(real64) :: expected_force, normalized_state, eta_a, eta_b
    real(real64), dimension(16) :: params
    real(real64), parameter, dimension(2) :: velocities = &
        [0.2d0, -0.2d0]

    ! Initialization
    rst = .true.
    call mdl%initialize(3)
    nrm = 10.0d0
    x = 0.0d0
    muc = 0.4d0
    mus = 0.8d0
    bv = 0.15d0
    vs = 0.5d0

    ! Define the model
    do i = 1, 3
        check = mdl%set_element_stiffness(i, real(i + 1, real64))
        if (.not.check) then
            rst = .false.
            return
        end if
        check = mdl%set_element_damping(i, 0.05d0 * real(i, real64))
        if (.not.check) then
            rst = .false.
            return
        end if
        check = mdl%set_element_scaling(i, 0.1d0 * real(i, real64))
        if (.not.check) then
            rst = .false.
            return
        end if
    end do
    mdl%coulomb_coefficient = muc
    mdl%static_coefficient = mus
    mdl%stribeck_velocity = vs
    mdl%viscous_damping = bv
    mdl%transition_sharpness = 100.0d0
    mdl%sliding_margin = 0.95d0
    mdl%reversal_sharpness = 10.0d0

    splus = nrm * (muc + (mus - muc) * &
        exp(-(velocities(1) / vs)**2))
    do i = 1, 3
        z_scale(i) = mdl%get_element_scaling(i) * splus / &
            mdl%get_element_stiffness(i)
    end do
    z = [0.4d0 * z_scale(1), 0.97d0 * z_scale(2), &
        -0.5d0 * z_scale(3)]

    do j = 1, size(velocities)
        v = velocities(j)
        expected_force = 0.0d0
        do i = 1, 3
            normalized_state = z(i) / z_scale(i)
            eta_a = 1.0d0 - 0.5d0 * tanh(mdl%transition_sharpness * &
                (normalized_state + mdl%sliding_margin)) + 0.5d0 * &
                tanh(mdl%transition_sharpness * &
                (normalized_state - mdl%sliding_margin))
            eta_b = 0.5d0 + 0.5d0 * tanh(mdl%reversal_sharpness * &
                normalized_state * v / vs)
            expected_rate(i) = v - eta_a * eta_b * abs(v) * &
                normalized_state
            expected_force = expected_force + &
                mdl%get_element_stiffness(i) * z(i) + &
                mdl%get_element_damping(i) * expected_rate(i)
        end do
        expected_force = expected_force + bv * v

        call mdl%state(0.0d0, x, v, nrm, z, dzdt)
        f = mdl%evaluate(0.0d0, x, v, nrm, z)
        do i = 1, 3
            if (.not.assert(dzdt(i), expected_rate(i))) then
                rst = .false.
                print *, "TEST FAILED: test_gmsm state rate", j, i
            end if
        end do
        if (.not.assert(f, expected_force)) then
            rst = .false.
            print *, "TEST FAILED: test_gmsm force", j
        end if
    end do

    ! A velocity reversal resets the slipping weight toward presliding.
    v = -velocities(1)
    dzdt(1) = mdl%element_state(1, 0.0d0, x, v, nrm, z_scale(1))
    if (abs(dzdt(1) - v) > 1.0d-3) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm reversal"
    end if

    ! The smoothed state rate has no jump at the presliding/sliding boundary.
    v = velocities(1)
    dzdt(1) = mdl%element_state(1, 0.0d0, x, v, nrm, &
        z_scale(1) * (mdl%sliding_margin - 1.0d-8))
    expected_rate(1) = mdl%element_state(1, 0.0d0, x, v, nrm, &
        z_scale(1) * (mdl%sliding_margin + 1.0d-8))
    if (abs(dzdt(1) - expected_rate(1)) > 1.0d-5) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm smooth transition"
    end if

    ! A stationary input must not move any internal state, even outside z_s.
    v = 0.0d0
    z = 2.0d0 * z_scale
    call mdl%state(0.0d0, x, v, nrm, z, dzdt)
    if (any(dzdt /= 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm zero velocity"
    end if

    ! Check the public parameter-array mapping for the new S-GMS parameters.
    call mdl%to_array(params)
    call mdl_copy%initialize(3)
    call mdl_copy%from_array(params)
    if (.not.assert(mdl_copy%transition_sharpness, &
        mdl%transition_sharpness) .or. &
        .not.assert(mdl_copy%sliding_margin, mdl%sliding_margin) .or. &
        .not.assert(mdl_copy%reversal_sharpness, &
        mdl%reversal_sharpness)) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm parameter mapping"
    end if

    if (.not.mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_gmsm -2"
    end if
end function

! ------------------------------------------------------------------------------
end module