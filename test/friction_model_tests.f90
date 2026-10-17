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
    use friction_data_io
    use fortran_test_helper
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    implicit none
contains
! ------------------------------------------------------------------------------
function test_coulomb() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: normal, f
    type(coulomb_model) :: mdl

    ! Initialization
    rst = .true.
    normal = 12.0d0
    mdl%friction_coefficient = 0.35d0

    ! Test
    f = mdl%evaluate(0.0d0, 0.0d0, 0.2d0, normal)
    if (.not.assert(f, 4.2d0)) then
        rst = .false.
        print *, "TEST FAILED: test_coulomb -1"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, -0.2d0, normal)
    if (.not.assert(f, -4.2d0)) then
        rst = .false.
        print *, "TEST FAILED: test_coulomb -2"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, normal)
    if (.not.assert(f, 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_coulomb -3"
    end if
    if (mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_coulomb -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_lugre() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    integer(int32) :: i
    real(real64) :: mus, mud, k, b, bv, vs, a, s, dsdt(1), v, n, f, &
        fans, dsans, g, a1, a2, Fc, Fs
    real(real64), parameter, dimension(2) :: velocities = &
        [0.2d0, -0.2d0]
    type(lugre_model) :: mdl

    ! Initialization
    rst = .true.
    mus = 0.8d0
    mud = 0.4d0
    k = 100.0d0
    b = 1.5d0
    bv = 0.2d0
    vs = 0.5d0
    a = 2.0d0
    s = 0.01d0
    n = 10.0d0
    mdl%static_coefficient = mus
    mdl%coulomb_coefficient = mud
    mdl%stribeck_velocity = vs
    mdl%shape_parameter = a
    mdl%stiffness = k
    mdl%damping = b
    mdl%viscous_damping = bv

    do i = 1, size(velocities)
        v = velocities(i)
        Fc = mdl%coulomb_coefficient * n
        Fs = mdl%static_coefficient * n
        a1 = Fc / mdl%stiffness
        a2 = (Fs - Fc) / mdl%stiffness
        g = a1 + a2 / &
            (1.0d0 + (abs(v) / mdl%stribeck_velocity)**mdl%shape_parameter)
        dsans = v - abs(v) * s / g
        fans = k * s + b * dsans + bv * v

        call mdl%state(0.0d0, 0.0d0, v, n, [s], dsdt)
        f = mdl%evaluate(0.0d0, 0.0d0, v, n, [s])
        if (.not.assert(dsans, dsdt(1))) then
            rst = .false.
            print *, "TEST FAILED: test_lugre state", i
        end if
        if (.not.assert(fans, f)) then
            rst = .false.
            print *, "TEST FAILED: test_lugre force", i
        end if
    end do
    v = 0.0d0
    call mdl%state(0.0d0, 0.0d0, v, n, [s], dsdt)
    if (.not.assert(dsdt(1), 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_lugre zero velocity"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, v, n, [s])
    if (.not.assert(f, k * s)) then
        rst = .false.
        print *, "TEST FAILED: test_lugre zero-velocity force"
    end if
    if (.not.mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_lugre state flag"
    end if
end function

! ------------------------------------------------------------------------------
function test_maxwell() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: f, normal
    type(maxwell_model) :: mdl

    ! Initialization
    rst = .true.
    normal = 10.0d0
    mdl%stiffness = 100.0d0
    mdl%friction_coefficient = 0.3d0

    ! Test
    f = mdl%evaluate(0.0d0, 0.01d0, 0.0d0, normal)
    if (.not.assert(f, 1.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell elastic loading"
    end if
    f = mdl%evaluate(0.0d0, 0.2d0, 0.0d0, normal)
    if (.not.assert(f, 3.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell positive saturation"
    end if
    f = mdl%evaluate(0.0d0, 0.19d0, 0.0d0, normal)
    if (.not.assert(f, 2.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell unloading"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, normal)
    if (.not.assert(f, -3.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell reversal"
    end if
    call mdl%reset()
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, normal)
    if (.not.assert(f, 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell reset"
    end if
    if (mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_maxwell state flag"
    end if
end function

! ------------------------------------------------------------------------------
function test_stribeck() result(rst)
    logical :: rst
    real(real64) :: nrm, dry_force, f, velocity
    type(stribeck_model) :: mdl

    rst = .true.
    nrm = 10.0d0
    velocity = 0.25d0
    mdl%static_friction_coefficient = 0.8d0
    mdl%coulomb_friction_coefficient = 0.4d0
    mdl%stribeck_velocity = 0.5d0
    mdl%viscous_damping = 2.0d0
    dry_force = nrm * (mdl%coulomb_friction_coefficient + &
        (mdl%static_friction_coefficient - &
        mdl%coulomb_friction_coefficient) * &
        exp(-(velocity / mdl%stribeck_velocity)**2))

    f = mdl%evaluate(0.0d0, 0.0d0, velocity, nrm)
    if (.not.assert(f, dry_force + mdl%viscous_damping * velocity)) then
        rst = .false.
        print *, "TEST FAILED: test_stribeck positive velocity"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, -velocity, nrm)
    if (.not.assert(f, -dry_force - mdl%viscous_damping * velocity)) then
        rst = .false.
        print *, "TEST FAILED: test_stribeck negative velocity"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, nrm)
    if (.not.assert(f, 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_stribeck zero velocity"
    end if
    if (mdl%has_internal_state()) then
        rst = .false.
        print *, "TEST FAILED: test_stribeck state flag"
    end if
end function

! ------------------------------------------------------------------------------
function test_modified_stribeck() result(rst)
    logical :: rst
    real(real64) :: f, nrm
    type(modified_stribeck_model) :: mdl

    rst = .true.
    nrm = 10.0d0
    mdl%static_friction_coefficient = 0.8d0
    mdl%coulomb_friction_coefficient = 0.4d0
    mdl%stribeck_velocity = 0.5d0
    mdl%viscous_damping = 2.0d0
    mdl%stiffness = 100.0d0

    f = mdl%evaluate(0.0d0, 0.03d0, 0.0d0, nrm)
    if (.not.assert(f, 3.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck presliding"
    end if
    f = mdl%evaluate(0.0d0, 0.2d0, 0.0d0, nrm)
    if (.not.assert(f, 8.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck positive saturation"
    end if
    f = mdl%evaluate(0.0d0, 0.18d0, 0.0d0, nrm)
    if (.not.assert(f, 6.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck unloading"
    end if
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, nrm)
    if (.not.assert(f, -8.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck reversal"
    end if
    call mdl%reset()
    f = mdl%evaluate(0.0d0, 0.0d0, 0.25d0, nrm)
    if (.not.assert(f, mdl%viscous_damping * 0.25d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck viscous force"
    end if
    call mdl%reset()
    f = mdl%evaluate(0.0d0, 0.0d0, -0.25d0, nrm)
    if (.not.assert(f, -mdl%viscous_damping * 0.25d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck reverse viscous force"
    end if
    call mdl%reset()
    f = mdl%evaluate(0.0d0, 0.0d0, 0.0d0, nrm)
    if (.not.assert(f, 0.0d0)) then
        rst = .false.
        print *, "TEST FAILED: test_modified_stribeck reset"
    end if
end function

! ------------------------------------------------------------------------------
function test_friction_data_io() result(rst)
    logical :: rst
    integer(int32) :: i, unit, ios
    character(len=128) :: header_line
    type(friction_data) :: data, roundtrip, semicolon_data
    character(len=*), parameter :: input_file = 'friction_data_io_input.csv'
    character(len=*), parameter :: output_file = 'friction_data_io_output.csv'
    character(len=*), parameter :: semicolon_file = &
        'friction_data_io_semicolon.csv'
    character(len=*), parameter :: large_file = 'friction_data_io_large.csv'
    character(len=*), parameter :: custom_header = &
        't|x|v|normal|friction'

    rst = .true.
    open(newunit=unit, file=input_file, status='replace', action='write', &
        iostat=ios)
    if (ios /= 0) then
        rst = .false.
        print *, "TEST FAILED: test_friction_data_io fixture open"
        return
    end if
    write(unit, '(A)') 'Time,Position,Velocity,Normal Force,Friction Force'
    write(unit, '(A)') '0.0,1.0,,3.0,4.0'
    write(unit, '(A)') ''
    write(unit, '(A)') '0.1,2.0,3.0,,5.0'
    write(unit, '(A)') '0.2,2.5,3.5,4.5'
    write(unit, '(A)') '0.3,invalid,3.6,4.6,5.6'
    close(unit)

    call read_friction_data(input_file, data, io_status=ios)
    if (ios /= 0) then
        rst = .false.
        print *, "TEST FAILED: test_friction_data_io read CSV"
    else
        if (size(data%time) /= 4) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io blank row handling"
        else
            if (.not.assert(data%time(4), 0.3d0) .or. &
                .not.ieee_is_nan(data%velocity(1)) .or. &
                .not.ieee_is_nan(data%normal_force(2)) .or. &
                .not.ieee_is_nan(data%friction_force(3)) .or. &
                .not.ieee_is_nan(data%position(4))) then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io missing data"
            end if
        end if
    end if

    call write_friction_data(output_file, data, io_status=ios)
    if (ios /= 0) then
        rst = .false.
        print *, "TEST FAILED: test_friction_data_io write CSV"
    else
        open(newunit=unit, file=output_file, status='old', action='read', &
            iostat=ios)
        if (ios /= 0) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io output open"
        else
            read(unit, '(A)', iostat=ios) header_line
            close(unit)
            if (ios /= 0 .or. trim(header_line) /= &
                '"Time","Position","Velocity","Normal Force","Friction Force"') then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io default header"
            end if
        end if
        call read_friction_data(output_file, roundtrip, io_status=ios)
        if (ios /= 0) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io round trip read"
        else if (size(roundtrip%time) /= size(data%time)) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io round trip rows"
        else if (.not.ieee_is_nan(roundtrip%velocity(1)) .or. &
            .not.assert(roundtrip%time(2), data%time(2))) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io round trip values"
        end if
    end if

    open(newunit=unit, file=semicolon_file, status='replace', &
        action='write', iostat=ios)
    if (ios /= 0) then
        rst = .false.
        print *, "TEST FAILED: test_friction_data_io semicolon fixture"
    else
        write(unit, '(A)') '1;2;3;4;5'
        write(unit, '(A)') ''
        write(unit, '(A)') '6;;8;9;10'
        close(unit)
        call read_friction_data(semicolon_file, semicolon_data, &
            delimiter=';', has_header=.false., io_status=ios)
        if (ios /= 0) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io headerless read"
        else if (size(semicolon_data%time) /= 2) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io headerless rows"
        else if (.not.assert(semicolon_data%friction_force(1), 5.0d0) .or. &
            .not.ieee_is_nan(semicolon_data%position(2))) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io custom delimiter"
        end if

        call write_friction_data(output_file, semicolon_data, &
            delimiter='|', header=custom_header, io_status=ios)
        if (ios /= 0) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io custom write"
        else
            open(newunit=unit, file=output_file, status='old', &
                action='read', iostat=ios)
            if (ios /= 0) then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io custom output open"
            else
                read(unit, '(A)', iostat=ios) header_line
                close(unit)
                if (ios /= 0 .or. trim(header_line) /= custom_header) then
                    rst = .false.
                    print *, "TEST FAILED: test_friction_data_io custom header"
                end if
            end if
            call read_friction_data(output_file, roundtrip, delimiter='|', &
                io_status=ios)
            if (ios /= 0) then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io custom round trip"
            else if (size(roundtrip%time) /= 2) then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io custom row count"
            else if (.not.assert(roundtrip%time(2), 6.0d0) .or. &
                .not.ieee_is_nan(roundtrip%position(2)) .or. &
                .not.assert(roundtrip%friction_force(2), 10.0d0)) then
                rst = .false.
                print *, "TEST FAILED: test_friction_data_io custom values"
            end if
        end if
    end if

    open(newunit=unit, file=large_file, status='replace', action='write', &
        iostat=ios)
    if (ios /= 0) then
        rst = .false.
        print *, "TEST FAILED: test_friction_data_io large fixture"
    else
        do i = 1, 2500
            write(unit, '(I0,A,I0,A,I0,A,I0,A,I0)') i, ',', i + 1, ',', &
                i + 2, ',', i + 3, ',', i + 4
        end do
        close(unit)
        call read_friction_data(large_file, roundtrip, has_header=.false., &
            io_status=ios)
        if (ios /= 0) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io large read"
        else if (size(roundtrip%time) /= 2500) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io storage growth count"
        else if (.not.assert(roundtrip%time(2500), 2500.0d0) .or. &
            .not.assert(roundtrip%friction_force(2500), 2504.0d0)) then
            rst = .false.
            print *, "TEST FAILED: test_friction_data_io storage growth values"
        end if
    end if

    call delete_test_file(input_file)
    call delete_test_file(output_file)
    call delete_test_file(semicolon_file)
    call delete_test_file(large_file)
end function

! ------------------------------------------------------------------------------
subroutine delete_test_file(filename)
    character(len=*), intent(in) :: filename
    integer(int32) :: unit, ios

    open(newunit=unit, file=filename, status='old', iostat=ios)
    if (ios == 0) close(unit, status='delete')
end subroutine

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