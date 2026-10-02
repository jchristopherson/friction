module friction_data_io
    use iso_fortran_env
    use friction_data_handling
    use friction_errors
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: friction_data
    public :: read_friction_data
    public :: write_friction_data

    character(len=*), parameter :: DEFAULT_HEADER = &
        '"Time","Position","Velocity","Normal Force","Friction Force"'
    integer(int32), parameter :: DATA_COLUMN_COUNT = 5
    integer(int32), parameter :: RECORD_LENGTH = 65536
contains
! ------------------------------------------------------------------------------
subroutine read_friction_data(filename, data, delimiter, has_header, io_status)
    !! Reads time-series data from a delimited text file into a friction_data
    !! container. The five fields are interpreted in this order: time,
    !! position, velocity, normal force, and friction force. Blank records are
    !! ignored; empty or nonnumeric fields and omitted trailing fields are set
    !! to quiet NaN. Fields after the first five are ignored. A numeric field
    !! may be surrounded by double quotes and whitespace.
    !!
    !! By default, the first nonblank record is treated as a header and the
    !! delimiter is a comma. The optional io_status argument allows callers to
    !! handle file and allocation errors; malformed numeric fields are treated
    !! as missing data and do not set an error status. If io_status is omitted,
    !! an I/O or allocation error terminates execution with
    !! FRICTION_INVALID_OPERATION_ERROR.
    !!
    !! @param filename The path of the input text file.
    !! @param data The output container. Its five allocatable arrays are
    !! allocated to the number of nonblank data records read.
    !! @param delimiter Optional single-character field delimiter; defaults to
    !! a comma.
    !! @param has_header Optional flag indicating whether the first nonblank
    !! record is a header; defaults to true.
    !! @param io_status Optional output status. It is zero on success, receives
    !! the native I/O status on file errors, or FRICTION_MEMORY_ERROR on
    !! allocation failure.
    character(len=*), intent(in) :: filename
        !! The path of the input text file.
    type(friction_data), intent(out) :: data
        !! The output friction-data container.
    character(len=1), intent(in), optional :: delimiter
        !! The field delimiter; defaults to a comma.
    logical, intent(in), optional :: has_header
        !! True if the first nonblank record is a header; defaults to true.
    integer(int32), intent(out), optional :: io_status
        !! Optional error status; zero indicates success.

    integer(int32) :: i, row_count, unit, ios, close_status
    real(real64), dimension(DATA_COLUMN_COUNT) :: values
    character(len=RECORD_LENGTH) :: line
    character(len=1) :: separator
    logical :: skip_header, first_record

    separator = ','
    if (present(delimiter)) separator = delimiter
    skip_header = .true.
    if (present(has_header)) skip_header = has_header
    if (present(io_status)) io_status = 0

    open(newunit=unit, file=trim(filename), status='old', action='read', &
        iostat=ios)
    if (ios /= 0) then
        call report_io_error(ios, io_status)
        return
    end if

    row_count = 0
    first_record = .true.
    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) then
            close(unit)
            call report_io_error(ios, io_status)
            return
        end if
        if (len_trim(line) == 0) cycle
        if (first_record) then
            first_record = .false.
            if (skip_header) cycle
        end if
        row_count = row_count + 1
    end do

    rewind(unit, iostat=ios)
    if (ios /= 0) then
        close(unit)
        call report_io_error(ios, io_status)
        return
    end if

    allocate(data%time(row_count), data%position(row_count), &
        data%velocity(row_count), data%normal_force(row_count), &
        data%friction_force(row_count), stat=ios)
    if (ios /= 0) then
        close(unit)
        call report_io_error(FRICTION_MEMORY_ERROR, io_status)
        return
    end if

    first_record = .true.
    i = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) then
            close(unit)
            call report_io_error(ios, io_status)
            return
        end if
        if (len_trim(line) == 0) cycle
        if (first_record) then
            first_record = .false.
            if (skip_header) cycle
        end if

        call parse_data_row(trim(line), separator, values)
        i = i + 1
        data%time(i) = values(1)
        data%position(i) = values(2)
        data%velocity(i) = values(3)
        data%normal_force(i) = values(4)
        data%friction_force(i) = values(5)
    end do
    close(unit, iostat=close_status)
    if (close_status /= 0) call report_io_error(close_status, io_status)
end subroutine

! ------------------------------------------------------------------------------
subroutine write_friction_data(filename, data, delimiter, header, io_status)
    !! Writes the five time-series arrays in a friction_data container to a
    !! delimited text file. The arrays are written in this order: time,
    !! position, velocity, normal force, and friction force. The file is
    !! replaced if it already exists, and the header is always written as the
    !! first record.
    !!
    !! The delimiter defaults to a comma. If header is omitted, the first
    !! record is exactly:
    !! "Time","Position","Velocity","Normal Force","Friction Force"
    !! A supplied header is written verbatim as one record. Numeric values are
    !! written in scientific notation with double-precision output precision.
    !! All five arrays must be allocated and have the same length.
    !!
    !! The optional io_status argument allows callers to handle file and
    !! container errors. It is zero on success, receives the I/O status on
    !! file errors, FRICTION_INVALID_OPERATION_ERROR if an array is unallocated,
    !! or FRICTION_ARRAY_SIZE_ERROR if the array lengths differ. If io_status is
    !! omitted, an error terminates execution with
    !! FRICTION_INVALID_OPERATION_ERROR.
    !!
    !! @param filename The path of the output text file.
    !! @param data The friction-data container to write.
    !! @param delimiter Optional single-character field delimiter; defaults to
    !! a comma.
    !! @param header Optional header record written verbatim; defaults to the
    !! standard five-column header.
    !! @param io_status Optional output status; zero indicates success.
    character(len=*), intent(in) :: filename
        !! The path of the output text file.
    type(friction_data), intent(in) :: data
        !! The friction-data container to write.
    character(len=1), intent(in), optional :: delimiter
        !! The field delimiter; defaults to a comma.
    character(len=*), intent(in), optional :: header
        !! The header record; defaults to the standard five-column header.
    integer(int32), intent(out), optional :: io_status
        !! Optional error status; zero indicates success.

    integer(int32) :: i, unit, ios, n
    character(len=1) :: separator
    character(len=:), allocatable :: header_line

    if (present(io_status)) io_status = 0
    if (.not.allocated(data%time) .or. .not.allocated(data%position) .or. &
        .not.allocated(data%velocity) .or. &
        .not.allocated(data%normal_force) .or. &
        .not.allocated(data%friction_force)) then
        call report_io_error(FRICTION_INVALID_OPERATION_ERROR, io_status)
        return
    end if

    n = size(data%time)
    if (size(data%position) /= n .or. size(data%velocity) /= n .or. &
        size(data%normal_force) /= n .or. &
        size(data%friction_force) /= n) then
        call report_io_error(FRICTION_ARRAY_SIZE_ERROR, io_status)
        return
    end if

    separator = ','
    if (present(delimiter)) separator = delimiter
    header_line = DEFAULT_HEADER
    if (present(header)) header_line = header

    open(newunit=unit, file=trim(filename), status='replace', &
        action='write', iostat=ios)
    if (ios /= 0) then
        call report_io_error(ios, io_status)
        return
    end if

    write(unit, '(A)', iostat=ios) header_line
    if (ios /= 0) then
        close(unit)
        call report_io_error(ios, io_status)
        return
    end if

    do i = 1, n
        write(unit, '(ES24.16E3,4(A,ES24.16E3))', iostat=ios) &
            data%time(i), separator, data%position(i), separator, &
            data%velocity(i), separator, data%normal_force(i), separator, &
            data%friction_force(i)
        if (ios /= 0) then
            close(unit)
            call report_io_error(ios, io_status)
            return
        end if
    end do

    close(unit, iostat=ios)
    if (ios /= 0) call report_io_error(ios, io_status)
end subroutine

! ------------------------------------------------------------------------------
subroutine parse_data_row(line, delimiter, values)
    character(len=*), intent(in) :: line
    character(len=1), intent(in) :: delimiter
    real(real64), intent(out), dimension(DATA_COLUMN_COUNT) :: values

    integer(int32) :: i, first, column, last

    values = quiet_nan()
    first = 1
    column = 1
    last = len(line)
    do i = 1, last
        if (line(i:i) /= delimiter) cycle
        if (column <= DATA_COLUMN_COUNT) then
            values(column) = parse_real_field(line(first:i-1))
        end if
        column = column + 1
        first = i + 1
    end do
    if (column <= DATA_COLUMN_COUNT .and. first <= last) then
        values(column) = parse_real_field(line(first:last))
    end if
end subroutine

! ------------------------------------------------------------------------------
function parse_real_field(field) result(value)
    character(len=*), intent(in) :: field
    real(real64) :: value

    integer(int32) :: ios, n
    character(len=:), allocatable :: token

    token = trim(adjustl(field))
    n = len(token)
    if (n >= 2) then
        if (token(1:1) == '"' .and. token(n:n) == '"') then
            token = trim(adjustl(token(2:n-1)))
        end if
    end if
    if (len(token) == 0) then
        value = quiet_nan()
        return
    end if

    read(token, *, iostat=ios) value
    if (ios /= 0) value = quiet_nan()
end function

! ------------------------------------------------------------------------------
pure function quiet_nan() result(value)
    real(real64) :: value
    value = ieee_value(0.0d0, ieee_quiet_nan)
end function

! ------------------------------------------------------------------------------
subroutine report_io_error(error_code, io_status)
    integer(int32), intent(in) :: error_code
    integer(int32), intent(out), optional :: io_status

    if (present(io_status)) then
        io_status = error_code
    else
        error stop FRICTION_INVALID_OPERATION_ERROR
    end if
end subroutine
end module