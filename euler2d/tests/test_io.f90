program test_io
    ! Unit tests for the 2D Euler output layer: the Tecplot density field,
    ! key=value run logs, and the non-finite/non-positive cell counter.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use precision_select
    use e2d_globals, only: n, m, beta, ubar, e2d_allocate, e2d_deallocate, &
                           e2d_set_gama
    use e2d_io, only: e2d_write_field, e2d_write_log, e2d_count_nonfinite
    use test_util
    implicit none

    character(len=*), parameter :: field_path = 'e2d_io_test.dat'
    character(len=*), parameter :: log_path = 'e2d_io_test.log'

    integer :: ierr
    character(len=256) :: msg
    character(len=512) :: line
    integer(integer_precision) :: j, k, nx, nz
    real(precision) :: x, z, rho
    logical :: ok

    nz = 4_integer_precision
    nx = 4_integer_precision
    call e2d_set_gama(1.4_precision)
    call e2d_allocate(nx, nz)

    do k = 0, m
        do j = 0, n
            ubar(j, k, 1) = 1.0_precision + real(j + k, precision)/100.0_precision
        end do
    end do

    ! ---------------------------------------------------------------
    ! 1. Field file: header lines and one point per grid node.
    ! ---------------------------------------------------------------
    call e2d_write_field(field_path, 0.0_precision, 0.0_precision, ierr, msg)
    call check_true(ierr == 0, 'write_field status')
    call check_true(file_exists(field_path), 'field file exists')
    call check_true(line_contains(field_path, 1, 'TITLE'), 'field line 1 has TITLE')
    call check_true(line_contains(field_path, 2, 'VARIABLES'), &
                    'field line 2 has VARIABLES')
    call check_true(line_contains(field_path, 3, 'ZONE'), 'field line 3 has ZONE')
    call check_true(data_line_count(field_path) == (n + 1)*(m + 1), &
                    'field has (n+1)*(m+1) data lines')

    ! ---------------------------------------------------------------
    ! 2. Physical coordinates honour the grid origins.
    ! ---------------------------------------------------------------
    call e2d_write_field(field_path, -0.1_precision, -0.5_precision, ierr, msg)
    call check_true(ierr == 0, 'write_field shifted status')
    call first_data_line(field_path, line, ok)
    call check_true(ok, 'field data line found')
    read (line, *, iostat=ierr) x, z, rho
    call check_true(ierr == 0, 'field data line parses')
    call check_close(x, -0.1_precision, 1.0e-14_precision, 'field x origin')
    call check_close(z, -0.5_precision, 1.0e-14_precision, 'field z origin')
    call check_close(rho, 1.0_precision, 1.0e-14_precision, 'field first density')

    ! ---------------------------------------------------------------
    ! 3. Log file: key=value metadata.
    ! ---------------------------------------------------------------
    call e2d_write_log(log_path, 'riemann1', 'QR', 'HCS4', &
                       16_integer_precision, 0.3_precision, &
                       0.001_precision, 1.4_precision, &
                       10_integer_precision, 0_integer_precision, &
                       0_integer_precision, 0_integer_precision, &
                       0.5_precision, ierr, msg)
    call check_true(ierr == 0, 'write_log status')
    call check_true(file_contains(log_path, 'case=riemann1'), 'log case')
    call check_true(file_contains(log_path, 'scheme=QR'), 'log scheme')
    call check_true(file_contains(log_path, 'operator=HCS4'), 'log operator')
    call check_true(file_contains(log_path, 'mesh=16'), 'log mesh')
    call check_true(file_contains(log_path, 'cfl='), 'log cfl')
    call check_true(file_contains(log_path, 'tend='), 'log tend')
    call check_true(file_contains(log_path, 'gama='), 'log gama')
    call check_true(file_contains(log_path, 'steps=10'), 'log steps')
    call check_true(file_contains(log_path, 'interp_corrections=0'), &
                    'log interpolation corrections')
    call check_true(file_contains(log_path, 'cell_corrections=0'), &
                    'log cell corrections')
    call check_true(file_contains(log_path, 'nan_count=0'), 'log nan_count')
    call check_true(file_contains(log_path, 'cpu_seconds='), 'log cpu_seconds')

    ! ---------------------------------------------------------------
    ! 4. Non-finite counting: one NaN density, then restored.
    ! ---------------------------------------------------------------
    ubar(:, :, 1) = 1.0_precision
    ubar(:, :, 2) = 0.0_precision
    ubar(:, :, 3) = 0.0_precision
    ubar(:, :, 4) = 1.0_precision/beta
    call check_true(e2d_count_nonfinite() == 0, 'count finite state')
    ubar(0_integer_precision, 0_integer_precision, 1) = &
        ieee_value(0.0_precision, ieee_quiet_nan)
    call check_true(e2d_count_nonfinite() == 1, 'count one NaN density')
    ubar(0_integer_precision, 0_integer_precision, 1) = 1.0_precision
    call check_true(e2d_count_nonfinite() == 0, 'count restored state')

    call delete_file(field_path)
    call delete_file(log_path)
    call e2d_deallocate()
    call finish()

contains

    logical function file_exists(name)
        ! Whether the named file can be found.
        character(len=*), intent(in) :: name
        logical :: exists

        inquire (file=trim(name), exist=exists)
        file_exists = exists
    end function file_exists

    logical function file_contains(name, token)
        ! Whether any line of the file contains the token.
        character(len=*), intent(in) :: name, token
        integer :: unit, ios
        character(len=512) :: line

        file_contains = .false.
        open (newunit=unit, file=trim(name), status='old', action='read', &
              iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, trim(token)) > 0) then
                file_contains = .true.
                exit
            end if
        end do
        close (unit)
    end function file_contains

    logical function line_contains(name, which, token)
        ! Whether line `which` of the file contains the token.
        character(len=*), intent(in) :: name, token
        integer, intent(in) :: which
        integer :: unit, ios, count
        character(len=512) :: line

        line_contains = .false.
        open (newunit=unit, file=trim(name), status='old', action='read', &
              iostat=ios)
        if (ios /= 0) return
        count = 0
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            count = count + 1
            if (count == which) then
                line_contains = index(line, trim(token)) > 0
                exit
            end if
        end do
        close (unit)
    end function line_contains

    integer function data_line_count(name)
        ! Number of point lines that follow the ZONE header.
        character(len=*), intent(in) :: name
        integer :: unit, ios
        character(len=512) :: line
        logical :: in_data

        data_line_count = 0
        in_data = .false.
        open (newunit=unit, file=trim(name), status='old', action='read', &
              iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (in_data) then
                data_line_count = data_line_count + 1
            else if (index(line, 'ZONE') > 0) then
                in_data = .true.
            end if
        end do
        close (unit)
    end function data_line_count

    subroutine first_data_line(name, line, ok)
        ! Return the first line after the ZONE header.
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: line
        logical, intent(out) :: ok
        integer :: unit, ios
        character(len=512) :: buffer
        logical :: in_data

        line = ''
        ok = .false.
        in_data = .false.
        open (newunit=unit, file=trim(name), status='old', action='read', &
              iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(a)', iostat=ios) buffer
            if (ios /= 0) exit
            if (in_data) then
                line = buffer
                ok = .true.
                exit
            else if (index(buffer, 'ZONE') > 0) then
                in_data = .true.
            end if
        end do
        close (unit)
    end subroutine first_data_line

    subroutine delete_file(name)
        ! Remove the named file if it exists.
        character(len=*), intent(in) :: name
        integer :: unit, ios
        logical :: exists

        inquire (file=trim(name), exist=exists)
        if (.not. exists) return
        open (newunit=unit, file=trim(name), status='old', iostat=ios)
        if (ios /= 0) return
        close (unit, status='delete')
    end subroutine delete_file

end program test_io
