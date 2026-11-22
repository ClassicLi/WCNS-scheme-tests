program accuracy_2d
    ! Driver for the 2D Euler accuracy tests.
    !
    ! Runs the selected case (vortex or sine) on nbase*2^(k-1) grids and
    ! writes the density error/order table.
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use scheme_params, only: scheme, read_scheme_namelist, validate_scheme_params
    use scheme_dispatch, only: scheme_output_name
    use accuracy_control
    use accuracy_order
    use e2d_globals, only: n, m, ubar, e2d_allocate, e2d_deallocate
    use e2d_cases, only: e2d_case_length, e2d_initial, e2d_exact_density
    use e2d_solver, only: e2d_boundary, e2d_solve_level
    use e2d_io, only: e2d_write_sample
    implicit none

    integer :: ierr, u_unit, k, nx
    character(len=256) :: msg
    character(len=32) :: op
    character(len=64) :: scheme_name
    character(len=256) :: out_name
    character(len=512) :: csv_path, log_path, sample_path
    real(precision) :: length, cpu, e1, e2, einf
    real(precision), allocatable :: rho(:, :)

#ifdef PRECISION_DP
    character(len=*), parameter :: precision_label = 'dp'
#else
    character(len=*), parameter :: precision_label = 'qp'
#endif

    ! ---------------------------------------------------------------
    ! Configuration
    ! ---------------------------------------------------------------
    call accuracy_control_init_2d('input.nml', ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    open (newunit=u_unit, file='input.nml', status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write (*, '(a)') 'error: cannot open input.nml'
        error stop 1
    end if
    call read_scheme_namelist(u_unit, ierr)
    close (u_unit)
    if (ierr /= 0) then
        write (*, '(a)') 'error: cannot read &scheme_nml'
        error stop 1
    end if

    call validate_scheme_params(ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call accuracy_resolve_operator(scheme, op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e2d_case_length(case_2d, length, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call accuracy_ensure_output_dirs(output_dir)

    scheme_name = trim(scheme_output_name(scheme))
    out_name = 'accuracy_2d_' // trim(case_2d) // '_' // trim(scheme_name) // '_' // trim(op)

    ! ---------------------------------------------------------------
    ! Refinement loop
    ! ---------------------------------------------------------------
    call order_reset()

    do k = 1, int(levels, 4)
        nx = int(nbase, 4)*2**(k - 1)

        call e2d_allocate(nx, length)
        call e2d_initial(case_2d, sine_rho0, sine_amp, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if
        call e2d_boundary()

        call e2d_solve_level(op, cfl, tend, cpu)

        allocate (rho(0:n, 0:m))
        call e2d_exact_density(case_2d, tend, sine_rho0, sine_amp, rho, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if

        call error_norms_2d(n, m, ubar(0:n, 0:m, 1), rho, e1, e2, einf)

        if (.not. ieee_is_finite(e1) .or. .not. ieee_is_finite(e2) .or. &
            .not. ieee_is_finite(einf)) then
            write (*, '(a,i0)') 'error: non-finite error at n = ', n
            error stop 1
        end if

        call order_record(k, int(n, integer_precision), e1, e2, einf, cpu)

        if (write_sample .and. k == int(levels, 4)) then
            sample_path = trim(output_dir) // '/figures/' // trim(out_name) // '.dat'
            call e2d_write_sample(sample_path)
        end if

        deallocate (rho)
        call e2d_deallocate()
    end do

    ! ---------------------------------------------------------------
    ! Output
    ! ---------------------------------------------------------------
    csv_path = trim(output_dir) // '/csv/' // trim(out_name) // '.csv'
    log_path = trim(output_dir) // '/logs/' // trim(out_name) // '.log'

    call order_write_csv(csv_path, scheme_name, op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call write_log(log_path, precision_label)

contains

    subroutine write_log(path, prec)
        ! Write run metadata as key=value lines.
        character(len=*), intent(in) :: path, prec

        integer :: unit

        open (newunit=unit, file=trim(path), status='replace', action='write')
        write (unit, '(a)') 'case=' // trim(case_2d)
        write (unit, '(a)') 'scheme=' // trim(scheme_name)
        write (unit, '(a)') 'operator=' // trim(op)
        write (unit, '(a)') 'precision=' // trim(prec)
        write (unit, '(a,i0)') 'nbase=', nbase
        write (unit, '(a,i0)') 'levels=', levels
        write (unit, '(a,es24.16)') 'cfl=', cfl
        write (unit, '(a,es24.16)') 'tend=', tend
        write (unit, '(a,es24.16)') 'sine_rho0=', sine_rho0
        write (unit, '(a,es24.16)') 'sine_amp=', sine_amp
        close (unit)
    end subroutine write_log

end program accuracy_2d
