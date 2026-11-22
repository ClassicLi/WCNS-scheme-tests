program euler1d
    ! Driver for the 1D Euler Riemann-problem test cases.
    use precision_select
    use scheme_params, only: scheme, read_scheme_namelist, validate_scheme_params
    use scheme_dispatch, only: scheme_output_name
    use e1d_globals, only: halo
    use e1d_control
    use e1d_initial, only: e1d_case_length, e1d_initial_state
    use e1d_solver, only: e1d_solve
    use e1d_io, only: e1d_write_profile, e1d_write_log, e1d_write_sample, &
                      e1d_count_nonfinite
    implicit none

    integer :: ierr, u_unit, n, steps, ncorrect, ncorr_interp, nan_count
    character(len=256) :: msg, csv_path, log_path, sample_path, base
    character(len=64) :: scheme_name
    character(len=32) :: op
    real(precision) :: length, h, t_final, cpu
    real(precision), allocatable :: u(:, :)

#ifdef PRECISION_DP
    character(len=*), parameter :: precision_label = 'dp'
#else
    character(len=*), parameter :: precision_label = 'qp'
#endif

    call e1d_control_init('input.nml', ierr, msg)
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

    call e1d_resolve_operator(trim(scheme), op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e1d_ensure_output_dirs(output_dir)

    n = int(mesh_num)
    length = e1d_case_length(case_1d)
    h = length/real(n, precision)

    allocate (u(-halo:n + halo, 3))
    call e1d_initial_state(n, case_1d, u)
    call e1d_solve(n, h, cfl, tend, op, u, steps, ncorrect, t_final, cpu, &
                   ncorr_interp)

    if (abs(t_final - tend) > 100.0_precision*epsilon(1.0_precision)*max(1.0_precision, tend)) then
        write (*, '(a)') 'error: solver did not reach tend'
        error stop 1
    end if

    scheme_name = trim(scheme_output_name(scheme))
    base = trim(case_1d) // '_' // trim(scheme_name) // '_' // trim(op)
    csv_path = trim(output_dir) // '/csv/' // trim(base) // '.csv'
    log_path = trim(output_dir) // '/logs/' // trim(base) // '.log'

    call e1d_write_profile(csv_path, n, length, u, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e1d_write_log(log_path, case_1d, scheme_name, op, precision_label, n, &
                       cfl, tend, steps, ncorrect, ncorr_interp, cpu, u, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    if (write_sample) then
        call execute_command_line('mkdir -p "' // trim(output_dir) // '/dat"')
        sample_path = trim(output_dir) // '/dat/' // trim(base) // '.dat'
        call e1d_write_sample(sample_path, n, length, u, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if
    end if

    nan_count = e1d_count_nonfinite(n, u)
    if (nan_count > 0) then
        write (*, '(a,i0)') 'error: non-finite cells: ', nan_count
        error stop 1
    end if

end program euler1d
