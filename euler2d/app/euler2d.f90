program euler2d
    ! Driver for the seven 2D Euler test cases.
    !
    ! Reads input.nml, allocates the grid, applies the initial and boundary
    ! states, advances to the requested output times, and writes the fields
    ! and a key=value log.  Solver and case logic live in their own modules.
    use precision_select
    use scheme_params, only: scheme
    use e2d_control
    use e2d_globals, only: e2d_allocate, e2d_deallocate
    use e2d_cases, only: e2d_initial, e2d_boundary
    use e2d_solver, only: e2d_advance
    use e2d_io, only: e2d_write_field, e2d_write_log, e2d_count_nonfinite
    implicit none

    integer :: ierr
    integer(integer_precision) :: steps, interp_corr, cell_corr, nan_count
    integer :: seg_steps, seg_interp, seg_cell
    character(len=256) :: msg
    character(len=512) :: base, out_path(2), log_path
    character(len=32) :: op
    character(len=16) :: nz_str, nx_str
    real(precision) :: cpu_seconds, seg_cpu

    call e2d_control_init('input.nml', ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e2d_resolve_operator(trim(scheme), op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e2d_ensure_output_dirs(trim(output_dir))

    call e2d_allocate(nx, nz)

    call e2d_initial(case_2d, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e2d_boundary(case_2d, 0.0_precision)

    write (nz_str, '(i0)') nz
    write (nx_str, '(i0)') nx
    base = trim(case_2d) // '_' // trim(nz_str) // 'x' // trim(nx_str)
    out_path(1) = trim(output_dir) // '/' // trim(base) // '_t1.dat'
    out_path(2) = trim(output_dir) // '/' // trim(base) // '_t2.dat'
    log_path = trim(output_dir) // '/' // trim(case_2d) // '.log'

    if (write_sample) then
        call e2d_write_field(trim(output_dir) // '/' // trim(base) // '_t0.dat', &
                             x_origin, z_origin, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if
    end if

    steps = 0
    interp_corr = 0
    cell_corr = 0
    cpu_seconds = 0.0_precision

    if (tend2 > 0.0_precision) then
        call e2d_advance(case_2d, op, cfl, 0.0_precision, tend2, seg_steps, &
                         seg_interp, seg_cell, seg_cpu)
        steps = steps + seg_steps
        interp_corr = interp_corr + seg_interp
        cell_corr = cell_corr + seg_cell
        cpu_seconds = cpu_seconds + seg_cpu

        call e2d_write_field(out_path(1), x_origin, z_origin, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if

        call e2d_advance(case_2d, op, cfl, tend2, tend, seg_steps, &
                         seg_interp, seg_cell, seg_cpu)
        steps = steps + seg_steps
        interp_corr = interp_corr + seg_interp
        cell_corr = cell_corr + seg_cell
        cpu_seconds = cpu_seconds + seg_cpu

        call e2d_write_field(out_path(2), x_origin, z_origin, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if
    else
        call e2d_advance(case_2d, op, cfl, 0.0_precision, tend, seg_steps, &
                         seg_interp, seg_cell, seg_cpu)
        steps = steps + seg_steps
        interp_corr = interp_corr + seg_interp
        cell_corr = cell_corr + seg_cell
        cpu_seconds = cpu_seconds + seg_cpu

        call e2d_write_field(out_path(1), x_origin, z_origin, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if
    end if

    nan_count = e2d_count_nonfinite()

    call e2d_write_log(log_path, case_2d, scheme, op, mesh_num, cfl, tend, &
                       gamma, steps, interp_corr, cell_corr, nan_count, &
                       cpu_seconds, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call e2d_deallocate()

    if (nan_count > 0) then
        write (*, '(a,i0)') 'error: non-finite cells: ', nan_count
        error stop 1
    end if

end program euler2d
