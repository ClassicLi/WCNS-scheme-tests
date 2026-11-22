program euler1d_exact
    ! Driver for the Sod/Lax exact Riemann solution.
    use precision_select
    use e1d_control, only: e1d_exact_control_init, e1d_ensure_output_dirs, &
                           case_1d, mesh_num, tend, output_dir
    use e1d_exact, only: e1d_exact_profile
    use e1d_io, only: e1d_write_exact
    implicit none

    integer :: ierr, n
    character(len=256) :: msg, path
    real(precision), allocatable :: rho(:), u(:), p(:)

    call e1d_exact_control_init('input.nml', ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    n = int(mesh_num)
    allocate (rho(0:n), u(0:n), p(0:n))
    call e1d_exact_profile(n, case_1d, tend, rho, u, p)

    call e1d_ensure_output_dirs(output_dir)
    path = trim(output_dir) // '/csv/exact_' // trim(case_1d) // '.csv'
    call e1d_write_exact(path, n, rho, u, p, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

end program euler1d_exact
