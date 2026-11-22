program adr_wcns
    ! ADR driver for WCNS interpolation with a selectable difference operator.
    use precision_select
    use adr_control
    use adr_dft, only: dft_mode, calc_phi
    use adr_time, only: rk3_step
    use adr_io, only: write_csv, write_log, write_tecplot
    use scheme_params
    use scheme_dispatch, only: scheme_output_name
    use wcns_params
    use wcns_dispatch, only: wcns_rhs
    implicit none

    integer :: ierr, u
    character(len=256) :: msg
    character(len=128) :: out_name
    character(len=512) :: csv_path, log_path, sample_path
    integer(integer_precision) :: n, j, n_mode
    real(precision), allocatable :: u0(:), u1(:), varphi(:)
    complex(precision), allocatable :: phi(:)
    complex(precision) :: hat0, hatt

    ! Read and validate &adr_nml.
    call adr_control_init('input.nml', ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    ! Read &scheme_nml and &difference_nml from the same input file.
    open(newunit=u, file='input.nml', status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot open input.nml'
        error stop 1
    end if
    call read_scheme_namelist(u, ierr)
    if (ierr /= 0) then
        close(u)
        write(*, '(a)') 'error: cannot read &scheme_nml'
        error stop 1
    end if
    call read_difference_namelist(u, ierr)
    close(u)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot read &difference_nml'
        error stop 1
    end if

    call validate_scheme_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call validate_wcns_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    out_name = trim(scheme_output_name(scheme)) // '_' // trim(difference_operator)

    ! Allocate per-mode output arrays and the real-space work arrays.
    n_mode = mesh_num / 2
    allocate(varphi(0:n_mode), phi(0:n_mode))
    allocate(u0(-7:mesh_num+7), u1(-7:mesh_num+7))

    call execute_command_line('mkdir -p "' // trim(output_dir) // '/csv" "' // &
        trim(output_dir) // '/logs" "' // trim(output_dir) // '/samples"')

    ! Sweep every resolved Fourier mode.
    do n = 0, n_mode
        do j = 0, mesh_num
            u0(j) = cos(real(n, precision) * real(j, precision) * 2.0_precision * PI &
                        / real(mesh_num, precision))
        end do

        call rk3_step(mesh_num, 7_integer_precision, dt, u0(0:mesh_num), wcns_rhs, u1(0:mesh_num))

        call dft_mode(u0(0:mesh_num), n, mesh_num, hat0)
        call dft_mode(u1(0:mesh_num), n, mesh_num, hatt)
        call calc_phi(hatt, hat0, dt, mesh_num, phi(n))
        varphi(n) = 2.0_precision * PI * real(n, precision) / real(mesh_num, precision)
    end do

    ! Write the CSV curve and the run metadata.
    csv_path = trim(output_dir) // '/csv/' // trim(out_name) // '.csv'
    log_path = trim(output_dir) // '/logs/' // trim(out_name) // '.log'
    call write_csv(csv_path, varphi, phi, n_mode)
    call write_log(log_path, out_name, difference_operator, mesh_num, dt, varphi, phi, n_mode)

    if (write_tecplot_sample) then
        sample_path = trim(output_dir) // '/samples/' // trim(out_name) // '.dat'
        call write_tecplot(sample_path, out_name, varphi, phi, n_mode)
    end if

end program adr_wcns
