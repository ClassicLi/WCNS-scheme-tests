program adr_weno
    ! ADR driver for the WENO reconstruction family.
    use precision_select
    use adr_control
    use adr_dft, only: dft_mode, calc_phi
    use adr_time, only: rk3_step
    use adr_io, only: write_csv, write_log, write_tecplot
    use weno_params
    use weno_dispatch, only: weno_rhs
    implicit none

    integer :: ierr, u
    character(len=256) :: msg
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

    ! Read &weno_nml from the same input file.
    open(newunit=u, file='input.nml', status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot open input.nml'
        error stop 1
    end if
    call read_weno_namelist(u, ierr)
    close(u)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot read &weno_nml'
        error stop 1
    end if

    call validate_weno_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    ! Allocate per-mode output arrays and the real-space work arrays.
    n_mode = mesh_num / 2
    allocate(varphi(0:n_mode), phi(0:n_mode))
    allocate(u0(-9:mesh_num+9), u1(-9:mesh_num+9))

    call execute_command_line('mkdir -p "' // trim(output_dir) // '/csv" "' // &
        trim(output_dir) // '/logs" "' // trim(output_dir) // '/samples"')

    ! Sweep every resolved Fourier mode.
    do n = 0, n_mode
        do j = 0, mesh_num
            u0(j) = cos(real(n, precision) * real(j, precision) * 2.0_precision * PI &
                        / real(mesh_num, precision))
        end do

        call rk3_step(mesh_num, 9_integer_precision, dt, u0(0:mesh_num), weno_rhs, u1(0:mesh_num))

        call dft_mode(u0(0:mesh_num), n, mesh_num, hat0)
        call dft_mode(u1(0:mesh_num), n, mesh_num, hatt)
        call calc_phi(hatt, hat0, dt, mesh_num, phi(n))
        varphi(n) = 2.0_precision * PI * real(n, precision) / real(mesh_num, precision)
    end do

    ! Write the CSV curve and the run metadata.
    csv_path = trim(output_dir) // '/csv/' // trim(scheme) // '.csv'
    log_path = trim(output_dir) // '/logs/' // trim(scheme) // '.log'
    call write_csv(csv_path, varphi, phi, n_mode)
    call write_log(log_path, scheme, '', mesh_num, dt, varphi, phi, n_mode)

    if (write_tecplot_sample) then
        sample_path = trim(output_dir) // '/samples/' // trim(scheme) // '.dat'
        call write_tecplot(sample_path, scheme, varphi, phi, n_mode)
    end if

end program adr_weno
