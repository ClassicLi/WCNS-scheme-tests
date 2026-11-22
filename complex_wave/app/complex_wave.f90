program complex_wave
    ! Driver for the complex-wave test.
    !
    ! Reads the three namelist groups, builds the reference four-feature
    ! initial profile, advances it for the requested number of wave periods
    ! with RK4, and writes one snapshot per period (including t = 0).
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use cw_control
    use cw_params
    use scheme_params, only: scheme, read_scheme_namelist, validate_scheme_params
    use scheme_dispatch, only: scheme_output_name
    use cw_periodic, only: boundary_cycle
    use cw_initial, only: initial_complex_wave
    use cw_time, only: rk4_step
    use cw_wcns_rhs, only: wcns_rhs
    use cw_io, only: ensure_output_dirs, write_csv_snapshot, &
                     write_tecplot_snapshots, write_run_log
    implicit none

    ! The interpolation template reaches six points on either side.
    integer(integer_precision), parameter :: halo = 7_integer_precision

#ifdef PRECISION_DP
    character(len=*), parameter :: precision_label = 'dp'
#else
    character(len=*), parameter :: precision_label = 'qp'
#endif

    integer :: ierr, u_unit
    character(len=256) :: msg
    character(len=64) :: run_name, tag
    character(len=512) :: csv_path, dat_path, log_path

    integer(integer_precision) :: nx, k, steps, nan_count
    real(precision), allocatable :: u(:), u1(:), snaps(:, :), times(:)
    real(precision) :: time, target, dt, period, umax_abs, h
    real(precision) :: dt_min, dt_max, umin, umax
    logical :: last

    ! ---------------------------------------------------------------
    ! Configuration
    ! ---------------------------------------------------------------

    ! Read &complex_wave_nml (mesh size, periods, CFL, output directory).
    call cw_control_init('input.nml', ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    ! Read &scheme_nml and &difference_nml from the same file.  The two
    ! readers rewind independently, so the order does not matter.
    open(newunit=u_unit, file='input.nml', status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot open input.nml'
        error stop 1
    end if

    call read_scheme_namelist(u_unit, ierr)
    if (ierr /= 0) then
        close(u_unit)
        write(*, '(a)') 'error: cannot read &scheme_nml'
        error stop 1
    end if

    call read_difference_namelist(u_unit, ierr)
    close(u_unit)
    if (ierr /= 0) then
        write(*, '(a)') 'error: cannot read &difference_nml'
        error stop 1
    end if

    call validate_scheme_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call validate_difference_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    ! ---------------------------------------------------------------
    ! Work arrays
    ! ---------------------------------------------------------------

    nx = mesh_num
    h = 2.0_precision/real(nx, precision)

    allocate(u(-halo:nx+halo))
    allocate(u1(0:nx))
    allocate(snaps(0:nx, 0:periods))
    allocate(times(0:periods))

    call ensure_output_dirs(output_dir)

    ! ---------------------------------------------------------------
    ! Initial condition and snapshot t = 0
    ! ---------------------------------------------------------------

    call initial_complex_wave(nx, u(0:nx))
    call boundary_cycle(nx, u, halo)

    snaps(:, 0) = u(0:nx)
    times(0) = 0.0_precision

    umin = huge(1.0_precision)
    umax = -huge(1.0_precision)
    nan_count = 0
    call update_stats(nx, u(0:nx), umin, umax, nan_count)

    ! ---------------------------------------------------------------
    ! Advance one period at a time
    ! ---------------------------------------------------------------

    time = 0.0_precision
    steps = 0
    dt_min = huge(1.0_precision)
    dt_max = 0.0_precision
    period = cw_period()

    do k = 1, periods
        target = real(k, precision)*period

        do
            ! Reference CFL step based on the current maximum amplitude.
            umax_abs = maxval(abs(u(0:nx)))
            if (.not. ieee_is_finite(umax_abs) .or. umax_abs <= 0.0_precision) then
                write(*, '(a)') 'error: non-finite or zero solution amplitude'
                error stop 1
            end if

            dt = cfl*h/umax_abs
            if (.not. ieee_is_finite(dt) .or. dt <= 0.0_precision) then
                write(*, '(a)') 'error: nonpositive time step'
                error stop 1
            end if

            ! Clip the last step so a snapshot lands exactly on k*T.
            if (time + dt >= target) then
                dt = target - time
                last = .true.
            else
                last = .false.
            end if

            call rk4_step(nx, halo, dt, u(0:nx), wcns_rhs, u1)
            u(0:nx) = u1
            time = time + dt
            steps = steps + 1
            dt_min = min(dt_min, dt)
            dt_max = max(dt_max, dt)
            call update_stats(nx, u(0:nx), umin, umax, nan_count)

            ! Abort immediately if this step produced NaN/Inf.  Without this
            ! check a non-finite value appearing on the final clipped step
            ! would still be written out with a zero exit status.
            if (nan_count > 0) then
                write(*, '(a)') 'error: non-finite solution; aborting before output'
                error stop 1
            end if

            if (last) exit
        end do

        ! Remove accumulated round-off in the reported time, then snapshot.
        time = target
        snaps(:, k) = u(0:nx)
        times(k) = time
    end do

    ! ---------------------------------------------------------------
    ! Output
    ! ---------------------------------------------------------------

    run_name = trim(scheme_output_name(scheme)) // '_' // trim(difference_operator)

    do k = 0, periods
        write(tag, '(i0)') k
        csv_path = trim(output_dir) // '/csv/' // trim(run_name) // '_t' // &
                   trim(tag) // '.csv'
        call write_csv_snapshot(csv_path, nx, snaps(:, k))
    end do

    dat_path = trim(output_dir) // '/dat/' // trim(run_name) // '.dat'
    call write_tecplot_snapshots(dat_path, trim(run_name), nx, periods, times, snaps)

    log_path = trim(output_dir) // '/logs/' // trim(run_name) // '.log'
    call write_run_log(log_path, precision_label, trim(scheme), &
                       trim(difference_operator), nx, periods, cfl, dt_min, &
                       dt_max, steps, umin, umax, nan_count)

    write(*, '(a)') 'complex_wave: wrote ' // trim(run_name)

contains

    subroutine update_stats(nx, values, umin, umax, nan_count)
        ! Update extrema and count non-finite entries over the interior.
        integer(integer_precision), intent(in) :: nx
        real(precision), intent(in) :: values(0:nx)
        real(precision), intent(inout) :: umin, umax
        integer(integer_precision), intent(inout) :: nan_count

        integer(integer_precision) :: j

        do j = 0, nx
            if (.not. ieee_is_finite(values(j))) then
                nan_count = nan_count + 1
            else
                umin = min(umin, values(j))
                umax = max(umax, values(j))
            end if
        end do
    end subroutine update_stats

end program complex_wave
