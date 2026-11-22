module e1d_io
    ! CSV, log, and optional Tecplot output for the 1D Euler test cases.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_select
    use e1d_globals, only: halo
    use e1d_flux, only: e1d_pressure
    implicit none

contains

    subroutine e1d_write_profile(path, n, length, u, ierr, msg)
        ! Write x, rho, u, p on the cell centres.
        character(len=*), intent(in) :: path
        integer, intent(in) :: n
        real(precision), intent(in) :: length
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios, j
        real(precision) :: x, rho, vel, p

        ierr = 0
        msg = ''

        open (newunit=unit, file=trim(path), status='replace', action='write', &
              iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot write profile: ' // trim(path)
            return
        end if

        write (unit, '(a)') 'x,rho,u,p'
        do j = 0, n
            x = real(j, precision)*length/real(n, precision)
            rho = u(j, 1)
            vel = u(j, 2)/rho
            p = e1d_pressure(u(j, :))
            write (unit, '(es24.16,a,es24.16,a,es24.16,a,es24.16)') &
                x, ',', rho, ',', vel, ',', p
        end do
        close (unit)
    end subroutine e1d_write_profile

    subroutine e1d_write_exact(path, n, rho, u, p, ierr, msg)
        ! Write the exact solution on x_j = j/n.
        character(len=*), intent(in) :: path
        integer, intent(in) :: n
        real(precision), intent(in) :: rho(0:n), u(0:n), p(0:n)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios, j
        real(precision) :: x

        ierr = 0
        msg = ''

        open (newunit=unit, file=trim(path), status='replace', action='write', &
              iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot write exact profile: ' // trim(path)
            return
        end if

        write (unit, '(a)') 'x,rho,u,p'
        do j = 0, n
            x = real(j, precision)/real(n, precision)
            write (unit, '(es24.16,a,es24.16,a,es24.16,a,es24.16)') &
                x, ',', rho(j), ',', u(j), ',', p(j)
        end do
        close (unit)
    end subroutine e1d_write_exact

    function e1d_count_nonfinite(n, u) result(count)
        ! Number of cells with a non-finite or non-positive state.
        integer, intent(in) :: n
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        integer :: count

        real(precision) :: rho, p
        integer :: j

        count = 0
        do j = 0, n
            rho = u(j, 1)
            if (.not. ieee_is_finite(u(j, 1)) .or. &
                .not. ieee_is_finite(u(j, 2)) .or. &
                .not. ieee_is_finite(u(j, 3)) .or. rho <= 0.0_precision) then
                count = count + 1
            else
                p = e1d_pressure(u(j, :))
                if (.not. ieee_is_finite(p) .or. p <= 0.0_precision) count = count + 1
            end if
        end do
    end function e1d_count_nonfinite

    subroutine e1d_write_log(path, case_name, scheme_name, op, prec, n, cfl, tend, &
                             steps, ncorrect, ncorr_interp, cpu, u, ierr, msg)
        ! Write run metadata and state extrema as key=value lines.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: case_name, scheme_name, op, prec
        integer, intent(in) :: n, steps, ncorrect, ncorr_interp
        real(precision), intent(in) :: cfl, tend, cpu
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios, j
        real(precision) :: rho, p
        real(precision) :: min_rho, max_rho, min_p, max_p

        ierr = 0
        msg = ''

        min_rho = huge(1.0_precision)
        max_rho = -huge(1.0_precision)
        min_p = huge(1.0_precision)
        max_p = -huge(1.0_precision)
        do j = 0, n
            rho = u(j, 1)
            if (ieee_is_finite(rho) .and. rho > 0.0_precision) then
                p = e1d_pressure(u(j, :))
            else
                p = huge(1.0_precision)
            end if
            if (ieee_is_finite(rho)) then
                min_rho = min(min_rho, rho)
                max_rho = max(max_rho, rho)
            end if
            if (ieee_is_finite(p)) then
                min_p = min(min_p, p)
                max_p = max(max_p, p)
            end if
        end do

        open (newunit=unit, file=trim(path), status='replace', action='write', &
              iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot write log: ' // trim(path)
            return
        end if

        write (unit, '(a)') 'case=' // trim(case_name)
        write (unit, '(a)') 'scheme=' // trim(scheme_name)
        write (unit, '(a)') 'operator=' // trim(op)
        write (unit, '(a)') 'precision=' // trim(prec)
        write (unit, '(a,i0)') 'mesh_num=', n
        write (unit, '(a,es24.16)') 'cfl=', cfl
        write (unit, '(a,es24.16)') 'tend=', tend
        write (unit, '(a,i0)') 'steps=', steps
        write (unit, '(a,i0)') 'corrections=', ncorrect
        write (unit, '(a,i0)') 'interp_corrections=', ncorr_interp
        write (unit, '(a,es24.16)') 'cpu_seconds=', cpu
        write (unit, '(a,es24.16)') 'min_rho=', min_rho
        write (unit, '(a,es24.16)') 'max_rho=', max_rho
        write (unit, '(a,es24.16)') 'min_p=', min_p
        write (unit, '(a,es24.16)') 'max_p=', max_p
        write (unit, '(a,i0)') 'nan_count=', e1d_count_nonfinite(n, u)
        close (unit)
    end subroutine e1d_write_log

    subroutine e1d_write_sample(path, n, length, u, ierr, msg)
        ! Write a single-zone Tecplot sample.
        character(len=*), intent(in) :: path
        integer, intent(in) :: n
        real(precision), intent(in) :: length
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios, j
        real(precision) :: x, rho, vel, p

        ierr = 0
        msg = ''

        open (newunit=unit, file=trim(path), status='replace', action='write', &
              iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot write sample: ' // trim(path)
            return
        end if

        write (unit, '(a)') 'title="euler1d"'
        write (unit, '(a)') 'variables="x","rho","u","p"'
        write (unit, '(a,i0,a)') 'zone i=', n + 1, ',f=point'
        do j = 0, n
            x = real(j, precision)*length/real(n, precision)
            rho = u(j, 1)
            vel = u(j, 2)/rho
            p = e1d_pressure(u(j, :))
            write (unit, '(4(es24.16,1x))') x, rho, vel, p
        end do
        close (unit)
    end subroutine e1d_write_sample

end module e1d_io
