module e2d_io
    ! Tecplot density output, key=value run logs, and non-finite counting
    ! for the 2D Euler test cases.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_select
    use e2d_globals, only: n, m, h, beta, ubar
    implicit none

contains

    subroutine e2d_write_field(path, x_origin, z_origin, ierr, msg)
        ! Write x, z, rho on the (n+1) x (m+1) grid as a point-packed zone.
        !
        ! The header matches the Tecplot ASCII layout used in the thesis:
        ! TITLE, VARIABLES, then a single ZONE with explicit I and J extents.
        character(len=*), intent(in) :: path
        real(precision), intent(in) :: x_origin, z_origin
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios
        integer(integer_precision) :: j, k

        ierr = 0
        msg = ''

        open (newunit=unit, file=trim(path), status='replace', action='write', &
              iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot write field: ' // trim(path)
            return
        end if

        write (unit, '(a)') 'TITLE="2D Euler"'
        write (unit, '(a)') 'VARIABLES="x","z","rho"'
        write (unit, '(a,i0,a,i0,a)') 'ZONE I=', n + 1, ', J=', m + 1, &
                                      ', DATAPACKING=POINT'
        do k = 0, m
            do j = 0, n
                write (unit, '(es24.16,1x,es24.16,1x,es24.16)') &
                    x_origin + real(j, precision)*h, &
                    z_origin + real(k, precision)*h, ubar(j, k, 1)
            end do
        end do
        close (unit)
    end subroutine e2d_write_field

    function e2d_count_nonfinite() result(count)
        ! Number of interior cells with a non-finite or non-positive state.
        !
        ! A cell is counted when any conserved component used here is not
        ! finite, when rho <= 0, or when the recovered pressure is not finite
        ! or not positive.
        integer :: count

        real(precision) :: rho, u, w, ee, p
        integer(integer_precision) :: j, k
        logical :: finite_state

        count = 0
        do k = 0, m
            do j = 0, n
                rho = ubar(j, k, 1)
                u = ubar(j, k, 2)
                w = ubar(j, k, 3)
                ee = ubar(j, k, 4)
                ! Test finiteness before any comparison or arithmetic:
                ! comparing a NaN with <= raises the invalid exception under
                ! -ffpe-trap=invalid, which would abort the counter instead
                ! of classifying the cell.
                finite_state = ieee_is_finite(rho) .and. &
                               ieee_is_finite(u) .and. &
                               ieee_is_finite(w) .and. &
                               ieee_is_finite(ee)
                if (.not. finite_state) then
                    count = count + 1
                else if (rho <= 0.0_precision) then
                    count = count + 1
                else
                    p = beta*(ee - 0.5_precision*rho* &
                              ((u/rho)**2 + (w/rho)**2))
                    if (.not. ieee_is_finite(p)) then
                        count = count + 1
                    else if (p <= 0.0_precision) then
                        count = count + 1
                    end if
                end if
            end do
        end do
    end function e2d_count_nonfinite

    subroutine e2d_write_log(path, case_name, scheme_name, op, mesh_num, cfl, &
                             tend, gama, steps, interp_corr, cell_corr, &
                             nan_count, cpu_seconds, ierr, msg)
        ! Write run metadata as key=value lines.
        character(len=*), intent(in) :: path, case_name, scheme_name, op
        integer(integer_precision), intent(in) :: mesh_num, steps, &
                                                  interp_corr, cell_corr, &
                                                  nan_count
        real(precision), intent(in) :: cfl, tend, gama, cpu_seconds
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios

        ierr = 0
        msg = ''

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
        write (unit, '(a,i0)') 'mesh=', mesh_num
        write (unit, '(a,es24.16)') 'cfl=', cfl
        write (unit, '(a,es24.16)') 'tend=', tend
        write (unit, '(a,es24.16)') 'gama=', gama
        write (unit, '(a,i0)') 'steps=', steps
        write (unit, '(a,i0)') 'interp_corrections=', interp_corr
        write (unit, '(a,i0)') 'cell_corrections=', cell_corr
        write (unit, '(a,i0)') 'nan_count=', nan_count
        write (unit, '(a,es24.16)') 'cpu_seconds=', cpu_seconds
        close (unit)
    end subroutine e2d_write_log

end module e2d_io
