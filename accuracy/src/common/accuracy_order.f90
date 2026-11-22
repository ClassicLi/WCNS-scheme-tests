module accuracy_order
    ! Error norms, convergence orders, and CSV output shared by both drivers.
    use precision_select
    implicit none

    ! Maximum number of refinement levels stored in one run.
    integer(integer_precision), parameter :: max_levels = 20_integer_precision

    integer(integer_precision) :: stored_levels = 0_integer_precision
    integer(integer_precision) :: stored_n(max_levels) = 0_integer_precision
    real(precision) :: stored_l1(max_levels) = 0.0_precision
    real(precision) :: stored_l2(max_levels) = 0.0_precision
    real(precision) :: stored_linf(max_levels) = 0.0_precision
    real(precision) :: stored_cpu(max_levels) = 0.0_precision

contains

    subroutine order_reset()
        ! Clear the stored table before a new run.
        stored_levels = 0_integer_precision
        stored_n = 0_integer_precision
        stored_l1 = 0.0_precision
        stored_l2 = 0.0_precision
        stored_linf = 0.0_precision
        stored_cpu = 0.0_precision
    end subroutine order_reset

    subroutine order_record(level, n, l1, l2, linf, cpu_time)
        ! Store one refinement level.  Levels are 1-based and consecutive.
        integer, intent(in) :: level
        integer(integer_precision), intent(in) :: n
        real(precision), intent(in) :: l1, l2, linf, cpu_time

        if (level < 1 .or. level > max_levels) then
            write (*, '(a,i0)') 'order_record: level out of range: ', level
            error stop 1
        end if

        stored_n(level) = n
        stored_l1(level) = l1
        stored_l2(level) = l2
        stored_linf(level) = linf
        stored_cpu(level) = cpu_time

        if (level > stored_levels) stored_levels = level
    end subroutine order_record

    function order_from_errors(e_prev, e_curr) result(p)
        ! Observed order p = log(e_prev/e_curr)/log(2).
        real(precision), intent(in) :: e_prev, e_curr
        real(precision) :: p

        if (e_prev <= 0.0_precision .or. e_curr <= 0.0_precision) then
            p = 0.0_precision
        else
            p = log(e_prev/e_curr)/log(2.0_precision)
        end if
    end function order_from_errors

    subroutine error_norms_1d(n, a, b, l1, l2, linf)
        ! L1, L2, and Linf differences between two 1D fields on 0:n.
        integer, intent(in) :: n
        real(precision), intent(in) :: a(0:n), b(0:n)
        real(precision), intent(out) :: l1, l2, linf

        real(precision) :: diff
        integer :: j

        l1 = 0.0_precision
        l2 = 0.0_precision
        linf = 0.0_precision

        do j = 0, n
            diff = abs(a(j) - b(j))
            l1 = l1 + diff
            l2 = l2 + diff*diff
            if (diff > linf) linf = diff
        end do

        l1 = l1/real(n + 1, precision)
        l2 = sqrt(l2/real(n + 1, precision))
    end subroutine error_norms_1d

    subroutine error_norms_2d(n, m, a, b, l1, l2, linf)
        ! L1, L2, and Linf differences between two 2D fields on 0:n, 0:m.
        integer, intent(in) :: n, m
        real(precision), intent(in) :: a(0:n, 0:m), b(0:n, 0:m)
        real(precision), intent(out) :: l1, l2, linf

        real(precision) :: diff
        integer :: j, k

        l1 = 0.0_precision
        l2 = 0.0_precision
        linf = 0.0_precision

        do k = 0, m
            do j = 0, n
                diff = abs(a(j, k) - b(j, k))
                l1 = l1 + diff
                l2 = l2 + diff*diff
                if (diff > linf) linf = diff
            end do
        end do

        l1 = l1/(real(n + 1, precision)*real(m + 1, precision))
        l2 = sqrt(l2/(real(n + 1, precision)*real(m + 1, precision)))
    end subroutine error_norms_2d

    subroutine order_write_csv(path, scheme_name, operator_name, ierr, msg)
        ! Write the stored table as CSV, adding the observed orders.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: scheme_name, operator_name
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: unit, ios, i

        ierr = 0
        msg = ''

        open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot open output file: ' // trim(path)
            return
        end if

        write (unit, '(a)') '# scheme=' // trim(scheme_name) // ' operator=' // &
            trim(operator_name)
        write (unit, '(a)') 'n,L1,order_L1,L2,order_L2,Linf,order_Linf,cpu_time'

        do i = 1, stored_levels
            if (i == 1) then
                write (unit, '(i0,7(a,es32.20e4))') stored_n(i), &
                    ',', stored_l1(i), ',', 0.0_precision, &
                    ',', stored_l2(i), ',', 0.0_precision, &
                    ',', stored_linf(i), ',', 0.0_precision, &
                    ',', stored_cpu(i)
            else
                write (unit, '(i0,7(a,es32.20e4))') stored_n(i), &
                    ',', stored_l1(i), ',', order_from_errors(stored_l1(i - 1), stored_l1(i)), &
                    ',', stored_l2(i), ',', order_from_errors(stored_l2(i - 1), stored_l2(i)), &
                    ',', stored_linf(i), ',', order_from_errors(stored_linf(i - 1), stored_linf(i)), &
                    ',', stored_cpu(i)
            end if
        end do

        close (unit)
    end subroutine order_write_csv

end module accuracy_order
