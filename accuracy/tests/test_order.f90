program test_order
    ! Unit tests for error norms, convergence order, and CSV output.
    use precision_select
    use accuracy_order
    use test_util
    implicit none

    real(precision) :: a(0:4), b(0:4)
    real(precision) :: e1, e2, einf
    integer :: i, ierr
    character(len=256) :: msg

    call check_close(order_from_errors(1.0_precision, 1.0_precision/8.0_precision), &
                     3.0_precision, test_tol, 'third order')
    call check_close(order_from_errors(1.0_precision, 1.0_precision/32.0_precision), &
                     5.0_precision, test_tol, 'fifth order')
    call check_close(order_from_errors(0.0_precision, 1.0_precision), &
                     0.0_precision, 0.0_precision, 'zero error gives zero order')

    do i = 0, 4
        a(i) = 1.0_precision
        b(i) = 1.0_precision
    end do
    b(0) = 3.0_precision

    call error_norms_1d(4, a, b, e1, e2, einf)
    call check_close(e1, 2.0_precision/5.0_precision, test_tol, 'L1')
    call check_close(e2, sqrt(4.0_precision/5.0_precision), test_tol, 'L2')
    call check_close(einf, 2.0_precision, test_tol, 'Linf')

    call order_reset()
    call order_record(1, 10_integer_precision, 1.0_precision, 1.0_precision, &
                      1.0_precision, 0.0_precision)
    call order_record(2, 20_integer_precision, 1.0_precision/8.0_precision, &
                      1.0_precision/8.0_precision, 1.0_precision/8.0_precision, &
                      0.0_precision)
    call order_write_csv('accuracy_order_test.csv', 'Z3', 'HCS4', ierr, msg)
    call check_true(ierr == 0, 'CSV written')

    call finish()
end program test_order
