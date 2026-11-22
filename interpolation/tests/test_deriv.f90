program test_deriv
    use precision_select
    use interp_deriv
    use test_util
    implicit none

    real(precision) :: u(-3:4)
    real(precision) :: coef(0:6)

    coef = (/ 0.3_precision, -1.2_precision, 0.7_precision, 0.5_precision, &
              -0.2_precision, 0.1_precision, 0.05_precision /)

    call fill_poly(u, coef, 2)
    call check_close(get_derivative_1st_tr_2(u), deriv_exact(coef, 1), test_tol, &
        '1st derivative tr2 exact on degree 2')

    call fill_poly(u, coef, 1)
    call check_close(get_derivative_1st_tr_1_left(u), deriv_exact(coef, 1), test_tol, &
        '1st derivative tr1 left exact on degree 1')
    call check_close(get_derivative_1st_tr_1_right(u), deriv_exact(coef, 1), test_tol, &
        '1st derivative tr1 right exact on degree 1')

    call fill_poly(u, coef, 3)
    call check_close(get_derivative_3rd_tr_2(u), deriv_exact(coef, 3), test_tol, &
        '3rd derivative tr2 exact on degree 3')
    call check_close(get_derivative_3rd_tr_1_left(u), deriv_exact(coef, 3), test_tol, &
        '3rd derivative tr1 left exact on degree 3')
    call check_close(get_derivative_3rd_tr_1_right(u), deriv_exact(coef, 3), test_tol, &
        '3rd derivative tr1 right exact on degree 3')

    call fill_poly(u, coef, 6)
    call check_close(get_derivative_5th_tr_2(u), deriv_exact(coef, 5), test_tol, &
        '5th derivative tr2 exact on degree 6')

    call fill_poly(u, coef, 5)
    call check_close(get_derivative_5th_tr_1_left(u), deriv_exact(coef, 5), test_tol, &
        '5th derivative tr1 left exact on degree 5')
    call check_close(get_derivative_5th_tr_1_right(u), deriv_exact(coef, 5), test_tol, &
        '5th derivative tr1 right exact on degree 5')

    call finish()

contains

    subroutine fill_poly(u, coef, degree)
        real(precision), intent(out) :: u(-3:4)
        real(precision), intent(in) :: coef(0:6)
        integer, intent(in) :: degree
        real(precision) :: c(0:6)
        integer :: l

        c = coef
        if (degree < 6) c(degree + 1:6) = 0.0_precision
        do l = -3, 4
            u(l) = poly_eval(c, real(l, precision))
        end do
    end subroutine fill_poly

    function poly_eval(c, x) result(y)
        real(precision), intent(in) :: c(0:6), x
        real(precision) :: y
        integer :: k

        y = 0.0_precision
        do k = 6, 0, -1
            y = y * x + c(k)
        end do
    end function poly_eval

    function deriv_exact(c, m) result(d)
        real(precision), intent(in) :: c(0:6)
        integer, intent(in) :: m
        real(precision) :: d, fact
        integer :: k

        fact = 1.0_precision
        do k = 1, m
            fact = fact * real(k, precision)
        end do
        d = fact * c(m)
    end function deriv_exact

end program test_deriv
