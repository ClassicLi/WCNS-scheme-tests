program test_linear
    use precision_select
    use interp_linear
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-3:4), u_half, u_ref
    real(precision) :: base(0:5), c(0:5)

    base = (/ 0.3_precision, -1.2_precision, 0.7_precision, 0.5_precision, -0.2_precision, 0.1_precision /)

    ! linear_2nd is exact through degree 1.
    c = base
    c(2:5) = 0.0_precision
    call fill_poly(u, c, h)
    call linear_2nd(u, u_half)
    call check_close(u_half, poly_eval(c, 0.5_precision * h), test_tol, &
        'linear_2nd exact on degree-1 polynomial')

    ! 3rd-order stencils are exact through degree 2.
    c = base
    c(3:5) = 0.0_precision
    call fill_poly(u, c, h)
    u_ref = poly_eval(c, 0.5_precision * h)
    call linear_3rd_minus1_to_plus1(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'linear_3rd_minus1_to_plus1 exact')
    call linear_3rd_minus2_to_0(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'linear_3rd_minus2_to_0 exact')
    call linear_3rd_0_to_plus2(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'linear_3rd_0_to_plus2 exact')
    call poly3(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'poly3 exact on degree-2 polynomial')

    ! 4th-order stencils are exact through degree 3.
    c = base
    c(4:5) = 0.0_precision
    call fill_poly(u, c, h)
    u_ref = poly_eval(c, 0.5_precision * h)
    call linear_4th_minus2_to_plus1(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'linear_4th_minus2_to_plus1 exact')
    call linear_4th_minus1_to_plus2(u, u_half)
    call check_close(u_half, u_ref, test_tol, 'linear_4th_minus1_to_plus2 exact')

    ! 5th-order stencil is exact through degree 4.
    c = base
    c(5) = 0.0_precision
    call fill_poly(u, c, h)
    call linear_5th_minus2_to_plus2(u, u_half)
    call check_close(u_half, poly_eval(c, 0.5_precision * h), test_tol, &
        'linear_5th_minus2_to_plus2 exact on degree-4 polynomial')

    ! 6th-order stencil is exact through degree 5.
    c = base
    call fill_poly(u, c, h)
    call linear_6th_minus2_to_plus3(u, u_half)
    call check_close(u_half, poly_eval(c, 0.5_precision * h), test_tol, &
        'linear_6th_minus2_to_plus3 exact on degree-5 polynomial')

    ! poly3 is an alias of the symmetric 3rd-order stencil.
    call linear_3rd_minus1_to_plus1(u, u_ref)
    call poly3(u, u_half)
    call check_close(u_half, u_ref, 0.0_precision, 'poly3 equals linear_3rd_minus1_to_plus1')

    call finish()

contains

    subroutine fill_poly(u, coef, h)
        real(precision), intent(out) :: u(-3:4)
        real(precision), intent(in) :: coef(0:5), h
        integer :: l

        do l = -3, 4
            u(l) = poly_eval(coef, real(l, precision) * h)
        end do
    end subroutine fill_poly

    function poly_eval(coef, x) result(y)
        real(precision), intent(in) :: coef(0:5), x
        real(precision) :: y
        integer :: k

        y = 0.0_precision
        do k = 5, 0, -1
            y = y * x + coef(k)
        end do
    end function poly_eval

end program test_linear
