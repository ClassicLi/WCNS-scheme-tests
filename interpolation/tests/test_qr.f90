program test_qr
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_qr
    use scheme_diagnostics_mod
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-3:4), u_half, exact
    type(scheme_diagnostics) :: diag
    integer :: l

    u = 2.0_precision
    call quasi_rational(u, u_half)
    call check_close(u_half, 2.0_precision, 1.0e-30_precision, 'QR constant data')

    do l = -3, 4
        u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
    end do
    call quasi_rational(u, u_half)
    exact = 0.3_precision - 1.2_precision * 0.5_precision * h
    call check_close(u_half, exact, test_tol, 'QR exact on linear data')

    u = 0.0_precision
    call quasi_rational(u, u_half)
    call check_true(ieee_is_finite(u_half), 'QR flat zero is finite')
    call check_close(u_half, 0.0_precision, 0.0_precision, 'QR flat zero')

    call quasi_rational(u, u_half, diag)
    call check_close(diag%beta(1), 0.0_precision, 0.0_precision, 'QR leaves beta zero')
    call check_close(diag%omega(1), 0.0_precision, 0.0_precision, 'QR leaves omega zero')

    call check_true(qr_order(40, 80) >= 2.7_precision, 'QR order 40->80 >= 2.7')
    call check_true(qr_order(80, 160) >= 2.7_precision, 'QR order 80->160 >= 2.7')
    call check_true(qr_order(160, 320) >= 2.7_precision, 'QR order 160->320 >= 2.7')

    call finish()

contains

    function qr_l1_error(n) result(err)
        integer, intent(in) :: n
        real(precision) :: err
        real(precision) :: hh, xj, ex
        real(precision) :: uu(-3:4), uh
        integer :: j, ll

        hh = 2.0_precision / real(n, precision)
        err = 0.0_precision
        do j = 0, n - 1
            xj = -1.0_precision + real(j, precision) * hh
            do ll = -3, 4
                uu(ll) = cos(xj + real(ll, precision) * hh)
            end do
            call quasi_rational(uu, uh)
            ex = cos(xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function qr_l1_error

    function qr_order(n1, n2) result(p)
        integer, intent(in) :: n1, n2
        real(precision) :: p

        p = log(qr_l1_error(n1) / qr_l1_error(n2)) / log(2.0_precision)
    end function qr_order

end program test_qr
