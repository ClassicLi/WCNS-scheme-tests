program test_moz
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_framework, only: MOZ
    use scheme_diagnostics_mod
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-3:4), u_half, exact
    type(scheme_diagnostics) :: diag
    integer :: l

    ! Linear data must be reproduced exactly.
    do l = -3, 4
        u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
    end do
    exact = 0.3_precision - 1.2_precision * 0.5_precision * h
    call MOZ(u, u_half)
    call check_close(u_half, exact, test_tol, 'MOZ linear exact')

    ! Flat zero data must stay finite and zero.
    u = 0.0_precision
    call MOZ(u, u_half)
    call check_true(ieee_is_finite(u_half), 'MOZ flat finite')
    call check_close(u_half, 0.0_precision, 0.0_precision, 'MOZ flat zero')

    ! Fifth-order convergence.
    call check_true(moz_order(1, 40, 80) >= 4.5_precision, 'MOZ cos order 40->80 >= 4.5')
    call check_true(moz_order(1, 80, 160) >= 4.5_precision, 'MOZ cos order 80->160 >= 4.5')
    call check_true(moz_order(1, 160, 320) >= 4.5_precision, 'MOZ cos order 160->320 >= 4.5')

    ! Second-order critical point.
    call check_true(moz_order(2, 40, 80) >= 4.4_precision, 'MOZ critical order 40->80 >= 4.4')
    call check_true(moz_order(2, 80, 160) >= 4.4_precision, 'MOZ critical order 80->160 >= 4.4')
    call check_true(moz_order(2, 160, 320) >= 4.4_precision, 'MOZ critical order 160->320 >= 4.4')

    ! Diagnostics.
    do l = -3, 4
        u(l) = cos(0.3_precision + real(l, precision) * h)
    end do
    call MOZ(u, u_half, diag)
    call check_true(all(diag%beta(1:3) >= 0.0_precision), 'MOZ beta nonnegative')
    call check_close(sum(diag%omega(1:3)), 1.0_precision, test_tol, 'MOZ omega sum')
    call check_true(diag%tau >= 0.0_precision, 'MOZ tau nonnegative')
    call check_true(diag%h_noscale >= 0.0_precision, 'MOZ h_noscale nonnegative')

    ! Extreme scaling must remain finite.
    do l = -3, 4
        u(l) = 1.0e30_precision * cos(0.3_precision + real(l, precision) * h)
    end do
    call MOZ(u, u_half)
    call check_true(ieee_is_finite(u_half), 'MOZ 1e30 scaling finite')

    call finish()

contains

    function test_value(kind, x) result(f)
        integer, intent(in) :: kind
        real(precision), intent(in) :: x
        real(precision) :: f

        select case (kind)
        case (1)
            f = cos(x)
        case (2)
            f = exp(0.75_precision * (x - 1.0_precision)) * x * x * x
        end select
    end function test_value

    function moz_l1_error(kind, n) result(err)
        integer, intent(in) :: kind, n
        real(precision) :: err
        real(precision) :: hh, xj, ex
        real(precision) :: uu(-3:4), uh
        integer :: j, ll

        hh = 2.0_precision / real(n, precision)
        err = 0.0_precision
        do j = 0, n - 1
            xj = -1.0_precision + real(j, precision) * hh
            do ll = -3, 4
                uu(ll) = test_value(kind, xj + real(ll, precision) * hh)
            end do
            call MOZ(uu, uh)
            ex = test_value(kind, xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function moz_l1_error

    function moz_order(kind, n1, n2) result(p)
        integer, intent(in) :: kind, n1, n2
        real(precision) :: p

        p = log(moz_l1_error(kind, n1) / moz_l1_error(kind, n2)) / log(2.0_precision)
    end function moz_order

end program test_moz
