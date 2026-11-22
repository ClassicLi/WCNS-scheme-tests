program test_moc322
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_framework, only: MOC322
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
    call MOC322(u, u_half)
    call check_close(u_half, exact, test_tol, 'MOC322 linear exact')

    ! Flat zero data must stay finite and zero.
    u = 0.0_precision
    call MOC322(u, u_half)
    call check_true(ieee_is_finite(u_half), 'MOC322 flat finite')
    call check_close(u_half, 0.0_precision, 0.0_precision, 'MOC322 flat zero')

    ! Third-order convergence on a smooth function.
    call check_true(moc322_order(1, 40, 80) >= 2.7_precision, 'MOC322 cos order 40->80 >= 2.7')
    call check_true(moc322_order(1, 80, 160) >= 2.7_precision, 'MOC322 cos order 80->160 >= 2.7')
    call check_true(moc322_order(1, 160, 320) >= 2.7_precision, 'MOC322 cos order 160->320 >= 2.7')

    ! MOC322 should retain third order at a first-order critical point.
    call check_true(moc322_order(2, 40, 80) >= 2.6_precision, 'MOC322 critical order 40->80 >= 2.6')
    call check_true(moc322_order(2, 80, 160) >= 2.6_precision, 'MOC322 critical order 80->160 >= 2.6')
    call check_true(moc322_order(2, 160, 320) >= 2.6_precision, 'MOC322 critical order 160->320 >= 2.6')

    ! Diagnostics.
    do l = -3, 4
        u(l) = cos(0.3_precision + real(l, precision) * h)
    end do
    call MOC322(u, u_half, diag)
    call check_true(diag%beta_L >= 0.0_precision, 'MOC322 beta_L nonnegative')
    call check_true(all(diag%beta(1:2) >= 0.0_precision), 'MOC322 sub betas nonnegative')
    call check_true(diag%omega_L >= 0.0_precision .and. diag%omega_L <= 1.0_precision, 'MOC322 omega_L in [0,1]')
    call check_true(all(diag%omega(1:2) >= 0.0_precision) .and. all(diag%omega(1:2) <= 1.0_precision), &
        'MOC322 sub omegas in [0,1]')
    call check_close(diag%omega_L + sum(diag%omega(1:2)), 1.0_precision, test_tol, 'MOC322 weights sum')

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
            f = exp(0.75_precision * (x - 1.0_precision)) * x * x
        end select
    end function test_value

    function moc322_l1_error(kind, n) result(err)
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
            call MOC322(uu, uh)
            ex = test_value(kind, xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function moc322_l1_error

    function moc322_order(kind, n1, n2) result(p)
        integer, intent(in) :: kind, n1, n2
        real(precision) :: p

        p = log(moc322_l1_error(kind, n1) / moc322_l1_error(kind, n2)) / log(2.0_precision)
    end function moc322_order

end program test_moc322
