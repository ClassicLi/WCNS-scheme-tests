program test_moc33
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_framework, only: MOC33
    use scheme_params, only: moc33_linear
    use scheme_diagnostics_mod
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision), parameter :: pi = acos(-1.0_precision)
    real(precision) :: u(-3:4), u_half, exact
    type(scheme_diagnostics) :: diag
    integer :: l

    call check_choice('LINEAR')
    call check_choice('GA')

    call finish()

contains

    subroutine check_choice(choice)
        character(len=*), intent(in) :: choice

        moc33_linear = choice

        do l = -3, 4
            u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
        end do
        exact = 0.3_precision - 1.2_precision * 0.5_precision * h
        call MOC33(u, u_half)
        call check_close(u_half, exact, test_tol, trim(choice) // ' linear exact')

        call check_true(moc33_order(choice, 1, 40, 80) >= 2.7_precision, trim(choice) // ' cos order 40->80 >= 2.7')
        call check_true(moc33_order(choice, 1, 80, 160) >= 2.7_precision, trim(choice) // ' cos order 80->160 >= 2.7')
        call check_true(moc33_order(choice, 1, 160, 320) >= 2.7_precision, trim(choice) // ' cos order 160->320 >= 2.7')

        call check_true(moc33_order(choice, 2, 40, 80) >= 2.6_precision, trim(choice) // ' critical order 40->80 >= 2.6')
        call check_true(moc33_order(choice, 2, 80, 160) >= 2.6_precision, trim(choice) // ' critical order 80->160 >= 2.6')
        call check_true(moc33_order(choice, 2, 160, 320) >= 2.6_precision, trim(choice) // ' critical order 160->320 >= 2.6')

        u = 0.0_precision
        call MOC33(u, u_half)
        call check_true(ieee_is_finite(u_half), trim(choice) // ' flat finite')
        call check_close(u_half, 0.0_precision, 0.0_precision, trim(choice) // ' flat zero')

        do l = -3, 4
            u(l) = cos(0.3_precision + real(l, precision) * h)
        end do
        call MOC33(u, u_half, diag)
        call check_true(diag%beta_L >= 0.0_precision .and. diag%beta_N >= 0.0_precision, &
            trim(choice) // ' betas nonnegative')
        call check_close(diag%omega_L + diag%omega_N, 1.0_precision, test_tol, &
            trim(choice) // ' weights sum')

        call check_switch(choice)
    end subroutine check_choice

    subroutine check_switch(choice)
        character(len=*), intent(in) :: choice
        integer, parameter :: n = 80
        real(precision) :: hh, xj, best
        integer :: j, jbest, ll
        type(scheme_diagnostics) :: d

        hh = 2.0_precision / real(n, precision)
        best = huge(1.0_precision)
        jbest = 0
        do j = 0, n - 1
            xj = -1.0_precision + real(j, precision) * hh
            if (abs(xj - 0.375_precision) < best) then
                best = abs(xj - 0.375_precision)
                jbest = j
            end if
        end do

        xj = -1.0_precision + real(jbest, precision) * hh
        do ll = -3, 4
            u(ll) = udf4(xj + real(ll, precision) * hh)
        end do
        call MOC33(u, u_half, d)
        call check_true(d%omega_L < 0.1_precision, trim(choice) // ' switches to QR at udf_4 jump')
    end subroutine check_switch

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

    function udf4(x) result(f)
        real(precision), intent(in) :: x
        real(precision) :: f

        if (x <= 0.375_precision) then
            f = x * x * x * exp(0.75_precision * (x - 1.0_precision))
        else
            f = 1.0_precision + sin(pi * x)
        end if
    end function udf4

    function moc33_l1_error(choice, kind, n) result(err)
        character(len=*), intent(in) :: choice
        integer, intent(in) :: kind, n
        real(precision) :: err
        real(precision) :: hh, xj, ex
        real(precision) :: uu(-3:4), uh
        integer :: j, ll

        moc33_linear = choice
        hh = 2.0_precision / real(n, precision)
        err = 0.0_precision
        do j = 0, n - 1
            xj = -1.0_precision + real(j, precision) * hh
            do ll = -3, 4
                uu(ll) = test_value(kind, xj + real(ll, precision) * hh)
            end do
            call MOC33(uu, uh)
            ex = test_value(kind, xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function moc33_l1_error

    function moc33_order(choice, kind, n1, n2) result(p)
        character(len=*), intent(in) :: choice
        integer, intent(in) :: kind, n1, n2
        real(precision) :: p

        p = log(moc33_l1_error(choice, kind, n1) / moc33_l1_error(choice, kind, n2)) / log(2.0_precision)
    end function moc33_order

end program test_moc33
