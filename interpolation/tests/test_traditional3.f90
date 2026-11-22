program test_traditional3
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_traditional
    use scheme_diagnostics_mod
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-3:4), u_half, exact
    type(scheme_diagnostics) :: diag
    integer :: kind, l

    ! Linear data must be reproduced exactly.
    do kind = 1, 3
        do l = -3, 4
            u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
        end do
        exact = 0.3_precision - 1.2_precision * 0.5_precision * h
        call trad3_apply(kind, u, u_half)
        call check_close(u_half, exact, test_tol, kind_label('linear exact', kind))
    end do

    ! Flat zero data must stay finite and zero.
    u = 0.0_precision
    do kind = 1, 3
        call trad3_apply(kind, u, u_half)
        call check_true(ieee_is_finite(u_half), kind_label('flat finite', kind))
        call check_close(u_half, 0.0_precision, 0.0_precision, kind_label('flat zero', kind))
    end do

    ! MUSCL is second order; JS3 and Z3 are third order.
    call check_true(trad3_order(1, 40, 80) >= 1.7_precision, 'MUSCL order 40->80 >= 1.7')
    call check_true(trad3_order(1, 80, 160) >= 1.7_precision, 'MUSCL order 80->160 >= 1.7')
    call check_true(trad3_order(1, 160, 320) >= 1.7_precision, 'MUSCL order 160->320 >= 1.7')

    do kind = 2, 3
        call check_true(trad3_order(kind, 40, 80) >= 2.7_precision, kind_label('order 40->80 >= 2.7', kind))
        call check_true(trad3_order(kind, 80, 160) >= 2.7_precision, kind_label('order 80->160 >= 2.7', kind))
        call check_true(trad3_order(kind, 160, 320) >= 2.7_precision, kind_label('order 160->320 >= 2.7', kind))
    end do

    ! Z3 diagnostics.
    do l = -3, 4
        u(l) = cos(0.3_precision + real(l, precision) * h)
    end do
    call Z3_interp(u, u_half, diag)
    call check_true(diag%beta(1) >= 0.0_precision .and. diag%beta(2) >= 0.0_precision, 'Z3 beta nonnegative')
    call check_close(diag%omega(1) + diag%omega(2), 1.0_precision, test_tol, 'Z3 omega sum')
    call check_true(diag%tau >= 0.0_precision, 'Z3 tau nonnegative')

    ! MUSCL monotonicity on a step.
    u = 0.0_precision
    u(0) = 1.0_precision
    u(1) = 1.0_precision
    call muscl_2nd(u, u_half)
    call check_true(u_half >= min(u(0), u(1)) .and. u_half <= max(u(0), u(1)), &
        'MUSCL step within local bounds')

    u = 0.0_precision
    u(-1) = 1.0_precision
    u(0) = 0.0_precision
    u(1) = 1.0_precision
    call muscl_2nd(u, u_half)
    call check_true(u_half >= min(u(0), u(1)) .and. u_half <= max(u(0), u(1)), &
        'MUSCL dip within local bounds')

    call finish()

contains

    subroutine trad3_apply(kind, uu, uh)
        integer, intent(in) :: kind
        real(precision), intent(in) :: uu(-3:4)
        real(precision), intent(out) :: uh

        select case (kind)
        case (1)
            call muscl_2nd(uu, uh)
        case (2)
            call JS3_interp(uu, uh)
        case (3)
            call Z3_interp(uu, uh)
        end select
    end subroutine trad3_apply

    function trad3_l1_error(kind, n) result(err)
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
                uu(ll) = cos(xj + real(ll, precision) * hh)
            end do
            call trad3_apply(kind, uu, uh)
            ex = cos(xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function trad3_l1_error

    function trad3_order(kind, n1, n2) result(p)
        integer, intent(in) :: kind, n1, n2
        real(precision) :: p

        p = log(trad3_l1_error(kind, n1) / trad3_l1_error(kind, n2)) / log(2.0_precision)
    end function trad3_order

    function kind_label(prefix, k) result(s)
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: k
        character(len=96) :: s

        write(s, '(a,i0)') trim(prefix) // ' kind=', k
    end function kind_label

end program test_traditional3
