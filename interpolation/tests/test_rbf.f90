program test_rbf
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use interp_rbf
    use interp_linear, only: linear_3rd_minus1_to_plus1
    use scheme_params, only: rbf_floor_value
    use scheme_diagnostics_mod
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
#ifdef PRECISION_DP
    ! GA at 160->320 is roundoff-limited in double precision.
    real(precision), parameter :: rbf5_order_high = 3.5_precision
#else
    real(precision), parameter :: rbf5_order_high = 4.5_precision
#endif
    real(precision) :: u(-3:4), u_half, u_lin, exact
    type(scheme_diagnostics) :: diag
    integer :: kind, l

    ! Linear data must be reproduced exactly by every variant.
    do kind = 1, 6
        do l = -3, 4
            u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
        end do
        exact = 0.3_precision - 1.2_precision * 0.5_precision * h
        call rbf_apply(kind, u, u_half)
        call check_close(u_half, exact, test_tol, kind_label('RBF linear exact', kind))
    end do

    ! Flat zero data must stay finite and zero for every variant.
    u = 0.0_precision
    do kind = 1, 6
        call rbf_apply(kind, u, u_half)
        call check_true(ieee_is_finite(u_half), kind_label('RBF flat finite', kind))
        call check_close(u_half, 0.0_precision, 0.0_precision, kind_label('RBF flat zero', kind))
    end do

    ! Three-point variants: third-order convergence.
    do kind = 1, 3
        call check_true(rbf_order(kind, 40, 80) >= 2.7_precision, kind_label('3-point order 40->80', kind))
        call check_true(rbf_order(kind, 80, 160) >= 2.7_precision, kind_label('3-point order 80->160', kind))
        call check_true(rbf_order(kind, 160, 320) >= 2.7_precision, kind_label('3-point order 160->320', kind))
    end do

    ! Five-point variants: fifth-order convergence.
    do kind = 4, 6
        call check_true(rbf_order(kind, 40, 80) >= 4.5_precision, kind_label('5-point order 40->80', kind))
        call check_true(rbf_order(kind, 80, 160) >= 4.5_precision, kind_label('5-point order 80->160', kind))
        call check_true(rbf_order(kind, 160, 320) >= rbf5_order_high, kind_label('5-point order 160->320', kind))
    end do

    ! Diagnostics: lambda_h2 is finite.
    do l = -3, 4
        u(l) = cos(0.3_precision + real(l, precision) * h)
    end do
    call RBF_GA_Perturbed_3(u, u_half, diag)
    call check_true(ieee_is_finite(diag%lambda_h2), 'RBF lambda_h2 finite')

    ! Floor semantics: a zero floor disables the perturbation, a large floor
    ! keeps the perturbation active.
    rbf_floor_value = 0.0_precision
    call RBF_GA_Perturbed_3(u, u_half)
    call linear_3rd_minus1_to_plus1(u, u_lin)
    call check_close(u_half, u_lin, test_tol, 'zero floor -> linear interpolation')

    rbf_floor_value = 1.0e6_precision
    call RBF_GA_Perturbed_3(u, u_half)
    call check_true(abs(u_half - u_lin) > 1.0e-12_precision, 'large floor keeps perturbation active')
    rbf_floor_value = 30.0_precision

    ! Extreme scaling must remain finite.
    do l = -3, 4
        u(l) = 1.0e30_precision * cos(0.3_precision + real(l, precision) * h)
    end do
    call RBF_IQ_Perturbed_5(u, u_half)
    call check_true(ieee_is_finite(u_half), 'RBF 1e30 scaling finite')

    do l = -3, 4
        u(l) = 1.0e-30_precision * cos(0.3_precision + real(l, precision) * h)
    end do
    call RBF_IQ_Perturbed_5(u, u_half)
    call check_true(ieee_is_finite(u_half), 'RBF 1e-30 scaling finite')

    call finish()

contains

    subroutine rbf_apply(kind, uu, uh)
        integer, intent(in) :: kind
        real(precision), intent(in) :: uu(-3:4)
        real(precision), intent(out) :: uh

        select case (kind)
        case (1)
            call RBF_GA_Perturbed_3(uu, uh)
        case (2)
            call RBF_MQ_Perturbed_3(uu, uh)
        case (3)
            call RBF_IQ_Perturbed_3(uu, uh)
        case (4)
            call RBF_GA_Perturbed_5(uu, uh)
        case (5)
            call RBF_MQ_Perturbed_5(uu, uh)
        case (6)
            call RBF_IQ_Perturbed_5(uu, uh)
        end select
    end subroutine rbf_apply

    function rbf_l1_error(kind, n) result(err)
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
            call rbf_apply(kind, uu, uh)
            ex = cos(xj + 0.5_precision * hh)
            err = err + abs(uh - ex)
        end do
        err = err / real(n, precision)
    end function rbf_l1_error

    function rbf_order(kind, n1, n2) result(p)
        integer, intent(in) :: kind, n1, n2
        real(precision) :: p

        p = log(rbf_l1_error(kind, n1) / rbf_l1_error(kind, n2)) / log(2.0_precision)
    end function rbf_order

    function kind_label(prefix, k) result(s)
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: k
        character(len=96) :: s

        write(s, '(a,i0)') trim(prefix) // ' kind=', k
    end function kind_label

end program test_rbf
