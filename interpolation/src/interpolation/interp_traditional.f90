module interp_traditional
    use precision_select
    use interp_linear, only: linear_3rd_minus2_to_0, linear_3rd_minus1_to_plus1, linear_3rd_0_to_plus2
    use scheme_params, only: js3_eps, js3_p, z3_eps, z3_q, js5_eps, js5_p, z5_eps, z5_q
    use scheme_diagnostics_mod, only: scheme_diagnostics
    implicit none

contains

    ! Second-order MUSCL interpolation with the minmod limiter.
    subroutine muscl_2nd(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        if (present(diag)) diag = scheme_diagnostics()

        u_half = u(0) + 0.5_precision * minmod_limiter(u(1) - u(0), u(0) - u(-1))
    end subroutine muscl_2nd

    ! Third-order JS-type nonlinear weighted interpolation.
    subroutine JS3_interp(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C(1:2) = (/ 0.25_precision, 0.75_precision /)
        real(precision) :: ta(1:2), beta(1:2), alpha(1:2), omega(1:2)
        integer :: k

        ta(1) = -0.5_precision * u(-1) + 1.5_precision * u(0)
        ta(2) = 0.5_precision * u(0) + 0.5_precision * u(1)

        beta(1) = (u(0) - u(-1))**2
        beta(2) = (u(1) - u(0))**2

        do k = 1, 2
            alpha(k) = C(k) / (beta(k) + js3_eps)**js3_p
        end do

        alpha = alpha / sum(alpha)
        omega = alpha
        u_half = omega(1) * ta(1) + omega(2) * ta(2)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta(1:2) = beta
            diag%omega(1:2) = omega
        end if
    end subroutine JS3_interp

    ! Third-order Z-type nonlinear weighted interpolation.
    subroutine Z3_interp(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C(1:2) = (/ 0.25_precision, 0.75_precision /)
        real(precision) :: ta(1:2), beta(1:2), alpha(1:2), omega(1:2)
        real(precision) :: tau
        integer :: k

        ta(1) = -0.5_precision * u(-1) + 1.5_precision * u(0)
        ta(2) = 0.5_precision * u(0) + 0.5_precision * u(1)

        beta(1) = (u(0) - u(-1))**2
        beta(2) = (u(1) - u(0))**2
        tau = abs(beta(1) - beta(2))

        do k = 1, 2
            alpha(k) = C(k) * (1.0_precision + (tau / (beta(k) + z3_eps))**z3_q)
        end do

        alpha = alpha / sum(alpha)
        omega = alpha
        u_half = omega(1) * ta(1) + omega(2) * ta(2)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta(1:2) = beta
            diag%omega(1:2) = omega
            diag%tau = tau
        end if
    end subroutine Z3_interp

    ! Fifth-order JS-type nonlinear weighted interpolation.
    subroutine JS5_interp(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C(0:2) = (/ 0.0625_precision, 0.625_precision, 0.3125_precision /)
        real(precision) :: ta(0:2), beta(0:2), alpha(0:2)
        integer :: k

        call linear_3rd_minus2_to_0(u, ta(0))
        call linear_3rd_minus1_to_plus1(u, ta(1))
        call linear_3rd_0_to_plus2(u, ta(2))

        beta(0) = 0.25_precision * (-u(-2) + real(4, precision) * u(-1) - real(3, precision) * u(0))**2 &
                + (u(-2) - real(2, precision) * u(-1) + u(0))**2
        beta(1) = 0.25_precision * (u(-1) - u(1))**2 &
                + (u(-1) - real(2, precision) * u(0) + u(1))**2
        beta(2) = 0.25_precision * (real(3, precision) * u(0) - real(4, precision) * u(1) + u(2))**2 &
                + (u(0) - real(2, precision) * u(1) + u(2))**2

        do k = 0, 2
            alpha(k) = C(k) / (beta(k) + js5_eps)**js5_p
        end do

        alpha = alpha / sum(alpha)
        u_half = alpha(0) * ta(0) + alpha(1) * ta(1) + alpha(2) * ta(2)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta(1:3) = beta(0:2)
            diag%omega(1:3) = alpha(0:2)
        end if
    end subroutine JS5_interp

    ! Fifth-order Z-type nonlinear weighted interpolation.
    subroutine Z5_interp(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C(0:2) = (/ 0.0625_precision, 0.625_precision, 0.3125_precision /)
        real(precision) :: ta(0:2), beta(0:2), alpha(0:2)
        real(precision) :: tau
        integer :: k

        call linear_3rd_minus2_to_0(u, ta(0))
        call linear_3rd_minus1_to_plus1(u, ta(1))
        call linear_3rd_0_to_plus2(u, ta(2))

        beta(0) = 0.25_precision * (-u(-2) + real(4, precision) * u(-1) - real(3, precision) * u(0))**2 &
                + (u(-2) - real(2, precision) * u(-1) + u(0))**2
        beta(1) = 0.25_precision * (u(-1) - u(1))**2 &
                + (u(-1) - real(2, precision) * u(0) + u(1))**2
        beta(2) = 0.25_precision * (real(3, precision) * u(0) - real(4, precision) * u(1) + u(2))**2 &
                + (u(0) - real(2, precision) * u(1) + u(2))**2

        tau = abs(beta(0) - beta(2))

        do k = 0, 2
            alpha(k) = C(k) * (1.0_precision + (tau / (beta(k) + z5_eps))**z5_q)
        end do

        alpha = alpha / sum(alpha)
        u_half = alpha(0) * ta(0) + alpha(1) * ta(1) + alpha(2) * ta(2)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta(1:3) = beta(0:2)
            diag%omega(1:3) = alpha(0:2)
            diag%tau = tau
        end if
    end subroutine Z5_interp

    function sgn(x) result(s)
        real(precision), intent(in) :: x
        real(precision) :: s

        if (x >= 0.0_precision) then
            s = 1.0_precision
        else
            s = -1.0_precision
        end if
    end function sgn

    function minmod_limiter(x, y) result(m)
        real(precision), intent(in) :: x, y
        real(precision) :: m

        m = 0.5_precision * (sgn(x) + sgn(y)) * min(abs(x), abs(y))
    end function minmod_limiter

end module interp_traditional
