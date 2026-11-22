module interp_framework
    use precision_select
    use interp_linear, only: poly3, linear_3rd_minus1_to_plus1, linear_5th_minus2_to_plus2, &
                             linear_3rd_minus2_to_0, linear_3rd_0_to_plus2
    use interp_rbf, only: RBF_GA_Perturbed_3, RBF_MQ_Perturbed_3, RBF_IQ_Perturbed_3, &
                          RBF_GA_Perturbed_5, RBF_MQ_Perturbed_5, RBF_IQ_Perturbed_5
    use interp_qr, only: quasi_rational
    use scheme_params, only: moc322_eps, moc322_M_q, moc33_eps, moc33_p, moc33_linear, &
                             moc53_eps, moc53_N_q, moc53_p, moc53_linear, &
                             moz_eps, moz_N_q, moz_M_q
    use scheme_diagnostics_mod, only: scheme_diagnostics
    implicit none

    private
    public :: MOC322, MOC33, MOC53, MOZ

contains

    ! MOC322: a 3rd-order main candidate plus two 2nd-order candidates.
    ! The smoothness indicators follow the MOC formulation in the thesis:
    ! beta_00 is the global indicator, beta_L/beta_R are the two sub-template
    ! indicators.  The alpha/weight form is the MOC form, equivalent to the
    ! Z-type expression used by the original code.
    subroutine MOC322(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C_L = 0.25_precision
        real(precision), parameter :: C_R = 0.75_precision

        real(precision) :: eta(-1:1)
        real(precision) :: r, gain
        real(precision) :: beta_00, beta_L, beta_R
        real(precision) :: alpha_00, alpha_L, alpha_R, sum_alpha
        real(precision) :: omega_00, omega_L, omega_R
        real(precision) :: u_00, u_L, u_R
        real(precision) :: eps

        eps = moc322_eps

        eta(-1) = abs(u(0) - real(2, precision) * u(-1) + u(-2)) + eps
        eta(0)  = abs(u(1) - real(2, precision) * u(0) + u(-1)) + eps
        eta(1)  = abs(u(2) - real(2, precision) * u(1) + u(0)) + eps

        r = ((eta(1) - eta(-1)) / (eta(1) + eta(-1)))**2

        beta_00 = r * eta(0)
        gain = (1.0_precision - r) / (4.0_precision * r + eps)
        beta_L = abs(u(0) - u(-1)) + gain * eta(0)
        beta_R = abs(u(1) - u(0)) + gain * eta(0)

        alpha_00 = 1.0_precision / (beta_00 + eps)**moc322_M_q
        alpha_L  = C_L / (beta_L + eps)**moc322_M_q
        alpha_R  = C_R / (beta_R + eps)**moc322_M_q
        sum_alpha = alpha_00 + alpha_L + alpha_R

        omega_00 = alpha_00 / sum_alpha
        omega_L = alpha_L / sum_alpha
        omega_R = alpha_R / sum_alpha

        call poly3(u, u_00)
        u_L = -0.5_precision * u(-1) + 1.5_precision * u(0)
        u_R = 0.5_precision * u(0) + 0.5_precision * u(1)

        u_half = omega_00 * u_00 + omega_L * u_L + omega_R * u_R

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta_L = beta_00
            diag%beta(1) = beta_L
            diag%beta(2) = beta_R
            diag%omega_L = omega_00
            diag%omega(1) = omega_L
            diag%omega(2) = omega_R
        end if
    end subroutine MOC322

    ! MOC33: a three-point linear/RBF candidate plus a three-point QR
    ! candidate.  The nonlinear candidate is fixed to QR.
    subroutine MOC33(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: eta(-1:1)
        real(precision) :: kappa
        real(precision) :: beta_00, beta_H(0:1), beta_N
        real(precision) :: alpha_L, alpha_N, sum_p
        real(precision) :: omega_L, omega_N
        real(precision) :: u_L, u_N
        real(precision) :: eps
        integer :: k

        eps = moc33_eps

        eta(-1) = abs(u(0) - real(2, precision) * u(-1) + u(-2)) + eps
        eta(0)  = abs(u(1) - real(2, precision) * u(0) + u(-1)) + eps
        eta(1)  = abs(u(2) - real(2, precision) * u(1) + u(0)) + eps

        kappa = abs((eta(1) - eta(-1)) / (eta(1) + eta(-1)))

        beta_00 = kappa * eta(0)
        do k = 0, 1
            beta_H(k) = abs(u(k) - u(k - 1)) + 0.5_precision * (1.0_precision - kappa) * eta(0)
        end do
        beta_N = beta_H(0) * beta_H(1) / (beta_H(0) + beta_H(1) + eps)

        alpha_L = 1.0_precision / (beta_00 + eps)
        alpha_N = 1.0_precision / (beta_N + eps)
        sum_p = alpha_L**moc33_p + alpha_N**moc33_p
        omega_L = alpha_L**moc33_p / sum_p
        omega_N = alpha_N**moc33_p / sum_p

        select case (trim(moc33_linear))
        case ('GA')
            call RBF_GA_Perturbed_3(u, u_L)
        case ('MQ')
            call RBF_MQ_Perturbed_3(u, u_L)
        case ('IQ')
            call RBF_IQ_Perturbed_3(u, u_L)
        case ('LINEAR')
            call linear_3rd_minus1_to_plus1(u, u_L)
        case default
            write(*, '(a)') 'invalid moc33_linear: ' // trim(moc33_linear)
            error stop 1
        end select

        call quasi_rational(u, u_N)

        u_half = omega_L * u_L + omega_N * u_N

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta_L = beta_00
            diag%beta_N = beta_N
            diag%omega_L = omega_L
            diag%omega_N = omega_N
            diag%beta(1) = beta_00
            diag%beta(2) = beta_N
            diag%omega(1) = omega_L
            diag%omega(2) = omega_N
        end if
    end subroutine MOC33

    ! MOC53: a five-point linear/RBF candidate plus a three-point QR
    ! candidate.  The nonlinear candidate is fixed to QR.
    subroutine MOC53(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: tau, h_noscale
        real(precision) :: beta(0:2)
        real(precision) :: beta_00, beta_10
        real(precision) :: alpha_L, alpha_N, sum_p
        real(precision) :: omega_L, omega_N
        real(precision) :: u_L, u_N

        call calc_tau5_improved(u, tau, h_noscale, moc53_N_q, moc53_eps)
        call calc_five_point_betas(u, beta)

        beta_00 = tau
        beta_10 = beta(0) * beta(1) * beta(2) &
                / (beta(0) * beta(1) + beta(1) * beta(2) + beta(0) * beta(2) + moc53_eps)

        alpha_L = 1.0_precision / (beta_00 + moc53_eps)
        alpha_N = 1.0_precision / (beta_10 + moc53_eps)
        sum_p = alpha_L**moc53_p + alpha_N**moc53_p
        omega_L = alpha_L**moc53_p / sum_p
        omega_N = alpha_N**moc53_p / sum_p

        select case (trim(moc53_linear))
        case ('GA')
            call RBF_GA_Perturbed_5(u, u_L)
        case ('MQ')
            call RBF_MQ_Perturbed_5(u, u_L)
        case ('IQ')
            call RBF_IQ_Perturbed_5(u, u_L)
        case ('LINEAR')
            call linear_5th_minus2_to_plus2(u, u_L)
        case default
            write(*, '(a)') 'invalid moc53_linear: ' // trim(moc53_linear)
            error stop 1
        end select

        call quasi_rational(u, u_N)

        u_half = omega_L * u_L + omega_N * u_N

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta_L = beta_00
            diag%beta_N = beta_10
            diag%omega_L = omega_L
            diag%omega_N = omega_N
            diag%h_noscale = h_noscale
            diag%beta(1) = beta_00
            diag%beta(2) = beta_10
            diag%omega(1) = omega_L
            diag%omega(2) = omega_N
        end if
    end subroutine MOC53

    ! MOZ: five-point Z-type weighting with the improved global indicator.
    subroutine MOZ(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision), parameter :: C(0:2) = (/ 0.0625_precision, 0.625_precision, 0.3125_precision /)
        real(precision) :: tau, h_noscale
        real(precision) :: beta(0:2), alpha(0:2), ta(0:2)
        integer :: k

        call calc_tau5_improved(u, tau, h_noscale, moz_N_q, moz_eps)
        call calc_five_point_betas(u, beta)

        do k = 0, 2
            alpha(k) = C(k) * (1.0_precision + (tau / (beta(k) + moz_eps))**moz_M_q)
        end do

        alpha = alpha / sum(alpha)

        call linear_3rd_minus2_to_0(u, ta(0))
        call linear_3rd_minus1_to_plus1(u, ta(1))
        call linear_3rd_0_to_plus2(u, ta(2))

        u_half = alpha(0) * ta(0) + alpha(1) * ta(1) + alpha(2) * ta(2)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%beta(1:3) = beta(0:2)
            diag%omega(1:3) = alpha(0:2)
            diag%tau = tau
            diag%h_noscale = h_noscale
        end if
    end subroutine MOZ

    ! Improved global smoothness indicator of MOZ/MOC53.
    subroutine calc_tau5_improved(u, tau, h_noscale, N_q, eps)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: tau, h_noscale
        integer(integer_precision), intent(in) :: N_q
        real(precision), intent(in) :: eps

        real(precision) :: dL, dR, num, den

        dL = u(1) - real(3, precision) * u(0) + real(3, precision) * u(-1) - u(-2)
        dR = u(2) - real(3, precision) * u(1) + real(3, precision) * u(0) - u(-1)

        num = dL * dL - dR * dR
        den = dL * dL + dR * dR

        h_noscale = abs(num / (den + eps))
        tau = h_noscale**N_q * abs(num)
    end subroutine calc_tau5_improved

    ! Three three-point sub-template smoothness indicators of the
    ! five-point Z5/MOZ/MOC53 stencil.
    subroutine calc_five_point_betas(u, beta)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: beta(0:2)

        beta(0) = 0.25_precision * (-u(-2) + real(4, precision) * u(-1) - real(3, precision) * u(0))**2 &
                + (u(-2) - real(2, precision) * u(-1) + u(0))**2
        beta(1) = 0.25_precision * (u(-1) - u(1))**2 &
                + (u(-1) - real(2, precision) * u(0) + u(1))**2
        beta(2) = 0.25_precision * (real(3, precision) * u(0) - real(4, precision) * u(1) + u(2))**2 &
                + (u(0) - real(2, precision) * u(1) + u(2))**2
    end subroutine calc_five_point_betas

end module interp_framework
