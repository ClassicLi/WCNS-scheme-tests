module weno_schemes
    ! WENO reconstruction routines for the ADR test.
    use precision_select
    use weno_params, only: js3_eps, js3_p, js5_eps, js5_p, z5_eps, z5_q, &
                           weno_c_c4, weno_c_c5, weno_c_eps, llr_eps, rbf_eps, &
                           teno_q, teno_ct, thinc_beta
    implicit none

contains

    ! Five-point WENO optimal linear reconstruction at j+1/2.
    subroutine weno5_linear_reconstruct(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision), parameter :: a = 2.0_precision
        real(precision), parameter :: b = -13.0_precision
        real(precision), parameter :: c = 47.0_precision
        real(precision), parameter :: d = 27.0_precision
        real(precision), parameter :: e = -3.0_precision
        real(precision), parameter :: denominator = 60.0_precision

        flux_half = (a*u(-2) + b*u(-1) + c*u(0) + d*u(1) + e*u(2))/denominator
    end subroutine weno5_linear_reconstruct

    ! Three-point linear reconstruction at j+1/2.
    subroutine weno3_linear(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: ideal_w(0:1)
        real(precision) :: tilde_f(0:1)

        ideal_w = [1.0_precision/3.0_precision, 2.0_precision/3.0_precision]

        tilde_f(0) = -0.5_precision * u(-1) + 1.5_precision * u(0)
        tilde_f(1) = 0.5_precision * u(0) + 0.5_precision * u(1)

        flux_half = ideal_w(0) * tilde_f(0) + ideal_w(1) * tilde_f(1)
    end subroutine weno3_linear

    ! Three-point JS-type nonlinear weighted reconstruction.
    subroutine weno3_js(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: beta(0:1)
        real(precision) :: tilde_f(0:1)
        real(precision) :: ideal_w(0:1)
        real(precision) :: alpha(0:1)
        real(precision) :: omega(0:1)

        ideal_w = [1.0_precision/3.0_precision, 2.0_precision/3.0_precision]

        beta(0) = (u(-1) - u(0))**2.0_precision
        beta(1) = (u(0) - u(1))**2.0_precision

        tilde_f(0) = -0.5_precision * u(-1) + 1.5_precision * u(0)
        tilde_f(1) = 0.5_precision * u(0) + 0.5_precision * u(1)

        alpha(0) = ideal_w(0) / ((beta(0) + js3_eps)**js3_p)
        alpha(1) = ideal_w(1) / ((beta(1) + js3_eps)**js3_p)

        omega(0) = alpha(0) / (alpha(0) + alpha(1))
        omega(1) = alpha(1) / (alpha(0) + alpha(1))

        flux_half = omega(0) * tilde_f(0) + omega(1) * tilde_f(1)
    end subroutine weno3_js

    ! Five-point JS-type nonlinear weighted reconstruction.
    subroutine weno5_js(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: tilde_a(0:2)
        real(precision) :: ideal_weight(0:2)
        real(precision) :: js_weight(0:2)
        real(precision) :: beta(0:2)
        real(precision) :: alpha(0:2)
        real(precision) :: sum_alpha
        integer(integer_precision) :: loop

        ideal_weight = [1.0_precision/10.0_precision, 6.0_precision/10.0_precision, &
                        3.0_precision/10.0_precision]

        tilde_a(0) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0)) / 6.0_precision
        tilde_a(1) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1)) / 6.0_precision
        tilde_a(2) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2)) / 6.0_precision

        beta(0) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)
        beta(1) = 0.25_precision*((u(1) - u(-1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(2) = 0.25_precision*((-3.0_precision*u(0) + 4.0_precision*u(1) - u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)

        do loop = 0, 2
            alpha(loop) = ideal_weight(loop) / ((beta(loop) + js5_eps)**js5_p)
        end do

        sum_alpha = sum(alpha)
        js_weight = alpha / sum_alpha

        flux_half = js_weight(0)*tilde_a(0) + js_weight(1)*tilde_a(1) + js_weight(2)*tilde_a(2)
    end subroutine weno5_js

    ! Five-point WENO-M reconstruction with mapped weights.
    subroutine weno5_m(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: tilde_a(0:2)
        real(precision) :: ideal_weight(0:2)
        real(precision) :: js_weight(0:2)
        real(precision) :: beta(0:2)
        real(precision) :: alpha(0:2)
        real(precision) :: sum_alpha
        real(precision) :: m_weight(0:2)
        integer(integer_precision) :: loop

        ideal_weight = [1.0_precision/10.0_precision, 6.0_precision/10.0_precision, &
                        3.0_precision/10.0_precision]

        tilde_a(0) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0)) / 6.0_precision
        tilde_a(1) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1)) / 6.0_precision
        tilde_a(2) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2)) / 6.0_precision

        beta(0) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)
        beta(1) = 0.25_precision*((u(1) - u(-1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(2) = 0.25_precision*((-3.0_precision*u(0) + 4.0_precision*u(1) - u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)

        do loop = 0, 2
            alpha(loop) = ideal_weight(loop) / ((beta(loop) + js5_eps)**js5_p)
        end do

        sum_alpha = sum(alpha)
        js_weight = alpha / sum_alpha

        call m_weight_5_func(js_weight, ideal_weight, m_weight)

        flux_half = m_weight(0)*tilde_a(0) + m_weight(1)*tilde_a(1) + m_weight(2)*tilde_a(2)
    end subroutine weno5_m

    ! Map JS weights to WENO-M weights.
    subroutine m_weight_5_func(omega_js, ideal_weight, omega_m)
        real(precision), intent(in) :: omega_js(0:2)
        real(precision), intent(in) :: ideal_weight(0:2)
        real(precision), intent(out) :: omega_m(0:2)

        integer(integer_precision) :: loop

        do loop = 0, 2
            omega_m(loop) = omega_js(loop) * (ideal_weight(loop) + ideal_weight(loop)*ideal_weight(loop) &
                - 3.0_precision*ideal_weight(loop)*omega_js(loop) + omega_js(loop)*omega_js(loop)) &
                / (ideal_weight(loop)*ideal_weight(loop) + omega_js(loop) &
                - 2.0_precision*ideal_weight(loop)*omega_js(loop))
        end do
    end subroutine m_weight_5_func

    ! Five-point Z-type nonlinear weighted reconstruction.
    subroutine weno5_z(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: tilde_a(0:2)
        real(precision) :: ideal_weight(0:2)
        real(precision) :: z_weight(0:2)
        real(precision) :: beta(0:2)
        real(precision) :: alpha(0:2)
        real(precision) :: sum_alpha
        real(precision) :: tau
        integer(integer_precision) :: loop

        ideal_weight = [1.0_precision/10.0_precision, 6.0_precision/10.0_precision, &
                        3.0_precision/10.0_precision]

        tilde_a(0) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0)) / 6.0_precision
        tilde_a(1) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1)) / 6.0_precision
        tilde_a(2) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2)) / 6.0_precision

        beta(0) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)
        beta(1) = 0.25_precision*((u(1) - u(-1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(2) = 0.25_precision*((-3.0_precision*u(0) + 4.0_precision*u(1) - u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)

        tau = abs(beta(2) - beta(0))

        do loop = 0, 2
            alpha(loop) = ideal_weight(loop) * (1.0_precision + (tau/(beta(loop) + z5_eps))**z5_q)
        end do

        sum_alpha = sum(alpha)
        z_weight = alpha / sum_alpha

        flux_half = z_weight(0)*tilde_a(0) + z_weight(1)*tilde_a(1) + z_weight(2)*tilde_a(2)
    end subroutine weno5_z

    ! Three-point trigonometric reconstruction (experimental normalization).
    subroutine weno3_trig(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: betah, cos_h, cos_2h

        betah = h / (sin(2.0_precision*h) * (2.0_precision + 4.0_precision*cos(h)))
        cos_h = cos(h)
        cos_2h = cos(2.0_precision*h)

        flux_half = 2.0_precision * betah * (2.0_precision*u(0) + u(1) &
                    + (-u(-1) + 2.0_precision*u(0) + u(1))*cos_h + u(0)*cos_2h)
    end subroutine weno3_trig

    ! Auxiliary function eta_+(a) used by LLR and DLR.
    function eta_plus(a) result(e)
        real(precision), intent(in) :: a
        real(precision) :: e

        e = -(log(1.0_precision - a) + a)/a/a
        if (a == 0.0_precision) e = 0.5_precision
    end function eta_plus

    ! Three-point LLR reconstruction.
    subroutine weno3_llr(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: d1, d2, a

        d1 = (u(0) - u(-1)) / h
        d2 = (u(1) - u(0)) / h

        a = (d2 - d1) / (d2 + llr_eps)

        flux_half = u(0) + d1 * h * eta_plus(a)

        if (a >= 1.0_precision) then
            flux_half = -1.0_precision/6.0_precision*u(-1) + 5.0_precision/6.0_precision*u(0) &
                        + 1.0_precision/3.0_precision*u(1)
        end if
    end subroutine weno3_llr

    ! Baeza-type third-order reconstruction.
    subroutine baeza_weno3(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: beta(0:2)
        real(precision), parameter :: vareps = 1.0e-40_precision
        real(precision) :: tilde_f(0:1)
        real(precision) :: ideal_w(0:1)
        real(precision) :: tiled_omega(0:1)
        real(precision) :: corrector
        real(precision) :: jj, tau, dd
        real(precision) :: omega(0:1)

        ideal_w = [1.0_precision/3.0_precision, 2.0_precision/3.0_precision]

        tilde_f(0) = -0.5_precision*u(-1) + 1.5_precision*u(0)
        tilde_f(1) = 0.5_precision*u(0) + 0.5_precision*u(1)

        beta(0) = (u(0) - u(-1))**2.0_precision
        beta(1) = (u(1) - u(0))**2.0_precision
        beta(2) = (u(2) - u(1))**2.0_precision

        jj = beta(0)*(beta(1) + beta(2)) + beta(2)*(beta(0) + beta(1))
        dd = (-u(-1) + 3.0_precision*u(0) - 3.0_precision*u(1) + u(2))**2.0_precision
        tau = dd*(beta(0) + beta(1) + beta(2))

        corrector = jj/(jj + tau + vareps)

        tiled_omega(0) = (beta(1) + vareps)/(beta(0) + beta(1) + 2.0_precision*vareps)
        tiled_omega(1) = (beta(0) + vareps)/(beta(0) + beta(1) + 2.0_precision*vareps)

        omega(0) = corrector*ideal_w(0) + (1.0_precision - corrector)*tiled_omega(0)
        omega(1) = corrector*ideal_w(1) + (1.0_precision - corrector)*tiled_omega(1)

        flux_half = omega(0)*tilde_f(0) + omega(1)*tilde_f(1)
    end subroutine baeza_weno3

    ! WENO-C reconstruction with 3/4/5-point candidates.
    subroutine weno_c_deheilan(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: A_3_2(0:2, 1:2), beta_3_r(0:2)
        real(precision) :: A_4_3(1:2, 1:3), beta_4_r(1:2)
        real(precision) :: A_5(1:4), beta_5
        real(precision) :: ideal_3(0:2), ideal_4(1:2)
        real(precision) :: flux_sub3(0:2), flux_sub4(1:2), flux_sub5
        real(precision) :: tilde_flux_4, tilde_flux_3
        real(precision) :: omega_3(0:2), alpha_3(0:2), sum_alpha_3
        real(precision) :: omega_4(1:2), alpha_4(1:2), sum_alpha_4
        real(precision) :: delta_4(1:2)
        real(precision) :: tilde_beta3, tilde_beta4, tilde_beta5
        real(precision) :: gama(1:3), tilde_alpha(1:3)
        integer(integer_precision) :: loop

        ideal_3 = [0.3_precision, 0.6_precision, 0.1_precision]
        ideal_4 = [0.6_precision, 0.4_precision]

        flux_sub3(0) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2))/6.0_precision
        flux_sub3(1) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1))/6.0_precision
        flux_sub3(2) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0))/6.0_precision

        flux_sub4(1) = (-u(-1) + 7.0_precision*u(0) + 7.0_precision*u(1) - u(2))/12.0_precision
        flux_sub4(2) = (u(-2) - 5.0_precision*u(-1) + 13.0_precision*u(0) + 3.0_precision*u(1))/12.0_precision

        flux_sub5 = ideal_3(0)*flux_sub3(0) + ideal_3(1)*flux_sub3(1) + ideal_3(2)*flux_sub3(2)

        A_3_2(0,1) = (3.0_precision*u(0) - 4.0_precision*u(1) + u(2))/2.0_precision
        A_3_2(0,2) = (u(0) - 2.0_precision*u(1) + u(2))/2.0_precision
        A_3_2(1,1) = (u(-1) - u(1))/2.0_precision
        A_3_2(1,2) = (u(-1) - 2.0_precision*u(0) + u(1))/2.0_precision
        A_3_2(2,1) = (u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))/2.0_precision
        A_3_2(2,2) = (u(-2) - 2.0_precision*u(-1) + u(0))/2.0_precision

        do loop = 0, 2
            beta_3_r(loop) = A_3_2(loop,1)**2.0_precision &
                + 13.0_precision*A_3_2(loop,2)**2.0_precision/3.0_precision
        end do

        A_4_3(1,1) = (-19.0_precision*u(-1) - 33.0_precision*u(0) + 63.0_precision*u(1) &
                      - 11.0_precision*u(2))/60.0_precision
        A_4_3(1,2) = (u(-1) - 2.0_precision*u(0) + u(1))/2.0_precision
        A_4_3(1,3) = (-u(-1) + 3.0_precision*u(0) - 3.0_precision*u(1) + u(2))/6.0_precision
        A_4_3(2,1) = (11.0_precision*u(-2) - 63.0_precision*u(-1) + 33.0_precision*u(0) &
                      + 19.0_precision*u(1))/60.0_precision
        A_4_3(2,2) = (u(-1) - 2.0_precision*u(0) + u(1))/2.0_precision
        A_4_3(2,3) = (-u(-2) + 3.0_precision*u(-1) - 3.0_precision*u(0) + u(1))/6.0_precision

        do loop = 1, 2
            beta_4_r(loop) = (A_4_3(loop,1) + A_4_3(loop,3)/10.0_precision)**2.0_precision &
                + 13.0_precision*A_4_3(loop,2)**2.0_precision/3.0_precision &
                + 781.0_precision*A_4_3(loop,3)**2.0_precision/20.0_precision
        end do

        A_5(1) = (11.0_precision*u(-2) - 82.0_precision*u(-1) + 82.0_precision*u(1) &
                  - 11.0_precision*u(2))/120.0_precision
        A_5(2) = (-3.0_precision*u(-2) + 40.0_precision*u(-1) - 74.0_precision*u(0) &
                  + 40.0_precision*u(1) - 3.0_precision*u(2))/56.0_precision
        A_5(3) = (-u(-2) + 2.0_precision*u(-1) - 2.0_precision*u(1) + u(2))/12.0_precision
        A_5(4) = (u(-2) - 4.0_precision*u(-1) + 6.0_precision*u(0) - 4.0_precision*u(1) &
                  + u(2))/24.0_precision

        beta_5 = (A_5(1) + A_5(3)/10.0_precision)**2.0_precision &
            + (13.0_precision/3.0_precision)*(A_5(2) + 123.0_precision/455.0_precision*A_5(4))**2.0_precision &
            + (781.0_precision/20.0_precision)*A_5(3)**2.0_precision &
            + (1421461.0_precision/2275.0_precision)*A_5(4)**2.0_precision

        do loop = 0, 2
            alpha_3(loop) = ideal_3(loop)/(beta_3_r(loop) + weno_c_eps)**2.0_precision
        end do
        sum_alpha_3 = sum(alpha_3)
        do loop = 0, 2
            omega_3(loop) = alpha_3(loop)/sum_alpha_3
        end do
        tilde_beta3 = sum(omega_3*beta_3_r)
        tilde_flux_3 = sum(omega_3*flux_sub3)

        do loop = 1, 2
            alpha_4(loop) = ideal_4(loop)/(beta_4_r(loop) + weno_c_eps)**2.0_precision
        end do
        sum_alpha_4 = alpha_4(1) + alpha_4(2)
        do loop = 1, 2
            omega_4(loop) = alpha_4(loop)/sum_alpha_4
        end do
        tilde_flux_4 = omega_4(1)*flux_sub4(1) + omega_4(2)*flux_sub4(2)

        do loop = 1, 2
            delta_4(loop) = (weno_c_c4/2.0_precision)*(abs(beta_4_r(loop) - beta_3_r(loop-1)) &
                + abs(beta_4_r(loop) - beta_3_r(loop)))
        end do
        tilde_beta4 = omega_4(1)*delta_4(1) + omega_4(2)*delta_4(2)

        tilde_beta5 = (weno_c_c5/3.0_precision)*(abs(beta_5 - beta_4_r(1)) + abs(beta_5 - beta_4_r(2)) &
            + abs(6.0_precision*beta_5 - (beta_3_r(0) + 4.0_precision*beta_3_r(1) &
            + beta_3_r(2)))/6.0_precision)

        tilde_alpha(1) = 1.0_precision/(tilde_beta3 + weno_c_eps)**2.0_precision
        tilde_alpha(2) = 1.0_precision/(tilde_beta4 + weno_c_eps)**2.0_precision
        tilde_alpha(3) = 1.0_precision/(tilde_beta5 + weno_c_eps)**2.0_precision

        gama = tilde_alpha/sum(tilde_alpha)

        flux_half = gama(1)*tilde_flux_3 + gama(2)*tilde_flux_4 + gama(3)*flux_sub5
    end subroutine weno_c_deheilan

    ! Fifth-order DLR reconstruction (experimental; see formula_audit).
    subroutine weno5_dlr(u, flux_half, h)
        real(precision), intent(in) :: u(-4:4)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: d_minus, d_plus, dd_minus, dd_plus
        real(precision) :: a, b, c, d
        real(precision) :: a0, a1, a2, qa, q0
        real(precision), parameter :: tol = 1.0_precision

        call dlr_negative_nonlinear(u(-4:3), h, d_minus, dd_minus)
        call dlr_negative_nonlinear(u(-3:4), h, d_plus, dd_plus)

        a0 = get_dlr_a0(h, d_minus, d_plus, dd_minus, dd_plus)

        if (abs(a0) < tol .and. abs(d_minus) > h) then
            dd_plus = dd_plus + h*h*h
        else if (abs(a0) < tol .and. abs(d_minus) <= h) then
            d_plus = d_plus + h*h*h*h
        end if

        a1 = get_dlr_a1(h, d_minus, d_plus, dd_minus, dd_plus)
        a2 = get_dlr_a2(h, d_minus, d_plus, dd_minus, dd_plus)

        a = a1 + a2
        b = a1 - a2

        qa = dlr_q(a, d_minus, d_plus, dd_minus, dd_plus, h)

        d = d_minus/2.0_precision - qa
        c = d_minus/2.0_precision + qa

        q0 = 2.0_precision*(d_plus*(a - 1.0_precision)**2.0_precision &
             + d_minus*(2.0_precision*a - 1.0_precision) - h*dd_minus)

        if (abs(q0) < tol) then
            dd_minus = dd_minus - h*h*h
        end if

        flux_half = u(0) + c*h*eta_plus(a) + d*h*eta_plus(b)
    end subroutine weno5_dlr

    ! DLR negative-side nonlinear derivative estimates.
    subroutine dlr_negative_nonlinear(u, h, d, dd)
        real(precision), intent(in) :: u(-4:3)
        real(precision), intent(in) :: h
        real(precision), intent(out) :: d, dd

        real(precision) :: tilde_d(-1:1)
        real(precision) :: tilde_dd(-1:1)
        real(precision), parameter :: ideal_w_d(-1:1) = [1.0_precision/24.0_precision, &
            1.0_precision/24.0_precision, 11.0_precision/12.0_precision]
        real(precision), parameter :: ideal_w_dd(-1:1) = [1.0_precision/16.0_precision, &
            1.0_precision/16.0_precision, 7.0_precision/8.0_precision]
        real(precision) :: beta(-1:1), alpha(-1:1)
        real(precision) :: omega_d(-1:1), omega_dd(-1:1)
        real(precision) :: denominator_omega_d, denominator_omega_dd
        real(precision) :: vareps
        integer(integer_precision) :: loop

        vareps = 1.0e-8_precision

        tilde_d(-1) = u(-3) - 3.0_precision*u(-2) + 2.0_precision*u(-1)
        tilde_d(0) = -u(-1) + u(0)
        tilde_d(1) = -2.0_precision*u(0) + 3.0_precision*u(1) - u(2)
        tilde_d = tilde_d / h

        tilde_dd(-1) = -3.0_precision*u(-4) + 11.0_precision*u(-3) - 13.0_precision*u(-2) &
                       + 5.0_precision*u(-1)
        tilde_dd(0) = u(-2) - u(-1) - u(0) + u(1)
        tilde_dd(1) = 5.0_precision*u(0) - 13.0_precision*u(1) + 11.0_precision*u(2) &
                      - 3.0_precision*u(3)
        tilde_dd = tilde_dd / h / h / 2.0_precision

        do loop = -1, 1
            beta(loop) = (h*tilde_d(loop))**2.0_precision + (h*h*tilde_dd(loop))**2.0_precision
        end do

        do loop = -1, 1
            alpha(loop) = h*h / (beta(loop) + vareps)
        end do

        denominator_omega_d = ideal_w_d(-1)*alpha(-1) + ideal_w_d(0)*alpha(0) + ideal_w_d(1)*alpha(1)
        denominator_omega_dd = ideal_w_dd(-1)*alpha(-1) + ideal_w_dd(0)*alpha(0) + ideal_w_dd(1)*alpha(1)

        do loop = -1, 1
            omega_d(loop) = ideal_w_d(loop)*alpha(loop)/denominator_omega_d
        end do

        do loop = -1, 1
            omega_dd(loop) = ideal_w_dd(loop)*alpha(loop)/denominator_omega_dd
        end do

        d = omega_d(-1)*tilde_d(-1) + omega_d(0)*tilde_d(0) + omega_d(1)*tilde_d(1)
        dd = omega_dd(-1)*tilde_dd(-1) + omega_dd(0)*tilde_dd(0) + omega_dd(1)*tilde_dd(1)
    end subroutine dlr_negative_nonlinear

    ! DLR q coefficient.
    function dlr_q(z, d1, d2, dd1, dd2, h) result(q)
        real(precision), intent(in) :: z, h, d1, d2, dd1, dd2
        real(precision) :: q
        real(precision) :: denominator, numerator

        numerator = d1*d2*(1.0_precision - z*z) - 2.0_precision*h*d2*dd1*(1.0_precision - z) &
                    + d1*(h*dd1 - d1)
        denominator = d2*(1.0_precision - z)**2.0_precision + d1*(2.0_precision*z - 1.0_precision) &
                      - h*dd1
        q = 0.5_precision*numerator/denominator
    end function dlr_q

    ! Auxiliary function eta_-(a) used by DLR.
    function eta_minus(a) result(e)
        real(precision), intent(in) :: a
        real(precision) :: e

        e = ((a - 1.0_precision)*log(1.0_precision - a) - a)/a/a
        if (a == 0.0_precision) e = -0.5_precision
        if (a == 1.0_precision) e = -1.0_precision
    end function eta_minus

    ! DLR auxiliary coefficient a0.
    function get_dlr_a0(h, d1, d2, dd1, dd2) result(a0)
        real(precision), intent(in) :: h, d1, d2, dd1, dd2
        real(precision) :: a0

        a0 = 2.0_precision*(d2*(d1 - d2) + h*d1*dd2)
    end function get_dlr_a0

    ! DLR auxiliary coefficient a1.
    function get_dlr_a1(h, d1, d2, dd1, dd2) result(a1)
        real(precision), intent(in) :: h, d1, d2, dd1, dd2
        real(precision) :: a1

        a1 = (dd1*dd2*h*h + (d1*dd2 + d2*dd1)*h + 2.0_precision*d2*(d1 - d2)) &
             / get_dlr_a0(h, d1, d2, dd1, dd2)
    end function get_dlr_a1

    ! DLR auxiliary coefficient a2.
    function get_dlr_a2(h, d1, d2, dd1, dd2) result(a2)
        real(precision), intent(in) :: h, d1, d2, dd1, dd2
        real(precision) :: a2
        real(precision) :: temp1, temp2, temp3, temp4

        temp1 = dd1*dd1*dd2*dd2*h*h*h*h + 2.0_precision*(d2*dd1*dd1*dd2 - d1*dd1*dd2*dd2)*h*h*h
        temp2 = (d1*d1*dd2*dd2 + 2.0_precision*d1*d2*dd1*dd2 + d2*d2*dd1*dd1)*h*h
        temp3 = 4.0_precision*(d1*d1*d1*dd2 - d1*d1*d2*dd2 + d1*d2*d2*dd1 - d2*d2*d2*dd1)*h
        temp4 = 4.0_precision*d1*d2*(d1 - d2)*(d1 - d2)

        a2 = sqrt(temp1 + temp2 + temp3 + temp4)/get_dlr_a0(h, d1, d2, dd1, dd2)
    end function get_dlr_a2

    ! Three-point GA-RBF reconstruction with data-defined epsilon.
    subroutine weno3_ga_rbf_eps_define(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: eps2h2, temp1, temp2
        real(precision), parameter :: smv = 1.0e-30_precision

        eps2h2 = (-u(-1) + 2.0_precision*u(0) - u(1)) &
                 / (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1) + smv)

        temp1 = 2.0_precision*eps2h2*exp(-eps2h2)/(1.0_precision - exp(-4.0_precision*eps2h2))
        temp2 = temp1

        flux_half = temp1*u(0) + temp2*u(1)
    end subroutine weno3_ga_rbf_eps_define

    ! Three-point GA-RBF reconstruction with fixed epsilon.
    subroutine weno3_ga_rbf_eps_undefine(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: eps, eps2h2, temp1, temp2

        eps = rbf_eps
        eps2h2 = eps*eps*h*h

        temp1 = 2.0_precision*eps2h2*exp(-eps2h2)/(1.0_precision - exp(-4.0_precision*eps2h2))
        temp2 = temp1

        flux_half = temp1*u(0) + temp2*u(1)
    end subroutine weno3_ga_rbf_eps_undefine

    ! Five-point linear reconstruction plus GA-RBF perturbation.
    subroutine rbf_5_linear(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: eps, eps2h2
        real(precision) :: coe(-2:2), true_coe(-2:2)
        real(precision) :: coe_ga(-2:2), true_coe_ga(-2:2)
        real(precision) :: temp1(-2:2), temp2(-2:2)
        integer(integer_precision) :: loop

        coe = [-2.0_precision, 13.0_precision, -47.0_precision, -27.0_precision, 3.0_precision]
        true_coe = coe/(-60.0_precision)

        coe_ga = [-12.0_precision, 18.0_precision, 18.0_precision, -42.0_precision, 18.0_precision]
        true_coe_ga = coe_ga/(-60.0_precision)

        call fulin_define_eps(eps, u(-2:3))
        eps2h2 = eps*eps*h*h

        do loop = -2, 2
            temp1(loop) = true_coe(loop)*u(loop)
            temp2(loop) = true_coe_ga(loop)*u(loop)*eps2h2
        end do

        flux_half = sum(temp1) + sum(temp2)
    end subroutine rbf_5_linear

    ! Data-defined RBF epsilon from the local smoothness measure.
    subroutine fulin_define_eps(eps, u)
        real(precision), intent(out) :: eps
        real(precision), intent(in) :: u(-2:3)

        real(precision) :: eta(-1:1)
        real(precision) :: mm, theta, gm
        real(precision) :: alpha1, alpha2
        integer(integer_precision) :: loop

        alpha1 = 3.0_precision
        alpha2 = 50.0_precision

        do loop = -1, 1
            eta(loop) = get_eta(u(loop-1:loop+1))
        end do

        theta = min(eta(-1), eta(0), eta(1))
        mm = 1.0_precision - min(1.0_precision, theta/0.4_precision)
        gm = (1.0_precision - mm)**4.0_precision*(1.0_precision + 4.0_precision*mm)

        eps = alpha1 + alpha2*(1.0_precision - gm)
    end subroutine fulin_define_eps

    ! Local smoothness ratio used by fulin_define_eps.
    function get_eta(u) result(e)
        real(precision), intent(in) :: u(-1:1)
        real(precision) :: e
        real(precision), parameter :: smv = 1.0e-40_precision
        real(precision) :: denominator, numerator

        denominator = (u(0) - u(-1))**2.0_precision + (u(1) - u(0))**2.0_precision + smv
        numerator = abs(2.0_precision*(u(1) - u(0))*(u(0) - u(-1))) + smv
        e = numerator/denominator
    end function get_eta

    ! Original TENO 3-3-3 reconstruction.
    subroutine teno_original_333(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: beta(0:2), tau5
        real(precision) :: gama(0:2), sum_gama
        real(precision) :: chi(0:2), delta_k(0:2)
        real(precision) :: weight(0:2)
        real(precision) :: tilde_f(0:2)
        real(precision), parameter :: smv = 10.0e-40_precision
        real(precision), parameter :: ideal_weight(0:2) = [3.0_precision/5.0_precision, &
            3.0_precision/10.0_precision, 1.0_precision/10.0_precision]
        integer(integer_precision) :: loop

        beta(0) = 0.25_precision*((u(-1) - u(1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(1) = 0.25_precision*((3.0_precision*u(0) - 4.0_precision*u(1) + u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)
        beta(2) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)

        tau5 = abs(beta(1) - beta(2))

        tilde_f(0) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1))/6.0_precision
        tilde_f(1) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2))/6.0_precision
        tilde_f(2) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0))/6.0_precision

        do loop = 0, 2
            gama(loop) = (1.0_precision + tau5/(beta(loop) + smv))**teno_q
        end do

        sum_gama = sum(gama)

        do loop = 0, 2
            chi(loop) = gama(loop)/sum_gama
        end do

        do loop = 0, 2
            delta_k(loop) = teno_cutoff(chi(loop))
        end do

        do loop = 0, 2
            weight(loop) = delta_k(loop)*ideal_weight(loop) &
                           / sum(delta_k*ideal_weight)
        end do

        flux_half = weight(0)*tilde_f(0) + weight(1)*tilde_f(1) + weight(2)*tilde_f(2)
    end subroutine teno_original_333

    ! TENO cutoff function.
    function teno_cutoff(x) result(c)
        real(precision), intent(in) :: x
        real(precision) :: c

        if (x < teno_ct) then
            c = 0.0_precision
        else
            c = 1.0_precision
        end if
    end function teno_cutoff

    ! TENO 3-3-3 with RBF perturbation (unvalidated).
    subroutine teno_original_333_rbf(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: beta(0:2), tau5
        real(precision) :: gama(0:2), sum_gama
        real(precision) :: chi(0:2), delta_k(0:2)
        real(precision) :: weight(0:2)
        real(precision) :: tilde_f(0:2), rbf_append(0:2)
        real(precision) :: eps, eps2h2
        real(precision), parameter :: smv = 10.0e-40_precision
        real(precision), parameter :: ideal_weight(0:2) = [3.0_precision/5.0_precision, &
            3.0_precision/10.0_precision, 1.0_precision/10.0_precision]
        integer(integer_precision) :: loop

        beta(0) = 0.25_precision*((u(-1) - u(1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(1) = 0.25_precision*((3.0_precision*u(0) - 4.0_precision*u(1) + u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)
        beta(2) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)

        tau5 = abs(beta(1) - beta(2))

        tilde_f(0) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1))/6.0_precision
        tilde_f(1) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2))/6.0_precision
        tilde_f(2) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0))/6.0_precision

        call fulin_define_eps(eps, u(-2:3))

        rbf_append(0) = eps2h2*(-2.0_precision*u(-1) - 2.0_precision*u(0) &
                                + 4.0_precision*u(1))/6.0_precision
        rbf_append(1) = eps2h2*(-4.0_precision*u(0) + 2.0_precision*u(1) &
                                + 2.0_precision*u(2))/(-6.0_precision)
        rbf_append(2) = eps2h2*(-18.0_precision*u(-1) + 18.0_precision*u(0))/(-6.0_precision)

        tilde_f = tilde_f + rbf_append

        do loop = 0, 2
            gama(loop) = (1.0_precision + tau5/(beta(loop) + smv))**teno_q
        end do

        sum_gama = sum(gama)

        do loop = 0, 2
            chi(loop) = gama(loop)/sum_gama
        end do

        do loop = 0, 2
            delta_k(loop) = teno_cutoff(chi(loop))
        end do

        do loop = 0, 2
            weight(loop) = delta_k(loop)*ideal_weight(loop) &
                           / sum(delta_k*ideal_weight)
        end do

        flux_half = weight(0)*tilde_f(0) + weight(1)*tilde_f(1) + weight(2)*tilde_f(2)
    end subroutine teno_original_333_rbf

    ! Modified TENO semi-ENO reconstruction.
    subroutine teno_mody_semi_eno(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: beta(0:2), tau5
        real(precision) :: gama(0:2), sum_gama
        real(precision) :: chi(0:2), delta_k(0:2)
        real(precision) :: tilde_f(1:6)
        real(precision), parameter :: smv = 10.0e-40_precision
        integer(integer_precision) :: loop

        beta(0) = 0.25_precision*((u(-1) - u(1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(1) = 0.25_precision*((3.0_precision*u(0) - 4.0_precision*u(1) + u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)
        beta(2) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)

        tau5 = abs(beta(1) - beta(2))

        tilde_f(1) = (-2.0_precision*u(-2) + 13.0_precision*u(-1) - 47.0_precision*u(0) &
                      - 27.0_precision*u(1) + 3.0_precision*u(2))/(-60.0_precision)
        tilde_f(2) = (u(-2) - 5.0_precision*u(-1) + 13.0_precision*u(0) &
                      + 3.0_precision*u(1))/12.0_precision
        tilde_f(3) = (u(-1) - 7.0_precision*u(0) - 7.0_precision*u(1) + u(2))/(-12.0_precision)
        tilde_f(4) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1))/6.0_precision
        tilde_f(5) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2))/6.0_precision
        tilde_f(6) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0))/6.0_precision

        do loop = 0, 2
            gama(loop) = (1.0_precision + tau5/(beta(loop) + smv))**teno_q
        end do

        sum_gama = sum(gama)
        do loop = 0, 2
            chi(loop) = gama(loop)/sum_gama
        end do
        do loop = 0, 2
            delta_k(loop) = teno_cutoff(chi(loop))
        end do

        if (delta_k(0)*delta_k(1)*delta_k(2) == 1.0_precision) then
            flux_half = tilde_f(1)
        else if (delta_k(0)*delta_k(2) == 1.0_precision .and. delta_k(1) == 0.0_precision) then
            flux_half = tilde_f(2)
        else if (delta_k(0)*delta_k(1) == 1.0_precision .and. delta_k(2) == 0.0_precision) then
            flux_half = tilde_f(3)
        else if (delta_k(1)*delta_k(2) == 0.0_precision .and. delta_k(0) == 1.0_precision) then
            flux_half = tilde_f(4)
        else if (delta_k(0)*delta_k(2) == 0.0_precision .and. delta_k(1) == 1.0_precision) then
            flux_half = tilde_f(5)
        else if (delta_k(0)*delta_k(1) == 0.0_precision .and. delta_k(2) == 1.0_precision) then
            flux_half = tilde_f(6)
        else
            write(*, '(a)') 'something wrong!'
        end if
    end subroutine teno_mody_semi_eno

    ! Modified TENO semi-ENO with GA-RBF perturbation.
    subroutine teno_mody_semi_eno_rbf_ga(u, flux_half, h)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        real(precision) :: beta(0:2), tau5
        real(precision) :: gama(0:2), sum_gama
        real(precision) :: chi(0:2), delta_k(0:2)
        real(precision) :: rbf_append(1:6)
        real(precision) :: tilde_f(1:6)
        real(precision) :: eps, eps2h2
        real(precision), parameter :: smv = 10.0e-40_precision
        integer(integer_precision) :: loop

        beta(0) = 0.25_precision*((u(-1) - u(1))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-1) - 2.0_precision*u(0) + u(1))**2.0_precision)
        beta(1) = 0.25_precision*((3.0_precision*u(0) - 4.0_precision*u(1) + u(2))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(0) - 2.0_precision*u(1) + u(2))**2.0_precision)
        beta(2) = 0.25_precision*((u(-2) - 4.0_precision*u(-1) + 3.0_precision*u(0))**2.0_precision) &
                  + (13.0_precision/12.0_precision)*((u(-2) - 2.0_precision*u(-1) + u(0))**2.0_precision)

        tau5 = abs(beta(1) - beta(2))

        tilde_f(1) = (-2.0_precision*u(-2) + 13.0_precision*u(-1) - 47.0_precision*u(0) &
                      - 27.0_precision*u(1) + 3.0_precision*u(2))/(-60.0_precision)
        tilde_f(2) = (u(-2) - 5.0_precision*u(-1) + 13.0_precision*u(0) &
                      + 3.0_precision*u(1))/12.0_precision
        tilde_f(3) = (u(-1) - 7.0_precision*u(0) - 7.0_precision*u(1) + u(2))/(-12.0_precision)
        tilde_f(4) = (-u(-1) + 5.0_precision*u(0) + 2.0_precision*u(1))/6.0_precision
        tilde_f(5) = (2.0_precision*u(0) + 5.0_precision*u(1) - u(2))/6.0_precision
        tilde_f(6) = (2.0_precision*u(-2) - 7.0_precision*u(-1) + 11.0_precision*u(0))/6.0_precision

        do loop = 0, 2
            gama(loop) = (1.0_precision + tau5/(beta(loop) + smv))**teno_q
        end do

        sum_gama = sum(gama)
        do loop = 0, 2
            chi(loop) = gama(loop)/sum_gama
        end do
        do loop = 0, 2
            delta_k(loop) = teno_cutoff(chi(loop))
        end do

        call fulin_define_eps(eps, u(-2:3))
        eps2h2 = eps*eps*h*h

        rbf_append(1) = eps2h2*(18.0_precision*u(2) - 42.0_precision*u(1) &
                                + 18.0_precision*u(0) + 18.0_precision*u(-1) &
                                - 12.0_precision*u(-2))/(-60.0_precision)
        rbf_append(2) = eps2h2*(9.0_precision*u(1) - 15.0_precision*u(0) &
                                + 3.0_precision*u(-1) + 3.0_precision*u(-2))/12.0_precision
        rbf_append(3) = eps2h2*(4.0_precision*u(-1) - 4.0_precision*u(0) &
                                - 4.0_precision*u(1) + 4.0_precision*u(2))/(-12.0_precision)
        rbf_append(4) = eps2h2*(-2.0_precision*u(-1) - 2.0_precision*u(0) &
                                + 4.0_precision*u(1))/6.0_precision
        rbf_append(5) = eps2h2*(-4.0_precision*u(0) + 2.0_precision*u(1) &
                                + 2.0_precision*u(2))/(-6.0_precision)
        rbf_append(6) = eps2h2*(-18.0_precision*u(-1) + 18.0_precision*u(0))/(-6.0_precision)

        if (delta_k(0)*delta_k(1)*delta_k(2) == 1.0_precision) then
            flux_half = tilde_f(1) + rbf_append(1)
        else if (delta_k(0)*delta_k(2) == 1.0_precision .and. delta_k(1) == 0.0_precision) then
            flux_half = tilde_f(2) + rbf_append(2)
        else if (delta_k(0)*delta_k(1) == 1.0_precision .and. delta_k(2) == 0.0_precision) then
            flux_half = tilde_f(3) + rbf_append(3)
        else if (delta_k(1)*delta_k(2) == 0.0_precision .and. delta_k(0) == 1.0_precision) then
            flux_half = tilde_f(4) + rbf_append(4)
        else if (delta_k(0)*delta_k(2) == 0.0_precision .and. delta_k(1) == 1.0_precision) then
            flux_half = tilde_f(5) + rbf_append(5)
        else if (delta_k(0)*delta_k(1) == 0.0_precision .and. delta_k(2) == 1.0_precision) then
            flux_half = tilde_f(6) + rbf_append(6)
        else
            write(*, '(a)') 'something wrong!'
        end if
    end subroutine teno_mody_semi_eno_rbf_ga

    ! THINC interface position.
    function get_x_hat(minus, mid, plus, beta) result(hat_out)
        real(precision), intent(in) :: minus, mid, plus, beta
        real(precision) :: hat_out
        real(precision) :: q_min, q_max, theta, hat

        q_min = min(minus, plus)
        q_max = max(minus, plus) - q_min
        theta = sign(plus - minus)

        hat = 2.0_precision*(mid - q_min - 0.5_precision*q_max)*beta/(q_max*theta)

        hat_out = -0.5_precision*(1.0_precision + log(exp(hat - 2.0_precision*beta) + exp(hat)))
    end function get_x_hat

    ! THINC reconstruction (unvalidated: original expression is suspect).
    subroutine get_thinc(u, flux_half)
        real(precision), intent(in) :: u(-2:3)
        real(precision), intent(out) :: flux_half

        real(precision) :: x_hat

        x_hat = get_x_hat(u(-1), u(0), u(1), thinc_beta)

        flux_half = min(u(-1), u(1)) + ((max(u(-1), u(1)) - min(u(-1), u(1)))/2.0_precision) &
                    * (1.0_precision + sign(u(1)) - u(-1))*tanh(1.0_precision*(1.0_precision - x_hat))
    end subroutine get_thinc

    ! Sign function returning -1, 0, or +1.
    function sign(a) result(s)
        real(precision), intent(in) :: a
        real(precision) :: s

        if (a == 0.0_precision) then
            s = 0.0_precision
        else if (a > 0.0_precision) then
            s = 1.0_precision
        else
            s = -1.0_precision
        end if
    end function sign

end module weno_schemes
