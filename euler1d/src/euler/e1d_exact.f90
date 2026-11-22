module e1d_exact
    ! Exact Riemann solution for the 1D Euler equations with gamma = 1.4.
    !
    ! The solver supports the four wave patterns (two shocks, two
    ! rarefactions, left shock/right rarefaction, left rarefaction/right
    ! shock) and the vacuum case.  The initial discontinuity sits at x = 0.5,
    ! so every sample uses xi = (x - 0.5)/t.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_select
    use e1d_globals, only: gama
    use e1d_initial, only: e1d_case_states
    implicit none

contains

    function e1d_f_pstar(p, p_k, rho_k) result(f)
        ! Pressure function f_K(p) used by the Riemann solver.
        real(precision), intent(in) :: p, p_k, rho_k
        real(precision) :: f

        real(precision) :: c_k

        c_k = sqrt(gama*p_k/rho_k)
        if (p > p_k) then
            f = (p - p_k)/(rho_k*c_k*sqrt((gama + 1.0_precision)/(2.0_precision*gama) &
                *p/p_k + (gama - 1.0_precision)/(2.0_precision*gama)))
        else
            f = 2.0_precision*c_k/(gama - 1.0_precision) &
                *((p/p_k)**((gama - 1.0_precision)/(2.0_precision*gama)) - 1.0_precision)
        end if
    end function e1d_f_pstar

    function e1d_f_pstar_deriv(p, p_k, rho_k) result(df)
        ! Derivative of f_K(p).
        real(precision), intent(in) :: p, p_k, rho_k
        real(precision) :: df

        real(precision) :: c_k, aa, bb, dd, ee

        c_k = sqrt(gama*p_k/rho_k)
        if (p > p_k) then
            aa = (gama + 1.0_precision)/(2.0_precision*gama)
            bb = (gama - 1.0_precision)/(2.0_precision*gama)
            dd = aa*p/p_k + bb
            df = (dd - (p - p_k)*aa/(2.0_precision*p_k)) &
                 /(rho_k*c_k*dd**1.5_precision)
        else
            ee = (gama - 1.0_precision)/(2.0_precision*gama)
            if (p <= 0.0_precision) then
                df = huge(1.0_precision)
            else
                df = 2.0_precision*c_k/(gama - 1.0_precision)*ee/p_k &
                     *(p/p_k)**(ee - 1.0_precision)
            end if
        end if
    end function e1d_f_pstar_deriv

    subroutine e1d_riemann_star(rho_l, u_l, p_l, rho_r, u_r, p_r, &
                                p_star, u_star, rho_star_l, rho_star_r, ierr, msg, &
                                max_iter)
        ! Star pressure, velocity, and densities from the Riemann problem.
        !
        ! max_iter is an optional iteration cap for tests; the default is 200.
        real(precision), intent(in) :: rho_l, u_l, p_l, rho_r, u_r, p_r
        real(precision), intent(out) :: p_star, u_star, rho_star_l, rho_star_r
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg
        integer, intent(in), optional :: max_iter

        real(precision) :: c_l, c_r, du, p_lo, p_hi, p_new, fsum, dfsum, tol
        real(precision) :: f_l, f_r
        integer :: iter, maxit
        logical :: bracketed, converged

        ierr = 0
        msg = ''
        maxit = 200
        if (present(max_iter)) maxit = max_iter

        c_l = sqrt(gama*p_l/rho_l)
        c_r = sqrt(gama*p_r/rho_r)
        du = u_l - u_r

        ! Two rarefactions that leave a vacuum in the centre.
        if (u_r - u_l >= 2.0_precision*c_l/(gama - 1.0_precision) &
            + 2.0_precision*c_r/(gama - 1.0_precision)) then
            p_star = 0.0_precision
            u_star = 0.5_precision*(u_l + u_r)
            rho_star_l = 0.0_precision
            rho_star_r = 0.0_precision
            return
        end if

        ! Bracket the root in [0, p_hi].
        bracketed = .false.
        p_lo = 0.0_precision
        p_hi = max(p_l, p_r, 1.0_precision)
        do iter = 1, 100
            fsum = e1d_f_pstar(p_hi, p_l, rho_l) + e1d_f_pstar(p_hi, p_r, rho_r) - du
            if (fsum >= 0.0_precision) then
                bracketed = .true.
                exit
            end if
            p_hi = 2.0_precision*p_hi
        end do

        if (.not. bracketed) then
            ierr = 2
            msg = 'Riemann solver failed to bracket the star pressure'
            return
        end if

        ! Newton iteration with a bisection fallback.
        converged = .false.
        p_star = max(min(p_l, p_r), 1.0e-12_precision)
        do iter = 1, maxit
            tol = 100.0_precision*epsilon(1.0_precision)*max(1.0_precision, p_star)
            fsum = e1d_f_pstar(p_star, p_l, rho_l) &
                 + e1d_f_pstar(p_star, p_r, rho_r) - du
            dfsum = e1d_f_pstar_deriv(p_star, p_l, rho_l) &
                  + e1d_f_pstar_deriv(p_star, p_r, rho_r)

            if (dfsum <= 0.0_precision .or. .not. ieee_is_finite(dfsum)) then
                p_new = 0.5_precision*(p_lo + p_hi)
            else
                p_new = p_star - fsum/dfsum
            end if

            if (p_new <= p_lo .or. p_new >= p_hi .or. .not. ieee_is_finite(p_new)) then
                p_new = 0.5_precision*(p_lo + p_hi)
            end if

            if (fsum > 0.0_precision) then
                p_hi = p_star
            else
                p_lo = p_star
            end if

            if (abs(p_new - p_star) <= tol) then
                p_star = p_new
                converged = .true.
                exit
            end if
            p_star = p_new
        end do

        if (.not. converged) then
            ierr = 1
            msg = 'Riemann solver did not converge'
            return
        end if

        f_l = e1d_f_pstar(p_star, p_l, rho_l)
        f_r = e1d_f_pstar(p_star, p_r, rho_r)
        u_star = 0.5_precision*(u_l + u_r - f_l + f_r)

        if (p_star > p_l) then
            rho_star_l = rho_l*((gama + 1.0_precision)*p_star/p_l + gama - 1.0_precision) &
                       /((gama - 1.0_precision)*p_star/p_l + gama + 1.0_precision)
        else
            rho_star_l = rho_l*(p_star/p_l)**(1.0_precision/gama)
        end if

        if (p_star > p_r) then
            rho_star_r = rho_r*((gama + 1.0_precision)*p_star/p_r + gama - 1.0_precision) &
                       /((gama - 1.0_precision)*p_star/p_r + gama + 1.0_precision)
        else
            rho_star_r = rho_r*(p_star/p_r)**(1.0_precision/gama)
        end if
    end subroutine e1d_riemann_star

    subroutine e1d_riemann_sample(rho_l, u_l, p_l, rho_r, u_r, p_r, t, x, rho, u, p)
        ! Sample the exact solution at (x, t).
        real(precision), intent(in) :: rho_l, u_l, p_l, rho_r, u_r, p_r, t, x
        real(precision), intent(out) :: rho, u, p

        real(precision) :: xi, c_l, c_r, ps, us, rsl, rsr
        real(precision) :: csl, csr, al, ar
        real(precision) :: zl_head, zl_tail, zl_shock
        real(precision) :: zr_head, zr_tail, zr_shock
        integer :: ierr
        character(len=256) :: msg

        if (t <= 0.0_precision) then
            if (x < 0.5_precision) then
                rho = rho_l
                u = u_l
                p = p_l
            else
                rho = rho_r
                u = u_r
                p = p_r
            end if
            return
        end if

        c_l = sqrt(gama*p_l/rho_l)
        c_r = sqrt(gama*p_r/rho_r)
        xi = (x - 0.5_precision)/t

        ! Vacuum: two rarefaction fans separated by an empty region.
        if (u_r - u_l >= 2.0_precision*c_l/(gama - 1.0_precision) &
            + 2.0_precision*c_r/(gama - 1.0_precision)) then
            zl_head = u_l - c_l
            zl_tail = u_l + 2.0_precision*c_l/(gama - 1.0_precision)
            zr_head = u_r + c_r
            zr_tail = u_r - 2.0_precision*c_r/(gama - 1.0_precision)
            if (xi < zl_head) then
                rho = rho_l
                u = u_l
                p = p_l
            else if (xi < zl_tail) then
                call left_fan(u_l, p_l, xi, c_l, rho, u, p)
            else if (xi < zr_tail) then
                rho = 0.0_precision
                u = 0.0_precision
                p = 0.0_precision
            else if (xi < zr_head) then
                call right_fan(u_r, p_r, xi, c_r, rho, u, p)
            else
                rho = rho_r
                u = u_r
                p = p_r
            end if
            return
        end if

        call e1d_riemann_star(rho_l, u_l, p_l, rho_r, u_r, p_r, &
                              ps, us, rsl, rsr, ierr, msg)
        if (ierr /= 0) then
            write (*, '(a)') 'error: ' // trim(msg)
            error stop 1
        end if

        csl = c_l*(ps/p_l)**((gama - 1.0_precision)/(2.0_precision*gama))
        csr = c_r*(ps/p_r)**((gama - 1.0_precision)/(2.0_precision*gama))
        al = sqrt((gama + 1.0_precision)/(2.0_precision*gama)*ps/p_l &
                  + (gama - 1.0_precision)/(2.0_precision*gama))
        ar = sqrt((gama + 1.0_precision)/(2.0_precision*gama)*ps/p_r &
                  + (gama - 1.0_precision)/(2.0_precision*gama))

        zl_head = u_l - c_l
        zl_tail = us - csl
        zl_shock = u_l - al*c_l
        zr_head = u_r + c_r
        zr_tail = us + csr
        zr_shock = u_r + ar*c_r

        if (ps > p_l) then
            if (ps > p_r) then
                ! Two shocks.
                if (xi < zl_shock) then
                    call left_state(rho_l, u_l, p_l, rho, u, p)
                else if (xi < us) then
                    call star_state(rsl, us, ps, rho, u, p)
                else if (xi < zr_shock) then
                    call star_state(rsr, us, ps, rho, u, p)
                else
                    call right_state(rho_r, u_r, p_r, rho, u, p)
                end if
            else
                ! Left shock, right rarefaction.
                if (xi < zl_shock) then
                    call left_state(rho_l, u_l, p_l, rho, u, p)
                else if (xi < us) then
                    call star_state(rsl, us, ps, rho, u, p)
                else if (xi < zr_tail) then
                    call star_state(rsr, us, ps, rho, u, p)
                else if (xi < zr_head) then
                    call right_fan(u_r, p_r, xi, c_r, rho, u, p)
                else
                    call right_state(rho_r, u_r, p_r, rho, u, p)
                end if
            end if
        else
            if (ps > p_r) then
                ! Left rarefaction, right shock.
                if (xi < zl_head) then
                    call left_state(rho_l, u_l, p_l, rho, u, p)
                else if (xi < zl_tail) then
                    call left_fan(u_l, p_l, xi, c_l, rho, u, p)
                else if (xi < us) then
                    call star_state(rsl, us, ps, rho, u, p)
                else if (xi < zr_shock) then
                    call star_state(rsr, us, ps, rho, u, p)
                else
                    call right_state(rho_r, u_r, p_r, rho, u, p)
                end if
            else
                ! Two rarefactions, non-vacuum.
                if (xi < zl_head) then
                    call left_state(rho_l, u_l, p_l, rho, u, p)
                else if (xi < zl_tail) then
                    call left_fan(u_l, p_l, xi, c_l, rho, u, p)
                else if (xi < us) then
                    call star_state(rsl, us, ps, rho, u, p)
                else if (xi < zr_tail) then
                    call star_state(rsr, us, ps, rho, u, p)
                else if (xi < zr_head) then
                    call right_fan(u_r, p_r, xi, c_r, rho, u, p)
                else
                    call right_state(rho_r, u_r, p_r, rho, u, p)
                end if
            end if
        end if
    end subroutine e1d_riemann_sample

    subroutine left_fan(u_l, p_l, xi, c_l, rho, u, p)
        ! Left rarefaction fan.
        real(precision), intent(in) :: u_l, p_l, xi, c_l
        real(precision), intent(out) :: rho, u, p

        real(precision) :: c_a

        c_a = (gama - 1.0_precision)/(gama + 1.0_precision)*(u_l - xi) &
            + 2.0_precision/(gama + 1.0_precision)*c_l
        u = xi + c_a
        p = p_l*(c_a/c_l)**(2.0_precision*gama/(gama - 1.0_precision))
        rho = gama*p/(c_a*c_a)
    end subroutine left_fan

    subroutine right_fan(u_r, p_r, xi, c_r, rho, u, p)
        ! Right rarefaction fan.
        real(precision), intent(in) :: u_r, p_r, xi, c_r
        real(precision), intent(out) :: rho, u, p

        real(precision) :: c_a

        c_a = (gama - 1.0_precision)/(gama + 1.0_precision)*(xi - u_r) &
            + 2.0_precision/(gama + 1.0_precision)*c_r
        u = xi - c_a
        p = p_r*(c_a/c_r)**(2.0_precision*gama/(gama - 1.0_precision))
        rho = gama*p/(c_a*c_a)
    end subroutine right_fan

    subroutine left_state(rho_l, u_l, p_l, rho, u, p)
        real(precision), intent(in) :: rho_l, u_l, p_l
        real(precision), intent(out) :: rho, u, p

        rho = rho_l
        u = u_l
        p = p_l
    end subroutine left_state

    subroutine right_state(rho_r, u_r, p_r, rho, u, p)
        real(precision), intent(in) :: rho_r, u_r, p_r
        real(precision), intent(out) :: rho, u, p

        rho = rho_r
        u = u_r
        p = p_r
    end subroutine right_state

    subroutine star_state(rho_s, u_s, p_s, rho, u, p)
        real(precision), intent(in) :: rho_s, u_s, p_s
        real(precision), intent(out) :: rho, u, p

        rho = rho_s
        u = u_s
        p = p_s
    end subroutine star_state

    subroutine e1d_exact_profile(n, case_name, t, rho, u, p)
        ! Exact profile on x_j = j/n, j = 0..n.
        integer, intent(in) :: n
        character(len=*), intent(in) :: case_name
        real(precision), intent(in) :: t
        real(precision), intent(out) :: rho(0:n), u(0:n), p(0:n)

        real(precision) :: rho_l, u_l, p_l, rho_r, u_r, p_r, x
        integer :: j

        call e1d_case_states(case_name, rho_l, u_l, p_l, rho_r, u_r, p_r)
        do j = 0, n
            x = real(j, precision)/real(n, precision)
            call e1d_riemann_sample(rho_l, u_l, p_l, rho_r, u_r, p_r, t, x, &
                                    rho(j), u(j), p(j))
        end do
    end subroutine e1d_exact_profile

end module e1d_exact
