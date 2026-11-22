module e1d_flux
    ! Thermodynamics, characteristic projection, Steger-Warming flux split,
    ! and characteristic interpolation for the 1D Euler equations.
    !
    ! The projection is the 3x3 reduction of the 4x4 matrices used by
    ! accuracy/src/euler2d/e2d_flux.f90 after removing the transverse
    ! velocity component.
    use precision_select
    use e1d_globals, only: halo, gama, beta
    use scheme_dispatch, only: interpolate_dispatch
    implicit none

    real(precision), parameter :: sw_eps2 = 0.12_precision**2
    real(precision), parameter :: sw_sig = 0.5_precision/gama

contains

    pure function e1d_pressure(u) result(p)
        ! Ideal-gas pressure from the conservative state (rho, rho*u, rho*E).
        real(precision), intent(in) :: u(3)
        real(precision) :: p

        p = beta*(u(3) - 0.5_precision*u(2)*u(2)/u(1))
    end function e1d_pressure

    pure function e1d_sound_speed(u) result(c)
        ! Physical sound speed c = sqrt(gama*p/rho).
        real(precision), intent(in) :: u(3)
        real(precision) :: c

        c = sqrt(gama*e1d_pressure(u)/u(1))
    end function e1d_sound_speed

    subroutine e1d_physical_flux(u, f)
        ! Physical Euler flux E(U).
        real(precision), intent(in) :: u(3)
        real(precision), intent(out) :: f(3)

        real(precision) :: rho, vel, p

        rho = u(1)
        vel = u(2)/rho
        p = e1d_pressure(u)

        f(1) = u(2)
        f(2) = u(2)*vel + p
        f(3) = (u(3) + p)*vel
    end subroutine e1d_physical_flux

    subroutine e1d_eigen_matrix(u, ml, mr)
        ! Left and right characteristic matrices for the 3x3 system.
        real(precision), intent(in) :: u(3)
        real(precision), intent(out) :: ml(3, 3), mr(3, 3)

        real(precision) :: rho, vel, p, c, c2, alfa, scale

        rho = u(1)
        vel = u(2)/rho
        p = e1d_pressure(u)
        c2 = gama*p/rho
        c = sqrt(c2)
        alfa = 0.5_precision*vel*vel

        ml(1, 1) = alfa*beta - c2
        ml(1, 2) = -beta*vel
        ml(1, 3) = beta
        ml(2, 1) = alfa*beta - c*vel
        ml(2, 2) = c - beta*vel
        ml(2, 3) = beta
        ml(3, 1) = alfa*beta + c*vel
        ml(3, 2) = -c - beta*vel
        ml(3, 3) = beta

        scale = 0.5_precision/c2
        mr(:, 1) = (/-2.0_precision, -2.0_precision*vel, -2.0_precision*alfa/)
        mr(:, 2) = (/1.0_precision, vel + c, alfa + c*vel + c2/beta/)
        mr(:, 3) = (/1.0_precision, vel - c, alfa - c*vel + c2/beta/)
        mr = mr*scale
    end subroutine e1d_eigen_matrix

    subroutine e1d_sw(ul, ur, fhalf)
        ! Modified Steger-Warming split at a half node.
        !
        ! fp uses the positive part of the left state and fn the negative
        ! part of the right state; the physical density multiplies each part.
        real(precision), intent(in) :: ul(3), ur(3)
        real(precision), intent(out) :: fhalf(3)

        real(precision) :: fp(3), fn(3)

        call e1d_sw_state(ul, 1.0_precision, fp)
        call e1d_sw_state(ur, -1.0_precision, fn)
        fhalf = ul(1)*fp + ur(1)*fn
    end subroutine e1d_sw

    subroutine e1d_sw_state(u, sgn, f)
        ! Steger-Warming split of one state; sgn = +1 for F+, -1 for F-.
        real(precision), intent(in) :: u(3), sgn
        real(precision), intent(out) :: f(3)

        real(precision) :: rho, vel, p, c2, c, alfa
        real(precision) :: lam1, lam4, lam5, diss, disp

        rho = u(1)
        vel = u(2)/rho
        p = e1d_pressure(u)
        c2 = gama*p/rho
        c = sqrt(c2)
        alfa = 0.5_precision*vel*vel

        lam1 = vel
        lam4 = vel + c
        lam5 = vel - c
        lam1 = (lam1 + sgn*sqrt(lam1*lam1 + sw_eps2))/2.0_precision
        lam4 = (lam4 + sgn*sqrt(lam4*lam4 + sw_eps2))/2.0_precision
        lam5 = (lam5 + sgn*sqrt(lam5*lam5 + sw_eps2))/2.0_precision

        diss = -sw_sig*(2.0_precision*lam1 - lam4 - lam5)
        disp = sw_sig*(lam4 - lam5)*c

        f(1) = diss + lam1
        f(2) = vel*diss + vel*lam1 + disp
        f(3) = (alfa + c2/beta)*diss + (alfa + p/(rho*beta))*lam1 + vel*disp
    end subroutine e1d_sw_state

    subroutine e1d_interpolate_character(ustencil, uhalf, correction)
        ! Interpolate one 7-point stencil and apply the positivity guard.
        real(precision), intent(in) :: ustencil(-3:4, 3)
        real(precision), intent(out) :: uhalf(3)
        integer, intent(out) :: correction

        call e1d_reconstruct_character(ustencil, uhalf)
        call e1d_apply_positivity(ustencil, uhalf, correction)
    end subroutine e1d_interpolate_character

    subroutine e1d_reconstruct_character(ustencil, uhalf)
        ! Characteristic interpolation without the positivity guard.
        real(precision), intent(in) :: ustencil(-3:4, 3)
        real(precision), intent(out) :: uhalf(3)

        real(precision) :: u(3), ml(3, 3), mr(3, 3)
        real(precision) :: q(-3:4, 3), qhalf(3)
        integer :: p

        u = ustencil(0, :)
        call e1d_eigen_matrix(u, ml, mr)

        do p = 1, 3
            q(:, p) = ml(p, 1)*ustencil(:, 1) + ml(p, 2)*ustencil(:, 2) &
                    + ml(p, 3)*ustencil(:, 3)
            call interpolate_dispatch(q(:, p), qhalf(p))
        end do

        uhalf = mr(:, 1)*qhalf(1) + mr(:, 2)*qhalf(2) + mr(:, 3)*qhalf(3)
    end subroutine e1d_reconstruct_character

    subroutine e1d_apply_positivity(ustencil, uhalf, correction)
        ! Return the centre state if the reconstruction is below 0.1% of the
        ! stencil minimum density or pressure.
        real(precision), intent(in) :: ustencil(-3:4, 3)
        real(precision), intent(inout) :: uhalf(3)
        integer, intent(out) :: correction

        real(precision) :: rho_min, p_min, rho_h, p_h
        integer :: j

        rho_min = minval(ustencil(-2:2, 1))
        p_min = huge(1.0_precision)
        do j = -2, 2
            p_min = min(p_min, e1d_pressure(ustencil(j, :)))
        end do

        correction = 0
        rho_h = uhalf(1)
        if (rho_h <= 0.001_precision*rho_min) then
            correction = 1
        else
            p_h = e1d_pressure(uhalf)
            if (p_h <= 0.001_precision*p_min) correction = 1
        end if

        if (correction /= 0) uhalf = ustencil(0, :)
    end subroutine e1d_apply_positivity

    subroutine e1d_interflux(n, u, eflux, ncorrection)
        ! Assemble all half-node fluxes on j = -3 .. n+2.
        !
        ! ncorrection optionally returns how many reconstructions used the
        ! positivity fallback.
        integer, intent(in) :: n
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        real(precision), intent(out) :: eflux(-3:n + 2, 3)
        integer, intent(out), optional :: ncorrection

        real(precision) :: ul(3), ur(3)
        integer :: j, correction, ncorr

        ncorr = 0
        do j = -3, n + 2
            call e1d_interpolate_character(u(j - 3:j + 4, :), ul, correction)
            ncorr = ncorr + correction
            call e1d_interpolate_character(u(j + 4:j - 3:-1, :), ur, correction)
            ncorr = ncorr + correction
            call e1d_sw(ul, ur, eflux(j, :))
        end do

        if (present(ncorrection)) ncorrection = ncorr
    end subroutine e1d_interflux

end module e1d_flux
