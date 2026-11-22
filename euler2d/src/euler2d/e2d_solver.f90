module e2d_solver
    ! Spatial residual, CFL time step, positivity correction, and fixed
    ! four-stage Runge-Kutta for the 2D Euler test cases.
    !
    ! All state lives in the e2d_globals working field ubar; e2d_rhs fills
    ! the half-node fluxes and returns -dE/dx - dG/dz on the interior nodes.
    use precision_select
    use e2d_globals, only: mt, n, m, h, gama, beta, ubar, &
                           eint, gint, eflux, gflux
    use e2d_flux, only: u_to_prim_e_g, interflux
    use e2d_cases, only: e2d_boundary, e2d_source
    use difference_operators, only: difference_apply
    implicit none

contains

    subroutine e2d_rhs(case_name, op, rhs, ncorr)
        ! Interior residual rhs = -dE/dx - dG/dz for the selected operator.
        !
        ! difference_apply already returns the negative of the derivative,
        ! so the two directional residuals are simply added.  ncorr returns
        ! how many interpolation positivity fallbacks fired.
        character(len=*), intent(in) :: case_name
        character(len=*), intent(in) :: op
        real(precision), intent(out) :: rhs(-mt:n + mt, -mt:m + mt, 4)
        integer, intent(out) :: ncorr

        real(precision) :: kx, ky
        integer :: j, k, p

        call u_to_prim_e_g()
        call interflux(ncorr)

        rhs = 0.0_precision
        do p = 1, 4
            do k = 0, m
                do j = 0, n
                    call difference_apply(op, h, eflux(j - 3:j + 2, k, p), &
                                          eint(j - 2:j + 2, k, p), kx)
                    call difference_apply(op, h, gflux(j, k - 3:k + 2, p), &
                                          gint(j, k - 2:k + 2, p), ky)
                    rhs(j, k, p) = kx + ky
                end do
            end do
        end do

        call e2d_source(case_name, rhs)
    end subroutine e2d_rhs

    function e2d_timestep(cfl) result(dt)
        ! CFL time step: dt = cfl*h/max(sqrt(u^2+w^2) + sqrt(gama*p/rho)).
        ! The second term is the sound speed used by the flux split.
        real(precision), intent(in) :: cfl
        real(precision) :: dt

        real(precision) :: ro, uu, ww, ps, temp, maxtemp
        integer :: j, k

        maxtemp = -huge(1.0_precision)
        do k = 0, m
            do j = 0, n
                ro = ubar(j, k, 1)
                uu = ubar(j, k, 2)/ro
                ww = ubar(j, k, 3)/ro
                ps = beta*(ubar(j, k, 4) - 0.5_precision*(uu*uu + ww*ww)*ro)
                temp = sqrt(uu*uu + ww*ww) + sqrt(gama*ps/ro)
                if (temp > maxtemp) maxtemp = temp
            end do
        end do

        if (maxtemp <= 0.0_precision) error stop 'e2d_timestep: non-positive speed'
        dt = cfl*h/maxtemp
    end function e2d_timestep

    subroutine e2d_apply_correction(u0, ncorr)
        ! Replace non-physical interior cells by the average of their four
        ! neighbours taken from the saved state u0.
        !
        ! The density test comes first so that a zero or negative density is
        ! never used as a divisor when the pressure is evaluated.  ncorr
        ! returns the number of corrected cells.
        real(precision), intent(in) :: u0(-mt:n + mt, -mt:m + mt, 4)
        integer, intent(out) :: ncorr

        real(precision) :: ro, uu, ww, ps
        integer :: j, k

        ncorr = 0
        do k = 0, m
            do j = 0, n
                ro = ubar(j, k, 1)
                if (ro < 1.0e-8_precision) then
                    ubar(j, k, :) = (u0(j - 1, k, :) + u0(j + 1, k, :) &
                                   + u0(j, k - 1, :) + u0(j, k + 1, :))/4.0_precision
                    ncorr = ncorr + 1
                else
                    uu = ubar(j, k, 2)/ro
                    ww = ubar(j, k, 3)/ro
                    ps = beta*(ubar(j, k, 4) - 0.5_precision*(uu*uu + ww*ww)*ro)
                    if (ps < 1.0e-8_precision) then
                        ubar(j, k, :) = (u0(j - 1, k, :) + u0(j + 1, k, :) &
                                       + u0(j, k - 1, :) + u0(j, k + 1, :))/4.0_precision
                        ncorr = ncorr + 1
                    end if
                end if
            end do
        end do
    end subroutine e2d_apply_correction

    subroutine e2d_advance(case_name, op, cfl, t_start, t_end, steps, &
                           interp_corr, cell_corr, cpu_seconds)
        ! Advance ubar from t_start to t_end with fixed four-stage RK4.
        !
        ! The halo is refreshed before every residual evaluation and the
        ! positivity correction is applied after every RK stage.  interp_corr
        ! and cell_corr accumulate the interpolation and cell corrections,
        ! steps counts the accepted time steps, and cpu_seconds reports the
        ! wall-clock time spent in the loop.
        character(len=*), intent(in) :: case_name
        character(len=*), intent(in) :: op
        real(precision), intent(in) :: cfl, t_start, t_end
        integer, intent(out) :: steps, interp_corr, cell_corr
        real(precision), intent(out) :: cpu_seconds

        real(precision), allocatable :: u0(:, :, :), k1(:, :, :)
        real(precision), allocatable :: k2(:, :, :), k3(:, :, :), k4(:, :, :)
        real(precision) :: dt, time, t0, t1
        integer :: nci, nc
        logical :: last

        allocate (u0(-mt:n + mt, -mt:m + mt, 4))
        allocate (k1(-mt:n + mt, -mt:m + mt, 4))
        allocate (k2(-mt:n + mt, -mt:m + mt, 4))
        allocate (k3(-mt:n + mt, -mt:m + mt, 4))
        allocate (k4(-mt:n + mt, -mt:m + mt, 4))

        call e2d_boundary(case_name, t_start)
        u0 = ubar
        steps = 0
        interp_corr = 0
        cell_corr = 0
        time = t_start
        call cpu_time(t0)

        do
            dt = e2d_timestep(cfl)
            last = .false.
            if (time + dt >= t_end) then
                dt = t_end - time
                last = .true.
            end if
            time = time + dt
            steps = steps + 1

            ! Stage 1.
            call e2d_boundary(case_name, time)
            call e2d_rhs(case_name, op, k1, nci)
            interp_corr = interp_corr + nci
            ubar = u0 + 0.5_precision*dt*k1
            call e2d_apply_correction(u0, nc)
            cell_corr = cell_corr + nc

            ! Stage 2.
            call e2d_boundary(case_name, time)
            call e2d_rhs(case_name, op, k2, nci)
            interp_corr = interp_corr + nci
            ubar = u0 + 0.5_precision*dt*k2
            call e2d_apply_correction(u0, nc)
            cell_corr = cell_corr + nc

            ! Stage 3.
            call e2d_boundary(case_name, time)
            call e2d_rhs(case_name, op, k3, nci)
            interp_corr = interp_corr + nci
            ubar = u0 + dt*k3
            call e2d_apply_correction(u0, nc)
            cell_corr = cell_corr + nc

            ! Stage 4.
            call e2d_boundary(case_name, time)
            call e2d_rhs(case_name, op, k4, nci)
            interp_corr = interp_corr + nci
            ubar = u0 + (k1 + 2.0_precision*k2 + 2.0_precision*k3 + k4) &
                  *dt/6.0_precision
            call e2d_apply_correction(u0, nc)
            cell_corr = cell_corr + nc

            call e2d_boundary(case_name, time)
            u0 = ubar
            if (last) exit
        end do

        call cpu_time(t1)
        cpu_seconds = t1 - t0

        deallocate (u0, k1, k2, k3, k4)
    end subroutine e2d_advance

end module e2d_solver
