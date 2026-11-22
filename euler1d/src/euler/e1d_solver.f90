module e1d_solver
    ! Boundary conditions, positivity correction, CFL time step, residual,
    ! and four-stage Runge-Kutta for the 1D Euler test cases.
    use precision_select
    use e1d_globals, only: halo
    use e1d_flux, only: e1d_interflux, e1d_physical_flux, e1d_pressure
    use difference_operators, only: difference_apply
    implicit none

contains

    subroutine e1d_boundary(n, u)
        ! Zero-order extrapolation into the halo.
        integer, intent(in) :: n
        real(precision), intent(inout) :: u(-halo:n + halo, 3)
        integer :: j

        do j = 1, halo
            u(-j, :) = u(j - 1, :)
            u(n + j, :) = u(n - j + 1, :)
        end do
    end subroutine e1d_boundary

    subroutine e1d_correction(n, u_old, u_new, ncorrect)
        ! Replace non-physical cells by the average of their old neighbours.
        !
        ! The density test comes first so that a zero or negative density is
        ! never used as a divisor when the pressure is evaluated.
        integer, intent(in) :: n
        real(precision), intent(in) :: u_old(-halo:n + halo, 3)
        real(precision), intent(inout) :: u_new(-halo:n + halo, 3)
        integer, intent(out) :: ncorrect

        real(precision) :: rho, p
        integer :: j

        ncorrect = 0
        do j = 0, n
            rho = u_new(j, 1)
            if (rho < 1.0e-8_precision) then
                u_new(j, :) = 0.5_precision*(u_old(j - 1, :) + u_old(j + 1, :))
                ncorrect = ncorrect + 1
            else
                p = e1d_pressure(u_new(j, :))
                if (p < 1.0e-8_precision) then
                    u_new(j, :) = 0.5_precision*(u_old(j - 1, :) + u_old(j + 1, :))
                    ncorrect = ncorrect + 1
                end if
            end if
        end do
    end subroutine e1d_correction

    function e1d_timestep(n, h, u, cfl) result(dt)
        ! CFL time step: dt = cfl*h/max(|u| + sqrt(p/rho)).
        integer, intent(in) :: n
        real(precision), intent(in) :: h, cfl
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        real(precision) :: dt

        real(precision) :: rho, vel, p, temp, maxtemp
        integer :: j

        maxtemp = -huge(1.0_precision)
        do j = 0, n
            rho = u(j, 1)
            p = e1d_pressure(u(j, :))
            if (rho <= 0.0_precision .or. p <= 0.0_precision) then
                write (*, '(a,i0)') 'e1d_timestep: non-physical state at j = ', j
                error stop 1
            end if
            vel = u(j, 2)/rho
            temp = abs(vel) + sqrt(p/rho)
            if (temp > maxtemp) maxtemp = temp
        end do

        if (maxtemp <= 0.0_precision) error stop 'e1d_timestep: non-positive speed'
        dt = cfl*h/maxtemp
    end function e1d_timestep

    subroutine e1d_rhs(n, h, u, op, rhs, ncorrection)
        ! Spatial right-hand side rhs = -dE/dx with the selected operator.
        !
        ! ncorrection optionally returns the number of interpolation
        ! positivity fallbacks used by this residual.
        integer, intent(in) :: n
        real(precision), intent(in) :: h
        real(precision), intent(in) :: u(-halo:n + halo, 3)
        character(len=*), intent(in) :: op
        real(precision), intent(out) :: rhs(-halo:n + halo, 3)
        integer, intent(out), optional :: ncorrection

        real(precision), allocatable :: eint(:, :), eflux(:, :)
        integer :: j, p, ncorr

        allocate (eint(-halo:n + halo, 3))
        allocate (eflux(-3:n + 2, 3))

        do j = -halo, n + halo
            call e1d_physical_flux(u(j, :), eint(j, :))
        end do
        call e1d_interflux(n, u, eflux, ncorr)

        rhs = 0.0_precision
        do p = 1, 3
            do j = 0, n
                call difference_apply(op, h, eflux(j - 3:j + 2, p), &
                                      eint(j - 2:j + 2, p), rhs(j, p))
            end do
        end do

        if (present(ncorrection)) ncorrection = ncorr

        deallocate (eint, eflux)
    end subroutine e1d_rhs

    subroutine e1d_solve(n, h, cfl, tend, op, u, steps, ncorrect, t_final, cpu, &
                         ncorr_interp)
        ! Advance the state to tend with four-stage RK4.
        !
        ! ncorr_interp optionally returns how many interpolation positivity
        ! fallbacks fired over all stages.
        integer, intent(in) :: n
        real(precision), intent(in) :: h, cfl, tend
        character(len=*), intent(in) :: op
        real(precision), intent(inout) :: u(-halo:n + halo, 3)
        integer, intent(out) :: steps, ncorrect
        real(precision), intent(out) :: t_final, cpu
        integer, intent(out), optional :: ncorr_interp

        real(precision), allocatable :: u0(:, :), k1(:, :), k2(:, :)
        real(precision), allocatable :: k3(:, :), k4(:, :), utemp(:, :)
        real(precision) :: dt, time, t0, t1
        logical :: last
        integer :: nc, nci, ncit

        allocate (u0(-halo:n + halo, 3))
        allocate (k1(-halo:n + halo, 3))
        allocate (k2(-halo:n + halo, 3))
        allocate (k3(-halo:n + halo, 3))
        allocate (k4(-halo:n + halo, 3))
        allocate (utemp(-halo:n + halo, 3))

        call e1d_boundary(n, u)
        steps = 0
        ncorrect = 0
        ncit = 0
        time = 0.0_precision
        call cpu_time(t0)

        do
            dt = e1d_timestep(n, h, u, cfl)
            last = .false.
            if (time + dt >= tend) then
                dt = tend - time
                last = .true.
            end if
            time = time + dt
            steps = steps + 1

            u0 = u

            ! Stage 1.
            call e1d_rhs(n, h, u0, op, k1, nci)
            ncit = ncit + nci
            utemp = u0 + 0.5_precision*dt*k1
            call e1d_correction(n, u0, utemp, nc)
            ncorrect = ncorrect + nc
            call e1d_boundary(n, utemp)

            ! Stage 2.
            call e1d_rhs(n, h, utemp, op, k2, nci)
            ncit = ncit + nci
            utemp = u0 + 0.5_precision*dt*k2
            call e1d_correction(n, u0, utemp, nc)
            ncorrect = ncorrect + nc
            call e1d_boundary(n, utemp)

            ! Stage 3.
            call e1d_rhs(n, h, utemp, op, k3, nci)
            ncit = ncit + nci
            utemp = u0 + dt*k3
            call e1d_correction(n, u0, utemp, nc)
            ncorrect = ncorrect + nc
            call e1d_boundary(n, utemp)

            ! Stage 4.
            call e1d_rhs(n, h, utemp, op, k4, nci)
            ncit = ncit + nci
            utemp = u0 + (k1 + 2.0_precision*k2 + 2.0_precision*k3 + k4) &
                    *dt/6.0_precision
            call e1d_correction(n, u0, utemp, nc)
            ncorrect = ncorrect + nc
            call e1d_boundary(n, utemp)

            u = utemp
            if (last) exit
        end do

        call cpu_time(t1)
        t_final = time
        cpu = t1 - t0
        if (present(ncorr_interp)) ncorr_interp = ncit

        deallocate (u0, k1, k2, k3, k4, utemp)
    end subroutine e1d_solve

end module e1d_solver
