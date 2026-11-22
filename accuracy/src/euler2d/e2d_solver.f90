module e2d_solver
    ! Boundary conditions, residual assembly, and RK4 stepping for 2D Euler.
    use precision_select
    use e2d_globals
    use e2d_flux, only: u_to_prim_e_g, interflux
    use difference_operators, only: difference_apply
    implicit none

contains

    subroutine e2d_boundary()
        ! Periodic halo fill in both directions.
        integer :: j, k, p

        do p = 1, mt
            do j = -mt, n + mt
                ubar(j, m + p, :) = ubar(j, p, :)
            end do
            do j = -mt, n + mt
                ubar(j, -p, :) = ubar(j, m - p, :)
            end do
            do k = -mt, m + mt
                ubar(-p, k, :) = ubar(n - p, k, :)
            end do
            do k = -mt, m + mt
                ubar(n + p, k, :) = ubar(p, k, :)
            end do
        end do
    end subroutine e2d_boundary

    subroutine e2d_rhs(op, rhs)
        ! rhs = -dE/dx - dG/dy using the selected difference operator.
        !
        ! In the hybrid compact operators the integer-node argument is the
        ! physical flux (Eint/Gint), not the conserved state ubar.
        character(len=*), intent(in) :: op
        real(precision), intent(out) :: rhs(-mt:n + mt, -mt:m + mt, 4)

        real(precision) :: kx, ky
        integer :: j, k, p

        call u_to_prim_e_g()
        call interflux()

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
    end subroutine e2d_rhs

    function e2d_timestep(cfl) result(dt)
        ! CFL-limited time step from the maximum local wave speed.
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
                temp = sqrt(uu*uu + ww*ww) + sqrt(ps/ro)
                if (temp > maxtemp) maxtemp = temp
            end do
        end do

        if (maxtemp <= 0.0_precision) error stop 'e2d_timestep: non-positive speed'
        dt = cfl*h/maxtemp
    end function e2d_timestep

    subroutine e2d_solve_level(op, cfl, tend, cpu_seconds)
        ! Advance the current state to tend with four-stage RK4.
        character(len=*), intent(in) :: op
        real(precision), intent(in) :: cfl, tend
        real(precision), intent(out) :: cpu_seconds

        real(precision), allocatable :: u0(:, :, :), k1(:, :, :), k2(:, :, :)
        real(precision), allocatable :: k3(:, :, :), k4(:, :, :)
        real(precision) :: dt, time, t0, t1
        logical :: last

        allocate (u0(-mt:n + mt, -mt:m + mt, 4))
        allocate (k1(-mt:n + mt, -mt:m + mt, 4))
        allocate (k2(-mt:n + mt, -mt:m + mt, 4))
        allocate (k3(-mt:n + mt, -mt:m + mt, 4))
        allocate (k4(-mt:n + mt, -mt:m + mt, 4))

        u0 = ubar
        time = 0.0_precision
        call cpu_time(t0)

        do
            dt = e2d_timestep(cfl)
            last = .false.
            if (time + dt >= tend) then
                dt = tend - time
                last = .true.
            end if
            time = time + dt

            call e2d_boundary()
            call e2d_rhs(op, k1)
            ubar = u0 + 0.5_precision*dt*k1

            call e2d_boundary()
            call e2d_rhs(op, k2)
            ubar = u0 + 0.5_precision*dt*k2

            call e2d_boundary()
            call e2d_rhs(op, k3)
            ubar = u0 + dt*k3

            call e2d_boundary()
            call e2d_rhs(op, k4)
            ubar = u0 + dt*(k1 + 2.0_precision*k2 + 2.0_precision*k3 + k4) &
                   /6.0_precision

            call e2d_boundary()
            u0 = ubar

            if (last) exit
        end do

        call cpu_time(t1)
        cpu_seconds = t1 - t0

        deallocate (u0, k1, k2, k3, k4)
    end subroutine e2d_solve_level

end module e2d_solver
