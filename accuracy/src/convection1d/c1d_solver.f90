module c1d_solver
    ! Spatial right-hand side and RK4 time stepping for the 1D test.
    !
    ! The residual uses the shared WCNS interpolation standard parts and the
    ! selected difference operator.  The half-node fluxes live at flux(k) for
    ! x = -1 + (k + 1/2) h.
    use precision_select
    use scheme_dispatch, only: interpolate_dispatch
    use difference_operators, only: difference_apply
    use c1d_initial_mod, only: c1d_initial
    implicit none

    integer(integer_precision), parameter :: halo = 7_integer_precision

contains

    subroutine c1d_boundary(n, u)
        ! Periodic halo fill.
        integer, intent(in) :: n
        real(precision), intent(inout) :: u(-halo:n + halo)

        integer :: i

        do i = 1, halo
            u(-i) = u(n - i)
            u(n + i) = u(i)
        end do
    end subroutine c1d_boundary

    subroutine c1d_rhs(n, u, op, f)
        ! f(0:n) = -dE/dx using WCNS interpolation and the selected operator.
        integer, intent(in) :: n
        real(precision), intent(in) :: u(-halo:n + halo)
        character(len=*), intent(in) :: op
        real(precision), intent(out) :: f(0:n)

        real(precision) :: flux(-3:n + 2)
        real(precision) :: h
        integer :: j, k

        h = 2.0_precision/real(n, precision)

        do k = -3, n + 2
            call interpolate_dispatch(u(k - 3:k + 4), flux(k))
        end do

        do j = 0, n
            call difference_apply(op, h, flux(j - 3:j + 2), u(j - 2:j + 2), f(j))
        end do
    end subroutine c1d_rhs

    subroutine c1d_solve_level(n, initial_name, op, cfl, tend, uend, cpu_seconds)
        ! Advance one grid level to tend and return u(0:n) and CPU time.
        integer, intent(in) :: n
        character(len=*), intent(in) :: initial_name, op
        real(precision), intent(in) :: cfl, tend
        real(precision), intent(out) :: uend(0:n)
        real(precision), intent(out) :: cpu_seconds

        real(precision), allocatable :: u(:), u1(:)
        real(precision), allocatable :: k1(:), k2(:), k3(:), k4(:)
        real(precision) :: h, time, dt, umax, t0, t1
        logical :: last

        h = 2.0_precision/real(n, precision)

        allocate (u(-halo:n + halo), u1(-halo:n + halo))
        allocate (k1(0:n), k2(0:n), k3(0:n), k4(0:n))

        call c1d_initial(n, initial_name, u)
        call c1d_boundary(n, u)

        time = 0.0_precision
        call cpu_time(t0)

        do
            umax = maxval(abs(u(0:n)))
            if (umax <= 0.0_precision) error stop 'c1d_solve_level: zero solution'

            dt = cfl*h/umax
            last = .false.
            if (time + dt >= tend) then
                dt = tend - time
                last = .true.
            end if
            time = time + dt

            call c1d_rhs(n, u, op, k1)
            u1(0:n) = u(0:n) + dt*0.5_precision*k1
            call c1d_boundary(n, u1)

            call c1d_rhs(n, u1, op, k2)
            u1(0:n) = u(0:n) + dt*0.5_precision*k2
            call c1d_boundary(n, u1)

            call c1d_rhs(n, u1, op, k3)
            u1(0:n) = u(0:n) + dt*k3
            call c1d_boundary(n, u1)

            call c1d_rhs(n, u1, op, k4)
            u(0:n) = u(0:n) + dt*(k1 + 2.0_precision*k2 + 2.0_precision*k3 + k4) &
                     /6.0_precision
            call c1d_boundary(n, u)

            if (last) exit
        end do

        call cpu_time(t1)
        cpu_seconds = t1 - t0

        uend = u(0:n)
        deallocate (u, u1, k1, k2, k3, k4)
    end subroutine c1d_solve_level

end module c1d_solver
