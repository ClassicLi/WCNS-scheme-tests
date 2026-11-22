module cw_time
    ! Four-stage fourth-order Runge-Kutta time advance with periodic halos.
    !
    ! The scheme is the classical RK4 used by the reference program:
    !   k1 = RHS(u)
    !   k2 = RHS(u + dt/2 k1)
    !   k3 = RHS(u + dt/2 k2)
    !   k4 = RHS(u + dt   k3)
    !   u1 = u + dt (k1 + 2 k2 + 2 k3 + k4) / 6
    use precision_select
    use cw_periodic, only: boundary_cycle
    implicit none

    abstract interface
        subroutine rhs_iface(nx, halo, u, f)
            ! Spatial right-hand side used by rk4_step.
            import :: precision, integer_precision
            integer(integer_precision), intent(in) :: nx, halo
            real(precision), intent(in) :: u(-halo:nx+halo)
            real(precision), intent(out) :: f(-halo:nx+halo)
        end subroutine rhs_iface
    end interface

contains

    subroutine rk4_step(nx, halo, dt, u0, rhs_op, u1)
        ! Advance one step and return only the interior u1(0:nx).
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: dt
        real(precision), intent(in) :: u0(0:nx)
        procedure(rhs_iface) :: rhs_op
        real(precision), intent(out) :: u1(0:nx)

        real(precision) :: un(-halo:nx+halo)
        real(precision) :: utemp(-halo:nx+halo)
        real(precision) :: k1(-halo:nx+halo), k2(-halo:nx+halo)
        real(precision) :: k3(-halo:nx+halo), k4(-halo:nx+halo)

        ! Stage 0: copy the interior and refresh the periodic halo.
        un = 0.0_precision
        un(0:nx) = u0
        call boundary_cycle(nx, un, halo)

        ! Stage 1.
        call rhs_op(nx, halo, un, k1)
        utemp = un + 0.5_precision*dt*k1
        call boundary_cycle(nx, utemp, halo)

        ! Stage 2.
        call rhs_op(nx, halo, utemp, k2)
        utemp = un + 0.5_precision*dt*k2
        call boundary_cycle(nx, utemp, halo)

        ! Stage 3.
        call rhs_op(nx, halo, utemp, k3)
        utemp = un + dt*k3
        call boundary_cycle(nx, utemp, halo)

        ! Stage 4 and final weighted average.
        call rhs_op(nx, halo, utemp, k4)
        u1 = u0 + dt*(k1(0:nx) + 2.0_precision*k2(0:nx) &
                    + 2.0_precision*k3(0:nx) + k4(0:nx))/6.0_precision
    end subroutine rk4_step

end module cw_time
