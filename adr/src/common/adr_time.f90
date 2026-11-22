module adr_time
    ! Three-stage third-order Runge-Kutta step with periodic halo refresh.
    use precision_select
    use adr_periodic, only: boundary_cycle
    implicit none

    abstract interface
        subroutine rhs_iface(nx, halo, u, f)
            ! Spatial right-hand side used by rk3_step.
            import :: precision, integer_precision
            integer(integer_precision), intent(in) :: nx, halo
            real(precision), intent(in) :: u(-halo:nx+halo)
            real(precision), intent(out) :: f(-halo:nx+halo)
        end subroutine rhs_iface
    end interface

contains

    subroutine rk3_step(nx, halo, dt, u0, rhs_op, u1)
        ! Advance one step with the three-stage third-order RK scheme.
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: dt
        real(precision), intent(in) :: u0(0:nx)
        procedure(rhs_iface) :: rhs_op
        real(precision), intent(out) :: u1(0:nx)

        real(precision) :: un(-halo:nx+halo)
        real(precision) :: utemp(-halo:nx+halo)
        real(precision) :: f(-halo:nx+halo)

        un(0:nx) = u0

        call boundary_cycle(nx, un, halo)
        call rhs_op(nx, halo, un, f)
        utemp = un + dt * f

        call boundary_cycle(nx, utemp, halo)
        call rhs_op(nx, halo, utemp, f)
        utemp = 0.75_precision * un + 0.25_precision * utemp + 0.25_precision * dt * f

        call boundary_cycle(nx, utemp, halo)
        call rhs_op(nx, halo, utemp, f)
        utemp = (un + 2.0_precision * utemp + 2.0_precision * dt * f) / 3.0_precision

        u1 = utemp(0:nx)
    end subroutine rk3_step

end module adr_time
