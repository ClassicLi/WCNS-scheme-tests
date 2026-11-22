program test_rk3
    ! Unit test for the three-stage third-order RK step.
    use precision_select
    use adr_time, only: rk3_step
    use test_util
    implicit none

    real(precision) :: u0(0:1), u1(0:1)
    real(precision), parameter :: dt = 0.1_precision

    u0 = 1.0_precision
    call rk3_step(1_integer_precision, 1_integer_precision, dt, u0, rhs_decay, u1)
    call check_close(u1(0), exp(-dt), 1.0e-5_precision, 'rk3 exponential decay')

    call finish()

contains

    subroutine rhs_decay(nx, halo, u, f)
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: u(-halo:nx+halo)
        real(precision), intent(out) :: f(-halo:nx+halo)

        f = -u
    end subroutine rhs_decay

end program test_rk3
