program test_rk4
    ! Unit test for the four-stage fourth-order Runge-Kutta step.
    use precision_select
    use cw_time, only: rk4_step
    use test_util
    implicit none

    real(precision) :: u0(0:1), u1(0:1), e1, e2
    real(precision), parameter :: tiny_dt = 0.01_precision

    ! One small step of u' = -u must match exp(-dt) to high accuracy.
    u0 = 1.0_precision
    call rk4_step(1_integer_precision, 1_integer_precision, tiny_dt, u0, rhs_decay, u1)
    call check_close(u1(0), exp(-tiny_dt), 1.0e-12_precision, 'rk4 one step')

    ! A single RK4 step has local truncation error O(dt^5), so halving dt
    ! must reduce the one-step error by about 2^5 = 32.
    u0 = 1.0_precision
    call rk4_step(1_integer_precision, 1_integer_precision, 0.1_precision, u0, rhs_decay, u1)
    e1 = abs(u1(0) - exp(-0.1_precision))
    u0 = 1.0_precision
    call rk4_step(1_integer_precision, 1_integer_precision, 0.05_precision, u0, rhs_decay, u1)
    e2 = abs(u1(0) - exp(-0.05_precision))
    call check_true(e1/e2 > 30.0_precision .and. e1/e2 < 34.0_precision, 'rk4 local order 5')

    call finish()
contains

    subroutine rhs_decay(nx, halo, u, f)
        ! Linear decay right-hand side, independent of the halo.
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: u(-halo:nx+halo)
        real(precision), intent(out) :: f(-halo:nx+halo)

        f = -u
    end subroutine rhs_decay

end program test_rk4
