program test_solver
    ! Unit tests for the 1D Euler boundary, correction, time step, and RK4.
    use precision_select
    use e1d_globals, only: halo
    use e1d_flux, only: e1d_pressure
    use e1d_initial, only: e1d_initial_state
    use e1d_solver
    use scheme_params, only: scheme
    use test_util
    implicit none

    integer, parameter :: n = 20
    real(precision), parameter :: h = 1.0_precision/real(n, precision)

    real(precision) :: u(-halo:n + halo, 3), rhs(-halo:n + halo, 3)
    real(precision) :: u_old(-halo:n + halo, 3), u_new(-halo:n + halo, 3)
    real(precision) :: uniform(-halo:n + halo, 3)
    real(precision) :: t_final, cpu
    integer :: steps, ncorrect, nc
    integer :: j
    logical :: ok

    do j = -halo, n + halo
        uniform(j, :) = (/1.0_precision, 0.3_precision, 2.5_precision/)
    end do

    ! A uniform flow has zero residual.
    u = uniform
    call e1d_rhs(n, h, u, 'WCNSE2', rhs)
    call check_true(maxval(abs(rhs(0:n, :))) < 1.0e-12_precision, 'uniform rhs zero')

    ! Zero-order extrapolation into the halo.
    u = 0.0_precision
    do j = 0, n
        u(j, 1) = real(j, precision)
    end do
    call e1d_boundary(n, u)
    call check_close(u(-1, 1), 0.0_precision, 0.0_precision, 'left halo')
    call check_close(u(-7, 1), 6.0_precision, 0.0_precision, 'left halo far')
    call check_close(u(n + 1, 1), real(n, precision), 0.0_precision, 'right halo')
    call check_close(u(n + 7, 1), real(n - 6, precision), 0.0_precision, 'right halo far')

    ! CFL time step uses sqrt(p/rho), matching e2d_timestep.
    call check_close(e1d_timestep(n, h, uniform, 0.1_precision), &
                     0.1_precision*h/(0.3_precision + sqrt(0.982_precision)), &
                     1.0e-14_precision, 'timestep formula')

    ! Positivity correction replaces a non-physical cell by its neighbours.
    u_old = uniform
    u_new = uniform
    u_new(5, 1) = -1.0_precision
    call e1d_correction(n, u_old, u_new, nc)
    call check_true(nc == 1, 'correction count')
    call check_close(u_new(5, 1), 1.0_precision, 1.0e-14_precision, 'correction average')

    ! A uniform flow is preserved by RK4, and the last step lands on tend.
    u = uniform
    call e1d_solve(n, h, 0.1_precision, 1.0e-4_precision, 'WCNSE2', u, &
                   steps, ncorrect, t_final, cpu)
    call check_close(t_final, 1.0e-4_precision, 1.0e-12_precision, 'final time clipped')
    call check_true(maxval(abs(u(0:n, :) - uniform(0:n, :))) < 1.0e-10_precision, &
                    'uniform state preserved')

    ! A coarse Sod run stays physical.
    call e1d_initial_state(n, 'sod', u)
    scheme = 'QR'
    call e1d_solve(n, h, 0.1_precision, 0.15_precision, 'HCS4', u, &
                   steps, ncorrect, t_final, cpu)
    call check_true(all(u(0:n, 1) > 0.0_precision), 'sod positive density')
    ok = .true.
    do j = 0, n
        ok = ok .and. (e1d_pressure(u(j, :)) > 0.0_precision)
    end do
    call check_true(ok, 'sod positive pressure')

    call finish()

end program test_solver
