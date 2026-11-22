program test_solver
    ! Unit tests for the 2D Euler CFL time step, residual, positivity
    ! correction, and the fixed four-stage RK4 advance.  A uniform state is
    ! used so that the residual and the advance have known answers.
    use precision_select
    use e2d_globals, only: mt, n, m, h, beta, ubar, e2d_allocate, &
                           e2d_deallocate, e2d_set_gama
    use e2d_cases, only: e2d_boundary
    use e2d_solver, only: e2d_rhs, e2d_timestep, e2d_apply_correction, &
                         e2d_advance
    use test_util
    implicit none

    integer(integer_precision), parameter :: nz = 16_integer_precision

    real(precision), allocatable :: rhs(:, :, :), u0(:, :, :), initial(:, :, :)
    real(precision) :: rho, uu, ww, ps, energy
    real(precision) :: dt, dt_expected, err_max, cpu_seconds
    integer :: steps, interp_corr, cell_corr, ncorr, p

    call e2d_set_gama(1.4_precision)
    call e2d_allocate(nz, nz)

    ! Uniform state: rho = 1, rest velocity, p = 1.
    rho = 1.0_precision
    uu = 0.0_precision
    ww = 0.0_precision
    ps = 1.0_precision
    energy = ps/beta + 0.5_precision*rho*(uu*uu + ww*ww)

    ubar(:, :, 1) = rho
    ubar(:, :, 2) = rho*uu
    ubar(:, :, 3) = rho*ww
    ubar(:, :, 4) = energy

    allocate (rhs(-mt:n + mt, -mt:m + mt, 4))
    allocate (u0(-mt:n + mt, -mt:m + mt, 4))
    allocate (initial(-mt:n + mt, -mt:m + mt, 4))

    ! ---------------------------------------------------------------
    ! 1. CFL time step on the uniform state.
    ! ---------------------------------------------------------------
    dt = e2d_timestep(0.3_precision)
    dt_expected = 0.3_precision*h/sqrt(1.4_precision)
    call check_close(dt, dt_expected, 1.0e-14_precision, 'CFL timestep')

    ! ---------------------------------------------------------------
    ! 2. Uniform state has a zero residual in all four components.
    ! ---------------------------------------------------------------
    call e2d_boundary('riemann1', 0.0_precision)
    call e2d_rhs('riemann1', 'HCS4', rhs, ncorr)

    call check_true(maxval(abs(rhs(0:n, 0:m, 1))) < 1.0e-12_precision, &
                    'uniform rhs rho')
    call check_true(maxval(abs(rhs(0:n, 0:m, 2))) < 1.0e-12_precision, &
                    'uniform rhs x-momentum')
    call check_true(maxval(abs(rhs(0:n, 0:m, 3))) < 1.0e-12_precision, &
                    'uniform rhs z-momentum')
    call check_true(maxval(abs(rhs(0:n, 0:m, 4))) < 1.0e-12_precision, &
                    'uniform rhs energy')

    ! ---------------------------------------------------------------
    ! 3. Positivity correction replaces one non-physical cell by the
    !    average of its four neighbours from the saved state.
    ! ---------------------------------------------------------------
    u0 = ubar
    ubar(0, 0, 1) = -1.0_precision
    call e2d_apply_correction(u0, ncorr)

    call check_true(ncorr == 1, 'correction count')
    call check_true(ubar(0, 0, 1) > 0.0_precision, 'corrected density positive')
    call check_close(ubar(0, 0, 1), rho, 1.0e-14_precision, &
                     'correction averages neighbours')

    ! ---------------------------------------------------------------
    ! 4. A uniform state is preserved by the RK4 advance.
    ! ---------------------------------------------------------------
    ubar(:, :, 1) = rho
    ubar(:, :, 2) = rho*uu
    ubar(:, :, 3) = rho*ww
    ubar(:, :, 4) = energy
    initial = ubar

    call e2d_advance('riemann1', 'HCS6', 0.1_precision, 0.0_precision, &
                     1.0e-6_precision, steps, interp_corr, cell_corr, cpu_seconds)

    call check_true(steps >= 1, 'advance takes at least one step')
    err_max = maxval(abs(ubar(0:n, 0:m, :) - initial(0:n, 0:m, :)))
    call check_true(err_max < 1.0e-12_precision, 'uniform state preserved by RK4')

    call e2d_deallocate()
    call finish()

end program test_solver
