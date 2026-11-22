program test_2d_core
    ! Unit tests for the 2D residual, periodic boundary, and one RK4 step.
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use e2d_globals
    use e2d_cases, only: e2d_initial
    use e2d_solver, only: e2d_boundary, e2d_rhs, e2d_solve_level
    use test_util
    implicit none

    integer :: ierr, j, k
    character(len=256) :: msg
    real(precision) :: rhs(-mt:8 + mt, -mt:8 + mt, 4), cpu

    call e2d_allocate(8, 10.0_precision)
    call e2d_initial('vortex', 0.0_precision, 0.0_precision, ierr, msg)

    ! ---------------------------------------------------------------
    ! Constant state must give a zero residual
    ! ---------------------------------------------------------------
    ubar = 0.0_precision
    ubar(:, :, 1) = 1.0_precision
    ubar(:, :, 4) = 1.0_precision/beta

    call e2d_boundary()
    call e2d_rhs('HCS4', rhs)
    call check_close(maxval(abs(rhs(0:n, 0:m, :))), 0.0_precision, &
                     1.0e-12_precision, 'constant rhs')

    ! ---------------------------------------------------------------
    ! Periodic halo
    ! ---------------------------------------------------------------
    ubar = 0.0_precision
    do k = -mt, m + mt
        do j = -mt, n + mt
            ubar(j, k, 1) = real(j, precision)
        end do
    end do
    call e2d_boundary()
    call check_close(ubar(-1, 0, 1), ubar(n - 1, 0, 1), 1.0e-12_precision, &
                     'periodic x')
    call check_close(ubar(0, -1, 1), ubar(0, m - 1, 1), 1.0e-12_precision, &
                     'periodic y')

    ! ---------------------------------------------------------------
    ! One RK4 step must remain finite
    ! ---------------------------------------------------------------
    call e2d_initial('vortex', 0.0_precision, 0.0_precision, ierr, msg)
    call e2d_boundary()
    call e2d_solve_level('HCS4', 0.1_precision, 0.01_precision, cpu)
    call check_true(all(ieee_is_finite(ubar(0:n, 0:m, :))), 'one step finite')

    call e2d_deallocate()
    call finish()
end program test_2d_core
