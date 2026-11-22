program test_1d_core
    ! Unit tests for the 1D initial profiles, RHS, and one-level solver.
    use precision_select
    use scheme_params, only: scheme
    use accuracy_order, only: error_norms_1d
    use c1d_initial_mod, only: c1d_profile, c1d_initial
    use c1d_exact_mod, only: c1d_exact
    use c1d_solver, only: c1d_boundary, c1d_rhs, c1d_solve_level
    use test_util
    implicit none

    integer, parameter :: n = 20
    real(precision), parameter :: pi = acos(-1.0_precision)
    real(precision) :: u(-7:n + 7), uex(0:n), f(0:n)
    real(precision) :: uend(0:40), uex40(0:40)
    real(precision) :: e1, e2, einf, cpu
    integer :: j

    scheme = 'Z3'

    call c1d_initial(n, 'sin2', u)
    call check_close(u(15), 1.0_precision, test_tol, 'initial peak')
    call check_close(c1d_profile(0.0_precision, 'gauss'), 1.0_precision, test_tol, 'gauss')
    call check_close(c1d_profile(0.0_precision, 'rational'), 1.0_precision, test_tol, 'rational')
    call check_close(c1d_profile(0.5_precision, 'sin'), 1.0_precision, test_tol, 'sin')
    call check_close(c1d_profile(0.5_precision, 'sin3'), 1.0_precision, test_tol, 'sin3')
    call check_close(c1d_profile(0.5_precision, 'sin_phase'), &
                     sin(pi*0.5_precision - sin(pi*0.5_precision)/pi), test_tol, 'sin_phase')

    call c1d_exact(n, 'sin2', 0.0_precision, uex)
    call check_close(uex(15), u(15), test_tol, 'exact at t=0')

    u = 1.0_precision
    call c1d_boundary(n, u)
    call c1d_rhs(n, u, 'HCS4', f)
    do j = 0, n
        call check_close(f(j), 0.0_precision, 1.0e-12_precision, 'constant rhs')
    end do

    call c1d_solve_level(40, 'sin2', 'HCS4', 0.1_precision, 0.2_precision, uend, cpu)
    call c1d_exact(40, 'sin2', 0.2_precision, uex40)
    call error_norms_1d(40, uend, uex40, e1, e2, einf)
    call check_true(e2 < 5.0e-2_precision, '1d smooth error')

    call finish()
end program test_1d_core
