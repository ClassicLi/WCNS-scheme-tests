program test_periodic_dft
    ! Unit tests for the periodic halo, DFT, and Phi extraction.
    use precision_select
    use adr_periodic, only: boundary_cycle
    use adr_dft, only: dft_mode, calc_phi
    use test_util
    implicit none

    real(precision), parameter :: PI = acos(-1.0_precision)
    real(precision) :: u(-2:6), v(0:17)
    complex(precision) :: II, hat, hat0, hat_t, phi
    real(precision) :: omega, dt
    integer(integer_precision) :: n_num, n0, j

    II = cmplx(0.0_precision, 1.0_precision, precision)

    u = 0.0_precision
    u(0:4) = [1.0_precision, 2.0_precision, 3.0_precision, 4.0_precision, 5.0_precision]
    call boundary_cycle(4_integer_precision, u, 2_integer_precision)
    call check_true(u(5) == 2.0_precision .and. u(6) == 3.0_precision .and. &
                    u(-1) == 4.0_precision .and. u(-2) == 3.0_precision, 'periodic halo')

    n_num = 17_integer_precision
    n0 = 3_integer_precision
    v(0:n_num) = [( cos(real(n0*j, precision)*2.0_precision*PI/real(n_num, precision)), &
                    j=0, n_num )]
    call dft_mode(v, n0, n_num, hat)
    call check_close(real(hat, precision), 0.5_precision, test_tol, 'dft cos real')
    call check_close(aimag(hat), 0.0_precision, test_tol, 'dft cos imag')

    hat0 = cmplx(1.0_precision, 0.0_precision, precision)
    omega = 0.5_precision
    dt = 1.0e-3_precision
    hat_t = hat0 * exp(-II*omega*dt)
    call calc_phi(hat_t, hat0, dt, n_num, phi)
    call check_close(real(phi, precision), 2.0_precision*PI*omega/real(n_num, precision), &
                     test_tol, 'phi real')
    call check_close(aimag(phi), 0.0_precision, test_tol, 'phi imag')

    call finish()
end program test_periodic_dft
