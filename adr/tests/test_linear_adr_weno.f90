program test_linear_adr_weno
    ! Analytic checks for the WENO linear reconstruction and flux difference.
    use precision_select
    use weno_params, only: scheme
    use weno_dispatch, only: weno_reconstruct, weno_rhs
    use test_util
    implicit none

    real(precision), parameter :: PI = acos(-1.0_precision)
    integer(integer_precision), parameter :: nx = 33, n = 4, halo = 9
    real(precision) :: u_cos(-halo:nx+halo), u_sin(-halo:nx+halo)
    real(precision) :: rhs_cos(-halo:nx+halo), rhs_sin(-halo:nx+halo)
    real(precision) :: u_poly(-4:4), u_half_cos, u_half_sin, theta, h
    complex(precision) :: II, K, expected, actual
    integer(integer_precision) :: j

    II = cmplx(0.0_precision, 1.0_precision, precision)
    theta = 2.0_precision * PI * real(n, precision) / real(nx, precision)
    h = 2.0_precision * PI / real(nx, precision)

    do j = -halo, nx + halo
        u_cos(j) = cos(theta * real(j, precision))
        u_sin(j) = sin(theta * real(j, precision))
    end do

    scheme = 'WENO3_LINEAR'
    call weno_reconstruct(u_cos(-4:4), u_half_cos, h)
    call weno_reconstruct(u_sin(-4:4), u_half_sin, h)
    K = cmplx(u_half_cos, u_half_sin, precision) / exp(II*theta/2.0_precision)
    call weno_rhs(nx, halo, u_cos, rhs_cos)
    call weno_rhs(nx, halo, u_sin, rhs_sin)
    expected = -2.0_precision * II * K * sin(theta/2.0_precision) / h
    actual = cmplx(rhs_cos(0), rhs_sin(0), precision)
    call check_close(real(actual, precision), real(expected, precision), test_tol, &
                     'WENO3 linear rhs real')
    call check_close(aimag(actual), aimag(expected), test_tol, 'WENO3 linear rhs imag')

    scheme = 'WENO5_LINEAR'
    call weno_reconstruct(u_cos(-4:4), u_half_cos, h)
    call weno_reconstruct(u_sin(-4:4), u_half_sin, h)
    K = cmplx(u_half_cos, u_half_sin, precision) / exp(II*theta/2.0_precision)
    call weno_rhs(nx, halo, u_cos, rhs_cos)
    call weno_rhs(nx, halo, u_sin, rhs_sin)
    expected = -2.0_precision * II * K * sin(theta/2.0_precision) / h
    actual = cmplx(rhs_cos(0), rhs_sin(0), precision)
    call check_close(real(actual, precision), real(expected, precision), test_tol, &
                     'WENO5 linear rhs real')
    call check_close(aimag(actual), aimag(expected), test_tol, 'WENO5 linear rhs imag')

    scheme = 'WENO3_LINEAR'
    u_poly = 0.0_precision
    do j = -2, 3
        u_poly(j) = real(j, precision)
    end do
    call weno_reconstruct(u_poly, u_half_cos, h)
    call check_close(u_half_cos, 0.5_precision, test_tol, 'WENO3 linear polynomial exact')

    scheme = 'WENO5_LINEAR'
    u_poly = 0.0_precision
    do j = -2, 2
        u_poly(j) = real(j, precision)
    end do
    call weno_reconstruct(u_poly, u_half_cos, h)
    call check_close(u_half_cos, 0.5_precision, test_tol, 'WENO5 linear polynomial exact')

    call finish()
end program test_linear_adr_weno
