program test_linear_adr_wcns
    ! Analytic checks for WCNS linear interpolation plus each difference operator.
    use precision_select
    use scheme_params, only: scheme
    use scheme_dispatch, only: interpolate_dispatch
    use wcns_params, only: difference_operator
    use wcns_dispatch, only: wcns_rhs
    use test_util
    implicit none

    real(precision), parameter :: PI = acos(-1.0_precision)
    integer(integer_precision), parameter :: nx = 33, n = 4, halo = 7
    real(precision) :: u_cos(-halo:nx+halo), u_sin(-halo:nx+halo)
    real(precision) :: rhs_cos(-halo:nx+halo), rhs_sin(-halo:nx+halo)
    real(precision) :: u_half_cos, u_half_sin, theta, h
    complex(precision) :: II, K, expected, actual
    character(len=8) :: ops(4)
    integer :: m
    integer(integer_precision) :: j

    II = cmplx(0.0_precision, 1.0_precision, precision)
    theta = 2.0_precision*PI*real(n, precision)/real(nx, precision)
    h = 2.0_precision*PI/real(nx, precision)

    ops = [character(len=8) :: 'WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4']

    do j = -halo, nx + halo
        u_cos(j) = cos(theta*real(j, precision))
        u_sin(j) = sin(theta*real(j, precision))
    end do

    scheme = 'LINEAR5_M2P2'

    do m = 1, 4
        difference_operator = trim(ops(m))

        call interpolate_dispatch(u_cos(-3:4), u_half_cos)
        call interpolate_dispatch(u_sin(-3:4), u_half_sin)
        K = cmplx(u_half_cos, u_half_sin, precision)/exp(II*theta/2.0_precision)

        call wcns_rhs(nx, halo, u_cos, rhs_cos)
        call wcns_rhs(nx, halo, u_sin, rhs_sin)

        select case (trim(ops(m)))
        case ('WCNSE2')
            expected = K*(-2.0_precision*II*sin(theta/2.0_precision)/h)
        case ('WCNSE4')
            expected = K*(-((9.0_precision/8.0_precision)*(2.0_precision*II*sin(theta/2.0_precision)) &
                - (1.0_precision/24.0_precision)*(2.0_precision*II*sin(1.5_precision*theta)))/h)
        case ('WCNSE6')
            expected = K*(-((75.0_precision/64.0_precision)*(2.0_precision*II*sin(theta/2.0_precision)) &
                - (25.0_precision/384.0_precision)*(2.0_precision*II*sin(1.5_precision*theta)) &
                + (3.0_precision/640.0_precision)*(2.0_precision*II*sin(2.5_precision*theta)))/h)
        case ('HCS4')
            ! The u_node term uses integer-node u, so it is not multiplied by K.
            expected = -((4.0_precision/3.0_precision)*K*(2.0_precision*II*sin(theta/2.0_precision)) &
                - (1.0_precision/6.0_precision)*(2.0_precision*II*sin(theta)))/h
        end select

        actual = cmplx(rhs_cos(0), rhs_sin(0), precision)
        call check_close(real(actual, precision), real(expected, precision), test_tol, &
                         trim(ops(m)) // ' rhs real')
        call check_close(aimag(actual), aimag(expected), test_tol, trim(ops(m)) // ' rhs imag')
    end do

    call finish()
end program test_linear_adr_wcns
