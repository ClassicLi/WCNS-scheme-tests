program test_rhs
    ! Analytic checks for the assembled WCNS spatial right-hand side.
    !
    ! The interpolation symbol is measured from the shared dispatch itself,
    ! so these checks isolate the operator assembly: half-node interpolation
    ! followed by one of the five difference operators.
    use precision_select
    use scheme_params, only: scheme
    use scheme_dispatch, only: interpolate_dispatch
    use cw_params, only: difference_operator
    use cw_wcns_rhs, only: wcns_rhs
    use test_util
    implicit none

    real(precision), parameter :: PI = acos(-1.0_precision)
    integer(integer_precision), parameter :: nx = 33, n = 4, halo = 7
    real(precision) :: u_cos(-halo:nx+halo), u_sin(-halo:nx+halo)
    real(precision) :: u_const(-halo:nx+halo)
    real(precision) :: rhs_cos(-halo:nx+halo), rhs_sin(-halo:nx+halo)
    real(precision) :: rhs_const(-halo:nx+halo)
    real(precision) :: u_half_cos, u_half_sin, theta, h
    complex(precision) :: II, K, expected, actual
    character(len=8) :: ops(5)
    integer :: m
    integer(integer_precision) :: j

    II = cmplx(0.0_precision, 1.0_precision, precision)
    theta = 2.0_precision*PI*real(n, precision)/real(nx, precision)
    h = 2.0_precision/real(nx, precision)

    ops = [character(len=8) :: 'WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6']

    do j = -halo, nx + halo
        u_cos(j) = cos(theta*real(j, precision))
        u_sin(j) = sin(theta*real(j, precision))
    end do

    scheme = 'LINEAR5_M2P2'

    do m = 1, 5
        difference_operator = trim(ops(m))

        ! Measure the interpolation symbol from the same shared dispatch.
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
            ! The integer-node term is not multiplied by the interpolation symbol.
            expected = -((4.0_precision/3.0_precision)*K*(2.0_precision*II*sin(theta/2.0_precision)) &
                - (1.0_precision/6.0_precision)*(2.0_precision*II*sin(theta)))/h
        case ('HCS6')
            expected = -((64.0_precision/45.0_precision)*K*(2.0_precision*II*sin(theta/2.0_precision)) &
                - (2.0_precision/9.0_precision)*(2.0_precision*II*sin(theta)) &
                + (1.0_precision/180.0_precision)*(2.0_precision*II*sin(2.0_precision*theta)))/h
        end select

        actual = cmplx(rhs_cos(0), rhs_sin(0), precision)
        call check_close(real(actual, precision), real(expected, precision), test_tol, &
                         trim(ops(m)) // ' rhs real')
        call check_close(aimag(actual), aimag(expected), test_tol, &
                         trim(ops(m)) // ' rhs imag')
    end do

    ! Every interpolation must preserve a constant, so its derivative is zero.
    scheme = 'Z5'
    u_const = 3.5_precision
    do m = 1, 5
        difference_operator = trim(ops(m))
        call wcns_rhs(nx, halo, u_const, rhs_const)
        call check_true(maxval(abs(rhs_const(0:nx))) < test_tol, &
                        trim(ops(m)) // ' constant preserved')
    end do

    call finish()
end program test_rhs
