program test_difference_operators
    ! Fourier-symbol checks for the five WCNS difference operators.
    !
    ! The half-node flux and integer-node values are sampled from cos/sin
    ! modes.  The imaginary part of the operator symbol must match the
    ! analytic discrete symbol, and the real part must vanish.
    use precision_select
    use difference_operators, only: difference_apply
    use test_util
    implicit none

    real(precision), parameter :: theta = 0.7_precision, h = 0.1_precision
    character(len=8) :: ops(5)
    real(precision) :: flux_cos(1:6), flux_sin(1:6)
    real(precision) :: u_cos(1:5), u_sin(1:5)
    real(precision) :: d_cos, d_sin, expected_im
    integer :: k, m

    ops = [character(len=8) :: 'WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6']

    do m = 1, 5
        ! flux_half(1:6) spans half nodes j-3/2 .. j+5/2.
        do k = 1, 6
            flux_cos(k) = cos(theta*(real(k, precision) - 3.5_precision))
            flux_sin(k) = sin(theta*(real(k, precision) - 3.5_precision))
        end do
        ! u_node(1:5) spans integer nodes j-2 .. j+2.
        do k = 1, 5
            u_cos(k) = cos(theta*real(k - 3, precision))
            u_sin(k) = sin(theta*real(k - 3, precision))
        end do

        call difference_apply(trim(ops(m)), h, flux_cos, u_cos, d_cos)
        call difference_apply(trim(ops(m)), h, flux_sin, u_sin, d_sin)

        select case (trim(ops(m)))
        case ('WCNSE2')
            expected_im = -(2.0_precision*sin(theta/2.0_precision))/h
        case ('WCNSE4')
            expected_im = -((9.0_precision/8.0_precision)*(2.0_precision*sin(theta/2.0_precision)) &
                - (1.0_precision/24.0_precision)*(2.0_precision*sin(1.5_precision*theta)))/h
        case ('WCNSE6')
            expected_im = -((75.0_precision/64.0_precision)*(2.0_precision*sin(theta/2.0_precision)) &
                - (25.0_precision/384.0_precision)*(2.0_precision*sin(1.5_precision*theta)) &
                + (3.0_precision/640.0_precision)*(2.0_precision*sin(2.5_precision*theta)))/h
        case ('HCS4')
            expected_im = -((4.0_precision/3.0_precision)*(2.0_precision*sin(theta/2.0_precision)) &
                - (1.0_precision/6.0_precision)*(2.0_precision*sin(theta)))/h
        case ('HCS6')
            expected_im = -((64.0_precision/45.0_precision)*(2.0_precision*sin(theta/2.0_precision)) &
                - (2.0_precision/9.0_precision)*(2.0_precision*sin(theta)) &
                + (1.0_precision/180.0_precision)*(2.0_precision*sin(2.0_precision*theta)))/h
        end select

        call check_close(d_cos, 0.0_precision, test_tol, trim(ops(m)) // ' cos real')
        call check_close(d_sin, expected_im, test_tol, trim(ops(m)) // ' sin imag')
    end do

    call finish()
end program test_difference_operators
