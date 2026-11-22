program test_dispatch
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use scheme_params, only: scheme, moc33_linear, moc53_linear
    use scheme_dispatch
    use scheme_diagnostics_mod
    use test_util
    implicit none

    character(len=32), parameter :: names(23) = (/ &
        'MUSCL        ', 'JS3          ', 'Z3           ', 'JS5          ', 'Z5           ', &
        'QR           ', 'MOC322       ', 'MOC33        ', 'MOC53        ', 'MOZ          ', &
        'LINEAR3_M1P1 ', 'LINEAR3_M2P0 ', 'LINEAR3_0P2  ', 'LINEAR4_M2P1 ', 'LINEAR4_M1P2 ', &
        'LINEAR5_M2P2 ', 'LINEAR6_M2P3 ', 'RBF_GA3      ', 'RBF_MQ3      ', 'RBF_IQ3      ', &
        'RBF_GA5      ', 'RBF_MQ5      ', 'RBF_IQ5      ' /)
    integer, parameter :: expected_order(23) = (/ &
        2, 3, 3, 5, 5, 3, 3, 3, 5, 5, 3, 3, 3, 4, 4, 5, 6, 3, 3, 3, 5, 5, 5 /)

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-3:4), u_half, exact
    type(scheme_diagnostics) :: diag
    character(len=64) :: out_name
    integer :: i, l

    do i = 1, size(names)
        call check_true(scheme_known(trim(names(i))), 'scheme_known ' // trim(names(i)))
        call check_true(scheme_order(trim(names(i))) == expected_order(i), 'scheme_order ' // trim(names(i)))
    end do
    call check_true(.not. scheme_known('NOPE'), 'unknown scheme rejected')

    do i = 1, size(names)
        scheme = trim(names(i))
        do l = -3, 4
            u(l) = 0.3_precision - 1.2_precision * real(l, precision) * h
        end do
        exact = 0.3_precision - 1.2_precision * 0.5_precision * h
        call interpolate_dispatch(u, u_half)
        call check_true(ieee_is_finite(u_half), 'dispatch finite ' // trim(names(i)))
        call check_close(u_half, exact, test_tol, 'dispatch linear exact ' // trim(names(i)))
    end do

    scheme = 'Z5'
    do l = -3, 4
        u(l) = cos(0.3_precision + real(l, precision) * h)
    end do
    call interpolate_dispatch(u, u_half, diag)
    call check_true(diag%omega(1) /= 0.0_precision, 'dispatch fills diag when present')
    call interpolate_dispatch(u, u_half)
    call check_true(ieee_is_finite(u_half), 'dispatch works without diag')

    call check_true(trim(scheme_output_name('RBF_GA5')) == 'RBF_GA5', 'output name RBF_GA5')
    call check_true(trim(scheme_output_name('RBF_MQ3')) == 'RBF_MQ3', 'output name RBF_MQ3')
    call check_true(trim(scheme_output_name('MOZ')) == 'MOZ', 'output name MOZ')

    moc53_linear = 'GA'
    out_name = scheme_output_name('MOC53')
    call check_true(trim(out_name) == 'MOC53_GA', 'output name MOC53_GA')

    moc33_linear = 'LINEAR'
    out_name = scheme_output_name('MOC33')
    call check_true(trim(out_name) == 'MOC33_LINEAR', 'output name MOC33_LINEAR')

    call finish()

end program test_dispatch
