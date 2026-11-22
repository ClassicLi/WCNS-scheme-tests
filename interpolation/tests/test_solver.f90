program test_solver
    use precision_select
    use control
    use mesh
    use solver
    use scheme_params, only: scheme
    use test_util
    implicit none

    integer(integer_precision) :: i, n
    real(precision) :: expected_l2

    call check_close(get_dx(10_integer_precision), 0.2_precision, test_tol, 'get_dx(10)')

    call calc_init(10_integer_precision)
    call check_close(integer_point(0)%x, -1.0_precision, test_tol, 'grid left end')
    call check_close(integer_point(10)%x, 1.0_precision, test_tol, 'grid right end')
    call check_close(half_point(0)%x, -0.9_precision, test_tol, 'first half point')
    deallocate(integer_point, half_point)

    ! Mock error/order check: constant exact value, linear error profile.
    iter_times = 2_integer_precision
    call allocate_errors(2)
    allocate(half_point(0:2))
    do i = 0, 2
        half_point(i)%u_exact = 1.0_precision
        half_point(i)%u_interpolate = 1.0_precision - 0.1_precision * real(i + 1, precision)
        half_point(i)%u_linear = 0.0_precision
        half_point(i)%omega_1 = 0.0_precision
        half_point(i)%omega_2 = 0.0_precision
        half_point(i)%omega_3 = 0.0_precision
    end do

    scheme = 'LINEAR3_M1P1'
    call calc_error(3_integer_precision, 1_integer_precision)
    expected_l2 = sqrt(0.14_precision / 3.0_precision)
    call check_close(L1_error(1), 0.2_precision, test_tol, 'mock L1')
    call check_close(L2_error(1), expected_l2, test_tol, 'mock L2')
    call check_close(L_infty_error(1), 0.3_precision, test_tol, 'mock Linf')
    call check_close(Error_nonLinear(1), 0.0_precision, 0.0_precision, 'linear scheme En is zero')

    L1_error(2) = L1_error(1) / 4.0_precision
    L2_error(2) = L2_error(1) / 4.0_precision
    L_infty_error(2) = L_infty_error(1) / 4.0_precision
    Error_nonLinear(2) = 0.0_precision
    Error_nonLinear_inf(2) = 0.0_precision
    call calc_order()
    call check_close(L1_order(2), 2.0_precision, test_tol, 'mock L1 order')
    call check_close(L2_order(2), 2.0_precision, test_tol, 'mock L2 order')
    call check_close(L_infty_order(2), 2.0_precision, test_tol, 'mock Linf order')
    call check_close(nonLinear_order(2), 0.0_precision, 0.0_precision, 'zero nonlinear order handled')
    deallocate(half_point)

    ! End-to-end: a linear fifth-order scheme should converge fifth order and
    ! have zero nonlinear error against its own background.
    ini_mesh_num = 20_integer_precision
    iter_times = 4_integer_precision
    test_func_type = 3_integer_precision
    udf_id = 3_integer_precision
    nonlinear_baseline = 'same_order'

    scheme = 'LINEAR5_M2P2'
    call allocate_errors(4)
    do i = 1, 4
        n = 2**(i - 1) * ini_mesh_num
        call calc_init(n)
        call interpolate(n)
        call calc_error(n, i)
        deallocate(integer_point, half_point)
    end do
    call calc_order()
    call check_true(L1_order(4) >= 4.5_precision, 'LINEAR5 L1 order >= 4.5')
    call check_true(abs(Error_nonLinear(4)) < 1.0e-18_precision, 'LINEAR5 En is zero')

    ! Z5 has a nonzero nonlinear error that decreases with refinement.
    scheme = 'Z5'
    call allocate_errors(4)
    do i = 1, 4
        n = 2**(i - 1) * ini_mesh_num
        call calc_init(n)
        call interpolate(n)
        call calc_error(n, i)
        deallocate(integer_point, half_point)
    end do
    call calc_order()
    call check_true(Error_nonLinear(4) > 0.0_precision, 'Z5 En positive')
#ifdef PRECISION_DP
    ! In double precision the nonlinear error reaches its roundoff floor by N=80.
    call check_true(Error_nonLinear(3) < 1.0e-15_precision, 'Z5 En at roundoff floor (N=80)')
    call check_true(Error_nonLinear(4) <= 1.0e-15_precision, 'Z5 En stays at roundoff floor (N=160)')
#else
    call check_true(Error_nonLinear(4) < Error_nonLinear(3), 'Z5 En decreases')
#endif

    ! delta_omega with exact ideal weights is zero.
    scheme = 'Z3'
    allocate(half_point(0:1))
    do i = 0, 1
        half_point(i)%omega_1 = 0.25_precision
        half_point(i)%omega_2 = 0.75_precision
        half_point(i)%omega_3 = 0.0_precision
    end do
    call calc_delta_omega(2_integer_precision)
    call check_close(half_point(0)%delta_omega_mml, 0.0_precision, 0.0_precision, 'delta_omega 2 mml')
    call check_close(half_point(0)%delta_omega_bjw, 0.0_precision, 0.0_precision, 'delta_omega 2 bjw')
    deallocate(half_point)

    scheme = 'Z5'
    allocate(half_point(0:1))
    do i = 0, 1
        half_point(i)%omega_1 = 0.0625_precision
        half_point(i)%omega_2 = 0.625_precision
        half_point(i)%omega_3 = 0.3125_precision
    end do
    call calc_delta_omega(2_integer_precision)
    call check_close(half_point(0)%delta_omega_mml, 0.0_precision, 0.0_precision, 'delta_omega 3 mml')
    call check_close(half_point(0)%delta_omega_bjw, 0.0_precision, 0.0_precision, 'delta_omega 3 bjw')
    deallocate(half_point)

    test_func_type = 3_integer_precision
    udf_id = 3_integer_precision
    call check_true(trim(func_label()) == 'udf3', 'func_label udf3')
    test_func_type = 1_integer_precision
    n_order_critical_point = 1_integer_precision
    call check_true(trim(func_label()) == 'cp1', 'func_label cp1')

    call finish()

end program test_solver
