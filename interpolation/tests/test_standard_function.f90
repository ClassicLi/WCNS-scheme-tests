program test_standard_function
    use precision_select
    use control
    use standard_function
    use test_util
    implicit none

    integer :: u, ierr
    character(len=128) :: msg

    call check_close(critical_point_func(1.0_precision, 0), 1.0_precision, test_tol, &
        'critical point n=0 at x=1')
    call check_close(critical_point_func(1.0_precision, 1), 1.0_precision, test_tol, &
        'critical point n=1 at x=1')
    call check_close(critical_point_func(1.0_precision, 2), 1.0_precision, test_tol, &
        'critical point n=2 at x=1')
    call check_close(critical_point_func(0.0_precision, 0), 0.0_precision, 0.0_precision, &
        'critical point n=0 at x=0')
    call check_close(critical_point_func(0.0_precision, 1), 0.0_precision, 0.0_precision, &
        'critical point n=1 at x=0')
    call check_close(critical_point_func(0.0_precision, 2), 0.0_precision, 0.0_precision, &
        'critical point n=2 at x=0')

    call check_close(wave_func(-0.3_precision), 1.0_precision, test_tol, 'wave_func plateau')
    call check_close(wave_func(0.05_precision), 0.5_precision, test_tol, 'wave_func tent')
    call check_close(wave_func(0.3_precision), 0.0_precision, 0.0_precision, 'wave_func outside bumps')

    call check_close(udf_func(0.0_precision, 4), 0.0_precision, test_tol, 'udf_4 at x=0')
    call check_close(udf_func(0.5_precision, 4), 2.0_precision, test_tol, 'udf_4 at x=0.5')

    call check_close(omega_offset_func(0.5_precision), 0.25_precision * exp(-0.375_precision), &
        test_tol, 'omega offset at x=0.5')

    test_func_type = 3
    udf_id = 4
    call check_close(standard_func(0.5_precision), 2.0_precision, test_tol, 'standard_func UDF dispatch')

    test_func_type = 1
    n_order_critical_point = 1
    call check_close(standard_func(1.0_precision), 1.0_precision, test_tol, &
        'standard_func critical point dispatch')

    open(newunit=u, status='scratch', form='formatted')
    write(u, '(a)') '&test_nml'
    write(u, '(a)') '  iter_times = 3, udf_id = 5'
    write(u, '(a)') '/'
    rewind(u)
    call read_test_namelist(u)
    close(u)
    call check_true(iter_times == 3_integer_precision, 'test namelist iter_times')
    call check_true(udf_id == 5_integer_precision, 'test namelist udf_id')

    test_func_type = 99
    call validate_test_params(ierr, msg)
    call check_true(ierr /= 0, 'invalid test_func_type rejected')

    test_func_type = 3
    call validate_test_params(ierr, msg)
    call check_true(ierr == 0, 'valid test params accepted')

    call finish()

end program test_standard_function
