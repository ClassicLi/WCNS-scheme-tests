program test_io
    use precision_select
    use control, only: ini_mesh_num, iter_times, test_func_type, udf_id, &
        n_order_critical_point, nonlinear_baseline, output_fig_flag => output_fig
    use mesh
    use solver
    use scheme_params, only: scheme
    use diagnostics, only: diag_point
    use io
    use test_util
    implicit none

    ini_mesh_num = 10_integer_precision
    iter_times = 2_integer_precision
    test_func_type = 3_integer_precision
    udf_id = 3_integer_precision
    output_fig_flag = .true.
    nonlinear_baseline = 'same_order'

    scheme = 'LINEAR3_M1P1'
    call procedure()
    call output_order()
    call output_csv()

    call calc_init(10_integer_precision)
    call interpolate(10_integer_precision)
    call output_fig(10_integer_precision)

    call check_true(file_exists('LINEAR3_M1P1_log.dat'), 'log file exists')
    call check_true(file_contains('LINEAR3_M1P1_log.dat', 'L1_error'), 'log header contains L1_error')
    call check_true(count_lines('LINEAR3_M1P1_log.dat') == 3, 'log has header plus two rows')

    call check_close(half_point(0)%delta_omega_mml, 0.0_precision, 0.0_precision, &
        'linear scheme delta_omega zero')
    call check_true(file_exists('rst10LINEAR3_M1P1.dat'), 'field file exists')
    call check_true(file_contains('rst10LINEAR3_M1P1.dat', 'beta_L'), 'field has beta_L')
    call check_true(file_contains('rst10LINEAR3_M1P1.dat', 'omega_L'), 'field has omega_L')
    call check_true(file_contains('rst10LINEAR3_M1P1.dat', 'lambda_h2'), 'field has lambda_h2')
    call check_true(.not. file_contains('rst10LINEAR3_M1P1.dat', 'A_1'), 'field drops A_1')

    call check_true(file_exists('results/convergence/LINEAR3_M1P1_udf3.csv'), 'CSV exists')
    call check_true(count_lines('results/convergence/LINEAR3_M1P1_udf3.csv') == 3, 'CSV has two rows')

    scheme = 'Z5'
    call calc_init(10_integer_precision)
    call interpolate(10_integer_precision)
    call output_fig(10_integer_precision)
    call check_true(file_exists('rst10Z5.dat'), 'Z5 field file exists')
    call check_true(diag_point(0)%beta(1) > 0.0_precision, 'Z5 beta diagnostics nonzero')
    call check_true(diag_point(0)%omega(1) > 0.0_precision, 'Z5 omega diagnostics nonzero')
    call check_true(half_point(0)%delta_omega_mml > 0.0_precision, 'Z5 delta_omega mml nonzero')
    call check_true(half_point(0)%delta_omega_bjw > 0.0_precision, 'Z5 delta_omega bjw nonzero')

    call finish()

contains

    logical function file_exists(name)
        character(len=*), intent(in) :: name
        logical :: exists

        inquire(file=trim(name), exist=exists)
        file_exists = exists
    end function file_exists

    integer function count_lines(name)
        character(len=*), intent(in) :: name
        integer :: unit, ios, n
        character(len=512) :: line

        count_lines = 0
        open(newunit=unit, file=trim(name), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        n = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            n = n + 1
        end do
        close(unit)
        count_lines = n
    end function count_lines

    logical function file_contains(name, token)
        character(len=*), intent(in) :: name, token
        integer :: unit, ios
        character(len=512) :: line
        logical :: found

        found = .false.
        open(newunit=unit, file=trim(name), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            file_contains = .false.
            return
        end if
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, trim(token)) > 0) then
                found = .true.
                exit
            end if
        end do
        close(unit)
        file_contains = found
    end function file_contains

end program test_io
