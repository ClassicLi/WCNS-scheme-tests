module solver
    use precision_select
    use control, only: iter_times, ini_mesh_num, test_func_type, &
        n_order_critical_point, udf_id, nonlinear_baseline
    use mesh
    use standard_function, only: standard_func
    use scheme_params, only: scheme
    use scheme_dispatch, only: interpolate_dispatch, scheme_order
    use interp_linear, only: linear_2nd, poly3, linear_5th_minus2_to_plus2, linear_6th_minus2_to_plus3
    use scheme_diagnostics_mod, only: scheme_diagnostics
    use diagnostics, only: init_diagnostics, store_diagnostics
    implicit none

    real(precision), parameter :: x_ini = -1.0_precision
    real(precision), parameter :: x_end = 1.0_precision

contains

    function get_dx(Mesh_points) result(h)
        integer(integer_precision), intent(in) :: Mesh_points
        real(precision) :: h

        h = (x_end - x_ini) / real(Mesh_points, precision)
    end function get_dx

    subroutine calc_init(Mesh_points)
        integer(integer_precision), intent(in) :: Mesh_points
        integer(integer_precision) :: i
        real(precision) :: h

        if (allocated(integer_point)) deallocate(integer_point)
        if (allocated(half_point)) deallocate(half_point)
        allocate(integer_point(-5:Mesh_points + 5))
        allocate(half_point(0:Mesh_points - 1))
        call init_diagnostics(int(Mesh_points))

        h = get_dx(Mesh_points)

        do i = -5, Mesh_points + 5
            integer_point(i)%index = real(i, precision)
            integer_point(i)%x = x_ini + real(i, precision) * h
            integer_point(i)%u_exact = standard_func(integer_point(i)%x)
        end do

        do i = 0, Mesh_points - 1
            half_point(i)%index = real(i, precision) + 0.5_precision
            half_point(i)%x = x_ini + real(i, precision) * h + 0.5_precision * h
            half_point(i)%u_exact = standard_func(half_point(i)%x)
            half_point(i)%u_linear = 0.0_precision
            half_point(i)%u_interpolate = 0.0_precision
            half_point(i)%omega_1 = 0.0_precision
            half_point(i)%omega_2 = 0.0_precision
            half_point(i)%omega_3 = 0.0_precision
            half_point(i)%delta_omega_mml = 0.0_precision
            half_point(i)%delta_omega_bjw = 0.0_precision
        end do
    end subroutine calc_init

    subroutine interpolate(Mesh_points)
        integer(integer_precision), intent(in) :: Mesh_points
        integer(integer_precision) :: i
        type(scheme_diagnostics) :: diag

        do i = 0, Mesh_points - 1
            call interpolate_dispatch(integer_point(i - 3:i + 4)%u_exact, &
                half_point(i)%u_interpolate, diag)
            half_point(i)%omega_1 = diag%omega(1)
            half_point(i)%omega_2 = diag%omega(2)
            half_point(i)%omega_3 = diag%omega(3)
            call store_diagnostics(int(i), diag)
        end do
    end subroutine interpolate

    subroutine calc_error(Mesh_points, iter)
        integer(integer_precision), intent(in) :: Mesh_points
        integer(integer_precision), intent(in) :: iter
        integer(integer_precision) :: i
        integer :: p
        real(precision) :: u_bg

        p = scheme_order(scheme)

        do i = 0, Mesh_points - 1
            if (is_linear_scheme()) then
                half_point(i)%u_linear = half_point(i)%u_interpolate
            else if (trim(nonlinear_baseline) == 'poly3') then
                call poly3(integer_point(i - 3:i + 4)%u_exact, half_point(i)%u_linear)
            else
                select case (p)
                case (2)
                    call linear_2nd(integer_point(i - 3:i + 4)%u_exact, half_point(i)%u_linear)
                case (3)
                    call poly3(integer_point(i - 3:i + 4)%u_exact, half_point(i)%u_linear)
                case (5)
                    call linear_5th_minus2_to_plus2(integer_point(i - 3:i + 4)%u_exact, &
                        half_point(i)%u_linear)
                case (6)
                    call linear_6th_minus2_to_plus3(integer_point(i - 3:i + 4)%u_exact, &
                        half_point(i)%u_linear)
                case default
                    half_point(i)%u_linear = half_point(i)%u_interpolate
                end select
            end if
        end do

        do i = 0, Mesh_points - 1
            u_bg = half_point(i)%u_linear
            L1_error(iter) = L1_error(iter) + abs(half_point(i)%u_exact - half_point(i)%u_interpolate)
            L2_error(iter) = L2_error(iter) &
                + (half_point(i)%u_exact - half_point(i)%u_interpolate)**2
            Error_nonLinear(iter) = Error_nonLinear(iter) &
                + abs(u_bg - half_point(i)%u_interpolate)
        end do

        L1_error(iter) = L1_error(iter) / real(Mesh_points, precision)
        L2_error(iter) = sqrt(L2_error(iter) / real(Mesh_points, precision))
        Error_nonLinear(iter) = Error_nonLinear(iter) / real(Mesh_points, precision)

        do i = 0, Mesh_points - 1
            u_bg = half_point(i)%u_linear
            L_infty_error(iter) = max(L_infty_error(iter), &
                abs(half_point(i)%u_exact - half_point(i)%u_interpolate))
            Error_nonLinear_inf(iter) = max(Error_nonLinear_inf(iter), &
                abs(u_bg - half_point(i)%u_interpolate))
        end do
    end subroutine calc_error

    subroutine calc_order()
        integer(integer_precision) :: i
        real(precision) :: ln2

        ln2 = log(2.0_precision)

        do i = 2, iter_times
            L1_order(i) = safe_order(L1_error(i - 1), L1_error(i), ln2)
            L2_order(i) = safe_order(L2_error(i - 1), L2_error(i), ln2)
            L_infty_order(i) = safe_order(L_infty_error(i - 1), L_infty_error(i), ln2)
            nonLinear_order(i) = safe_order(Error_nonLinear(i - 1), Error_nonLinear(i), ln2)
            nonLinear_order_inf(i) = safe_order(Error_nonLinear_inf(i - 1), &
                Error_nonLinear_inf(i), ln2)
        end do
    end subroutine calc_order

    subroutine procedure()
        integer(integer_precision) :: i
        integer(integer_precision) :: mesh_points

        call allocate_errors(int(iter_times))
        time_used = 0

        call system_clock(time_ini)

        do i = 1, iter_times
            mesh_points = 2**(i - 1) * ini_mesh_num
            call calc_init(mesh_points)
            call interpolate(mesh_points)
            call calc_error(mesh_points, i)
            call calc_delta_omega(mesh_points)
            deallocate(integer_point)
            deallocate(half_point)
        end do

        call system_clock(time_end)
        time_used = time_end - time_ini

        call calc_order()
    end subroutine procedure

    subroutine calc_delta_omega(Mesh_points)
        integer(integer_precision), intent(in) :: Mesh_points
        integer(integer_precision) :: i
        integer :: nw
        real(precision) :: w1, w2, w3

        nw = n_weights()

        if (nw == 2) then
            do i = 0, Mesh_points - 1
                w1 = half_point(i)%omega_1
                w2 = half_point(i)%omega_2
                half_point(i)%delta_omega_mml = abs(w1 - 0.25_precision) + abs(w2 - 0.75_precision)
                half_point(i)%delta_omega_bjw = sqrt((w1 - 0.25_precision)**2 &
                    + (w2 - 0.75_precision)**2)
            end do
        else if (nw == 3) then
            do i = 0, Mesh_points - 1
                w1 = half_point(i)%omega_1
                w2 = half_point(i)%omega_2
                w3 = half_point(i)%omega_3
                half_point(i)%delta_omega_mml = abs(w1 - 0.0625_precision) &
                    + abs(w2 - 0.625_precision) + abs(w3 - 0.3125_precision)
                half_point(i)%delta_omega_bjw = sqrt((w1 - 0.0625_precision)**2 &
                    + (w2 - 0.625_precision)**2 + (w3 - 0.3125_precision)**2)
            end do
        else
            do i = 0, Mesh_points - 1
                half_point(i)%delta_omega_mml = 0.0_precision
                half_point(i)%delta_omega_bjw = 0.0_precision
            end do
        end if
    end subroutine calc_delta_omega

    function func_label() result(s)
        character(len=32) :: s

        select case (test_func_type)
        case (1_integer_precision)
            write(s, '(a,i0)') 'cp', n_order_critical_point
        case (2_integer_precision)
            s = 'wave'
        case (3_integer_precision)
            write(s, '(a,i0)') 'udf', udf_id
        case (4_integer_precision)
            s = 'omega'
        case default
            s = 'unknown'
        end select
    end function func_label

    logical function is_linear_scheme()
        is_linear_scheme = index(trim(scheme), 'LINEAR') == 1
    end function is_linear_scheme

    integer function n_weights()
        select case (trim(scheme))
        case ('JS3', 'Z3', 'MOC322', 'MOC33', 'MOC53')
            n_weights = 2
        case ('JS5', 'Z5', 'MOZ')
            n_weights = 3
        case default
            n_weights = 0
        end select
    end function n_weights

    function safe_order(e1, e2, ln2) result(p)
        real(precision), intent(in) :: e1, e2, ln2
        real(precision) :: p

        if (e1 > 0.0_precision .and. e2 > 0.0_precision) then
            p = log(e1 / e2) / ln2
        else
            p = 0.0_precision
        end if
    end function safe_order

end module solver
