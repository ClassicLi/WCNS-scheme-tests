module standard_function
    use precision_select
    use control
    implicit none

    real(precision), parameter :: PI = acos(-1.0_precision)

contains

    function critical_point_func(x, n) result(f)
        real(precision), intent(in) :: x
        integer, intent(in) :: n
        real(precision) :: f

        f = exp(0.75_precision * (x - 1.0_precision)) * x**(n + 1)
    end function critical_point_func

    function gaussian(x, beta, z) result(val)
        real(precision), intent(in) :: x, beta, z
        real(precision) :: val

        val = exp(-beta * (x - z)**2)
    end function gaussian

    function bump(x, alpha, a) result(val)
        real(precision), intent(in) :: x, alpha, a
        real(precision) :: val

        val = sqrt(max(0.0_precision, 1.0_precision - alpha * alpha * (x - a)**2))
    end function bump

    function wave_func(x) result(f)
        real(precision), intent(in) :: x
        real(precision) :: f
        real(precision) :: beta
        real(precision), parameter :: z = -0.7_precision
        real(precision), parameter :: delta = 0.005_precision
        real(precision), parameter :: a = 0.5_precision
        real(precision), parameter :: alpha = 10.0_precision

        beta = log(2.0_precision) / (36.0_precision * delta**2)

        if (x >= -0.8_precision .and. x <= -0.6_precision) then
            f = (gaussian(x, beta, z - delta) + 4.0_precision * gaussian(x, beta, z) &
               + gaussian(x, beta, z + delta)) / 6.0_precision
        else if (x >= -0.4_precision .and. x <= -0.2_precision) then
            f = 1.0_precision
        else if (x >= 0.0_precision .and. x <= 0.2_precision) then
            f = 1.0_precision - abs(10.0_precision * (x - 0.1_precision))
        else if (x >= 0.4_precision .and. x <= 0.6_precision) then
            f = (bump(x, alpha, a - delta) + 4.0_precision * bump(x, alpha, a) &
               + bump(x, alpha, a + delta)) / 6.0_precision
        else
            f = 0.0_precision
        end if
    end function wave_func

    function udf_func(x, id) result(f)
        real(precision), intent(in) :: x
        integer, intent(in) :: id
        real(precision) :: f

        select case (id)
        case (1)
            if (x <= 0.0_precision) then
                f = real(1, precision) - (x + 0.5_precision)**2
            else
                f = -real(1, precision) + (x - 0.5_precision)**2
            end if
        case (2)
            if (x <= 0.0_precision) then
                f = -sin(PI * x)
            else
                f = -sin(PI * x) + 1.0_precision
            end if
        case (3)
            f = cos(x)
        case (4)
            if (x <= 0.375_precision) then
                f = x * x * x * exp(0.75_precision * (x - 1.0_precision))
            else
                f = 1.0_precision + sin(PI * x)
            end if
        case (5)
            if (x >= -0.25_precision .and. x <= 0.25_precision) then
                f = x * x * x * exp(0.75_precision * (x - 1.0_precision))
            else if (x > 0.25_precision) then
                f = 1.0_precision + sin(PI * x)
            else
                f = sin(PI * (x + 0.5_precision))**4 + 0.1_precision
            end if
        case (6)
            if (x <= 0.0_precision) then
                f = -sin(PI * x) - 0.5_precision * x * x * x
            else
                f = -sin(PI * x) - 0.5_precision * x * x * x + 1.0_precision
            end if
        case default
            f = 0.0_precision
        end select
    end function udf_func

    function omega_offset_func(x) result(f)
        real(precision), intent(in) :: x
        real(precision) :: f

        if (x <= 0.525_precision) then
            f = x * x * exp(0.75_precision * (x - 1.0_precision))
        else
            f = 2.0_precision + x * x * exp(0.75_precision * (x - 1.0_precision))
        end if
    end function omega_offset_func

    function standard_func(x) result(f)
        real(precision), intent(in) :: x
        real(precision) :: f

        select case (test_func_type)
        case (1_integer_precision)
            f = critical_point_func(x, int(n_order_critical_point))
        case (2_integer_precision)
            f = wave_func(x)
        case (3_integer_precision)
            f = udf_func(x, int(udf_id))
        case (4_integer_precision)
            f = omega_offset_func(x)
        case default
            f = 0.0_precision
        end select
    end function standard_func

end module standard_function
