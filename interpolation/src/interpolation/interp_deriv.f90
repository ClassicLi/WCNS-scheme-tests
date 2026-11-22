module interp_deriv
    use precision_select
    implicit none

    ! Derivative approximations used by the RBF shape-parameter formulas.
    ! All routines assume unit grid spacing and use the stencil u(-3:4).

contains

    ! First derivative, central, second-order truncation.
    function get_derivative_1st_tr_2(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = 0.5_precision * (u(1) - u(-1))
    end function get_derivative_1st_tr_2

    ! First derivative, left-biased, first-order truncation.
    function get_derivative_1st_tr_1_left(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = u(0) - u(-1)
    end function get_derivative_1st_tr_1_left

    ! First derivative, right-biased, first-order truncation.
    function get_derivative_1st_tr_1_right(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = u(1) - u(0)
    end function get_derivative_1st_tr_1_right

    ! Third derivative, central, second-order truncation.
    function get_derivative_3rd_tr_2(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = 0.5_precision * (-u(-2) + real(2, precision) * u(-1) &
            - real(2, precision) * u(1) + u(2))
    end function get_derivative_3rd_tr_2

    ! Third derivative, left-biased, first-order truncation.
    function get_derivative_3rd_tr_1_left(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = -u(-2) + real(3, precision) * u(-1) - real(3, precision) * u(0) + u(1)
    end function get_derivative_3rd_tr_1_left

    ! Third derivative, right-biased, first-order truncation.
    function get_derivative_3rd_tr_1_right(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d

        d = -u(-1) + real(3, precision) * u(0) - real(3, precision) * u(1) + u(2)
    end function get_derivative_3rd_tr_1_right

    ! Fifth derivative, central, second-order truncation.
    function get_derivative_5th_tr_2(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d
        real(precision) :: coe(-3:3)
        integer :: k

        coe = (/ -1.0_precision, 4.0_precision, -5.0_precision, 0.0_precision, &
                  5.0_precision, -4.0_precision, 1.0_precision /)
        d = 0.0_precision
        do k = -3, 3
            d = d + coe(k) * u(k)
        end do
        d = 0.5_precision * d
    end function get_derivative_5th_tr_2

    ! Fifth derivative, left-biased, first-order truncation.
    function get_derivative_5th_tr_1_left(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d
        real(precision) :: coe(-3:2)
        integer :: k

        coe = (/ -1.0_precision, 5.0_precision, -10.0_precision, &
                 10.0_precision, -5.0_precision, 1.0_precision /)
        d = 0.0_precision
        do k = -3, 2
            d = d + coe(k) * u(k)
        end do
    end function get_derivative_5th_tr_1_left

    ! Fifth derivative, right-biased, first-order truncation.
    function get_derivative_5th_tr_1_right(u) result(d)
        real(precision), intent(in) :: u(-3:4)
        real(precision) :: d
        real(precision) :: coe(-2:3)
        integer :: k

        coe = (/ -1.0_precision, 5.0_precision, -10.0_precision, &
                 10.0_precision, -5.0_precision, 1.0_precision /)
        d = 0.0_precision
        do k = -2, 3
            d = d + coe(k) * u(k)
        end do
    end function get_derivative_5th_tr_1_right

end module interp_deriv
