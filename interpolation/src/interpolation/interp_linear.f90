module interp_linear
    use precision_select
    implicit none

    ! All routines interpolate from the stencil u(-3:4) to the half-point
    ! between u(0) and u(1).  They are internal helpers; diagnostics are
    ! attached only to the scheme routines exposed by scheme_dispatch.

contains

    subroutine linear_2nd(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = 0.5_precision * (u(0) + u(1))
    end subroutine linear_2nd

    subroutine linear_3rd_minus1_to_plus1(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = u(0) + (real(3, precision) * u(1) - real(2, precision) * u(0) - u(-1)) &
               / real(8, precision)
    end subroutine linear_3rd_minus1_to_plus1

    subroutine linear_3rd_minus2_to_0(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (real(3, precision) * u(-2) - real(10, precision) * u(-1) &
               + real(15, precision) * u(0)) / real(8, precision)
    end subroutine linear_3rd_minus2_to_0

    subroutine linear_3rd_0_to_plus2(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (real(3, precision) * u(0) + real(6, precision) * u(1) - u(2)) &
               / real(8, precision)
    end subroutine linear_3rd_0_to_plus2

    subroutine linear_4th_minus2_to_plus1(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (u(-2) - real(5, precision) * u(-1) + real(15, precision) * u(0) &
               + real(5, precision) * u(1)) / real(16, precision)
    end subroutine linear_4th_minus2_to_plus1

    subroutine linear_4th_minus1_to_plus2(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (-u(-1) + real(9, precision) * u(0) + real(9, precision) * u(1) - u(2)) &
               / real(16, precision)
    end subroutine linear_4th_minus1_to_plus2

    subroutine linear_5th_minus2_to_plus2(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (real(3, precision) * u(-2) - real(20, precision) * u(-1) &
               + real(90, precision) * u(0) + real(60, precision) * u(1) &
               - real(5, precision) * u(2)) / real(128, precision)
    end subroutine linear_5th_minus2_to_plus2

    subroutine linear_6th_minus2_to_plus3(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        u_half = (real(3, precision) * u(-2) - real(25, precision) * u(-1) &
               + real(150, precision) * u(0) + real(150, precision) * u(1) &
               - real(25, precision) * u(2) + real(3, precision) * u(3)) &
               / real(256, precision)
    end subroutine linear_6th_minus2_to_plus3

    ! The 3rd-order background linear interpolation used by the nonlinear
    ! error definition for 3rd-order schemes.
    subroutine poly3(u, u_half)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half

        call linear_3rd_minus1_to_plus1(u, u_half)
    end subroutine poly3

end module interp_linear
