module precision_select
    implicit none

    ! Single precision: about 6 digits, exponent range 10^-37 to 10^37-1
    integer, parameter :: sp = selected_real_kind(6, 37)
    ! Double precision: about 15 digits, exponent range 10^-307 to 10^307-1
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! Quadruple precision: about 33 digits, exponent range 10^-4931 to 10^4931-1
    integer, parameter :: qp = selected_real_kind(33, 4931)

    integer, parameter :: i1 = selected_int_kind(2)
    integer, parameter :: i2 = selected_int_kind(4)
    integer, parameter :: i4 = selected_int_kind(9)
    integer, parameter :: i8 = selected_int_kind(18)

    ! Default quad precision; build with PRECISION=dp (-DPRECISION_DP) to select double.
#ifdef PRECISION_DP
    integer, parameter :: precision = dp
#else
    integer, parameter :: precision = qp
#endif
    integer, parameter :: integer_precision = i8

end module precision_select
