module e1d_globals
    ! Shared constants for the 1D Euler test cases.
    !
    ! The halo is wide enough for the 7-point interpolation stencil and for
    ! the farthest integer-node read of the hybrid difference operators.
    use precision_select
    implicit none

    integer(integer_precision), parameter :: halo = 7_integer_precision
    real(precision), parameter :: gama = 1.4_precision
    real(precision), parameter :: beta = gama - 1.0_precision

end module e1d_globals
