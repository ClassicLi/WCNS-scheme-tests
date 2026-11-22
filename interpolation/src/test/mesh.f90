module mesh
    use precision_select
    implicit none

    type :: points
        real(precision) :: index
        real(precision) :: x
        real(precision) :: u_exact
    end type points

    type :: half_points
        real(precision) :: index
        real(precision) :: x
        real(precision) :: u_exact
        real(precision) :: u_linear
        real(precision) :: u_interpolate
        real(precision) :: omega_1
        real(precision) :: omega_2
        real(precision) :: omega_3
        real(precision) :: delta_omega_mml
        real(precision) :: delta_omega_bjw
    end type half_points

    type(points), allocatable :: integer_point(:)
    type(half_points), allocatable :: half_point(:)

    real(precision), allocatable :: L1_error(:)
    real(precision), allocatable :: L2_error(:)
    real(precision), allocatable :: L_infty_error(:)
    real(precision), allocatable :: L1_order(:)
    real(precision), allocatable :: L2_order(:)
    real(precision), allocatable :: L_infty_order(:)
    real(precision), allocatable :: Error_nonLinear(:)
    real(precision), allocatable :: nonLinear_order(:)
    real(precision), allocatable :: Error_nonLinear_inf(:)
    real(precision), allocatable :: nonLinear_order_inf(:)

    integer(integer_precision) :: time_ini
    integer(integer_precision) :: time_end
    integer(integer_precision) :: time_used

contains

    subroutine allocate_errors(n)
        integer, intent(in) :: n

        if (allocated(L1_error)) then
            deallocate(L1_error, L2_error, L_infty_error)
            deallocate(L1_order, L2_order, L_infty_order)
            deallocate(Error_nonLinear, nonLinear_order)
            deallocate(Error_nonLinear_inf, nonLinear_order_inf)
        end if

        allocate(L1_error(n), L2_error(n), L_infty_error(n))
        allocate(L1_order(n), L2_order(n), L_infty_order(n))
        allocate(Error_nonLinear(n), nonLinear_order(n))
        allocate(Error_nonLinear_inf(n), nonLinear_order_inf(n))

        L1_error = 0.0_precision
        L2_error = 0.0_precision
        L_infty_error = 0.0_precision
        L1_order = 0.0_precision
        L2_order = 0.0_precision
        L_infty_order = 0.0_precision
        Error_nonLinear = 0.0_precision
        nonLinear_order = 0.0_precision
        Error_nonLinear_inf = 0.0_precision
        nonLinear_order_inf = 0.0_precision
    end subroutine allocate_errors

end module mesh
