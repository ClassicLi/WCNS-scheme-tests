program test_periodic
    ! Unit test for the periodic halo refresh.
    use precision_select
    use cw_periodic, only: boundary_cycle
    use test_util
    implicit none

    integer(integer_precision), parameter :: nx = 8, halo = 7
    real(precision) :: u(-halo:nx+halo)
    integer(integer_precision) :: i, j

    ! The grid stores nx + 1 points, u(0) through u(nx), where u(nx)
    ! duplicates u(0) to close the periodic interval.
    do j = 0, nx
        u(j) = real(j, precision)
    end do

    call boundary_cycle(nx, u, halo)

    do i = 1, halo
        call check_close(u(-i), real(nx - i, precision), 0.0_precision, 'left halo')
        call check_close(u(nx + i), real(i, precision), 0.0_precision, 'right halo')
    end do

    call finish()
end program test_periodic
