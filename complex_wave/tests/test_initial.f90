program test_initial
    ! Unit test for the four-segment complex-wave initial condition.
    use precision_select
    use cw_initial, only: initial_complex_wave
    use test_util
    implicit none

    integer(integer_precision), parameter :: nx = 200
    real(precision) :: u(0:nx), u201(0:201)
    integer(integer_precision) :: k

    call initial_complex_wave(nx, u)

    ! The square block sets the global maximum to exactly one.
    call check_close(maxval(u(0:nx)), 1.0_precision, 1.0e-12_precision, 'max amplitude')

    ! Four ten-percent segments separated by zero gaps.
    do k = 0, 19
        call check_close(u(k), 0.0_precision, 0.0_precision, 'left zero')
    end do
    do k = 41, 59
        call check_close(u(k), 0.0_precision, 0.0_precision, 'gap zero')
    end do
    do k = 81, 99
        call check_close(u(k), 0.0_precision, 0.0_precision, 'gap zero')
    end do
    do k = 121, 139
        call check_close(u(k), 0.0_precision, 0.0_precision, 'gap zero')
    end do

    call check_true(all(abs(u(60:80) - 1.0_precision) <= 0.0_precision), 'square block')
    call check_close(u(30), 0.9936_precision, 1.0e-3_precision, 'gaussian peak')
    call check_close(u(110), 1.0_precision, 1.0e-12_precision, 'triangle apex')
    call check_close(u(150), 0.9996_precision, 1.0e-3_precision, 'semicircle peak')
    call check_true(maxloc(u(20:40), dim=1) == 11, 'gaussian centered at k=30')
    call check_true(maxloc(u(140:160), dim=1) == 11, 'semicircle centered at k=150')

    ! Non-multiple-of-10 mesh: segments must stay in bounds and not overlap.
    call initial_complex_wave(201_integer_precision, u201)
    call check_true(all(u201 >= -1.0e-12_precision), 'mesh=201 lower bound')
    call check_true(all(u201 <= 1.0_precision + 1.0e-12_precision), 'mesh=201 upper bound')
    call check_close(maxval(u201), 1.0_precision, 1.0e-12_precision, 'mesh=201 max amplitude')

    call finish()
end program test_initial
