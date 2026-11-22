module cw_initial
    ! Complex initial condition used by the reference complex-wave test.
    !
    ! The profile combines four features on [-1, 1]:
    !   - a smoothed Gaussian bump centered at x = -0.7,
    !   - a square block of unit height on [-0.4, -0.2],
    !   - a triangular wave with apex at x = 0.1,
    !   - a smoothed semicircle centered at x = 0.5.
    ! The Gaussian and semicircle are smoothed with the three-point stencil
    ! (f(x-delta) + 4 f(x) + f(x+delta)) / 6 using delta = 0.005, exactly as
    ! in the reference program (at mesh 200 this is delta = h/2).
    use precision_select
    implicit none

    ! Half-width of the smoothing stencil, fixed as in the reference program.
    real(precision), parameter :: init_delta = 0.005_precision

contains

    subroutine initial_complex_wave(nx, u)
        ! Fill the interior u(0:nx) with the four-feature profile.
        integer(integer_precision), intent(in) :: nx
        real(precision), intent(out) :: u(0:nx)

        real(precision), parameter :: z_gauss = -0.7_precision
        real(precision), parameter :: a_circle = 0.5_precision
        real(precision) :: h, x
        integer(integer_precision) :: k

        h = 2.0_precision / real(nx, precision)
        u = 0.0_precision

        ! Smooth Gaussian bump on the first ten-percent segment.
        do k = nx/10, 2*nx/10
            x = -1.0_precision + real(k, precision)*h
            u(k) = (gaussian_profile(x, z_gauss - init_delta) &
                  + 4.0_precision*gaussian_profile(x, z_gauss) &
                  + gaussian_profile(x, z_gauss + init_delta))/6.0_precision
        end do

        ! Square block on the third ten-percent segment.
        u(3*nx/10:4*nx/10) = 1.0_precision

        ! Triangular wave on the fifth ten-percent segment.  The max(0, .)
        ! only matters when nx is not a multiple of 10: the segment endpoints
        ! can then fall slightly outside x in [0, 0.2], where the reference
        ! formula would give a small negative value.
        do k = 5*nx/10, 6*nx/10
            x = -1.0_precision + real(k, precision)*h
            u(k) = max(0.0_precision, 1.0_precision - abs(10.0_precision*x - 1.0_precision))
        end do

        ! Smoothed semicircle on the seventh ten-percent segment.
        do k = 7*nx/10, 8*nx/10
            x = -1.0_precision + real(k, precision)*h
            u(k) = (semicircle_profile(x, a_circle - init_delta) &
                  + 4.0_precision*semicircle_profile(x, a_circle) &
                  + semicircle_profile(x, a_circle + init_delta))/6.0_precision
        end do
    end subroutine initial_complex_wave

    real(precision) function gaussian_profile(x, center) result(value)
        ! Gaussian bump with a half-height radius of 0.03.
        real(precision), intent(in) :: x, center
        real(precision), parameter :: beta = log(2.0_precision)/(0.03_precision**2)

        value = exp(-beta*(x - center)**2)
    end function gaussian_profile

    real(precision) function semicircle_profile(x, center) result(value)
        ! Unit semicircle of half-width 0.1; zero outside its support.
        real(precision), intent(in) :: x, center
        real(precision), parameter :: alfa = 10.0_precision

        value = 1.0_precision - (alfa*(x - center))**2
        if (value >= 0.0_precision) then
            value = sqrt(value)
        else
            value = 0.0_precision
        end if
    end function semicircle_profile

end module cw_initial
