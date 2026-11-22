module c1d_exact_mod
    ! Exact solution for the 1D periodic convection test.
    use precision_select
    use c1d_initial_mod, only: c1d_profile
    implicit none

contains

    subroutine c1d_exact(n, initial_name, time, u)
        ! u(x, t) = u0(x - t), wrapped periodically into [-1, 1].
        integer, intent(in) :: n
        character(len=*), intent(in) :: initial_name
        real(precision), intent(in) :: time
        real(precision), intent(out) :: u(0:n)

        real(precision) :: h, x, xw
        integer :: k

        h = 2.0_precision/real(n, precision)
        do k = 0, n
            x = -1.0_precision + real(k, precision)*h - time
            xw = modulo(x + 1.0_precision, 2.0_precision) - 1.0_precision
            u(k) = c1d_profile(xw, initial_name)
        end do
    end subroutine c1d_exact

end module c1d_exact_mod
