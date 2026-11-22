module adr_dft
    ! Discrete Fourier transform and modified-wavenumber extraction.
    use precision_select
    use adr_control, only: PI, II
    implicit none

contains

    subroutine dft_mode(u, n, n_num, u_hat)
        ! Compute the n-th Fourier coefficient of u(0:n_num-1).
        integer(integer_precision), intent(in) :: n, n_num
        real(precision), intent(in) :: u(0:n_num)
        complex(precision), intent(out) :: u_hat

        real(precision) :: phyn
        integer(integer_precision) :: j

        phyn = 2.0_precision * PI * real(n, precision) / real(n_num, precision)
        u_hat = cmplx(0.0_precision, 0.0_precision, precision)

        do j = 0, n_num - 1
            u_hat = u_hat + u(j) * exp(-II * real(j, precision) * phyn)
        end do

        u_hat = u_hat / real(n_num, precision)
    end subroutine dft_mode

    subroutine calc_phi(hat_t, hat_0, dt, n_num, phi)
        ! phi = 2*pi*i*log(hat_t/hat_0)/(dt*n_num).
        complex(precision), intent(in) :: hat_t, hat_0
        real(precision), intent(in) :: dt
        integer(integer_precision), intent(in) :: n_num
        complex(precision), intent(out) :: phi

        phi = 2.0_precision * PI * II * log(hat_t / hat_0) / dt / real(n_num, precision)
    end subroutine calc_phi

end module adr_dft
