module io
    use precision_select
    use control, only: iter_times, ini_mesh_num, output_fields => output_fig
    use mesh
    use solver, only: func_label, x_ini, x_end, calc_delta_omega
    use standard_function, only: standard_func
    use scheme_params, only: scheme
    use scheme_dispatch, only: scheme_output_name
    use diagnostics, only: diag_point
    implicit none

contains

    subroutine output_exact_standard_func()
        integer, parameter :: n = 10000
        integer :: unit, i
        real(precision) :: x

        open(newunit=unit, file='u_exact.dat', status='replace', action='write')
        write(unit, '(a)') 'VARIABLES = x, u_exact'
        write(unit, '(a,i0,a)') 'ZONE T="u_exact", I=', n + 1, ', DATAPACKING=BLOCK'

        do i = 0, n
            x = x_ini + (x_end - x_ini) * real(i, precision) / real(n, precision)
            write(unit, '(es32.20e4)') x
        end do
        do i = 0, n
            x = x_ini + (x_end - x_ini) * real(i, precision) / real(n, precision)
            write(unit, '(es32.20e4)') standard_func(x)
        end do

        close(unit)
    end subroutine output_exact_standard_func

    subroutine output_order()
        integer :: unit, i
        character(len=64) :: method

        method = trim(scheme_output_name(scheme))

        open(newunit=unit, file=trim(method)//'_log.dat', status='replace', action='write')
        write(unit, '(a)') 'N L1_error L1_order L2_error L2_order Linf_error Linf_order ' &
            // 'En En_order En_inf En_inf_order'
        do i = 1, iter_times
            write(unit, '(i0,10(1x,es32.20e4))') 2**(i - 1) * int(ini_mesh_num), &
                L1_error(i), L1_order(i), L2_error(i), L2_order(i), &
                L_infty_error(i), L_infty_order(i), Error_nonLinear(i), &
                nonLinear_order(i), Error_nonLinear_inf(i), nonLinear_order_inf(i)
        end do
        close(unit)
    end subroutine output_order

    subroutine output_csv()
        integer :: unit, i
        character(len=64) :: method

        method = trim(scheme_output_name(scheme))

        open(newunit=unit, file='results/convergence/'//trim(method)//'_'//trim(func_label())//'.csv', &
            status='replace', action='write')
        write(unit, '(a)') 'N,L1_error,L1_order,L2_error,L2_order,Linf_error,Linf_order,' &
            // 'En,En_order,En_inf,En_inf_order'
        do i = 1, iter_times
            write(unit, '(i0,10(a,es32.20e4))') 2**(i - 1) * int(ini_mesh_num), &
                ',', L1_error(i), ',', L1_order(i), ',', L2_error(i), ',', L2_order(i), &
                ',', L_infty_error(i), ',', L_infty_order(i), ',', Error_nonLinear(i), &
                ',', nonLinear_order(i), ',', Error_nonLinear_inf(i), &
                ',', nonLinear_order_inf(i)
        end do
        close(unit)
    end subroutine output_csv

    subroutine output_fig(Mesh_points)
        integer(integer_precision), intent(in) :: Mesh_points
        integer :: unit, i, n
        character(len=16) :: nstr
        character(len=64) :: method, fname
        real(precision), allocatable :: col(:)

        if (.not. output_fields) return

        call calc_delta_omega(Mesh_points)

        n = int(Mesh_points)
        method = trim(scheme_output_name(scheme))
        write(nstr, '(i0)') n
        fname = 'rst'//trim(nstr)//trim(method)//'.dat'

        open(newunit=unit, file=trim(fname), status='replace', action='write')
        write(unit, '(a)') 'VARIABLES =index, x, u_exact, u_interp, error, ' &
            // 'beta_1, beta_2, beta_3, omega_1, omega_2, omega_3, tau, ' &
            // 'beta_L, beta_N, omega_L, omega_N, lambda_h2, h_noscale, ' &
            // 'delta_omega_mml, delta_omega_bjw'
        write(unit, '(a)') 'ZONE T="N: '//trim(nstr)//'", I='//trim(nstr)//', DATAPACKING=BLOCK'

        allocate(col(0:n - 1))

        do i = 0, n - 1
            col(i) = half_point(i)%index
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = half_point(i)%x
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = half_point(i)%u_exact
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = half_point(i)%u_interpolate
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = abs(half_point(i)%u_exact - half_point(i)%u_interpolate)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%beta(1)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%beta(2)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%beta(3)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%omega(1)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%omega(2)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%omega(3)
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%tau
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%beta_L
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%beta_N
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%omega_L
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%omega_N
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%lambda_h2
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = diag_point(i)%h_noscale
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = half_point(i)%delta_omega_mml
        end do
        call write_column(unit, n, col)

        do i = 0, n - 1
            col(i) = half_point(i)%delta_omega_bjw
        end do
        call write_column(unit, n, col)

        deallocate(col)
        close(unit)
    end subroutine output_fig

    subroutine write_column(unit, n, values)
        integer, intent(in) :: unit, n
        real(precision), intent(in) :: values(0:n - 1)
        integer :: i, j

        do i = 0, n - 1, 6
            write(unit, '(6(es32.20e4,2x))') (values(j), j=i, min(i + 5, n - 1))
        end do
    end subroutine write_column

end module io
