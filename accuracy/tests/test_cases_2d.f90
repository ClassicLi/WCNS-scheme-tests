program test_cases_2d
    ! Unit tests for the 2D isentropic-vortex and sine-density cases.
    use precision_select
    use e2d_globals
    use e2d_cases
    use test_util
    implicit none

    integer :: ierr, j, k
    character(len=256) :: msg
    real(precision) :: rho(0:10, 0:10), length

    ! ---------------------------------------------------------------
    ! Sine density wave on [0, 2*pi]^2
    ! ---------------------------------------------------------------
    call e2d_allocate(8, 2.0_precision*acos(-1.0_precision))
    call e2d_initial('sine', 1.0_precision, 0.2_precision, ierr, msg)
    call check_true(ierr == 0, 'sine initial accepted')

    call e2d_exact_density('sine', 0.0_precision, 1.0_precision, 0.2_precision, &
                           rho(0:8, 0:8), ierr, msg)
    do k = 0, 8
        do j = 0, 8
            call check_close(rho(j, k), ubar(j, k, 1), 1.0e-12_precision, &
                             'sine exact t=0')
        end do
    end do

    call e2d_exact_density('sine', 0.5_precision, 1.0_precision, 0.2_precision, &
                           rho(0:8, 0:8), ierr, msg)
    call check_close(rho(2, 3), 1.0_precision + 0.2_precision* &
                     sin(2.0_precision*h + 3.0_precision*h - 1.0_precision), &
                     1.0e-12_precision, 'sine phase')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! Isentropic vortex on [0, 10]^2
    ! ---------------------------------------------------------------
    call e2d_allocate(10, 10.0_precision)
    call e2d_initial('vortex', 0.0_precision, 0.0_precision, ierr, msg)
    call check_true(ierr == 0, 'vortex initial accepted')

    call e2d_exact_density('vortex', 0.0_precision, 0.0_precision, 0.0_precision, &
                           rho, ierr, msg)
    call check_close(rho(4, 4), ubar(4, 4, 1), 1.0e-12_precision, 'vortex exact t=0')

    call e2d_exact_density('vortex', 1.0_precision, 0.0_precision, 0.0_precision, &
                           rho, ierr, msg)
    call check_close(rho(4, 4), ubar(3, 4, 1), 1.0e-12_precision, &
                     'vortex advects one cell')

    call e2d_case_length('nope', length, ierr, msg)
    call check_true(ierr /= 0, 'unknown case rejected')
    call e2d_deallocate()

    call finish()
end program test_cases_2d
