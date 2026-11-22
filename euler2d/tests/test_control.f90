program test_control
    ! Unit tests for 2D Euler namelist parsing, case defaults, and operator pairing.
    use precision_select
    use e2d_globals, only: gama
    use e2d_control
    use test_util
    implicit none

    integer :: ierr
    character(len=256) :: msg
    character(len=32) :: op

    ! ---------------------------------------------------------------
    ! Defaults and case metadata
    ! ---------------------------------------------------------------
    call write_nml('riemann1', 0, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'riemann1 accepted')
    call check_true(nz == 400 .and. nx == 400, 'riemann1 defaults 400x400')
    call check_close(tend, 0.8_precision, 1.0e-14_precision, 'riemann1 tend')
    call check_close(gama, 1.4_precision, 1.0e-14_precision, 'riemann1 gamma')

    call write_nml('shock_bubble', 0, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'shock_bubble accepted')
    call check_true(nz == 800 .and. nx == 1600, 'shock_bubble 800x1600')
    call check_close(tend2, 0.2_precision, 1.0e-14_precision, 'shock_bubble tend2')

    call write_nml('rayleigh_taylor', 0, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'rt accepted')
    call check_true(nz == 400 .and. nx == 100, 'rt 400x100')
    call check_close(gama, 5.0_precision/3.0_precision, 1.0e-14_precision, 'rt gamma')

    ! ---------------------------------------------------------------
    ! Invalid configuration
    ! ---------------------------------------------------------------
    call write_nml('sedov', 64, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'unknown case rejected')

    call write_nml('riemann1', 7, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'mesh below 8 rejected')

    call write_nml('rayleigh_taylor', 10, 0.3_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'rt mesh not divisible by 4 rejected')

    call write_nml('riemann1', 64, 0.0_precision, -1.0_precision, -1.0_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'zero cfl rejected')

    call write_nml('riemann1', 64, 0.3_precision, -1.0_precision, 0.1_precision, &
                   0.0_precision, 'QR', 'AUTO')
    call e2d_control_init('e2d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'tend2 on non-shock_bubble rejected')

    ! ---------------------------------------------------------------
    ! AUTO operator pairing
    ! ---------------------------------------------------------------
    difference_operator = 'AUTO'
    call e2d_resolve_operator('MUSCL', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'WCNSE2', 'AUTO MUSCL -> WCNSE2')
    call e2d_resolve_operator('QR', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS4', 'AUTO QR -> HCS4')
    call e2d_resolve_operator('LINEAR4_M2P1', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'WCNSE4', 'AUTO 4th -> WCNSE4')
    call e2d_resolve_operator('Z5', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'AUTO Z5 -> HCS6')
    call e2d_resolve_operator('LINEAR6_M2P3', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'WCNSE6', 'AUTO 6th -> WCNSE6')

    ! Explicit mismatch warns but continues.
    difference_operator = 'HCS6'
    call e2d_resolve_operator('QR', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'explicit mismatch continues')

    difference_operator = 'E6'
    call e2d_resolve_operator('QR', op, ierr, msg)
    call check_true(ierr /= 0, 'legacy E6 name rejected')

    call delete_nml()
    call finish()

contains

    subroutine write_nml(case_name, mesh, cfl_in, tend_in, tend2_in, gamma_in, &
                         scheme_in, operator)
        ! Write a complete input.nml for one test case.
        character(len=*), intent(in) :: case_name
        integer, intent(in) :: mesh
        real(precision), intent(in) :: cfl_in, tend_in, tend2_in, gamma_in
        character(len=*), intent(in) :: scheme_in, operator
        integer :: u

        open (newunit=u, file='e2d_params_test.nml', status='replace', action='write')
        write (u, '(a)') '&euler2d_nml'
        write (u, '(a)') "  case_2d = '"//trim(case_name)//"',"
        write (u, '(a,i0)') '  mesh_num = ', mesh
        write (u, '(a,es24.16)') '  cfl = ', cfl_in
        write (u, '(a,es24.16)') '  tend = ', tend_in
        write (u, '(a,es24.16)') '  tend2 = ', tend2_in
        write (u, '(a,es24.16)') '  gamma = ', gamma_in
        write (u, '(a)') "  output_dir = 'results',"
        write (u, '(a)') '  write_sample = .false.'
        write (u, '(a)') '/'
        write (u, '(a)') '&scheme_nml'
        write (u, '(a)') "  scheme = '"//trim(scheme_in)//"'"
        write (u, '(a)') '/'
        write (u, '(a)') '&difference_nml'
        write (u, '(a)') "  difference_operator = '"//trim(operator)//"'"
        write (u, '(a)') '/'
        close (u)
    end subroutine write_nml

    subroutine delete_nml()
        ! Remove the temporary namelist file.
        integer :: u

        open (newunit=u, file='e2d_params_test.nml', status='old')
        close (u, status='delete')
    end subroutine delete_nml

end program test_control
