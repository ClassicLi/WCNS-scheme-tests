program test_control
    ! Unit tests for euler1d namelist parsing, validation, and operator pairing.
    use precision_select
    use e1d_control
    use test_util
    implicit none

    integer :: ierr
    character(len=256) :: msg
    character(len=32) :: op

    ! ---------------------------------------------------------------
    ! Solver configuration
    ! ---------------------------------------------------------------
    call write_nml('sod', 64, 0.2_precision, -1.0_precision, 'out', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'sod configuration accepted')
    call check_true(trim(case_1d) == 'sod', 'case read')
    call check_true(mesh_num == 64, 'mesh read')
    call check_close(cfl, 0.2_precision, 1.0e-14_precision, 'cfl read')
    call check_close(tend, 0.15_precision, 1.0e-14_precision, 'tend defaulted')
    call check_true(trim(output_dir) == 'out', 'output_dir read')
    call check_true(write_sample, 'write_sample read')
    call delete_nml()

    call write_nml('shu_osher', 200, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'shu configuration accepted')
    call check_close(tend, 1.8_precision, 1.0e-14_precision, 'shu tend defaulted')
    call delete_nml()

    call write_nml('sedov', 64, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'unknown case rejected')
    call delete_nml()

    call write_nml('sod', 7, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'mesh below 8 rejected')
    call delete_nml()

    call write_nml('sod', 64, 0.0_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'zero cfl rejected')
    call delete_nml()

    call write_nml('sod', 64, 0.1_precision, -1.0_precision, '', 'AUTO')
    call e1d_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'empty output_dir rejected')
    call delete_nml()

    ! ---------------------------------------------------------------
    ! Operator pairing
    ! ---------------------------------------------------------------
    difference_operator = 'AUTO'
    call e1d_resolve_operator('QR', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS4', 'AUTO QR -> HCS4')

    call e1d_resolve_operator('Z5', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'AUTO Z5 -> HCS6')

    call e1d_resolve_operator('MUSCL', op, ierr, msg)
    call check_true(ierr /= 0, 'AUTO MUSCL rejected')

    difference_operator = 'HCS6'
    call e1d_resolve_operator('Z3', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'explicit HCS6 accepted')

    difference_operator = 'NOPE'
    call e1d_resolve_operator('Z3', op, ierr, msg)
    call check_true(ierr /= 0, 'unknown operator rejected')

    ! ---------------------------------------------------------------
    ! Exact-solution configuration (same namelist group)
    ! ---------------------------------------------------------------
    call write_nml('sod', 64, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_exact_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'exact sod accepted')
    call check_close(tend, 0.15_precision, 1.0e-14_precision, 'exact tend defaulted')
    call delete_nml()

    call write_nml('shu_osher', 64, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_exact_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'exact shu_osher rejected')
    call delete_nml()

    call write_nml('sod', 1, 0.1_precision, -1.0_precision, 'results', 'AUTO')
    call e1d_exact_control_init('euler1d_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'exact mesh 1 rejected')
    call delete_nml()

    call finish()

contains

    subroutine write_nml(case_name, mesh, cfl_in, tend_in, outdir, operator)
        ! Write a complete input.nml for one test case.
        character(len=*), intent(in) :: case_name
        integer, intent(in) :: mesh
        real(precision), intent(in) :: cfl_in, tend_in
        character(len=*), intent(in) :: outdir, operator
        integer :: u

        open (newunit=u, file='euler1d_params_test.nml', status='replace', action='write')
        write (u, '(a)') '&euler1d_nml'
        write (u, '(a)') "  case_1d = '"//trim(case_name)//"',"
        write (u, '(a,i0)') '  mesh_num = ', mesh
        write (u, '(a,es24.16)') '  cfl = ', cfl_in
        write (u, '(a,es24.16)') '  tend = ', tend_in
        write (u, '(a)') "  output_dir = '"//trim(outdir)//"',"
        write (u, '(a)') '  write_sample = .true.'
        write (u, '(a)') '/'
        write (u, '(a)') '&difference_nml'
        write (u, '(a)') "  difference_operator = '"//trim(operator)//"'"
        write (u, '(a)') '/'
        close (u)
    end subroutine write_nml

    subroutine delete_nml()
        ! Remove the temporary namelist file.
        integer :: u

        open (newunit=u, file='euler1d_params_test.nml', status='old')
        close (u, status='delete')
    end subroutine delete_nml

end program test_control
