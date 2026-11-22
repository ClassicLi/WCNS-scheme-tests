program test_params
    ! Unit tests for ADR, WENO, and WCNS parameter validation.
    use precision_select
    use adr_control
    use weno_params
    use wcns_params
    use test_util
    implicit none

    integer :: ierr
    character(len=256) :: msg

#ifdef PRECISION_DP
    call check_close(adr_default_dt(), 1.0e-10_precision, 0.0_precision, 'default dt dp')
#else
    call check_close(adr_default_dt(), 1.0e-20_precision, 0.0_precision, 'default dt qp')
#endif

    call write_nml('  mesh_num = 33, dt = 0.0d0')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'valid odd mesh accepted')
    call check_true(mesh_num == 33, 'mesh_num read')
    call check_close(dt, adr_default_dt(), 0.0_precision, 'dt auto')
    call delete_nml()

    call write_nml('  mesh_num = 33, dt = 1.0d-7')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'explicit dt accepted')
    call check_close(dt, 1.0e-7_precision, 0.0_precision, 'explicit dt preserved')
    call delete_nml()

    call write_nml('  mesh_num = 34, dt = 0.0d0')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'even mesh rejected')
    call check_true(index(msg, 'odd') > 0, 'even mesh message mentions odd')
    call delete_nml()

    call write_nml('  mesh_num = 3, dt = 0.0d0')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'mesh below 5 rejected')
    call delete_nml()

    call write_nml('  mesh_num = 33, dt = -1.0d0')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'negative dt rejected')
    call delete_nml()

    call write_nml('  dt = 0.0d0')
    call adr_control_init('adr_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'omitted mesh_num accepted')
    call check_true(mesh_num == 827, 'default mesh_num 827')
    call delete_nml()

    call check_true(weno_scheme_allowed('WENO5_LINEAR'), 'linear WENO scheme allowed')
    call check_true(.not. weno_scheme_allowed('NOT_A_SCHEME'), 'unknown WENO scheme rejected')

    scheme = 'NOT_A_SCHEME'
    call validate_weno_params(ierr, msg)
    call check_true(ierr /= 0, 'validate rejects unknown WENO scheme')

    scheme = 'WENO5_LINEAR'
    call validate_weno_params(ierr, msg)
    call check_true(ierr == 0, 'validate accepts linear WENO scheme')

    call check_true(difference_operator_allowed('HCS4'), 'HCS4 operator allowed')
    call check_true(difference_operator_allowed('HCS6'), 'HCS6 operator allowed')
    call check_true(.not. difference_operator_allowed('NOPE'), 'unknown operator rejected')

    difference_operator = 'NOPE'
    call validate_wcns_params(ierr, msg)
    call check_true(ierr /= 0, 'validate rejects unknown operator')

    difference_operator = 'HCS4'
    call validate_wcns_params(ierr, msg)
    call check_true(ierr == 0, 'validate accepts HCS4')

    call finish()

contains

    subroutine write_nml(body)
        character(len=*), intent(in) :: body
        integer :: u

        open(newunit=u, file='adr_params_test.nml', status='replace', action='write')
        write(u, '(a)') '&adr_nml'
        write(u, '(a)') trim(body)
        write(u, '(a)') '/'
        close(u)
    end subroutine write_nml

    subroutine delete_nml()
        integer :: u

        open(newunit=u, file='adr_params_test.nml', status='old')
        close(u, status='delete')
    end subroutine delete_nml

end program test_params
