program test_control
    ! Unit tests for accuracy configuration parsing, validation, and
    ! automatic difference-operator pairing.
    use precision_select
    use accuracy_control
    use test_util
    implicit none

    integer :: ierr
    character(len=256) :: msg
    character(len=32) :: op

    ! ---------------------------------------------------------------
    ! 1D configuration
    ! ---------------------------------------------------------------
    call write_nml('sin2', 3)
    call accuracy_control_init_1d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr == 0, '1d control accepted')
    call check_true(nbase == 20 .and. levels == 3, '1d control read')
    call check_close(cfl, 0.2_precision, 0.0_precision, '1d cfl read')
    call check_close(tend, 1.5_precision, 0.0_precision, '1d tend read')
    call check_true(trim(initial_1d) == 'sin2', '1d initial read')
    call delete_nml()

    call write_nml('gauss', 3)
    call accuracy_control_init_1d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'gauss initial accepted')
    call check_true(trim(initial_1d) == 'gauss', 'gauss initial read')
    call delete_nml()

    call write_nml('NOPE', 3)
    call accuracy_control_init_1d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'unknown 1D initial rejected')
    call delete_nml()

    call write_nml('sin2', 21)
    call accuracy_control_init_1d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'levels > 20 rejected')
    call delete_nml()

    call write_nml('sin2', 1)
    call accuracy_control_init_1d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr == 0, 'levels = 1 accepted')
    call delete_nml()

    ! ---------------------------------------------------------------
    ! 2D configuration
    ! ---------------------------------------------------------------
    call write_nml('sin2', 3)
    call accuracy_control_init_2d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr == 0, '2d control accepted')
    call check_true(trim(case_2d) == 'vortex', '2d case read')
    call check_close(sine_rho0, 1.0_precision, 0.0_precision, 'sine rho0 read')
    call check_close(sine_amp, 0.2_precision, 0.0_precision, 'sine amp read')
    call delete_nml()

    ! A sine amplitude at or above the base density can drive rho to zero.
    call write_nml('sin2', 3, 1.0_precision, 'sine')
    call accuracy_control_init_2d('accuracy_params_test.nml', ierr, msg)
    call check_true(ierr /= 0, 'sine_amp >= sine_rho0 rejected')
    call check_true(index(msg, 'sine_amp') > 0, 'density message mentions sine_amp')
    call delete_nml()

    ! ---------------------------------------------------------------
    ! Operator pairing
    ! ---------------------------------------------------------------
    difference_operator = 'AUTO'

    call accuracy_resolve_operator('Z3', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS4', 'AUTO Z3 -> HCS4')

    call accuracy_resolve_operator('Z5', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'AUTO Z5 -> HCS6')

    call accuracy_resolve_operator('QR', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS4', 'AUTO QR -> HCS4')

    call accuracy_resolve_operator('MUSCL', op, ierr, msg)
    call check_true(ierr /= 0, 'AUTO rejects unsupported order')

    difference_operator = 'HCS6'
    call accuracy_resolve_operator('Z3', op, ierr, msg)
    call check_true(ierr == 0 .and. trim(op) == 'HCS6', 'explicit HCS6 accepted')

    difference_operator = 'NOPE'
    call accuracy_resolve_operator('Z3', op, ierr, msg)
    call check_true(ierr /= 0, 'unknown operator rejected')

    call finish()

contains

    subroutine write_nml(initial_name, nlev, amp, case_name)
        character(len=*), intent(in) :: initial_name
        integer, intent(in) :: nlev
        real(precision), intent(in), optional :: amp
        character(len=*), intent(in), optional :: case_name
        real(precision) :: sine_amplitude
        character(len=16) :: cname
        integer :: u

        sine_amplitude = 0.2_precision
        if (present(amp)) sine_amplitude = amp
        cname = 'vortex'
        if (present(case_name)) cname = case_name

        open (newunit=u, file='accuracy_params_test.nml', status='replace', action='write')
        write (u, '(a)') '&accuracy_nml'
        write (u, '(a)') "  write_sample = .false., output_dir = 'results'"
        write (u, '(a)') '/'
        write (u, '(a)') '&accuracy_1d_nml'
        write (u, '(a)') "  initial_1d = '"//trim(initial_name)//"', nbase = 20, levels = "// &
            trim(adjustl(int_char(nlev)))//", cfl = 0.2d0, tend = 1.5d0"
        write (u, '(a)') '/'
        write (u, '(a)') '&accuracy_2d_nml'
        write (u, '(a)') "  case_2d = '"//trim(cname)//"',"
        write (u, '(a,es24.16)') "  nbase = 16, levels = 2, cfl = 0.1d0, tend = 1.0d0, "// &
            "sine_rho0 = 1.0d0, sine_amp = ", sine_amplitude
        write (u, '(a)') '/'
        write (u, '(a)') '&scheme_nml'
        write (u, '(a)') "  scheme = 'Z3'"
        write (u, '(a)') '/'
        write (u, '(a)') '&difference_nml'
        write (u, '(a)') "  difference_operator = 'AUTO'"
        write (u, '(a)') '/'
        close (u)
    end subroutine write_nml

    function int_char(value) result(str)
        integer, intent(in) :: value
        character(len=16) :: str

        write (str, '(i0)') value
    end function int_char

    subroutine delete_nml()
        integer :: u

        open (newunit=u, file='accuracy_params_test.nml', status='old')
        close (u, status='delete')
    end subroutine delete_nml

end program test_control
