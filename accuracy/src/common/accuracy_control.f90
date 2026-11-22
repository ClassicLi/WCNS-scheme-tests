module accuracy_control
    ! Configuration, validation, and operator pairing for the accuracy tests.
    !
    ! The two drivers share one input.nml.  Each driver reads &accuracy_nml,
    ! its own case group, and &difference_nml; the interpolation scheme comes
    ! from &scheme_nml through scheme_params.
    use precision_select
    use scheme_dispatch, only: scheme_order
    use accuracy_order, only: max_levels
    implicit none

    ! &accuracy_nml
    logical :: write_sample = .false.
    character(len=256) :: output_dir = 'results'

    ! &accuracy_1d_nml
    character(len=16) :: initial_1d = 'sin2'
    integer(integer_precision) :: nbase = 10_integer_precision
    integer(integer_precision) :: levels = 8_integer_precision
    real(precision) :: cfl = 0.1_precision
    real(precision) :: tend = 2.0_precision

    ! &accuracy_2d_nml
    character(len=16) :: case_2d = 'vortex'
    real(precision) :: sine_rho0 = 1.0_precision
    real(precision) :: sine_amp = 0.2_precision

    ! &difference_nml
    character(len=32) :: difference_operator = 'AUTO'

    namelist /accuracy_nml/ write_sample, output_dir
    namelist /accuracy_1d_nml/ initial_1d, nbase, levels, cfl, tend
    namelist /accuracy_2d_nml/ case_2d, nbase, levels, cfl, tend, sine_rho0, sine_amp
    namelist /difference_nml/ difference_operator

contains

    subroutine accuracy_control_init_1d(input_file, ierr, msg)
        ! Reset defaults, read the 1D groups, and validate them.
        character(len=*), intent(in) :: input_file
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: u, ios

        call reset_defaults()
        ierr = 0
        msg = ''

        open (newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot open input file: ' // trim(input_file)
            return
        end if

        rewind (u)
        read (u, nml=accuracy_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &accuracy_nml'
            close (u)
            return
        end if

        rewind (u)
        read (u, nml=accuracy_1d_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            msg = 'cannot read &accuracy_1d_nml'
            close (u)
            return
        end if

        rewind (u)
        read (u, nml=difference_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 4
            msg = 'cannot read &difference_nml'
            close (u)
            return
        end if

        close (u)
        call accuracy_validate_1d(ierr, msg)
    end subroutine accuracy_control_init_1d

    subroutine accuracy_control_init_2d(input_file, ierr, msg)
        ! Reset defaults, read the 2D groups, and validate them.
        character(len=*), intent(in) :: input_file
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: u, ios

        call reset_defaults()
        ierr = 0
        msg = ''

        open (newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot open input file: ' // trim(input_file)
            return
        end if

        rewind (u)
        read (u, nml=accuracy_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &accuracy_nml'
            close (u)
            return
        end if

        rewind (u)
        read (u, nml=accuracy_2d_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            msg = 'cannot read &accuracy_2d_nml'
            close (u)
            return
        end if

        rewind (u)
        read (u, nml=difference_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 4
            msg = 'cannot read &difference_nml'
            close (u)
            return
        end if

        close (u)
        call accuracy_validate_2d(ierr, msg)
    end subroutine accuracy_control_init_2d

    subroutine reset_defaults()
        ! Restore every namelist variable so repeated calls are deterministic.
        write_sample = .false.
        output_dir = 'results'

        initial_1d = 'sin2'
        nbase = 10_integer_precision
        levels = 8_integer_precision
        cfl = 0.1_precision
        tend = 2.0_precision

        case_2d = 'vortex'
        sine_rho0 = 1.0_precision
        sine_amp = 0.2_precision

        difference_operator = 'AUTO'
    end subroutine reset_defaults

    subroutine accuracy_validate_1d(ierr, msg)
        ! Reject values that would break the 1D driver.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. accuracy_initial_allowed(initial_1d)) then
            ierr = 10
            msg = 'unknown initial_1d: ' // trim(initial_1d)
            return
        end if

        call validate_common(ierr, msg)
        if (ierr /= 0) return

        ! The 1D reference problem is fixed to [-1, 1] with Tend = 2.
        if (abs(tend - 2.0_precision) > 1.0e-12_precision) then
            write (*, '(a)') 'warning: 1D reference problem uses tend = 2'
        end if
    end subroutine accuracy_validate_1d

    subroutine accuracy_validate_2d(ierr, msg)
        ! Reject values that would break the 2D driver.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (trim(case_2d) /= 'vortex' .and. trim(case_2d) /= 'sine') then
            ierr = 20
            msg = 'unknown case_2d: ' // trim(case_2d)
            return
        end if

        call validate_common(ierr, msg)
        if (ierr /= 0) return

        if (sine_rho0 <= 0.0_precision) then
            ierr = 21
            msg = 'sine_rho0 must be positive'
            return
        end if

        if (sine_amp < 0.0_precision) then
            ierr = 22
            msg = 'sine_amp must be non-negative'
            return
        end if

        ! The sine density must stay positive on the whole grid.
        if (trim(case_2d) == 'sine' .and. sine_amp >= sine_rho0) then
            ierr = 23
            msg = 'sine_amp must be smaller than sine_rho0 (density must stay positive)'
            return
        end if

        ! The 2D reference problems use Tend = 1.
        if (abs(tend - 1.0_precision) > 1.0e-12_precision) then
            write (*, '(a)') 'warning: 2D reference problems use tend = 1'
        end if
    end subroutine accuracy_validate_2d

    subroutine validate_common(ierr, msg)
        ! Checks shared by the 1D and 2D drivers.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (nbase < 4_integer_precision) then
            ierr = 30
            msg = 'nbase must be >= 4'
            return
        end if

        if (levels < 1_integer_precision .or. levels > max_levels) then
            ierr = 31
            msg = 'levels must be between 1 and 20'
            return
        end if

        if (cfl <= 0.0_precision) then
            ierr = 32
            msg = 'cfl must be positive'
            return
        end if

        if (tend <= 0.0_precision) then
            ierr = 33
            msg = 'tend must be positive'
            return
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 34
            msg = 'output_dir must not be empty'
            return
        end if
    end subroutine validate_common

    function accuracy_initial_allowed(name) result(ok)
        ! Return true for a known 1D initial-profile ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('sin2', 'sin', 'sin3', 'sin_phase', 'gauss', 'sqrt', 'rational')
            ok = .true.
        case default
            ok = .false.
        end select
    end function accuracy_initial_allowed

    function accuracy_operator_allowed(name) result(ok)
        ! Return true for a known difference-operator ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6')
            ok = .true.
        case default
            ok = .false.
        end select
    end function accuracy_operator_allowed

    subroutine accuracy_resolve_operator(scheme_name, op, ierr, msg)
        ! Resolve 'AUTO' from the scheme order, or validate an explicit operator.
        character(len=*), intent(in) :: scheme_name
        character(len=32), intent(out) :: op
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: p
        character(len=16) :: pstr

        ierr = 0
        msg = ''
        op = ''

        if (trim(difference_operator) == 'AUTO') then
            p = scheme_order(scheme_name)
            select case (p)
            case (3)
                op = 'HCS4'
            case (5)
                op = 'HCS6'
            case default
                ierr = 40
                msg = 'AUTO has no difference operator for scheme ' // &
                      trim(scheme_name) // ' (order 0 or unsupported)'
                return
            end select
        else
            if (.not. accuracy_operator_allowed(difference_operator)) then
                ierr = 41
                msg = 'unknown difference operator: ' // trim(difference_operator)
                return
            end if
            op = trim(difference_operator)

            p = scheme_order(scheme_name)
            write (pstr, '(i0)') p
            if ((p == 3 .and. trim(op) == 'HCS6') .or. &
                (p == 5 .and. trim(op) == 'HCS4')) then
                write (*, '(a)') 'warning: ' // trim(scheme_name) // ' (order ' // &
                    trim(pstr) // ') with ' // trim(op) // ' is not the default pairing'
            end if
        end if
    end subroutine accuracy_resolve_operator

    subroutine accuracy_ensure_output_dirs(root)
        ! Create the result subdirectories used by both drivers.
        character(len=*), intent(in) :: root

        call execute_command_line('mkdir -p "' // trim(root) // '/csv" "' // &
            trim(root) // '/logs" "' // trim(root) // '/figures"')
    end subroutine accuracy_ensure_output_dirs

end module accuracy_control
