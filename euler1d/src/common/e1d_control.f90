module e1d_control
    ! Configuration, validation, and operator pairing for the 1D Euler tests.
    !
    ! Both drivers read the same &euler1d_nml.  The exact-solution driver
    ! additionally rejects cases other than sod and lax, and ignores cfl.
    use precision_select
    use e1d_initial, only: e1d_case_allowed, e1d_case_tend
    use scheme_dispatch, only: scheme_order
    implicit none

    ! &euler1d_nml
    character(len=16) :: case_1d = 'sod'
    integer(integer_precision) :: mesh_num = 200_integer_precision
    real(precision) :: cfl = 0.1_precision
    real(precision) :: tend = -1.0_precision
    character(len=256) :: output_dir = 'results'
    logical :: write_sample = .false.

    ! &difference_nml
    character(len=32) :: difference_operator = 'AUTO'

    namelist /euler1d_nml/ case_1d, mesh_num, cfl, tend, output_dir, write_sample
    namelist /difference_nml/ difference_operator

contains

    subroutine e1d_control_init(input_file, ierr, msg)
        ! Reset defaults, read the solver groups, and validate them.
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
        read (u, nml=euler1d_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &euler1d_nml'
            close (u)
            return
        end if

        rewind (u)
        read (u, nml=difference_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 3
            msg = 'cannot read &difference_nml'
            close (u)
            return
        end if

        close (u)
        call e1d_validate_solver(ierr, msg)
    end subroutine e1d_control_init

    subroutine e1d_exact_control_init(input_file, ierr, msg)
        ! Read the shared group and validate the exact-solution subset.
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
        read (u, nml=euler1d_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &euler1d_nml'
            close (u)
            return
        end if

        close (u)
        call e1d_validate_exact(ierr, msg)
    end subroutine e1d_exact_control_init

    subroutine reset_defaults()
        ! Restore every namelist variable so repeated calls are deterministic.
        case_1d = 'sod'
        mesh_num = 200_integer_precision
        cfl = 0.1_precision
        tend = -1.0_precision
        output_dir = 'results'
        write_sample = .false.
        difference_operator = 'AUTO'
    end subroutine reset_defaults

    subroutine e1d_validate_solver(ierr, msg)
        ! Reject values that would break the numerical driver.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. e1d_case_allowed(case_1d)) then
            ierr = 10
            msg = 'unknown case_1d: ' // trim(case_1d)
            return
        end if

        if (mesh_num < 8_integer_precision) then
            ierr = 11
            msg = 'mesh_num must be at least 8'
            return
        end if

        if (cfl <= 0.0_precision) then
            ierr = 12
            msg = 'cfl must be positive'
            return
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 13
            msg = 'output_dir must not be empty'
            return
        end if

        if (tend <= 0.0_precision) tend = e1d_case_tend(case_1d)
    end subroutine e1d_validate_solver

    subroutine e1d_validate_exact(ierr, msg)
        ! Reject cases that have no closed-form Riemann solution here.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (trim(case_1d) /= 'sod' .and. trim(case_1d) /= 'lax') then
            ierr = 20
            msg = 'exact solution only supports sod and lax, got: ' // trim(case_1d)
            return
        end if

        if (mesh_num < 2_integer_precision) then
            ierr = 21
            msg = 'mesh_num must be at least 2'
            return
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 22
            msg = 'output_dir must not be empty'
            return
        end if

        if (tend <= 0.0_precision) tend = e1d_case_tend(case_1d)
    end subroutine e1d_validate_exact

    function e1d_operator_allowed(name) result(ok)
        ! Return true for a known difference-operator ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6')
            ok = .true.
        case default
            ok = .false.
        end select
    end function e1d_operator_allowed

    subroutine e1d_resolve_operator(scheme_name, op, ierr, msg)
        ! Resolve AUTO from the scheme order, or validate an explicit operator.
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
            if (.not. e1d_operator_allowed(difference_operator)) then
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
    end subroutine e1d_resolve_operator

    subroutine e1d_ensure_output_dirs(root)
        ! Create the result subdirectories used by both drivers.
        character(len=*), intent(in) :: root

        call execute_command_line('mkdir -p "' // trim(root) // '/csv" "' // &
            trim(root) // '/logs" "' // trim(root) // '/figures"')
    end subroutine e1d_ensure_output_dirs

end module e1d_control
