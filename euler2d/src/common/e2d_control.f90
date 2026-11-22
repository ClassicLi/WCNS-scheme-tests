module e2d_control
    ! Configuration, case defaults, and operator pairing for the 2D Euler tests.
    !
    ! Reads &euler2d_nml, &scheme_nml, and &difference_nml, resolves the case
    ! defaults, and validates the resulting mesh and time parameters.
    use precision_select
    use e2d_cases, only: e2d_case_meta, e2d_case_allowed, e2d_case_meta_get
    use e2d_globals, only: e2d_set_gama
    use scheme_params, only: read_scheme_namelist, validate_scheme_params
    use scheme_dispatch, only: scheme_order
    implicit none

    ! &euler2d_nml
    character(len=32) :: case_2d = 'riemann1'
    integer(integer_precision) :: mesh_num = 0_integer_precision
    real(precision) :: cfl = 0.3_precision
    real(precision) :: tend = -1.0_precision
    real(precision) :: tend2 = -1.0_precision
    real(precision) :: gamma = 0.0_precision
    character(len=256) :: output_dir = 'results'
    logical :: write_sample = .false.

    ! &difference_nml
    character(len=32) :: difference_operator = 'AUTO'

    ! Derived grid size and physical-coordinate origin.
    integer(integer_precision) :: nz = 0_integer_precision
    integer(integer_precision) :: nx = 0_integer_precision
    real(precision) :: x_origin = 0.0_precision
    real(precision) :: z_origin = 0.0_precision

    namelist /euler2d_nml/ case_2d, mesh_num, cfl, tend, tend2, gamma, &
                           output_dir, write_sample
    namelist /difference_nml/ difference_operator

contains

    subroutine e2d_control_init(input_file, ierr, msg)
        ! Reset defaults, read the namelist groups, resolve the case, and validate.
        character(len=*), intent(in) :: input_file
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        type(e2d_case_meta) :: meta
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
        read (u, nml=euler2d_nml, iostat=ios)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &euler2d_nml'
            close (u)
            return
        end if

        call read_scheme_namelist(u, ios)
        if (ios /= 0) then
            ierr = 3
            msg = 'cannot read &scheme_nml'
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

        call validate_scheme_params(ierr, msg)
        if (ierr /= 0) then
            msg = 'invalid scheme parameters: ' // trim(msg)
            return
        end if

        if (.not. e2d_case_allowed(case_2d)) then
            ierr = 10
            msg = 'unknown case_2d: ' // trim(case_2d)
            return
        end if

        call e2d_case_meta_get(case_2d, meta, ierr, msg)
        if (ierr /= 0) return

        ! A zero or negative value selects the case default.
        if (mesh_num == 0) mesh_num = meta%default_nz
        if (tend <= 0.0_precision) tend = meta%default_tend
        if (tend2 < 0.0_precision) tend2 = meta%default_tend2
        if (gamma <= 0.0_precision) gamma = meta%default_gama

        call e2d_set_gama(gamma)

        call validate_config(meta, ierr, msg)
        if (ierr /= 0) return

        nz = mesh_num
        nx = nz*meta%ratio_num/meta%ratio_den
        x_origin = meta%x_origin
        z_origin = meta%z_origin
    end subroutine e2d_control_init

    subroutine reset_defaults()
        ! Restore every namelist variable and derived value for repeated calls.
        case_2d = 'riemann1'
        mesh_num = 0_integer_precision
        cfl = 0.3_precision
        tend = -1.0_precision
        tend2 = -1.0_precision
        gamma = 0.0_precision
        output_dir = 'results'
        write_sample = .false.
        difference_operator = 'AUTO'
        nz = 0_integer_precision
        nx = 0_integer_precision
        x_origin = 0.0_precision
        z_origin = 0.0_precision
    end subroutine reset_defaults

    subroutine validate_config(meta, ierr, msg)
        ! Reject values that would break the numerical driver.
        type(e2d_case_meta), intent(in) :: meta
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (mesh_num < 8_integer_precision) then
            ierr = 11
            msg = 'mesh_num must be at least 8'
            return
        end if

        if (mod(mesh_num, meta%ratio_den) /= 0) then
            ierr = 12
            msg = 'mesh_num must be divisible by the case ratio denominator'
            return
        end if

        if (cfl <= 0.0_precision) then
            ierr = 13
            msg = 'cfl must be positive'
            return
        end if

        if (tend2 > 0.0_precision) then
            if (trim(case_2d) /= 'shock_bubble') then
                ierr = 14
                msg = 'tend2 is only supported for shock_bubble'
                return
            end if
            if (tend2 >= tend) then
                ierr = 15
                msg = 'tend2 must be less than tend'
                return
            end if
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 16
            msg = 'output_dir must not be empty'
            return
        end if
    end subroutine validate_config

    function e2d_operator_allowed(name) result(ok)
        ! Return true for a known difference-operator ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6')
            ok = .true.
        case default
            ok = .false.
        end select
    end function e2d_operator_allowed

    subroutine e2d_resolve_operator(scheme_name, op, ierr, msg)
        ! Resolve AUTO from the scheme order, or validate an explicit operator.
        character(len=*), intent(in) :: scheme_name
        character(len=32), intent(out) :: op
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: p
        character(len=16) :: pstr
        character(len=32) :: recommended

        ierr = 0
        msg = ''
        op = ''

        if (trim(difference_operator) == 'AUTO') then
            p = scheme_order(scheme_name)
            recommended = auto_operator(p)
            if (len_trim(recommended) == 0) then
                ierr = 40
                msg = 'AUTO has no difference operator for scheme ' // &
                      trim(scheme_name) // ' (order 0 or unsupported)'
                return
            end if
            op = recommended
        else
            if (.not. e2d_operator_allowed(difference_operator)) then
                ierr = 41
                msg = 'unknown difference operator: ' // trim(difference_operator)
                return
            end if
            op = trim(difference_operator)

            p = scheme_order(scheme_name)
            recommended = auto_operator(p)
            if (len_trim(recommended) > 0 .and. trim(op) /= trim(recommended)) then
                write (pstr, '(i0)') p
                write (*, '(a)') 'warning: ' // trim(scheme_name) // ' (order ' // &
                    trim(pstr) // ') with ' // trim(op) // ' is not the default pairing'
            end if
        end if
    end subroutine e2d_resolve_operator

    function auto_operator(p) result(op)
        ! Map an interpolation order to its recommended difference operator.
        integer, intent(in) :: p
        character(len=32) :: op

        select case (p)
        case (2)
            op = 'WCNSE2'
        case (3)
            op = 'HCS4'
        case (4)
            op = 'WCNSE4'
        case (5)
            op = 'HCS6'
        case (6)
            op = 'WCNSE6'
        case default
            op = ''
        end select
    end function auto_operator

    subroutine e2d_ensure_output_dirs(root)
        ! Create the output directory (and any missing parents).
        character(len=*), intent(in) :: root

        call execute_command_line('mkdir -p "' // trim(root) // '"')
    end subroutine e2d_ensure_output_dirs

end module e2d_control
