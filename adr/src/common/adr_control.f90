module adr_control
    ! Global ADR parameters, namelist input, and input validation.
    use precision_select
    implicit none

    ! Mathematical constants used by the ADR analysis.
    real(precision), parameter :: PI = acos(-1.0_precision)
    complex(precision), parameter :: II = cmplx(0.0_precision, 1.0_precision, precision)

    ! Problem defaults; overridden by &adr_nml.
    integer(integer_precision) :: mesh_num = 827_integer_precision
    real(precision) :: dt = 0.0_precision
    character(len=256) :: output_dir = 'results'
    logical :: write_tecplot_sample = .false.

    namelist /adr_nml/ mesh_num, dt, output_dir, write_tecplot_sample

contains

    function adr_default_dt() result(dt0)
        ! Precision-dependent default time step (qp: 1e-20, dp: 1e-10).
        real(precision) :: dt0

#ifdef PRECISION_DP
        dt0 = 1.0e-10_precision
#else
        dt0 = 1.0e-20_precision
#endif
    end function adr_default_dt

    subroutine adr_control_init(input_file, ierr, msg)
        ! Reset defaults, read &adr_nml, validate, and resolve dt = 0.
        character(len=*), intent(in) :: input_file
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: u, ios

        ierr = 0
        msg = ''

        mesh_num = 827_integer_precision
        dt = 0.0_precision
        output_dir = 'results'
        write_tecplot_sample = .false.

        open(newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot open input file: ' // trim(input_file)
            return
        end if

        read(u, nml=adr_nml, iostat=ios)
        close(u)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &adr_nml'
            return
        end if

        call adr_control_validate(ierr, msg)
        if (ierr /= 0) return

        if (dt == 0.0_precision) dt = adr_default_dt()
    end subroutine adr_control_init

    subroutine adr_control_validate(ierr, msg)
        ! Check the odd mesh, nonnegative dt, and nonempty output directory.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (mesh_num < 5_integer_precision .or. mod(mesh_num, 2_integer_precision) == 0) then
            ierr = 3
            msg = 'mesh_num must be odd and >= 5'
            return
        end if

        if (dt < 0.0_precision) then
            ierr = 4
            msg = 'dt must be nonnegative'
            return
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 5
            msg = 'output_dir must not be empty'
            return
        end if
    end subroutine adr_control_validate

end module adr_control
