module cw_control
    ! Global settings for the complex-wave test.
    !
    ! The physics is fixed to the reference problem: the linear advection
    ! equation u_t + u_x = 0 on [-1, 1] with unit wave speed, so one wave
    ! period is T = 2.  Mesh size, number of periods, CFL number, and the
    ! output directory are read from the &complex_wave_nml group.
    use precision_select
    implicit none

    ! Defaults; every field can be overridden by &complex_wave_nml.
    integer(integer_precision) :: mesh_num = 200_integer_precision
    integer(integer_precision) :: periods = 4_integer_precision
    real(precision) :: cfl = 0.3_precision
    character(len=256) :: output_dir = 'results'

    ! Fixed physical parameters of the reference problem.
    real(precision), parameter :: domain_length = 2.0_precision
    real(precision), parameter :: wave_speed = 1.0_precision

    namelist /complex_wave_nml/ mesh_num, periods, cfl, output_dir

contains

    subroutine cw_control_init(input_file, ierr, msg)
        ! Reset the defaults, read &complex_wave_nml, and validate the result.
        character(len=*), intent(in) :: input_file
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: u, ios

        ierr = 0
        msg = ''

        ! Reset first so repeated calls in one process are deterministic.
        mesh_num = 200_integer_precision
        periods = 4_integer_precision
        cfl = 0.3_precision
        output_dir = 'results'

        open(newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            msg = 'cannot open input file: ' // trim(input_file)
            return
        end if

        read(u, nml=complex_wave_nml, iostat=ios)
        close(u)
        if (ios /= 0) then
            ierr = 2
            msg = 'cannot read &complex_wave_nml'
            return
        end if

        call cw_control_validate(ierr, msg)
    end subroutine cw_control_init

    subroutine cw_control_validate(ierr, msg)
        ! Reject values that would break the initial condition or time loop.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        ! The initial condition uses four ten-percent segments; below 20
        ! those segments cannot be resolved on the grid.
        if (mesh_num < 20_integer_precision) then
            ierr = 3
            msg = 'mesh_num must be >= 20'
            return
        end if

        ! At least one period is required, and a negative count is meaningless.
        if (periods < 1_integer_precision) then
            ierr = 4
            msg = 'periods must be a positive integer'
            return
        end if

        ! cfl = 0 would produce dt = 0 and stall the time loop.
        if (cfl <= 0.0_precision) then
            ierr = 5
            msg = 'cfl must be positive'
            return
        end if

        if (len_trim(output_dir) == 0) then
            ierr = 6
            msg = 'output_dir must not be empty'
            return
        end if
    end subroutine cw_control_validate

    function cw_period() result(T)
        ! One wave period on the periodic domain.
        real(precision) :: T

        T = domain_length / wave_speed
    end function cw_period

    function cw_tend() result(Tend)
        ! Final time after the requested number of periods.
        real(precision) :: Tend

        Tend = real(periods, precision) * cw_period()
    end function cw_tend

end module cw_control
