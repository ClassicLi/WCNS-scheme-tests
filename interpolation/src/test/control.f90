module control
    use precision_select
    implicit none

    integer(integer_precision) :: ini_mesh_num = 10_integer_precision
    integer(integer_precision) :: iter_times = 10_integer_precision
    integer(integer_precision) :: test_func_type = 3_integer_precision
    integer(integer_precision) :: n_order_critical_point = 1_integer_precision
    integer(integer_precision) :: udf_id = 4_integer_precision
    character(len=16) :: nonlinear_baseline = 'same_order'
    logical :: output_fig = .true.

    namelist /test_nml/ ini_mesh_num, iter_times, test_func_type, &
        n_order_critical_point, udf_id, nonlinear_baseline, output_fig

contains

    subroutine read_test_namelist(unit, ierr)
        integer, intent(in) :: unit
        integer, intent(out), optional :: ierr
        integer :: ios

        rewind(unit)
        read(unit, nml=test_nml, iostat=ios)
        if (present(ierr)) ierr = ios
    end subroutine read_test_namelist

    subroutine validate_test_params(ierr, msg)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (ini_mesh_num < 1_integer_precision .or. iter_times < 1_integer_precision) then
            ierr = 1
            msg = 'ini_mesh_num and iter_times must be positive'
            return
        end if

        if (test_func_type < 1_integer_precision .or. test_func_type > 4_integer_precision) then
            ierr = 1
            msg = 'test_func_type must be 1..4'
            return
        end if

        if (n_order_critical_point < 0_integer_precision .or. &
            n_order_critical_point > 2_integer_precision) then
            ierr = 1
            msg = 'n_order_critical_point must be 0..2'
            return
        end if

        if (udf_id < 1_integer_precision .or. udf_id > 6_integer_precision) then
            ierr = 1
            msg = 'udf_id must be 1..6'
            return
        end if

        select case (trim(nonlinear_baseline))
        case ('same_order', 'poly3')
        case default
            ierr = 1
            msg = 'nonlinear_baseline must be same_order or poly3'
            return
        end select
    end subroutine validate_test_params

end module control
