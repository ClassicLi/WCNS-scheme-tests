program main
    use precision_select
    use control, only: ini_mesh_num, output_fig_flag => output_fig, &
                       read_test_namelist, validate_test_params
    use scheme_params, only: scheme, read_scheme_namelist, validate_scheme_params
    use solver, only: procedure, time_used, func_label, calc_init, interpolate
    use io, only: output_exact_standard_func, output_order, output_csv, output_fig
    use scheme_dispatch, only: scheme_output_name
    implicit none

    integer :: unit, ios, ierr
    character(len=128) :: msg
    character(len=64) :: method

    open(newunit=unit, file='input.nml', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*, '(a)') 'Error: input.nml not found'
        error stop 1
    end if

    call read_test_namelist(unit, ios)
    if (ios /= 0) then
        write(*, '(a)') 'Error reading &test_nml'
        error stop 1
    end if

    call read_scheme_namelist(unit, ios)
    if (ios /= 0) then
        write(*, '(a)') 'Error reading &scheme_nml'
        error stop 1
    end if
    close(unit)

    call validate_test_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'Error: ' // trim(msg)
        error stop 1
    end if

    call validate_scheme_params(ierr, msg)
    if (ierr /= 0) then
        write(*, '(a)') 'Error: ' // trim(msg)
        error stop 1
    end if

    method = trim(scheme_output_name(scheme))

    call procedure()
    call output_exact_standard_func()
    call output_order()
    call output_csv()

    if (output_fig_flag) then
        call calc_init(ini_mesh_num)
        call interpolate(ini_mesh_num)
        call output_fig(ini_mesh_num)
    end if

    write(*, '(a)') 'scheme: ' // trim(method)
    write(*, '(a,es12.4)') 'clock ticks used: ', real(time_used, precision)
    write(*, '(a)') 'log: ' // trim(method) // '_log.dat'
    write(*, '(a)') 'csv: results/convergence/' // trim(method) // '_' // trim(func_label()) // '.csv'
end program main
