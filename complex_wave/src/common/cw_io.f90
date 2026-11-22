module cw_io
    ! Snapshot, Tecplot, and log output for the complex-wave driver.
    !
    ! Each run writes one CSV per snapshot, one multi-zone Tecplot file with
    ! all snapshots, and one key=value log with run metadata and extrema.
    use precision_select
    implicit none

contains

    subroutine ensure_output_dirs(output_dir)
        ! Create the result subdirectories used by this driver.
        character(len=*), intent(in) :: output_dir

        call execute_command_line('mkdir -p "' // trim(output_dir) // '/csv" "' // &
            trim(output_dir) // '/dat" "' // trim(output_dir) // '/logs" "' // &
            trim(output_dir) // '/figures"')
    end subroutine ensure_output_dirs

    subroutine write_csv_snapshot(path, nx, u)
        ! Write one snapshot as x,u columns.
        character(len=*), intent(in) :: path
        integer(integer_precision), intent(in) :: nx
        real(precision), intent(in) :: u(0:nx)

        integer :: unit
        integer(integer_precision) :: j
        real(precision) :: h, x

        h = 2.0_precision/real(nx, precision)

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(a)') 'x,u'
        do j = 0, nx
            x = -1.0_precision + real(j, precision)*h
            write(unit, '(es24.16,a,es24.16)') x, ',', u(j)
        end do
        close(unit)
    end subroutine write_csv_snapshot

    subroutine write_tecplot_snapshots(path, title, nx, n_snap, times, snaps)
        ! Write all snapshots as point-packed Tecplot zones.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: title
        integer(integer_precision), intent(in) :: nx, n_snap
        real(precision), intent(in) :: times(0:n_snap)
        real(precision), intent(in) :: snaps(0:nx, 0:n_snap)

        integer :: unit
        integer(integer_precision) :: j, k
        real(precision) :: h, x

        h = 2.0_precision/real(nx, precision)

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(a)') 'title="complex wave: ' // trim(title) // '"'
        write(unit, '(a)') 'variables="X","U"'
        do k = 0, n_snap
            write(unit, '(a,i0,a,es12.5,a)') 'Zone I=', nx + 1, &
                ', datapacking=POINT, T="t=', times(k), '"'
            do j = 0, nx
                x = -1.0_precision + real(j, precision)*h
                write(unit, '(es24.16,1x,es24.16)') x, snaps(j, k)
            end do
        end do
        close(unit)
    end subroutine write_tecplot_snapshots

    subroutine write_run_log(path, precision_label, scheme_name, operator_name, &
                             mesh_num, periods, cfl, dt_min, dt_max, steps, &
                             umin, umax, nan_count)
        ! Write run metadata and solution extrema as key=value lines.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: precision_label
        character(len=*), intent(in) :: scheme_name
        character(len=*), intent(in) :: operator_name
        integer(integer_precision), intent(in) :: mesh_num, periods
        real(precision), intent(in) :: cfl
        real(precision), intent(in) :: dt_min, dt_max
        integer(integer_precision), intent(in) :: steps
        real(precision), intent(in) :: umin, umax
        integer(integer_precision), intent(in) :: nan_count

        integer :: unit

        open(newunit=unit, file=trim(path), status='replace', action='write')
        write(unit, '(a)') 'version=1'
        write(unit, '(a)') 'precision=' // trim(precision_label)
        write(unit, '(a)') 'scheme=' // trim(scheme_name)
        write(unit, '(a)') 'operator=' // trim(operator_name)
        write(unit, '(a,i0)') 'mesh_num=', mesh_num
        write(unit, '(a,i0)') 'periods=', periods
        write(unit, '(a,es24.16)') 'cfl=', cfl
        write(unit, '(a,es24.16)') 'dt_min=', dt_min
        write(unit, '(a,es24.16)') 'dt_max=', dt_max
        write(unit, '(a,i0)') 'steps=', steps
        write(unit, '(a,es24.16)') 'umin=', umin
        write(unit, '(a,es24.16)') 'umax=', umax
        write(unit, '(a,i0)') 'nan_count=', nan_count
        close(unit)
    end subroutine write_run_log

end module cw_io
