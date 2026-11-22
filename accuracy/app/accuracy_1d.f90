program accuracy_1d
    ! Driver for the 1D periodic convection accuracy test.
    !
    ! Reads configuration and the selected initial profile, refines the mesh
    ! from nbase to nbase*2^(levels-1), and writes the error/order table.
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use scheme_params, only: scheme, read_scheme_namelist, validate_scheme_params
    use scheme_dispatch, only: scheme_output_name
    use accuracy_control
    use accuracy_order
    use c1d_exact_mod, only: c1d_exact
    use c1d_solver, only: c1d_solve_level
    implicit none

    integer :: ierr, u_unit
    character(len=256) :: msg
    character(len=32) :: op
    character(len=64) :: initial, scheme_name
    character(len=256) :: out_name
    character(len=512) :: csv_path, log_path
    integer :: k, n
    real(precision) :: cpu, e1, e2, einf
    real(precision), allocatable :: uend(:), uex(:)

#ifdef PRECISION_DP
    character(len=*), parameter :: precision_label = 'dp'
#else
    character(len=*), parameter :: precision_label = 'qp'
#endif

    ! ---------------------------------------------------------------
    ! Configuration
    ! ---------------------------------------------------------------
    call accuracy_control_init_1d('input.nml', ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    open (newunit=u_unit, file='input.nml', status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write (*, '(a)') 'error: cannot open input.nml'
        error stop 1
    end if
    call read_scheme_namelist(u_unit, ierr)
    close (u_unit)
    if (ierr /= 0) then
        write (*, '(a)') 'error: cannot read &scheme_nml'
        error stop 1
    end if

    call validate_scheme_params(ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call accuracy_resolve_operator(scheme, op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call accuracy_ensure_output_dirs(output_dir)

    initial = trim(initial_1d)
    scheme_name = trim(scheme_output_name(scheme))

    ! ---------------------------------------------------------------
    ! Refinement loop
    ! ---------------------------------------------------------------
    call order_reset()

    do k = 1, int(levels, 4)
        n = int(nbase, 4)*2**(k - 1)

        allocate (uend(0:n), uex(0:n))

        call c1d_solve_level(n, initial, op, cfl, tend, uend, cpu)
        call c1d_exact(n, initial, tend, uex)
        call error_norms_1d(n, uend, uex, e1, e2, einf)

        if (.not. ieee_is_finite(e1) .or. .not. ieee_is_finite(e2) .or. &
            .not. ieee_is_finite(einf)) then
            write (*, '(a,i0)') 'error: non-finite error at n = ', n
            error stop 1
        end if

        call order_record(k, int(n, integer_precision), e1, e2, einf, cpu)
        deallocate (uend, uex)
    end do

    ! ---------------------------------------------------------------
    ! Output
    ! ---------------------------------------------------------------
    out_name = 'accuracy_1d_' // trim(initial) // '_' // trim(scheme_name) // '_' // trim(op)
    csv_path = trim(output_dir) // '/csv/' // trim(out_name) // '.csv'
    log_path = trim(output_dir) // '/logs/' // trim(out_name) // '.log'

    call order_write_csv(csv_path, scheme_name, op, ierr, msg)
    if (ierr /= 0) then
        write (*, '(a)') 'error: ' // trim(msg)
        error stop 1
    end if

    call write_log(log_path, initial, scheme_name, op, precision_label)

contains

    subroutine write_log(path, initial_name, method_name, operator_name, prec)
        ! Write run metadata as key=value lines.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: initial_name, method_name, operator_name, prec

        integer :: unit

        open (newunit=unit, file=trim(path), status='replace', action='write')
        write (unit, '(a)') 'initial=' // trim(initial_name)
        write (unit, '(a)') 'scheme=' // trim(method_name)
        write (unit, '(a)') 'operator=' // trim(operator_name)
        write (unit, '(a)') 'precision=' // trim(prec)
        write (unit, '(a,i0)') 'nbase=', nbase
        write (unit, '(a,i0)') 'levels=', levels
        write (unit, '(a,es24.16)') 'cfl=', cfl
        write (unit, '(a,es24.16)') 'tend=', tend
        close (unit)
    end subroutine write_log

end program accuracy_1d
