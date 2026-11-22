module weno_params
    ! Scheme selection and tunable parameters for the WENO ADR driver.
    use precision_select
    implicit none

    character(len=32) :: scheme = 'WENO5_LINEAR'

    ! JS- and Z-type weighting parameters (defaults match the WCNS standard set).
    real(precision) :: js3_eps = 1.0e-6_precision
    real(precision) :: js3_p = 2.0_precision
    real(precision) :: js5_eps = 1.0e-6_precision
    real(precision) :: js5_p = 2.0_precision
    real(precision) :: z5_eps = 1.0e-40_precision
    real(precision) :: z5_q = 2.0_precision

    ! WENO-C parameters.
    real(precision) :: weno_c_c4 = 0.2_precision
    real(precision) :: weno_c_c5 = 0.1_precision
    real(precision) :: weno_c_eps = 1.0e-6_precision

    ! TENO parameters.
    real(precision) :: teno_q = 6.0_precision
    real(precision) :: teno_ct = 1.0e-4_precision

    ! THINC, RBF, and LLR parameters.
    real(precision) :: thinc_beta = 2.0_precision
    real(precision) :: rbf_eps = 75.0_precision
    real(precision) :: llr_eps = 1.0e-20_precision

    namelist /weno_nml/ scheme, js3_eps, js3_p, js5_eps, js5_p, z5_eps, z5_q, &
        weno_c_c4, weno_c_c5, weno_c_eps, teno_q, teno_ct, thinc_beta, rbf_eps, llr_eps

contains

    subroutine read_weno_namelist(unit, ierr)
        ! Read the &weno_nml group from an open input file.
        integer, intent(in) :: unit
        integer, intent(out), optional :: ierr

        integer :: ios

        rewind(unit)
        read(unit, nml=weno_nml, iostat=ios)
        if (present(ierr)) ierr = ios
    end subroutine read_weno_namelist

    function weno_scheme_allowed(name) result(ok)
        ! Return true for a known WENO scheme ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('WENO5_LINEAR', 'WENO3_LINEAR', 'WENO3_JS', 'WENO5_JS', 'WENO5_M', 'WENO5_Z', &
              'WENO_C', 'BAEZA_WENO3', 'WENO3_LLR', 'WENO3_TRIG', &
              'WENO5_DLR', 'WENO3_GA_RBF_DEFINE', 'WENO3_GA_RBF_UNDEFINE', 'RBF5_LINEAR', &
              'TENO_333', 'TENO_333_RBF', 'TENO_MOD_SEMI_ENO', 'TENO_MOD_SEMI_ENO_RBF', 'THINC')
            ok = .true.
        case default
            ok = .false.
        end select
    end function weno_scheme_allowed

    subroutine validate_weno_params(ierr, msg)
        ! Check the scheme ID, epsilon values, and exponents.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. weno_scheme_allowed(scheme)) then
            ierr = 1
            msg = 'unknown WENO scheme: ' // trim(scheme)
            return
        end if

        if (js3_eps <= 0.0_precision .or. js5_eps <= 0.0_precision .or. &
            z5_eps <= 0.0_precision .or. weno_c_eps <= 0.0_precision .or. &
            rbf_eps <= 0.0_precision .or. llr_eps <= 0.0_precision) then
            ierr = 2
            msg = 'WENO epsilon parameters must be positive'
            return
        end if

        if (js3_p <= 0.0_precision .or. js5_p <= 0.0_precision .or. &
            z5_q <= 0.0_precision .or. teno_q <= 0.0_precision) then
            ierr = 3
            msg = 'WENO exponents must be positive'
            return
        end if
    end subroutine validate_weno_params

end module weno_params
