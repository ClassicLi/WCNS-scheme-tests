module scheme_params
    use precision_select
    implicit none

    character(len=32) :: scheme = 'LINEAR3_M1P1'

    ! JS- and Z-type nonlinear weighted interpolation
    real(precision) :: js3_eps = 1.0e-6_precision
    real(precision) :: js3_p = 2.0_precision
    real(precision) :: z3_eps = 1.0e-40_precision
    real(precision) :: z3_q = 2.0_precision
    real(precision) :: js5_eps = 1.0e-6_precision
    real(precision) :: js5_p = 2.0_precision
    real(precision) :: z5_eps = 1.0e-40_precision
    real(precision) :: z5_q = 2.0_precision

    ! Quasi-rational interpolation
    real(precision) :: qr_smv = 1.0e-80_precision

    ! RBF perturbed interpolation
    integer(integer_precision) :: rbf_size = 5_integer_precision
    character(len=8) :: rbf_kernel = 'GA'
    character(len=8) :: rbf_num_deriv = 'right'
    character(len=8) :: rbf_denom_deriv = 'right'
    real(precision) :: rbf_floor_value = 30.0_precision

    ! MOC322: main 3rd-order candidate plus two 2nd-order candidates
    real(precision) :: moc322_eps = 1.0e-80_precision
    integer(integer_precision) :: moc322_M_q = 2_integer_precision

    ! MOC33: three-point linear/RBF candidate plus three-point QR candidate
    real(precision) :: moc33_eps = 1.0e-80_precision
    integer(integer_precision) :: moc33_p = 2_integer_precision
    character(len=8) :: moc33_linear = 'GA'

    ! MOC53: five-point linear/RBF candidate plus three-point QR candidate
    real(precision) :: moc53_eps = 1.0e-80_precision
    integer(integer_precision) :: moc53_N_q = 4_integer_precision
    integer(integer_precision) :: moc53_p = 2_integer_precision
    character(len=8) :: moc53_linear = 'LINEAR'

    ! MOZ: improved five-point Z-type weighting
    real(precision) :: moz_eps = 1.0e-80_precision
    integer(integer_precision) :: moz_N_q = 4_integer_precision
    integer(integer_precision) :: moz_M_q = 2_integer_precision

    namelist /scheme_nml/ scheme, &
        js3_eps, js3_p, z3_eps, z3_q, js5_eps, js5_p, z5_eps, z5_q, &
        qr_smv, &
        rbf_size, rbf_kernel, rbf_num_deriv, rbf_denom_deriv, rbf_floor_value, &
        moc322_eps, moc322_M_q, &
        moc33_eps, moc33_p, moc33_linear, &
        moc53_eps, moc53_N_q, moc53_p, moc53_linear, &
        moz_eps, moz_N_q, moz_M_q

contains

    subroutine read_scheme_namelist(unit, ierr)
        integer, intent(in) :: unit
        integer, intent(out), optional :: ierr
        integer :: ios

        rewind(unit)
        read(unit, nml=scheme_nml, iostat=ios)
        if (present(ierr)) ierr = ios
    end subroutine read_scheme_namelist

    function scheme_allowed(name) result(ok)
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('MUSCL', 'JS3', 'Z3', 'JS5', 'Z5', 'QR', &
              'MOC322', 'MOC33', 'MOC53', 'MOZ', &
              'LINEAR3_M1P1', 'LINEAR3_M2P0', 'LINEAR3_0P2', &
              'LINEAR4_M2P1', 'LINEAR4_M1P2', &
              'LINEAR5_M2P2', 'LINEAR6_M2P3', &
              'RBF_GA3', 'RBF_MQ3', 'RBF_IQ3', &
              'RBF_GA5', 'RBF_MQ5', 'RBF_IQ5')
            ok = .true.
        case default
            ok = .false.
        end select
    end function scheme_allowed

    subroutine validate_scheme_params(ierr, msg)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. scheme_allowed(scheme)) then
            ierr = 1
            msg = 'unknown scheme: ' // trim(scheme)
            return
        end if

        select case (trim(rbf_kernel))
        case ('GA', 'MQ', 'IQ')
        case default
            ierr = 1
            msg = 'invalid rbf_kernel: ' // trim(rbf_kernel)
            return
        end select

        if (.not. deriv_type_valid(rbf_num_deriv)) then
            ierr = 1
            msg = 'invalid rbf_num_deriv: ' // trim(rbf_num_deriv)
            return
        end if

        if (.not. deriv_type_valid(rbf_denom_deriv)) then
            ierr = 1
            msg = 'invalid rbf_denom_deriv: ' // trim(rbf_denom_deriv)
            return
        end if

        if (rbf_size /= 3 .and. rbf_size /= 5) then
            ierr = 1
            msg = 'rbf_size must be 3 or 5'
            return
        end if

        if (.not. linear_choice_valid(moc33_linear)) then
            ierr = 1
            msg = 'invalid moc33_linear: ' // trim(moc33_linear)
            return
        end if

        if (.not. linear_choice_valid(moc53_linear)) then
            ierr = 1
            msg = 'invalid moc53_linear: ' // trim(moc53_linear)
            return
        end if

        if (moc322_M_q < 1 .or. moc33_p < 1 .or. moc53_N_q < 1 .or. &
            moc53_p < 1 .or. moz_N_q < 1 .or. moz_M_q < 1) then
            ierr = 1
            msg = 'integer scheme exponents must be positive'
            return
        end if

        if (.not. (js3_p > 0.0_precision) .or. .not. (js5_p > 0.0_precision) .or. &
            .not. (z3_q > 0.0_precision) .or. .not. (z5_q > 0.0_precision)) then
            ierr = 1
            msg = 'real scheme exponents must be positive'
            return
        end if

        if (js3_eps <= 0.0_precision .or. z3_eps <= 0.0_precision .or. &
            js5_eps <= 0.0_precision .or. z5_eps <= 0.0_precision .or. &
            qr_smv <= 0.0_precision .or. rbf_floor_value < 0.0_precision .or. &
            moc322_eps <= 0.0_precision .or. moc33_eps <= 0.0_precision .or. &
            moc53_eps <= 0.0_precision .or. moz_eps <= 0.0_precision) then
            ierr = 1
            msg = 'epsilon parameters must be positive'
            return
        end if
    end subroutine validate_scheme_params

    logical function deriv_type_valid(name)
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ('left', 'mid', 'right')
            deriv_type_valid = .true.
        case default
            deriv_type_valid = .false.
        end select
    end function deriv_type_valid

    logical function linear_choice_valid(name)
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ('GA', 'MQ', 'IQ', 'LINEAR')
            linear_choice_valid = .true.
        case default
            linear_choice_valid = .false.
        end select
    end function linear_choice_valid

end module scheme_params
