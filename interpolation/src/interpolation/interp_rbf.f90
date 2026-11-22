module interp_rbf
    use precision_select
    use interp_linear, only: linear_3rd_minus1_to_plus1, linear_5th_minus2_to_plus2
    use interp_deriv
    use scheme_params, only: rbf_num_deriv, rbf_denom_deriv, rbf_floor_value
    use scheme_diagnostics_mod, only: scheme_diagnostics
    implicit none

    real(precision), parameter :: rbf_eps = 1.0e-80_precision

contains

    subroutine RBF_GA_Perturbed_3(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-1:1)

        call linear_3rd_minus1_to_plus1(u, u_half)

        lambda_h2 = -real(1, precision) / real(6, precision) &
            * select_deriv_3rd(u, rbf_num_deriv) / (select_deriv_1st(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-1) = -real(3, precision) / real(32, precision) * lambda_h2
        temp(0) = -real(3, precision) / real(16, precision) * lambda_h2
        temp(1) = real(9, precision) / real(32, precision) * lambda_h2

        u_half = u_half + temp(-1) * u(-1) + temp(0) * u(0) + temp(1) * u(1)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_GA_Perturbed_3

    subroutine RBF_MQ_Perturbed_3(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-1:1)

        call linear_3rd_minus1_to_plus1(u, u_half)

        lambda_h2 = -real(1, precision) / real(3, precision) &
            * select_deriv_3rd(u, rbf_num_deriv) / (select_deriv_1st(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-1) = 0.0_precision
        temp(0) = -real(3, precision) / real(16, precision) * lambda_h2
        temp(1) = real(3, precision) / real(16, precision) * lambda_h2

        u_half = u_half + temp(-1) * u(-1) + temp(0) * u(0) + temp(1) * u(1)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_MQ_Perturbed_3

    subroutine RBF_IQ_Perturbed_3(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-1:1)

        call linear_3rd_minus1_to_plus1(u, u_half)

        lambda_h2 = -real(1, precision) / real(12, precision) &
            * select_deriv_3rd(u, rbf_num_deriv) / (select_deriv_1st(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-1) = -real(9, precision) / real(80, precision) * lambda_h2
        temp(0) = -real(21, precision) / real(40, precision) * lambda_h2
        temp(1) = real(51, precision) / real(80, precision) * lambda_h2

        u_half = u_half + temp(-1) * u(-1) + temp(0) * u(0) + temp(1) * u(1)

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_IQ_Perturbed_3

    subroutine RBF_GA_Perturbed_5(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-2:2)
        integer :: loop

        call linear_5th_minus2_to_plus2(u, u_half)

        lambda_h2 = -real(1, precision) / real(20, precision) &
            * select_deriv_5th(u, rbf_num_deriv) / (select_deriv_3rd(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-2) = real(45, precision) / real(512, precision) * lambda_h2
        temp(-1) = -real(15, precision) / real(128, precision) * lambda_h2
        temp(0) = -real(45, precision) / real(256, precision) * lambda_h2
        temp(1) = real(45, precision) / real(128, precision) * lambda_h2
        temp(2) = -real(75, precision) / real(512, precision) * lambda_h2

        do loop = -2, 2
            u_half = u_half + temp(loop) * u(loop)
        end do

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_GA_Perturbed_5

    subroutine RBF_MQ_Perturbed_5(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-2:2)
        integer :: loop

        call linear_5th_minus2_to_plus2(u, u_half)

        lambda_h2 = -real(1, precision) / real(40, precision) &
            * select_deriv_5th(u, rbf_num_deriv) / (select_deriv_3rd(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-2) = real(1125, precision) / real(7168, precision) * lambda_h2
        temp(-1) = -real(285, precision) / real(1792, precision) * lambda_h2
        temp(0) = -real(1665, precision) / real(3584, precision) * lambda_h2
        temp(1) = real(1395, precision) / real(1792, precision) * lambda_h2
        temp(2) = -real(2235, precision) / real(7168, precision) * lambda_h2

        do loop = -2, 2
            u_half = u_half + temp(loop) * u(loop)
        end do

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_MQ_Perturbed_5

    subroutine RBF_IQ_Perturbed_5(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        real(precision) :: lambda_h2
        real(precision) :: temp(-2:2)
        integer :: loop

        call linear_5th_minus2_to_plus2(u, u_half)

        lambda_h2 = -real(3, precision) / real(220, precision) &
            * select_deriv_5th(u, rbf_num_deriv) / (select_deriv_3rd(u, rbf_denom_deriv) + rbf_eps)

        if (abs(lambda_h2) >= rbf_floor_value) lambda_h2 = 0.0_precision

        temp(-2) = real(22775, precision) / real(76288, precision) * lambda_h2
        temp(-1) = -real(6385, precision) / real(19072, precision) * lambda_h2
        temp(0) = -real(30015, precision) / real(38144, precision) * lambda_h2
        temp(1) = real(26395, precision) / real(19072, precision) * lambda_h2
        temp(2) = -real(42785, precision) / real(76288, precision) * lambda_h2

        do loop = -2, 2
            u_half = u_half + temp(loop) * u(loop)
        end do

        if (present(diag)) then
            diag = scheme_diagnostics()
            diag%lambda_h2 = lambda_h2
        end if
    end subroutine RBF_IQ_Perturbed_5

    function select_deriv_1st(u, deriv_type) result(d)
        real(precision), intent(in) :: u(-3:4)
        character(len=*), intent(in) :: deriv_type
        real(precision) :: d

        select case (trim(deriv_type))
        case ('left')
            d = get_derivative_1st_tr_1_left(u)
        case ('mid')
            d = get_derivative_1st_tr_2(u)
        case ('right')
            d = get_derivative_1st_tr_1_right(u)
        case default
            d = get_derivative_1st_tr_1_right(u)
        end select
    end function select_deriv_1st

    function select_deriv_3rd(u, deriv_type) result(d)
        real(precision), intent(in) :: u(-3:4)
        character(len=*), intent(in) :: deriv_type
        real(precision) :: d

        select case (trim(deriv_type))
        case ('left')
            d = get_derivative_3rd_tr_1_left(u)
        case ('mid')
            d = get_derivative_3rd_tr_2(u)
        case ('right')
            d = get_derivative_3rd_tr_1_right(u)
        case default
            d = get_derivative_3rd_tr_1_right(u)
        end select
    end function select_deriv_3rd

    function select_deriv_5th(u, deriv_type) result(d)
        real(precision), intent(in) :: u(-3:4)
        character(len=*), intent(in) :: deriv_type
        real(precision) :: d

        select case (trim(deriv_type))
        case ('left')
            d = get_derivative_5th_tr_1_left(u)
        case ('mid')
            d = get_derivative_5th_tr_2(u)
        case ('right')
            d = get_derivative_5th_tr_1_right(u)
        case default
            d = get_derivative_5th_tr_1_right(u)
        end select
    end function select_deriv_5th

end module interp_rbf
