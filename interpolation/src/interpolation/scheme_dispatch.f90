module scheme_dispatch
    use precision_select
    use scheme_params, only: scheme, scheme_allowed, moc33_linear, moc53_linear
    use scheme_diagnostics_mod, only: scheme_diagnostics
    use interp_linear, only: linear_3rd_minus1_to_plus1, linear_3rd_minus2_to_0, &
                             linear_3rd_0_to_plus2, linear_4th_minus2_to_plus1, &
                             linear_4th_minus1_to_plus2, linear_5th_minus2_to_plus2, &
                             linear_6th_minus2_to_plus3
    use interp_qr, only: quasi_rational
    use interp_rbf, only: RBF_GA_Perturbed_3, RBF_MQ_Perturbed_3, RBF_IQ_Perturbed_3, &
                          RBF_GA_Perturbed_5, RBF_MQ_Perturbed_5, RBF_IQ_Perturbed_5
    use interp_traditional, only: muscl_2nd, JS3_interp, Z3_interp, JS5_interp, Z5_interp
    use interp_framework, only: MOC322, MOC33, MOC53, MOZ
    implicit none

contains

    subroutine interpolate_dispatch(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag

        if (present(diag)) diag = scheme_diagnostics()

        select case (trim(scheme))
        case ('MUSCL')
            call muscl_2nd(u, u_half, diag)
        case ('JS3')
            call JS3_interp(u, u_half, diag)
        case ('Z3')
            call Z3_interp(u, u_half, diag)
        case ('JS5')
            call JS5_interp(u, u_half, diag)
        case ('Z5')
            call Z5_interp(u, u_half, diag)
        case ('QR')
            call quasi_rational(u, u_half, diag)
        case ('MOC322')
            call MOC322(u, u_half, diag)
        case ('MOC33')
            call MOC33(u, u_half, diag)
        case ('MOC53')
            call MOC53(u, u_half, diag)
        case ('MOZ')
            call MOZ(u, u_half, diag)
        case ('LINEAR3_M1P1')
            call linear_3rd_minus1_to_plus1(u, u_half)
        case ('LINEAR3_M2P0')
            call linear_3rd_minus2_to_0(u, u_half)
        case ('LINEAR3_0P2')
            call linear_3rd_0_to_plus2(u, u_half)
        case ('LINEAR4_M2P1')
            call linear_4th_minus2_to_plus1(u, u_half)
        case ('LINEAR4_M1P2')
            call linear_4th_minus1_to_plus2(u, u_half)
        case ('LINEAR5_M2P2')
            call linear_5th_minus2_to_plus2(u, u_half)
        case ('LINEAR6_M2P3')
            call linear_6th_minus2_to_plus3(u, u_half)
        case ('RBF_GA3')
            call RBF_GA_Perturbed_3(u, u_half, diag)
        case ('RBF_MQ3')
            call RBF_MQ_Perturbed_3(u, u_half, diag)
        case ('RBF_IQ3')
            call RBF_IQ_Perturbed_3(u, u_half, diag)
        case ('RBF_GA5')
            call RBF_GA_Perturbed_5(u, u_half, diag)
        case ('RBF_MQ5')
            call RBF_MQ_Perturbed_5(u, u_half, diag)
        case ('RBF_IQ5')
            call RBF_IQ_Perturbed_5(u, u_half, diag)
        case default
            write(*, '(a)') 'unknown scheme: ' // trim(scheme)
            error stop 1
        end select
    end subroutine interpolate_dispatch

    function scheme_known(name) result(ok)
        character(len=*), intent(in) :: name
        logical :: ok

        ok = scheme_allowed(name)
    end function scheme_known

    function scheme_order(name) result(p)
        character(len=*), intent(in) :: name
        integer :: p

        select case (trim(name))
        case ('MUSCL')
            p = 2
        case ('JS3', 'Z3', 'QR', 'MOC322', 'MOC33', &
              'LINEAR3_M1P1', 'LINEAR3_M2P0', 'LINEAR3_0P2', &
              'RBF_GA3', 'RBF_MQ3', 'RBF_IQ3')
            p = 3
        case ('LINEAR4_M2P1', 'LINEAR4_M1P2')
            p = 4
        case ('JS5', 'Z5', 'MOC53', 'MOZ', 'LINEAR5_M2P2', &
              'RBF_GA5', 'RBF_MQ5', 'RBF_IQ5')
            p = 5
        case ('LINEAR6_M2P3')
            p = 6
        case default
            p = 0
        end select
    end function scheme_order

    function scheme_output_name(name) result(str)
        character(len=*), intent(in) :: name
        character(len=64) :: str

        select case (trim(name))
        case ('MOC33')
            str = trim(name) // '_' // trim(moc33_linear)
        case ('MOC53')
            str = trim(name) // '_' // trim(moc53_linear)
        case default
            str = trim(name)
        end select
    end function scheme_output_name

end module scheme_dispatch
