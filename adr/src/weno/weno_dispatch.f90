module weno_dispatch
    ! Runtime dispatch for WENO reconstruction and the WENO ADR residual.
    use precision_select
    use adr_control, only: PI
    use weno_params, only: scheme
    use weno_schemes, only: weno5_linear_reconstruct, weno3_linear, &
                            weno3_js, weno5_js, weno5_m, weno5_z, &
                            weno3_trig, weno3_llr, baeza_weno3, weno_c_deheilan, &
                            weno5_dlr, weno3_ga_rbf_eps_define, weno3_ga_rbf_eps_undefine, &
                            rbf_5_linear, teno_original_333, teno_original_333_rbf, &
                            teno_mody_semi_eno, teno_mody_semi_eno_rbf_ga, get_thinc
    implicit none

contains

    subroutine weno_reconstruct(u, flux_half, h)
        ! Reconstruct the half-node flux for the selected WENO scheme.
        real(precision), intent(in) :: u(-4:4)
        real(precision), intent(out) :: flux_half
        real(precision), intent(in) :: h

        select case (trim(scheme))
        case ('WENO5_LINEAR')
            call weno5_linear_reconstruct(u(-2:3), flux_half)
        case ('WENO3_LINEAR')
            call weno3_linear(u(-2:3), flux_half)
        case ('WENO3_JS')
            call weno3_js(u(-2:3), flux_half)
        case ('WENO5_JS')
            call weno5_js(u(-2:3), flux_half)
        case ('WENO5_M')
            call weno5_m(u(-2:3), flux_half)
        case ('WENO5_Z')
            call weno5_z(u(-2:3), flux_half)
        case ('WENO3_TRIG')
            call weno3_trig(u(-2:3), flux_half, h)
        case ('WENO3_LLR')
            call weno3_llr(u(-2:3), flux_half, h)
        case ('BAEZA_WENO3')
            call baeza_weno3(u(-2:3), flux_half)
        case ('WENO_C')
            call weno_c_deheilan(u(-2:3), flux_half)
        case ('WENO5_DLR')
            call weno5_dlr(u, flux_half, h)
        case ('WENO3_GA_RBF_DEFINE')
            call weno3_ga_rbf_eps_define(u(-2:3), flux_half, h)
        case ('WENO3_GA_RBF_UNDEFINE')
            call weno3_ga_rbf_eps_undefine(u(-2:3), flux_half, h)
        case ('RBF5_LINEAR')
            call rbf_5_linear(u(-2:3), flux_half, h)
        case ('TENO_333')
            call teno_original_333(u(-2:3), flux_half, h)
        case ('TENO_333_RBF')
            call teno_original_333_rbf(u(-2:3), flux_half, h)
        case ('TENO_MOD_SEMI_ENO')
            call teno_mody_semi_eno(u(-2:3), flux_half, h)
        case ('TENO_MOD_SEMI_ENO_RBF')
            call teno_mody_semi_eno_rbf_ga(u(-2:3), flux_half, h)
        case ('THINC')
            call get_thinc(u(-2:3), flux_half)
        case default
            write(*, '(a)') 'unknown WENO scheme: ' // trim(scheme)
            error stop 1
        end select
    end subroutine weno_reconstruct

    subroutine weno_rhs(nx, halo, u, f)
        ! Second-order flux difference of the WENO-reconstructed flux.
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: u(-halo:nx+halo)
        real(precision), intent(out) :: f(-halo:nx+halo)

        real(precision) :: flux(-3:nx+2)
        real(precision) :: h
        integer(integer_precision) :: k, j

        h = 2.0_precision * PI / real(nx, precision)

        flux = 0.0_precision
        do k = -3, nx + 2
            call weno_reconstruct(u(k-4:k+4), flux(k), h)
        end do

        f = 0.0_precision
        do j = 0, nx
            f(j) = -(flux(j) - flux(j-1)) / h
        end do
    end subroutine weno_rhs

end module weno_dispatch
