module interp_qr
    use precision_select
    use scheme_params, only: qr_smv
    use scheme_diagnostics_mod, only: scheme_diagnostics
    implicit none

contains

    ! Quasi-rational interpolation on the three-point stencil.
    ! QR does not produce nonlinear weights, so diag is reset to zero.
    subroutine quasi_rational(u, u_half, diag)
        real(precision), intent(in) :: u(-3:4)
        real(precision), intent(out) :: u_half
        type(scheme_diagnostics), intent(out), optional :: diag
        real(precision) :: numerator, denominator

        if (present(diag)) diag = scheme_diagnostics()

        numerator = real(3, precision) * abs(u(0) - u(-1)) * (u(1) - u(0)) &
                  + abs(u(1) - u(0)) * (u(0) - u(-1))
        denominator = abs(u(1) - u(0)) + real(3, precision) * abs(u(0) - u(-1))

        u_half = u(0) + numerator / (real(2, precision) * denominator + qr_smv)
    end subroutine quasi_rational

end module interp_qr
