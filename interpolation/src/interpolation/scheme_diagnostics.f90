module scheme_diagnostics_mod
    use precision_select
    implicit none

    ! Optional diagnostics returned by a scheme when the caller requests them.
    ! Fields that are not produced by the active scheme stay zero.
    type :: scheme_diagnostics
        real(precision) :: beta(3) = 0.0_precision
        real(precision) :: omega(3) = 0.0_precision
        real(precision) :: tau = 0.0_precision
        real(precision) :: beta_L = 0.0_precision
        real(precision) :: beta_N = 0.0_precision
        real(precision) :: omega_L = 0.0_precision
        real(precision) :: omega_N = 0.0_precision
        real(precision) :: lambda_h2 = 0.0_precision
        real(precision) :: h_noscale = 0.0_precision
    end type scheme_diagnostics

end module scheme_diagnostics_mod
