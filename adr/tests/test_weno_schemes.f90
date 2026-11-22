program test_weno_schemes
    ! Constant, linear, and finiteness checks for the retained WENO schemes.
    use, intrinsic :: ieee_arithmetic
    use precision_select
    use weno_params, only: scheme
    use weno_dispatch, only: weno_reconstruct
    use test_util
    implicit none

    real(precision), parameter :: h = 0.1_precision
    real(precision) :: u(-4:4), flux
    integer :: j, k
    character(len=24) :: names(9)
    character(len=24) :: names_finite(5)

    names = [character(len=24) :: 'WENO3_JS', 'WENO5_JS', 'WENO5_M', 'WENO5_Z', &
             'WENO_C', 'BAEZA_WENO3', 'WENO3_LLR', 'TENO_333', 'TENO_MOD_SEMI_ENO']

    do k = 1, size(names)
        scheme = trim(names(k))

        u = 2.5_precision
        call weno_reconstruct(u, flux, h)
        call check_close(flux, 2.5_precision, test_tol, trim(names(k)) // ' constant')

        do j = -4, 4
            u(j) = real(j, precision)
        end do
        call weno_reconstruct(u, flux, h)
        call check_close(flux, 0.5_precision, test_tol, trim(names(k)) // ' linear exact')

        do j = -4, 4
            u(j) = cos(real(j, precision))
        end do
        call weno_reconstruct(u, flux, h)
        call check_true(ieee_is_finite(flux), trim(names(k)) // ' finite')
    end do

    ! These variants use trigonometric or RBF normalizations that are not
    ! polynomial-exact at finite h; only require finite results.
    names_finite = [character(len=24) :: 'WENO3_TRIG', 'WENO3_GA_RBF_DEFINE', &
                    'WENO3_GA_RBF_UNDEFINE', 'RBF5_LINEAR', 'TENO_MOD_SEMI_ENO_RBF']

    do k = 1, size(names_finite)
        scheme = trim(names_finite(k))
        do j = -4, 4
            u(j) = cos(real(j, precision))
        end do
        call weno_reconstruct(u, flux, h)
        call check_true(ieee_is_finite(flux), trim(names_finite(k)) // ' finite')
    end do

    ! WENO5_DLR is kept as an experimental subroutine. The original file marks
    ! the 5th-order DLR as problematic and it can produce NaN on ADR data, so
    ! this test only exercises the dispatch branch and makes no numeric claim.
    scheme = 'WENO5_DLR'
    do j = -4, 4
        u(j) = cos(real(j, precision))
    end do
    call weno_reconstruct(u, flux, h)

    ! TENO_333_RBF keeps the original uninitialized eps2h2 (unvalidated), and
    ! THINC's original expression is dimensionally suspect; exercise both
    ! branches without a numeric claim.
    scheme = 'TENO_333_RBF'
    call weno_reconstruct(u, flux, h)
    scheme = 'THINC'
    call weno_reconstruct(u, flux, h)

    call finish()
end program test_weno_schemes
