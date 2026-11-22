program test_params
    use precision_select
    use scheme_diagnostics_mod
    use scheme_params
    use test_util
    implicit none
    integer :: u, ierr
    character(len=128) :: msg

#ifdef PRECISION_DP
    call check_true(storage_size(real(0, precision)) >= 64, 'dp is at least 64-bit')
#else
    call check_true(storage_size(real(0, precision)) >= 128, 'qp is at least 128-bit')
#endif
    call check_close(real(moc322_M_q, precision), 2.0_precision, 0.0_precision, 'moc322_M_q default 2')
    call check_close(real(moc53_N_q, precision), 4.0_precision, 0.0_precision, 'moc53_N_q default 4')
    call check_close(js5_p, 2.0_precision, 0.0_precision, 'js5_p default 2')
    call check_true(trim(scheme) == 'LINEAR3_M1P1', 'default scheme')

    open(newunit=u, status='scratch', form='formatted')
    write(u, '(a)') "&scheme_nml"
    write(u, '(a)') "  scheme = 'Z5', z5_q = 3.0"
    write(u, '(a)') "/"
    rewind(u)
    call read_scheme_namelist(u)
    close(u)
    call check_true(trim(scheme) == 'Z5', 'namelist scheme override')
    call check_close(z5_q, 3.0_precision, 1.0e-14_precision, 'namelist z5_q override')

    call validate_scheme_params(ierr, msg)
    call check_true(ierr == 0, 'valid scheme params')

    scheme = 'NOT_A_SCHEME'
    call validate_scheme_params(ierr, msg)
    call check_true(ierr /= 0, 'invalid scheme rejected')

    scheme = 'Z3'
    z3_q = -2.0_precision
    call validate_scheme_params(ierr, msg)
    call check_true(ierr /= 0, 'negative z3_q rejected')
    z3_q = 2.0_precision
    call validate_scheme_params(ierr, msg)
    call check_true(ierr == 0, 'restored z3_q accepted')

    z5_q = 0.0_precision
    call validate_scheme_params(ierr, msg)
    call check_true(ierr /= 0, 'zero z5_q rejected')
    z5_q = 2.0_precision

    js3_p = -1.0_precision
    call validate_scheme_params(ierr, msg)
    call check_true(ierr /= 0, 'negative js3_p rejected')
    js3_p = 2.0_precision

    js5_p = 0.0_precision
    call validate_scheme_params(ierr, msg)
    call check_true(ierr /= 0, 'zero js5_p rejected')
    js5_p = 1.0_precision

    call validate_scheme_params(ierr, msg)
    call check_true(ierr == 0, 'restored real exponents accepted')

    call finish()
end program test_params
