program test_precision
    use precision_select
    use test_util
    implicit none

#ifdef PRECISION_DP
    call check_true(precision == dp, 'PRECISION_DP selects dp')
    call check_true(storage_size(real(0, precision)) >= 64, 'dp is at least 64-bit')
#else
    call check_true(precision == qp, 'default selects qp')
    call check_true(storage_size(real(0, precision)) >= 128, 'qp is at least 128-bit')
#endif
    call finish()
end program test_precision
