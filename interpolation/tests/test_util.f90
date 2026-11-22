module test_util
    use precision_select
    implicit none

    integer :: n_pass = 0
    integer :: n_fail = 0

#ifdef PRECISION_DP
    real(precision), parameter :: test_tol = 1.0e-10_precision
#else
    real(precision), parameter :: test_tol = 1.0e-20_precision
#endif

contains

    subroutine check_true(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name

        if (cond) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*, '(a)') 'FAIL: ' // trim(name)
        end if
    end subroutine check_true

    subroutine check_close(a, b, tol, name)
        real(precision), intent(in) :: a, b, tol
        character(len=*), intent(in) :: name

        if (abs(a - b) <= tol) then
            n_pass = n_pass + 1
        else
            n_fail = n_fail + 1
            write(*, '(a,2(1x,es24.16))') 'FAIL: ' // trim(name), a, b
        end if
    end subroutine check_close

    subroutine check_order(observed, expected, tol, name)
        real(precision), intent(in) :: observed, expected, tol
        character(len=*), intent(in) :: name

        call check_close(observed, expected, tol, name)
    end subroutine check_order

    subroutine finish()
        write(*, '(a,i0,a,i0)') 'PASS: ', n_pass, ' checks, FAIL: ', n_fail
        if (n_fail > 0) error stop 1
    end subroutine finish

end module test_util
