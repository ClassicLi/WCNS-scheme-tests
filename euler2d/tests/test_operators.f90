program test_operators
    ! Polynomial-exactness tests for the selectable WCNS difference operators.
    !
    ! difference_apply returns the residual -dE/dx, so a polynomial f of
    ! degree p is differentiated exactly when the reported value matches
    ! -p*x0**(p-1).  Each operator must be exact through its design order
    ! and inexact at the next degree.
    use precision_select
    use difference_operators, only: difference_apply
    use test_util
    implicit none

    integer :: p

    ! WCNSE2 is exact for degrees 1 and 2, but not degree 3.
    call check_poly('WCNSE2', 1, .true.)
    call check_poly('WCNSE2', 2, .true.)
    call check_poly('WCNSE2', 3, .false.)

    ! Fourth-order operators are exact through degree 4, but not degree 5.
    do p = 1, 4
        call check_poly('WCNSE4', p, .true.)
        call check_poly('HCS4', p, .true.)
    end do
    call check_poly('WCNSE4', 5, .false.)
    call check_poly('HCS4', 5, .false.)

    ! Sixth-order operators are exact through degree 6, but not degree 7.
    do p = 1, 6
        call check_poly('WCNSE6', p, .true.)
        call check_poly('HCS6', p, .true.)
    end do
    call check_poly('WCNSE6', 7, .false.)
    call check_poly('HCS6', 7, .false.)

    call finish()

contains

    subroutine check_poly(op, power, exact)
        ! Check one operator on f(x) = x**power sampled at nodes and half nodes.
        character(len=*), intent(in) :: op
        integer, intent(in) :: power
        logical, intent(in) :: exact
        real(precision) :: fh(1:6), un(1:5), dfdx, expected, err
        real(precision), parameter :: x0 = 0.3_precision
        real(precision), parameter :: hh = 0.1_precision
        integer :: k

        do k = 1, 6
            fh(k) = (x0 + (real(k, precision) - 3.5_precision)*hh)**power
        end do
        do k = 1, 5
            un(k) = (x0 + real(k - 3, precision)*hh)**power
        end do
        call difference_apply(op, hh, fh, un, dfdx)
        ! The operator returns the residual -dE/dx.
        expected = -real(power, precision)*x0**(power - 1)
        err = abs(dfdx - expected)
        if (exact) then
            call check_true(err < 1.0e-10_precision, trim(op)//' exact p='//achar(48 + power))
        else
            call check_true(err > 1.0e-8_precision, trim(op)//' not exact p='//achar(48 + power))
        end if
    end subroutine check_poly

end program test_operators
