module difference_operators
    ! Selectable WCNS difference operators for the complex-wave test.
    !
    ! The WCNS interpolation supplies half-node fluxes E_{j+1/2}.  These
    ! operators turn them (and, for the hybrid forms, the integer-node
    ! values u_j) into the spatial right-hand side K_j = -dE/dx at node j.
    !
    ! Array layout used by every operator:
    !   flux_half(1:6) = E_{j-5/2}, E_{j-3/2}, E_{j-1/2}, E_{j+1/2},
    !                    E_{j+3/2}, E_{j+5/2}
    !   u_node(1:5)    = u_{j-2}, u_{j-1}, u_j, u_{j+1}, u_{j+2}
    ! so flux_half(4) - flux_half(3) is E_{j+1/2} - E_{j-1/2}, and
    ! u_node(4) - u_node(2) is u_{j+1} - u_{j-1}.
    use precision_select
    implicit none

contains

    subroutine difference_apply(operator, h, flux_half, u_node, dfdx)
        ! Dispatch one of WCNSE2, WCNSE4, WCNSE6, HCS4, or HCS6.
        character(len=*), intent(in) :: operator
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(in) :: u_node(1:5)
        real(precision), intent(out) :: dfdx

        select case (trim(operator))
        case ('WCNSE2')
            call wcnse2(h, flux_half, dfdx)
        case ('WCNSE4')
            call wcnse4(h, flux_half, dfdx)
        case ('WCNSE6')
            call wcnse6(h, flux_half, dfdx)
        case ('HCS4')
            call hcs4(h, flux_half, u_node, dfdx)
        case ('HCS6')
            call hcs6(h, flux_half, u_node, dfdx)
        case default
            write(*, '(a)') 'unknown difference operator: ' // trim(operator)
            error stop 1
        end select
    end subroutine difference_apply

    subroutine wcnse2(h, flux_half, dfdx)
        ! Second-order flux difference using one pair of half nodes.
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(out) :: dfdx

        dfdx = -(flux_half(4) - flux_half(3))/h
    end subroutine wcnse2

    subroutine wcnse4(h, flux_half, dfdx)
        ! Fourth-order flux difference using two pairs of half nodes.
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(out) :: dfdx

        dfdx = -((9.0_precision/8.0_precision)*(flux_half(4) - flux_half(3)) &
                 - (1.0_precision/24.0_precision)*(flux_half(5) - flux_half(2)))/h
    end subroutine wcnse4

    subroutine wcnse6(h, flux_half, dfdx)
        ! Sixth-order flux difference using three pairs of half nodes.
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(out) :: dfdx

        dfdx = -((75.0_precision/64.0_precision)*(flux_half(4) - flux_half(3)) &
                 - (25.0_precision/384.0_precision)*(flux_half(5) - flux_half(2)) &
                 + (3.0_precision/640.0_precision)*(flux_half(6) - flux_half(1)))/h
    end subroutine wcnse6

    subroutine hcs4(h, flux_half, u_node, dfdx)
        ! Hybrid compact fourth-order operator: one pair of half nodes plus
        ! one pair of integer nodes.
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(in) :: u_node(1:5)
        real(precision), intent(out) :: dfdx

        real(precision), parameter :: a = 4.0_precision/3.0_precision
        real(precision), parameter :: b = -1.0_precision/6.0_precision

        dfdx = a*(flux_half(4) - flux_half(3)) + b*(u_node(4) - u_node(2))
        dfdx = -dfdx/h
    end subroutine hcs4

    subroutine hcs6(h, flux_half, u_node, dfdx)
        ! Hybrid compact sixth-order operator: one pair of half nodes plus
        ! two pairs of integer nodes.  The coefficients follow from matching
        ! the Taylor series through O(h^6):
        !   a + 2b + 4c = 1,  a/24 + b/3 + 8c/3 = 0,
        !   a/1920 + b/60 + 8c/15 = 0
        ! giving a = 64/45, b = -2/9, c = 1/180.
        real(precision), intent(in) :: h
        real(precision), intent(in) :: flux_half(1:6)
        real(precision), intent(in) :: u_node(1:5)
        real(precision), intent(out) :: dfdx

        real(precision), parameter :: a = 64.0_precision/45.0_precision
        real(precision), parameter :: b = -2.0_precision/9.0_precision
        real(precision), parameter :: c = 1.0_precision/180.0_precision

        dfdx = a*(flux_half(4) - flux_half(3)) &
             + b*(u_node(4) - u_node(2)) &
             + c*(u_node(5) - u_node(1))
        dfdx = -dfdx/h
    end subroutine hcs6

end module difference_operators
