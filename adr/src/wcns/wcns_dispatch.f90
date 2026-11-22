module wcns_dispatch
    ! WCNS residual: half-node interpolation plus the selected difference operator.
    use precision_select
    use adr_control, only: PI
    use scheme_dispatch, only: interpolate_dispatch
    use wcns_params, only: difference_operator
    use difference_operators, only: difference_apply
    implicit none

contains

    subroutine wcns_rhs(nx, halo, u, f)
        ! Interpolate to half nodes, then apply the selected difference operator.
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: u(-halo:nx+halo)
        real(precision), intent(out) :: f(-halo:nx+halo)

        real(precision) :: flux(-3:nx+2)
        real(precision) :: h
        integer(integer_precision) :: k, j

        h = 2.0_precision * PI / real(nx, precision)

        flux = 0.0_precision
        do k = -3, nx + 2
            call interpolate_dispatch(u(k-3:k+4), flux(k))
        end do

        f = 0.0_precision
        do j = 0, nx
            call difference_apply(difference_operator, h, flux(j-3:j+2), u(j-2:j+2), f(j))
        end do
    end subroutine wcns_rhs

end module wcns_dispatch
