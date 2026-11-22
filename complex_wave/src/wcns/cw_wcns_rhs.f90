module cw_wcns_rhs
    ! Spatial right-hand side of the complex-wave test.
    !
    ! The residual is assembled in two steps:
    !   1. WCNS interpolation gives the half-node fluxes E_{j+1/2};
    !   2. the selected difference operator gives K_j = -dE/dx.
    ! For u_t + u_x = 0 with unit wave speed, K_j is du/dt at node j.
    use precision_select
    use scheme_dispatch, only: interpolate_dispatch
    use cw_params, only: difference_operator
    use difference_operators, only: difference_apply
    implicit none

contains

    subroutine wcns_rhs(nx, halo, u, f)
        ! Compute f(0:nx) = -dE/dx from the interpolated half-node fluxes.
        integer(integer_precision), intent(in) :: nx, halo
        real(precision), intent(in) :: u(-halo:nx+halo)
        real(precision), intent(out) :: f(-halo:nx+halo)

        real(precision) :: flux(-3:nx+2)
        real(precision) :: h
        integer(integer_precision) :: k, j

        ! The physical domain is [-1, 1], so h = 2 / nx.
        h = 2.0_precision/real(nx, precision)

        ! Interpolate every half node needed by the difference stencil.
        ! flux(k) lives at x = (k + 1/2) h.
        flux = 0.0_precision
        do k = -3, nx + 2
            call interpolate_dispatch(u(k-3:k+4), flux(k))
        end do

        ! Apply the selected difference operator at every interior node.
        f = 0.0_precision
        do j = 0, nx
            call difference_apply(difference_operator, h, flux(j-3:j+2), &
                                  u(j-2:j+2), f(j))
        end do
    end subroutine wcns_rhs

end module cw_wcns_rhs
