module e2d_globals
    ! Shared mesh, gas constants, and conservative fields for 2D Euler.
    !
    ! The halo covers the widest interpolation stencil (-mt:mt) and the
    ! integer-node range read by the hybrid difference operators.
    use precision_select
    implicit none

    integer(integer_precision), parameter :: mt = 7_integer_precision

    integer(integer_precision) :: n = 0_integer_precision
    integer(integer_precision) :: m = 0_integer_precision
    real(precision) :: h = 0.0_precision
    real(precision) :: gama = 1.4_precision
    real(precision) :: beta = 0.4_precision
    real(precision), parameter :: pi = acos(-1.0_precision)

    ! Cell-centred conservative state.
    real(precision), allocatable :: ubar(:, :, :)
    real(precision), allocatable :: eint(:, :, :)
    real(precision), allocatable :: gint(:, :, :)

    ! Half-node fluxes, wide enough for every difference operator.
    real(precision), allocatable :: eflux(:, :, :)
    real(precision), allocatable :: gflux(:, :, :)

contains

    subroutine e2d_set_gama(g)
        ! Set the ratio of specific heats and the derived pressure coefficient.
        real(precision), intent(in) :: g

        gama = g
        beta = gama - 1.0_precision
    end subroutine e2d_set_gama

    subroutine e2d_allocate(nx, nz)
        ! Allocate every field for an (nx+1) x (nz+1) grid on [0, ratio] x [0, 1].
        integer(integer_precision), intent(in) :: nx, nz

        n = nx
        m = nz
        h = 1.0_precision/real(m, precision)

        allocate (ubar(-mt:n + mt, -mt:m + mt, 4))
        allocate (eint(-mt:n + mt, -mt:m + mt, 4))
        allocate (gint(-mt:n + mt, -mt:m + mt, 4))
        allocate (eflux(-3:n + 2, 0:m, 4))
        allocate (gflux(0:n, -3:m + 2, 4))

        ubar = 0.0_precision
        eint = 0.0_precision
        gint = 0.0_precision
        eflux = 0.0_precision
        gflux = 0.0_precision
    end subroutine e2d_allocate

    subroutine e2d_deallocate()
        ! Release every field allocated by e2d_allocate.
        if (allocated(ubar)) deallocate (ubar)
        if (allocated(eint)) deallocate (eint)
        if (allocated(gint)) deallocate (gint)
        if (allocated(eflux)) deallocate (eflux)
        if (allocated(gflux)) deallocate (gflux)
    end subroutine e2d_deallocate

end module e2d_globals
