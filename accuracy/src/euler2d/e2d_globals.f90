module e2d_globals
    ! Global fields for the 2D Euler accuracy tests.
    !
    ! The mesh size is a runtime value so one process can refine through all
    ! levels.  Every large field is allocatable; the halo is wide enough for
    ! the interpolation and difference stencils at j = -3 .. n+2.
    use precision_select
    implicit none

    integer(integer_precision), parameter :: mt = 7_integer_precision

    integer :: n = 0, m = 0
    real(precision) :: h = 0.0_precision
    real(precision) :: domain_length = 10.0_precision

    ! Cell-centred conservative state and integer-node fluxes.
    real(precision), allocatable :: ubar(:, :, :)
    real(precision), allocatable :: eint(:, :, :)
    real(precision), allocatable :: gint(:, :, :)

    ! Half-node fluxes, wide enough for every difference operator.
    real(precision), allocatable :: eflux(:, :, :)
    real(precision), allocatable :: gflux(:, :, :)

    real(precision), parameter :: gama = 1.4_precision
    real(precision), parameter :: beta = gama - 1.0_precision
    real(precision), parameter :: eps = 1.0e-40_precision

contains

    subroutine e2d_allocate(nx, length)
        ! Allocate every field for an (nx+1) x (nx+1) grid on [0, length]^2.
        integer, intent(in) :: nx
        real(precision), intent(in) :: length

        if (allocated(ubar)) call e2d_deallocate()

        n = nx
        m = nx
        domain_length = length
        h = length/real(n, precision)

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
