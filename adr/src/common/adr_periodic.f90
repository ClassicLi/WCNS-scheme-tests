module adr_periodic
    ! Periodic halo handling for the ADR grid.
    use precision_select
    implicit none

contains

    subroutine boundary_cycle(n_max, u, n_out)
        ! Fill n_out halo points on both sides from the periodic interior.
        integer(integer_precision), intent(in) :: n_max, n_out
        real(precision), intent(inout) :: u(-n_out:n_max+n_out)

        integer(integer_precision) :: i

        do i = 1, n_out
            u(n_max + i) = u(i)
            u(-i) = u(n_max - i)
        end do
    end subroutine boundary_cycle

end module adr_periodic
