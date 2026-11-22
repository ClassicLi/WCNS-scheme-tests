module diagnostics
    use precision_select
    use scheme_diagnostics_mod, only: scheme_diagnostics
    implicit none

    type(scheme_diagnostics), allocatable :: diag_point(:)

contains

    subroutine init_diagnostics(n)
        integer, intent(in) :: n
        integer :: i

        if (allocated(diag_point)) deallocate(diag_point)
        allocate(diag_point(0:n - 1))
        do i = 0, n - 1
            diag_point(i) = scheme_diagnostics()
        end do
    end subroutine init_diagnostics

    subroutine store_diagnostics(i, diag)
        integer, intent(in) :: i
        type(scheme_diagnostics), intent(in) :: diag

        diag_point(i) = diag
    end subroutine store_diagnostics

end module diagnostics
