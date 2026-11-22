module c1d_initial_mod
    ! Initial profiles for the 1D convection accuracy test.
    !
    ! Every selectable profile is a smooth function of x on the periodic
    ! domain [-1, 1].  The exact solution is the same profile shifted by t.
    use precision_select
    implicit none

    real(precision), parameter :: pi = acos(-1.0_precision)
    integer(integer_precision), parameter :: halo = 7_integer_precision

contains

    function c1d_profile(x, initial_name) result(value)
        ! Value of the selected initial profile at x.
        real(precision), intent(in) :: x
        character(len=*), intent(in) :: initial_name
        real(precision) :: value

        select case (trim(initial_name))
        case ('sin2')
            value = sin(pi*x)**2
        case ('sin')
            value = sin(pi*x)
        case ('sin3')
            value = sin(pi*x)**3
        case ('sin_phase')
            value = sin(pi*x - sin(pi*x)/pi)
        case ('gauss')
            value = exp(-20.0_precision*x*x)
        case ('sqrt')
            value = sqrt(1.0_precision + 20.0_precision*x*x)
        case ('rational')
            value = 1.0_precision/(1.0_precision + 20.0_precision*x*x)
        case default
            write (*, '(a)') 'unknown initial_1d: ' // trim(initial_name)
            error stop 1
        end select
    end function c1d_profile

    subroutine c1d_initial(n, initial_name, u)
        ! Fill the periodic grid and its halo with the selected profile.
        integer, intent(in) :: n
        character(len=*), intent(in) :: initial_name
        real(precision), intent(out) :: u(-halo:n + halo)

        real(precision) :: h
        integer :: k

        h = 2.0_precision/real(n, precision)
        do k = -halo, n + halo
            u(k) = c1d_profile(-1.0_precision + real(k, precision)*h, initial_name)
        end do
    end subroutine c1d_initial

end module c1d_initial_mod
