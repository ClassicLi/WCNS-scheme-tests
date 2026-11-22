module e1d_initial
    ! Case definitions and initial states for the 1D Euler test cases.
    !
    ! Three classical Riemann problems are supported:
    !   sod        [0, 1],  t = 0.15
    !   lax        [0, 1],  t = 0.15
    !   shu_osher  [0, 10], t = 1.8
    use precision_select
    use e1d_globals, only: halo, beta
    implicit none

contains

    function e1d_case_allowed(name) result(ok)
        ! Return true for a supported case name.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('sod', 'lax', 'shu_osher')
            ok = .true.
        case default
            ok = .false.
        end select
    end function e1d_case_allowed

    function e1d_case_length(name) result(length)
        ! Physical domain length of a case.
        character(len=*), intent(in) :: name
        real(precision) :: length

        select case (trim(name))
        case ('sod', 'lax')
            length = 1.0_precision
        case ('shu_osher')
            length = 10.0_precision
        case default
            write (*, '(a)') 'unknown case: ' // trim(name)
            error stop 1
        end select
    end function e1d_case_length

    function e1d_case_tend(name) result(tend)
        ! Default final time of a case.
        character(len=*), intent(in) :: name
        real(precision) :: tend

        select case (trim(name))
        case ('sod', 'lax')
            tend = 0.15_precision
        case ('shu_osher')
            tend = 1.8_precision
        case default
            write (*, '(a)') 'unknown case: ' // trim(name)
            error stop 1
        end select
    end function e1d_case_tend

    subroutine e1d_case_states(name, rho_l, u_l, p_l, rho_r, u_r, p_r)
        ! Left and right states of a case.
        !
        ! For shu_osher the right state is the base state of the smooth
        ! segment; the sinusoidal density is applied in e1d_initial_state.
        character(len=*), intent(in) :: name
        real(precision), intent(out) :: rho_l, u_l, p_l
        real(precision), intent(out) :: rho_r, u_r, p_r

        select case (trim(name))
        case ('sod')
            rho_l = 1.0_precision
            u_l = 0.0_precision
            p_l = 1.0_precision
            rho_r = 0.125_precision
            u_r = 0.0_precision
            p_r = 0.1_precision
        case ('lax')
            rho_l = 0.445_precision
            u_l = 0.698_precision
            p_l = 3.528_precision
            rho_r = 0.5_precision
            u_r = 0.0_precision
            p_r = 0.571_precision
        case ('shu_osher')
            rho_l = 3.857143_precision
            u_l = 2.629369_precision
            p_l = 10.3333_precision
            rho_r = 1.0_precision
            u_r = 0.0_precision
            p_r = 1.0_precision
        case default
            write (*, '(a)') 'unknown case: ' // trim(name)
            error stop 1
        end select
    end subroutine e1d_case_states

    subroutine e1d_initial_state(n, name, u)
        ! Fill the cell-centred conservative state, including the halo.
        !
        ! The halo uses the same piecewise definition as the interior; the
        ! solver later replaces it with zero-order extrapolation.
        integer, intent(in) :: n
        character(len=*), intent(in) :: name
        real(precision), intent(out) :: u(-halo:n + halo, 3)

        real(precision) :: rho_l, u_l, p_l, rho_r, u_r, p_r
        real(precision) :: h, x, rho, vel, pre
        integer :: j

        call e1d_case_states(name, rho_l, u_l, p_l, rho_r, u_r, p_r)

        select case (trim(name))
        case ('sod', 'lax')
            do j = -halo, n + halo
                if (j <= n/2) then
                    rho = rho_l
                    vel = u_l
                    pre = p_l
                else
                    rho = rho_r
                    vel = u_r
                    pre = p_r
                end if
                u(j, 1) = rho
                u(j, 2) = rho*vel
                u(j, 3) = pre/beta + 0.5_precision*rho*vel*vel
            end do
        case ('shu_osher')
            h = 10.0_precision/real(n, precision)
            do j = -halo, n + halo
                if (j < n/10) then
                    rho = rho_l
                    vel = u_l
                    pre = p_l
                else
                    x = real(j, precision)*h
                    rho = 1.0_precision + 0.2_precision*sin(5.0_precision*x)
                    vel = u_r
                    pre = p_r
                end if
                u(j, 1) = rho
                u(j, 2) = rho*vel
                u(j, 3) = pre/beta + 0.5_precision*rho*vel*vel
            end do
        case default
            write (*, '(a)') 'unknown case: ' // trim(name)
            error stop 1
        end select
    end subroutine e1d_initial_state

end module e1d_initial
