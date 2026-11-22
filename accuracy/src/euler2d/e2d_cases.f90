module e2d_cases
    ! Initial conditions and exact densities for the 2D accuracy cases.
    !
    ! Both cases are periodic and advected by a uniform velocity field:
    !   vortex: isentropic vortex on [0, 10]^2, centre (5, 5);
    !   sine  : rho = rho0 + amp*sin(x + y), u = w = 1, p = 1.
    use precision_select
    use e2d_globals, only: n, m, h, ubar, gama, beta
    implicit none

    real(precision), parameter :: pi = acos(-1.0_precision)

contains

    subroutine e2d_case_length(case_name, length, ierr, msg)
        ! Return the physical side length of the selected case.
        character(len=*), intent(in) :: case_name
        real(precision), intent(out) :: length
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''
        length = 0.0_precision

        select case (trim(case_name))
        case ('vortex')
            length = 10.0_precision
        case ('sine')
            length = 2.0_precision*pi
        case default
            ierr = 1
            msg = 'unknown case_2d: ' // trim(case_name)
        end select
    end subroutine e2d_case_length

    subroutine e2d_initial(case_name, rho0, amp, ierr, msg)
        ! Fill ubar(0:n, 0:m, 1:4); the halo is filled by e2d_boundary.
        character(len=*), intent(in) :: case_name
        real(precision), intent(in) :: rho0, amp
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        real(precision) :: x, y, rho, uu, ww, ps
        integer :: j, k

        ierr = 0
        msg = ''

        select case (trim(case_name))
        case ('vortex')
            call vortex_initial()
        case ('sine')
            do k = 0, m
                y = real(k, precision)*h
                do j = 0, n
                    x = real(j, precision)*h
                    rho = rho0 + amp*sin(x + y)
                    uu = 1.0_precision
                    ww = 1.0_precision
                    ps = 1.0_precision
                    ubar(j, k, 1) = rho
                    ubar(j, k, 2) = rho*uu
                    ubar(j, k, 3) = rho*ww
                    ubar(j, k, 4) = ps/beta + 0.5_precision*rho*(uu*uu + ww*ww)
                end do
            end do
        case default
            ierr = 1
            msg = 'unknown case_2d: ' // trim(case_name)
        end select
    end subroutine e2d_initial

    subroutine vortex_initial()
        ! Reference isentropic vortex at t = 0.
        real(precision), parameter :: EKC = 0.3_precision
        real(precision), parameter :: RC = 0.5_precision
        real(precision), parameter :: RFA = 0.35_precision
        real(precision) :: xbar, ybar, rb2, rb, taob2, tempxp
        real(precision) :: rho, uu, ww, ps, temp
        integer :: j, k

        do k = 0, m
            ybar = real(k, precision)*h - 5.0_precision
            do j = 0, n
                xbar = real(j, precision)*h - 5.0_precision
                rb2 = xbar*xbar + ybar*ybar
                rb = sqrt(rb2)
                taob2 = (rb/RC)**2
                tempxp = exp(RFA*(1.0_precision - taob2))
                uu = 1.0_precision + EKC/RC*tempxp*ybar
                ww = -EKC/RC*tempxp*xbar
                temp = 1.0_precision - (gama - 1.0_precision)*EKC*EKC*tempxp*tempxp &
                       /(4.0_precision*RFA*gama)
                rho = temp**(1.0_precision/(gama - 1.0_precision))
                ps = rho**gama
                ubar(j, k, 1) = rho
                ubar(j, k, 2) = rho*uu
                ubar(j, k, 3) = rho*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*rho*(uu*uu + ww*ww)
            end do
        end do
    end subroutine vortex_initial

    subroutine e2d_exact_density(case_name, time, rho0, amp, rho, ierr, msg)
        ! Exact density at the requested time on the (n+1) x (m+1) grid.
        character(len=*), intent(in) :: case_name
        real(precision), intent(in) :: time, rho0, amp
        real(precision), intent(out) :: rho(0:n, 0:m)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        real(precision) :: x, y
        integer :: j, k

        ierr = 0
        msg = ''

        select case (trim(case_name))
        case ('vortex')
            do k = 0, m
                do j = 0, n
                    rho(j, k) = vortex_density(real(j, precision)*h - 5.0_precision - time, &
                                               real(k, precision)*h - 5.0_precision)
                end do
            end do
        case ('sine')
            do k = 0, m
                y = real(k, precision)*h
                do j = 0, n
                    x = real(j, precision)*h
                    rho(j, k) = rho0 + amp*sin(x + y - 2.0_precision*time)
                end do
            end do
        case default
            ierr = 1
            msg = 'unknown case_2d: ' // trim(case_name)
        end select
    end subroutine e2d_exact_density

    function vortex_density(xbar, ybar) result(rho)
        ! Isentropic-vortex density at a point relative to the vortex centre.
        real(precision), intent(in) :: xbar, ybar
        real(precision) :: rho

        real(precision), parameter :: EKC = 0.3_precision
        real(precision), parameter :: RC = 0.5_precision
        real(precision), parameter :: RFA = 0.35_precision
        real(precision) :: rb2, rb, taob2, tempxp, temp

        rb2 = xbar*xbar + ybar*ybar
        rb = sqrt(rb2)
        taob2 = (rb/RC)**2
        tempxp = exp(RFA*(1.0_precision - taob2))
        temp = 1.0_precision - (gama - 1.0_precision)*EKC*EKC*tempxp*tempxp &
               /(4.0_precision*RFA*gama)
        rho = temp**(1.0_precision/(gama - 1.0_precision))
    end function vortex_density

end module e2d_cases
