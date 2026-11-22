module e2d_cases
    ! Case metadata, initial states, halo fills, and source terms for the
    ! seven 2D Euler test cases.
    !
    ! Each case fixes the mesh aspect ratio, the default resolution and end
    ! time, the ratio of specific heats, and the physical-coordinate origin.
    ! The initial-state and boundary routines mirror the reference programs;
    ! only the global names (M -> m, N -> n, MT -> mt, Ubar -> ubar) and the
    ! runtime gas constants differ.
    use precision_select
    use e2d_globals, only: n, m, h, mt, gama, beta, pi, ubar
    implicit none

    private
    public :: e2d_case_meta, e2d_case_allowed, e2d_case_meta_get
    public :: e2d_initial, e2d_boundary, e2d_source

    type :: e2d_case_meta
        character(len=32) :: name = ''
        integer(integer_precision) :: ratio_num = 1
        integer(integer_precision) :: ratio_den = 1
        integer(integer_precision) :: default_nz = 0
        real(precision) :: default_tend = 0.0_precision
        real(precision) :: default_tend2 = -1.0_precision
        real(precision) :: default_gama = 1.4_precision
        real(precision) :: x_origin = 0.0_precision
        real(precision) :: z_origin = 0.0_precision
    end type e2d_case_meta

    ! Fields: name, ratio_num, ratio_den, default_nz, default_tend,
    !         default_tend2, default_gama, x_origin, z_origin.
    type(e2d_case_meta), parameter :: case_table(7) = [ &
        e2d_case_meta('double_mach', 4_integer_precision, 1_integer_precision, &
                      200_integer_precision, 0.2_precision, -1.0_precision, &
                      1.4_precision, 0.0_precision, 0.0_precision), &
        e2d_case_meta('kelvin_helmholtz', 1_integer_precision, 1_integer_precision, &
                      400_integer_precision, 0.5_precision, -1.0_precision, &
                      1.4_precision, 0.0_precision, 0.0_precision), &
        e2d_case_meta('rayleigh_taylor', 1_integer_precision, 4_integer_precision, &
                      400_integer_precision, 1.95_precision, -1.0_precision, &
                      5.0_precision/3.0_precision, 0.0_precision, 0.0_precision), &
        e2d_case_meta('riemann1', 1_integer_precision, 1_integer_precision, &
                      400_integer_precision, 0.8_precision, -1.0_precision, &
                      1.4_precision, 0.0_precision, 0.0_precision), &
        e2d_case_meta('riemann2', 1_integer_precision, 1_integer_precision, &
                      1600_integer_precision, 0.3_precision, -1.0_precision, &
                      1.4_precision, 0.0_precision, 0.0_precision), &
        e2d_case_meta('shock_bubble', 2_integer_precision, 1_integer_precision, &
                      800_integer_precision, 0.4_precision, 0.2_precision, &
                      1.4_precision, -0.1_precision, -0.5_precision), &
        e2d_case_meta('shock_vortex', 2_integer_precision, 1_integer_precision, &
                      180_integer_precision, 0.6_precision, -1.0_precision, &
                      1.4_precision, -0.25_precision, -0.5_precision) ]

contains

    function e2d_case_allowed(name) result(ok)
        ! Return true when the name matches a known case.
        character(len=*), intent(in) :: name
        logical :: ok

        integer :: i

        ok = .false.
        do i = 1, size(case_table)
            if (trim(name) == trim(case_table(i)%name)) then
                ok = .true.
                return
            end if
        end do
    end function e2d_case_allowed

    subroutine e2d_case_meta_get(name, meta, ierr, msg)
        ! Copy the metadata entry for one case, or report an unknown name.
        character(len=*), intent(in) :: name
        type(e2d_case_meta), intent(out) :: meta
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        integer :: i

        meta = e2d_case_meta()
        ierr = 1
        msg = 'unknown case_2d: ' // trim(name)

        do i = 1, size(case_table)
            if (trim(name) == trim(case_table(i)%name)) then
                meta = case_table(i)
                ierr = 0
                msg = ''
                return
            end if
        end do
    end subroutine e2d_case_meta_get

    subroutine e2d_initial(case_name, ierr, msg)
        ! Fill the interior ubar(0:n,0:m,:) with the selected case state.
        ! The halo is left untouched; e2d_boundary fills it.
        character(len=*), intent(in) :: case_name
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        select case (trim(case_name))
        case ('double_mach')
            call case_double_mach_initial()
        case ('kelvin_helmholtz')
            call case_kelvin_helmholtz_initial()
        case ('rayleigh_taylor')
            call case_rayleigh_taylor_initial()
        case ('riemann1')
            call case_riemann1_initial()
        case ('riemann2')
            call case_riemann2_initial()
        case ('shock_bubble')
            call case_shock_bubble_initial()
        case ('shock_vortex')
            call case_shock_vortex_initial()
        case default
            error stop 1
        end select
    end subroutine e2d_initial

    subroutine e2d_boundary(case_name, t)
        ! Fill every halo cell with the case-specific boundary rule.
        character(len=*), intent(in) :: case_name
        real(precision), intent(in) :: t

        select case (trim(case_name))
        case ('double_mach')
            call case_double_mach_boundary(t)
        case ('kelvin_helmholtz')
            call case_kelvin_helmholtz_boundary()
        case ('rayleigh_taylor')
            call case_rayleigh_taylor_boundary()
        case ('riemann1')
            call case_riemann1_boundary()
        case ('riemann2')
            call case_riemann2_boundary()
        case ('shock_bubble')
            call case_shock_bubble_boundary()
        case ('shock_vortex')
            call case_shock_vortex_boundary()
        case default
            error stop 1
        end select
    end subroutine e2d_boundary

    subroutine e2d_source(case_name, rhs)
        ! Add the Rayleigh-Taylor gravity source; the other cases have none.
        character(len=*), intent(in) :: case_name
        real(precision), intent(inout) :: rhs(-mt:n + mt, -mt:m + mt, 4)

        select case (trim(case_name))
        case ('double_mach', 'kelvin_helmholtz', 'riemann1', 'riemann2', &
              'shock_bubble', 'shock_vortex')
            ! No source term.
        case ('rayleigh_taylor')
            call case_rayleigh_taylor_source(rhs)
        case default
            error stop 1
        end select
    end subroutine e2d_source

    ! -----------------------------------------------------------------
    ! double_mach
    ! -----------------------------------------------------------------

    subroutine case_double_mach_initial()
        ! Oblique shock at the lower-left wall; the shocked state fills the
        ! region above the line z = (x - 1/6)*sqrt(3).
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z
        real(precision), parameter :: eps = 1.0e-8_precision
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                x = j*h
                z = k*h
                if (z - (x - 1.0_precision/6.0_precision)*sqrt(3.0_precision) < eps) then
                    ro = 1.4_precision
                    uu = 0.0_precision
                    ww = 0.0_precision
                    ps = 1.0_precision
                else
                    ro = 8.0_precision
                    uu = 7.14471_precision
                    ww = -4.125_precision
                    ps = 116.5_precision
                end if
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_double_mach_initial

    subroutine case_double_mach_boundary(t)
        ! Moving-shock top boundary, inflow or mirror bottom, inflow left,
        ! and zero-order extrapolation right.
        real(precision), intent(in) :: t
        real(precision) :: den, vel_u, vel_w, pre
        real(precision), parameter :: denL = 8.0_precision
        real(precision), parameter :: vel_uL = 7.14471_precision
        real(precision), parameter :: vel_wL = -4.125_precision
        real(precision), parameter :: preL = 116.5_precision
        real(precision), parameter :: eps = 1.0e-8_precision
        real(precision) :: x, xs
        integer(integer_precision) :: j, k, p

        ! Top: the shock front moves as xs(t).
        xs = 1.0_precision/6.0_precision + (1.0_precision + 20.0_precision*t) &
             /sqrt(3.0_precision)
        do j = -mt, n + mt
            x = j*h + 0.5_precision*h
            if (x - xs > eps) then
                den = 1.4_precision
                vel_u = 0.0_precision
                vel_w = 0.0_precision
                pre = 1.0_precision
            else
                den = 8.0_precision
                vel_u = 7.14471_precision
                vel_w = -4.125_precision
                pre = 116.5_precision
            end if
            ubar(j, m + 1:m + mt, 1) = den
            ubar(j, m + 1:m + mt, 2) = den*vel_u
            ubar(j, m + 1:m + mt, 3) = den*vel_w
            ubar(j, m + 1:m + mt, 4) = pre/beta + 0.5_precision*den*(vel_u**2 + vel_w**2)
        end do

        ! Bottom: inflow near the left corner, reflecting elsewhere.
        do j = -mt, n + mt
            x = j*h + 0.5_precision*h
            if (x - 1.0_precision/6.0_precision < eps) then
                ubar(j, -mt:-1, 1) = denL
                ubar(j, -mt:-1, 2) = denL*vel_uL
                ubar(j, -mt:-1, 3) = denL*vel_wL
                ubar(j, -mt:-1, 4) = preL/beta + 0.5_precision*denL*(vel_uL**2 + vel_wL**2)
            else
                do p = 1, mt
                    ubar(j, -p, :) = ubar(j, p - 1, :)
                    ubar(j, -p, 3) = -ubar(j, p - 1, 3)
                end do
            end if
        end do

        ! Left inflow and right zero-order extrapolation.
        do k = -mt, m + mt
            ubar(-mt:-1, k, 1) = denL
            ubar(-mt:-1, k, 2) = denL*vel_uL
            ubar(-mt:-1, k, 3) = denL*vel_wL
            ubar(-mt:-1, k, 4) = preL/beta + 0.5_precision*denL*(vel_uL**2 + vel_wL**2)
            do p = 1, mt
                ubar(n + p, k, :) = ubar(n, k, :)
            end do
        end do
    end subroutine case_double_mach_boundary

    ! -----------------------------------------------------------------
    ! kelvin_helmholtz
    ! -----------------------------------------------------------------

    subroutine case_kelvin_helmholtz_initial()
        ! Two density layers in z with a sinusoidal vertical velocity.
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                x = j*h
                z = k*h
                if (z >= 0.25_precision .and. z <= 0.75_precision) then
                    ro = 2.0_precision
                    uu = -0.5_precision
                    ww = 0.01_precision*sin(2.0_precision*pi*x)
                    ps = 2.5_precision
                else
                    ro = 1.0_precision
                    uu = 0.5_precision
                    ww = 0.01_precision*sin(2.0_precision*pi*x)
                    ps = 2.5_precision
                end if
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_kelvin_helmholtz_initial

    subroutine case_kelvin_helmholtz_boundary()
        ! Periodic halo by four mirror copies.
        integer(integer_precision) :: p

        do p = 1, mt
            ubar(:, -p, :) = ubar(:, p - 1, :)
            ubar(:, m + p, :) = ubar(:, m - p + 1, :)
            ubar(-p, :, :) = ubar(p - 1, :, :)
            ubar(n + p, :, :) = ubar(n - p + 1, :, :)
        end do
    end subroutine case_kelvin_helmholtz_boundary

    ! -----------------------------------------------------------------
    ! rayleigh_taylor
    ! -----------------------------------------------------------------

    subroutine case_rayleigh_taylor_initial()
        ! Heavy fluid below a perturbed interface; hydrostatic pressure.
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z, vc
        integer(integer_precision) :: j, k

        do j = 0, n
            x = j*h
            do k = 0, m/2
                z = k*h
                ro = 2.0_precision
                ps = 2.0_precision*z + 1.0_precision
                vc = sqrt(gama*ps/ro)
                uu = 0.0_precision
                ww = -0.025_precision*vc*cos(8.0_precision*pi*x)
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
            do k = m/2 + 1, m
                z = k*h
                ro = 1.0_precision
                ps = z + 1.5_precision
                vc = sqrt(gama*ps/ro)
                uu = 0.0_precision
                ww = -0.025_precision*vc*cos(8.0_precision*pi*x)
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_rayleigh_taylor_initial

    subroutine case_rayleigh_taylor_boundary()
        ! Fixed top and bottom states; reflecting x walls.
        real(precision) :: ro, ps
        integer(integer_precision) :: k, p

        ro = 1.0_precision
        ps = 2.5_precision
        ubar(0:n, m + 1:m + mt, 1) = ro
        ubar(0:n, m + 1:m + mt, 2) = 0.0_precision
        ubar(0:n, m + 1:m + mt, 3) = 0.0_precision
        ubar(0:n, m + 1:m + mt, 4) = ps/beta

        ro = 2.0_precision
        ps = 1.0_precision
        ubar(0:n, -mt:-1, 1) = ro
        ubar(0:n, -mt:-1, 2) = 0.0_precision
        ubar(0:n, -mt:-1, 3) = 0.0_precision
        ubar(0:n, -mt:-1, 4) = ps/beta

        do p = 1, mt
            do k = -mt, m + mt
                ubar(-p, k, :) = ubar(p, k, :)
                ubar(-p, k, 2) = -ubar(p, k, 2)
                ubar(n + p, k, :) = ubar(n - p, k, :)
                ubar(n + p, k, 2) = -ubar(n - p, k, 2)
            end do
        end do
    end subroutine case_rayleigh_taylor_boundary

    subroutine case_rayleigh_taylor_source(rhs)
        ! Gravity acts on the z momentum and energy.
        real(precision), intent(inout) :: rhs(-mt:n + mt, -mt:m + mt, 4)
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                rhs(j, k, 3) = rhs(j, k, 3) + ubar(j, k, 1)
                rhs(j, k, 4) = rhs(j, k, 4) + ubar(j, k, 3)
            end do
        end do
    end subroutine case_rayleigh_taylor_source

    ! -----------------------------------------------------------------
    ! riemann1
    ! -----------------------------------------------------------------

    subroutine case_riemann1_initial()
        ! Four constant quadrants (2D Riemann problem, configuration 1).
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                x = j*h
                z = k*h
                if (x >= 0.8_precision .and. x <= 1.0_precision .and. &
                    z >= 0.8_precision .and. z <= 1.0_precision) then
                    ro = 1.5_precision
                    uu = 0.0_precision
                    ww = 0.0_precision
                    ps = 1.5_precision
                else if (x >= 0.0_precision .and. x < 0.8_precision .and. &
                         z >= 0.8_precision .and. z <= 1.0_precision) then
                    ro = 0.5323_precision
                    uu = 1.206_precision
                    ww = 0.0_precision
                    ps = 0.3_precision
                else if (x >= 0.0_precision .and. x < 0.8_precision .and. &
                         z >= 0.0_precision .and. z < 0.8_precision) then
                    ro = 0.138_precision
                    uu = 1.206_precision
                    ww = 1.206_precision
                    ps = 0.029_precision
                else
                    ro = 0.5323_precision
                    uu = 0.0_precision
                    ww = 1.206_precision
                    ps = 0.3_precision
                end if
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_riemann1_initial

    subroutine case_riemann1_boundary()
        ! Periodic halo by four mirror copies.
        integer(integer_precision) :: p

        do p = 1, mt
            ubar(:, -p, :) = ubar(:, p - 1, :)
            ubar(:, m + p, :) = ubar(:, m - p + 1, :)
            ubar(-p, :, :) = ubar(p - 1, :, :)
            ubar(n + p, :, :) = ubar(n - p + 1, :, :)
        end do
    end subroutine case_riemann1_boundary

    ! -----------------------------------------------------------------
    ! riemann2
    ! -----------------------------------------------------------------

    subroutine case_riemann2_initial()
        ! Four constant quadrants (2D Riemann problem, configuration 12).
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                x = j*h
                z = k*h
                if (x >= 0.5_precision .and. x <= 1.0_precision .and. &
                    z >= 0.5_precision .and. z <= 1.0_precision) then
                    ro = 1.0_precision
                    uu = 0.75_precision
                    ww = -0.5_precision
                    ps = 1.0_precision
                else if (x >= 0.0_precision .and. x < 0.5_precision .and. &
                         z >= 0.5_precision .and. z <= 1.0_precision) then
                    ro = 2.0_precision
                    uu = 0.75_precision
                    ww = 0.5_precision
                    ps = 1.0_precision
                else if (x >= 0.0_precision .and. x < 0.5_precision .and. &
                         z >= 0.0_precision .and. z < 0.5_precision) then
                    ro = 1.0_precision
                    uu = -0.75_precision
                    ww = 0.5_precision
                    ps = 1.0_precision
                else
                    ro = 3.0_precision
                    uu = -0.75_precision
                    ww = -0.5_precision
                    ps = 1.0_precision
                end if
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_riemann2_initial

    subroutine case_riemann2_boundary()
        ! Periodic halo by four mirror copies.
        integer(integer_precision) :: p

        do p = 1, mt
            ubar(:, -p, :) = ubar(:, p - 1, :)
            ubar(:, m + p, :) = ubar(:, m - p + 1, :)
            ubar(-p, :, :) = ubar(p - 1, :, :)
            ubar(n + p, :, :) = ubar(n - p + 1, :, :)
        end do
    end subroutine case_riemann2_boundary

    ! -----------------------------------------------------------------
    ! shock_bubble
    ! -----------------------------------------------------------------

    subroutine case_shock_bubble_initial()
        ! A strong shock in a tube with a light cylindrical bubble.
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n
                x = j*h - 0.1_precision
                z = k*h - 0.5_precision
                if (x >= -0.1_precision .and. x <= 0.0_precision) then
                    ro = 3.81_precision
                    uu = 2.85_precision
                    ww = 0.0_precision
                    ps = 10.0_precision
                else if (sqrt((x - 0.3_precision)**2 + z*z) <= 0.2_precision) then
                    ro = 0.1_precision
                    uu = 0.0_precision
                    ww = 0.0_precision
                    ps = 1.0_precision
                else
                    ro = 1.0_precision
                    uu = 0.0_precision
                    ww = 0.0_precision
                    ps = 1.0_precision
                end if
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_shock_bubble_initial

    subroutine case_shock_bubble_boundary()
        ! Reflecting top and bottom, inflow left, zero-order extrapolation
        ! right.
        real(precision) :: den, vel_u, vel_w, pre
        integer(integer_precision) :: j, k, p

        do p = 1, mt
            do j = -mt, n + mt
                ubar(j, -p, :) = ubar(j, p, :)
                ubar(j, -p, 3) = -ubar(j, p, 3)
                ubar(j, m + p, :) = ubar(j, m - p, :)
                ubar(j, m + p, 3) = -ubar(j, m - p, 3)
            end do
        end do

        den = 3.81_precision
        vel_u = 2.85_precision
        vel_w = 0.0_precision
        pre = 10.0_precision
        do k = -mt, m + mt
            ubar(-mt:-1, k, 1) = den
            ubar(-mt:-1, k, 2) = den*vel_u
            ubar(-mt:-1, k, 3) = den*vel_w
            ubar(-mt:-1, k, 4) = pre/beta + 0.5_precision*den*(vel_u**2 + vel_w**2)
        end do

        do k = -mt, m + mt
            do p = 1, mt
                ubar(n + p, k, :) = ubar(n, k, :)
            end do
        end do
    end subroutine case_shock_bubble_boundary

    ! -----------------------------------------------------------------
    ! shock_vortex
    ! -----------------------------------------------------------------

    subroutine case_shock_vortex_initial()
        ! Stationary vortex left of a moving shock.
        real(precision) :: ro, uu, ww, ps
        real(precision) :: x, z, tao2, tt
        real(precision), parameter :: eps = 0.3_precision
        real(precision), parameter :: alfa = 0.204_precision
        real(precision), parameter :: rc = 0.05_precision
        integer(integer_precision) :: j, k

        do k = 0, m
            do j = 0, n/4
                x = j*h - 0.25_precision
                z = k*h - 0.50_precision
                tao2 = (x**2 + z**2)/rc**2
                tt = 1.0_precision - eps**2/(14.0_precision*alfa) &
                     *exp(2.0_precision*alfa*(1.0_precision - tao2))
                ro = tt**2.5_precision
                ps = tt**3.5_precision
                uu = 1.30154_precision + eps/rc*exp(alfa*(1.0_precision - tao2))*z
                ww = 0.0_precision - eps/rc*exp(alfa*(1.0_precision - tao2))*x
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
            do j = n/4 + 1, n
                x = j*h - 0.25_precision
                z = k*h - 0.50_precision
                tao2 = (x**2 + z**2)/rc**2
                tt = 1.06494_precision - eps**2/(14.0_precision*alfa) &
                     *exp(2.0_precision*alfa*(1.0_precision - tao2))
                ro = tt**2.5_precision
                ps = tt**3.5_precision
                uu = 1.1133_precision + eps/rc*exp(alfa*(1.0_precision - tao2))*z
                ww = 0.0_precision - eps/rc*exp(alfa*(1.0_precision - tao2))*x
                ubar(j, k, 1) = ro
                ubar(j, k, 2) = ro*uu
                ubar(j, k, 3) = ro*ww
                ubar(j, k, 4) = ps/beta + 0.5_precision*ro*(uu**2 + ww**2)
            end do
        end do
    end subroutine case_shock_vortex_initial

    subroutine case_shock_vortex_boundary()
        ! Periodic halo by four mirror copies.
        integer(integer_precision) :: p

        do p = 1, mt
            ubar(:, -p, :) = ubar(:, p - 1, :)
            ubar(:, m + p, :) = ubar(:, m - p + 1, :)
            ubar(-p, :, :) = ubar(p - 1, :, :)
            ubar(n + p, :, :) = ubar(n - p + 1, :, :)
        end do
    end subroutine case_shock_vortex_boundary

end module e2d_cases
