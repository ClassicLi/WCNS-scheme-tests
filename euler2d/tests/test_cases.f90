program test_cases
    ! Unit tests for the seven 2D Euler cases: initial states, halo fills,
    ! and the Rayleigh-Taylor gravity source.
    use precision_select
    use e2d_globals, only: mt, n, m, beta, ubar, e2d_allocate, &
                           e2d_deallocate, e2d_set_gama
    use e2d_cases, only: e2d_case_meta, e2d_case_meta_get, e2d_initial, &
                         e2d_boundary, e2d_source
    use test_util
    implicit none

    type(e2d_case_meta) :: meta
    integer :: ierr
    character(len=256) :: msg
    character(len=32), parameter :: other_cases(6) = [character(len=32) :: &
        'double_mach', 'kelvin_helmholtz', 'riemann1', 'riemann2', &
        'shock_bubble', 'shock_vortex']
    real(precision) :: rho, uu, ww, ps, err
    real(precision), allocatable :: rhs(:, :, :), rhs_ref(:, :, :)
    integer(integer_precision) :: p, nx, nz
    integer :: c

    nz = 16_integer_precision

    ! ---------------------------------------------------------------
    ! Initial states: one known interior point per case.
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('double_mach', meta, ierr, msg)
    call check_true(ierr == 0, 'double_mach metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('double_mach', ierr, msg)
    call check_true(ierr == 0, 'double_mach initial status')
    call primitive_at(0_integer_precision, 0_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 8.0_precision, 0.0_precision, 'double_mach rho at (0,0)')
    call primitive_at(2*m, m, rho, uu, ww, ps)
    call check_close(rho, 1.4_precision, 0.0_precision, 'double_mach rho at (2m,m)')
    call e2d_deallocate()

    call e2d_case_meta_get('kelvin_helmholtz', meta, ierr, msg)
    call check_true(ierr == 0, 'kelvin_helmholtz metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('kelvin_helmholtz', ierr, msg)
    call check_true(ierr == 0, 'kelvin_helmholtz initial status')
    call primitive_at(0_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 0.0_precision, 'kelvin_helmholtz rho below band')
    call primitive_at(0_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 2.0_precision, 0.0_precision, 'kelvin_helmholtz rho in band')
    call primitive_at(0_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(ww, 0.01_precision*sin(0.0_precision), 1.0e-14_precision, &
                     'kelvin_helmholtz w at x=0')
    call e2d_deallocate()

    call e2d_case_meta_get('rayleigh_taylor', meta, ierr, msg)
    call check_true(ierr == 0, 'rayleigh_taylor metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('rayleigh_taylor', ierr, msg)
    call check_true(ierr == 0, 'rayleigh_taylor initial status')
    call primitive_at(0_integer_precision, 4_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 2.0_precision, 0.0_precision, 'rayleigh_taylor rho lower half')
    call check_close(ps, 1.5_precision, 1.0e-14_precision, 'rayleigh_taylor p lower half')
    call primitive_at(0_integer_precision, 12_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 0.0_precision, 'rayleigh_taylor rho upper half')
    call check_close(ps, 2.25_precision, 1.0e-14_precision, 'rayleigh_taylor p upper half')
    call e2d_deallocate()

    call e2d_case_meta_get('riemann1', meta, ierr, msg)
    call check_true(ierr == 0, 'riemann1 metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('riemann1', ierr, msg)
    call check_true(ierr == 0, 'riemann1 initial status')
    call primitive_at(2_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 0.138_precision, 1.0e-14_precision, 'riemann1 rho (2,2)')
    call check_close(uu, 1.206_precision, 1.0e-14_precision, 'riemann1 u (2,2)')
    call check_close(ww, 1.206_precision, 1.0e-14_precision, 'riemann1 w (2,2)')
    call check_close(ps, 0.029_precision, 1.0e-14_precision, 'riemann1 p (2,2)')
    call primitive_at(14_integer_precision, 14_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.5_precision, 1.0e-14_precision, 'riemann1 rho (14,14)')
    call check_close(uu, 0.0_precision, 1.0e-14_precision, 'riemann1 u (14,14)')
    call check_close(ww, 0.0_precision, 1.0e-14_precision, 'riemann1 w (14,14)')
    call check_close(ps, 1.5_precision, 1.0e-14_precision, 'riemann1 p (14,14)')
    call primitive_at(2_integer_precision, 14_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 0.5323_precision, 1.0e-14_precision, 'riemann1 rho (2,14)')
    call check_close(uu, 1.206_precision, 1.0e-14_precision, 'riemann1 u (2,14)')
    call check_close(ww, 0.0_precision, 1.0e-14_precision, 'riemann1 w (2,14)')
    call check_close(ps, 0.3_precision, 1.0e-14_precision, 'riemann1 p (2,14)')
    call primitive_at(14_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 0.5323_precision, 1.0e-14_precision, 'riemann1 rho (14,2)')
    call check_close(uu, 0.0_precision, 1.0e-14_precision, 'riemann1 u (14,2)')
    call check_close(ww, 1.206_precision, 1.0e-14_precision, 'riemann1 w (14,2)')
    call check_close(ps, 0.3_precision, 1.0e-14_precision, 'riemann1 p (14,2)')
    call e2d_deallocate()

    call e2d_case_meta_get('riemann2', meta, ierr, msg)
    call check_true(ierr == 0, 'riemann2 metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('riemann2', ierr, msg)
    call check_true(ierr == 0, 'riemann2 initial status')
    call primitive_at(2_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 1.0e-14_precision, 'riemann2 rho (2,2)')
    call check_close(uu, -0.75_precision, 1.0e-14_precision, 'riemann2 u (2,2)')
    call check_close(ww, 0.5_precision, 1.0e-14_precision, 'riemann2 w (2,2)')
    call primitive_at(14_integer_precision, 2_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 3.0_precision, 1.0e-14_precision, 'riemann2 rho (14,2)')
    call check_close(uu, -0.75_precision, 1.0e-14_precision, 'riemann2 u (14,2)')
    call check_close(ww, -0.5_precision, 1.0e-14_precision, 'riemann2 w (14,2)')
    call primitive_at(2_integer_precision, 14_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 2.0_precision, 1.0e-14_precision, 'riemann2 rho (2,14)')
    call check_close(uu, 0.75_precision, 1.0e-14_precision, 'riemann2 u (2,14)')
    call check_close(ww, 0.5_precision, 1.0e-14_precision, 'riemann2 w (2,14)')
    call primitive_at(14_integer_precision, 14_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 1.0e-14_precision, 'riemann2 rho (14,14)')
    call check_close(uu, 0.75_precision, 1.0e-14_precision, 'riemann2 u (14,14)')
    call check_close(ww, -0.5_precision, 1.0e-14_precision, 'riemann2 w (14,14)')
    call e2d_deallocate()

    call e2d_case_meta_get('shock_bubble', meta, ierr, msg)
    call check_true(ierr == 0, 'shock_bubble metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('shock_bubble', ierr, msg)
    call check_true(ierr == 0, 'shock_bubble initial status')
    call primitive_at(0_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 3.81_precision, 1.0e-14_precision, 'shock_bubble rho inflow')
    call check_close(uu, 2.85_precision, 1.0e-14_precision, 'shock_bubble u inflow')
    call check_close(ww, 0.0_precision, 1.0e-14_precision, 'shock_bubble w inflow')
    call check_close(ps, 10.0_precision, 1.0e-14_precision, 'shock_bubble p inflow')
    call primitive_at(8_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 0.1_precision, 1.0e-14_precision, 'shock_bubble rho bubble')
    call check_close(ps, 1.0_precision, 1.0e-14_precision, 'shock_bubble p bubble')
    call primitive_at(16_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 1.0e-14_precision, 'shock_bubble rho ambient')
    call check_close(ps, 1.0_precision, 1.0e-14_precision, 'shock_bubble p ambient')
    call e2d_deallocate()

    call e2d_case_meta_get('shock_vortex', meta, ierr, msg)
    call check_true(ierr == 0, 'shock_vortex metadata')
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('shock_vortex', ierr, msg)
    call check_true(ierr == 0, 'shock_vortex initial status')
    call primitive_at(0_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 1.0e-4_precision, 'shock_vortex rho left')
    call check_close(ps, 1.0_precision, 1.0e-4_precision, 'shock_vortex p left')
    call check_close(uu, 1.30154_precision, 1.0e-4_precision, 'shock_vortex u left')
    ! At x = -0.25, z = 0 the swirl is nonzero:
    ! w = -(0.3/0.05)*exp(0.204*(1 - 25))*(-0.25) = 1.5*exp(-4.896).
    call check_close(ww, 1.5_precision*exp(-4.896_precision), 1.0e-14_precision, &
                     'shock_vortex w left')
    call primitive_at(16_integer_precision, 8_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.06494_precision**2.5_precision, 1.0e-4_precision, &
                     'shock_vortex rho right')
    call check_close(ps, 1.06494_precision**3.5_precision, 1.0e-4_precision, &
                     'shock_vortex p right')
    call check_close(uu, 1.1133_precision, 1.0e-4_precision, 'shock_vortex u right')
    call check_close(ww, 0.0_precision, 1.0e-4_precision, 'shock_vortex w right')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! Mirror boundaries: KH, R1, R2, SV (all four halo bands).
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('kelvin_helmholtz', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call check_mirror('kelvin_helmholtz')
    call e2d_deallocate()

    call e2d_case_meta_get('riemann1', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call check_mirror('riemann1')
    call e2d_deallocate()

    call e2d_case_meta_get('riemann2', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call check_mirror('riemann2')
    call e2d_deallocate()

    call e2d_case_meta_get('shock_vortex', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call check_mirror('shock_vortex')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! shock_bubble: reflecting z walls, inflow left, extrapolated right.
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('shock_bubble', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call fill_pattern()
    call e2d_boundary('shock_bubble', 0.0_precision)

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(0:n, -p, 3) + ubar(0:n, p, 3))))
    end do
    call check_true(err == 0.0_precision, 'shock_bubble bottom z-momentum')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(0:n, m + p, 3) + ubar(0:n, m - p, 3))))
    end do
    call check_true(err == 0.0_precision, 'shock_bubble top z-momentum')

    err = maxval(abs(ubar(-mt:-1, -mt:m + mt, 1) - 3.81_precision))
    err = max(err, maxval(abs(ubar(-mt:-1, -mt:m + mt, 2) - 3.81_precision*2.85_precision)))
    err = max(err, maxval(abs(ubar(-mt:-1, -mt:m + mt, 3))))
    err = max(err, maxval(abs(ubar(-mt:-1, -mt:m + mt, 4) - &
                           (10.0_precision/beta + 0.5_precision*3.81_precision* &
                            (2.85_precision**2 + 0.0_precision)))))
    call check_true(err == 0.0_precision, 'shock_bubble left inflow')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(n + p, -mt:m + mt, :) - &
                                   ubar(n, -mt:m + mt, :))))
    end do
    call check_true(err == 0.0_precision, 'shock_bubble right extrapolation')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! double_mach: moving shock top and bottom inflow / mirror / extrapolation.
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('double_mach', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call fill_pattern()
    call e2d_boundary('double_mach', 0.0_precision)

    call check_close(ubar(0_integer_precision, m + 1_integer_precision, 1), &
                     8.0_precision, 0.0_precision, 'double_mach top shock rho at x=0')
    call check_close(ubar(n, m + 1_integer_precision, 1), &
                     1.4_precision, 0.0_precision, 'double_mach top post-shock rho')

    err = maxval(abs(ubar(0_integer_precision, -mt:-1, 1) - 8.0_precision))
    err = max(err, maxval(abs(ubar(0_integer_precision, -mt:-1, 2) - &
                               8.0_precision*7.14471_precision)))
    err = max(err, maxval(abs(ubar(0_integer_precision, -mt:-1, 3) - &
                               8.0_precision*(-4.125_precision))))
    err = max(err, maxval(abs(ubar(0_integer_precision, -mt:-1, 4) - &
                               (116.5_precision/beta + 0.5_precision*8.0_precision* &
                                (7.14471_precision**2 + 4.125_precision**2)))))
    call check_true(err == 0.0_precision, 'double_mach bottom inflow')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, abs(ubar(16_integer_precision, -p, 3) + &
                           ubar(16_integer_precision, p - 1, 3)))
    end do
    call check_true(err == 0.0_precision, 'double_mach bottom mirror')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(n + p, -mt:m + mt, :) - &
                                   ubar(n, -mt:m + mt, :))))
    end do
    call check_true(err == 0.0_precision, 'double_mach right extrapolation')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! rayleigh_taylor: fixed top and bottom, reflecting x walls.
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('rayleigh_taylor', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call fill_pattern()
    call e2d_boundary('rayleigh_taylor', 0.0_precision)

    call primitive_at(0_integer_precision, m + 1_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 1.0_precision, 0.0_precision, 'rayleigh_taylor top rho')
    call check_close(ps, 2.5_precision, 1.0e-14_precision, 'rayleigh_taylor top p')
    call primitive_at(0_integer_precision, -1_integer_precision, rho, uu, ww, ps)
    call check_close(rho, 2.0_precision, 0.0_precision, 'rayleigh_taylor bottom rho')
    call check_close(ps, 1.0_precision, 1.0e-14_precision, 'rayleigh_taylor bottom p')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(-p, -mt:m + mt, 2) + &
                                   ubar(p, -mt:m + mt, 2))))
    end do
    call check_true(err == 0.0_precision, 'rayleigh_taylor left reflection')

    err = 0.0_precision
    do p = 1, mt
        err = max(err, maxval(abs(ubar(n + p, -mt:m + mt, 2) + &
                                   ubar(n - p, -mt:m + mt, 2))))
    end do
    call check_true(err == 0.0_precision, 'rayleigh_taylor right reflection')
    call e2d_deallocate()

    ! ---------------------------------------------------------------
    ! Source terms: only rayleigh_taylor adds gravity.
    ! ---------------------------------------------------------------
    call e2d_case_meta_get('rayleigh_taylor', meta, ierr, msg)
    call e2d_set_gama(meta%default_gama)
    nx = nz*meta%ratio_num/meta%ratio_den
    call e2d_allocate(nx, nz)
    call e2d_initial('rayleigh_taylor', ierr, msg)
    allocate (rhs(-mt:n + mt, -mt:m + mt, 4))
    allocate (rhs_ref(-mt:n + mt, -mt:m + mt, 4))
    rhs = 0.0_precision
    call e2d_source('rayleigh_taylor', rhs)
    call check_close(rhs(0_integer_precision, 0_integer_precision, 3), &
                     ubar(0_integer_precision, 0_integer_precision, 1), &
                     0.0_precision, 'rayleigh_taylor source rho')
    call check_close(rhs(0_integer_precision, 0_integer_precision, 4), &
                     ubar(0_integer_precision, 0_integer_precision, 3), &
                     0.0_precision, 'rayleigh_taylor source rho*w')
    call check_close(rhs(0_integer_precision, 0_integer_precision, 1), &
                     0.0_precision, 0.0_precision, 'rayleigh_taylor source mass')
    call check_close(rhs(0_integer_precision, 0_integer_precision, 2), &
                     0.0_precision, 0.0_precision, 'rayleigh_taylor source x-momentum')

    rhs = -3.5_precision
    do c = 1, size(other_cases)
        rhs_ref = rhs
        call e2d_source(trim(other_cases(c)), rhs)
        call check_true(maxval(abs(rhs - rhs_ref)) == 0.0_precision, &
                        'source no-op for '//trim(other_cases(c)))
    end do
    deallocate (rhs)
    deallocate (rhs_ref)
    call e2d_deallocate()

    call finish()

contains

    subroutine primitive_at(jj, kk, ro, vx, vz, pr)
        ! Recover primitive variables from the conservative state at (jj,kk).
        integer(integer_precision), intent(in) :: jj, kk
        real(precision), intent(out) :: ro, vx, vz, pr

        ro = ubar(jj, kk, 1)
        vx = ubar(jj, kk, 2)/ro
        vz = ubar(jj, kk, 3)/ro
        pr = beta*(ubar(jj, kk, 4) - 0.5_precision*ro*(vx*vx + vz*vz))
    end subroutine primitive_at

    subroutine fill_pattern()
        ! Give every cell, including halos, a distinct value per index.
        integer(integer_precision) :: jj, kk, ll

        do ll = 1, 4
            do kk = -mt, m + mt
                do jj = -mt, n + mt
                    ubar(jj, kk, ll) = real(ll, precision)*1000.0_precision + &
                                       real(jj, precision)*100.0_precision + &
                                       real(kk, precision)
                end do
            end do
        end do
    end subroutine fill_pattern

    subroutine check_mirror(name)
        ! Verify the four mirror halo relations used by periodic cases.
        character(len=*), intent(in) :: name
        integer(integer_precision) :: pp
        real(precision) :: e1, e2, e3, e4

        call fill_pattern()
        call e2d_boundary(name, 0.0_precision)

        e1 = 0.0_precision
        e2 = 0.0_precision
        e3 = 0.0_precision
        e4 = 0.0_precision
        do pp = 1, mt
            e1 = max(e1, maxval(abs(ubar(:, -pp, :) - ubar(:, pp - 1, :))))
            e2 = max(e2, maxval(abs(ubar(:, m + pp, :) - ubar(:, m - pp + 1, :))))
            e3 = max(e3, maxval(abs(ubar(-pp, :, :) - ubar(pp - 1, :, :))))
            e4 = max(e4, maxval(abs(ubar(n + pp, :, :) - ubar(n - pp + 1, :, :))))
        end do

        call check_true(e1 == 0.0_precision, trim(name)//' bottom mirror')
        call check_true(e2 == 0.0_precision, trim(name)//' top mirror')
        call check_true(e3 == 0.0_precision, trim(name)//' left mirror')
        call check_true(e4 == 0.0_precision, trim(name)//' right mirror')
    end subroutine check_mirror

end program test_cases
