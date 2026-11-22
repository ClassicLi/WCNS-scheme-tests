program test_flux
    ! Unit tests for the 1D Euler thermodynamics, characteristic matrices,
    ! Steger-Warming split, and characteristic interpolation.
    use precision_select
    use e1d_globals, only: halo
    use e1d_flux
    use test_util
    implicit none

    real(precision) :: u(3), f(3), fh(3)
    real(precision) :: ml(3, 3), mr(3, 3), ident(3, 3)
    real(precision) :: stencil(-3:4, 3), uh(3)
    real(precision) :: uf(-halo:16 + halo, 3), eflux(-3:16 + 2, 3)
    real(precision) :: states(3, 3)
    integer :: correction, s, j, ncorr

    states(:, 1) = (/1.0_precision, 0.3_precision, 2.5_precision/)
    states(:, 2) = (/0.125_precision, 0.0_precision, 0.25_precision/)
    states(:, 3) = (/2.0_precision, -1.2_precision, 10.0_precision/)

    u = states(:, 1)
    call check_close(e1d_pressure(u), 0.982_precision, test_tol, 'pressure')
    call check_close(e1d_sound_speed(u), sqrt(1.4_precision*0.982_precision), &
                     test_tol, 'sound speed')

    call e1d_physical_flux(u, f)
    call check_close(f(1), 0.3_precision, test_tol, 'flux rho')
    call check_close(f(2), 1.072_precision, test_tol, 'flux momentum')
    call check_close(f(3), 1.0446_precision, test_tol, 'flux energy')

    ! Characteristic matrices: mr * ml = I and U = mr * (ml * U).
    ident = 0.0_precision
    ident(1, 1) = 1.0_precision
    ident(2, 2) = 1.0_precision
    ident(3, 3) = 1.0_precision

    do s = 1, 3
        u = states(:, s)
        call e1d_eigen_matrix(u, ml, mr)
        call check_true(maxval(abs(matmul(mr, ml) - ident)) < 100.0_precision*test_tol, &
                        'ml/mr inverse')
        call check_true(maxval(abs(matmul(mr, matmul(ml, u)) - u)) < &
                        100.0_precision*test_tol, 'characteristic round trip')
    end do

    ! Steger-Warming reduces to the physical flux when both states agree.
    u = states(:, 1)
    call e1d_physical_flux(u, f)
    call e1d_sw(u, u, fh)
    call check_true(maxval(abs(fh - f)) < 100.0_precision*test_tol, 'SW consistency')

    ! A uniform stencil is interpolated back to the same state.
    do j = -3, 4
        stencil(j, :) = u
    end do
    call e1d_interpolate_character(stencil, uh, correction)
    call check_true(correction == 0, 'uniform no correction')
    call check_true(maxval(abs(uh - u)) < 100.0_precision*test_tol, &
                    'uniform interpolation')

    ! The same check through the full half-node flux assembly.
    do j = -halo, 16 + halo
        uf(j, :) = u
    end do
    call e1d_interflux(16, uf, eflux, ncorr)
    call check_true(ncorr == 0, 'uniform interflux no correction')
    do j = -3, 16 + 2
        call check_true(maxval(abs(eflux(j, :) - f)) < 100.0_precision*test_tol, &
                        'uniform interflux')
    end do

    ! The positivity guard returns the centre state when the reconstructed
    ! density or pressure falls below 0.1% of the stencil minimum.
    do j = -3, 4
        stencil(j, :) = u
    end do
    uh = (/1.0e-6_precision, 0.0_precision, 2.5_precision/)
    call e1d_apply_positivity(stencil, uh, correction)
    call check_true(correction == 1, 'positivity rho triggers')
    call check_true(maxval(abs(uh - u)) < 100.0_precision*test_tol, &
                    'positivity rho fallback')

    uh = (/1.0_precision, 0.0_precision, 1.0e-6_precision/)
    call e1d_apply_positivity(stencil, uh, correction)
    call check_true(correction == 1, 'positivity p triggers')
    call check_true(maxval(abs(uh - u)) < 100.0_precision*test_tol, &
                    'positivity p fallback')

    ! A tiny but physical uniform state must not trigger the fallback.
    u = (/1.0e-12_precision, 0.0_precision, 1.0e-12_precision/)
    do j = -3, 4
        stencil(j, :) = u
    end do
    call e1d_interpolate_character(stencil, uh, correction)
    call check_true(correction == 0, 'tiny uniform no correction')

    call finish()

end program test_flux
