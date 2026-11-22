program test_flux
    ! Unit tests for the 2D Euler characteristic matrices, modified
    ! Steger-Warming split, integer-node fluxes, and the characteristic
    ! half-node interpolation pipeline.
    use precision_select
    use e2d_globals, only: e2d_allocate, e2d_deallocate, ubar, eint, gint, &
                           eflux, gflux, beta, n, m
    use e2d_flux, only: u_to_prim_e_g, interflux, sw_x, sw_z, &
                        eigen_matrix_x, eigen_matrix_z
    use test_util
    implicit none

    real(precision) :: u(4), flux(4)
    real(precision) :: ml(4, 4), mr(4, 4), ident(4, 4)
    real(precision) :: ro, uu, ww, ps, energy, err_max
    integer :: j, k, l, ncorr

    call e2d_allocate(4_integer_precision, 4_integer_precision)

    u = (/1.0_precision, 0.3_precision, -0.2_precision, 2.5_precision/)

    ! ---------------------------------------------------------------
    ! Characteristic matrices: ml * mr = I.
    ! ---------------------------------------------------------------
    ident = 0.0_precision
    do l = 1, 4
        ident(l, l) = 1.0_precision
    end do

    call eigen_matrix_x(u, ml, mr)
    call check_true(maxval(abs(matmul(ml, mr) - ident)) < 1.0e-12_precision, &
                    'eigen_matrix_x inverse')

    call eigen_matrix_z(u, ml, mr)
    call check_true(maxval(abs(matmul(ml, mr) - ident)) < 1.0e-12_precision, &
                    'eigen_matrix_z inverse')

    ! ---------------------------------------------------------------
    ! Steger-Warming split with equal states reduces to Eint/Gint.
    ! ---------------------------------------------------------------
    do l = 1, 4
        ubar(:, :, l) = u(l)
    end do
    call u_to_prim_e_g()

    call sw_x(u, u, flux)
    call check_true(maxval(abs(flux - eint(0, 0, :))) < 1.0e-12_precision, &
                    'sw_x physical flux')

    call sw_z(u, u, flux)
    call check_true(maxval(abs(flux - gint(0, 0, :))) < 1.0e-12_precision, &
                    'sw_z physical flux')

    ! ---------------------------------------------------------------
    ! Primitive/energy relations at the integer node.
    ! ---------------------------------------------------------------
    ro = 1.0_precision
    uu = 0.3_precision
    ww = -0.2_precision
    ps = 1.0_precision
    energy = ps/beta + 0.5_precision*(uu*uu + ww*ww)*ro

    ! Fill the whole grid with this valid state: u_to_prim_e_g reads every
    ! halo cell, so leaving zero-density cells would divide by zero.
    ubar(:, :, 1) = ro
    ubar(:, :, 2) = ro*uu
    ubar(:, :, 3) = ro*ww
    ubar(:, :, 4) = energy
    call u_to_prim_e_g()

    call check_close(eint(0, 0, 2), ro*uu*uu + ps, 1.0e-14_precision, &
                     'eint x-momentum flux')
    call check_close(eint(0, 0, 4), (energy + ps)*uu, 1.0e-14_precision, &
                     'eint energy flux')

    ! ---------------------------------------------------------------
    ! A uniform state survives the full half-node interpolation pipeline.
    ! ---------------------------------------------------------------
    do l = 1, 4
        ubar(:, :, l) = u(l)
    end do
    call u_to_prim_e_g()
    call interflux(ncorr)
    call check_true(ncorr == 0, 'uniform interflux no correction')

    err_max = 0.0_precision
    do k = 0, m
        do j = -3, n + 2
            err_max = max(err_max, maxval(abs(eflux(j, k, :) - eint(0, 0, :))))
        end do
    end do
    call check_true(err_max < 1.0e-10_precision, 'uniform eflux')

    err_max = 0.0_precision
    do k = -3, m + 2
        do j = 0, n
            err_max = max(err_max, maxval(abs(gflux(j, k, :) - gint(0, 0, :))))
        end do
    end do
    call check_true(err_max < 1.0e-10_precision, 'uniform gflux')

    call e2d_deallocate()
    call finish()

end program test_flux
