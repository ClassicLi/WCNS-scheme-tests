module e2d_flux
    ! Spatial flux assembly for the 2D Euler accuracy tests.
    !
    ! The pipeline is kept identical to the reference program:
    !   1. E5_X / E5_Z project a 7-point stencil to characteristic variables,
    !      interpolate each characteristic field to the half node, and project
    !      back to a conserved half-node state;
    !   2. InterFlux reconstructs left/right states and applies the modified
    !      Steger-Warming split to obtain Eflux/Gflux;
    !   3. e2d_solver applies the selected difference operator to the half-node
    !      fluxes and the integer-node fluxes Eint/Gint.
    !
    ! Only the interpolation call inside E5_X / E5_Z is replaced by the shared
    ! scheme_dispatch standard part.
    use precision_select
    use e2d_globals
    use scheme_dispatch, only: interpolate_dispatch
    implicit none

contains

    subroutine u_to_prim_e_g()
        ! Integer-node fluxes Eint/Gint from the conservative state.
        integer :: j, k
        real(precision) :: ro, uu, ww, ps

        do k = -mt, m + mt
            do j = -mt, n + mt
                ro = ubar(j, k, 1)
                uu = ubar(j, k, 2)/ro
                ww = ubar(j, k, 3)/ro
                ps = beta*(ubar(j, k, 4) - 0.5_precision*(uu*uu + ww*ww)*ro)

                eint(j, k, 1) = ubar(j, k, 2)
                eint(j, k, 2) = ubar(j, k, 2)*uu + ps
                eint(j, k, 3) = ubar(j, k, 2)*ww
                eint(j, k, 4) = (ubar(j, k, 4) + ps)*uu

                gint(j, k, 1) = ubar(j, k, 3)
                gint(j, k, 2) = ubar(j, k, 3)*uu
                gint(j, k, 3) = ubar(j, k, 3)*ww + ps
                gint(j, k, 4) = (ubar(j, k, 4) + ps)*ww
            end do
        end do
    end subroutine u_to_prim_e_g

    subroutine interflux()
        ! Half-node fluxes from left and right characteristic reconstructions.
        integer :: j, k
        real(precision) :: ul(4), ur(4)
        real(precision) :: qplate_l(-3:4, 4), qplate_r(-3:4, 4)
        integer :: correction

        do k = 0, m
            do j = -3, n + 2
                qplate_l = ubar(j - 3:j + 4, k, :)
                qplate_r = ubar(j + 4:j - 3:-1, k, :)
                call e5_x(qplate_l, ul, correction)
                call e5_x(qplate_r, ur, correction)
                call sw_x(ul, ur, eflux(j, k, :))
            end do
        end do

        do k = -3, m + 2
            do j = 0, n
                qplate_l = ubar(j, k - 3:k + 4, :)
                qplate_r = ubar(j, k + 4:k - 3:-1, :)
                call e5_z(qplate_l, ul, correction)
                call e5_z(qplate_r, ur, correction)
                call sw_z(ul, ur, gflux(j, k, :))
            end do
        end do
    end subroutine interflux

    subroutine e5_x(u, uhalf, correction)
        ! Characteristic interpolation of one x-direction stencil.
        real(precision), intent(in) :: u(-3:4, 4)
        real(precision), intent(out) :: uhalf(4)
        integer, intent(out) :: correction

        real(precision) :: q(-3:4, 4), qhalf(4)
        real(precision) :: ml(4, 4), mr(4, 4)
        real(precision) :: ro, uu, ww, ps
        real(precision) :: prim(-3:4, 4), primhalf(4)
        real(precision) :: denminvalue, preminvalue
        integer :: p, j

        call eigen_matrix_x(u(0, :), ml, mr)

        do p = 1, 4
            q(:, p) = ml(p, 1)*u(:, 1) + ml(p, 2)*u(:, 2) &
                    + ml(p, 3)*u(:, 3) + ml(p, 4)*u(:, 4)
            call interpolate_dispatch(q(:, p), qhalf(p))
        end do

        uhalf(:) = mr(:, 1)*qhalf(1) + mr(:, 2)*qhalf(2) &
                 + mr(:, 3)*qhalf(3) + mr(:, 4)*qhalf(4)

        ! Positivity fallback: if the reconstructed density or pressure is
        ! below 0.1% of the stencil minimum, use the centre state.
        do j = -2, 3
            ro = u(j, 1)
            uu = u(j, 2)/ro
            ww = u(j, 3)/ro
            ps = beta*(u(j, 4) - 0.5_precision*(uu*uu + ww*ww)*ro)
            prim(j, 1) = ro
            prim(j, 2) = uu
            prim(j, 3) = ww
            prim(j, 4) = ps
        end do

        ro = uhalf(1)
        uu = uhalf(2)/ro
        ww = uhalf(3)/ro
        ps = beta*(uhalf(4) - 0.5_precision*(uu*uu + ww*ww)*ro)
        primhalf(1) = ro
        primhalf(2) = uu
        primhalf(3) = ww
        primhalf(4) = ps

        correction = 0
        denminvalue = 0.001_precision*min(prim(-2, 1), prim(-1, 1), prim(0, 1), &
                                          prim(1, 1), prim(2, 1))
        preminvalue = 0.001_precision*min(prim(-2, 4), prim(-1, 4), prim(0, 4), &
                                          prim(1, 4), prim(2, 4))
        if (primhalf(1) <= denminvalue .or. primhalf(4) <= preminvalue) then
            uhalf(:) = u(0, :)
            correction = 1
        end if
    end subroutine e5_x

    subroutine e5_z(u, uhalf, correction)
        ! Characteristic interpolation of one z-direction stencil.
        real(precision), intent(in) :: u(-3:4, 4)
        real(precision), intent(out) :: uhalf(4)
        integer, intent(out) :: correction

        real(precision) :: q(-3:4, 4), qhalf(4)
        real(precision) :: ml(4, 4), mr(4, 4)
        real(precision) :: ro, uu, ww, ps
        real(precision) :: prim(-3:4, 4), primhalf(4)
        real(precision) :: denminvalue, preminvalue
        integer :: p, j

        call eigen_matrix_z(u(0, :), ml, mr)

        do p = 1, 4
            q(:, p) = ml(p, 1)*u(:, 1) + ml(p, 2)*u(:, 2) &
                    + ml(p, 3)*u(:, 3) + ml(p, 4)*u(:, 4)
            call interpolate_dispatch(q(:, p), qhalf(p))
        end do

        uhalf(:) = mr(:, 1)*qhalf(1) + mr(:, 2)*qhalf(2) &
                 + mr(:, 3)*qhalf(3) + mr(:, 4)*qhalf(4)

        do j = -2, 3
            ro = u(j, 1)
            uu = u(j, 2)/ro
            ww = u(j, 3)/ro
            ps = beta*(u(j, 4) - 0.5_precision*(uu*uu + ww*ww)*ro)
            prim(j, 1) = ro
            prim(j, 2) = uu
            prim(j, 3) = ww
            prim(j, 4) = ps
        end do

        ro = uhalf(1)
        uu = uhalf(2)/ro
        ww = uhalf(3)/ro
        ps = beta*(uhalf(4) - 0.5_precision*(uu*uu + ww*ww)*ro)
        primhalf(1) = ro
        primhalf(2) = uu
        primhalf(3) = ww
        primhalf(4) = ps

        correction = 0
        denminvalue = 0.001_precision*min(prim(-2, 1), prim(-1, 1), prim(0, 1), &
                                          prim(1, 1), prim(2, 1))
        preminvalue = 0.001_precision*min(prim(-2, 4), prim(-1, 4), prim(0, 4), &
                                          prim(1, 4), prim(2, 4))
        if (primhalf(1) <= denminvalue .or. primhalf(4) <= preminvalue) then
            uhalf(:) = u(0, :)
            correction = 1
        end if
    end subroutine e5_z

    subroutine sw_x(ul, ur, eflux_half)
        ! Modified Steger-Warming split in the x direction.
        real(precision), intent(in) :: ul(4), ur(4)
        real(precision), intent(out) :: eflux_half(4)

        real(precision), parameter :: eps2 = 0.12_precision**2
        real(precision), parameter :: sig = 0.5_precision/gama
        real(precision) :: ro, uu, ww, ps, alfa, c2, c
        real(precision) :: lam1, lam4, lam5, diss, disp
        real(precision) :: fp(4), fn(4)

        ! Left state, positive eigenvalues.
        ro = ul(1)
        uu = ul(2)/ro
        ww = ul(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(ul(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        lam1 = uu
        lam4 = uu + c
        lam5 = uu - c
        lam1 = (lam1 + sqrt(lam1**2 + eps2))/2.0_precision
        lam4 = (lam4 + sqrt(lam4**2 + eps2))/2.0_precision
        lam5 = (lam5 + sqrt(lam5**2 + eps2))/2.0_precision

        diss = -sig*(2.0_precision*lam1 - lam4 - lam5)
        disp = sig*(lam4 - lam5)*c
        fp(1) = diss + lam1
        fp(2) = uu*diss + uu*lam1 + disp
        fp(3) = ww*diss + ww*lam1
        fp(4) = (alfa + c2/beta)*diss + (alfa + ps/ro/beta)*lam1 + uu*disp

        ! Right state, negative eigenvalues.
        ro = ur(1)
        uu = ur(2)/ro
        ww = ur(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(ur(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        lam1 = uu
        lam4 = uu + c
        lam5 = uu - c
        lam1 = (lam1 - sqrt(lam1**2 + eps2))/2.0_precision
        lam4 = (lam4 - sqrt(lam4**2 + eps2))/2.0_precision
        lam5 = (lam5 - sqrt(lam5**2 + eps2))/2.0_precision

        diss = -sig*(2.0_precision*lam1 - lam4 - lam5)
        disp = sig*(lam4 - lam5)*c
        fn(1) = diss + lam1
        fn(2) = uu*diss + uu*lam1 + disp
        fn(3) = ww*diss + ww*lam1
        fn(4) = (alfa + c2/beta)*diss + (alfa + ps/ro/beta)*lam1 + uu*disp

        eflux_half = ul(1)*fp + ur(1)*fn
    end subroutine sw_x

    subroutine sw_z(ul, ur, gflux_half)
        ! Modified Steger-Warming split in the z direction.
        real(precision), intent(in) :: ul(4), ur(4)
        real(precision), intent(out) :: gflux_half(4)

        real(precision), parameter :: eps2 = 0.12_precision**2
        real(precision), parameter :: sig = 0.5_precision/gama
        real(precision) :: ro, uu, ww, ps, alfa, c2, c
        real(precision) :: lam1, lam4, lam5, diss, disp
        real(precision) :: fp(4), fn(4)

        ! Left state, positive eigenvalues.
        ro = ul(1)
        uu = ul(2)/ro
        ww = ul(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(ul(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        lam1 = ww
        lam4 = ww + c
        lam5 = ww - c
        lam1 = (lam1 + sqrt(lam1**2 + eps2))/2.0_precision
        lam4 = (lam4 + sqrt(lam4**2 + eps2))/2.0_precision
        lam5 = (lam5 + sqrt(lam5**2 + eps2))/2.0_precision

        diss = -sig*(2.0_precision*lam1 - lam4 - lam5)
        disp = sig*(lam4 - lam5)*c
        fp(1) = diss + lam1
        fp(2) = uu*diss + uu*lam1
        fp(3) = ww*diss + ww*lam1 + disp
        fp(4) = (alfa + c2/beta)*diss + (alfa + ps/ro/beta)*lam1 + ww*disp

        ! Right state, negative eigenvalues.
        ro = ur(1)
        uu = ur(2)/ro
        ww = ur(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(ur(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        lam1 = ww
        lam4 = ww + c
        lam5 = ww - c
        lam1 = (lam1 - sqrt(lam1**2 + eps2))/2.0_precision
        lam4 = (lam4 - sqrt(lam4**2 + eps2))/2.0_precision
        lam5 = (lam5 - sqrt(lam5**2 + eps2))/2.0_precision

        diss = -sig*(2.0_precision*lam1 - lam4 - lam5)
        disp = sig*(lam4 - lam5)*c
        fn(1) = diss + lam1
        fn(2) = uu*diss + uu*lam1
        fn(3) = ww*diss + ww*lam1 + disp
        fn(4) = (alfa + c2/beta)*diss + (alfa + ps/ro/beta)*lam1 + ww*disp

        gflux_half = ul(1)*fp + ur(1)*fn
    end subroutine sw_z

    subroutine eigen_matrix_x(u, ml, mr)
        ! Left/right characteristic matrices for the x direction.
        real(precision), intent(in) :: u(4)
        real(precision), intent(out) :: ml(4, 4), mr(4, 4)

        real(precision) :: ro, uu, ww, ps, c, c2, alfa

        ro = u(1)
        uu = u(2)/ro
        ww = u(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(u(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        ml(1, 1) = alfa*beta - c2
        ml(1, 2) = -beta*uu
        ml(1, 3) = -beta*ww
        ml(1, 4) = beta
        ml(2, 1) = -ww
        ml(2, 2) = 0.0_precision
        ml(2, 3) = 1.0_precision
        ml(2, 4) = 0.0_precision
        ml(3, 1) = alfa*beta - c*uu
        ml(3, 2) = c - beta*uu
        ml(3, 3) = -beta*ww
        ml(3, 4) = beta
        ml(4, 1) = alfa*beta + c*uu
        ml(4, 2) = -c - beta*uu
        ml(4, 3) = -beta*ww
        ml(4, 4) = beta

        mr(1, 1) = -2.0_precision
        mr(2, 1) = -2.0_precision*uu
        mr(3, 1) = -2.0_precision*ww
        mr(4, 1) = -2.0_precision*alfa
        mr(1, 2) = 0.0_precision
        mr(2, 2) = 0.0_precision
        mr(3, 2) = 2.0_precision*c2
        mr(4, 2) = ww*2.0_precision*c2
        mr(1, 3) = 1.0_precision
        mr(2, 3) = uu + c
        mr(3, 3) = ww
        mr(4, 3) = alfa + c*uu + c2/beta
        mr(1, 4) = 1.0_precision
        mr(2, 4) = uu - c
        mr(3, 4) = ww
        mr(4, 4) = alfa - c*uu + c2/beta
        mr = mr*0.5_precision/c2
    end subroutine eigen_matrix_x

    subroutine eigen_matrix_z(u, ml, mr)
        ! Left/right characteristic matrices for the z direction.
        real(precision), intent(in) :: u(4)
        real(precision), intent(out) :: ml(4, 4), mr(4, 4)

        real(precision) :: ro, uu, ww, ps, c, c2, alfa

        ro = u(1)
        uu = u(2)/ro
        ww = u(3)/ro
        alfa = 0.5_precision*(uu*uu + ww*ww)
        ps = beta*(u(4) - alfa*ro)
        c2 = gama*ps/ro
        c = sqrt(c2)

        ml(1, 1) = alfa*beta - c2
        ml(1, 2) = -beta*uu
        ml(1, 3) = -beta*ww
        ml(1, 4) = beta
        ml(2, 1) = -uu
        ml(2, 2) = 1.0_precision
        ml(2, 3) = 0.0_precision
        ml(2, 4) = 0.0_precision
        ml(3, 1) = alfa*beta - c*ww
        ml(3, 2) = -beta*uu
        ml(3, 3) = c - beta*ww
        ml(3, 4) = beta
        ml(4, 1) = alfa*beta + c*ww
        ml(4, 2) = -beta*uu
        ml(4, 3) = -c - beta*ww
        ml(4, 4) = beta

        mr(1, 1) = -2.0_precision
        mr(2, 1) = -2.0_precision*uu
        mr(3, 1) = -2.0_precision*ww
        mr(4, 1) = -2.0_precision*alfa
        mr(1, 2) = 0.0_precision
        mr(2, 2) = 2.0_precision*c2
        mr(3, 2) = 0.0_precision
        mr(4, 2) = uu*2.0_precision*c2
        mr(1, 3) = 1.0_precision
        mr(2, 3) = uu
        mr(3, 3) = ww + c
        mr(4, 3) = alfa + c*ww + c2/beta
        mr(1, 4) = 1.0_precision
        mr(2, 4) = uu
        mr(3, 4) = ww - c
        mr(4, 4) = alfa - c*ww + c2/beta
        mr = mr*0.5_precision/c2
    end subroutine eigen_matrix_z

end module e2d_flux
