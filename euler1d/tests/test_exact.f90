program test_exact
    ! Unit tests for the Sod/Lax exact Riemann solution.
    use precision_select
    use e1d_initial, only: e1d_case_states
    use e1d_exact
    use test_util
    implicit none

    real(precision) :: rl, ul, pl, rr, ur, pr
    real(precision) :: ps, us, rsl, rsr
    real(precision) :: rho, u, p
    integer :: ierr
    character(len=256) :: msg

    ! Star states of the two supported cases.
    call e1d_case_states('sod', rl, ul, pl, rr, ur, pr)
    call e1d_riemann_star(rl, ul, pl, rr, ur, pr, ps, us, rsl, rsr, ierr, msg)
    call check_true(ierr == 0, 'sod star converged')
    call check_close(ps, 0.30313017805064674_precision, 1.0e-8_precision, 'sod p_star')
    call check_close(us, 0.9274526200489498_precision, 1.0e-8_precision, 'sod u_star')
    call check_close(rsl, 0.4263194281784951_precision, 1.0e-8_precision, &
                     'sod rho_star_l')
    call check_close(rsr, 0.2655737117053071_precision, 1.0e-8_precision, &
                     'sod rho_star_r')

    ! A deliberately truncated iteration must report failure, so the
    ! caller does not use an unconverged star state.
    call e1d_riemann_star(rl, ul, pl, rr, ur, pr, ps, us, rsl, rsr, ierr, msg, &
                          max_iter=1)
    call check_true(ierr /= 0, 'star reports non-convergence')

    call e1d_case_states('lax', rl, ul, pl, rr, ur, pr)
    call e1d_riemann_star(rl, ul, pl, rr, ur, pr, ps, us, rsl, rsr, ierr, msg)
    call check_true(ierr == 0, 'lax star converged')
    call check_close(ps, 2.4660979192073564_precision, 1.0e-8_precision, 'lax p_star')
    call check_close(us, 1.528723026632884_precision, 1.0e-8_precision, 'lax u_star')
    call check_close(rsl, 0.34456847418960945_precision, 1.0e-8_precision, &
                     'lax rho_star_l')
    call check_close(rsr, 1.3040845320261996_precision, 1.0e-8_precision, &
                     'lax rho_star_r')

    ! Wave coordinate: the initial discontinuity is at x = 0.5.
    call e1d_riemann_sample(1.0_precision, 0.0_precision, 1.0_precision, &
                            0.125_precision, 0.0_precision, 0.1_precision, &
                            0.15_precision, 0.2_precision, rho, u, p)
    call check_close(rho, 1.0_precision, 1.0e-12_precision, 'sod x=0.2 left rho')
    call e1d_riemann_sample(1.0_precision, 0.0_precision, 1.0_precision, &
                            0.125_precision, 0.0_precision, 0.1_precision, &
                            0.15_precision, 0.9_precision, rho, u, p)
    call check_close(rho, 0.125_precision, 1.0e-12_precision, 'sod x=0.9 right rho')

    ! t = 0 reduces to the initial piecewise state.
    call e1d_riemann_sample(1.0_precision, 0.0_precision, 1.0_precision, &
                            0.125_precision, 0.0_precision, 0.1_precision, &
                            0.0_precision, 0.25_precision, rho, u, p)
    call check_close(rho, 1.0_precision, 1.0e-14_precision, 't=0 left')
    call e1d_riemann_sample(1.0_precision, 0.0_precision, 1.0_precision, &
                            0.125_precision, 0.0_precision, 0.1_precision, &
                            0.0_precision, 0.75_precision, rho, u, p)
    call check_close(rho, 0.125_precision, 1.0e-14_precision, 't=0 right')

    ! Double shock.
    call e1d_riemann_star(1.0_precision, 1.0_precision, 1.0_precision, &
                          1.0_precision, -1.0_precision, 1.0_precision, &
                          ps, us, rsl, rsr, ierr, msg)
    call check_true(ierr == 0 .and. ps > 1.0_precision, 'double shock p*')
    call check_close(us, 0.0_precision, 1.0e-10_precision, 'double shock u*')

    ! Double rarefaction, non-vacuum.
    call e1d_riemann_star(1.0_precision, -1.0_precision, 1.0_precision, &
                          1.0_precision, 1.0_precision, 1.0_precision, &
                          ps, us, rsl, rsr, ierr, msg)
    call check_true(ierr == 0 .and. ps > 0.0_precision .and. ps < 1.0_precision, &
                    'double rarefaction p*')

    ! Vacuum.
    call e1d_riemann_star(1.0_precision, -6.0_precision, 1.0_precision, &
                          1.0_precision, 6.0_precision, 1.0_precision, &
                          ps, us, rsl, rsr, ierr, msg)
    call check_close(ps, 0.0_precision, 1.0e-14_precision, 'vacuum p*')
    call e1d_riemann_sample(1.0_precision, -6.0_precision, 1.0_precision, &
                            1.0_precision, 6.0_precision, 1.0_precision, &
                            0.1_precision, 0.5_precision, rho, u, p)
    call check_true(rho == 0.0_precision .and. p == 0.0_precision, 'vacuum center')

    ! Left shock / right rarefaction.
    call e1d_riemann_star(1.0_precision, 0.0_precision, 0.1_precision, &
                          1.0_precision, 0.0_precision, 1.0_precision, &
                          ps, us, rsl, rsr, ierr, msg)
    call check_true(ierr == 0 .and. ps > 0.1_precision .and. ps < 1.0_precision, &
                    'left shock p*')

    call finish()

end program test_exact
