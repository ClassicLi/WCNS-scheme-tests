program test_initial
    ! Unit tests for the 1D Euler case definitions and initial states.
    use precision_select
    use e1d_globals, only: halo
    use e1d_initial
    use test_util
    implicit none

    real(precision) :: rl, ul, pl, rr, ur, pr
    real(precision) :: u8(-halo:8 + halo, 3)
    real(precision) :: u10(-halo:10 + halo, 3)

    ! Case selection.
    call check_true(e1d_case_allowed('sod'), 'sod allowed')
    call check_true(e1d_case_allowed('lax'), 'lax allowed')
    call check_true(e1d_case_allowed('shu_osher'), 'shu_osher allowed')
    call check_true(.not. e1d_case_allowed('sedov'), 'sedov rejected')
    call check_true(.not. e1d_case_allowed('tt'), 'tt rejected')

    ! Domain length and final time.
    call check_close(e1d_case_length('sod'), 1.0_precision, 0.0_precision, 'sod length')
    call check_close(e1d_case_length('shu_osher'), 10.0_precision, 0.0_precision, &
                     'shu length')
    call check_close(e1d_case_tend('sod'), 0.15_precision, 0.0_precision, 'sod tend')
    call check_close(e1d_case_tend('shu_osher'), 1.8_precision, 0.0_precision, &
                     'shu tend')

    ! Sod and Lax left/right states.
    call e1d_case_states('sod', rl, ul, pl, rr, ur, pr)
    call check_close(rl, 1.0_precision, 0.0_precision, 'sod rho_l')
    call check_close(rr, 0.125_precision, 0.0_precision, 'sod rho_r')
    call check_close(pr, 0.1_precision, 0.0_precision, 'sod p_r')

    call e1d_case_states('lax', rl, ul, pl, rr, ur, pr)
    call check_close(rl, 0.445_precision, 0.0_precision, 'lax rho_l')
    call check_close(ul, 0.698_precision, 0.0_precision, 'lax u_l')
    call check_close(pl, 3.528_precision, 0.0_precision, 'lax p_l')
    call check_close(rr, 0.5_precision, 0.0_precision, 'lax rho_r')
    call check_close(pr, 0.571_precision, 0.0_precision, 'lax p_r')

    ! Sod grid: discontinuity between j = n/2 and j = n/2 + 1.
    call e1d_initial_state(8, 'sod', u8)
    call check_close(u8(0, 1), 1.0_precision, 0.0_precision, 'sod j=0 rho')
    call check_close(u8(4, 1), 1.0_precision, 0.0_precision, 'sod j=n/2 rho')
    call check_close(u8(5, 1), 0.125_precision, 0.0_precision, 'sod j=n/2+1 rho')
    call check_close(u8(-7, 1), 1.0_precision, 0.0_precision, 'sod left halo')
    call check_close(u8(15, 1), 0.125_precision, 0.0_precision, 'sod right halo')
    call check_close(u8(0, 3), 2.5_precision, 1.0e-14_precision, 'sod E = p/beta')

    ! Shu-Osher grid: the right segment overwrites j = n/10.
    call e1d_initial_state(10, 'shu_osher', u10)
    call check_close(u10(0, 1), 3.857143_precision, 0.0_precision, 'shu left rho')
    call check_close(u10(1, 1), 1.0_precision + 0.2_precision*sin(5.0_precision), &
                     1.0e-12_precision, 'shu j=n/10 is right state')
    call check_close(u10(5, 1), 1.0_precision + 0.2_precision*sin(25.0_precision), &
                     1.0e-12_precision, 'shu sinusoid')
    call check_close(u10(1, 3), 2.5_precision, 1.0e-14_precision, 'shu E = 1/0.4')

    call finish()

end program test_initial
