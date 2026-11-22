program invalid_choice_driver
    use precision_select
    use interp_framework, only: MOC33
    use scheme_params, only: moc33_linear
    implicit none

    real(precision) :: u(-3:4), u_half
    integer :: i

    do i = -3, 4
        u(i) = 1.0_precision
    end do
    u(0) = 0.0_precision
    moc33_linear = 'BOGUS'
    call MOC33(u, u_half)
    write(*, '(a,es12.4)') 'unexpected result: ', u_half
end program invalid_choice_driver
