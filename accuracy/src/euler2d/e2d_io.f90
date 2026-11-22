module e2d_io
    ! Optional Tecplot sample output for the 2D accuracy tests.
    use precision_select
    use e2d_globals, only: n, m, h, ubar
    implicit none

contains

    subroutine e2d_write_sample(path)
        ! Write the density field as a point-packed Tecplot zone.
        character(len=*), intent(in) :: path

        integer :: unit, j, k

        open (newunit=unit, file=trim(path), status='replace', action='write')
        write (unit, '(a)') 'title="2D accuracy sample"'
        write (unit, '(a)') 'variables="x","y","rho"'
        write (unit, '(a,i0,a,i0,a)') 'Zone I=', n + 1, ' J=', m + 1, &
            ', datapacking=POINT'
        do k = 0, m
            do j = 0, n
                write (unit, '(3(es24.16,1x))') real(j, precision)*h, &
                    real(k, precision)*h, ubar(j, k, 1)
            end do
        end do
        close (unit)
    end subroutine e2d_write_sample

end module e2d_io
