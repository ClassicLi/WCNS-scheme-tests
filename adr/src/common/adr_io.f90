module adr_io
    ! CSV, log, and optional Tecplot output for the ADR drivers.
    use, intrinsic :: ieee_arithmetic
    use precision_select
    implicit none

contains

    subroutine write_csv(path, varphi, phi, n)
        ! Write varphi, Re(Phi), Im(Phi) as comma-separated values.
        character(len=*), intent(in) :: path
        integer(integer_precision), intent(in) :: n
        real(precision), intent(in) :: varphi(0:n)
        complex(precision), intent(in) :: phi(0:n)

        integer :: u
        integer(integer_precision) :: j

        open(newunit=u, file=trim(path), status='replace', action='write')
        write(u, '(a)') 'varphi,phi_re,phi_im'
        do j = 0, n
            write(u, '(es24.16,a,es24.16,a,es24.16)') varphi(j), ',', &
                real(phi(j), precision), ',', aimag(phi(j))
        end do
        close(u)
    end subroutine write_csv

    subroutine write_log(path, scheme_name, operator_name, n_num, dt, varphi, phi, n)
        ! Write run metadata, curve extrema, and a NaN count.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: scheme_name
        character(len=*), intent(in) :: operator_name
        integer(integer_precision), intent(in) :: n_num
        real(precision), intent(in) :: dt
        real(precision), intent(in) :: varphi(0:n)
        complex(precision), intent(in) :: phi(0:n)
        integer(integer_precision), intent(in) :: n

        integer :: u
        integer(integer_precision) :: j, n_nan
        real(precision) :: re_min, re_max, im_min, im_max, v_min, v_max

        n_nan = 0
        do j = 0, n
            if (ieee_is_nan(real(phi(j), precision)) .or. ieee_is_nan(aimag(phi(j)))) then
                n_nan = n_nan + 1
            end if
        end do

        re_min = minval(real(phi(0:n), precision))
        re_max = maxval(real(phi(0:n), precision))
        im_min = minval(aimag(phi(0:n)))
        im_max = maxval(aimag(phi(0:n)))
        v_min = minval(varphi(0:n))
        v_max = maxval(varphi(0:n))

        open(newunit=u, file=trim(path), status='replace', action='write')
        write(u, '(a)') 'version=1'
        write(u, '(a)') 'scheme=' // trim(scheme_name)
        write(u, '(a)') 'operator=' // trim(operator_name)
        write(u, '(a,i0)') 'mesh_num=', n_num
        write(u, '(a,es24.16)') 'dt=', dt
        write(u, '(a,es24.16)') 'varphi_min=', v_min
        write(u, '(a,es24.16)') 'varphi_max=', v_max
        write(u, '(a,es24.16)') 'phi_re_min=', re_min
        write(u, '(a,es24.16)') 'phi_re_max=', re_max
        write(u, '(a,es24.16)') 'phi_im_min=', im_min
        write(u, '(a,es24.16)') 'phi_im_max=', im_max
        write(u, '(a,i0)') 'nan_count=', n_nan
        close(u)
    end subroutine write_log

    subroutine write_tecplot(path, name, varphi, phi, n)
        ! Write a block-format Tecplot sample of the ADR curve.
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: name
        real(precision), intent(in) :: varphi(0:n)
        complex(precision), intent(in) :: phi(0:n)
        integer(integer_precision), intent(in) :: n

        integer :: u
        integer(integer_precision) :: j

        open(newunit=u, file=trim(path), status='replace', action='write')
        write(u, '(a)') 'Title="ADR ' // trim(name) // '"'
        write(u, '(a)') 'Variables = varphi, phi_re, phi_im'
        write(u, '(a)') 'ZONE T="' // trim(name) // '"'
        write(u, '(a,i0,a)') 'I=', n + 1, ',J=1,K=1,f=block'
        do j = 0, n
            write(u, '(es24.16)') varphi(j)
        end do
        do j = 0, n
            write(u, '(es24.16)') real(phi(j), precision)
        end do
        do j = 0, n
            write(u, '(es24.16)') aimag(phi(j))
        end do
        close(u)
    end subroutine write_tecplot

end module adr_io
