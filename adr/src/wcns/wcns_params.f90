module wcns_params
    ! Difference-operator selection and validation for the WCNS ADR driver.
    use precision_select
    implicit none

    ! One of WCNSE2, WCNSE4, WCNSE6, HCS4, HCS6.
    character(len=32) :: difference_operator = 'HCS4'

    namelist /difference_nml/ difference_operator

contains

    subroutine read_difference_namelist(unit, ierr)
        ! Read the &difference_nml group from an open input file.
        integer, intent(in) :: unit
        integer, intent(out), optional :: ierr

        integer :: ios

        rewind(unit)
        read(unit, nml=difference_nml, iostat=ios)
        if (present(ierr)) ierr = ios
    end subroutine read_difference_namelist

    function difference_operator_allowed(name) result(ok)
        ! Return true for a known difference-operator ID.
        character(len=*), intent(in) :: name
        logical :: ok

        select case (trim(name))
        case ('WCNSE2', 'WCNSE4', 'WCNSE6', 'HCS4', 'HCS6')
            ok = .true.
        case default
            ok = .false.
        end select
    end function difference_operator_allowed

    subroutine validate_wcns_params(ierr, msg)
        ! Check that the selected difference operator is known.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. difference_operator_allowed(difference_operator)) then
            ierr = 1
            msg = 'unknown difference operator: ' // trim(difference_operator)
        end if
    end subroutine validate_wcns_params

end module wcns_params
