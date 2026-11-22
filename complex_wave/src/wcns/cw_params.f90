module cw_params
    ! Difference-operator selection for the complex-wave test.
    !
    ! The half-node fluxes come from a WCNS interpolation scheme; this module
    ! selects how those fluxes (and, for the hybrid operators, the integer-node
    ! values) are combined into the spatial right-hand side.
    use precision_select
    implicit none

    ! One of WCNSE2, WCNSE4, WCNSE6, HCS4, HCS6.
    character(len=32) :: difference_operator = 'HCS6'

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

    subroutine validate_difference_params(ierr, msg)
        ! Check that the selected difference operator is known.
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: msg

        ierr = 0
        msg = ''

        if (.not. difference_operator_allowed(difference_operator)) then
            ierr = 1
            msg = 'unknown difference operator: ' // trim(difference_operator)
        end if
    end subroutine validate_difference_params

end module cw_params
