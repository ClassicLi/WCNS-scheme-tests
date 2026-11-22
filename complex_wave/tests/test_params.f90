program test_params
    ! Unit test for the complex-wave control parameters and operator validation.
    use precision_select
    use cw_control, only: mesh_num, periods, cfl, output_dir, cw_control_validate
    use cw_params, only: difference_operator, difference_operator_allowed, &
                         validate_difference_params
    use test_util
    implicit none

    integer :: ierr
    character(len=256) :: msg

    ! Valid defaults must pass every check.
    mesh_num = 200
    periods = 4
    cfl = 0.3_precision
    output_dir = 'results'
    call cw_control_validate(ierr, msg)
    call check_true(ierr == 0, 'defaults valid')

    ! The period count must be a positive integer.
    periods = 0
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'periods=0 rejected')
    periods = -1
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'periods<0 rejected')
    periods = 4

    ! A nonpositive CFL would give a zero or negative time step.
    cfl = 0.0_precision
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'cfl=0 rejected')
    cfl = -0.1_precision
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'cfl<0 rejected')
    cfl = 0.3_precision

    ! The initial condition uses ten-percent segments, so 20 is the floor.
    mesh_num = 19
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'mesh<20 rejected')
    mesh_num = 201
    call cw_control_validate(ierr, msg)
    call check_true(ierr == 0, 'mesh=201 accepted')
    mesh_num = 200

    ! The output directory must be named.
    output_dir = ''
    call cw_control_validate(ierr, msg)
    call check_true(ierr /= 0, 'empty output_dir rejected')
    output_dir = 'results'

    ! Operator names are case-sensitive and must be in the allowed set.
    call check_true(difference_operator_allowed('HCS6'), 'HCS6 allowed')
    call check_true(difference_operator_allowed('WCNSE6'), 'WCNSE6 allowed')
    call check_true(.not. difference_operator_allowed('NOPE'), 'NOPE rejected')
    difference_operator = 'NOPE'
    call validate_difference_params(ierr, msg)
    call check_true(ierr /= 0, 'validate unknown operator')

    call finish()
end program test_params
