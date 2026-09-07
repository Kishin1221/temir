subroutine make_full_displacement(nnode, num_bc_set, bc_node_set, num_node_in_set, shift_index, is_fixed, &
                                  fixed_disp_magn, reduced_disp, disp)

    use parameters
    implicit none

    ! Input arguments
    integer, intent(in) :: nnode, num_bc_set, bc_node_set(nbc_max,1000), num_node_in_set(nbc_max)
    integer, intent(in) :: shift_index(2*nnode_max)
    logical, intent(in) :: is_fixed(2*nnode_max)    
    real(8), intent(in) :: fixed_disp_magn(nbc_max,3), reduced_disp(2*nnode_max)

    ! Output arguments
    real(8), intent(out) :: disp(2*nnode_max)

    ! Local variables in make_full_displacement
    integer :: i, k

    !!! Expand reduced_disp to disp
    disp = 0.0d0

    do i = 1, 2*nnode
        if (is_fixed(i) .eqv. .false.) then
            disp(i) = reduced_disp(i - shift_index(i))
        end if
    end do

    do k = 1, num_bc_set
        if (fixed_disp_magn(k, 1) /= 0.0d0) then
            do i = 1, num_node_in_set(k)
                disp(2 * bc_node_set(k, i) - 1) = fixed_disp_magn(k, 1)
            end do
        end if

        if (fixed_disp_magn(k, 2) /= 0.0d0) then
            do i = 1, num_node_in_set(k)
                disp(2 * bc_node_set(k, i)) = fixed_disp_magn(k, 2)
            end do
        end if
    end do

end subroutine make_full_displacement

subroutine reaction(nnode, is_fixed, stiffness, disp, reaction_force)

    use parameters
    implicit none

    ! Input arguments
    integer, intent(in) :: nnode
    logical, intent(in) :: is_fixed(2*nnode_max)
    real(8), intent(in) :: stiffness(2*nnode_max, 2*nnode_max), disp(2*nnode_max)

    ! Output arguments
    real(8), intent(out) :: reaction_force(2*nnode_max)

    ! Local variables in reaction_force
    integer :: i, j

    !!! Vamos a calcular
    reaction_force = 0.0d0

    do i = 1, 2*nnode
        if (is_fixed(i)) then
            do j = 1, 2*nnode
                reaction_force(i) = reaction_force(i) + stiffness(i, j) * disp(j)
            end do
        end if
    end do

end subroutine reaction

subroutine strain_stress(nelem, connect, D, B, disp, strain, stress)

    use parameters
    implicit none

    ! Input arguments
    integer, intent(in) :: nelem, connect(nelem_max, 3)
    real(8), intent(in) :: D(3,3), B(nelem_max,3,6)
    real(8), intent(in) :: disp(2*nnode_max)

    ! Output arguments
    real(8), intent(out) :: strain(nelem_max, 3), stress(nelem_max, 3)

    ! Local variables in strain_stress
    integer :: ielem, i, j, k
    real(8) :: local_disp(6)

    !!! Vamos a calcular
    do ielem = 1, nelem
        local_disp = 0.0d0
        strain(ielem, :) = 0.0d0
        stress(ielem, :) = 0.0d0

        do i =1, 3
            local_disp(2*i - 1) = disp(2*connect(ielem, i) -1)                  ! x-axis
            local_disp(2*i) = disp(2*connect(ielem, i))                         ! y-axis
        end do

        do j = 1, 3
            do k = 1, 6
                strain(ielem, j) = strain(ielem, j) + B(ielem, j, k) * local_disp(k)
            end do
        end do

        do j = 1, 3
            do k = 1, 3
                stress(ielem, j) = stress(ielem, j) + D(j, k) * strain(ielem, k)
            end do
        end do
    end do

end subroutine strain_stress



