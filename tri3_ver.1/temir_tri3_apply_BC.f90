subroutine apply_BC(nnode, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn, &
                    stiffness, shift_index, num_not_fixed, is_fixed, reduced_stiffness, force, reduced_force)

    use parameters
    implicit none

    !!!!! Get BC and global_stiffness matrix, and assemble in equations.

    ! Input arguments
    integer, intent(in) :: nnode, num_bc_set, bc_node_set(nbc_max,1000), num_node_in_set(nbc_max)
    integer, intent(in) :: fixed_disp_vector(nbc_max,3) 
    real(8), intent(in) :: fixed_disp_magn(nbc_max,3), point_load_magn(nbc_max,2)
    real(8), intent(in) :: stiffness(2*nnode_max, 2*nnode_max)

    ! Output arguments
    integer, intent(out) :: shift_index(2*nnode_max), num_not_fixed
    logical, intent(out) :: is_fixed(2*nnode_max)
    real(8), intent(out) :: reduced_stiffness(2*nnode_max,2*nnode_max)
    real(8), intent(out) :: force(2*nnode_max), reduced_force(2*nnode_max)

    ! Local variables in apply_BC
    integer :: i, j, k, trans_dof

    !!! Find fixed degrees of freedom 
    ! The index of is_fixed is corresponding to degrees of freedom of all nodes.
    is_fixed = .false.                                      ! Initialize : assume that all degrees of freedom is not fixed
    do k = 1, num_bc_set
        do i = 1, 3
            if (fixed_disp_vector(k, i) == 1) then          ! When displacement of x is fixed
                do j = 1, num_node_in_set(k)
                    is_fixed(2 * bc_node_set(k, j) - 1) = .true.                ! .true. means this degree of freedom is fixed
                end do
                
            else if(fixed_disp_vector(k, i) == 2) then      ! When displacement of y is fixed
                do j = 1, num_node_in_set(k)
                    is_fixed(2 * bc_node_set(k, j)) = .true.                    ! .true. means this degree of freedom is fixed
                end do

            end if
        end do
    end do

    ! Count the number of DOF fixed by BC of fixed_disp so far
    ! shift_index is defined to original index, not to after reductions. 英語が怪しいな
    shift_index = 0                                         ! Initialize
    do i = 1, 2*nnode
        if (is_fixed(i) .eqv. .true.) then
            if (i == 1) then
                shift_index(i) = 1                          ! When i = 1, shift_index(i-1) does not exist.
            
            else 
                shift_index(i) = shift_index(i-1) + 1       ! shift_index incrase
        
            end if
        
        else
            if (i == 1) then
                shift_index(i) = 0                          ! When i = 1, shift_index(i-1) does not exist.

            else
                shift_index(i) = shift_index(i-1)
            
            end if

        end if
    end do

    ! Count the number of free DOF
    num_not_fixed = 0
    do i = 1, 2*nnode
        if (.not. is_fixed(i)) then
            num_not_fixed = num_not_fixed + 1
        end if
    end do

    ! Reduction of stiffness matrix
    reduced_stiffness = 0.0d0                               ! Initialize
    do i = 1, 2*nnode
        do j = 1, 2* nnode
            if (.not. is_fixed(i) .and. .not. is_fixed(j)) then                                     ! If both of row and colmun are not fixed
                reduced_stiffness(i - shift_index(i), j - shift_index(j)) = stiffness(i, j)         ! Reduction. Store the entry of stiffness matrix
            
            else
                ! Do not store

            end if
        end do
    end do

    ! Apply BC of fixed_disp to global force vector
    force = 0.0d0
    do k = 1, num_bc_set
        if (point_load_magn(k, 1) /= 0.0d0) then            !!! Caution : point_load_magn is not expected to include round-off error.
            do i = 1, num_node_in_set(k)
                force(2 * bc_node_set(k, i) - 1) = force(2 * bc_node_set(k, i) - 1) + point_load_magn(k, 1)
            end do
        end if

        if (point_load_magn(k, 2) /= 0.0d0) then
            do i =1, num_node_in_set(k)
                force(2 * bc_node_set(k, i)) = force(2 * bc_node_set(k, i)) + point_load_magn(k, 2)
            end do
        end if
    end do

    ! Transposition
    do k = 1, num_bc_set                                    !!! Caution : fixed_disp_magn is not expected to include round-off error.
        if (fixed_disp_magn(k, 1) /= 0.0d0) then            ! If a displacement given for x-axis is not 0, need transposition.
            do i = 1, num_node_in_set(k)
                trans_dof = 2 * bc_node_set(k, i) - 1
                do j = 1, 2*nnode
                    force(j) = force(j) - stiffness(j, trans_dof) * fixed_disp_magn(k, 1)
                end do
            end do
        end if

    if (fixed_disp_magn(k, 2) /= 0.0d0) then                ! If a displacement given for y-axis is not 0, need transposition. 
            do i = 1, num_node_in_set(k)
                trans_dof = 2 * bc_node_set(k, i)
                do j = 1, 2*nnode
                    force(j) = force(j) -stiffness(j, trans_dof) * fixed_disp_magn(k, 2)
                end do
            end do
        end if
    end do

    ! Reduction of force vector
    reduced_force = 0.0d0
    do j = 1, 2*nnode
        if (is_fixed(j) .eqv. .false.) then                    ! If this DOF is not fixed
            reduced_force(j - shift_index(j)) = force(j)    ! Reduction
        
        else                                                ! If this DOF is fixed
            ! Do not store
        
        end if
    end do

end subroutine apply_BC