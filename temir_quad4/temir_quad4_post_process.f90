subroutine make_full_displacement(bc_node_set, num_node_in_set, shift_index, is_fixed, fixed_disp_magn, &
                                    reduced_disp, nnode, num_bc_set, disp)

    use parameters
    implicit none

    ! Input arguments
    logical, allocatable, intent(in) :: is_fixed(:)
    integer, allocatable, intent(in) :: bc_node_set(:,:), num_node_in_set(:), shift_index(:) 
    real(8), allocatable, intent(in) :: fixed_disp_magn(:,:), reduced_disp(:)
    integer, intent(in) :: nnode, num_bc_set
    
    ! Output arguments
    real(8), allocatable, intent(out) :: disp(:)

    ! Local variables
    integer :: i, k

    ! Declare matrix size
    allocate (disp(2*nnode))
    
    ! Initialize
    disp = 0.0d0

    !!! Expand reduced_disp to disp
    ! Copy disp values of not-fixed degree of freedom
    do i =1, 2*nnode
        if (is_fixed(i) .eqv. .false.) then
            disp(i) = reduced_disp(i - shift_index(i))
        end if
    end do

    ! Copy disp value of degree of freedom which is fixed by fixed_disp
    do k = 1, num_bc_set
        if (fixed_disp_magn(k, 1) /= 0.0d0) then
            do i = 1, num_node_in_set(k)
                disp(2 * bc_node_set(k, i) - 1) = fixed_disp_magn(k, 1)
            end do
        end if
        
        if (fixed_disp_magn(k, 2) /= 0.0d0) then
            do i =1, num_node_in_set(k)
                disp(2 * bc_node_set(k, i)) = fixed_disp_magn(k,2)
            end do
        end if
    end do
    
end subroutine make_full_displacement

subroutine reaction(is_fixed, stiffness, disp, nnode, reaction_force)

    use parameters
    implicit none

    ! Input arguments
    logical, allocatable, intent(in) :: is_fixed(:)
    real(8), allocatable, intent(in) :: stiffness(:,:), disp(:)
    integer, intent(in) :: nnode

    ! Output arguments
    real(8), allocatable, intent(out) :: reaction_force(:)

    ! Local variables
    integer :: i, j

    ! Declare matrix size
    allocate (reaction_force(2*nnode))

    ! Initialize
    reaction_force = 0.0d0


    !!! Vamos a calcular
    do i =1, 2*nnode
        if (is_fixed(i)) then
            do j = 1, 2*nnode
                reaction_force(i) = reaction_force(i) + stiffness(i, j) *disp(j)
            end do
        end if
    end do

end subroutine reaction

subroutine strain_stress()

    use parameters
    implicit none

    ! Inpput arguments
    integer, allocatable, intent(in) :: connect(:,:)
    real(8), allocatable, intent(in) :: B(:,:,:,:), disp(:)
    integer, intent(in) :: nelem
    real(8), intent(in) :: D(3,3)
    
    ! Output arguments
    real(8), allocatable, intent(out) :: strain(:,:,:), stress(:,:,:)
    
    ! Local variables
    integer :: ielem
    real(8) :: local_disp(8)


    ! Declare matrices size
    allocate (strain(nelem, ngauss, 3))
    allocate (stress(nelem, ngauss, 3))




end subroutine strain_stress