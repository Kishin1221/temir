subroutine make_elem_stiffness(B, detJ, nelem, D, t, elem_stiffness)

    use parameters
    implicit none

    !!!!! Program to calculate element stiffness matrix.
    !!! K_e = t * (Be^T)De(Be) * (jacobian at gauss point) * (weight at gauss point)

    ! Input arguments
    real(8), allocatable, intent(in) :: B(:,:,:,:), detJ(:,:)
    integer, intent(in) :: nelem
    real(8), intent(in) :: D(3,3), t

    ! Output arguments
    real(8), allocatable, intent(out) :: elem_stiffness(:,:,:)

    ! Local variables
    integer :: ielem, igp, i, j, p, q

    ! Declare matrix size
    allocate (elem_stiffness(nelem, 8, 8))

    ! Initialize
    elem_stiffness = 0.0d0


    !!! Calculate entries of elem_stiffness
    do ielem = 1, nelem                                     ! Loop of elements
        do igp = 1, 4                                       ! Loop of gauss point
            do j = 1, 8                                     ! Loop of free index for 3rd degree of elem_stiffness 
                do i = 1, 8                                 ! Loop of free index for 2nd degree of elem_stiffness
                    do q = 1, 3                             ! Loop of 2nd dummy index 
                        do p = 1, 3                         ! Loop of 1st dummy index
                            elem_stiffness(ielem, i, j) = elem_stiffness(ielem, i, j) + B(ielem, igp, p, i) * D(p, q) * B(ielem, igp, q, j) * detJ(ielem, igp) * gauss_wt(igp)
                        end do
                    end do
                end do
            end do
        end do
    end do

    elem_stiffness = t * elem_stiffness                     ! Thickness

end subroutine make_elem_stiffness


subroutine assemble_stiffness(connect, elem_stiffness, nelem, nnode, stiffness)

    use parameters
    implicit none

    !!!!! Program to make global stiffness matrix by assembling element stiffness matrix. 
    !!! Decide direction of degree of freedum by judging mod2 of index.

    ! Input arguments
    integer, allocatable, intent(in) :: connect(:,:)
    real(8), allocatable, intent(in) :: elem_stiffness(:,:,:)
    integer, intent(in) :: nelem, nnode

    ! Output arguments
    real(8), allocatable, intent(out) :: stiffness(:,:)

    ! Local variables
    integer :: idof_u, idof_f
    integer :: ielem, i, j

    ! Declare matrix size
    allocate (stiffness(2*nnode, 2*nnode))

    ! Initialize
    stiffness = 0.0d0


    ! Assemble elem_stiffness to stiffness
    do ielem = 1, nelem
        do j = 1, 8
            do i = 1, 8
                if (mod(i, 2) == 1) then                                        ! The 1st index "i" is odd number
                    idof_u = 2 * connect(ielem, (i+1)/2) - 1
                
                else if (mod(i, 2) == 0) then                                   ! The 1st index "i" is even number
                    idof_u = 2 * connect(ielem, i/2)
                
                end if

                if (mod(j, 2) == 1) then                                        ! The 2nd index "j" is odd number
                    idof_f = 2 * connect(ielem, (j+1)/2) - 1
                
                else if (mod(j, 2) == 0) then                                   ! The 2nd index "j" is even number
                    idof_f = 2 * connect(ielem, j/2)
                
                end if

                stiffness(idof_u, idof_f) = stiffness(idof_u, idof_f) + elem_stiffness(ielem, i, j)
            
            end do
        end do
    end do

end subroutine assemble_stiffness