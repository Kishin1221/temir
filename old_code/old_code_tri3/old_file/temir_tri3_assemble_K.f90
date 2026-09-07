subroutine assemble_stiffness(nnode, nelem, connect, Ke, stiffness)
    use parameters
    implicit none

    !!!!! 
    
    ! Declear variables
    integer, intent(in) :: nnode, nelem, connect(nelem_max, 3)
    real(8), intent(in) :: elem_stiffness(nelem_max,6,6)
    real(8), intent(out) :: stiffness(2*nnode_max, 2*nnode_max)

    ! Declear local variables
    integer :: ielem, i, j, idof_f, idof_u


    ! Assemble Ke to stiffness
    stiffness = 0.0d0

    do ielem = 1, nelem
        do i = 1, 6
            do j = 1, 6
                if (mod(i, 2) == 1) then                    ! i is odd number
                    idof_u = 2 * connect(ielem, (i+1)/2) -1
                else if (mod(i, 2) == 0) then               ! i is even number
                    idof_u = 2 * connect(ielem, i/2)
                end if

                if (mod(j, 2) == 1) then                    ! j is odd number
                    idof_f = 2 * connect(ielem, (j+1)/2) -1
                else if (mod(j, 2) == 0) then               ! j is even number
                    idof_f = 2 * connect(ielem, j/2)
                end if

                stiffness(idof_u, idof_f) = stiffness(idof_u, idof_f) + elem_stiffness(ielem, i, j)
            end do
        end do
    end do

end subroutine assemble_stiffness