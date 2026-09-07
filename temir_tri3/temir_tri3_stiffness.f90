subroutine make_elem_stiffness(nelem, D, B, area, t, elem_stiffness)

    use parameters
    implicit none

    !!!!! Get D and B matrix from main and make element_stifness matrix.

    ! Input arguments
    integer, intent(in) :: nelem
    real(8), intent(in) :: D(3,3), B(nelem_max,3,6), area(nelem_max), t
    
    ! Output arguments
    real(8), intent(out) :: elem_stiffness(nelem_max,6,6)

    ! Local variables in make_elem_stiffness
    integer :: ielem, i, j, k
    real(8) :: DB(3,6)

    !!! Calculate elem_stiffness matrix of ielem'th element
    do ielem = 1, nelem
        DB = 0.0d0                                          ! Initialize
        elem_stiffness(ielem, :, :) = 0.0d0                 ! Initialize

        do i = 1, 3
            do j = 1, 6
                do k = 1, 3
                    DB(i, j) = DB(i, j) + D(i, k) * B(ielem, k, j)    ! Calculate B*D of (B-transpose)DB first
                end do
            end do
        end do

        do i = 1, 6
            do j = 1, 6
                do k = 1, 3
                    elem_stiffness(ielem, i, j) = elem_stiffness(ielem, i, j) + (B(ielem, k, i) * DB(k, j))  ! Calculate (B-transpose)DB
                end do
                elem_stiffness(ielem, i, j) = elem_stiffness(ielem, i, j) * t * area(ielem)                         ! Maltiple thickness and area
            end do
        end do
    end do

end subroutine make_elem_stiffness


subroutine assemble_stiffness(nelem, connect, elem_stiffness, stiffness)

    use parameters
    implicit none

    !!!!! Get element_stiffness matrix and assemble in global_stiffness matrix.
    
    ! Input arguments
    integer, intent(in) :: nelem, connect(nelem_max, 3)
    real(8), intent(in) :: elem_stiffness(nelem_max,6,6)
    real(8), intent(out) :: stiffness(2*nnode_max, 2*nnode_max)

    ! Local variables in assemble_stiffness
    integer :: ielem, i, j, idof_f, idof_u

    ! Assemble elem_stiffness to stiffness
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