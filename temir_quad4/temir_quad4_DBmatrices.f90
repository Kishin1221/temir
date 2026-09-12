subroutine make_D(D, E, nu)

    use parameters
    implicit none

    ! Input arguments
    real(8), intent(in) :: E, nu

    ! Output argument
    real(8), intent(out) :: D(3,3)
    
    ! Local variable
    real(8) :: factor

    ! Initialize
    D = 0.0d0

    !!! This program only accept plane stress model now.
    ! Caluculate D matrix
    factor = E / (1-nu**2)

    D(1,1) = 1
    D(1,2) = nu
    D(2,1) = nu
    D(2,2) = 1
    D(3,3) = (1-nu) / 2

    D = D * factor

end subroutine make_D

subroutine derivative_shape_function(dNdxi, dNdeta)

    use parameters
    implicit none

    ! Output arguments
    real(8), intent(out) :: dNdxi(4,4), dNdeta(4,4)

    !Local variables
    integer :: igp, inode

    do igp = 1, ngauss                                      ! In fortran, row loop should be arrenged outside.
        do inode = 1, 4
            dNdxi(inode, igp)  = 0.25d0 * node_xi(inode)  * (1.0d0 + node_eta(inode) * gauss_pt(igp, 2))
            dNdeta(inode, igp) = 0.25d0 * node_eta(inode) * (1.0d0 + node_xi(inode)  * gauss_pt(igp, 1))
        end do
    end do


end subroutine derivative_shape_function


subroutine make_B(connect, coord, nelem, dNdxi, dNdeta, B, detJ)

    use parameters
    implicit none

    ! Input arguments
    integer, allocatable, intent(in) :: connect(:,:)
    real(8), allocatable, intent(in) :: coord(:,:)
    integer, intent(in) :: nelem
    real(8), intent(in) :: dNdxi(4,4), dNdeta(4,4)

    ! Output arguments
    real(8), allocatable, intent(out) :: B(:,:,:,:), detJ(:,:)

    ! Local variables
    integer :: ielem, inode, igp
    real(8) :: J(2,2), H(2,2), dNdx(4), dNdy(4)

    ! Declare matrix size
    allocate (B(nelem, ngauss, 3, 8))
    allocate (detJ(nelem, ngauss))
    B = 0.0d0

    !!! Start
    do ielem = 1, nelem
        do igp = 1, ngauss
            
            J = 0.0d0
            dNdx = 0.0d0
            dNdy = 0.0d0

            do inode = 1, 4
                J(1, 1) = J(1, 1) + coord(connect(ielem, inode), 1) * dNdxi(inode, igp)
                J(2, 1) = J(2, 1) + coord(connect(ielem, inode), 2) * dNdxi(inode, igp)
                J(1, 2) = J(1, 2) + coord(connect(ielem, inode), 1) * dNdeta(inode, igp)
                J(2, 2) = J(2, 2) + coord(connect(ielem, inode), 2) * dNdeta(inode, igp)
            end do 

            detJ(ielem, igp) = J(1,1) * J(2, 2) - J(1, 2) * J(2, 1)

            H(1,1) =  J(2,2) / detJ(ielem, igp)
            H(1,2) = -J(1,2) / detJ(ielem, igp)
            H(2,1) = -J(2,1) / detJ(ielem, igp)
            H(2,2) =  J(1,1) / detJ(ielem, igp)

            do inode = 1, 4
                dNdx(inode) = dNdxi(inode, igp) * H(1, 1) + dNdeta(inode, igp) * H(2, 1)
                dNdy(inode) = dNdxi(inode, igp) * H(1, 2) + dNdeta(inode, igp) * H(2, 2)

                B(ielem, igp, 1, 2*inode-1) = dNdx(inode)
                B(ielem, igp, 2, 2*inode-1) = 0.0d0
                B(ielem, igp, 3, 2*inode-1) = dNdy(inode)
                B(ielem, igp, 1, 2*inode)   = 0.0d0
                B(ielem, igp, 2, 2*inode)   = dNdy(inode)
                B(ielem, igp, 3, 2*inode)   = dNdx(inode)
            end do
        end do
    end do

end subroutine make_B