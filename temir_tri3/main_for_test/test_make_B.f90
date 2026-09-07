module parameters
    implicit none
    integer, parameter :: nnode_max = 1000
    integer, parameter :: nelem_max = 2000
end module parameters

program test_make_B
    use parameters
    implicit none

    integer :: i,j
    integer :: nnode, nelem
    integer :: connect(nelem_max, 3)
    real(8) :: coord(nnode_max, 2)
    real(8) :: B(nelem_max,3,6)

    nnode = 4
    nelem = 2
    
    connect(1,1) = 1
    connect(1,2) = 4
    connect(1,3) = 3
    connect(2,1) = 4
    connect(2,2) = 1
    connect(2,3) = 2

    coord(1,1) = 0
    coord(1,2) = 0
    coord(2,1) = 2
    coord(2,2) = 0
    coord(3,1) = 0
    coord(3,2) = 1
    coord(4,1) = 2
    coord(4,2) = 1

    call make_B(connect, coord, nelem, B)

    do i=1, nelem
        do j=1, 3
            print *, B(i,j,:)
        end do
        print *
    end do

end program test_make_B

