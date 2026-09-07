module parameters
    implicit none
    integer, parameter :: nnode_max = 1000
    integer, parameter :: nelem_max = 2000
end module parameters

program test_strain_stress
    use parameters
    implicit none
    
    integer :: i,j
    integer :: nelem
    integer :: connect(nelem_max, 3)
    real(8) :: coord(nnode_max, 2)
    real(8) :: E, nu
    real(8) :: D(3,3)
    real(8) :: B(nelem_max,3,6)
    real(8) :: disp(2*nnode_max)
    real(8)  :: strain(nelem_max, 3)
    real(8)  :: stress(nelem_max, 3)
    
    E = 100.d0
    nu = 1.0d0/3.0d0
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

    call make_D(D, E, nu)

    call make_B(connect, coord, nelem, B)

    disp(1) = 0
    disp(2) = 0
    disp(3) = 6.0d0/100.0d0
    disp(4) = 0
    disp(5) = 0
    disp(6) = -1.0d0/100.0d0
    disp(7) = 6.0d0/100.0d0
    disp(8) = -1.0d0/100.0d0

    call strain_stress(nelem, connect, D, B, disp, strain, stress)
      
    do i=1, nelem
        do j=1, 3
            print *, strain(i, j)
        end do
        print *
    end do

    do i=1, nelem
        do j=1, 3
            print *, stress(i, j)
        end do
        print *
    end do

end program test_strain_stress