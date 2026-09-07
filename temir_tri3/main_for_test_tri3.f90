!!! Module set 
module parameters
    implicit none

    ! Specify the limit of model size
    integer, parameter :: nnode_max = 1000          
    integer, parameter :: nelem_max = 2000
    integer, parameter :: nbc_max = 10

contains
    !!! This is the functions to convert values from marc_style_exponent natation to real(8) notation. 
    function expo2double(field) result(value)       
        character(len = 256), intent(in) :: field
        integer :: expo, pos
        real(8) :: base, value

        pos = scan(field, "+-", BACK=.true.)        ! Buscar sign
        read(field(:pos-1),*) base                  ! Extract base of value
        read(field(pos:),*) expo                    ! Extract exponent of value

        value =base * 10.0d0**expo                  ! Convert to double
    end function expo2double
end module parameters

!!! Main of the program 
program main
    use parameters
    implicit none

    !Declear variables for read_geometry
    integer :: datfile, nnode, nelem                ! datfile : input file, nnode : the number of nodes, nelem : the number of element
    integer :: connect(nelem_max, 3)                ! connect : connectivity of each element 
    real(8) :: coord(nnode_max, 2)                  ! coord : coordinate of each node
    real(8) :: E, nu, t                             ! E : Young modulous, nu = poisson ration, t = thickness of plane stress plate

    ! Declear variables for read_BC
    integer :: num_bc_set, bc_node_set(nbc_max,1000), fixed_disp_vector(nbc_max,3), num_node_in_set(nbc_max)
    real(8) :: fixed_disp_magn(nbc_max,3), point_load_magn(nbc_max,2)

    !Declear variable for make_D
    real(8) :: D(3,3)                               ! D : D-matrix 

    ! Declear variables for make_B
    real(8) :: B(nelem_max,3,6), area(nelem_max)    ! B : B-matrix

    ! Declear variables for make_elem_stiffness
    real(8) :: elem_stiffness(nelem_max,6,6)

    ! Declear variabls for stiffness
    real(8) :: stiffness(2*nnode_max, 2*nnode_max)

    ! Delear variables for apply_BC
    real(8) :: reduced_stiffness(2*nnode_max,2*nnode_max), force(2*nnode_max), reduced_force(2*nnode_max)

    ! Declear variables for this test_read_dat_1.f90
    integer :: i, j, k
    
    call read_geometry(datfile, connect, coord, nnode, nelem, E, nu, t)

    call read_BC(datfile, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn)

    call make_D(D, E, nu)
    
    call make_B(connect, coord, nelem, B, area)

    call make_elem_stiffness(nelem, D, B, area, t, elem_stiffness)

    call assemble_stiffness(nnode, nelem, connect, elem_stiffness, stiffness)

    call apply_BC(nnode, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn, &
                  stiffness, reduced_stiffness, force, reduced_force)

   !print *, "connectivity"
   !do i = 1, nelem
   !    do j = 1, 3
   !        print *, connect(i, j)
   !    end do
   !    print *
   !end do
   !
   !print *, "coordinates"
   !do i = 1, nnode
   !    do j = 1, 2
   !        print *, coord(i, j)
   !    end do
   !    print *
   !end do

   !print *, "Boundary Conditions"
   !do i = 1, num_bc_set
   !    print *, "The list of node to applied"
   !    print *, bc_node_set(i, 1:num_node_in_set(i))
   !    do j = 1, 2
   !        print *, "fixed_disp_magnitude"
   !        print *, fixed_disp_magn(i, j)
   !        print *, "fixed_disp_vector"
   !        print *, fixed_disp_vector(i, j)
   !        print *, "point load magnitude"
   !        print *, point_load_magn(i, j)
   !    end do
   !end do

   !print *, "E = "
   !print *, E
   !print *, "nu = "
   !print *, nu

   !print *, "D-matrix"
   !do i=1, 3
   !    do j=1, 3
   !        print *, D(i, j)
   !    end do
   !    print *
   !end do

   !print *, "B-matrix"
   !do i=1, nelem
   !    do j=1, 3
   !        print *, B(i,j,:)
   !    end do
   !    print *
   !end do

   !print *, "thickness"
   !print *, "t = ", t    

   !do k = 1, nelem
   !    print *, "K",k
   !    do i = 1, 6
   !        do j = 1, 6
   !            print *,elem_stiffness(k, i, j)*8.0/75000.0
   !        end do
   !        print *
   !    end do
   !end do 

   !print *, "stiffness matrix"
   !do i = 1, 2*nnode
   !    do j = 1, 2*nnode
   !        print *, stiffness(i, j)*8.0/75000.0
   !    end do
   !    print *
   !end do

    print *, "=== Reduced system: K * x = f ==="
    do i = 1, 2*nnode
        print *, "row", i
        do j = 1, 2*nnode
            write(*, '(F12.4)', advance="no") reduced_stiffness(i, j) *8.0d0/75000.d0
        end do
        write(*, '(A, F12.4)') "  |  f = ", reduced_force(i)
        print *
    end do

end program main