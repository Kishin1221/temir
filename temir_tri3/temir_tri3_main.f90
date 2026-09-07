!!!!!!! "temir" is a FEM program for elastic homogeneous material. !!!!!!!

!!!!! Program Specifications !!!!! 
! This version "temir_tri3" only accept marc element type 201.
! Numbers of elements and nodes must be continuous.
! Max size of model is specified in module. !!!!!
! Running command should be "~/path/temir_tri3_main.exe modeldatfil.dat". !

!!! To Avoid Error (The situations below are not expected in this program.) !!! 
! All boundary conditions defined in mentat must be used in job. !

!!! For compile
! gfortran temir_tri3_main.f90 temir_tri3_read_geometry.f90 temir_tri3_read_BC.f90 temir_tri3_DBmatrices.f90 temir_tri3_stiffness.f90 temir_tri3_apply_BC.f90 temir_tri3_solver.f90 temir_tri3_post_process.f90 temir_tri3_export_result.f90 -o temir_tri3.exe 


!!! Module set 
module parameters
    implicit none

    ! Specify the limit of model size
    integer, parameter :: nnode_max = 256
    integer, parameter :: nelem_max = 256
    integer, parameter :: nbc_max = 32

contains
    !!! This is the functions to convert values from marc_style_exponent natation to real(8) notation. 
    function expo2double(field) result(value)       
        character(len = *), intent(in) :: field
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

    ! Variables for read_geometry
    character(len=256) :: datfilename
    integer :: datfile, nnode, nelem                ! datfile : input file, nnode : the number of nodes, nelem : the number of element
    integer :: connect(nelem_max, 3)                ! connect : connectivity of each element 
    real(8) :: coord(nnode_max, 2)                  ! coord : coordinate of each node
    real(8) :: E, nu, t                             ! E : Young modulous, nu = poisson ration, t = thickness of plane stress plate

    ! Variables for read_BC
    integer :: num_bc_set, bc_node_set(nbc_max,1000), fixed_disp_vector(nbc_max,3), num_node_in_set(nbc_max)
    real(8) :: fixed_disp_magn(nbc_max,3), point_load_magn(nbc_max,2)

    ! Variables for make_D
    real(8) :: D(3,3)                               ! D : D-matrix 

    ! Variables for make_B
    real(8) :: B(nelem_max,3,6), area(nelem_max)    ! B : B-matrix

    ! Variables for make_elem_stiffness
    real(8) :: elem_stiffness(nelem_max,6,6)

    ! Variables for assemble_stiffness
    real(8) :: stiffness(2*nnode_max, 2*nnode_max)

    ! Variables for apply_BC
    integer :: shift_index(2*nnode_max), num_not_fixed
    logical :: is_fixed(2*nnode_max)
    real(8) :: reduced_stiffness(2*nnode_max,2*nnode_max)
    real(8) :: force(2*nnode_max), reduced_force(2*nnode_max)

    ! Variables for solver
    real(8) :: reduced_disp(2*nnode_max)

    ! Variables for make_full_displacement
    real(8) :: disp(2*nnode_max)

    ! Variables for reaction_force
    real(8) :: reaction_force(2*nnode_max)

    ! Variables for strain_stress
    real(8) :: strain(nelem_max,3), stress(nelem_max,3)

    ! Local variable
    integer :: i


    !!! Start Sbroutines    
    call read_geometry(datfilename, datfile, nnode, nelem, connect, coord, E, nu, t)

    call read_BC(datfile, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn)

    call make_D(D, E, nu)

    call make_B(nelem, connect, coord, B, area)

    call make_elem_stiffness(nelem, D, B, area, t, elem_stiffness)

    call assemble_stiffness(nelem, connect, elem_stiffness, stiffness)

    call apply_BC(nnode, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn, &
                  stiffness, shift_index, num_not_fixed, is_fixed, reduced_stiffness, force, reduced_force)

    call solver(num_not_fixed, reduced_stiffness, reduced_force, reduced_disp)

    call make_full_displacement(nnode, num_bc_set, bc_node_set, num_node_in_set, shift_index, is_fixed, &
                                fixed_disp_magn, reduced_disp, disp)

    call reaction(nnode, is_fixed, stiffness, disp, reaction_force)

    call strain_stress(nelem, connect, D, B, disp, strain, stress)

    call export_result(datfilename, nnode, nelem, coord, E, nu, t, disp, reaction_force, strain, stress)

end program