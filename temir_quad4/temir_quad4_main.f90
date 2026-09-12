!!!!!!! "temir" is a FEM program for elastic homogeneous material. !!!!!!!

!!!!! Program specifications !!!!!

!!! To avoid ERROR !!!

!!! For compile !!!


!!! Module set
module parameters
    implicit none
    
    !!! Parameters for gauss Integration (2×2 for quad)
    ! Nodes of master elements
    real(8), parameter :: node_xi(4) = [-1.0d0, 1.0d0, 1.0d0, -1.0d0]
    real(8), parameter :: node_eta(4) = [-1.0d0, -1.0d0, 1.0d0, 1.0d0]

    ! Nodes and weights of gauss point
    integer, parameter :: ngauss = 4
    real(8), parameter :: gauss_pt(4,2) = reshape([-0.5773502692d0, 0.5773502692d0, 0.5773502692d0, -0.5773502692d0, &
                                            -0.5773502692d0, -0.5773502692d0, 0.5773502692d0, 0.5773502692d0], shape=[4,2])
    real(8), parameter :: gauss_wt(4) = [1.0d0, 1.0d0, 1.0d0, 1.0d0]
   
contains
    !!! This is the function to convert values from marc_style_exponent notation to real(8) notation.
    function expo2double(field) result(value)
        character(len = *), intent(in) :: field
        integer :: expo, pos
        real(8) :: base, value

        pos = scan(field, "+-", BACK=.true.)
        read(field(:pos-1),*) base
        read(field(pos:),*) expo

        value = base * 10.0d0**expo

    end function expo2double
end module parameters


!!! Main of the program
program main
    use parameters
    implicit none

    ! Variables to open input datfile
    character(len=256) ::datfilename
    integer :: datfile

    ! Variables for read_geometry
    integer :: nelem, nnode
    real(8) :: E, nu, t
    integer, allocatable :: connect(:, :)
    real(8), allocatable :: coord(:, :)

    ! Open input datfile
    call get_command_argument(1, datfilename)
    print *, "Input file = ", trim(datfilename)
    open(newunit = datfile, file = trim(datfilename), status = "old", action = "read")

    ! Start subroutine
    call read_geometry(datfile, nnode, nelem, connect, coord, E, nu, t)

    call read_BC(datfile, nnode, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn, num_bc_set)

    call make_D(D, E, nu)

    call derivative_shape_function(dNdxi, dNdeta)

    call make_B(connect, coord, nelem, dNdxi, dNdeta, B, detJ)
    
