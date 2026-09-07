!!!!!!! "temir" is a FEM program for elastic homogeneous material. !!!!!!!

!!!!! Program specifications !!!!!

!!! To avoid ERROR !!!

!!! For compile !!!


!!! Module set
module parameters
    implicit none
   
contains
    !!! This is the function to convert values from marc_style_exponent notation to real(8) notation.
    function expodouble(field) result(value)
        character(len = *), intent(in) :: field
        integer :: expo, pos
        real(8) :: base, value

        pos = scan(field, "+-", BACK=.true.)
        read(field(:pos-1),*) base
        read(field(pos:),*) expo

        value = base * 10.0d0**expo

    end function expodouble
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
