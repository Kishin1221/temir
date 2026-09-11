subroutine read_geometry(datfile, nnode, nelem, connect, coord, E, nu, t)

    use parameters
    implicit none

    ! Input arguments
    integer, intent(in) :: datfile

    ! Output arguments
    integer, intent(out) :: nnode, nelem
    real(8), intent(out) :: E, nu, t
    integer, allocatable, intent(out) :: connect(:, :)
    real(8), allocatable, intent(out) :: coord(:, :)

    ! Loacl variables
    character(len=256) :: line, field
    integer :: elem_type
    integer :: ios_sizing, ios_nelem, ios_nnode, ios_connectivity, ios_coordinate, ios_isotropic, ios_geom
    integer :: i, j, l, m, n


    !!! Get model size
    do 
        read(datfile, "(A)", iostat = ios_sizing) line
        if (ios_sizing /= 0) then
            
            print *, "ERROR : Sizng Block Not Found !!!!!"
            error stop
        
        end if
        
        if (line(1:10) == "sizing") then                    ! Buscar "sizing"
            read(line(41:50), *, iostat = ios_nelem) nelem  ! Store the number of element 
            read(line(51:60), *, iostat = ios_nnode) nnode  ! Store the number of node

            ! Error section
            if (ios_nelem == 0 .and. ios_nnode == 0) then
                print *, "The number of element = ", nelem
                print *, "The number of node    = ", nnode
            else
                print *, "!!!!! ERROR : Model size cannot be read"
                print *, "line = ", line
                error stop
            end if

            exit
        end if
    end do

    ! Declare matrix size
    allocate(connect(nelem, 4))
    allocate(coord(nnode, 3))

    ! Store element connectivity
    do 
        read(datfile, "(A)", iostat = ios_connectivity) line
        if (ios_connectivity /= 0) then
            print *, "!!!!! ERROR : Connectivity Block Not Found !!!!!"
            error stop
        end if            

        ! declear matrix size here

        if (line(1:12) == "connectivity") then
            read(datfile, *)

            do m = 1, nelem
                read(datfile, "(A)") line
                read(line(16:20), *, iostat = ios_elem_type) elem_type

                if (elem_type == 3) then                    ! Confirm the type of the element == 3
                    read(line(6:10), *) i                   ! Store element number
                    
                    read(line(26:30),*) connect(i, 1)           ! Store number of 1st node 
                    read(line(36:40),*) connect(i, 2)           ! Store number of 2nd node
                    read(line(46:50),*) connect(i, 3)           ! Store number of 3rd node
                    read(line(56:60),*) connect(i, 4)           ! Store number of 4th node
                
                else
                    read(line(6:10), *) i
                    print *, "!!!!! ERROR : Unexpected element type !!!!!"
                    print *, " element number = ", i, ", element type = ", elem_type
                    error stop
                
                end if
            end do
            exit
        end if
    end do

    ! Store node coordinate
    do
        read(datfile, "(A)", iostat = ios_coordinate) line
        if (ios_coordinate /= 0) then
            print *, "!!!!! ERROR : Coordinate Block Not Found !!!!!"
            error stop
        end if

        if (line(1:11) == "coordinates") then               ! Buscar "coordinate"
            read(datfile,*)

            do n = 1, nnode
                read(datfile, "(A)") line
                read(line(6:10), *) j                       ! Store node number
                                                            ! The number of line "n" has possible to differ from node number "j".
                do l = 1, 3
                    read(line(20*l-9:20*l+10), *) field     ! Read value of coordinate
                    coord(j, l) = expo2double(field)        ! Store the value of coordinate
                end do
            end do
            exit
        end if
    end do
    

    !!! Get material properties
    do
        read(datfile, "(A)", iostat = ios_isotropic) line
        if (ios_isotropic /= 0) then
            print *, "ERROR : Isotropic Block Not Found !!!!!"
            error stop
        end if

        if (line(1:9) == "isotropic") then
            read(datfile, *)
            read(datfile, *)
            read(datfile, *)

            read(datfile, "(A)") line
            read(line(1:20), *) field
            E = expo2double(field)

            read(line(21:40), *) field
            nu = expo2double(field)
            exit
        end if
    end do


    !!! Get geometry info
    do
        read(datfile, "(A)", iostat = ios_geom) line
        if (ios_geom /= 0) then
            print *, "!!!!! ERROR : Geometr Block Not Found !!!!!"
            error stop
        end if

        if (line(1:16) == "geometry") then
            read(datfile, *)
            read(datfile, *)
            read(datfile, *)
    
            read(datfile, "(A)") line
            read(line(1:20), *) field
            t = expo2double(field)
            exit
        end if
    end do

end subroutine