subroutine read_geometry(datfilename, datfile, nnode, nelem, connect, coord, E, nu, t)
    use parameters
    implicit none

    !!!!! Get .dat file and extract model info and boundary conditions.

    ! Output arguments
    character(len=256), intent(out) :: datfilename
    integer, intent(out) :: datfile
    integer, intent(out) :: nnode, nelem
    integer, intent(out) :: connect(nelem_max, 3)
    real(8), intent(out) :: coord(nnode_max, 2)
    real(8), intent(out) :: E, nu, t

    ! Local variables in read_geometry
    character(len=256) :: line, field
    integer :: ios_nelem, ios_nnode, elem_type, ios_elem_type
    integer :: i, j, l, m, n

    ! Local variables to find error
    integer :: ios_sizing, ios_connectivity, ios_coordinate, ios_isotropic, ios_geom
 
    !!! Specify .dat file to input
    call get_command_argument(1,datfilename)                        ! Get the name of input file
    print *, "Input file = ",trim(datfilename)
    open(newunit=datfile, file = trim(datfilename), status="old", action="read")


    !!! Get geometry info of model
    ! Get model size
    do 
        read(datfile,"(A)",iostat = ios_sizing) line
        if (ios_sizing /= 0) then
            print *, "!!!!! ERROR : Sizing Block Not Found !!!!!"
            error stop
        end if

        if (line(1:10) == "sizing") then                            ! Buscar "sizing"
            read(line(41:50),*,iostat=ios_nelem) nelem              ! Store number of elements
            read(line(51:60),*,iostat=ios_nnode) nnode              ! Store number of nodes
            
            ! Error Section 
            !! いずれ節点数や要素数の制限も設けたい．
            if(ios_nelem == 0 .and. ios_nnode == 0) then
                print *, "number of element = ", nelem
                print *, "number of node = ", nnode
            else
                print *, "!!!!! ERROR : Size data read !!!!!"
                print *, "line = ",line
                error stop
            end if

            exit
        end if
    end do
  
    ! Store element connectivity
    do 
        read(datfile,"(A)", iostat = ios_connectivity) line
        if (ios_connectivity /= 0) then
            print *, "!!!!! ERROR : Connectivity Block Not Found !!!!!"
            error stop
        end if

        if(line(1:12) == "connectivity") then                       ! Buscar "connectivity"
            read(datfile,*)

            do m = 1, nelem
                read(datfile,"(A)") line
                read(line(16:20),*,iostat=ios_elem_type) elem_type
                    if(elem_type == 201) then                       ! Confirm type of elemnt == 201

                        read(line(6:10),*) i                        ! Store element number
                                                                    ! The number of line "m" has possible to differ from element number "i".
                        read(line(26:30),*) connect(i, 1)           ! Store number of 1st node 
                        read(line(36:40),*) connect(i, 2)           ! Store number of 2nd node
                        read(line(46:50),*) connect(i, 3)           ! Store number of 3rd node

                    else
                        read(line(6:10),*) i                        ! Store element number                        
                        print *, "!!!!! Unexpected type of element : element number = ", i, ", element type = ", elem_type
                        error stop
                    end if
            end do
            exit
        end if
    end do

    ! Store node coordinate
    do 
        read(datfile,"(A)", iostat = ios_coordinate) line
        if (ios_coordinate /= 0) then
            print *, "!!!!! ERROR : Coordinate Block Not Found !!!!!"
            error stop
        end if

        if(line(1:11) == "coordinates") then                        ! Buscar "coordinates"
            read(datfile,*)
            do n = 1, nnode 
                read(datfile,"(A)") line
                read(line(6:10),*) j                                ! Store node number
                                                                    ! The number of line "n" has possible to differ from node number "j".
                do l = 1, 2
                    read(line(20*l-9:20*l+10),*) field              ! Read value of node in each axis
                    coord(j, l) = expo2double(field)                ! Calcular value of node in each axis
                end do
            end do
            exit
        end if
    end do


    !!! Get material properties
    !!暫定的なやつ．仕様書がよくわからない．
    do 
        read(datfile,"(A)", iostat = ios_isotropic) line
        if (ios_isotropic /= 0) then
            print *, "!!!!! ERROR : Isotropic Block Not Found !!!!!"
            error stop
        end if

        if(line(1:9) == "isotropic") then
            read(datfile,*)
            read(datfile,*)
            read(datfile,*)
            
            read(datfile,"(A)") line
            read(line(1:20),*) field
            E = expo2double(field)                                  ! Convert notation by expo2double
            
            read(line(21:40),*) field
            nu = expo2double(field)                                 ! Convert notation by expo2double
            exit
        end if
    end do


    !!! Get geometry info
    do
        read(datfile,"(A)", iostat = ios_geom) line
        if(ios_geom /= 0) then
            print *, "!!!!! ERROR : Geometry Block Not Found !!!!!"
            error stop
        end if

        if (line(1:16) == "geometry") then
            read(datfile,*)
            read(datfile,*)
            read(datfile,*)

            read(datfile,"(A)") line
            read(line(1:20),*) field
            t = expo2double(field)
            exit
        end if
    end do

end subroutine read_geometry