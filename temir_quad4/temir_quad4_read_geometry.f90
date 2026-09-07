subroutine read_geometry(datfile, nnode, nelem, connect, coord, E, nu, t)

    use parameters
    implicit none

    ! Input arguments
    integer, intent(in) :: datfile

    ! Output arguments
    integer, intent(out) :: nnode, nelem
    real(8), intent(out) :: E, nu, t
    integer, allocatable :: connect(:, :)
    real(8), allocatable :: coord(:, :)

    ! Loacl variables
    character(len=256) :: line
    integer :: ios_sizing, ios_nelem, ios_nnode


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

    
    
end subroutine