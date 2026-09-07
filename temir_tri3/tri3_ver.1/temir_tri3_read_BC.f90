subroutine read_BC(datfile, num_bc_set, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn)
    use parameters
    implicit none

    !!!!! Read the same file as subroutine read_geometry, and make vector of External Force and Given Displacement.

    ! Input arguments
    integer, intent(in) :: datfile

    ! Output arguments
    integer, intent(out) :: num_bc_set, bc_node_set(nbc_max,1000)
    integer, intent(out) :: fixed_disp_vector(nbc_max,3), num_node_in_set(nbc_max)
    real(8), intent(out) :: fixed_disp_magn(nbc_max,3), point_load_magn(nbc_max,2)

    ! Local variables(common in read_bc)
    character(len=256) :: line
    character(len=20)  :: field, tokens(100), bc_set_name(nbc_max)
    integer :: bc_set_id
    integer :: i, j, l

    ! Local variables used in "define"
    character(len=256) :: bc_set_type
    integer :: ini_node, fin_node, temp_node, token_count, set_node_count
    logical :: find_to

    ! Declear variables used in fixed disp and point load
    integer :: num_stored_bc
    real(8) :: magnitude

    !Declear temporary variable used in fixed disp and point load
    character(len=20) :: temp_bc_name(nbc_max)
    integer :: temp_fixed_disp_vector(nbc_max,3)
    real(8) :: temp_fixed_disp_magn(nbc_max,3), temp_point_load_magn(nbc_max,2)

    ! Declear variables to find error
    integer :: ios_define, ios_temp_node, ios_token, ios_buscar

    
    !!! Go to the top of .dat input file
    rewind(datfile)
    bc_set_id = 1

    !!! Skip to first "define"
    do                                             
        read(datfile,"(A)", iostat = ios_define) line
        if (ios_define /= 0) then
            print *, "!!!!! ERROR : Blocks of define Not Found !!!!!"
            error stop
        end if

        if (line(1:12) == "define") then                    ! Find a block of "define"
            exit
        end if
    end do

    !!! Get "set" info of boundary conditions in define blocks
    do 
        read(line(21:40),*) bc_set_type                     ! Check the type of set of boundary condition
        
        if (bc_set_id > nbc_max) then
            print*, "!!!!! ERROR : Boundary conditions exist more than !!!!!", nbc_max
            error stop
        else
            ! No problem
        end if

        ! Confirm the type of set is "node". This program accept displacement or POINT LOAD as a boundary conditions.
        if (bc_set_type == "node") then
            ! No problem
        else
            print *, "!!!!! ERROR : THe type of boundary conditions are wrong !!!!! "
            error stop
        end if

        read(line(61:80),*) bc_set_name(bc_set_id)          ! Store the name of set of BC 
        set_node_count = 1

        ! Loop A1 !
        do 
            read(datfile,"(A)") line
            tokens = "@"
            read(line,*,iostat = ios_token) tokens                              ! Sprit line to block and memorize as array of "tokens"
            token_count = count(tokens /= "@")                                  ! Count th number of entry of "tokens"

            find_to = any(tokens(1:token_count) == "to")                        ! Continuas nodes are expressed with "to".
            if (find_to) then
                read(tokens(1),*) ini_node                                      ! Store initial node number of "to" list
                read(tokens(3),*) fin_node                                      ! Store final node number of "to" list
                num_node_in_set(bc_set_id) = fin_node - ini_node + 1

                do l = 1, num_node_in_set(bc_set_id)                                           
                    bc_node_set(bc_set_id, l) = ini_node + l - 1                ! Store node number to apply BC
                end do

            else 
                ! Loop A2 !
                do i = 1, token_count
                    if (tokens(i) == "c") then                                  ! "c" means bc_node_set continue to the next line.
                        exit
                    else
                        read(tokens(i),*,iostat=ios_temp_node) temp_node        ! Try to store token as Integer
                        if (ios_temp_node == 0) then                            ! Integer means bc_node_set continue to next token. 
                            bc_node_set(bc_set_id, set_node_count) = temp_node  ! Store node number to apply BC
                            set_node_count = set_node_count + 1
                        else                                                    ! This line should contain nothing except integer or "c".
                            print *, "!!!!! ERROR : Failed to read bc_set !!!!!"    
                            error stop
                        end if
                    end if
                end do
                num_node_in_set(bc_set_id) = set_node_count -1                  ! Store the number of nodes included in BC set
            
                if(tokens(token_count) == "c") then         ! If last token is "c", go to next line and repeat Loop A2.
                    ! Read next line
                else
                    exit                                    ! If last token is integer, this "define" block end.
                end if
            end if
        end do
        
        bc_set_id = bc_set_id + 1                                   
        read(datfile,"(A)") line
        if (line(1:12) == "define") then                    ! If next line start with "define", repeat Loop A1.
            ! Continuar do loop
        else                                                ! All of "define" has finished.
            exit
        end if
    end do

    num_bc_set = bc_set_id - 1                              ! This is used to specify times of loop to find bc_set_name
    bc_set_id = 1                                           ! Already used in upper section. Need to be reset.
    num_stored_bc = 0                                       ! This is used to judge whether all BC have been read.
    temp_fixed_disp_magn = 0.0d0                            ! Initialize
    temp_fixed_disp_vector = 0                              ! Initialize
    fixed_disp_magn = 0.0d0                                 ! Initialize
    fixed_disp_vector = 0                                   ! Initialize
    temp_point_load_magn = 0.0d0                            ! Initialize
    point_load_magn = 0.0d0                                 ! Initialize

    ! Loop B1 !
    Buscar_All_BC : do
    ! Try to find all BC by searching the characters of "fixed disp" or "point load".
    ! Each block of BC (both of the blocks of "fixed disp" and "point load") consist of 6 lines in total. 
        read(datfile,"(A)",iostat = ios_buscar) line
        if (ios_buscar /= 0) then
            print *, "!!!!! ERROR : The Boundary Conditions not found !!!!!"
            error stop
        end if

        !!! Read fixed disp
        if (line(1:12) == "fixed disp") then                ! Find a block of "fixed disp"
            read(datfile,*)                                 ! Skip the next line of the characters of "fixed disp". 
            !Now the next line is the 1st line of the 1st "fixed disp" block.

            ! Loop B2(fixed disp) !
            do
                num_stored_bc = num_stored_bc + 1           ! Count up the number of BC which are already read.
                read(datfile, *)
                read(datfile,"(A)") line
                token_count = len_trim(line) / 20           ! Count the number of vector to fixed by mesuring the length of line

                do i = 1, token_count
                    read(line(20*i-19:20*i), "(A)") field
                    magnitude = expo2double(field)
                    temp_fixed_disp_magn(num_stored_bc, i) = magnitude
                end do

                read(datfile, *)

                read(datfile, "(A)") line                   ! Read line of disp vector
                tokens = "@"                                ! Initialize
                read(line,*, iostat = ios_token) tokens

                do j = 1, token_count
                    read(tokens(j),*) temp_fixed_disp_vector(num_stored_bc, j)          ! Store vector of BC
                end do               
                
                read(datfile,*)
                
                read(datfile,"(A)") line
                temp_bc_name(num_stored_bc) = line(1:20)    ! Store the name of BC (This is used to map define block and BC block later)

                do l = 1, num_bc_set                                                            ! "l" is just a counter to search the same name.
                    if (temp_bc_name(num_stored_bc) == bc_set_name(l)) then                     ! Find the set of "define" which has the same name.
                        do i = 1, 3
                            fixed_disp_magn(l, i) = temp_fixed_disp_magn(num_stored_bc, i)      ! Assign
                            fixed_disp_vector(l, i) = temp_fixed_disp_vector(num_stored_bc, i)  ! Plug in 
                        end do
                        exit
                    end if
                end do                
        
                                read(datfile,"(A)") line
                backspace(datfile)                                  ! Read the next line
                if (line(1:12) == "point load") then                ! Judge whether fixed disp end or still exist below
                    exit
                end if

                if (num_stored_bc == num_bc_set) then               ! Judgement of end of BC to read
                    exit Buscar_All_BC
                else
                    ! Continuar reading BC
                end if                
            end do

        !!! Read point load
        else if (line(1:12) == "point load") then                   ! Find a block of "point load" 
            read(datfile,*)
            
            ! Loop B2(point load) !
            do
                num_stored_bc = num_stored_bc + 1
                read(datfile,*)
                read(datfile,"(A)") line        
                token_count = len_trim(line) / 20

                do i = 1, token_count
                    read(line(20*i-19:20*i), "(A)") field
                    magnitude = expo2double(field)
                    temp_point_load_magn(num_stored_bc, i) = magnitude
                end do                

                read(datfile,*)
                read(datfile,*)

                read(datfile,"(A)") line
                temp_bc_name(num_stored_bc) = line(1:20)            ! Get the name of this BC
                do l = 1, num_bc_set
                    if (temp_bc_name(num_stored_bc) == bc_set_name(l)) then                     ! Find the set of "define" which has the same name. 
                        do i = 1, 2
                            point_load_magn(l, i) = temp_point_load_magn(num_stored_bc, i)      ! Assign
                        end do
                        exit
                    end if
                end do

                ! Count up the number of BC which are already read.
                if (num_stored_bc == num_bc_set) then               ! Judgement of end of BC to read
                    exit Buscar_All_BC
                else
                    ! Continuar reading BC
                end if               
            end do
        end if
    end do Buscar_All_BC                                            ! End of boundary conditions

    ! Close input dat file
    close(datfile)

end subroutine read_BC                