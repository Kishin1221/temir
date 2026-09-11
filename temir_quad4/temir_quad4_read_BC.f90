subroutine read_BC(datfile, nnode, bc_node_set, num_node_in_set, fixed_disp_vector, fixed_disp_magn, point_load_magn, num_bc_set)

    use parameters
    implicit none

    !!!!! Read input datfile and Store boundary conditions. 
    
    !!! Algorithm
    ! According to the structure of .dat file, informations about boundary conditions are separated in different blocks.
    ! In early section, blocks of "define" exist, then blocks of "fixed disp" or "point load" follow.

    ! In this program, we first read "define" and get informations below.
        ! the name of BC "bc_set_name"
        ! node numbers to apply this BC "bc_node_set"
        ! The number of elements of bc_node_in_set "num_node_in_bc_set"

    ! Then we find "fixed disp" or "point load".


 
    ! Input arguments
    integer, intent(in) :: datfile, nnode

    ! Output arguments
    integer, allocatable, intent(out) :: bc_node_set(:,:), num_node_in_set(:), fixed_disp_vector(:,:) 
    real(8), allocatable, intent(out) :: fixed_disp_magn(:,:), point_load_magn(:,:)
    integer, intent(out) :: num_bc_set
    
    ! Local variables (common in read_BC)
    character(len=20), allocatable :: bc_set_name(:), temp_bc_name(:)
    character(len=256) :: line
    character(len=20) :: field, tokens(100)
    integer :: bc_set_id
    integer :: i, j, l 

    ! Local variables used in "define"
    character(len=20), allocatable :: expand_bc_set_name(:)
    integer, allocatable :: expand_num_node_in_set(:), expand_bc_node_set(:,:)
    character(len=20) :: bc_set_type
    integer :: ini_node, fin_node, temp_node, token_count, set_node_count
    logical :: find_to

    ! Local variables used in fixed disp and point load
    integer :: num_stored_bc
    real(8) :: magnitude

    ! Temporary variables used in fixed disp and point load
    integer, allocatable :: temp_fixed_disp_vector(:,:)
    real(8), allocatable :: temp_fixed_disp_magn(:,:), temp_point_load_magn(:,:)

    ! ios variables
    integer :: ios_define, ios_token, ios_temp_node, ios_buscar



    !!! Go to the top of datfile 
    rewind (datfile)
    bc_set_id = 1

    !!! Skip to first "define"
    do
        read (datfile, "(A)", iostat = ios_define) line
        if (ios_define /= 0) then
            print *, "!!!!! ERROR : Define Block Not Found !!!!!"
            error stop
        end if

        if (line(1:12) == "define") then                    ! Find a block of "define"
            exit 
        end if
    end do

    ! Declare initial sizes of allocatable matricies
    allocate (bc_set_name(0))                               ! Initialize
    allocate (num_node_in_set(0))                           ! Initialize
    allocate (bc_node_set(0, nnode))                        ! initialize
    bc_set_id = 1                                           ! Initialize

    !!! Get "set" info of boundary conditions in define block
    do                                                      ! This is the largest loop to read "define" block.
        read (line(21:40), *) bc_set_type

        ! Confirm that the type if BC is "node".
        if (bc_set_type == "node") then
            ! No problem
        
        else 
            print *, "!!!!! ERROR : The type of boundary conditions are wrong !!!!!"
            error stop
        
        end if

        ! Expand bc_set_name by one element, by claude
        allocate (expand_bc_set_name(bc_set_id))
        expand_bc_set_name(1:bc_set_id-1) = bc_set_name(1:bc_set_id-1)
        call move_alloc (expand_bc_set_name, bc_set_name)

        ! Expand num_node_in_set by one element , by claude
        allocate (expand_num_node_in_set(bc_set_id))
        expand_num_node_in_set(1:bc_set_id-1) = num_node_in_set(1:bc_set_id-1)
        call move_alloc (expand_num_node_in_set, num_node_in_set) 

        ! Expand bc_node_set by one row (1st dimension), 2nd dimension fixed at nnode
        allocate(expand_bc_node_set(bc_set_id, nnode))
        expand_bc_node_set(1:bc_set_id-1, 1:nnode) = bc_node_set(1:bc_set_id-1, 1:nnode)
        call move_alloc(expand_bc_node_set, bc_node_set)
        
        read (line(61:80), *) bc_set_name(bc_set_id)        ! Store the name of set of BC
        set_node_count = 1

        ! Loop A1 !
        do 
            read (datfile, "(A)") line
            tokens = "@"
            read (line, *, iostat = ios_token) tokens                           ! Split line by space and store to array of tokens
            token_count = count(tokens /= "@")                                  ! Count the number of "tokens"
            find_to = any(tokens(1:token_count) == "to")                        ! Buscar "to" in "tokens"

            if (find_to) then                                                   ! If "tokens" contains "to", it means BC is applied to continuous nodes.
                read (tokens(1), *) ini_node
                read (tokens(3), *) fin_node
                num_node_in_set(bc_set_id) = fin_node - ini_node + 1            ! Store the number of nodes to apply this BC

                do l = 1, num_node_in_set(bc_set_id)
                    bc_node_set(bc_set_id, l) = ini_node + l - 1                ! Store node numbers to apply this BC 
                end do
            
            else 
                !Loop A2 !
                do i = 1, token_count
                    if (tokens(i) == "c") then                                  ! "c" means bc_node_set still continues to the next line 
                        exit

                    else
                        read (tokens(i), *, iostat = ios_temp_node) temp_node   ! Try to store token as Integer

                        if (ios_temp_node == 0) then                            ! If Integer, bc_node_set on this line still continues  
                            bc_node_set(bc_set_id, set_node_count) = temp_node  ! Store node number to apply BC
                            set_node_count = set_node_count + 1
                        
                        else                                                    ! This line should contain nothing except integers or "c"
                            print *, "!!!!! ERROR : Failed to read bc_set !!!!!"
                            error stop

                        end if
                    end if
                end do
                num_node_in_set(bc_set_id) = set_node_count - 1                 ! Store the number of nodes included in BC set

                if (tokens(token_count) == "c") then        ! If the last token is "c", go to the next line and repeat Loop A2
                    ! Read next line

                else                                        ! If the last token is integer, this "define" block end
                    exit
                
                end if
            end if
        end do

        bc_set_id = bc_set_id + 1
        read (datfile, "(A)") line

        if (line(1:12) == "define") then                    ! If the next line start with "define", repeat Loop A1
            ! Continuar do loop
        
        else                                                ! All "define" block has read. Completed.  
            exit

        end if
    end do

    !!! Initialize
    num_bc_set = bc_set_id -1
    bc_set_id = 1
    num_stored_bc = 0

    !!! Decide sizes of allocatable array
    allocate (fixed_disp_vector(num_bc_set, 2))
    allocate (fixed_disp_magn(num_bc_set, 2))
    allocate (point_load_magn(num_bc_set, 2))
    fixed_disp_vector = 0    
    fixed_disp_magn = 0.0d0
    point_load_magn = 0.0d0

    allocate (temp_fixed_disp_vector(num_bc_set, 2))
    allocate (temp_fixed_disp_magn(num_bc_set, 2))
    allocate (temp_point_load_magn(num_bc_set, 2))
    temp_fixed_disp_vector = 0
    temp_fixed_disp_magn = 0.0d0
    temp_point_load_magn = 0.0d0

    allocate (temp_bc_name(num_bc_set))

    ! Loop B1 !
    Buscar_All_BC : do
    ! 
        read (datfile, "(A)", iostat = ios_buscar) line
        if (ios_buscar /= 0) then
            print *, "!!!!! ERROR : The Boundary Conditions Not Found !!!!!"
            error stop
        end if

        !!! Read fixed disp
        if (line(1:12) == "fixed disp") then                ! Find a block of "fixed disp"
            read (datfile, *)                               ! Skip the next line. It contain no neccessary informations.
            
            ! Loop B2(fixed disp) !
            do 
                num_stored_bc = num_stored_bc + 1
                read (datfile, *)
                read (datfile, "(A)") line
                token_count = len_trim(line) / 20           ! Count the number of vector to fix by mesuring the length of the line

                do i = 1, token_count
                    read (line(20*i-19:20*i), "(A)") field
                    magnitude = expo2double(field)
                    temp_fixed_disp_magn(num_stored_bc, i) = magnitude          ! Store mafnitude of fixed disp for each vector
                end do

                read (datfile, *)

                read (datfile, "(A)") line
                tokens = "@"                                ! Initialize
                read (line, *, iostat = ios_token) tokens

                do j = 1, token_count
                    read (tokens(j), *) temp_fixed_disp_vector(num_stored_bc, j)                    ! Store all vector to fix displacement
                end do

                read (datfile, *)

                read (datfile, "(A)") line
                temp_bc_name(num_stored_bc) = line(1:20)    ! Get the name of this BC (This is used to map define blocks and BC blocks)

                do l = 1, num_bc_set                                                                ! "l" is just a index to search the same name.
                    if (temp_bc_name(num_stored_bc) == bc_set_name(l)) then                         ! When find the same name
                        do i = 1, 2
                            fixed_disp_magn(l, i) = temp_fixed_disp_magn(num_stored_bc, i)          ! Assign
                            fixed_disp_vector(l, i) = temp_fixed_disp_vector(num_stored_bc, i)      ! Assign
                        end do
                        exit
                    end if
                end do

                read (datfile, "(A)") line
                backspace (datfile)

                if (line(1:12) == "point load") then
                    exit
                end if

                if (num_stored_bc == num_bc_set) then
                    exit Buscar_All_BC
                
                else 
                    ! Continuar reading BC
                
                end if
            end do
        
        !!! Read point load
        else if (line(1:12) == "point load") then
            read (datfile, *)

            ! Loop B2(point load) !
            do
                num_stored_bc = num_stored_bc + 1
                read (datfile, *)
                read (datfile, "(A)") line
                token_count = len_trim(line) / 20           ! Count the number of vector to give load by mesuring the length of the line

                do i =1, token_count
                    read (line(20*i-19:20*i), "(A)") field
                    magnitude = expo2double(field)
                    temp_point_load_magn(num_stored_bc, i) = magnitude          ! Store magnitude of load for each vector
                end do
                
                read (datfile, *)
                read (datfile, *)

                read (datfile, "(A)") line
                temp_bc_name(num_stored_bc) = line(1:20)    ! Get the name of this BC

                do l = 1, num_bc_set                                                                ! "l" is just an index to find the same name
                    if (temp_bc_name(num_stored_bc) == bc_set_name(l)) then                         ! Find the same name
                        do i = 1, 2
                            point_load_magn(l, i) = temp_point_load_magn(num_stored_bc, i)          ! Assign
                        end do
                    exit
                    end if
                end do

                read (datfile, "(A)") line
                backspace (datfile)

                if (line(1:12) == "fixed disp") then
                    exit
                end if
                
                if (num_stored_bc == num_bc_set) then
                    exit Buscar_All_BC
                
                else
                    ! Continuar reading BC

                end if
            end do
        end if
    end do Buscar_All_BC

    !!! Close input datfile
    close (datfile)
              
end subroutine read_BC