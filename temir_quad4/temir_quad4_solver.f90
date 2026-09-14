subroutine solver(reduced_stiffness, reduced_force, num_not_fixed, reduced_disp)

    use parameters
    implicit none

    ! Input arguments
    real(8), allocatable, intent(inout) :: reduced_stiffness(:,:), reduced_force(:)
    integer, intent(in) :: num_not_fixed

    ! Output arguments
    real(8), allocatable, intent(out) :: reduced_disp(:)

    ! Local variables
    integer, allocatable :: piv(:)
    integer :: i, j, ipiv
    integer :: best_row, temp_piv
    real(8) :: best_value, row_factor

    ! Declare matrices size
    allocate (piv(num_not_fixed))
    allocate (reduced_disp(num_not_fixed))

    do i = 1, num_not_fixed
        ! The equation No.i of reduced_stiffness is handled as No.piv(i) when calculate.  
        piv(i) = i                                          ! Initialize
    end do        

    !!! Forward Reduction
    do ipiv = 1, num_not_fixed
        !!! Pivoting
        ! Initialize
        best_row = ipiv
        best_value = reduced_stiffness(piv(ipiv), ipiv)

        do i = ipiv+1, num_not_fixed
            if(abs(reduced_stiffness(piv(i), ipiv)) > abs(best_value)) then
                best_row = i 
                best_value = reduced_stiffness(piv(i), ipiv)

            else
                ! Do nothing
            
            end if
        end do

        temp_piv = piv(ipiv)
        piv(ipiv) = piv(best_row)
        piv(best_row) = temp_piv

        ! Normalize ipiv-ipiv's entry of reduced_stiffness matrix
        do j = ipiv, num_not_fixed
            reduced_stiffness(piv(ipiv), j) = reduced_stiffness(piv(ipiv), j) / best_value
        end do

        ! Normalize ipiv's entry of reduced_force vector
        reduced_force(piv(ipiv)) = reduced_force(piv(ipiv)) / best_value 

        ! Subtracting
        do i = ipiv + 1, num_not_fixed
            row_factor = reduced_stiffness(piv(i), ipiv)

            do j = ipiv + 1, num_not_fixed
                reduced_stiffness(piv(i), j) = reduced_stiffness(piv(i), j) -  row_factor * reduced_stiffness(piv(ipiv), j) 
            end do

            reduced_force(piv(i)) = reduced_force(piv(i)) - row_factor * reduced_force(piv(ipiv))
        end do
    end do 
        
    !!! Back Substitution
    do ipiv = num_not_fixed, 1, -1
        reduced_disp(piv(ipiv)) = reduced_force(piv(ipiv))
        
        do j = ipiv + 1, num_not_fixed
            reduced_disp(piv(ipiv)) = reduced_disp(piv(ipiv)) - reduced_stiffness(piv(ipiv), j) * reduced_disp(j)
        end do    
    end do   

end subroutine solver