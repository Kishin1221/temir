subroutine solver(num_not_fixed, reduced_stiffness, reduced_force, reduced_disp)
    
    use parameters
    implicit none

    !!!!! Solve simulatenous equation by using Gauss elimination method.

    ! Input arguments
    integer, intent(in) :: num_not_fixed
    real(8), intent(inout) :: reduced_stiffness(2*nnode_max, 2*nnode_max)
    real(8), intent(inout) :: reduced_force(2*nnode_max)
    
    ! Output arguments
    real(8), intent(out) :: reduced_disp(2*nnode_max)

    ! Local variables in solver
    integer :: i, j, k
    integer :: piv(2*nnode_max), best_row, temp_piv
    real(8) :: best_value, row_factor

    do i = 1, num_not_fixed
        ! The equation No.i of reduced_stiffness is handled as No.piv(i) when calculate.  
        piv(i) = i                                          ! Initialize
    end do

    !!! Forward Reduction
    do k = 1, num_not_fixed
        ! Pivoting
        best_row = k
        best_value = reduced_stiffness(piv(k), k)           ! Initialize

        do i = k + 1, num_not_fixed
            if (abs(reduced_stiffness(piv(i), k)) > best_value) then
                best_row = i                                                    ! Used for pivoting
                best_value = reduced_stiffness(piv(i), k)                       ! Used for normalization
                
            else 
                ! Do nothing

            end if
        end do

        if (abs(best_value) < 1.0d-8) then                   ! DEBUG
            print *, "k = ", k, " best_value = ", best_value  ! DEBUG
        end if

        temp_piv = piv(k)                                   ! Evacuation
        piv(k) = piv(best_row)                              ! The best row is moved to k'th row
        piv(best_row) = temp_piv                            ! The k'th row 
        
        ! Normalize k-k's entry of reduced_stiffness matrix
        do j = k, num_not_fixed
            reduced_stiffness(piv(k), j) = reduced_stiffness(piv(k), j) / best_value
        end do

        ! Normalize k's entry of reduced_force vector
        reduced_force(piv(k)) = reduced_force(piv(k)) / best_value 

        ! Subtracting
        do i = k + 1, num_not_fixed
            row_factor = reduced_stiffness(piv(i), k)

            do j = k + 1, num_not_fixed
                reduced_stiffness(piv(i), j) = reduced_stiffness(piv(i), j) -  row_factor * reduced_stiffness(piv(k), j) 
            end do

            reduced_force(piv(i)) = reduced_force(piv(i)) - row_factor * reduced_force(piv(k))
        end do
    end do 
        
    !!! Back Substitution
    do k = num_not_fixed, 1, -1
        reduced_disp(piv(k)) = reduced_force(piv(k))
        
        do j = k + 1, num_not_fixed
            reduced_disp(piv(k)) = reduced_disp(piv(k)) - reduced_stiffness(piv(k), j) * reduced_disp(j)
        end do    
    end do
        
end subroutine solver