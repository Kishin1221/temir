subroutine strain_stress(nelem, connect, D, B, reduced_disp, disp, strain, stress)
    use parameters
    implicit none

    !!!!! Get displacement and D and B matrix, then calculate strain and stress.

    ! Declear variable get from main
    integer, intent(in) :: nelem, connect(nelem_max, 3)
    real(8), intent(in) :: D(3,3), B(nelem_max,3,6)
    real(8), intent(in) :: reduced_disp(2*nnode_max)

    ! Declear variables give to main
    real(8), intent(out) :: strain(nelem_max, 3)
    real(8), intent(out) :: stress(nelem_max, 3)

    ! Declear local variable in strain_stress
    integer :: ielem, i, j, k, l
    real(8) :: disp_local(6)

    !!! Vamos a Calcular
    do ielem = 1, nelem

        !Extract displacement of element No.e 
        do i = 1, 3
            disp_local(2*i-1) = disp(2*connect(ielem, i)-1)
            disp_local(2*i) = disp(2*connect(ielem, i))
        end do
        
        !Reset strain and stress
        do j = 1, 3
            strain(ielem, j) = 0.0d0
            stress(ielem, j) = 0.0d0
        end do

        !Calculate strain
        do j = 1, 3
            do k = 1, 6
                strain(ielem, j) = strain(ielem, j) + B(ielem, j, k) * disp_local(k)
            end do
        end do

        !Calculate stress
        do j = 1,3
            do l = 1, 3
                stress(ielem, j) = stress(ielem, j) + D(j, l) * strain(ielem, l)
            end do
        end do
    
    end do

end subroutine strain_stress