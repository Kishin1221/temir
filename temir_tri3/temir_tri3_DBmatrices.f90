subroutine make_D(D, E, nu)
    implicit none

    !!!!! Get model info and make D-matrix.

    !! 平面ひずみか平面応力かをdatから読み取って両方計算できたらいいね．将来．!!
    
    !!! Plane Stress

    ! Input arguments
    real(8), intent(in) :: E, nu            ! E is elastic modulus, nu is poisson ratio.

    ! Output arguments
    real(8), intent(out) :: D(3,3)

    ! Local variables in make_D
    real(8) :: factor

    ! Calculate values of entry
    factor = E/(1-nu**2)                    ! common factor               

    D = 0.0d0

    D(1,1) = 1                              ! Calcular each entry
    D(1,2) = nu
    D(2,1) = nu
    D(2,2) = 1
    D(3,3) = (1-nu)/2

    D = factor * D

end subroutine make_D

subroutine make_B(nelem, connect, coord, B, area)
    use parameters
    implicit none

    !!!!! Get model info and make B-matrix

    ! Input arfuments
    integer, intent(in) :: nelem, connect(nelem_max, 3)
    real(8), intent(in) :: coord(nnode_max, 2)

    ! Output arguments
    real(8), intent(out) :: B(nelem_max, 3, 6)
    real(8), intent(out) :: area(nelem_max)

    ! Local variables in make_B 
    integer :: ielem, i, m, n
    real(8) :: bi(nelem_max,3), ci(nelem_max,3)

    !!! Calculate entries of B-matrix
    do ielem = 1, nelem
        
        !Reset to 0
        do m = 1, 3
            do n = 1, 6
                B (ielem, m, n) = 0.0d0                                 ! Initialize
            end do
        end do

        area(ielem) = 0.0d0                                             ! Initialize

        bi(ielem,1) = coord(connect(ielem, 2),2) - coord(connect(ielem, 3),2)
        bi(ielem,2) = coord(connect(ielem, 3),2) - coord(connect(ielem, 1),2)
        bi(ielem,3) = coord(connect(ielem, 1),2) - coord(connect(ielem, 2),2)
        ci(ielem,1) = coord(connect(ielem, 3),1) - coord(connect(ielem, 2),1)
        ci(ielem,2) = coord(connect(ielem, 1),1) - coord(connect(ielem, 3),1)
        ci(ielem,3) = coord(connect(ielem, 2),1) - coord(connect(ielem, 1),1)  

        do i = 1,3
            area(ielem) = area(ielem) + (coord(connect(ielem, i),1) * bi(ielem, i))/2     ! Calcular area of element
        end do
        
        !Caluculate entry of B-matrix
        B(ielem, 1, 1) = bi(ielem, 1)/(2*area(ielem))
        B(ielem, 1, 3) = bi(ielem, 2)/(2*area(ielem))
        B(ielem, 1, 5) = bi(ielem, 3)/(2*area(ielem))
        B(ielem, 2, 2) = ci(ielem, 1)/(2*area(ielem))
        B(ielem, 2, 4) = ci(ielem, 2)/(2*area(ielem))
        B(ielem, 2, 6) = ci(ielem, 3)/(2*area(ielem))
        B(ielem, 3, 1) = B(ielem, 2, 2) 
        B(ielem, 3, 2) = B(ielem, 1, 1) 
        B(ielem, 3, 3) = B(ielem, 2, 4) 
        B(ielem, 3, 4) = B(ielem, 1, 3) 
        B(ielem, 3, 5) = B(ielem, 2, 6) 
        B(ielem, 3, 6) = B(ielem, 1, 5) 

    end do
    
end subroutine make_B