subroutine make_D(D, E, nu)
    implicit none

    !!!!! Get model info and make D-matrix.

    !! 平面ひずみか平面応力かをdatから読み取って両方計算できたらいいね．将来．!!
    
    !!! Plane Stress

    ! Declear variable
    real(8), intent(in) :: E, nu            ! E is elastic modulus, nu is poisson ratio.
    real(8) :: factor

    ! Declear matrices
    real(8), intent(out) :: D(3,3)

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