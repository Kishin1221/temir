subroutine make_D(D, E, nu)

    use parameters
    implicit none

    ! Input arguments
    real(8), intent(in) :: E, nu

    ! Output argument
    real(8), intent(out) :: D(3,3)
    
    ! Local variable
    real(8) :: factor

    ! Initialize
    D = 0.0d0

    !!! This program only accept plane stress model now.
    ! Caluculate D matrix
    factor = E / (1-nu**2)

    D(1,1) = 1
    D(1,2) = nu
    D(2,1) = nu
    D(2,2) = 1
    D(3,3) = (1-nu) / 2

    D = D * factor

end subroutine make_D

subroutine make_B()

    use parameters
    implicit none

end subroutine make_B