program test_make_D
    use parameters
    implicit none

    real(8) :: E, nu
    real(8) :: D(3,3)

    E = 100.d0
    nu = 1.0d0/3.0d0

    call make_D(D, E, nu)

    print *, D

end program test_make_D