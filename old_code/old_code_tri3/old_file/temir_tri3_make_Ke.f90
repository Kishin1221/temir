subroutine make_Ke(nelem, D, B, area, t, Ke)
    use parameters
    implicit none

    !!!!! Get D and B matrix from main and make element K stifness.

    ! Declear variables get from main
    integer, intent(in) :: nelem
    real(8), intent(in) :: D(3,3), B(nelem_max,3,6), area(nelem_max), t
    real(8), intent(out) :: Ke(nelem_max,6,6)

    ! Declear local variables in make_Ke
    integer :: ielem, i, j, k
    real(8) :: DB(nelem_max,3,6)

    !!! Calculate Ke matrix of ielem'th element
    do ielem = 1, nelem
        DB(ielem, :, :) = 0.0d0                             ! Initialize
        Ke(ielem, :, :) = 0.0d0                             ! Initialize

        do i = 1, 3
            do j = 1, 6
                do k = 1, 3
                    DB(ielem, i, j) = DB(ielem, i, j) + D(i, k) * B(ielem, k, j)            ! Calculate B*D of (B-transpose)DB first
                end do
            end do
        end do

        do i = 1, 6
            do j = 1, 6
                do k = 1, 3
                    Ke(ielem, i, j) = Ke(ielem, i, j) + (B(ielem, k, i) * DB(ielem, k, j))  ! Calculate (B-transpose)DB
                end do
                Ke(ielem, i, j) = Ke(ielem, i, j) * t * area(ielem)                         ! Maltiple thickness and area
            end do
        end do
    end do

end subroutine make_Ke