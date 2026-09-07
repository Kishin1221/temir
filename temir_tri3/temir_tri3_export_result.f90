subroutine export_result(datfilename, nnode, nelem, coord, E, nu, t, disp, reaction_force, strain, stress)

    use parameters
    implicit none

    ! Input aarguments
    character(len=256), intent(in) :: datfilename
    integer, intent(in) :: nnode, nelem
    real(8), intent(in) :: coord(nnode_max, 2), E, nu, t
    real(8), intent(in) :: disp(2*nnode_max), reaction_force(2*nnode_max)
    real(8), intent(in) :: strain(nelem_max, 3), stress(nelem_max, 3)

    ! Local variables in export_result
    character(len=256) :: jobname
    integer :: outfile, pos, ielem, inode


    ! Extract job name from input file name 
    pos = index(datfilename, ".dat")
    jobname = datfilename(1 : pos-1)

    ! Make outfile
    open (newunit = outfile, file = trim(jobname)//"_temir.out", status = "replace", action = "write")

    ! Write model data
    write(outfile, *) "Job : ", trim(jobname)
    write(outfile, *)

    write(outfile, *) "----- Sizing -----"
    write(outfile, *) "    Number of Elements = ", nelem 
    write(outfile, *) "    Number of Nodes    = ", nnode
    write(outfile, *)

    write(outfile, *) "----- Material -----"
    write(outfile, *) "    Young's Modulus(E)  = ", E
    write(outfile, *) "    Poisson's Ratio(nu) = ", nu
    write(outfile, *)

    write(outfile, *) "----- Plane stress -----"
    write(outfile, *) "    Thickness (t) = ", t
    write(outfile, *)

    ! Write result
    write(outfile, *) "----- Result of Element Values -----"
    write(outfile, "(A10, 6A16)") "Element", "strain_xx", "strain_yy", "strain_xy", "stress_xx", "stress_yy", "stress_xy"
    
    do ielem = 1, nelem
        write(outfile, "(I10, 6E16.7)") ielem, strain(ielem, 1), strain(ielem, 2), strain(ielem, 3), stress(ielem, 1), stress(ielem, 2), stress(ielem, 3)
    end do

    write(outfile, *)

    write(outfile, *) "----- Result of Nodal Value -----"
    write(outfile, "(A10, 6A16)") "Node", "Coordinate_x", "Coordinate_y", "Displacement_x", "Displacement_y", "React_force_x", "React_force_y"

    do inode = 1, nnode
        write(outfile, "(I10, 6E16.7)") inode, coord(inode, 1), coord(inode, 2), disp(2*inode - 1), disp(2*inode), reaction_force(2*inode - 1), reaction_force(2*inode)
    end do

    write(outfile, *)
    write(outfile, *) "outfile end"

    close(outfile)

    print *, "Output file generated successfully."
    
end subroutine export_result