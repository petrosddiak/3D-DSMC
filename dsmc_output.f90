module results_output
        implicit none
    contains
        subroutine output(updated_sampled_values,num_cells,num_cells_x,num_cells_y,num_cells_z,node_locs,number_of_sampled_variables)
            integer :: num_cells, num_cells_x, num_cells_y, num_cells_z, i, number_of_sampled_variables
            real*8, dimension((num_cells_x+1)*(num_cells_y+1)*(num_cells_z+1),3) :: node_locs
            real*8, dimension(num_cells,number_of_sampled_variables) :: updated_sampled_values
            character(len=1024) :: filename, name

            !write(filename,"(A20)") 'results_final.dat'
            write(filename,"(A40)") 'results_couette_10000.dat'
            open(unit=111, file=filename)

            write (111,'(A100)') 'Variables= "X", "Y", "Z", "Vel", "U", "V", "W", "N", "P", "T", "rho"'
            write (111,*) 'ZONE I=', num_cells_x+1,', J= ',num_cells_y+1,', K= ',num_cells_z+1, ', DATAPACKING=BLOCK, VARLOCATION=([4-11]=CELLCENTERED)'
            !write (111,*) 'ZONE DATAPACKING=POINT, I=', num_cells_x,', J= ',num_cells_y,' K= ',num_cells_z
            !do i=1,num_cells_x
            !    do j=1,num_cells_y
            !        do k=1,num_cells_z
            !            write(111,'(F12.8,F12.8,F12.8,F13.6,F13.6,F13.6,F13.6)') cell_centers(counter,4),cell_centers(counter,5),cell_centers(counter,6),&
            !            updated_sampled_values(counter,1),updated_sampled_values(counter,2),updated_sampled_values(counter,3),updated_sampled_values(counter,4)
            !            counter = counter + 1
            !        end do
            !    end do
            !end do

            do i=1,size(node_locs(:,1))
                write(111,'(F17.12)') node_locs(i,1)
            end do

            do i=1,size(node_locs(:,2))
                write(111,'(F17.12)') node_locs(i,2)
            end do

            do i=1,size(node_locs(:,3))
                write(111,'(F17.12)') node_locs(i,3)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,1)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,2)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,3)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,4)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,*) updated_sampled_values(i,5)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,*) updated_sampled_values(i,6)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,7)
            end do

            do i=1,num_cells_x*num_cells_y*num_cells_z
                write(111,'(F15.8)') updated_sampled_values(i,8)
            end do            

            close(111)

        end subroutine output
end module results_output
