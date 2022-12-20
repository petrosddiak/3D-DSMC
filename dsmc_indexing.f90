module indexing
        implicit none
    contains
        !NEEDS CARE, AS A CONSTANT (UNIFORM) SPACING IS ASSUMED BETWEEN CELLS
        subroutine pointers(position,num_particles,num_cells,x0,dx,y0,dy,z0,dz,indexes,num_cells_x,num_cells_y,num_cells_z,d_subcell_x,&
                            d_subcell_y,d_subcell_z,sub_cell_indexes)
            real*8 :: x0, dx, y0, dy, z0, dz, d_subcell_x, d_subcell_y, d_subcell_z
            integer :: num_particles, num_cells, i, j, k, l, m, bounce_back, num_cells_x, num_cells_y, num_cells_z, current_nop, summary, cell_counter
            real*8, dimension(num_particles,4) :: position
            integer*4, dimension(num_particles) :: local_particles, sub_cell_particles
            integer*4, dimension(num_particles,4) :: indexes, sub_cell_indexes
            character(len=1024) :: filename, name

            do i=1,num_particles
                indexes(i,2) = int((position(i,2)-x0)/dx)+1
                indexes(i,3) = int((position(i,3)-y0)/dy)+1
                indexes(i,4) = int((position(i,4)-z0)/dz)+1
            end do

            cell_counter = 1
            do i=1,num_cells_x
                do j=1,num_cells_y
                    do k=1,num_cells_z
                        current_nop = count((indexes(:,2).eq.i).and.(indexes(:,3).eq.j).and.(indexes(:,4).eq.k))
                        local_particles(1:current_nop) = pack(indexes(:,1),(indexes(:,2).eq.i).and.(indexes(:,3).eq.j).and.(indexes(:,4).eq.k))
                        do m=1,size(local_particles(1:current_nop))
                            sub_cell_indexes(local_particles(m),2) = int((position(local_particles(m),2)-(i-1)*dx)/d_subcell_x)+1
                            sub_cell_indexes(local_particles(m),3) = int((position(local_particles(m),3)-(j-1)*dy)/d_subcell_y)+1
                            sub_cell_indexes(local_particles(m),4) = int((position(local_particles(m),4)-(k-1)*dz)/d_subcell_z)+1
                        end do
                        cell_counter = cell_counter + 1
                    end do
                end do
            end do

            !write(filename,"(A25)") 'indexing.dat'
            !open(4,file=filename)
            !do i=1,num_particles
            !    write(4,*) indexes(i,:)
            !end do
            !close(4)

        end subroutine pointers
end module indexing