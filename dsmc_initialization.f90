module initialization
        implicit none
    contains
        subroutine init(xsize,ysize,zsize,dx,dy,dz,num_cells,num_particles,boltz,Tref,mol_mass,mfp,n_density,dref,ump, &
                        ump_ini,position,velocity,indexes,x0,y0,z0,num_cells_x,num_cells_y,num_cells_z,cell_centers,vel_x_free,&
                        vel_y_free,vel_z_free,T_ini,node_locs,velocity_wall,collision_model,flow_case,num_particles_left,&
                        num_particles_right,Vrel_max,d_subcell_x,d_subcell_y,d_subcell_z,sub_cells,sub_cell_indexes)
            integer*4 :: num_cells, num_cells_x, num_cells_y, num_cells_z, num_particles, i, j, k, l, m, n , o, seed, cell_counter, collision_model,flow_case,&
                        num_particles_left,num_particles_right,sub_cells, total_sub_cells, current_nop, current_subcell_nop
            real*8 :: xsize, ysize, zsize, dx, dy, dz, boltz, tref, mol_mass, n_density, dref, mfp, ump, ump_ini, t_ini, &
                         Kn, viscosity, density, number_density, pressure, u_wall, volume, x0, y0, z0, vel_x_free, vel_y_free, &
                         vel_z_free,velocity_wall,pressure_left,pressure_right,d_subcell_x,d_subcell_y,d_subcell_z
            real*8, dimension(num_particles,4) :: position, velocity
            integer*4, dimension(num_particles,4) :: indexes, sub_cell_indexes 
            real*8, dimension(num_particles,3) :: gaussian_random_number
            real*8, dimension(num_particles) :: rand1, rand2
            character(len=1024) :: filename, name
            real*8, parameter :: pi = 3.141592654
            real*8, dimension(num_cells) :: Vrel_max
            real*8, dimension(num_cells,6) :: cell_centers
            real*8, dimension((num_cells_x+1)*(num_cells_y+1)*(num_cells_z+1),3) :: node_locs
            integer*4, dimension(num_particles) :: local_particles
            integer*4, dimension(num_particles) :: sub_cell_particles

            do i=1,num_particles
                position(i,1) = i
                velocity(i,1) = i
                indexes(i,1) = i
                sub_cell_indexes(i,1) = i
            end do

            do i=1,num_cells
                Vrel_max(i) =  3*sqrt(boltz*Tref/mol_mass)
            end do

            if(flow_case .ne. 4) then

                call random_seed()
                call random_number(position(:,2:4))
                position(:,2) = position(:,2)*(xsize - x0) + x0
                position(:,3) = position(:,3)*(ysize - y0) + y0
                position(:,4) = position(:,4)*(zsize - z0) + z0

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,1) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:))) ! "1-random" instead of "random", to avoid getting 0 value at log

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,2) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:)))

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,3) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:)))

                velocity(:,2) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,1) + vel_x_free
                velocity(:,3) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,2) + vel_y_free
                velocity(:,4) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,3) + vel_z_free

                !velocity(:,2) = velocity(:,2) + velocity_wall*2*(position(:,3)/ysize-0.5)
                !velocity(:,3) = velocity(:,3) + 100

            else

                call random_seed()
                call random_number(position(:num_particles_left,2:4))
                position(:num_particles_left,2) = position(:num_particles_left,2)*(xsize*0.5 - x0) + x0
                position(:num_particles_left,3) = position(:num_particles_left,3)*(ysize - y0) + y0
                position(:num_particles_left,4) = position(:num_particles_left,4)*(zsize - z0) + z0

                call random_seed()
                call random_number(position(num_particles_left+1:,2:4))
                position(num_particles_left+1:,2) = position(num_particles_left+1:,2)*(xsize - xsize*0.5) + xsize*0.5
                position(num_particles_left+1:,3) = position(num_particles_left+1:,3)*(ysize - y0) + y0
                position(num_particles_left+1:,4) = position(num_particles_left+1:,4)*(zsize - z0) + z0

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,1) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:))) ! "1-random" instead of "random", to avoid getting 0 value at log

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,2) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:)))

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number(:,3) = sqrt(-2*log(1-rand1(:)))*cos(2*pi*(1-rand2(:)))

                velocity(:,2) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,1) + vel_x_free
                velocity(:,3) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,2) + vel_y_free
                velocity(:,4) = sqrt(2*boltz*tref/mol_mass)*gaussian_random_number(:,3) + vel_z_free
                
            end if

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
                        !print*, local_particles(1:current_nop)
                        do m=1,size(local_particles(1:current_nop))
                            sub_cell_indexes(local_particles(m),2) = int((position(local_particles(m),2)-(i-1)*dx)/d_subcell_x)+1
                            sub_cell_indexes(local_particles(m),3) = int((position(local_particles(m),3)-(j-1)*dy)/d_subcell_y)+1
                            sub_cell_indexes(local_particles(m),4) = int((position(local_particles(m),4)-(k-1)*dz)/d_subcell_z)+1
                            !print*, (i-1)*dx, d_subcell_x!(position(local_particles(m),2)-(i-1)*dx)/d_subcell_x
                        end do
                        cell_counter = cell_counter + 1
                    end do
                end do
            end do

            write(filename,"(A16)") 'cell_centers.dat'
            open(7000,file=filename)
            cell_counter = 1
            do i=1,num_cells_x
                do j=1,num_cells_y
                    do k=1,num_cells_z
                        cell_centers(cell_counter,1) = i  ! indices
                        cell_centers(cell_counter,2) = j  ! indices
                        cell_centers(cell_counter,3) = k  ! indices
                        cell_centers(cell_counter,4) = x0 + (i-0.5)*xsize/num_cells_x
                        cell_centers(cell_counter,5) = y0 + (j-0.5)*ysize/num_cells_y
                        cell_centers(cell_counter,6) = z0 + (k-0.5)*zsize/num_cells_z
                        write(7000,'(I3,I3,I3,F16.10,F16.10,F16.10)'), int(cell_centers(cell_counter,1)),int(cell_centers(cell_counter,2)),int(cell_centers(cell_counter,3)),cell_centers(cell_counter,4),cell_centers(cell_counter,5),cell_centers(cell_counter,6)
                        cell_counter = cell_counter + 1
                    end do
                end do
            end do
            close(7000)

            cell_counter = 1
            do k=1,num_cells_z+1
                do j=1,num_cells_y+1
                    do i=1,num_cells_x+1
                        node_locs(cell_counter,1) = x0 +(i-1)*xsize/num_cells_x
                        node_locs(cell_counter,2) = y0 +(j-1)*ysize/num_cells_y
                        node_locs(cell_counter,3) = z0 +(k-1)*zsize/num_cells_z
                        cell_counter = cell_counter + 1
                    end do
                end do
            end do

            write(filename,"(A25)") 'position_initial.dat'
            open(1,file=filename)
            do j=1,num_particles
                write(1,'(F14.6,F17.12,F17.12,F17.12)') position(j,:)
            end do
            close(1)

            write(filename,"(A20)") 'velocity_initial.dat'
            open(2,file=filename)
            do j=1,num_particles
                write(2,'(F14.6,F12.6,F12.6,F12.6)') velocity(j,:)
            end do
            close(2)

            write(filename,"(A25)") 'node_locations.dat'
            open(5,file=filename)
            do j=1,size(node_locs(:,1))
                write(5,*) node_locs(j,:)
            end do
            close(5)

            write(filename,"(A15)") 'indexing.dat'
            open(6,file=filename)
            do j=1,num_particles
                write(6,'(I4,I4,I4,I4)') indexes(j,:)
            end do
            close(6)

        end subroutine init
end module initialization
