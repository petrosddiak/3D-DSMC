module check
        implicit none
    contains
        subroutine parameter_check(num_particles,num_cells_x,num_cells_y,num_cells_z,dx,dy,dz,dt,volume,velocity,xsize,ysize,zsize,position)
            integer*4 :: num_particles,num_cells_x,num_cells_y,num_cells_z,i
            real*8, dimension(num_particles,4) :: velocity, position
            real*8 :: dx,dy,dz,dt,volume,xsize,ysize,zsize

            print*, 'maximum CFL is ', max(maxval(velocity(:,2)/(dx/dt)),maxval(velocity(:,3)/(dy/dt)),maxval(velocity(:,4)/(dz/dt)))
            do i=1,num_particles
                if(((velocity(i,2)/(dx/dt)) .gt. 1) .or. ((velocity(i,3)/(dy/dt)) .gt. 1) .or. ((velocity(i,4)/(dz/dt)) .gt. 1)) then
                    print*, '-------CFL EXCEEDS 1 SOMEWHERE-------'
                    !call exit(0)
                end if

                if((position(i,2).lt.0) .or. (position(i,2).gt.xsize) .or. (position(i,3).lt.0) .or. (position(i,3).gt.ysize) .or. (position(i,4).lt.0) .or. (position(i,4).gt.zsize)) then
                    print*, 'Particle #', i
                    print*, '-------PARTICLE OUT OF BOUNDS-------'
                    !call exit(0)
                end if
            end do

        end subroutine parameter_check
end module check