module boundary
        implicit none
    contains
        subroutine boundaries(xsize,ysize,zsize,num_particles,position,velocity,dt,velocity_wall,flow_case,x0,y0,z0)
            real*8 :: dt, xsize, ysize, zsize, velocity_wall, x0, y0, z0
            integer :: num_particles, i, j, bounce_back, flow_case
            real*8, dimension(num_particles,4) :: position, velocity
            real*8, dimension(num_particles) :: dt_ac
            character(len=1024) :: filename, name


            select case(flow_case)
                case(1)
                    ! Y DIRECTION
                    where(position(:,3).gt.ysize)                         ! checking whether a wall boundary was hit
                        dt_ac = (position(:,3)-ysize)/velocity(:,3)       ! time after collision
                        velocity(:,3) = -velocity(:,3)                    ! reversing velocity
                        velocity(:,2) = velocity(:,2) + velocity_wall     ! x-velocity component gets wall velocity update
                        position(:,3) = ysize + dt_ac*velocity(:,3)       ! updating molecule position
                    elsewhere(position(:,3).lt.0)
                        dt_ac = (position(:,3))/velocity(:,3)             ! time after collision
                        velocity(:,3) = -velocity(:,3)                    ! reversing velocity
                        position(:,3) = dt_ac*velocity(:,3)               ! updating molecule position
                    ! Z DIRECTION
                    elsewhere(position(:,4).gt.zsize)
                        dt_ac = (position(:,4)-zsize)/velocity(:,4)       ! time after collision
                        velocity(:,4) = -velocity(:,4)                    ! reversing velocity
                        position(:,4) = zsize + dt_ac*velocity(:,4)       ! updating molecule position
                    elsewhere(position(:,4).lt.0)
                        dt_ac = (position(:,4))/velocity(:,4)             ! time after collision
                        velocity(:,4) = -velocity(:,4)                    ! reversing velocity
                        position(:,4) = dt_ac*velocity(:,4)               ! updating molecule position
                    ! X DIRECTION
                    elsewhere(position(:,2).gt.xsize)
                        dt_ac = (position(:,2)-xsize)/velocity(:,2)       ! time after collision
                        velocity(:,2) = -velocity(:,2)                    ! reversing velocity
                        position(:,2) = xsize + dt_ac*velocity(:,2)       ! updating molecule position
                    elsewhere(position(:,2).lt.0)
                        dt_ac = (position(:,2))/velocity(:,2)             ! time after collision
                        velocity(:,2) = -velocity(:,2)                    ! reversing velocity
                        position(:,2) = dt_ac*velocity(:,2)               ! updating molecule position  
                    end where

                    if(any(position(:,3).gt.ysize) .or. any(position(:,3).lt.0)) then
                        print*, 'Particles over top wall: ', count(position(:,3).gt.ysize)
                        print*, 'Particles below bottom wall: ', count(position(:,3).lt.0)
                        print*, '---------------Some particles exceed wall limits---------------'
                        call exit(0)
                    end if

                case(2) ! COUETTE FLOW
                    ! X DIRECTION
                    where(position(:,2).gt.xsize)
                        position(:,2) = position(:,2) - xsize             ! updating molecule position
                    elsewhere(position(:,2).lt.0)
                        position(:,2) = position(:,2) + xsize
                    !! Y DIRECTION
                    !elsewhere(position(:,3).gt.ysize)
                    !    dt_ac = (position(:,3)-ysize)/velocity(:,3)       ! time after collision
                    !    velocity(:,3) = -velocity(:,3)                    ! reversing velocity
                    !    position(:,3) = ysize + dt_ac*velocity(:,3)       ! updating molecule position
                    !    velocity(:,2) = velocity(:,2) + velocity_wall
                    !elsewhere(position(:,3).lt.0)
                    !    dt_ac = (position(:,3))/velocity(:,3)             ! time after collision
                    !    velocity(:,3) = -velocity(:,3)                    ! reversing velocity
                    !    position(:,3) = y0 + dt_ac*velocity(:,3)          ! updating molecule position
                    !    velocity(:,2) = velocity(:,2) - velocity_wall
                    ! Z DIRECTION
                    elsewhere(position(:,4).gt.zsize)
                        position(:,4) = position(:,4) - zsize             ! updating molecule position
                    elsewhere(position(:,4).lt.0)
                        position(:,4) = position(:,4) + zsize
                    end where   

                case(3) ! THIN WING
                    ! Z DIRECTION
                    where(position(:,4).gt.zsize)
                        dt_ac = (position(:,4)-zsize)/velocity(:,4)       ! time after collision
                        velocity(:,4) = -velocity(:,4)                    ! reversing velocity
                        position(:,4) = zsize + dt_ac*velocity(:,4)       ! updating molecule position
                    elsewhere(position(:,4).lt.0)
                        dt_ac = (position(:,4))/velocity(:,4)             ! time after collision
                        velocity(:,4) = -velocity(:,4)                    ! reversing velocity
                        position(:,4) = dt_ac*velocity(:,4)               ! updating molecule position
                    ! X DIRECTION
                    elsewhere(position(:,2).gt.xsize)
                        position(:,2) = position(:,2) - xsize             ! updating molecule position
                    elsewhere(position(:,2).lt.0)
                        position(:,2) = position(:,2) + xsize             ! updating molecule position  
                    end where

                case(4) ! SHOCK TUBE
                    !! X DIRECTION
                    !where(position(:,2).lt.x0)
                    !    dt_ac = (position(:,2)-x0)/velocity(:,2)       ! time after collision
                    !    velocity(:,2) = -velocity(:,2)                 ! reversing velocity
                    !    position(:,2) = x0 + dt_ac*velocity(:,2)       ! updating molecule position
                    !elsewhere(position(:,2).gt.xsize)
                    !    dt_ac = (position(:,2)-xsize)/velocity(:,2)       ! time after collision
                    !    velocity(:,2) = -velocity(:,2)                    ! reversing velocity
                    !    position(:,2) = xsize + dt_ac*velocity(:,2)       ! updating molecule position
                    !! Y DIRECTION
                    !elsewhere(position(:,3).lt.y0)
                    !    dt_ac = (position(:,3)-y0)/velocity(:,3)       ! time after collision
                    !    velocity(:,3) = -velocity(:,3)                 ! reversing velocity
                    !    position(:,3) = y0 + dt_ac*velocity(:,3)       ! updating molecule position
                    !elsewhere(position(:,3).gt.ysize)
                    !    dt_ac = (position(:,3)-ysize)/velocity(:,3)       ! time after collision
                    !    velocity(:,3) = -velocity(:,3)                    ! reversing velocity
                    !    position(:,3) = ysize + dt_ac*velocity(:,3)       ! updating molecule position
                    !! Z DIRECTION
                    !elsewhere(position(:,4).lt.z0)
                    !    dt_ac = (position(:,4)-z0)/velocity(:,4)       ! time after collision
                    !    velocity(:,4) = -velocity(:,4)                 ! reversing velocity
                    !    position(:,4) = z0 + dt_ac*velocity(:,4)       ! updating molecule position
                    !elsewhere(position(:,4).gt.zsize)
                    !    dt_ac = (position(:,4)-zsize)/velocity(:,4)       ! time after collision
                    !    velocity(:,4) = -velocity(:,4)                    ! reversing velocity
                    !    position(:,4) = zsize + dt_ac*velocity(:,4)       ! updating molecule position
                    !end where
            end select

        end subroutine boundaries
end module boundary