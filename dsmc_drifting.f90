module drifting
        implicit none
    contains
        subroutine drift(xsize,ysize,x0,y0,z0,zsize,num_particles,position,velocity,dt,boltz,Twall,mol_mass,velocity_wall,flow_case,wall_bound_type,&
            x_wing_leading,y_wing,z_wing,wing_chord)
            real*8 :: nt, dt, xsize, ysize, zsize, velocity_wall, x_wall, y_wall, z_wall, x0, y0, z0, dt1, normal_comp, boltz, Twall, mol_mass, x_wing, x_wing_leading, y_wing, z_wing,&
                    wing_chord
            integer*4 :: num_particles, i, j, bounce_back, wall_bound_type, flow_case, component_indicator
            real*8, dimension(num_particles,4) :: position, velocity
            real*8, dimension(num_particles) :: dt_ac, x_old, y_old, z_old
            character(len=1024) :: filename, name

            x_old = position(:,2)
            y_old = position(:,3)
            z_old = position(:,4)

            position(:,2) = position(:,2) + dt*velocity(:,2)     ! advection in x direction
            position(:,3) = position(:,3) + dt*velocity(:,3)     ! advection in y direction
            position(:,4) = position(:,4) + dt*velocity(:,4)     ! advection in z direction

            do i=1,num_particles
                select case(flow_case)
                    case(2)
                        if((position(i,3)-ysize)*(y_old(i)-ysize).lt.0) then
                            component_indicator = 2
                            x_wall = (x_old(i)*(ysize-position(i,3))+position(i,2)*(y_old(i)-ysize))/(y_old(i)-position(i,3)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(ysize-position(i,3))+position(i,4)*(y_old(i)-ysize))/(y_old(i)-position(i,3)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(y_old(i)-ysize)/(y_old(i)-position(i,3))
                            normal_comp = -1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,3) = -velocity(i,3)                    ! reversing velocity
                                    velocity(i,2) = velocity(i,2) + velocity_wall     ! x-velocity component gets wall velocity update
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = ysize  + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if

                        if((position(i,3)-y0)*(y_old(i)-y0).lt.0) then
                            component_indicator = 2
                            x_wall = (x_old(i)*(y0-position(i,3))+position(i,2)*(y_old(i)-y0))/(y_old(i)-position(i,3)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(y0-position(i,3))+position(i,4)*(y_old(i)-y0))/(y_old(i)-position(i,3)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(y_old(i)-y0)/(y_old(i)-position(i,3))
                            normal_comp = 1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,3) = -velocity(i,3)                    ! reversing velocity
                                    velocity(i,2) = velocity(i,2) - velocity_wall     ! x-velocity component gets wall velocity update
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = y0     + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if

                    case(3)
                        if((position(i,3)-y_wing)*(y_old(i)-y_wing).lt.0) then
                            component_indicator = 2
                            x_wing = (x_old(i)*(y_wing-position(i,3))+position(i,2)*(y_old(i)-y_wing))/(y_old(i)-position(i,3)) !linear interp. to find x-component of location
                            z_wing = (z_old(i)*(y_wing-position(i,3))+position(i,4)*(y_old(i)-y_wing))/(y_old(i)-position(i,3)) !linear interp. to find z-component of location
                            if(x_wing.gt.x_wing_leading .and. x_wing.lt.(x_wing_leading+wing_chord)) then
                                dt1 = dt - dt*(y_old(i)-y_wing)/(y_old(i)-position(i,3))
                                normal_comp = -1.0
                                select case(wall_bound_type)
                                    case(1)
                                        velocity(i,3) = -velocity(i,3)                    ! reversing velocity
                                    case(2)
                                        call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                                end select
                                position(i,2) = x_wing + velocity(i,2)*dt1
                                position(i,3) = y_wing + velocity(i,3)*dt1
                                position(i,4) = z_wing + velocity(i,4)*dt1
                            end if
                        end if

                    case(4)
                        ! Y WALL
                        if((position(i,3)-ysize)*(y_old(i)-ysize).lt.0) then
                            component_indicator = 2
                            x_wall = (x_old(i)*(ysize-position(i,3))+position(i,2)*(y_old(i)-ysize))/(y_old(i)-position(i,3)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(ysize-position(i,3))+position(i,4)*(y_old(i)-ysize))/(y_old(i)-position(i,3)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(y_old(i)-ysize)/(y_old(i)-position(i,3))                                           !linear interp. to find dt
                            normal_comp = -1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,3) = -velocity(i,3)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = ysize  + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if

                        if((position(i,3)-y0)*(y_old(i)-y0).lt.0) then
                            component_indicator = 2
                            x_wall = (x_old(i)*(y0-position(i,3))+position(i,2)*(y_old(i)-y0))/(y_old(i)-position(i,3)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(y0-position(i,3))+position(i,4)*(y_old(i)-y0))/(y_old(i)-position(i,3)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(y_old(i)-y0)/(y_old(i)-position(i,3))
                            normal_comp = 1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,3) = -velocity(i,3)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = y0     + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if
                        ! X WALL
                        if((position(i,2)-xsize)*(x_old(i)-xsize).lt.0) then
                            component_indicator = 1
                            y_wall = (y_old(i)*(xsize-position(i,2))+position(i,3)*(x_old(i)-xsize))/(x_old(i)-position(i,2)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(xsize-position(i,2))+position(i,4)*(x_old(i)-xsize))/(x_old(i)-position(i,2)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(x_old(i)-xsize)/(x_old(i)-position(i,2))
                            normal_comp = -1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,2) = -velocity(i,2)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = xsize  + velocity(i,2)*dt1
                            position(i,3) = y_wall + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if

                        if((position(i,2)-x0)*(x_old(i)-x0).lt.0) then
                            component_indicator = 1
                            y_wall = (y_old(i)*(x0-position(i,2))+position(i,3)*(x_old(i)-x0))/(x_old(i)-position(i,2)) !linear interp. to find x-component of location
                            z_wall = (z_old(i)*(x0-position(i,2))+position(i,4)*(x_old(i)-x0))/(x_old(i)-position(i,2)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(x_old(i)-x0)/(x_old(i)-position(i,2))
                            normal_comp = 1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,2) = -velocity(i,2)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x0     + velocity(i,2)*dt1
                            position(i,3) = y_wall + velocity(i,3)*dt1
                            position(i,4) = z_wall + velocity(i,4)*dt1
                        end if
                        ! Z WALL
                        if((position(i,4)-zsize)*(z_old(i)-zsize).lt.0) then
                            component_indicator = 3
                            x_wall = (x_old(i)*(zsize-position(i,4))+position(i,2)*(z_old(i)-zsize))/(z_old(i)-position(i,4)) !linear interp. to find x-component of location
                            y_wall = (y_old(i)*(zsize-position(i,4))+position(i,3)*(z_old(i)-zsize))/(z_old(i)-position(i,4)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(z_old(i)-zsize)/(z_old(i)-position(i,4))
                            normal_comp = -1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,4) = -velocity(i,4)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = y_wall + velocity(i,3)*dt1
                            position(i,4) = zsize  + velocity(i,4)*dt1
                        end if

                        if((position(i,4)-z0)*(z_old(i)-z0).lt.0) then
                            component_indicator = 3
                            x_wall = (x_old(i)*(z0-position(i,4))+position(i,2)*(z_old(i)-z0))/(z_old(i)-position(i,4)) !linear interp. to find x-component of location
                            y_wall = (y_old(i)*(z0-position(i,4))+position(i,3)*(z_old(i)-z0))/(z_old(i)-position(i,4)) !linear interp. to find z-component of location
                            dt1 = dt - dt*(z_old(i)-z0)/(z_old(i)-position(i,4))
                            normal_comp = 1.0
                            select case(wall_bound_type)
                                case(1)
                                    velocity(i,4) = -velocity(i,4)                    ! reversing velocity
                                case(2)
                                    call diffuse_scatter(velocity(i,2:4),boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
                            end select
                            position(i,2) = x_wall + velocity(i,2)*dt1
                            position(i,3) = y_wall + velocity(i,3)*dt1
                            position(i,4) = z0     + velocity(i,4)*dt1
                        end if
                end select
            end do

        end subroutine drift

        subroutine diffuse_scatter(velocity,boltz,Twall,mol_mass,normal_comp,velocity_wall,flow_case,component_indicator)
            real*8 :: boltz, Twall, mol_mass, rand1, rand2, gaussian_random_number, rayleigh_random_number, velocity_wall, normal_comp
            real*8, dimension(3) :: velocity
            integer*4 :: flow_case, component_indicator
            real*8, parameter :: pi = 3.141592654
          
            call random_seed()

            if(component_indicator .eq. 1) then
                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2)) ! "1-random" instead of "random", to avoid getting 0 value at log
                velocity(2) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number

                call random_number(rand1)
                call random_number(rand2)
                rayleigh_random_number = sqrt(-2*log(1-rand1))
                velocity(1) = normal_comp*sqrt(boltz*Twall/mol_mass)*rayleigh_random_number

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2))
                velocity(3) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number

            else if(component_indicator .eq. 2) then
                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2)) ! "1-random" instead of "random", to avoid getting 0 value at log
                velocity(1) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number

                call random_number(rand1)
                call random_number(rand2)
                rayleigh_random_number = sqrt(-2*log(1-rand1))
                velocity(2) = normal_comp*sqrt(boltz*Twall/mol_mass)*rayleigh_random_number

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2))
                velocity(3) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number

            else if (component_indicator .eq. 3) then
                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2)) ! "1-random" instead of "random", to avoid getting 0 value at log
                velocity(1) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number

                call random_number(rand1)
                call random_number(rand2)
                rayleigh_random_number = sqrt(-2*log(1-rand1))
                velocity(3) = normal_comp*sqrt(boltz*Twall/mol_mass)*rayleigh_random_number

                call random_number(rand1)
                call random_number(rand2)
                gaussian_random_number = sqrt(-2*log(1-rand1))*cos(2*pi*(1-rand2))
                velocity(2) = sqrt(boltz*Twall/mol_mass)*gaussian_random_number
            
            end if

            if(flow_case.eq.2) then
                velocity(1) = velocity(1) + velocity_wall*(-normal_comp)       ! adding the wall velocity
            end if

        end subroutine diffuse_scatter

end module drifting