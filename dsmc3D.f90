program dsmc3d
    !use omp_lib
    use check
    use initialization
    use drifting
    use boundary
    use indexing
    use collision
    use sampling
    use results_output

    implicit none

    external convergence
    integer*4 :: max_par_percell, num_particles, max_num_par, num_cells_x, num_cells_y, num_cells_z, num_cells, step,&
                current_nop, i, j, k, l, m, n, ii, jj, kk, timestep, num_collisions, i_prop, j_prop, cell_num_collisions,&
                counter, num_timesteps, cell_counter, samples, sampling_frequency, flow_case, max_num_particles,&
                sampling_start, collision_model, output_option, wall_bound_type, num_particles_left, num_particles_right,trnt_t_step,&
                total_trnts, transient_samples, number_of_sampled_variables, residual_counter, sub_cells, total_sub_cells, sub_cell_nop,&
                j_prop_loc
    real*8, allocatable :: position(:,:), velocity(:,:), local_variables(:,:), Vrel_max(:)
    integer*4, allocatable :: indexes(:,:), local_particles(:), sub_cell_indexes(:,:), same_sub_cell_particles(:)
    integer*4, dimension(3) :: i_prop_indexes
    real*8 :: xsize, ysize, zsize, t, dt, vmax_cell, xmin, xmax, dx, dy, dz, lambda, volume, n_density, boltz, mol_mass,&
        density, tref, mfp, ump, Kn, ump_ini, t_ini, pressure, viscosity, number_density,&
        Rf, Sigma, Vrel, S_max, cell_volume, first_molecule, second_molecule, velocityx_mean, velocityy_mean,&
        velocityz_mean, velocity_wall, eff_num, coeff, timer_start, timer_finish, num_cand_collisions,vel_x_free,vel_y_free,&
        vel_z_free, alpha, t1, t2, diameter, dref, Vel_rel_ref, Tfree, omega, SigmaRef, Smax, Twall, x0, y0, z0, x_wing_leading,&
        y_wing, z_wing, wing_chord, weight, pressure_left, pressure_right, T_left, T_right, n_density_left, n_density_right, real_time,&
        d_subcell_x,d_subcell_y,d_subcell_z
    real*8, dimension(3) :: vel_rel_vec
    real*8, parameter :: pi = 3.141592654
    integer*4, dimension(8) :: nops
    logical :: conditionx
    character(len=8) :: actual_time
    character(len=1024) :: filename, name
    real*8, allocatable :: sampled_values(:,:), cell_centers(:,:), updated_sampled_values(:,:), node_locs(:,:), old_sampled_values(:,:), average_residual(:)

    !call random_seed()
    T_ini = 100.0d0; Tfree = 300.0d0; Twall = 300; T_left = 273; T_right = 273; Tref = 273
    Kn = 0.1
    pressure = 1000d0; boltz = 1.3806e-23; mol_mass = 6.63e-26    ! Mass of argon atom (kg)
    viscosity = 2.117e-05
    !density = 2.685e25     ! Number density of argon at STP (m^-3)
    number_density = pressure/(boltz*T_ini)
    ump_ini = sqrt(2.d0*boltz*Tref/mol_mass)
    vel_x_free = 0; vel_y_free = 0; vel_z_free = 0
    mfp = viscosity*16.0d0/(5.0d0*sqrt(pi)*number_density*ump_ini) !!hard sphere mean path, for others page 139 of rarefied gas dynamics
    xsize = 6.33e-7; ysize = 6.33e-7; zsize = 1  !minimum size is 10e-7, this can be changed in the initialization and output write to file commands
    y_wing = 0.5; x_wing_leading = 0.4; wing_chord = 0.2
    dt = 4.33e-12; velocity_wall = 50
    x0 = 0; y0 = 0; z0 = 0
    volume = xsize*ysize*zsize
    num_particles = 90000
    max_num_particles = 1000000
    num_cells_x = 50; num_cells_y = 50; num_cells_z = 1
    sub_cells = 8
    num_cells = num_cells_x*num_cells_y*num_cells_z
    total_sub_cells = sub_cells*num_cells
    dx = xsize/float(num_cells_x)
    dy = ysize/float(num_cells_y)
    dz = zsize/float(num_cells_z)
    d_subcell_x = dx/(sub_cells**(1.0/3.0))
    d_subcell_y = dy/(sub_cells**(1.0/3.0))
    d_subcell_z = dz/(sub_cells**(1.0/3.0))
    !weight = dx*dy*dz*number_density/(num_particles*volume)
    weight = 0
    eff_num = number_density*volume/num_particles !Number of actual molecules, that 1 particle represents
    omega = 0.81 ! viscosity coefficient index

    collision_model = 1
    select case(collision_model)
        case(1) ! HS
            dref = 3.66e-10
            diameter = dref
            SigmaRef = pi*dref**2
            Sigma = SigmaRef
            viscosity = 5/16*sqrt(mol_mass*boltz*Tref/pi)/(diameter**2)

        case(2) ! VHS
            dref = 4.17e-10
            SigmaRef = pi*dref**2
            Vel_rel_ref = sqrt(4*boltz*Tref/mol_mass)/(gamma(2.5-omega)**(1/(2*omega-1)))

        case(3) ! VSS
            dref = 4.11e-10
            SigmaRef = pi*dref**2
            Vel_rel_ref = sqrt(4*boltz*Tref/mol_mass)/(gamma(2.5-omega)**(1/(2*omega-1)))      
    end select

    cell_volume = volume/num_cells
    num_timesteps = 10000; total_trnts = 1; sampling_frequency = 10; sampling_start = 7000; output_option = 0; flow_case = 2; wall_bound_type = 2; number_of_sampled_variables = 8 
    !flow_case: 2 is couette, 3 is thin wing, 4 is shock tube --- wall_bound_type: 1 is specular, 2 is diffuse reflection
    alpha = 3.1

    if(flow_case .eq. 4) then
        eff_num = 1.2e10
        pressure_right = 44.2
        pressure_left = 11.9*pressure_right
        n_density_left = pressure_left/(boltz*T_left)
        n_density_right = pressure_right/(boltz*T_right)
        !n_density_left = 1e-4/mol_mass
        !n_density_right = 0.125e-4/mol_mass
        num_particles_left = n_density_left*xsize*ysize*zsize*0.5/eff_num
        num_particles_right = n_density_right*xsize*ysize*zsize*0.5/eff_num
        num_particles = num_particles_left + num_particles_right
        print*, 'NOP: ', num_particles_left, num_particles_right, n_density_left*boltz*T_left, n_density_right*boltz*T_right
    else if(flow_case .eq. 2) then
        eff_num = 1e11
        pressure = 100000
        number_density = pressure/(boltz*Tref)
        !num_particles = number_density*xsize*ysize*zsize/eff_num
        print*, 'NOP: ', num_particles
    end if
    
    coeff = 0.5*eff_num*pi*dref**2*dt/(volume/num_cells)
    if(coeff.lt.1e-4) then
        !call exit(0)
    end if

    print*, 'ndensityleft and right is ', n_density_left,n_density_right
    print*, 'coeff ', coeff, 'avg nop per cell ', num_particles/num_cells

    !print*, 0.5*pi*dref**2*dt/(volume/num_cells)
    !print*, Vrel_max, S_max, 'Eff num ', eff_num, 'coeff is ', coeff, dt
    print*, 'number density ', number_density
    !print*, 'umpini is ', ump_ini, 'dx/dt', dx/dt

    call cpu_time(timer_start)

    allocate(position(num_particles,4),velocity(num_particles,4),indexes(num_particles,4),sampled_values(num_cells,number_of_sampled_variables),cell_centers(num_cells,6),&
            updated_sampled_values(num_cells,number_of_sampled_variables),node_locs((num_cells_x+1)*(num_cells_y+1)*(num_cells_z+1),3),Vrel_max(num_cells),&
            old_sampled_values(num_cells,number_of_sampled_variables),average_residual(total_trnts*((num_timesteps-sampling_start)/sampling_frequency+1)),&
            sub_cell_indexes(num_particles,4),same_sub_cell_particles(num_particles))

    allocate(local_particles(num_particles))
    allocate(local_variables(num_particles,4))

    transient_samples = 0
    residual_counter = 0
    do trnt_t_step=1,total_trnts
        call init(xsize,ysize,zsize,dx,dy,dz,num_cells,num_particles,boltz,tref,mol_mass,mfp,n_density,dref,ump, &
            ump_ini,position,velocity,indexes,x0,y0,z0,num_cells_x,num_cells_y,num_cells_z,cell_centers,vel_x_free,vel_y_free,vel_z_free,T_ini,&
            node_locs,velocity_wall,collision_model,flow_case,num_particles_left,num_particles_right,Vrel_max,d_subcell_x,d_subcell_y,&
            d_subcell_z,sub_cells,sub_cell_indexes)

        !call parameter_check(num_particles,num_cells_x,num_cells_y,num_cells_z,dx,dy,dz,dt,volume,velocity,xsize,ysize,zsize,position)

        old_sampled_values(:,:) = 1
        real_time = 0
        transient_samples = transient_samples + 1
        samples = 0;

        do timestep=1,num_timesteps
            real_time = real_time + dt
            !write(*,1000), 'Timestep ',timestep,' of ',num_timesteps
            !1000 format (A8,1X,I4,1X,A4,1X,I4)

            call drift(xsize,ysize,x0,y0,z0,zsize,num_particles,position,velocity,dt,boltz,Twall,mol_mass,velocity_wall,flow_case,wall_bound_type,x_wing_leading,y_wing,z_wing,wing_chord)
            call boundaries(xsize,ysize,zsize,num_particles,position,velocity,dt,velocity_wall,flow_case,x0,y0,z0)
            call pointers(position,num_particles,num_cells,x0,dx,y0,dy,z0,dz,indexes,num_cells_x,num_cells_y,num_cells_z,d_subcell_x,d_subcell_y,d_subcell_z,sub_cell_indexes) 
            counter = 1

            do k=1,num_cells_z
                do j=1,num_cells_y
                    do i=1,num_cells_x
                        current_nop = count((indexes(:,2).eq.i).and.(indexes(:,3).eq.j).and.(indexes(:,4).eq.k)) ! counting cell's # of particles
                        if(current_nop.lt.2) then
                            cycle
                        end if

                        local_particles = pack(indexes(:,1),(indexes(:,2).eq.i).and.(indexes(:,3).eq.j).and.(indexes(:,4).eq.k)) !finding particle numbers

                        do m=1,4
                            local_variables(1:current_nop,m) = velocity(local_particles,m)
                        end do

                        Smax = SigmaRef*Vrel_max(counter)

                        num_cand_collisions = coeff*current_nop**2*Vrel_max(counter)
                        cell_num_collisions = 0
                        !print*, '**********************************************************************'
                        !print*, num_cand_collisions, current_nop
                        do l=1,int(num_cand_collisions)
                            !call random_seed()
                            call random_number(first_molecule)
                            i_prop = int(first_molecule*size(local_particles))+1
                            i_prop_indexes = (sub_cell_indexes(local_particles(i_prop),2:4))
                            sub_cell_nop = count((sub_cell_indexes(local_particles,2).eq.i_prop_indexes(1)).and.(sub_cell_indexes(local_particles,3)&
                                                .eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3)))
                            same_sub_cell_particles = pack(sub_cell_indexes(local_particles,1),((sub_cell_indexes(local_particles,2).eq.i_prop_indexes(1))&
                                                            .and.(sub_cell_indexes(local_particles,3).eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3))))
                            !print*, size(local_particles), sub_cell_nop
                            call random_number(second_molecule)
                            call random_number(Rf)
                            j_prop = int(second_molecule*size(same_sub_cell_particles))+1

                            !if(sub_cell_nop.le.1) then
                            !    print*, 'ton ipiame'
                            !end if
                            
                            conditionx = i_prop_indexes(1).eq.1 .or. conditionx.eq.sub_cells**(1.0/3.0)
                            if(sub_cell_nop.lt.2) then
                                call random_number(second_molecule)
                                if(conditionx) then
                                    if(i_prop_indexes(1).eq.1) then
                                        sub_cell_nop = count((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)+1)).and.(sub_cell_indexes(local_particles,3)&
                                                            .eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3)))
                                        if(sub_cell_nop.eq.0) exit
                                        same_sub_cell_particles = pack(sub_cell_indexes(local_particles,1),((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)+1))&
                                                                .and.(sub_cell_indexes(local_particles,3).eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3))))
                                        j_prop = int(second_molecule*size(same_sub_cell_particles))+1
                                        j_prop_loc = findloc(local_particles,same_sub_cell_particles(j_prop),1)
                                    else
                                        sub_cell_nop = count((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)-1)).and.(sub_cell_indexes(local_particles,3)&
                                                            .eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3)))
                                        if(sub_cell_nop.eq.0) exit
                                        same_sub_cell_particles = pack(sub_cell_indexes(local_particles,1),((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)-1))&
                                                                .and.(sub_cell_indexes(local_particles,3).eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3))))
                                        j_prop = int(second_molecule*size(same_sub_cell_particles))+1
                                        j_prop_loc = findloc(local_particles,same_sub_cell_particles(j_prop),1)
                                    end if

                                else
                                    sub_cell_nop = count((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)+1)).and.(sub_cell_indexes(local_particles,3)&
                                                            .eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3)))
                                    if(sub_cell_nop.eq.0) exit
                                    same_sub_cell_particles = pack(sub_cell_indexes(local_particles,1),((sub_cell_indexes(local_particles,2).eq.(i_prop_indexes(1)+1))&
                                                                .and.(sub_cell_indexes(local_particles,3).eq.i_prop_indexes(2)).and.(sub_cell_indexes(local_particles,4).eq.i_prop_indexes(3))))
                                    j_prop = int(second_molecule*size(same_sub_cell_particles))+1
                                    j_prop_loc = findloc(local_particles,same_sub_cell_particles(j_prop),1)
                                end if
                            end if

                            do while(local_particles(i_prop).eq.same_sub_cell_particles(j_prop))
                                call random_number(second_molecule)
                                j_prop = int(second_molecule*size(same_sub_cell_particles))+1
                            end do

                            j_prop_loc = findloc(local_particles,same_sub_cell_particles(j_prop),1)
                            !print*, local_particles
                            !print*, local_particles(i_prop), local_particles(j_prop_loc)
                            !print*, i_prop, j_prop_loc
                            !print*, sub_cell_indexes(:,1)
                            !print*, sub_cell_indexes(local_particles(i_prop),:)
                            !print*, sub_cell_indexes(local_particles(j_prop_loc),:)
                            !print*, '-------------------------------------------------------'

                            vel_rel_vec(1) = local_variables(i_prop,2) - local_variables(j_prop_loc,2)
                            vel_rel_vec(2) = local_variables(i_prop,3) - local_variables(j_prop_loc,3)
                            vel_rel_vec(3) = local_variables(i_prop,4) - local_variables(j_prop_loc,4)                        

                            Vrel = sqrt(vel_rel_vec(1)**2+vel_rel_vec(2)**2+vel_rel_vec(3)**2)
                            if(Vrel.gt.Vrel_max(counter)) then
                                Vrel_max(counter) = Vrel
                                Smax = SigmaRef*Vrel_max(counter)
                            end if
                           
                            if(collision_model.eq.2 .or. collision_model.eq.3) then
                                Sigma = SigmaRef*(Vel_rel_ref/Vrel)**(2*omega-1)
                            end if
                            
                            if((Vrel*Sigma)/(Smax).gt.Rf) then
                                call collide(local_variables(i_prop,2:4),local_variables(j_prop_loc,2:4),Vrel,vel_rel_vec,collision_model,alpha)
                                cell_num_collisions = cell_num_collisions + 1
                                velocity(local_particles(i_prop),2:4) = local_variables(i_prop,2:4)
                                velocity(local_particles(j_prop_loc),2:4) = local_variables(j_prop_loc,2:4)
                            end if
                            
                        end do
                        
                        !print*, '--------------------------------------------------------------'
                        !do l=1,current_nop
                        !    print*, local_particles(l),sub_cell_indexes(local_particles(l),:)
                        !end do

                        if(output_option) then
                            write(*,2000), 'Current # of particles is ', current_nop
                            2000 format (A26,1X,I6)

                            write(*,3000), 'Candidate collisions are' , int(num_cand_collisions)
                            3000 format (A24,1X,I8)

                            write(*,4000), 'A total of ', cell_num_collisions, ' collisions happened for cell  ', i,j,k, '  (', cell_num_collisions/num_cand_collisions*100, '%)'
                            4000 format (A11,I5,A30,I3,I3,I3,A3,F4.1,A2)

                            write(*,5000), 'Vrel_max is ', Vrel_max(counter), ' Vrel is ', Vrel
                            5000 format (A11,F8.2,A8,F8.2)
                        end if

                        num_collisions = num_collisions + cell_num_collisions

                        if (timestep.ge.sampling_start .and. modulo(timestep,sampling_frequency).eq.0) then
                            call sample(velocity(local_particles,:), current_nop, num_cells, counter, sampled_values, weight, dx, dy, dz, boltz, mol_mass, eff_num, k, j,&
                                number_of_sampled_variables)
                            do m=1,number_of_sampled_variables
                                updated_sampled_values(counter,m) = updated_sampled_values(counter,m) + sampled_values(counter,m)
                            end do
                        end if
                        counter = counter + 1
                    end do
                end do
            end do
            if (timestep.ge.sampling_start .and. modulo(timestep,sampling_frequency).eq.0) then
                samples = samples + 1
                residual_counter = residual_counter + 1
                call convergence(old_sampled_values,updated_sampled_values,average_residual,residual_counter,num_cells,number_of_sampled_variables,num_timesteps,total_trnts,&
                                sampling_frequency,sampling_start)
                old_sampled_values = updated_sampled_values
                !print*, average_residual
            end if
            call cpu_time(timer_finish)
            write(*,1000), 'Timestep ', timestep,' of ', num_timesteps, ', transient timestep ', trnt_t_step, ' of ',total_trnts, ', real time (s) ', real_time,&
                            ', elapsed time: ', (timer_finish)/60, ' min.',', remaining time: ', (timer_finish)/60*(num_timesteps*total_trnts/timestep-1), ' min. (',&
                            (timer_finish)/60*(num_timesteps*total_trnts/timestep-1)/60,'h)'
            1000 format (A8,1X,I5,1X,A4,1X,I5,1X,A21,I2,A4,I2,A15,E10.2,A16,F10.2,A5,A18,F8.2,A7,F4.1,A2)
            
        !call parameter_check(num_particles,num_cells_x,num_cells_y,num_cells_z,dx,dy,dz,dt,volume,velocity,xsize,ysize,zsize,position)
        end do

        call parameter_check(num_particles,num_cells_x,num_cells_y,num_cells_z,dx,dy,dz,dt,volume,velocity,xsize,ysize,zsize,position)

        write(filename,"(A15)") 'position.dat'
        open(1,file=filename)
        do j=1,num_particles
            write(1,'(F14.6,F17.12,F17.12,F17.12)') position(j,:)
        end do
        close(1)
        print*, '-------------------------------------------------------------------------------'
    end do
    call time(actual_time)
    write(*,"('Time at end of execution ', A8)") actual_time

    write(filename,"(A15)") 'residuals.dat'
    open(222,file=filename)
    do j=1,size(average_residual)
        write(222,'(I5,F13.9)') j, average_residual(j)
    end do
    close(222)

    print*, samples*transient_samples
    
    updated_sampled_values = updated_sampled_values/samples
    updated_sampled_values = updated_sampled_values/transient_samples

    deallocate(local_particles)
    deallocate(local_variables)

    call output(updated_sampled_values,num_cells,num_cells_x,num_cells_y,num_cells_z,node_locs,number_of_sampled_variables)

    call cpu_time(timer_finish)
    write(*,10000) "Total execution time is ", (timer_finish-timer_start)/60, " minutes"
    10000 format (A24,1X,F8.2,1X,A8)

end program dsmc3d

subroutine convergence(old_sampled_values,updated_sampled_values,average_residual,residual_counter,num_cells,number_of_sampled_variables,num_timesteps,total_trnts,&
                        sampling_frequency,sampling_start)

    implicit none
    integer :: i, num_cells, residual_counter, num_timesteps, total_trnts, sampling_frequency, sampling_start, number_of_sampled_variables
    real*8, dimension(num_cells,number_of_sampled_variables) :: old_sampled_values, updated_sampled_values
    real*8, dimension(num_cells) :: res_errors
    real*8, dimension(total_trnts*((num_timesteps-sampling_start)/sampling_frequency+1)) :: average_residual

    do i=1,num_cells
        res_errors(i) = abs((old_sampled_values(i,2)-updated_sampled_values(i,2))/(old_sampled_values(i,2)))
    end do
    
    average_residual(residual_counter) = sum(res_errors)/num_cells
    

end subroutine convergence