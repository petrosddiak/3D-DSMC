module sampling
        implicit none
    contains
        subroutine sample(local_variables, current_nop, num_cells, counter, sampled_values, weight, dx, dy, dz, boltz, mol_mass, eff_num, k, j, number_of_sampled_variables)
            integer*4 :: current_nop, num_cells, counter, i, k, j, number_of_sampled_variables
            real*8, dimension(current_nop,4) :: local_variables
            real*8, dimension(current_nop,3) :: velocity_fluct
            real*8 :: weight, dx, dy, dz, boltz, density, eff_num, n_density, Ru, pressure, T, mol_mass, fluct_sum
            real*8, dimension(num_cells,number_of_sampled_variables) :: sampled_values
            real*8, dimension(3) :: velocity_mean

            do i=1,3
                velocity_mean(i) = sum(local_variables(:,i+1))/current_nop !x,y,z mean velocities for current cell
            end do

            do i=1,3
                velocity_fluct(:,i) = velocity_mean(i) - local_variables(:,i+1) !fluctuating components of each particle in current cell
            end do

            n_density = eff_num*current_nop/(dx*dy*dz) !local number density
            density = n_density*mol_mass
            Ru = 0.2081 !kJ/(kg*K)

            do i=1,current_nop
                fluct_sum = fluct_sum + velocity_fluct(i,1)**2+velocity_fluct(i,2)**2+velocity_fluct(i,3)**2
            end do

            pressure = density*fluct_sum/3.0/current_nop 
            T = pressure/(density*Ru*1000)
            
            sampled_values(counter,1) = sqrt(velocity_mean(1)**2+velocity_mean(2)**2+velocity_mean(3)**2)
            sampled_values(counter,2) = velocity_mean(1)
            sampled_values(counter,3) = velocity_mean(2)
            sampled_values(counter,4) = velocity_mean(3)
            sampled_values(counter,5) = n_density
            sampled_values(counter,6) = pressure
            sampled_values(counter,7) = T
            sampled_values(counter,8) = density

        end subroutine sample
end module sampling