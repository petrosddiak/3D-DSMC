 module collision
        implicit none
    contains
        subroutine collide(velocity_i,velocity_j,Vrel,vel_rel_vec,collision_model,alpha)
            integer*4 :: k, collision_model
            real*8 :: theta, phi, chi, epsilon, R1, R2, Vrel, alpha
            real*8, dimension(3) :: velocity_i, velocity_j, Vcm, Vrel_post, Vel_post_i, Vel_post_j, vel_rel_vec
            character(len=1024) :: filename, name

            select case(collision_model)
                case(1) ! HS
                    call random_seed()
                    call random_number(R1)
                    call random_number(R2)
                    Vcm = (/(velocity_i(1) + velocity_j(1))*0.5,(velocity_i(2) + velocity_j(2))*0.5,(velocity_i(3) + velocity_j(3))*0.5/)
                    theta = acos(2*R1-1)
                    phi = 2*acos(-1.d0)*R2

                    Vrel_post(1) = Vrel*cos(theta)
                    Vrel_post(2) = Vrel*sin(theta)*cos(phi)
                    Vrel_post(3) = Vrel*sin(theta)*sin(phi)

                case(2) ! VHS
                    call random_seed()
                    call random_number(R1)
                    call random_number(R2)
                    Vcm = (/(velocity_i(1) + velocity_j(1))*0.5,(velocity_i(2) + velocity_j(2))*0.5,(velocity_i(3) + velocity_j(3))*0.5/)
                    theta = acos(2*R1-1)
                    phi = 2*acos(-1.d0)*R2

                    Vrel_post(1) = Vrel*cos(theta)
                    Vrel_post(2) = Vrel*sin(theta)*cos(phi)
                    Vrel_post(3) = Vrel*sin(theta)*sin(phi)

                case(3) ! VSS

                    call random_seed()
                    call random_number(R1)
                    call random_number(R2)
                    Vcm = (/(velocity_i(1) + velocity_j(1))*0.5,(velocity_i(2) + velocity_j(2))*0.5,(velocity_i(3) + velocity_j(3))*0.5/)
                    chi = acos(2*(R1)**(1/alpha)-1)
                    phi = 2*acos(-1.d0)*R2

                    Vrel_post(1) = cos(chi)*vel_rel_vec(1) + sin(chi)*sin(phi)*sqrt(vel_rel_vec(2)**2+vel_rel_vec(3)**2)             
                    Vrel_post(2) = cos(chi)*vel_rel_vec(2) + sin(chi)*(Vrel*vel_rel_vec(3)*cos(phi)-vel_rel_vec(1)*vel_rel_vec(2)*sin(phi))/sqrt(vel_rel_vec(2)**2+vel_rel_vec(3)**2)
                    Vrel_post(3) = cos(chi)*vel_rel_vec(3) - sin(chi)*(Vrel*vel_rel_vec(2)*cos(phi)+vel_rel_vec(1)*vel_rel_vec(3)*sin(phi))/sqrt(vel_rel_vec(2)**2+vel_rel_vec(3)**2)

            end select

            Vel_post_i(1) = Vcm(1) + 0.5*Vrel_post(1)
            Vel_post_i(2) = Vcm(2) + 0.5*Vrel_post(2)
            Vel_post_i(3) = Vcm(3) + 0.5*Vrel_post(3)

            Vel_post_j(1) = Vcm(1) - 0.5*Vrel_post(1)
            Vel_post_j(2) = Vcm(2) - 0.5*Vrel_post(2)
            Vel_post_j(3) = Vcm(3) - 0.5*Vrel_post(3)
            
            do k=1,3
                velocity_i(k) = Vel_post_i(k)
                velocity_j(k) = Vel_post_j(k)
            end do

        end subroutine collide
end module collision

