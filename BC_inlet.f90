! Subroutine to impose boundary conditions
subroutine BC_inlet ()
use global_parameters 
!use para
use some_functions
implicit none
real(kind=dp) :: press_mms, uvel_mms, vvel_mms, rho_mms



do j=1, jmax
    V_vec(0,j,4)  = 2*V_vec(1,j,1)-V_vec(2,j,1)
    V_vec(0,j,1)  = V_vec(1,j,1)*(V_vec(0,j,4)/V_vec(1,j,4))
    V_vec(0,j,2)  = -(V_vec(1,j,2)*facex(1,j,3)+V_vec(1,j,3)*facex(1,j,4))*facex(1,j,3) &
                    +(V_vec(1,j,2)*facex(1,j,5)+V_vec(1,j,3)*facex(1,j,6))*facex(1,j,5)
    V_vec(0,j,3)  = -(V_vec(1,j,2)*facex(1,j,3)+V_vec(1,j,3)*facex(1,j,4))*facex(1,j,4) &
                    +(V_vec(1,j,2)*facex(1,j,5)+V_vec(1,j,3)*facex(1,j,6))*facex(1,j,6)
!!
!    print*, V_vec(0,j,2), V_vec(0,j,3)
!    pause 
    V_vec(-1,j,4) = 2*V_vec(0,j,1)-V_vec(1,j,1)
    V_vec(-1,j,1) = V_vec(0,j,1)*(V_vec(-1,j,4)/V_vec(0,j,4))
    V_vec(-1,j,2) = -(V_vec(2,j,2)*facex(1,j,3)+V_vec(2,j,3)*facex(1,j,4))*facex(1,j,3) &
                    +(V_vec(2,j,2)*facex(1,j,5)+V_vec(2,j,3)*facex(1,j,6))*facex(1,j,5)
    V_vec(-1,j,3) = -(V_vec(2,j,2)*facex(1,j,3)+V_vec(2,j,3)*facex(1,j,4))*facex(1,j,4) &
                    +(V_vec(2,j,2)*facex(1,j,5)+V_vec(2,j,3)*facex(1,j,6))*facex(1,j,6)
end do

do j=1, jmax
    V_vec(imax+1,j,1)  = 2*V_vec(imax,j,1)-V_vec(imax-1,j,1)
    V_vec(imax+1,j,2)  = 2*V_vec(imax,j,2)-V_vec(imax-1,j,2)
    V_vec(imax+1,j,3)  = 2*V_vec(imax,j,3)-V_vec(imax-1,j,3)
    V_vec(imax+1,j,4)  = 2*V_vec(imax,j,4)-V_vec(imax-1,j,4)

    V_vec(imax+2,j,1)  = 2*V_vec(imax+1,j,1)-V_vec(imax,j,1)
    V_vec(imax+2,j,2)  = 2*V_vec(imax+1,j,2)-V_vec(imax,j,2)
    V_vec(imax+2,j,3)  = 2*V_vec(imax+1,j,3)-V_vec(imax,j,3)
    V_vec(imax+2,j,4)  = 2*V_vec(imax+1,j,4)-V_vec(imax,j,4)
end do


!do i=1, imax
!    V_vec(i,0,4)   = 2*V_vec(i,1,4)-V_vec(i,2,4)
!    V_vec(i,0,1)   = V_vec(i,1,1)*(V_vec(i,0,4)/V_vec(i,1,4))
!    V_vec(i,0,2)   = -(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,3) &
!                     +(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,5)
!    V_vec(i,0,3)   = -(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,4) &
!                     +(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,6)
!    
!
!    V_vec(i,-1,4)  = 2*V_vec(i,0,4)-V_vec(i,1,4)
!    V_vec(i,-1,1)  = V_vec(i,0,1)*(V_vec(i,-1,4)/V_vec(i,0,4))
!    V_vec(i,-1,2)  = -(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,3) &
!                     +(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,5)
!    V_vec(i,-1,3)  = -(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,4) &
!                     +(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,6)
!end do
!
!do i=1, imax
!    V_vec(i,jmax+1,4)  = 2*V_vec(i,jmax,4)-V_vec(i,jmax-1,4)
!    V_vec(i,jmax+1,1)  = V_vec(i,jmax,1)*(V_vec(i,jmax+1,4)/V_vec(i,jmax,4))
!    V_vec(i,jmax+1,2)  = -(V_vec(i,jmax,2)*facey(i,jmax+1,3)+V_vec(i,jmax,3)*facey(i,jmax+1,4))*facey(i,jmax+1,3) &
!                         +(V_vec(i,jmax,2)*facey(i,jmax+1,5)+V_vec(i,jmax,3)*facey(i,jmax+1,6))*facey(i,jmax+1,5)
!    V_vec(i,jmax+1,3)  = -(V_vec(i,jmax,2)*facey(i,jmax+1,3)+V_vec(i,jmax,3)*facey(i,jmax+1,4))*facey(i,jmax+1,4) &
!                         +(V_vec(i,jmax,2)*facey(i,jmax+1,5)+V_vec(i,jmax,3)*facey(i,jmax+1,6))*facey(i,jmax+1,6)
!    
!    V_vec(i,jmax+2,4)  = 2*V_vec(i,jmax+1,4)-V_vec(i,jmax,4)
!    V_vec(i,jmax+2,1)  = V_vec(i,jmax+1,1)*(V_vec(i,jmax+2,4)/V_vec(i,jmax+1,4))
!    V_vec(i,jmax+2,2)  = -(V_vec(i,jmax-1,2)*facey(i,jmax+1,3)+V_vec(i,jmax-1,3)*facey(i,jmax+1,4))*facey(i,jmax+1,3) &
!                         +(V_vec(i,jmax-1,2)*facey(i,jmax+1,5)+V_vec(i,jmax-1,3)*facey(i,jmax+1,6))*facey(i,jmax+1,5)
!    V_vec(i,jmax+2,3)  = -(V_vec(i,jmax-1,2)*facey(i,jmax+1,3)+V_vec(i,jmax-1,3)*facey(i,jmax+1,4))*facey(i,jmax+1,4) &
!                         +(V_vec(i,jmax-1,2)*facey(i,jmax+1,5)+V_vec(i,jmax-1,3)*facey(i,jmax+1,6))*facey(i,jmax+1,6)
!end do


do i=1, imax
if (facey(i,1,1) > 0.347) then 
    V_vec(i,0,4)   = 2*V_vec(i,1,4)-V_vec(i,2,4)
    V_vec(i,0,1)   = V_vec(i,1,1)*(V_vec(i,0,4)/V_vec(i,1,4))
    V_vec(i,0,2)   = -(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,3) &
                     +(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,5)
    V_vec(i,0,3)   = 0.0_dp !-(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,4) &
                     !+(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,6)

    V_vec(i,-1,4)  = 2*V_vec(i,0,4)-V_vec(i,1,4)
    V_vec(i,-1,1)  = V_vec(i,0,1)*(V_vec(i,-1,4)/V_vec(i,0,4))
    V_vec(i,-1,2)  = -(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,3) &
                     +(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,5)
    V_vec(i,-1,3)  = 0.0_dp !-(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,4) &
                     !+(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,6)
else
    V_vec(i,0,4)   = 2*V_vec(i,1,4)-V_vec(i,2,4)
    V_vec(i,0,1)   = V_vec(i,1,1)*(V_vec(i,0,4)/V_vec(i,1,4))
    V_vec(i,0,2)   = -(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,3) &
                     +(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,5)
    V_vec(i,0,3)   = -(V_vec(i,1,2)*facey(i,1,3)+V_vec(i,1,3)*facey(i,1,4))*facey(i,1,4) &
                     +(V_vec(i,1,2)*facey(i,1,5)+V_vec(i,1,3)*facey(i,1,6))*facey(i,1,6)

    V_vec(i,-1,4)  = 2*V_vec(i,0,4)-V_vec(i,1,4)
    V_vec(i,-1,1)  = V_vec(i,0,1)*(V_vec(i,-1,4)/V_vec(i,0,4))
    V_vec(i,-1,2)  = -(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,3) &
                     +(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,5)
    V_vec(i,-1,3)  = -(V_vec(i,2,2)*facey(i,1,3)+V_vec(i,2,3)*facey(i,1,4))*facey(i,1,4) &
                     +(V_vec(i,2,2)*facey(i,1,5)+V_vec(i,2,3)*facey(i,1,6))*facey(i,1,6)
endif
end do

do i=1, imax
if (facey(i,jmax+1,3) < -0.01) then
    V_vec(i,jmax+1,4)  = p_in
    V_vec(i,jmax+1,1)  = rho_in
    V_vec(i,jmax+1,2)  = vel_in
    V_vec(i,jmax+1,3)  = 0.0_dp
        
    V_vec(i,jmax+2,4)  = p_in
    V_vec(i,jmax+2,1)  = rho_in
    V_vec(i,jmax+2,2)  = vel_in
    V_vec(i,jmax+2,3)  = 0.0_dp
else
    V_vec(i,jmax+1,4)  = V_vec(i,jmax,4) !2*V_vec(i,jmax,4)-V_vec(i,jmax-1,4)
    V_vec(i,jmax+1,1)  = V_vec(i,jmax,1)*(V_vec(i,jmax+1,4)/V_vec(i,jmax,4))
    V_vec(i,jmax+1,2)  = -(V_vec(i,jmax,2)*facey(i,jmax+1,3)+V_vec(i,jmax,3)*facey(i,jmax+1,4))*facey(i,jmax+1,3) &
                         +(V_vec(i,jmax,2)*facey(i,jmax+1,5)+V_vec(i,jmax,3)*facey(i,jmax+1,6))*facey(i,jmax+1,5)
    V_vec(i,jmax+1,3)  = 0.0_dp !-(V_vec(i,jmax,2)*facey(i,jmax+1,3)+V_vec(i,jmax,3)*facey(i,jmax+1,4))*facey(i,jmax+1,4) &
                         !+(V_vec(i,jmax,2)*facey(i,jmax+1,5)+V_vec(i,jmax,3)*facey(i,jmax+1,6))*facey(i,jmax+1,6)
        
    V_vec(i,jmax+2,4)  = V_vec(i,jmax+1,4) !2*V_vec(i,jmax+1,4)-V_vec(i,jmax,4)
    V_vec(i,jmax+2,1)  = V_vec(i,jmax+1,1)*(V_vec(i,jmax+2,4)/V_vec(i,jmax+1,4))
    V_vec(i,jmax+2,2)  = -(V_vec(i,jmax-1,2)*facey(i,jmax+1,3)+V_vec(i,jmax-1,3)*facey(i,jmax+1,4))*facey(i,jmax+1,3) &
                         +(V_vec(i,jmax-1,2)*facey(i,jmax+1,5)+V_vec(i,jmax-1,3)*facey(i,jmax+1,6))*facey(i,jmax+1,5)
    V_vec(i,jmax+2,3)  = 0.0_dp !-(V_vec(i,jmax-1,2)*facey(i,jmax+1,3)+V_vec(i,jmax-1,3)*facey(i,jmax+1,4))*facey(i,jmax+1,4) &
                         !+(V_vec(i,jmax-1,2)*facey(i,jmax+1,5)+V_vec(i,jmax-1,3)*facey(i,jmax+1,6))*facey(i,jmax+1,6)
endif
end do

    V_vec(0,0,:)           = V_vec(1,1,:)
    V_vec(imax+1,jmax+1,:) = V_vec(imax,jmax,:)
    V_vec(imax+1,0,:)      = V_vec(imax,1,:)
    V_vec(0,jmax+1,:)      = V_vec(0,jmax,:)

return
end
  
