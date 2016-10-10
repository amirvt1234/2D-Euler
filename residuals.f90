subroutine residuals(alpha_RK)

use global_parameters
use some_functions
use fluxes

implicit none
real(kind=dp) :: alpha_RK


! calls subroutine to calculate dissipations
call BC_inlet()

!call write_output()
!pause

call flux_param(V_vec,flux_tot_x,flux_tot_y)
!res_norm      = f_norm(res(1:imax,1:jmax,:), imax, jmax)
source_mms = 0.0_dp
res = 0.0_dp
do j=1,jmax
    do i=1,imax

        res(i,j,1)  = flux_tot_x(i,j,1)*facex(i,j,7)-flux_tot_x(i+1,j,1)*facex(i+1,j,7) &
                     +flux_tot_y(i,j,1)*facey(i,j,7)-flux_tot_y(i,j+1,1)*facey(i,j+1,7) &
                     +source_mms(i,j,1)*cell(i,j,3)   ! residual of eq1
        res(i,j,2)  = flux_tot_x(i,j,2)*facex(i,j,7)-flux_tot_x(i+1,j,2)*facex(i+1,j,7) &
                     +flux_tot_y(i,j,2)*facey(i,j,7)-flux_tot_y(i,j+1,2)*facey(i,j+1,7) &
                     +source_mms(i,j,2) *cell(i,j,3)   ! residual of eq1
        res(i,j,3)  = flux_tot_x(i,j,3)*facex(i,j,7)-flux_tot_x(i+1,j,3)*facex(i+1,j,7) &
                     +flux_tot_y(i,j,3)*facey(i,j,7)-flux_tot_y(i,j+1,3)*facey(i,j+1,7) &
                     +source_mms(i,j,3) *cell(i,j,3)   ! residual of eq1
        res(i,j,4)  = flux_tot_x(i,j,4)*facex(i,j,7)-flux_tot_x(i+1,j,4)*facex(i+1,j,7) &
                     +flux_tot_y(i,j,4)*facey(i,j,7)-flux_tot_y(i,j+1,4)*facey(i,j+1,7) &
                     +source_mms(i,j,4)*cell(i,j,3)   ! residual of eq1

        dt   = CFL*(dabs(facex(i,j,1)-facex(i+1,j,1)))*(dabs(facey(i,j+1,2)-facey(i,j,2)))/  &
               ((dsqrt(dabs(gamma*V_vec(i,j,4)/V_vec(i,j,1)))+dabs(V_vec(i,j,2)))* &
                (dabs(facey(i,j+1,2)-facey(i,j,2)))                                   +&
                (dsqrt(dabs(gamma*V_vec(i,j,4)/V_vec(i,j,1)))+dabs(V_vec(i,j,3)))* &
                (dabs(facex(i+1,j,1)-facex(i,j,1)))                                  )
!        if (niter >60000) then
!            if (res(i,j,1) > 1.0*res_norm(2,1))  res(i,j,1)=1.0*res_norm(2,1)
!            if (res(i,j,2) > 1.0*res_norm(2,2))  res(i,j,2)=1.0*res_norm(2,2)
!            if (res(i,j,3) > 1.0*res_norm(2,3))  res(i,j,3)=1.0*res_norm(2,3)
!            if (res(i,j,4) > 1.0*res_norm(2,4))  res(i,j,4)=1.0*res_norm(2,4)
!        end if


        U_vec(i,j,1)    = U_vec_old(i,j,1)+alpha_RK*res(i,j,1)*(dt/cell(i,j,3))   ! conserved variable : u1
        U_vec(i,j,2)    = U_vec_old(i,j,2)+alpha_RK*res(i,j,2)*(dt/cell(i,j,3))   ! conserved variable : u1
        U_vec(i,j,3)    = U_vec_old(i,j,3)+alpha_RK*res(i,j,3)*(dt/cell(i,j,3))   ! conserved variable : u1
        U_vec(i,j,4)    = U_vec_old(i,j,4)+alpha_RK*res(i,j,4)*(dt/cell(i,j,3))

        ! calculate primitive variables from conserved variables
        V_vec(i,j,:) = UtoV(U_vec(i,j,:))
        if (V_vec(i,j,1) < 0.1 ) V_vec(i,j,1)=0.1
!        if (V_vec(i,j,2) < 0.10 ) V_vec(i,j,2)=0.10
        if (V_vec(i,j,4) < 10.0 ) V_vec(i,j,4)=10.0

!        if (V_vec(i,j,1) > 8.0 ) V_vec(i,j,1)=8.0
!        if (abs(V_vec(i,j,2)) > 3.0*vel_in ) V_vec(i,j,2)=3.0*vel_in*sign(1.0_dp,V_vec(i,j,2))
!        if (abs(V_vec(i,j,3)) > 3.0*vel_in ) V_vec(i,j,3)=3.0*vel_in*sign(1.0_dp,V_vec(i,j,3))
!        if (V_vec(i,j,4) > 4.0*10.0**7 ) V_vec(i,j,4)=4.0*10.0**7
        U_vec(i,j,:) =VtoU(V_vec(i,j,:))

    end do
end do



return
end 
