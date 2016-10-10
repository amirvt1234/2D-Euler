subroutine write_output_airfoil()

use global_parameters
use set_inputs
implicit none
real(kind=dp), dimension(1:imax,2) :: press_coefl, press_coefu
real(kind=dp) ::  wingl, lup, ldown, lift_force, drag_force, drag_coef, lift_coef, forcex, forcey

press_coefl =0.0_dp ;press_coefu =0.0_dp

do j=1,jmax+1
   do i=1,imax+1 
      V_vec_node(i,j,:) = (V_vec(i-1,j-1,:)+V_vec(i-1,j,:)+V_vec(i,j-1,:)+V_vec(i,j,:))/4.0_dp
   end do
end do

open(10,file='NS.dat',status='unknown',form='formatted')   !
do j = 1,jmax+1
   do i = 1,imax+1
      write(10,10) x(i,j),y(i,j), V_vec_node(i,j,1), V_vec_node(i,j,2), V_vec_node(i,j,3), V_vec_node(i,j,4)
   end do
end do

10 FORMAT (6(1X,F24.14))
close(10)


j=0
k=-1
wingl =0.0_dp; lup=0.0_dp; ldown = 0.0_dp
do i=1, imax
  if (x(i,1) < 0.152401) then
	if(y(i,1)<-0.0000000001_dp) then 
         j=j+1
	 press_coefl(j,1) = x(i,1) !ldown
	 press_coefl(j,2) = (V_vec_node(i,1,4)+0.5*V_vec_node(i,1,1)*  &
                            (V_vec_node(i,1,2)**2+V_vec_node(i,1,3)**2)-p_in)/(0.5*rho_in*vel_in**2)
         press_coefl(j,2) = (V_vec_node(i,1,4)-p_in)/(0.5*rho_in*vel_in**2)
	 ldown = ldown + dsqrt((x(i+1,1)-x(i,1))**2+(y(i+1,1)-y(i,1))**2)
	else
         k=k+1
         wingl = wingl + dsqrt(x(i,1)**2+y(i,1)**2)
	 press_coefu(k,1) = lup
	 press_coefu(k,2) = (V_vec_node(i,1,4)+0.5*V_vec_node(i,1,1)*  &
                            (V_vec_node(i,1,2)**2+V_vec_node(i,1,3)**2)-p_in)/(0.5*rho_in*vel_in**2)
        press_coefu(k,2) = (V_vec_node(i,1,4)-p_in)/(0.5*rho_in*vel_in**2)
	 if (k .ge. 1)  lup = lup + dsqrt((x(i+1,1)-x(i,1))**2+(y(i+1,1)-y(i,1))**2)
	end if
  end if
end do
do i=1,j
	press_coefl(i,1) = press_coefu(k-i,1)
end do
lift_force = 0.0_dp ; drag_force = 0.0_dp; forcey=0.0_dp; forcex=0.0_dp
do i=1, imax
  if (facey(i,1,1) < 0.152401) then 
        forcey = forcey + 0.5*(V_vec_node(i,1,4)+V_vec_node(i+1,1,4))*facey(i,1,4)*facey(i,1,7)
        forcex = forcex + 0.5*(V_vec_node(i,1,4)+V_vec_node(i+1,1,4))*facey(i,1,3)*facey(i,1,7)
  end if
end do
forcex = -forcex; forcey = -forcey
drag_force =  forcex*dcos(angle*pi/180.0)+forcey*dsin(angle*pi/180.0)
lift_force =  forcey*dcos(angle*pi/180.0)-forcex*dsin(angle*pi/180.0)
!pause
drag_coef = drag_force/(0.5*rho_in*vel_in**2*0.152400)
lift_coef = lift_force/(0.5*rho_in*vel_in**2*0.152400)

open(13,file='press_coefl.dat',status='unknown',form='formatted')   !
open(14,file='press_coefu.dat',status='unknown',form='formatted')   !
open(15,file='dragandlift.dat',status='unknown',form='formatted')   !
open(16,file='dragandlift_time.dat',status='unknown',form='formatted')   !
write(16,*) drag_coef, lift_coef
write(15,*) drag_coef, lift_coef

do i=1,j
      write(13,13)  press_coefl(i,1)/press_coefu(k,1), press_coefl(i,2)
end do
do i=1,k
      write(14,13)  press_coefu(i,1)/press_coefu(k,1), press_coefu(i,2)
end do

13 FORMAT (3(1X,F24.14))
!15 FORMAT (4(1X,F24.14))
close(13)
close(14)   
close(15)

return
end subroutine write_output_airfoil
