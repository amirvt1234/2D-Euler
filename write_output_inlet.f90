subroutine write_output_inlet()

use global_parameters
use some_functions
implicit none
real(kind=dp) :: pressure_lost, pexit


do j=1,jmax+1
   do i=1,imax+1 
      V_vec_node(i,j,:) = (V_vec(i-1,j-1,:)+V_vec(i-1,j,:)+V_vec(i,j-1,:)+V_vec(i,j,:))/4.0_dp
   end do
end do

open(10,file='NS.dat',status='unknown',form='formatted')   !
do j = 1,jmax+1
   do i = 1,imax+1
      write(10,10) x(i,j),y(i,j), V_vec_node(i,j,1), V_vec_node(i,j,2), &
                                  V_vec_node(i,j,3), V_vec_node(i,j,4)  !, &
!                                  V_vec_node(i,j,1)-V_vec_mms_c(i,j,1), &
!                                  V_vec_node(i,j,2)-V_vec_mms_c(i,j,2), &
!                                  V_vec_node(i,j,3)-V_vec_mms_c(i,j,3), &
!                                  V_vec_node(i,j,4)-V_vec_mms_c(i,j,4)
   end do
end do

10 FORMAT (10(1X,F24.14))
close(10)


!open(11,file='boundary.dat',status='unknown',form='formatted')   !
!do j = -1,jmax+2
!   do i = -1,imax+2
!      write(11,11)  cell_dummy_ghost(i,j,1), cell_dummy_ghost(i,j,2), V_vec(i,j,1), V_vec(i,j,2), &
!                                  V_vec(i,j,3), V_vec(i,j,4)  !, &
!!                                  V_vec_node(i,j,1)-V_vec_mms_c(i,j,1), &
!!                                  V_vec_node(i,j,2)-V_vec_mms_c(i,j,2), &
!!                                  V_vec_node(i,j,3)-V_vec_mms_c(i,j,3), &
!!                                  V_vec_node(i,j,4)-V_vec_mms_c(i,j,4)
!   end do
!end do
!
!11 FORMAT (10(1X,F24.14))
!close(11)
pexit = 0.0_dp
do j=1,jmax

    pexit = pexit + (0.5*(V_vec_node(imax+1,j,4)+V_vec_node(imax+1,j+1,4))+0.5*   &
             0.5*(V_vec_node(imax+1,j,1)+V_vec_node(imax+1,j+1,1))*       &   ! 
             ((0.5*(V_vec_node(imax+1,j,2)+V_vec_node(imax+1,j+1,2)))**2+   &
              (0.5*(V_vec_node(imax+1,j,3)+V_vec_node(imax+1,j+1,3)))**2)  )*facex(imax+1,j,7)

enddo

pexit = pexit/0.4_dp
!print*, pexit, (p_in+0.5*rho_in*vel_in**2)

!mms_error = f_norm(errord(1:imax,1:jmax,:),imax, jmax)
open(13,file='p_drop.dat',status='unknown',form='formatted')   !
open(14,file='p_drop_time_series.dat',status='unknown',form='formatted')   !

write(13,13) niter*1.0_dp, pexit-(p_in+0.5*rho_in*vel_in**2)
write(14,13) niter*1.0_dp, pexit-(p_in+0.5*rho_in*vel_in**2)
13 FORMAT (2(1X,F24.14))
close(13)

return
end subroutine write_output_inlet