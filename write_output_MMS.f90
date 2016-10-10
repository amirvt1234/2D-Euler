subroutine write_output_MMS()

use global_parameters
use some_functions
implicit none
real(kind=dp), dimension(3,dim) :: mms_error


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
do j=1,jmax
do i=1,imax
    errord(i,j,:) = VtoU(V_vec(i,j,:))-VtoU(V_vec_mms_c(i,j,:))
enddo
enddo

mms_error = f_norm(errord(1:imax,1:jmax,:),imax, jmax)
open(13,file='error_mms.dat',status='unknown',form='formatted')   !
open(14,file='error_mms_time_series.dat',status='unknown',form='formatted')   !

write(13,13) imax*1.0_dp, mms_error
write(14,13) imax*1.0_dp, mms_error
13 FORMAT (13(1X,F24.14))
close(13)

return
end subroutine write_output_MMS