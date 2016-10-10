!
!  CFD_ZAINALI_AMIR_HW2
!
program main

use global_parameters
use some_functions
use set_precision, only : dp
implicit none
real(kind=dp):: start_time, end_time
call cpu_time(start_time)

call read_grid()

!call BC_inlet()   ! Imposes Boundary conditions
!call residuals(0.5_dp) 
niter= 0


!call residuals(1.0_dp)
dt = 1._dp*10._dp**(-7)

do niter=1,nmax
	
    
!    do j=1,jmax
!        do i=1,imax
!            U_vec(i,j,:) =VtoU(V_vec(i,j,:))
!        end do
!    end do
    U_vec_old = U_vec
!    call residuals(1.0_dp/4.)
!    call residuals(1.0_dp/3.)
    
!    call residuals(1.0_dp/2.)
    call residuals(1.0_dp)
    


    if (niter > 0) then
    if (mod(niter,1000)==0 .OR. niter==1) then
        if (niter==1) then
            res_norm0 = f_norm(res(1:imax,1:jmax,:), imax, jmax)
        end if
        !res_norm_old  = res_norm
        res_norm      = f_norm(res(1:imax,1:jmax,:), imax, jmax)
        !delta_resnorm = abs(res_norm-res_norm_old) 
        open(1111,file='res.dat',status='unknown',form='formatted')
!        do j = 1,jmax+1
!           do i = 1,imax+1
                write(1111,133) res_norm/res_norm0
!           end do
!        end do
        
        !close(1111)
        print*, 'number of iterations=', niter
        print*, 'dt=', dt
        print*, 'L_{2norm}(U1)=', res_norm(2,1)/ res_norm0(2,1)
        print*, 'L_{2norm}(U2)=', res_norm(2,2)/ res_norm0(2,2)
        print*, 'L_{2norm}(U3)=', res_norm(2,3)/ res_norm0(2,3)
        print*, 'L_{2norm}(U4)=', res_norm(2,4)/ res_norm0(2,4)
        !pause
        call write_output_inlet()
        !if (maxval(res_norm(2,:)) < (10.**(-4)) .and.  maxval(delta_resnorm(2,:)) < (10.**(-14)) ) exit
        if (maxval(res_norm(2,:)/res_norm0(2,:)) < (10.**(-12))  ) exit
    end if
    end if
!    if (maxval(norm(3,:)) < (10.**(-4)) .and.  maxval(delta_norm(3,:)) < (10.**(-13)) ) exit

end do
133 FORMAT (13(1X,F24.14))
!call write_output_airfoil()
call cpu_time(end_time)
print*, '===================================='
print*, 'execution time=', end_time-start_time
print*, '===================================='             
end ! End 


