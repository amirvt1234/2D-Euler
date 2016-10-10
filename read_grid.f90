subroutine read_grid()
use global_parameters   !, only: x, y, imax, jmax, kmax, rho, u, v, p, i, j, k
use set_precision, only : dp
use set_inputs
use some_functions
implicit none


real(kind=dp) :: rmassconv, energyconv, ymtmconv, xmtmconv
real(kind=dp) :: press_mms, uvel_mms, vvel_mms, rho_mms

open(unit=12,file='Inlet.105x33.grd')
read(12,*) 
read(12,*) imax, jmax, kmax

imax = imax-1
jmax = jmax-1
kmax = kmax-1

!imax=8
!jmax=imax

allocate(x(1:imax+1,1:jmax+1), y(1:imax+1,1:jmax+1), &
         facex(1:imax+1,1:jmax+1,7), facey(1:imax+1,1:jmax+1,7) )   !!! nodes
allocate(flux_tot_x(1:imax+1,1:jmax+1,dim), flux_tot_y(1:imax+1,1:jmax+1,dim))
allocate(cell(1:imax,1:jmax,3), source_mms(1:imax,1:jmax,dim)) 
allocate(U_vec(-1:imax+2,-1:jmax+2,dim), U_vec_old(-1:imax+2,-1:jmax+2,dim), &
         V_vec(-1:imax+2,-1:jmax+2,dim), res(-1:imax+2,-1:jmax+2,dim))
allocate(V_vec_mms_c(1:imax,1:jmax,dim), V_vec_node(1:imax+1,1:jmax+1,dim) )
allocate(cell_dummy_ghost(-1:imax+2,-1:jmax+2,2))
allocate(errord(1:imax,1:jmax,dim))


 
read(12,*) (((x(i,j),i=1,imax+1),j=1,jmax+1),k=1,kmax+1),  &
           (((y(i,j),i=1,imax+1),j=1,jmax+1),k=1,kmax+1)



!do j=1,jmax+1
!    do i=1,imax+1
!        x(i,j) = 0.0_dp + (1.0_dp/(imax))*(i-1)
!        y(i,j) = 0.0_dp + (1.0_dp/(imax))*(j-1)
!    enddo
!enddo

close(12)




open(unit=122,file='dimension.out')
write(122,*) imax, jmax
close(122)
call initialize_constants

U_vec=0.0_dp
U_vec_old = 0.0_dp

!rho = 1.0_dp
!u   = 800.0_dp
!v   = 800.0_dp
!p   = 100000.0_dp


T_in     = 217.0_dp
p_in     = 12270.0_dp 
mach_in  = 4.0_dp                                      
rho_in   = p_in/(R*T_in)
vel_in     = mach_in*dsqrt(gamma*p_in/rho_in)

!angle    = 8.0_dp
!T_in     = 300.0_dp
!p_in     = 67243.5_dp !65855.8_dp !
!mach_in  = 0.82_dp    !0.84_dp                                      
!rho_in   = p_in/(R*T_in)
!vel_in     = mach_in*dsqrt(gamma*p_in/rho_in)

V_vec(:,:,1) = 1.0_dp
V_vec(:,:,2) = 800.0_dp
V_vec(:,:,3) = 800.0_dp
V_vec(:,:,4) = 100000.0_dp

V_vec(:,:,1) = rho_in
V_vec(:,:,2) = vel_in*dcos(angle*pi/180._dp)
V_vec(:,:,3) = vel_in*dsin(angle*pi/180._dp)
V_vec(:,:,4) = p_in !100000.0_dp


do j=1,jmax
    do i=1,imax
        U_vec(i,j,:) =VtoU(V_vec(i,j,:))
    end do
end do 




do j=1, jmax
    do i=1, imax 
        cell(i,j,1)  =(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))/4  ! x coordinate of cell
        cell(i,j,2)  =(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))/4  ! y coordinate of cell
        cell(i,j,3)  =0.5*dabs((x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i+1,j))-(y(i+1,j+1)-y(i,j))*(x(i,j+1)-x(i+1,j))) ! volume of cell

        source_mms(i,j,1) = rmassconv (ilength,cell(i,j,1),cell(i,j,2))
        source_mms(i,j,2) = xmtmconv  (ilength,cell(i,j,1),cell(i,j,2))
        source_mms(i,j,3) = ymtmconv  (ilength,cell(i,j,1),cell(i,j,2))
        source_mms(i,j,4) = energyconv(gamma,ilength,cell(i,j,1),cell(i,j,2))
      
    end do
end do
cell_dummy_ghost = 0.0_dp
cell_dummy_ghost(1:imax,1:jmax,1) = cell(1:imax, 1:jmax, 1)
cell_dummy_ghost(1:imax,1:jmax,2) = cell(1:imax, 1:jmax, 2)
cell_dummy_ghost(0,1:jmax,1) = 2*cell(1,1:jmax,1)-cell(2,1:jmax,1)
cell_dummy_ghost(0,1:jmax,2) = 2*cell(1,1:jmax,2)-cell(2,1:jmax,2)
cell_dummy_ghost(-1,1:jmax,1) = 2*cell_dummy_ghost(0,1:jmax,1)-cell(1,1:jmax,1)
cell_dummy_ghost(-1,1:jmax,2) = 2*cell_dummy_ghost(0,1:jmax,2)-cell(1,1:jmax,2)
cell_dummy_ghost(imax+1,1:jmax,1) = 2*cell(imax,1:jmax,1)-cell(imax-1,1:jmax,1)
cell_dummy_ghost(imax+1,1:jmax,2) = 2*cell(imax,1:jmax,2)-cell(imax-1,1:jmax,2)
cell_dummy_ghost(imax+2,1:jmax,1) = 2*cell_dummy_ghost(imax+1,1:jmax,1)-cell(imax,1:jmax,1)
cell_dummy_ghost(imax+2,1:jmax,2) = 2*cell_dummy_ghost(imax+1,1:jmax,2)-cell(imax,1:jmax,2)
cell_dummy_ghost(1:imax,0,1) = 2*cell(1:imax,1,1)-cell(1:imax,2,1)
cell_dummy_ghost(1:imax,0,2) = 2*cell(1:imax,1,2)-cell(1:imax,2,2)
cell_dummy_ghost(1:imax,-1,1) = 2*cell_dummy_ghost(1:imax,0,1)-cell(1:imax,1,1)
cell_dummy_ghost(1:imax,-1,2) = 2*cell_dummy_ghost(1:imax,0,2)-cell(1:imax,1,2)
cell_dummy_ghost(1:imax,jmax+1,1) = 2*cell(1:imax,jmax,1)-cell(1:imax,jmax-1,1)
cell_dummy_ghost(1:imax,jmax+1,2) = 2*cell(1:imax,jmax,2)-cell(1:imax,jmax-1,2)
cell_dummy_ghost(1:imax,jmax+2,1) = 2*cell_dummy_ghost(1:imax,jmax+1,1)-cell(1:imax,jmax,1)
cell_dummy_ghost(1:imax,jmax+2,2) = 2*cell_dummy_ghost(1:imax,jmax+1,2)-cell(1:imax,jmax,2)

do j=1, jmax
    do i=1, imax         
        V_vec_mms_c(i,j,1) = rho_mms   (ilength,cell(i,j,1),cell(i,j,2))
        V_vec_mms_c(i,j,2) = uvel_mms  (ilength,cell(i,j,1),cell(i,j,2))
        V_vec_mms_c(i,j,3) = vvel_mms  (ilength,cell(i,j,1),cell(i,j,2))
        V_vec_mms_c(i,j,4) = press_mms (ilength,cell(i,j,1),cell(i,j,2))        
    end do
end do

do j=1, jmax
    do i=1, imax+1
        facex     (i,j,1)=(x(i,j)+x(i,j+1))/2    ! x location
        facex     (i,j,2)=(y(i,j)+y(i,j+1))/2    ! y location
        facex     (i,j,7)=dsqrt((y(i,j+1)-y(i,j))**2+(x(i,j+1)-x(i,j))**2) ! area of face
        facex     (i,j,3)= (y(i,j+1)-y(i,j))/facex(i,j,7)  ! nx
        facex     (i,j,4)=-(x(i,j+1)-x(i,j))/facex(i,j,7) ! ny
        facex     (i,j,5)= (x(i,j+1)-x(i,j))/facex(i,j,7)  ! tx
        facex     (i,j,6)= (y(i,j+1)-y(i,j))/facex(i,j,7) ! ty
    end do
end do


do j=1, jmax+1
    do i=1, imax
        facey     (i,j,1)=(x(i,j)+x(i+1,j))/2 !x location
        facey     (i,j,2)=(y(i,j)+y(i+1,j))/2 !y location
        facey     (i,j,7)=dsqrt((y(i+1,j)-y(i,j))**2+(x(i+1,j)-x(i,j))**2) ! area of face
        facey     (i,j,3)= (y(i,j)-y(i+1,j))/facey(i,j,7) ! nx
        facey     (i,j,4)=-(x(i,j)-x(i+1,j))/facey(i,j,7) ! ny
        facey     (i,j,5)= (x(i,j)-x(i+1,j))/facey(i,j,7) ! nx
        facey     (i,j,6)= (y(i,j)-y(i+1,j))/facey(i,j,7) ! ny        
    end do
end do

!open(30,file='grido.dat',status='unknown')
!
!do j = 1, jmax+1
!    do i = 1, imax+1
!        write(30,19) x(i,j), y(i,j)
!    end do
!end do
!19 FORMAT (6(1X,F24.14))
!open(40,file='cell_c.dat',status='unknown')
!
!do j = 1, jmax
!    do i = 1, imax
!        write(40,20) cell(i,j,1), cell(i,j,2), cell(i,j,3), V_vec_mms_c(i,j,1), V_vec_mms_c(i,j,2), V_vec_mms_c(i,j,3), V_vec_mms_c(i,j,4)
!    end do
!end do
!
!close(40)
!20 FORMAT (7(1X,F24.14))
!
!open(50,file='facex.dat',status='unknown')
!
!do j=1, jmax
!    do i=1, imax+1
!        write(50,21) facex(i,j,1), facex(i,j,2), facex(i,j,3), facex(i,j,4) , facex(i,j,5), facex(i,j,6), facex(i,j,7) 
!    end do
!end do
!
!open(51,file='facey.dat',status='unknown')
!
!do j=1, jmax+1
!    do i=1, imax
!        write(51,21) facey(i,j,1), facey(i,j,2), facey(i,j,3), facey(i,j,4) , facey(i,j,5), facey(i,j,6), facey(i,j,7) 
!    end do
!end do
!
!close(50)
!close(51)
!21 FORMAT (7(F24.14))

end subroutine read_grid
