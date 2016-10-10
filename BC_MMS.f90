! Subroutine to impose boundary conditions
subroutine BC_MMS ()
use global_parameters 
!use para
use some_functions
implicit none
real(kind=dp) :: press_mms, uvel_mms, vvel_mms, rho_mms



do j=1, jmax
    V_vec(0,j,1)  = rho_mms   (ilength,cell_dummy_ghost(0,j,1),cell_dummy_ghost(0,j,2))
    V_vec(0,j,2)  = uvel_mms  (ilength,cell_dummy_ghost(0,j,1),cell_dummy_ghost(0,j,2))
    V_vec(0,j,3)  = vvel_mms  (ilength,cell_dummy_ghost(0,j,1),cell_dummy_ghost(0,j,2))
    V_vec(0,j,4)  = press_mms (ilength,cell_dummy_ghost(0,j,1),cell_dummy_ghost(0,j,2))

    V_vec(-1,j,1) = rho_mms   (ilength,cell_dummy_ghost(-1,j,1),cell_dummy_ghost(-1,j,2))
    V_vec(-1,j,2) = uvel_mms  (ilength,cell_dummy_ghost(-1,j,1),cell_dummy_ghost(-1,j,2))
    V_vec(-1,j,3) = vvel_mms  (ilength,cell_dummy_ghost(-1,j,1),cell_dummy_ghost(-1,j,2))
    V_vec(-1,j,4) = press_mms (ilength,cell_dummy_ghost(-1,j,1),cell_dummy_ghost(-1,j,2))
end do

do j=1, jmax
    V_vec(imax+1,j,1)  = rho_mms   (ilength,cell_dummy_ghost(imax+1,j,1),cell_dummy_ghost(imax+1,j,2))
    V_vec(imax+1,j,2)  = uvel_mms  (ilength,cell_dummy_ghost(imax+1,j,1),cell_dummy_ghost(imax+1,j,2))
    V_vec(imax+1,j,3)  = vvel_mms  (ilength,cell_dummy_ghost(imax+1,j,1),cell_dummy_ghost(imax+1,j,2))
    V_vec(imax+1,j,4)  = press_mms (ilength,cell_dummy_ghost(imax+1,j,1),cell_dummy_ghost(imax+1,j,2))

    V_vec(imax+2,j,1) = rho_mms   (ilength,cell_dummy_ghost(imax+2,j,1),cell_dummy_ghost(imax+2,j,2))
    V_vec(imax+2,j,2) = uvel_mms  (ilength,cell_dummy_ghost(imax+2,j,1),cell_dummy_ghost(imax+2,j,2))
    V_vec(imax+2,j,3) = vvel_mms  (ilength,cell_dummy_ghost(imax+2,j,1),cell_dummy_ghost(imax+2,j,2))
    V_vec(imax+2,j,4) = press_mms (ilength,cell_dummy_ghost(imax+2,j,1),cell_dummy_ghost(imax+2,j,2))
end do

do i=1, imax
    V_vec(i,0,1)   = rho_mms   (ilength,cell_dummy_ghost(i,0,1),cell_dummy_ghost(i,0,2))
    V_vec(i,0,2)   = uvel_mms  (ilength,cell_dummy_ghost(i,0,1),cell_dummy_ghost(i,0,2))
    V_vec(i,0,3)   = vvel_mms  (ilength,cell_dummy_ghost(i,0,1),cell_dummy_ghost(i,0,2))
    V_vec(i,0,4)   = press_mms (ilength,cell_dummy_ghost(i,0,1),cell_dummy_ghost(i,0,2))

    V_vec(i,-1,1)  = rho_mms   (ilength,cell_dummy_ghost(i,-1,1),cell_dummy_ghost(i,-1,2))
    V_vec(i,-1,2)  = uvel_mms  (ilength,cell_dummy_ghost(i,-1,1),cell_dummy_ghost(i,-1,2))
    V_vec(i,-1,3)  = vvel_mms  (ilength,cell_dummy_ghost(i,-1,1),cell_dummy_ghost(i,-1,2))
    V_vec(i,-1,4)  = press_mms (ilength,cell_dummy_ghost(i,-1,1),cell_dummy_ghost(i,-1,2))
end do

do i=1, imax
    V_vec(i,jmax+1,1)  = rho_mms   (ilength,cell_dummy_ghost(i,jmax+1,1),cell_dummy_ghost(i,jmax+1,2))
    V_vec(i,jmax+1,2)  = uvel_mms  (ilength,cell_dummy_ghost(i,jmax+1,1),cell_dummy_ghost(i,jmax+1,2))
    V_vec(i,jmax+1,3)  = vvel_mms  (ilength,cell_dummy_ghost(i,jmax+1,1),cell_dummy_ghost(i,jmax+1,2))
    V_vec(i,jmax+1,4)  = press_mms (ilength,cell_dummy_ghost(i,jmax+1,1),cell_dummy_ghost(i,jmax+1,2))

    V_vec(i,jmax+2,1)  = rho_mms   (ilength,cell_dummy_ghost(i,jmax+2,1),cell_dummy_ghost(i,jmax+2,2))
    V_vec(i,jmax+2,2)  = uvel_mms  (ilength,cell_dummy_ghost(i,jmax+2,1),cell_dummy_ghost(i,jmax+2,2))
    V_vec(i,jmax+2,3)  = vvel_mms  (ilength,cell_dummy_ghost(i,jmax+2,1),cell_dummy_ghost(i,jmax+2,2))
    V_vec(i,jmax+2,4)  = press_mms (ilength,cell_dummy_ghost(i,jmax+2,1),cell_dummy_ghost(i,jmax+2,2))
end do

    V_vec(0,0,:)           = V_vec(1,1,:)
    V_vec(imax+1,jmax+1,:) = V_vec(imax,jmax,:)
    V_vec(imax+1,0,:)      = V_vec(imax,1,:)
    V_vec(0,jmax+1,:)      = V_vec(0,jmax,:)

return
end
  
