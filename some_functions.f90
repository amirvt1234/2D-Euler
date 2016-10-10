module some_functions

implicit none


contains

function VtoU (V_vec_d)
use global_parameters, only: gamma, dim
use set_precision, only : dp
implicit none
real(kind=dp), dimension(dim), intent(in)   :: V_vec_d
real(kind=dp), dimension(dim) :: VtoU
VtoU(1)   = V_vec_d(1)         
VtoU(2)   = V_vec_d(1)*V_vec_d(2)
VtoU(3)   = V_vec_d(1)*V_vec_d(3)
VtoU(4)   = V_vec_d(4)/(gamma-1.)+V_vec_d(1)*(V_vec_d(2)**2+V_vec_d(3)**2)/2.
end function VtoU

function UtoV (U_vec_d)
use set_precision, only : dp
use global_parameters, only: gamma, dim
implicit none
real(kind=dp), intent(in), dimension(dim)   :: U_vec_d  ! loop counter
real(kind=dp), dimension(dim) :: UtoV
UtoV(1)     = U_vec_d(1)                       
UtoV(2)     = U_vec_d(2)/U_vec_d(1)
UtoV(3)     = U_vec_d(3)/U_vec_d(1) 
Utov(4)     = ((gamma-1.)*U_vec_d(4)-0.5*(gamma-1.)*(U_vec_d(2)**2+U_vec_d(3)**2)/U_vec_d(1))
end function UtoV

function f_norm(ff, n, m)
use global_parameters, only: dp,dim
implicit none
real(kind=dp), dimension(n,m,dim) :: ff
real(kind=dp), dimension(3,dim) :: f_norm
integer :: i, n, m
do i=1,dim
    f_norm(1,i)  = sum(dabs(ff(:,:,i)))/(n*m)
    f_norm(2,i)  = dsqrt(sum(ff(:,:,i)**2)/(n*m))
    f_norm(3,i)  = maxval(dabs(ff(:,:,i)))
end do
end function f_norm

function dfdu(f,v,nx,ny)
use global_parameters, only: dp,dim, gamma
implicit none
real(kind=dp), dimension(1:dim) :: f, v
real(kind=dp), dimension(dim,dim) :: dfdu
real(kind=dp) :: nx,ny, nvu, vel

nvu = nx*v(1)+ny*v(2)
vel = v(2)**2+v(3)**2

dfdu (1,1)  = nvu 
dfdu (1,2)  = nx*v(1)
dfdu (1,3)  = ny*v(1)
dfdu (1,4)  = 0.0_dp

dfdu (2,1)  = nvu*v(2) 
dfdu (2,2)  = v(1)*(2*nx*v(2)+ny*v(3))
dfdu (2,3)  = ny*v(1)*v(2)
dfdu (2,4)  = nx

dfdu (3,1)  = nvu*v(3) 
dfdu (3,2)  = nx*v(1)*v(3)
dfdu (3,3)  = v(1)*(nx*v(2)+2*ny*v(3))
dfdu (3,4)  = ny

dfdu (4,1)  = 0.5*vel*nvu
dfdu (4,2)  = v(1)*v(2)*nvu+nx*v(1)*(0.5*vel+(gamma*v(4)/(v(1)*(gamma-1))))
dfdu (4,3)  = v(1)*v(3)*nvu+ny*v(1)*(0.5*vel+(gamma*v(4)/(v(1)*(gamma-1))))
dfdu (4,4)  = gamma*nvu/(gamma-1)
 
end function dfdu


end module some_functions
