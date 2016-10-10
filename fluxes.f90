module fluxes
use global_parameters, only:dp, niter, i, j, gamma, jmax, imax, dim, facex, facey, res_norm, globalc, epsilon_flux, k_flux
use set_precision, only : dp
use some_functions, only: VtoU


implicit none
!save




contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine flux_param(V_vec_d,flux_tot_xd,flux_tot_yd)
implicit none
real(kind=dp), dimension(-1:imax+2,-1:jmax+2,1:dim), intent(in)  :: V_vec_d
real(kind=dp), dimension(1:imax+1,1:jmax+1,1:dim),  intent(out) :: flux_tot_xd,flux_tot_yd
real(kind=dp), dimension(-1:imax+2,-1:jmax+2,1:dim,2)   ::  V_vec_tilda_x, V_vec_tilda_y, psi_porm_x, psi_porm_y
real(kind=dp), dimension(-1:imax+2,-1:jmax+2,2) :: speed_s_lorr_x, speed_s_lorr_y 

psi_porm_x= 0.0
psi_porm_y= 0.0
do j=2,jmax
    do i=2,jmax
        psi_porm_x (i,j,:,:) = psi_calc(V_vec_d(i-2:i+1,j,:))
        psi_porm_y (i,j,:,:) = psi_calc(V_vec_d(i,j-2:j+1,:))
    end do
end do
!psi_porm_x= 1.0_dp
!psi_porm_y= 1.0_dp



do j=1,jmax
    do i=1,imax+1
        V_vec_tilda_x(i,j,:,:) = V_vec_tilda_calc(V_vec_d(i-2:i+1,j,:),psi_porm_x(i-1:i+1,j,:,:))
        if (V_vec_tilda_x(i,j,1,1) < 0.1 ) V_vec_tilda_x(i,j,1,1)=0.1
        if (V_vec_tilda_x(i,j,4,1) < 10.0 ) V_vec_tilda_x(i,j,4,1)=10.0
        if (V_vec_tilda_x(i,j,1,2) < 0.1 ) V_vec_tilda_x(i,j,1,2)=0.1
        if (V_vec_tilda_x(i,j,4,2) < 10.0 ) V_vec_tilda_x(i,j,4,2)=10.0
        speed_s_lorr_x(i,j,1)  = dsqrt(dabs(gamma*V_vec_tilda_x(i,j,4,1)/V_vec_tilda_x(i,j,1,1)))
        speed_s_lorr_x(i,j,2)  = dsqrt(dabs(gamma*V_vec_tilda_x(i,j,4,2)/V_vec_tilda_x(i,j,1,2)))
    
        flux_tot_xd(i,j,:) = flux_calc(speed_s_lorr_x(i,j,:), V_vec_tilda_x(i,j,1,:) &
                                      , V_vec_tilda_x(i,j,2,:), V_vec_tilda_x(i,j,3,:), V_vec_tilda_x(i,j,4,:)  &
                                      , facex(i,j,3), facex(i,j,4))
    end do
end do

do j=1,jmax+1
    do i=1,imax
        V_vec_tilda_y(i,j,:,:) = V_vec_tilda_calc(V_vec_d(i,j-2:j+1,:),psi_porm_x(i,j-1:j+1,:,:))
        if (V_vec_tilda_y(i,j,1,1) < 0.1 ) V_vec_tilda_y(i,j,1,1)=0.1
        if (V_vec_tilda_y(i,j,4,1) < 10.0 ) V_vec_tilda_y(i,j,4,1)=10.0
        if (V_vec_tilda_y(i,j,1,2) < 0.1 ) V_vec_tilda_y(i,j,1,2)=0.1
        if (V_vec_tilda_y(i,j,4,2) < 10.0 ) V_vec_tilda_y(i,j,4,2)=10.0
        speed_s_lorr_y(i,j,1)  = dsqrt(dabs(gamma*V_vec_tilda_y(i,j,4,1)/V_vec_tilda_y(i,j,1,1)))
        speed_s_lorr_y(i,j,2)  = dsqrt(dabs(gamma*V_vec_tilda_y(i,j,4,2)/V_vec_tilda_y(i,j,1,2)))
    
        flux_tot_yd(i,j,:) = flux_calc(speed_s_lorr_y(i,j,:), V_vec_tilda_y(i,j,1,:) &
                                      , V_vec_tilda_y(i,j,2,:), V_vec_tilda_y(i,j,3,:), V_vec_tilda_y(i,j,4,:) &
                                      , facey(i,j,3), facey(i,j,4))
    end do
end do

end subroutine flux_param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function V_vec_tilda_calc(V_vec_dd,psi_porm_d)
implicit none
real(kind=dp), dimension(4,dim), intent(in) ::V_vec_dd   ! loop counter
real(kind=dp), dimension(3,dim,2), intent(in)   ::   psi_porm_d
real(kind=dp), dimension(1:dim,2)   ::   V_vec_tilda_calc


if (niter == 1200000) then
k_flux = -1.0_dp
epsilon_flux = 1.0_dp
endif

!k_flux = 0.0_dp
!epsilon_flux = 0.0_dp
!psi_porm_d = 1.0_dp
 

V_vec_tilda_calc(:,1) = V_vec_dd(2,:) +0.25*epsilon_flux*                             &
                    ((1-k_flux)*psi_porm_d (1,:,1)*(V_vec_dd(2,:)-V_vec_dd(1,:)) + &
                       (1+k_flux)*psi_porm_d (2,:,2)*(V_vec_dd(3,:)-V_vec_dd(2,:))      ) 

V_vec_tilda_calc(:,2) = V_vec_dd(3,:) -0.25*epsilon_flux*                             &
                    ((1-k_flux)*psi_porm_d (3,:,2)*(V_vec_dd(4,:)-V_vec_dd(3,:)) + &
                       (1+k_flux)*psi_porm_d (2,:,1)*(V_vec_dd(3,:)-V_vec_dd(2,:))      )

end function V_vec_tilda_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function psi_calc(V_vec_psi)
implicit none
real(kind=dp), dimension(4,dim), intent(in) ::V_vec_psi   ! loop counter
real(kind=dp), dimension(1:dim,2)   ::   r_porm
real(kind=dp), dimension(1:dim,2)   ::   psi_calc
integer :: ii
do ii=1,dim
    r_porm   (ii,1) = (V_vec_psi(4,ii)-V_vec_psi(3,ii))/(sign(1.0_dp,V_vec_psi(3,ii)-V_vec_psi(2,ii))* &
		       max(dabs(V_vec_psi(3,ii)-V_vec_psi(2,ii)),(10._dp**(-3))))
    r_porm   (ii,2) = (V_vec_psi(2,ii)-V_vec_psi(1,ii))/(sign(1.0_dp,V_vec_psi(3,ii)-V_vec_psi(2,ii))* & 
		       max(dabs(V_vec_psi(3,ii)-V_vec_psi(2,ii)),(10._dp**(-3))))


    psi_calc (ii,1) = (r_porm(ii,1)+dabs(r_porm(ii,1)))/(1+dabs(r_porm(ii,1)))
    psi_calc (ii,2) = (r_porm(ii,2)+dabs(r_porm(ii,2)))/(1+dabs(r_porm(ii,2)))
end do
end function psi_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function flux_calc(speed_sd, rho_lorr, u_lorr, v_lorr, p_lorr, nx_f, ny_f)

implicit none
real(kind=dp), intent(in) :: nx_f, ny_f  
real(kind=dp), dimension(2), intent(in) :: speed_sd  ! loop counter
real(kind=dp), dimension(2), intent(in) :: rho_lorr, u_lorr, v_lorr, p_lorr
real(kind=dp), dimension(dim) :: flux_calc

!! Begin: Variables for FVS scheme
   real(kind=dp), dimension(dim) :: flux_conv,flux_pres
   real(kind=dp), dimension(2)   :: M_lorr, beta
   real(kind=dp), dimension(2) :: M_porm, alpha_porm, c_porm
   real(kind=dp), dimension(2) :: d_porm, pdb_porm
!! End
!! Begin: variables for FDS Scheme
   real(kind=dp)  ::   rhodbar, R_fds, udbar, vdbar, htdbar, a2dbar, adbar
   real(kind=dp)  ::   lambda1dbar,lambda2dbar,lambda3dbar, lambda4dbar, lambda_param
   real(kind=dp), dimension(dim) ::   deltaw
   real(kind=dp), dimension(dim) ::   r1dbar, r2dbar, r3dbar, r4dbar
   real(kind=dp), dimension(dim) ::   flux_l, flux_r
!! End
   real(kind=dp), dimension(dim,2) :: U_vec_lorr, V_vec_lorr

   character(len=3) :: flux_selector



flux_selector ='FVS'
!print*, flux_selector 
!pause
select case (flux_selector)
case ('FVS')

    M_lorr     (1)   = (u_lorr(1)*nx_f+v_lorr(1)*ny_f)/speed_sd(1)
    M_lorr     (2)   = (u_lorr(2)*nx_f+v_lorr(2)*ny_f)/speed_sd(2)
    M_porm     (1) =  0.25*(M_lorr(1)+1)**2 !M_porm(1) = M^{plus}(i)
    M_porm     (2) = -0.25*(M_lorr(2)-1)**2 !M_porm(2) = M^{minus}(i)
    beta       (1)   = -max(0,1-int(abs(M_lorr(1))))
    beta       (2)   = -max(0,1-int(abs(M_lorr(2))))
    alpha_porm (1) =  0.5*(1._dp+sign(1._dp,M_lorr(1)))    !alpha_porm(1) = alpha^{plus}(i)
    alpha_porm (2) =  0.5*(1._dp-sign(1._dp,M_lorr(2)))    !alpha_porm(2) = alpha^{minus}(i)
    c_porm     (1) =  alpha_porm (1)*(1+beta(1))*M_lorr(1)-beta(1)*M_porm(1) !...
    c_porm     (2) =  alpha_porm (2)*(1+beta(2))*M_lorr(2)-beta(2)*M_porm(2) !...
    pdb_porm   (1) =  M_porm(1)*(-M_lorr(1)+2)  ! p^{dbar}
    pdb_porm   (2) =  M_porm(2)*(-M_lorr(2)-2)  ! p^{dbar}
    d_porm     (1) =  alpha_porm (1)*(1+beta(1))-beta(1)*pdb_porm(1)
    d_porm     (2) =  alpha_porm (2)*(1+beta(2))-beta(2)*pdb_porm(2)

    flux_pres  (1) =  0.0_dp ! zero
    flux_pres  (2) =  d_porm(1)*p_lorr(1)*nx_f+d_porm(2)*p_lorr(2)*nx_f !
    flux_pres  (3) =  d_porm(1)*p_lorr(1)*ny_f+d_porm(2)*p_lorr(2)*ny_f !
    flux_pres  (4) =  0.0_dp ! zero

    flux_conv  (1) =  rho_lorr(1)*speed_sd(1)*c_porm(1) * 1._dp  +  &
                              rho_lorr(2)*speed_sd(2)*c_porm(2) * 1._dp     !
    flux_conv  (2) =  rho_lorr(1)*speed_sd(1)*c_porm(1) * u_lorr(1)  +  &
                              rho_lorr(2)*speed_sd(2)*c_porm(2) * u_lorr(2)       !
    flux_conv  (3) =  rho_lorr(1)*speed_sd(1)*c_porm(1) * v_lorr(1)  +  &
                              rho_lorr(2)*speed_sd(2)*c_porm(2) * v_lorr(2)       !
    flux_conv  (4) =  rho_lorr(1)*speed_sd(1)*c_porm(1) * (p_lorr(1)/rho_lorr(1)*gamma/(gamma-1)+ &
		             0.5*(u_lorr(1)**2+v_lorr(1)**2)) + &
                      rho_lorr(2)*speed_sd(2)*c_porm(2) * (p_lorr(2)/rho_lorr(2)*gamma/(gamma-1)+ & 
		             0.5*(u_lorr(2)**2+v_lorr(2)**2))       !

    flux_calc   (1) = flux_conv(1)+flux_pres(1)
    flux_calc   (2) = flux_conv(2)+flux_pres(2)
    flux_calc   (3) = flux_conv(3)+flux_pres(3)
    flux_calc   (4) = flux_conv(4)+flux_pres(4)
    
case ('FDS')
    V_vec_lorr (1,:) = rho_lorr
    V_vec_lorr (2,:) = u_lorr
    V_vec_lorr (3,:) = v_lorr
    V_vec_lorr (4,:) = p_lorr
    U_vec_lorr (:,1) = VtoU(V_vec_lorr(:,1))
    U_vec_lorr (:,2) = VtoU(V_vec_lorr(:,2))

    R_fds     = dsqrt(dabs(rho_lorr(2)/rho_lorr(1)))
    rhodbar   = R_fds*rho_lorr(1)
    udbar     = (R_fds*u_lorr(2)+u_lorr(1))/(R_fds+1)
    vdbar     = (R_fds*v_lorr(2)+v_lorr(1))/(R_fds+1)
    htdbar    = (R_fds*(p_lorr(2)/rho_lorr(2)*gamma/(gamma-1)+0.5*(u_lorr(2)**2+v_lorr(2)**2))+             &
                       (p_lorr(1)/rho_lorr(1)*gamma/(gamma-1)+0.5*(u_lorr(1)**2+v_lorr(1)**2)))/(R_fds+1)
    a2dbar    = dabs((gamma-1.)*(htdbar-0.5*(udbar**2+vdbar**2))) 
    adbar     = dsqrt(a2dbar)

    lambda_param = 0.1 ! where    0.0<lambda_param<0.1
    lambda1dbar  = udbar*nx_f+vdbar*ny_f
    if (abs(lambda1dbar ) < 2*lambda_param*adbar) lambda1dbar=lambda1dbar**2/(4*lambda_param*a2dbar)+lambda_param*adbar
    lambda2dbar  = udbar*nx_f+vdbar*ny_f
    if (abs(lambda2dbar ) < 2*lambda_param*adbar) lambda2dbar=lambda2dbar**2/(4*lambda_param*a2dbar)+lambda_param*adbar
    lambda3dbar  = udbar*nx_f+vdbar*ny_f+adbar
    if (abs(lambda3dbar ) < 2*lambda_param*adbar) lambda3dbar=lambda3dbar**2/(4*lambda_param*a2dbar)+lambda_param*adbar
    lambda4dbar  = udbar*nx_f+vdbar*ny_f-adbar
    if (abs(lambda4dbar ) < 2*lambda_param*adbar) lambda4dbar=lambda4dbar**2/(4*lambda_param*a2dbar)+lambda_param*adbar

    r1dbar(1)  = 1._dp
    r1dbar(2)  = udbar
    r1dbar(3)  = vdbar
    r1dbar(4)  = (udbar**2+vdbar**2)/2.

    r2dbar(1)  = 0.0_dp
    r2dbar(2)  = ny_f*rhodbar
    r2dbar(3)  = -nx_f*rhodbar
    r2dbar(4)  = rhodbar*(ny_f*udbar-nx_f*vdbar)

    r3dbar(1)  = (rhodbar/(2*adbar))*1.0
    r3dbar(2)  = (rhodbar/(2*adbar))*(udbar+nx_f*adbar)
    r3dbar(3)  = (rhodbar/(2*adbar))*(vdbar+ny_f*adbar)
    r3dbar(4)  = (rhodbar/(2*adbar))*(htdbar+adbar*udbar*nx_f+adbar*vdbar*ny_f)

    r4dbar(1) = (-rhodbar/(2*adbar))*1.0
    r4dbar(2) = (-rhodbar/(2*adbar))*(udbar-nx_f*adbar)
    r4dbar(3) = (-rhodbar/(2*adbar))*(vdbar-ny_f*adbar)
    r4dbar(4) = (-rhodbar/(2*adbar))*(htdbar-adbar*udbar*nx_f-adbar*vdbar*ny_f)

    deltaw (1)  = (rho_lorr(2)-rho_lorr(1)) - (p_lorr(2)-p_lorr(1))/a2dbar
    deltaw (2)  =  ny_f*(u_lorr(2)-u_lorr(1))-nx_f*(v_lorr(2)-v_lorr(1))
    deltaw (3)  = +nx_f*(u_lorr(2)-u_lorr(1))+ny_f*(v_lorr(2)-v_lorr(1)) + (p_lorr(2)-p_lorr(1))/(adbar*rhodbar)
    deltaw (4)  = +nx_f*(u_lorr(2)-u_lorr(1))+ny_f*(v_lorr(2)-v_lorr(1)) - (p_lorr(2)-p_lorr(1))/(adbar*rhodbar)

    flux_l = flux(V_vec_lorr(:,1),nx_f,ny_f)
    flux_r = flux(V_vec_lorr(:,2),nx_f,ny_f)
       
    flux_calc   (1) = 0.5*(flux_l(1)+flux_r(1))  -0.5*(abs(lambda1dbar)*deltaw(1)*r1dbar(1)  + &
                                                       abs(lambda2dbar)*deltaw(2)*r2dbar(1)  + &
                                                       abs(lambda3dbar)*deltaw(3)*r3dbar(1)  + & 
                                                       abs(lambda4dbar)*deltaw(4)*r4dbar(1)     )   

    flux_calc   (2) = 0.5*(flux_l(2)+flux_r(2))  -0.5*(abs(lambda1dbar)*deltaw(1)*r1dbar(2)  + &
                                                       abs(lambda2dbar)*deltaw(2)*r2dbar(2)  + &
                                                       abs(lambda3dbar)*deltaw(3)*r3dbar(2)  + & 
                                                       abs(lambda4dbar)*deltaw(4)*r4dbar(2)     )

    flux_calc   (3) = 0.5*(flux_l(3)+flux_r(3))  -0.5*(abs(lambda1dbar)*deltaw(1)*r1dbar(3)  + &
                                                       abs(lambda2dbar)*deltaw(2)*r2dbar(3)  + &
                                                       abs(lambda3dbar)*deltaw(3)*r3dbar(3)  + & 
                                                       abs(lambda4dbar)*deltaw(4)*r4dbar(3)     )

    flux_calc   (4) = 0.5*(flux_l(4)+flux_r(4))  -0.5*(abs(lambda1dbar)*deltaw(1)*r1dbar(4)  + &
                                                       abs(lambda2dbar)*deltaw(2)*r2dbar(4)  + &
                                                       abs(lambda3dbar)*deltaw(3)*r3dbar(4)  + & 
                                                       abs(lambda4dbar)*deltaw(4)*r4dbar(4)     )

end select

!!return
end function flux_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function flux (f_V, nx_ff, ny_ff)
use set_precision, only : dp
use global_parameters, only: gamma, dim
implicit none
real(kind=dp), intent(in) :: nx_ff, ny_ff
!integer, intent(in) :: ii  ! loop counter
real(kind=dp), dimension(dim), intent(in) :: f_V
real(kind=dp), dimension(dim) :: flux
flux(1) = f_V(1)*(f_V(2)*nx_ff+f_V(3)*ny_ff)          ! calculates fluxes
flux(2) = f_V(1)*f_V(2)*(f_V(2)*nx_ff+f_V(3)*ny_ff)+f_V(4)*nx_ff   !       (gamma-1.)*f_U(3)+0.5*(3.-gamma)*f_U(2)**2/f_U(1) 
flux(3) = f_V(1)*f_V(3)*(f_V(2)*nx_ff+f_V(3)*ny_ff)+f_V(4)*ny_ff
flux(4) = f_V(1)*((gamma*f_V(4)/(f_V(1)*(gamma-1))+0.5*(f_V(2)**2+f_V(3)**2)))*(f_V(2)*nx_ff+f_V(3)*ny_ff)
end function flux

end module fluxes
