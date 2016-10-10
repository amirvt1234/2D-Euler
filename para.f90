
module set_precision
implicit none
save
integer, parameter :: sngl = selected_real_kind( 6,  37)
integer, parameter :: dbl  = selected_real_kind(15, 307)
integer, parameter :: dp   = dbl

end module

module global_parameters
use set_precision, only : dp
implicit none
save
integer, parameter :: nmax = 10000000
integer, parameter :: dim = 4 ! for 1D = 3 2D=4
integer     :: imax, jmax, kmax
real(kind=dp), parameter ::  gamma=1.4_dp, R=287.0_dp
integer :: i, niter, j, k
real(kind=dp) :: dt=0.0_dp, CFL = 0.2_dp, angle
real(kind=dp) :: ilength=1.0_dp, jlength=1.0_dp, dx, dy
integer       :: case_select, globalc=0
real(kind=dp), dimension(3,dim) :: res_norm = 100.0_dp , res_norm0
real(kind=dp) :: epsilon_flux, k_flux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dp), dimension(:,:), allocatable :: x, y   !!! nodes
real(kind=dp), dimension(:,:,:), allocatable :: facex  !_x, facex_y, facex_nx, facex_ny, facex_area
real(kind=dp), dimension(:,:,:), allocatable :: facey  !_x, facey_y, facey_nx, facey_ny, facey_area
real(kind=dp), dimension(:,:,:),allocatable :: flux_tot_x, flux_tot_y 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dp), dimension(:,:,:),allocatable :: cell ! x & y position of cell center ...
real(kind=dp), dimension(:,:,:),allocatable :: cell_dummy_ghost
!real(kind=dp), dimension(:,:),allocatable :: source_mass, source_ymtm, source_xmtm, source_energy
real(kind=dp), dimension(:,:,:),allocatable :: source_mms


real(kind=dp), dimension(:,:,:),allocatable :: U_vec, U_vec_old
real(kind=dp), dimension(:,:,:),allocatable :: V_vec, V_vec_node
real(kind=dp), dimension(:,:,:),allocatable :: res 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dp), dimension(:,:,:),allocatable :: V_vec_mms_c
real(kind=dp), dimension(:,:,:),allocatable :: errord


real(kind=dp) :: T_in=217.0_dp, p_in=12270.0_dp, mach_in=4.0_dp, rho_in, vel_in

end module global_parameters
