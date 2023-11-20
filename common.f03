module common
  use amplitude_mod
  use amplitude_QCD_mod
  implicit none
  real*8,parameter  :: alphaS=0.119d0
  integer :: next,nfin,hel_picked

  type(amplitude) :: amplitudes
  type(amplitude_QCD) :: amps

  ! timing
  real*4 :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.
  real*8 :: amp2,weight
  real*8,dimension(:),allocatable :: amp2_hel
  real(kind=8),dimension(:,:),allocatable,public :: p
  real(kind=8),public :: jac,xbjrk(2)
  

  ! counting events
  integer(kind=4) :: passed=0
  integer(kind=4) :: all_evt=0

  ! setup
  real*8,parameter :: pT_min     = -1d0
  real*8,parameter :: DRjj_min   = -1d0
  real*8,parameter :: eta_max    = -1d0
  real*8,parameter :: sqrt_s_min = 30d0
  
end module common

