module common
  use amplitude_QCD_mod
  implicit none
  real*8,parameter  :: alphaS=0.119d0,alphaEW=7.547E-003
  integer :: next,nfin_glu,hel_picked

  type(amplitude_QCD) :: amps
  type(amplitude_QCD) :: amps_sf

  real*8 :: amp2,weight
  real*8,dimension(:),allocatable :: amp2_hel
  real(kind=8),dimension(:,:),allocatable,public :: p
  real(kind=8),public :: jac,xbjrk(2)

  ! timing
  real*4 :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.

  ! technical
  logical,parameter :: smooth_cuts=.false.
  logical,parameter :: include_pdf=.true.
  logical,parameter :: read_from_file=.false.
  
  ! counting events
  integer(kind=4) :: passed=0
  real(kind=8),dimension(0:5) :: passed_it1(0:5)=0
  integer(kind=4) :: all_evt=0
  integer(kind=4) :: num_error=0

  ! setup cuts. Set to '-1d0' means do not apply any cut on this variable
  real*8,parameter :: pT_min     = 30d0
  real*8,parameter :: DRjj_min   = 0.4d0  ! max allowed value: Drjj_min=1d0
  real*8,parameter :: eta_max    = 5d0
  real*8,parameter :: sqrt_s_min =-1d0
  
end module common

