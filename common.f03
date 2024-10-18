module common
  use amplitude_QCD_mod
  implicit none
  ! coupling constants
  real(kind=8),parameter  :: alphaS=0.119d0,alphaEW=0.00754677114d0

  ! timing
  real(kind=4) :: tBefore,tAfter,tTot_A,tTot_B
  real(kind=4) :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.

  ! technical
  logical,parameter :: smooth_cuts=.false.
  logical,parameter :: include_pdf=.true.
  
  ! setup cuts. Set to '-1d0' means do not apply any cut on this variable
  real(kind=8),parameter :: pT_min     = 30d0
  real(kind=8),parameter :: DRjj_min   = 0.4d0  ! max allowed value: Drjj_min=1d0
  real(kind=8),parameter :: eta_max    = 6d0
  real(kind=8),parameter :: sqrt_s_min = -1d0
  
end module common

