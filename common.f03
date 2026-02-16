module common
  use particles
  use simple_integrator_mod
  implicit none
  ! Model and process parameters
  type(physics_model) :: phys_model
  real(kind=8),parameter :: alphaS_MZ=0.119d0,alphaEW=0.007546771114d0
  real(kind=8),parameter :: sqrts=14000.d0
  real(kind=8) :: alphaS
  ! scale choice:
  ! 0: ren/fac scale = M_Z
  ! 1: ren/fac scale = H_T
  ! 2: ren/fac scale = H_T/2
  ! 3: ren/fac scale = min(jet-pT)
  ! 4: ren/fac scale = M(non-jet system)
  integer,parameter :: scale_choice=2
  logical,parameter :: use_lhapdf=.true. ! (if false, use internal NNPDF_2.3_nlo_qed)
!!$  character(len=128),parameter :: lhapdfset='NNPDF31_lo_as_0130'
  character(len=128),parameter :: lhapdfset='NNPDF23_nlo_as_0119_qed'

  ! setup cuts. Set to '-1d0' means do not apply any cut on this variable
  ! jets:
  real(kind=8),parameter :: pTj_min      = 30d0
  real(kind=8),parameter :: DRjj_min     = 0.4d0  ! max allowed value: Drjj_min=1d0
  real(kind=8),parameter :: etaj_max     = 6d0
  real(kind=8),parameter :: sqrt_sjj_min = -1d0

  ! photons:
  real(kind=8),parameter :: pTa_min      = 30d0
  real(kind=8),parameter :: DRaa_min     = 0.4d0  ! max allowed value: Draa_min=1d0
  real(kind=8),parameter :: etaa_max     = 6d0
  real(kind=8),parameter :: sqrt_saa_min = -1d0

  ! leptons:
  real(kind=8),parameter :: pTl_min      = -1d0
  real(kind=8),parameter :: DRll_min     = -1d0  ! max allowed value: Drll_min=1d0
  real(kind=8),parameter :: etal_max     = -1d0
  real(kind=8),parameter :: sqrt_sll_min = -1d0

  ! jets+photons:
  real(kind=8),parameter :: DRja_min     = 0.4d0  ! max allowed value: Drja_min=1d0
  real(kind=8),parameter :: sqrt_sja_min = -1d0

  ! jets+leptons:
  real(kind=8),parameter :: DRjl_min     = -1d0  ! max allowed value: Drjl_min=1d0
  real(kind=8),parameter :: sqrt_sjl_min = -1d0

  ! leptons+photons:
  real(kind=8),parameter :: DRla_min     = 0.4d0  ! max allowed value: Drla_min=1d0
  real(kind=8),parameter :: sqrt_sla_min = -1d0
  
  ! technical
  logical,parameter :: include_pdf=.true.
  logical,parameter :: use_colour_singlet_multichannel=.true.
  logical,parameter :: reduce_to_unique_matrix_elements=.true.
  logical,parameter :: decompose_same_flavour_into_two_diff_flavour=.true.
  logical,parameter :: use_cross_process_optimisation_of_currents=.true.
  logical :: keep_processes_separate
  real(kind=4) :: tBefore,tAfter,tTot_A,tTot_B
  real(kind=4) :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.,t_Proc_init=0.
  type(integrator) :: simple_integrator
  real(kind=8) :: scale_ren,scale_fac,scale_shower
end module common

