module common
  use particles, only: physics_model
  use run_parameters
  use simple_integrator_mod
  implicit none
  ! Model and process parameters
  type(physics_model) :: phys_model
  real(kind=8) :: alphaS

  ! technical
  logical,parameter :: use_colour_singlet_multichannel=.true.
  logical,parameter :: reduce_to_unique_matrix_elements=.true.
  logical,parameter :: decompose_same_flavour_amplitudes=.true.
  logical,parameter :: use_cross_process_optimisation_of_currents=.true.
  logical :: keep_processes_separate=.true.
  integer,parameter :: timing_none=0,timing_basic=1,timing_detailed=2
  integer :: timing_mode=timing_basic,timing_sample=100
  real(kind=8) :: tBefore,tAfter,tTot_A,tTot_B
  real(kind=8) :: t_Initialise=0d0,t_PS_init=0d0,t_Amp_init=0d0,t_PS=0d0,t_Amp=0d0,t_all=0d0,t_mat=0d0,&
       t_Proc_init=0d0,t_Amp_opt=0d0,t_weight=0d0,t_lib_check=0d0,t_other=0d0,&
       t_Int_init=0d0,t_Int_get=0d0,t_Int_fill=0d0,t_Evt_write=0d0,&
       t_Evt_wgt_assign=0d0,t_Evt_wgt_update=0d0,t_Int_loop=0d0,t_Finalise=0d0
  type(integrator) :: simple_integrator
  real(kind=8) :: scale_ren,scale_fac,scale_shower
end module common
