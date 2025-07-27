module common
  use amplitude_QCD_mod
  use phase_space_base
  use particles
  implicit none
  ! Model and process parameters
  type(physics_model) :: phys_model
  real(kind=8),parameter  :: alphaS=0.119d0,alphaEW=0.007546771114d0
  real(kind=8) :: sqrts

  ! timing
  real(kind=4) :: tBefore,tAfter,tTot_A,tTot_B
  real(kind=4) :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.

  ! technical
  logical,parameter :: include_pdf=.true.
  logical,parameter :: use_colour_singlet_multichannel=.true.
  type :: multichan_info
     ! if adding variables here, also update the finalize_multichan_info subroutine
     integer,dimension(:,:),allocatable :: channels,unique_channelgroup_list
     integer,dimension(:),allocatable :: unique_channel_list,map_proc_to_channelgroup,number_of_channels
     integer :: max_channels,n_unique_channels,n_unique_channelgroups
  end type multichan_info

  type phase_space_order_group
     ! if adding variables here, also update the finalize_phase_space_order_group subroutine
     type(amplitude_QCD) :: amps
     class(phase_space_type),allocatable :: phase_space
     type(multichan_info) :: multichan
     integer,dimension(:,:),allocatable :: processes,color_orders
     integer,dimension(:),allocatable :: iden_iproc,phase_space_orders
     integer :: nproc
     real(kind=8),dimension(:,:),allocatable :: val_procs,idenCOandMAPfactor
     integer,dimension(:,:,:),allocatable :: iden_processes
     integer(kind=4),dimension(:,:),allocatable :: spin
     integer(kind=8),dimension(:),allocatable :: iden
     logical,dimension(-6:7,2) :: ipdgs
     integer(kind=4) :: nhel,next
     integer,dimension(:),allocatable :: col_fac
     real(kind=8),dimension(:),allocatable :: amp2,amp2_hel
     integer(kind=4),dimension(:),allocatable :: hel,hel_fac
     integer(kind=4) :: passed=0,all_evt=0
     integer,dimension(:),allocatable :: include_hel
     ! cuts
     double precision,dimension(:),allocatable :: pT_min,eta_max
     double precision,dimension(:,:),allocatable :: DR_min,sqrt_s_min
  end type phase_space_order_group

contains
  subroutine finalize_multichan_info(mi)
    type(multichan_info),intent(inout) :: mi
    if (allocated(mi%channels)) deallocate(mi%channels)
    if (allocated(mi%unique_channelgroup_list)) deallocate(mi%unique_channelgroup_list)
    if (allocated(mi%unique_channel_list)) deallocate(mi%unique_channel_list)
    if (allocated(mi%map_proc_to_channelgroup)) deallocate(mi%map_proc_to_channelgroup)
    if (allocated(mi%number_of_channels)) deallocate(mi%number_of_channels)
  end subroutine finalize_multichan_info

  subroutine finalize_phase_space_order_group(pgl)
    type(phase_space_order_group),intent(inout) :: pgl
    call finalize_amplitude_QCD(pgl%amps)
    if(allocated(pgl%phase_space)) then
       call pgl%phase_space%cleanup()
       deallocate(pgl%phase_space)
    endif
    call finalize_multichan_info(pgl%multichan)
    if (allocated(pgl%processes)) deallocate(pgl%processes)
    if (allocated(pgl%color_orders)) deallocate(pgl%color_orders)
    if (allocated(pgl%iden_iproc)) deallocate(pgl%iden_iproc)
    if (allocated(pgl%phase_space_orders)) deallocate(pgl%phase_space_orders)
    if (allocated(pgl%val_procs)) deallocate(pgl%val_procs)
    if (allocated(pgl%idenCOandMAPfactor)) deallocate(pgl%idenCOandMAPfactor)
    if (allocated(pgl%iden_processes)) deallocate(pgl%iden_processes)
    if (allocated(pgl%spin)) deallocate(pgl%spin)
    if (allocated(pgl%iden)) deallocate(pgl%iden)
    if (allocated(pgl%col_fac)) deallocate(pgl%col_fac)
    if (allocated(pgl%amp2)) deallocate(pgl%amp2)
    if (allocated(pgl%amp2_hel)) deallocate(pgl%amp2_hel)
    if (allocated(pgl%hel)) deallocate(pgl%hel)
    if (allocated(pgl%hel_fac)) deallocate(pgl%hel_fac)
    if (allocated(pgl%include_hel)) deallocate(pgl%include_hel)
    if (allocated(pgl%pT_min)) deallocate(pgl%pT_min)
    if (allocated(pgl%eta_max)) deallocate(pgl%eta_max)
    if (allocated(pgl%DR_min)) deallocate(pgl%DR_min)
    if (allocated(pgl%sqrt_s_min)) deallocate(pgl%sqrt_s_min)
  end subroutine finalize_phase_space_order_group

end module common

