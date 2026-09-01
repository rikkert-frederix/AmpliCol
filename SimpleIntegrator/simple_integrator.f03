!===============================================================================
! SimpleIntegrator
!===============================================================================
!
! Purpose
! -------
! This module provides a compact adaptive Monte Carlo integrator with optional
! unweighted event generation.  It is intended for matrix-element, phase-space,
! or similar event-generation programs where the caller owns the physics
! function and the integrator owns:
!
!   * adaptive one-dimensional grids for every integration variable,
!   * distribution of requested points across channels and integrals,
!   * accumulated signed and absolute integral estimates,
!   * candidate-event storage and final unweighted-event weights.
!
! Public module and dependencies
! ------------------------------
! Use the public type through
!
!     use simple_integrator_mod
!     type(integrator) :: integ
!
! The module expects the helper modules in helper_modules.f03 and the random
! number function ran2() from ranmar.f (or a caller-provided compatible
! function returning a double precision number in (0,1)).  Several progress
! messages are written to stdout and to Fortran unit 99.  Standalone programs
! should open unit 99 before using the integrator if they want to control the
! log destination; otherwise many compilers create a default file for that unit.
!
! Concepts
! --------
! A "channel" is a separately adapted sampling map, usually a phase-space
! channel.  Each channel can contain one or more "integrals", for example
! subprocesses sharing the same channel.  The public channel and integral labels
! returned by get_points are one-based indices.
!
! Each point has ndim adapted coordinates and ndim_extra flat random
! coordinates.  Only the first ndim coordinates are included in adaptive grid
! weights and in compute_wgt_from_x.  The caller may use ndim_extra for random
! choices that should not affect grid adaptation.
!
! Public workflow
! ---------------
! 1. Initialise once:
!
!        call integ%init(nchannel, ndim, ndim_extra, nintegral, &
!             nevts_unw_req, niters, adaptation_classes=classes)
!
!    ndim, ndim_extra, and nintegral are arrays of length nchannel.
!    nevts_unw_req is the requested number of unweighted events.  niters is the
!    maximum number of adaptation/generation iterations.  The optional
!    adaptation_classes(:,channel) assigns each integral to a positive,
!    contiguous class.  Integrals in one class share grids; different classes
!    have independent grids while retaining independent integral estimates.
!
! 2. Repeatedly request points, evaluate the caller's integrand, and return the
!    values:
!
!        done = .false.
!        do while (.not. done)
!           call integ%get_points(npoints, ichan, iint)
!           do ip = 1, npoints
!              ! integ%x(:,ip) contains the random point.
!              ! integ%wgt(ip) is the grid Jacobian/volume factor.
!              ! Compute f_abs(ip) >= 0 and f(ip), including integ%wgt(ip).
!           end do
!           call integ%fill_points(npoints, f_abs, f, to_write, done, accepted, &
!                external_converged=my_external_test)
!           ! If to_write(ip) is true, write/store the corresponding event now.
!        end do
!
!    fill_points may be called with fewer points than were returned by
!    get_points; unused trailing points are discarded.  A new get_points call
!    must not be made until the previous batch has been returned with
!    fill_points.  If external_converged is present and false, satisfying the
!    requested statistical accuracy does not terminate the integration; the
!    next accuracy-only budget is doubled and its grids are refined.  The
!    iteration cap remains authoritative, so callers should report when their
!    external criterion is still false at that cap.
!
! 3. After done is true, retrieve final event weights:
!
!        call integ%assign_evnt_wgts(wgts)
!
!    wgts has shape (3, number_of_written_candidate_events).  Column i contains
!    the nominal event weight, the adjusted event weight including any
!    overweight correction, and the overweight excess.  Events rejected by the
!    final unweighting have zero weights.
!
! Integrand convention
! --------------------
! f_abs is the non-negative envelope used for maximum-weight estimates,
! unweighting, and the absolute-envelope accuracy target.  f is the signed
! contribution to the physical integral.  Event-generation grids adapt to
! f_abs; accuracy-only grids adapt to the mean of abs(f), after local
! cancellations.  The optional accepted array records
! points that pass the caller's cuts; accepted zero-weight points count towards
! the per-iteration statistics, while rejected points do not.
!
! compute_wgt_from_x can be used when an external multichannel combination needs
! the current adaptive-grid weight for a point already known in a specific
! channel.
!
! Tunable internal parameters
! ---------------------------
! The parameters below the type declarations control minimum statistics, grid
! sizes, event-generation startup, and allowed overweight fraction.  They are
! compile-time parameters in this standalone version.
!
! Limitations
! -----------
! This implementation is serial and keeps candidate events in memory until final
! weighting.  It does not currently read or write grids; read_all_grids and
! write_all_grids are placeholders.
!
module simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use random_number_interface, only: ran2
  implicit none
  private
  ! One adaptive sampling channel.  A channel owns one grid per adapted
  ! dimension and adaptation class.  Integral estimates may share a class.
  type :: channel
     integer :: ndim,nintegral,naux,current_integral,current_iter&
          &,number,max_iters,nevts_unw_req,ndim_extra,nadaptation,warmup_iterations
     integer(kind=8) :: npoints,npoints_iter
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,unc_iter
     real(kind=8),allocatable,dimension(:) :: aux_res,aux_unc,aux_res_iter,aux_unc_iter
     real(kind=8) :: overweight
     logical :: done,evgen_done,integration_only
     integer,allocatable,dimension(:) :: adaptation_class
     type(grid),allocatable,dimension(:,:,:) :: grids
     type(integral),allocatable,dimension(:) :: integrals
   contains
     procedure,private :: init => channel_init
     procedure,private :: add_point => channel_add_point
     procedure,private :: get_point => channel_get_point
     procedure,private :: update_result_iter => channel_update_result_iter
     procedure,private :: combine_iters => channel_combine_iters
     procedure,private :: print_result_iter => channel_print_result_iter
     procedure,private :: print_combined_result => channel_print_combined_result
     procedure,private :: init_next_iter => channel_init_next_iter
     procedure,private :: check_gen_evnts => channel_check_gen_evnts
     procedure,private :: update_grids => channel_update_grids
     procedure,private :: update_nevts_unw_req => channel_update_nevts_unw_req
     procedure,private :: recompute_wgt_from_x
  end type channel
  ! One physical integral inside a channel.  It accumulates iteration estimates
  ! and stores candidate events until the final unweighting decision.
  type :: integral
     real(kind=8) :: max_value,overweight
     real(kind=8),dimension(:),allocatable :: f_max
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,accum&
          &,accum2,unc_iter
     real(kind=8),allocatable,dimension(:) :: aux_res,aux_unc,aux_res_iter,aux_res2_iter&
          &,aux_accum,aux_accum2,aux_unc_iter
     integer :: ichan,nevts_unw_gen,evnt,nevnt_in_list,ndim &
          &,current_iter,max_iters,nevts_unw_req,adaptation_class,warmup_iterations
     integer(kind=8) :: npoints_iter,npoints,npoints_requested&
          &,npoints_nonzero,npoints_nonzero_total
     logical :: done,evgen_done,integration_only
     type(evnt),dimension(:),allocatable :: evnt_list
   contains
     procedure,private :: init => integral_init
     procedure,private :: add_point => integral_add_point
     procedure,private :: update_result_iter => integral_update_result_iter
     procedure,private :: combine_iters => integral_combine_iters
     procedure,private :: init_next_iter => integral_init_next_iter
     procedure,private :: compute_fmax => integral_compute_fmax
     procedure,private :: compute_fmax_next_iter => integral_compute_fmax_next_iter
     procedure,private :: unwgt => integral_unwgt
     procedure,private :: update_max_value,check_write_evnt,increase_size_evnt_list,compute_wgts,check_overweight
  end type integral
  ! One monotone one-dimensional adaptive grid.  current maps uniform random
  ! cells to physical integration coordinates; accum stores the adaptation data.
  type :: grid
     integer :: size,size_fill
     real(kind=8),allocatable,dimension(:) :: current,accum,current_for_fillcell
     integer(kind=8),allocatable,dimension(:) :: nhits
     logical :: use_mean_abs
   contains
     procedure,private :: init => grid_init
     procedure,private :: add_point => grid_add_point
     procedure,private :: get_x,get_wgt,massage_accum,find_cell&
          &,interpolate_current,find_cell_to_fill
     procedure,private :: update => grid_update
  end type grid
  ! Candidate event metadata saved during event generation.
  type :: evnt
     real(kind=8),allocatable,dimension(:) :: x,f_abs
     real(kind=8) :: wgt,rnd,overwgt
     integer :: iter,label
     logical :: unwgt
  end type evnt
  ! Public driver object.  Users call init, then alternate get_points and
  ! fill_points until done, then call assign_evnt_wgts if events were written.
  type,public :: integrator
     integer :: nchannel=0,current_channel=0,nevts_unw_req=0,npoints_gen=0,&
          allocation_cursor=0,warmup_iterations=0,event_label=0,&
          candidate_event_limit=0
     integer(kind=8) :: npoints_requested=0_8,invalid_points=0_8
     real(kind=8),dimension(2) :: res=0d0,unc=0d0
     real(kind=8) :: requested_accuracy=0d0
     logical :: event_quotas_initialised=.false.,integration_only=.false.,&
          force_external_refinement=.false.
     real(kind=8),allocatable,dimension(:,:),public :: x
     real(kind=8),allocatable,dimension(:),public :: wgt
     type(channel),allocatable,dimension(:) :: channels
     ! x and wgt are allocated by get_points and released by fill_points.
   contains
     procedure,public :: init,get_points,fill_points,compute_wgt_from_x,assign_evnt_wgts,get_channel_results&
          &,get_channel_aux_results
     procedure,private :: read_all_grids,write_all_grids&
          &,get_channel_and_integral,update_points_requested&
          &,print_results,compute_total_rate,update_nevts_unw_req&
          &,count_unweighted_evnts,init_next_iter&
          &,get_npoints_nonzero_iter,get_npoints_iter,finalise_iter,update_grids
  end type integrator
  integer,parameter :: importance_sampling_strategy=3
  ! Fraction of largest candidate weights ignored when estimating f_max.
  real(kind=8),parameter :: write_evnt_fraction=0.05d0
  integer,parameter :: min_points_per_channel=1024
  integer,parameter :: min_points_per_integral=128
  integer(kind=8),parameter :: accuracy_pilot_points_per_stratum=1024_8
  integer(kind=8),parameter :: min_accuracy_points_per_stratum=128_8
  real(kind=8),parameter :: accuracy_exploration_fraction=0.25d0
  integer(kind=8),parameter :: accuracy_quota_growth_limit=16_8
  real(kind=8),parameter :: required_accuracy_factor=10d0
  integer,parameter :: min_grid_size=8
  integer,parameter :: max_grid_size=2048
  integer(kind=8),parameter :: initial_event_buffer_size=16_8
  ! Bound all public-input-driven workspaces before allocation.  The grid
  ! bound is deliberately conservative: it charges every saved iteration for
  ! maximum-size coordinate and fill-cell arrays, even when the actual fill
  ! grid is much smaller.  Candidate storage has a separate global bound
  ! because each retained point owns two allocatable arrays.
  integer(kind=8),parameter :: max_integrator_extent=1000000_8
  integer(kind=8),parameter :: max_integrator_workspace_units=268435456_8
  integer(kind=8),parameter :: max_batch_workspace_units=67108864_8
  integer(kind=8),parameter :: max_candidate_storage_units=67108864_8
  integer(kind=8),parameter :: max_exact_point_budget=9007199254740991_8
  integer(kind=8),parameter :: grid_storage_units=4_8*int(max_grid_size+1,kind=8)+32_8
  integer(kind=8),parameter :: event_record_overhead_units=16_8
  real(kind=8),parameter :: allowed_overweight_factor=0.001d0
  integer,parameter :: final_n_iters_for_evnt_gen=8
  ! Raw sums include squares of every returned value.  This conservative bound
  ! leaves an enormous physical range while keeping both the square and any
  ! feasible accumulation representable.
  real(kind=8),parameter :: integrator_value_limit=0.25d0*huge(1d0)**0.25d0
  real(kind=8),parameter :: integrator_square_floor=sqrt(tiny(1d0))
  real(kind=8),parameter :: event_random_floor=sqrt(tiny(1d0))
contains

  integer(kind=8) function checked_count_product(first,second,context)
    implicit none
    integer(kind=8),intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0_8 .or. second.lt.0_8) then
       write (*,*) 'ERROR: negative count in ',trim(context),first,second
       stop 1
    endif
    if (first.ne.0_8) then
       if (second.gt.huge(checked_count_product)/first) then
          write (*,*) 'ERROR: 64-bit count overflow in ',trim(context),first,second
          stop 1
       endif
    endif
    checked_count_product=first*second
  end function checked_count_product

  integer(kind=8) function checked_count_sum(first,second,context)
    implicit none
    integer(kind=8),intent(in) :: first,second
    character(len=*),intent(in) :: context
    if (first.lt.0_8 .or. second.lt.0_8) then
       write (*,*) 'ERROR: negative count in ',trim(context),first,second
       stop 1
    endif
    if (second.gt.huge(checked_count_sum)-first) then
       write (*,*) 'ERROR: 64-bit count overflow in ',trim(context),first,second
       stop 1
    endif
    checked_count_sum=first+second
  end function checked_count_sum

  integer(kind=8) function checked_point_sum(first,second,context)
    implicit none
    integer(kind=8),intent(in) :: first,second
    character(len=*),intent(in) :: context
    checked_point_sum=checked_count_sum(first,second,context)
    if (checked_point_sum.gt.max_exact_point_budget) then
       write (*,*) 'ERROR: integration point count exceeds exact supported range in ',&
            trim(context),checked_point_sum,max_exact_point_budget
       stop 1
    endif
  end function checked_point_sum

  ! Initialise the public integrator object and all channel/integral state.
  subroutine init(this,nchannel,ndim,ndim_extra,nintegral,nevts_unw_req,niters,accuracy,naux,adaptation_classes)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: nchannel,nevts_unw_req,niters
    integer,dimension(:),intent(in) :: ndim,nintegral,ndim_extra
    real(kind=8),intent(in),optional :: accuracy
    integer,intent(in),optional :: naux
    integer,dimension(:,:),intent(in),optional :: adaptation_classes
    integer :: i,naux_local,iclass,max_class,ios
    integer(kind=8) :: initial_channel_points,total_integrals,min_channel_points,&
         min_integral_points,min_total_integral_points,workspace_units,grid_slots,&
         channel_units,integral_units,per_integral_units,channel_aux_units,&
         event_footprint,candidate_limit8
    character(len=256) :: allocation_message
    if (nchannel.lt.1) then
       write (*,*) 'ERROR: nchannel must be at least 1'
       stop 1
    endif
    if (int(nchannel,kind=8).gt.max_integrator_extent) then
       write (*,*) 'ERROR: nchannel exceeds the supported integrator extent:',nchannel
       stop 1
    endif
    if (size(ndim).ne.nchannel .or. size(ndim_extra).ne.nchannel .or. &
         size(nintegral).ne.nchannel) then
       write (*,*) 'ERROR: integrator channel metadata has incompatible shape',&
            nchannel,size(ndim),size(ndim_extra),size(nintegral)
       stop 1
    endif
    if (niters.lt.1) then
       write (*,*) 'ERROR: niters must be at least 1'
       stop 1
    endif
    if (int(niters,kind=8).gt.max_integrator_extent) then
       write (*,*) 'ERROR: niters exceeds the supported integrator extent:',niters
       stop 1
    endif
    if (nevts_unw_req.lt.1) then
       write (*,*) 'ERROR: nevts_unw_req must be at least 1'
       stop 1
    endif
    this%integration_only=.false.
    this%force_external_refinement=.false.
    this%requested_accuracy=0d0
    if (present(accuracy)) then
       if (.not.ieee_is_finite(accuracy)) then
          write (*,*) 'ERROR: requested accuracy must be between 0 and 1'
          stop 1
       endif
       if (accuracy.le.0d0 .or. accuracy.ge.1d0) then
          write (*,*) 'ERROR: requested accuracy must be between 0 and 1'
          stop 1
       endif
       this%integration_only=.true.
       this%requested_accuracy=accuracy
    endif
    if (any(ndim.lt.1)) then
       write (*,*) 'ERROR: all channels must have at least one adapted dimension'
       stop 1
    endif
    if (any(int(ndim,kind=8).gt.max_integrator_extent)) then
       write (*,*) 'ERROR: an adapted dimension exceeds the supported integrator extent'
       stop 1
    endif
    if (any(ndim_extra.lt.0)) then
       write (*,*) 'ERROR: ndim_extra cannot be negative'
       stop 1
    endif
    if (any(int(ndim_extra,kind=8).gt.max_integrator_extent)) then
       write (*,*) 'ERROR: a flat dimension exceeds the supported integrator extent'
       stop 1
    endif
    if (any(nintegral.lt.1)) then
       write (*,*) 'ERROR: all channels must have at least one integral'
       stop 1
    endif
    if (any(int(nintegral,kind=8).gt.max_integrator_extent)) then
       write (*,*) 'ERROR: a channel integral count exceeds the supported integrator extent'
       stop 1
    endif
    total_integrals=0_8
    do i=1,nchannel
       if (total_integrals.gt.huge(total_integrals)-int(nintegral(i),kind=8)) then
          write (*,*) 'ERROR: total integral count overflows 64-bit integer'
          stop 1
       endif
       total_integrals=total_integrals+int(nintegral(i),kind=8)
    enddo
    if (present(adaptation_classes)) then
       if (size(adaptation_classes,1).lt.maxval(nintegral) .or. &
            size(adaptation_classes,2).ne.nchannel) then
          write (*,*) 'ERROR: adaptation_classes has incompatible shape'
          stop 1
       endif
       do i=1,nchannel
          if (any(adaptation_classes(1:nintegral(i),i).lt.1)) then
             write (*,*) 'ERROR: adaptation classes must be positive for channel',i
             stop 1
          endif
          max_class=maxval(adaptation_classes(1:nintegral(i),i))
          if (max_class.gt.nintegral(i)) then
             write (*,*) 'ERROR: adaptation class exceeds the channel integral count',i,max_class
             stop 1
          endif
          do iclass=1,max_class
             if (count(adaptation_classes(1:nintegral(i),i).eq.iclass).eq.0) then
                write (*,*) 'ERROR: adaptation classes must be contiguous for channel',i
                stop 1
             endif
          enddo
       enddo
    endif
    naux_local=0
    if (present(naux)) naux_local=naux
    if (naux_local.lt.0) then
       write (*,*) 'ERROR: naux cannot be negative'
       stop 1
    endif
    if (int(naux_local,kind=8).gt.max_integrator_extent) then
       write (*,*) 'ERROR: naux exceeds the supported integrator extent:',naux_local
       stop 1
    endif

    ! Reject metadata combinations whose eventual saved-grid and accumulator
    ! state would exceed the documented process-local workspace ceiling.
    workspace_units=0_8
    do i=1,nchannel
       max_class=1
       if (present(adaptation_classes)) &
            max_class=maxval(adaptation_classes(1:nintegral(i),i))
       grid_slots=checked_count_product(int(ndim(i),kind=8),int(niters,kind=8)+1_8,&
            'integrator grid slots')
       grid_slots=checked_count_product(grid_slots,int(max_class,kind=8),&
            'integrator adaptation grid slots')
       channel_units=checked_count_product(grid_slots,grid_storage_units,&
            'integrator grid storage')
       workspace_units=checked_count_sum(workspace_units,channel_units,&
            'integrator workspace')
    enddo
    per_integral_units=int(niters,kind=8)+7_8*int(naux_local,kind=8)+64_8
    if (.not.this%integration_only) per_integral_units=checked_count_sum(per_integral_units,&
         checked_count_product(initial_event_buffer_size,event_record_overhead_units,&
         'initial candidate-event records'),'integral workspace')
    integral_units=checked_count_product(total_integrals,per_integral_units,&
         'integral state storage')
    workspace_units=checked_count_sum(workspace_units,integral_units,&
         'integrator workspace')
    channel_aux_units=checked_count_product(int(nchannel,kind=8),&
         4_8*int(naux_local,kind=8)+32_8,'channel state storage')
    workspace_units=checked_count_sum(workspace_units,channel_aux_units,&
         'integrator workspace')
    if (workspace_units.gt.max_integrator_workspace_units) then
       write (*,*) 'ERROR: integrator metadata exceeds the supported workspace:',&
            workspace_units,max_integrator_workspace_units
       stop 1
    endif

    event_footprint=checked_count_sum(int(niters,kind=8),int(maxval(ndim),kind=8),&
         'candidate-event footprint')
    event_footprint=checked_count_sum(event_footprint,event_record_overhead_units,&
         'candidate-event footprint')
    ! Buffer growth briefly owns both the old and the deep-copied allocatable
    ! payloads, so reserve a factor of two for that transient peak.
    event_footprint=checked_count_product(event_footprint,2_8,&
         'candidate-event growth footprint')
    candidate_limit8=max(1_8,max_candidate_storage_units/event_footprint)
    candidate_limit8=min(candidate_limit8,int(huge(this%candidate_event_limit),kind=8))
    if (.not.this%integration_only .and. int(nevts_unw_req,kind=8).gt.candidate_limit8) then
       write (*,*) 'ERROR: requested events exceed the retained-candidate workspace:',&
            nevts_unw_req,candidate_limit8
       stop 1
    endif
    if (allocated(this%channels)) deallocate(this%channels)
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%wgt)) deallocate(this%wgt)
    this%event_label=0
    this%candidate_event_limit=int(candidate_limit8)
    this%invalid_points=0_8
    this%nchannel=nchannel
    this%nevts_unw_req=nevts_unw_req
    allocate(this%channels(this%nchannel),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate integrator channels: ',trim(allocation_message)
       stop 1
    endif
    if (this%integration_only) then
       this%warmup_iterations=0
       ! Accuracy-only runs start with equal, event-independent pilot samples
       ! for every channel/integral leaf stratum.
       this%npoints_requested=checked_count_product(accuracy_pilot_points_per_stratum,total_integrals,&
            'accuracy pilot point count')
       do i=1,this%nchannel
          initial_channel_points=checked_count_product(accuracy_pilot_points_per_stratum,&
               int(nintegral(i),kind=8),'channel pilot point count')
          if (present(adaptation_classes)) then
             call this%channels(i)%init(ndim(i),ndim_extra(i),nintegral(i),initial_channel_points,niters,i,&
                  naux_local,this%integration_only,this%warmup_iterations,&
                  adaptation_classes(1:nintegral(i),i))
          else
             call this%channels(i)%init(ndim(i),ndim_extra(i),nintegral(i),initial_channel_points,&
                  niters,i,naux_local,this%integration_only,this%warmup_iterations)
          endif
       enddo
    else
       ! If we assume 1% unweighting efficiency, we expect ~10% time
       ! spent in iterations that do not produce events.
       this%warmup_iterations=5
       this%npoints_requested=int(nevts_unw_req/(0.1d0*2**this%warmup_iterations),kind=8)
       min_channel_points=int(min_points_per_channel,kind=8)
       min_integral_points=checked_count_product(int(min_points_per_integral,kind=8),&
            int(maxval(nintegral),kind=8),'minimum integral point count')
       do while (this%npoints_requested/int(this%nchannel,kind=8).lt.&
            max(min_channel_points,min_integral_points) &
            .and. this%warmup_iterations.gt.3)
          this%warmup_iterations=this%warmup_iterations-1
          this%npoints_requested=int(nevts_unw_req/(0.03d0*2**this%warmup_iterations),kind=8)
       enddo
       min_channel_points=checked_count_product(min_channel_points,int(this%nchannel,kind=8),&
            'minimum channel point budget')
       min_total_integral_points=checked_count_product(min_integral_points,int(this%nchannel,kind=8),&
            'minimum integral point budget')
       this%npoints_requested=max(this%npoints_requested,min_channel_points,min_total_integral_points)
       do i=1,this%nchannel
          if (present(adaptation_classes)) then
             call this%channels(i)%init(ndim(i),ndim_extra(i),nintegral(i),&
                  this%npoints_requested/int(nchannel,kind=8),niters,i,&
                  naux_local,this%integration_only,this%warmup_iterations,&
                  adaptation_classes(1:nintegral(i),i))
          else
             call this%channels(i)%init(ndim(i),ndim_extra(i),nintegral(i),&
                  this%npoints_requested/int(nchannel,kind=8),&
                  niters,i,naux_local,this%integration_only,this%warmup_iterations)
          endif
       enddo
    endif
    this%current_channel=0
    this%allocation_cursor=0
    this%npoints_gen=0
    this%event_quotas_initialised=.false.
    this%res=0d0
    this%unc=0d0
  end subroutine init

  ! Initialise one channel, including its first set of grids and integrals.
  subroutine channel_init(this,ndim,ndim_extra,nintegral,npoints,niters,ichan,naux,&
       integration_only,warmup_iterations,adaptation_classes)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: ndim,nintegral,niters,ichan,ndim_extra,naux,warmup_iterations
    logical,intent(in) :: integration_only
    integer(kind=8) :: npoints
    integer,dimension(:),intent(in),optional :: adaptation_classes
    integer :: i,iadapt,ios
    integer(kind=8) :: adaptation_points,class_count,integral_count
    character(len=256) :: allocation_message
    this%ndim=ndim
    this%ndim_extra=ndim_extra
    this%max_iters=niters
    this%nintegral=nintegral
    this%naux=naux
    this%integration_only=integration_only
    this%warmup_iterations=warmup_iterations
    this%nadaptation=1
    allocate(this%adaptation_class(this%nintegral),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate channel adaptation classes: ',&
            trim(allocation_message)
       stop 1
    endif
    this%adaptation_class=1
    if (present(adaptation_classes)) then
       if (size(adaptation_classes).ne.this%nintegral) then
          write (*,*) 'ERROR: channel adaptation-class count differs from integral count'
          stop 1
       endif
       this%adaptation_class=adaptation_classes
       this%nadaptation=maxval(this%adaptation_class)
    endif
    this%number=ichan
    this%nevts_unw_req=0
    allocate(this%grids(1:this%ndim,1:this%max_iters+1,1:this%nadaptation),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate channel grids: ',trim(allocation_message)
       stop 1
    endif
    allocate(this%integrals(1:this%nintegral),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate channel integrals: ',trim(allocation_message)
       stop 1
    endif
    do iadapt=1,this%nadaptation
       class_count=int(count(this%adaptation_class.eq.iadapt),kind=8)
       integral_count=int(this%nintegral,kind=8)
       adaptation_points=(npoints/integral_count)*class_count+&
            (mod(npoints,integral_count)*class_count)/integral_count
       adaptation_points=max(1_8,adaptation_points)
       do i=1,this%ndim
          call this%grids(i,1,iadapt)%init(adaptation_points)
       enddo
    enddo
    do i=1,this%nintegral
       call this%integrals(i)%init(ndim,npoints/this%nintegral,this%number,this%max_iters,naux,&
            this%adaptation_class(i),this%integration_only,this%warmup_iterations)
    enddo
    this%current_integral=0
    this%current_iter=0
    this%npoints=0_8
    this%res=0d0
    this%unc=0d0
    allocate(this%aux_res(naux),this%aux_unc(naux),this%aux_res_iter(naux),&
         this%aux_unc_iter(naux),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate channel auxiliary accumulators: ',&
            trim(allocation_message)
       stop 1
    endif
    this%aux_res=0d0
    this%aux_unc=0d0
    this%overweight=0d0
    call this%init_next_iter()
    this%evgen_done=.false.
  end subroutine channel_init

  ! Reset channel accumulators and advance to the next iteration.
  subroutine channel_init_next_iter(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    this%res_iter=0d0
    this%res2_iter=0d0
    this%unc_iter=0d0
    this%aux_res_iter=0d0
    this%aux_unc_iter=0d0
    this%npoints_iter=0_8
    this%done=.false.
    if (.not.this%integration_only .and. all(this%integrals%evgen_done)) then
       this%evgen_done=.true.
       this%done=.true.
    else
       this%evgen_done=.false.
    endif
    do i=1,this%nintegral
       call this%integrals(i)%init_next_iter(this)
    enddo
    this%current_iter=this%current_iter+1
  end subroutine channel_init_next_iter
  
  ! Reset one integral for a new iteration and choose the active f_max estimate.
  subroutine integral_init_next_iter(this,thischan)
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    this%current_iter=this%current_iter+1
    this%res_iter=0d0
    this%res2_iter=0d0
    this%accum=0d0
    this%accum2=0d0
    this%unc_iter=0d0
    this%aux_res_iter=0d0
    this%aux_res2_iter=0d0
    this%aux_accum=0d0
    this%aux_accum2=0d0
    this%aux_unc_iter=0d0
    this%npoints_iter=0_8
    this%npoints_nonzero=0_8
    this%evnt=0
    if (this%integration_only) then
       this%evgen_done=.false.
       this%done=.false.
    elseif (.not. this%evgen_done) then
       this%done=.false.
    endif
    if (this%current_iter.eq.1) then
       this%f_max(this%current_iter)=-1d0
    elseif (this%current_iter.le.this%warmup_iterations+1) then
       this%f_max(this%current_iter)=this%max_value
    else
       call this%compute_fmax_next_iter(thischan)
    endif
    this%max_value=0d0
  end subroutine integral_init_next_iter
  
  ! Initialise a grid either uniformly or by interpolating a previous grid.
  subroutine grid_init(this,npoints,current)
    implicit none
    class(grid),intent(inout) :: this
    integer(kind=8),intent(in) :: npoints
    real(kind=8),dimension(:),intent(in),optional :: current
    integer :: i,isize,ios
    character(len=256) :: allocation_message
    this%size=max_grid_size
    call determine_sizefill(npoints,this%size_fill)
    allocate(this%current(0:this%size),this%current_for_fillcell(0:this%size_fill),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate adaptive-grid coordinates: ',&
            trim(allocation_message)
       stop 1
    endif
    if (present(current)) then
       isize=size(current)-1
       if (isize.ne.this%size_fill) then
          call this%interpolate_current(isize,this%size_fill,current,this%current_for_fillcell)
       else
          this%current_for_fillcell=current
       endif
       call this%interpolate_current(isize,this%size,current,this%current)
    else
       do i=0,this%size_fill
          this%current_for_fillcell(i)=dble(i)/this%size_fill
       enddo
       do i=0,this%size
          this%current(i)=dble(i)/this%size
       enddo
    endif
    allocate(this%accum(0:this%size_fill),this%nhits(this%size_fill),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate adaptive-grid accumulators: ',&
            trim(allocation_message)
       stop 1
    endif
    this%accum=0d0
    this%nhits=0_8
    this%use_mean_abs=.false.
  end subroutine grid_init

  ! Choose the number of adaptation fill cells from the requested statistics.
  subroutine determine_sizefill(npoints,isize)
    implicit none
    integer(kind=8),intent(in) :: npoints
    integer,intent(out) :: isize
    isize=int(sqrt(dble(npoints))/10)
    isize=max(isize,min_grid_size)
    isize=min(isize,max_grid_size)
  end subroutine determine_sizefill
  
  ! Initialise one integral estimate and its candidate-event buffer.
  subroutine integral_init(this,ndim,npoints,ichan,niters,naux,adaptation_class,&
       integration_only,warmup_iterations)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: ndim,ichan,niters,naux,adaptation_class,warmup_iterations
    logical,intent(in) :: integration_only
    integer(kind=8) :: npoints,event_buffer_extent
    integer :: ios
    character(len=256) :: allocation_message
    this%ndim=ndim
    this%ichan=ichan
    this%adaptation_class=adaptation_class
    this%integration_only=integration_only
    this%warmup_iterations=warmup_iterations
    this%npoints=0_8
    this%npoints_iter=0_8
    this%npoints_nonzero=0_8
    this%npoints_requested=npoints
    this%max_iters=niters
    this%nevts_unw_req=0
    this%nevts_unw_gen=0
    this%evnt=0
    allocate(this%f_max(this%max_iters),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate event-envelope history: ',&
            trim(allocation_message)
       stop 1
    endif
    this%f_max=-1d0
    if (this%integration_only) then
       event_buffer_extent=1_8
    else
       event_buffer_extent=max(1_8,min(npoints,initial_event_buffer_size))
    endif
    allocate(this%evnt_list(event_buffer_extent),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate initial candidate-event buffer:',&
            event_buffer_extent,trim(allocation_message)
       stop 1
    endif
    this%evnt_list%unwgt=.false.
    this%evnt_list%overwgt=0d0
    this%nevnt_in_list=0
    this%evgen_done=.false.
    this%done=.false.
    this%current_iter=0
    this%npoints_nonzero_total=0_8
    this%overweight=0d0
    this%max_value=0d0
    this%res=0d0
    this%unc=0d0
    this%res_iter=0d0
    this%res2_iter=0d0
    this%unc_iter=0d0
    this%accum=0d0
    this%accum2=0d0
    allocate(this%aux_res(naux),this%aux_unc(naux),this%aux_res_iter(naux),this%aux_res2_iter(naux),&
         this%aux_accum(naux),this%aux_accum2(naux),this%aux_unc_iter(naux),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate integral auxiliary accumulators: ',&
            trim(allocation_message)
       stop 1
    endif
    this%aux_res=0d0
    this%aux_unc=0d0
    this%aux_res_iter=0d0
    this%aux_res2_iter=0d0
    this%aux_accum=0d0
    this%aux_accum2=0d0
    this%aux_unc_iter=0d0
  end subroutine integral_init
  
  ! Select an active channel/integral and generate a batch of random points.
  subroutine get_points(this,npoints,ichan,iint)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: npoints
    integer,intent(out) :: ichan,iint
    integer :: i,ntot,ios
    integer(kind=8) :: batch_units
    real(kind=8) :: wgt_chan
    character(len=256) :: allocation_message
    if (npoints.lt.1) then
       write (*,*) 'ERROR: get_points requires at least one point'
       stop 1
    endif
    if (this%npoints_gen.ne.0) then
       write (*,*) 'ERROR: previous points must be returned with fill_points before get_points is called again'
       stop 1
    endif
    if (all(this%channels%done .or. this%channels%evgen_done)) then
       write (*,*) 'ERROR: get_points called after integration is done'
       stop 1
    endif
    
    call this%get_channel_and_integral(ichan,iint,wgt_chan)
    this%current_channel=ichan
    if (this%integration_only) then
       if (int(npoints,kind=8).gt.this%channels(ichan)%integrals(iint)%npoints_requested-&
            this%channels(ichan)%integrals(iint)%npoints_iter) then
          write (*,*) 'ERROR: accuracy-only batch exceeds its leaf quota'
          stop 1
       endif
    endif

    if (this%channels(this%current_channel)%ndim.gt.&
         huge(ntot)-this%channels(this%current_channel)%ndim_extra) then
       write (*,*) 'ERROR: integration point dimension overflows default integer'
       stop 1
    endif
    ntot=this%channels(this%current_channel)%ndim+this%channels(this%current_channel)%ndim_extra
    batch_units=checked_count_product(int(ntot,kind=8),int(npoints,kind=8),&
         'integration point batch')
    batch_units=checked_count_sum(batch_units,int(npoints,kind=8),&
         'integration point batch')
    if (batch_units.gt.max_batch_workspace_units) then
       write (*,*) 'ERROR: integration point batch exceeds the supported workspace:',&
            batch_units,max_batch_workspace_units
       stop 1
    endif
    allocate(this%x(1:ntot,1:npoints),this%wgt(1:npoints),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate integration point batch: ',&
            trim(allocation_message)
       stop 1
    endif

    do i=1,npoints
       call this%channels(this%current_channel)%get_point(this%x(1,i),this%wgt(i))
    enddo
    this%wgt=this%wgt*wgt_chan
    this%npoints_gen=npoints

  end subroutine get_points
  
  ! Return evaluated values for the most recent batch and trigger iteration
  ! finalisation when all active channels/integrals have enough statistics.
  ! iteration_finished lets clients synchronize per-iteration diagnostics
  ! without exposing the integrator's private channel state.
  subroutine fill_points(this,npoints,f_abs,f,to_write,done,accepted,f_aux,iteration_finished,external_converged)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: npoints
    real(kind=8),dimension(npoints),intent(in) :: f,f_abs
    logical,dimension(npoints),intent(out) :: to_write
    logical,intent(out) :: done
    logical,dimension(npoints),intent(in),optional :: accepted
    real(kind=8),dimension(:,:),intent(in),optional :: f_aux
    logical,intent(out),optional :: iteration_finished
    logical,intent(in),optional :: external_converged
    integer :: i,ndim,ios
    logical :: convergence_ok,point_is_valid,point_is_accepted
    real(kind=8) :: safe_f,safe_f_abs,safe_wgt
    real(kind=8),allocatable :: safe_x(:),safe_aux(:)
    character(len=256) :: allocation_message
    done=.false.
    convergence_ok=.true.
    if (present(external_converged)) convergence_ok=external_converged
    if (present(iteration_finished)) iteration_finished=.false.
    if (this%npoints_gen.eq.0) then
       write (*,*) 'ERROR: fill_points called before get_points'
       stop 1
    endif
    if (npoints.lt.1) then
       write (*,*) 'ERROR: fill_points requires at least one point'
       stop 1
    endif
    if (npoints.gt.this%npoints_gen) then
       write (*,*) 'ERROR: too many points returned'
       stop 1
    endif
    if (present(f_aux)) then
       if (size(f_aux,1).ne.this%channels(this%current_channel)%naux .or. size(f_aux,2).ne.npoints) then
          write (*,*) 'ERROR: f_aux has incompatible shape'
          stop 1
       endif
    elseif (this%channels(this%current_channel)%naux.ne.0) then
       write (*,*) 'ERROR: fill_points requires f_aux after init with naux > 0'
       stop 1
    endif
    ndim=this%channels(this%current_channel)%ndim
    allocate(safe_x(ndim),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate point-validation workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    if (present(f_aux)) then
       allocate(safe_aux(size(f_aux,1)),stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'ERROR: could not allocate auxiliary-validation workspace: ',&
               trim(allocation_message)
          stop 1
       endif
    endif
    do i=1,npoints
       safe_x=this%x(1:ndim,i)
       safe_wgt=this%wgt(i)
       safe_f_abs=f_abs(i)
       safe_f=f(i)
       point_is_valid=all(ieee_is_finite(this%x(:,i)))
       if (point_is_valid) point_is_valid=all(this%x(:,i).ge.0d0) .and. &
            all(this%x(:,i).le.1d0)
       if (point_is_valid) point_is_valid=ieee_is_finite(safe_wgt)
       if (point_is_valid) point_is_valid=safe_wgt.gt.0d0 .and. &
            safe_wgt.le.integrator_value_limit
       if (point_is_valid) point_is_valid=squareable_integrator_value(safe_f_abs)
       if (point_is_valid) point_is_valid=squareable_integrator_value(safe_f)
       if (point_is_valid) point_is_valid=safe_f_abs.ge.0d0
       if (point_is_valid) point_is_valid=safe_f_abs+64d0*epsilon(1d0)*&
            max(1d0,safe_f_abs,abs(safe_f)).ge.abs(safe_f)
       if (present(f_aux)) then
          safe_aux=f_aux(:,i)
          point_is_valid=point_is_valid .and. all(squareable_integrator_value(safe_aux))
       endif
       if (.not.point_is_valid) then
          call report_invalid_integrand_point(this%invalid_points,this%current_channel,&
               this%channels(this%current_channel)%current_integral)
          safe_x=0.5d0
          safe_wgt=0d0
          safe_f_abs=0d0
          safe_f=0d0
          if (present(f_aux)) safe_aux=0d0
       endif
       if (present(accepted)) then
          point_is_accepted=accepted(i) .and. point_is_valid
       else
          point_is_accepted=point_is_valid .and. safe_f_abs.gt.0d0
       endif
       if (present(f_aux)) then
          call this%channels(this%current_channel)%add_point(safe_x,safe_wgt,safe_f_abs,safe_f,to_write(i),&
               point_is_accepted,this%event_label,this%candidate_event_limit,safe_aux)
       else
          call this%channels(this%current_channel)%add_point(safe_x,safe_wgt,safe_f_abs,safe_f,to_write(i),&
               point_is_accepted,this%event_label,this%candidate_event_limit)
       endif
    enddo
    deallocate(safe_x)
    if (allocated(safe_aux)) deallocate(safe_aux)
    this%npoints_gen=0
    if (all(this%channels%done)) then
       call this%finalise_iter(done,convergence_ok)
       if (present(iteration_finished)) iteration_finished=.true.
    endif
    deallocate(this%x)
    deallocate(this%wgt)
  end subroutine fill_points

  subroutine report_invalid_integrand_point(invalid_points,ichan,iint)
    implicit none
    integer,intent(in) :: ichan,iint
    integer(kind=8),intent(inout) :: invalid_points

    invalid_points=checked_count_sum(invalid_points,1_8,&
         'invalid-integrand point counter')
    if (invalid_points.le.10_8 .or. mod(invalid_points,1000_8).eq.0_8) then
       write(*,'(a,i0,a,2(i0,1x))') 'WARNING: rejected invalid integrator point, count=',&
            invalid_points,' channel/integral=',ichan,iint
       write(99,'(a,i0,a,2(i0,1x))') 'WARNING: rejected invalid integrator point, count=',&
            invalid_points,' channel/integral=',ichan,iint
    elseif (invalid_points.eq.11_8) then
       write(*,'(a)') 'WARNING: further invalid-integrator-point messages suppressed'
       write(99,'(a)') 'WARNING: further invalid-integrator-point messages suppressed'
    endif
  end subroutine report_invalid_integrand_point

  elemental pure logical function squareable_integrator_value(value)
    real(kind=8),intent(in) :: value
    squareable_integrator_value=.false.
    if (.not.ieee_is_finite(value)) return
    if (abs(value).gt.integrator_value_limit) return
    if (value.ne.0d0 .and. abs(value).lt.integrator_square_floor) return
    squareable_integrator_value=.true.
  end function squareable_integrator_value

  real(kind=8) function validated_uniform_random(stage)
    implicit none
    character(len=*),intent(in) :: stage
    validated_uniform_random=ran2()
    if (.not.ieee_is_finite(validated_uniform_random)) then
       write (*,*) 'ERROR: random-number generator returned a non-finite value during ',trim(stage)
       stop 1
    endif
    if (validated_uniform_random.lt.0d0 .or. validated_uniform_random.ge.1d0) then
       write (*,*) 'ERROR: random-number generator returned a value outside [0,1) during ',trim(stage),&
            ':',validated_uniform_random
       stop 1
    endif
  end function validated_uniform_random

  ! Finish one global iteration: combine rates, unweight candidates, update
  ! grids, and decide whether the requested event sample is complete.
  subroutine finalise_iter(this,done,external_converged)
    implicit none
    class(integrator),intent(inout) :: this
    logical,intent(out) :: done
    logical,intent(in) :: external_converged
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    character(len=19) :: formatted
    integer(kind=8) :: npoints_report
    call date_and_time(date, time, zone)
    write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
    call this%get_npoints_iter(npoints_report)
    write (*,*) ''
    write (*,'(a,x,i4,x,a,x,i10,x,a)') &
         'iteration',this%channels(1)%current_iter,'(',npoints_report, &
         'points) '//trim(formatted)//' :'
    write (99,*) ''
    write (99,'(a,x,i4,x,a,x,i10,x,a)') &
         'iteration',this%channels(1)%current_iter,'(',npoints_report, &
         'points) '//trim(formatted)//' :'
    call this%compute_total_rate()
    if (.not.this%integration_only) call this%count_unweighted_evnts()
    call this%print_results()
    this%force_external_refinement=this%integration_only .and. .not.external_converged
    call this%update_grids()
    call this%init_next_iter()
    this%force_external_refinement=.false.
    if (.not.this%integration_only .and. all(this%channels%evgen_done)) done=.true.
    if (this%integration_only) then
       if (this%res(1).gt.0d0) then
          if (this%unc(1)/this%res(1).lt.this%requested_accuracy .and. external_converged) done=.true.
       endif
    endif
    if (all(this%channels%done)) done=.true.
    call flush(99)
  end subroutine finalise_iter

  ! Update all active channel grids after an iteration has been finalised.
  subroutine update_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    logical :: accuracy_refining
    accuracy_refining=.false.
    if (this%integration_only .and. this%res(1).gt.0d0) then
       accuracy_refining=this%unc(1)/this%res(1).gt.2d0*this%requested_accuracy
    endif
    if (this%force_external_refinement) accuracy_refining=.true.
    do i=1,this%nchannel
       call this%channels(i)%update_grids(accuracy_refining)
    enddo
  end subroutine update_grids
  
  ! Count non-zero points from the current iteration across all integrals.
  subroutine get_npoints_nonzero_iter(this,npoints_nonzero)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i,j
    integer(kind=8) :: npoints_nonzero
    npoints_nonzero=0_8
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          npoints_nonzero=checked_point_sum(npoints_nonzero,&
               this%channels(i)%integrals(j)%npoints_nonzero,&
               'non-zero iteration point total')
       enddo
    enddo
  end subroutine get_npoints_nonzero_iter

  ! Count every evaluated point from the current iteration.  Accuracy-only
  ! quotas are defined in terms of evaluations, including zero-weight points.
  subroutine get_npoints_iter(this,npoints)
    implicit none
    class(integrator),intent(in) :: this
    integer(kind=8),intent(out) :: npoints
    integer :: i,j
    npoints=0
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          npoints=checked_point_sum(npoints,&
               this%channels(i)%integrals(j)%npoints_iter,&
               'iteration point total')
       enddo
    enddo
  end subroutine get_npoints_iter
  
  ! Start the next iteration for channels that have not reached max_iters.
  subroutine init_next_iter(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       if (this%channels(i)%current_iter.lt.this%channels(i)%max_iters) then
          call this%channels(i)%init_next_iter()
       endif
    enddo
    call this%update_points_requested()
  end subroutine init_next_iter
    
  ! Distribute requested events, test stored candidates, and mark completed
  ! channel/integral event-generation tasks.
  subroutine count_unweighted_evnts(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    ! Do not let provisional warm-up rates mark zero-quota channels complete.
    ! Once event generation starts, draw the integer channel/integral quotas
    ! exactly once.  Redrawing a one-event quota every iteration can move it to
    ! a different channel, clear the previously selected candidate, and reach
    ! the iteration cap with no event despite adequate candidates overall.
    if (this%channels(1)%current_iter.lt.this%warmup_iterations) return
    if (.not.this%event_quotas_initialised) then
       call this%update_nevts_unw_req
       this%event_quotas_initialised=.true.
    endif
    do i=1,this%nchannel
       call this%channels(i)%check_gen_evnts()
    enddo
  end subroutine count_unweighted_evnts

  ! Split the total requested event count over channels in proportion to their
  ! current absolute integral estimates.
  subroutine update_nevts_unw_req(this)
    use sort_array_mod
    implicit none
    class(integrator),intent(inout) :: this
    real(kind=8),dimension(this%nchannel) :: res
    integer :: nevts_to_distribute,i
    integer,dimension(this%nchannel) :: idx
    real(kind=8) :: total,quota_real
    res=this%channels%res(1)
    if (.not.all(ieee_is_finite(res))) then
       write (*,*) 'ERROR: invalid channel absolute rates during event allocation'
       stop 1
    endif
    if (any(res.lt.0d0)) then
       write (*,*) 'ERROR: invalid channel absolute rates during event allocation'
       stop 1
    endif
    call sort_indices_by_values(res,idx)
    nevts_to_distribute=this%nevts_unw_req
    total=sum(res)
    if (total.le.0d0) then
       do i=1,this%nchannel
          this%channels(i)%nevts_unw_req=0
          call this%channels(i)%update_nevts_unw_req()
       enddo
       return
    endif
    do i=1,this%nchannel-1
       call draw_stochastic_event_quota(nevts_to_distribute,res(idx(i)),total,&
            this%channels(idx(i))%nevts_unw_req,quota_real)
       nevts_to_distribute=nevts_to_distribute-this%channels(idx(i))%nevts_unw_req
       total=max(0d0,total-res(idx(i)))
    enddo
    this%channels(idx(this%nchannel))%nevts_unw_req=nevts_to_distribute
    do i=1,this%nchannel
       call this%channels(i)%update_nevts_unw_req()
    enddo
  end subroutine update_nevts_unw_req

  ! Split one channel's requested event count across its integrals.
  subroutine channel_update_nevts_unw_req(this)
    use sort_array_mod
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%nintegral) :: res
    integer :: nevts_to_distribute,i
    integer,dimension(this%nintegral) :: idx
    real(kind=8) :: total,quota_real
    res=this%integrals%res(1)
    if (.not.all(ieee_is_finite(res))) then
       write (*,*) 'ERROR: invalid integral absolute rates during event allocation'
       stop 1
    endif
    if (any(res.lt.0d0)) then
       write (*,*) 'ERROR: invalid integral absolute rates during event allocation'
       stop 1
    endif
    call sort_indices_by_values(res,idx)
    nevts_to_distribute=this%nevts_unw_req
    total=sum(res)
    if (total.le.0d0) then
       this%integrals%nevts_unw_req=0
       return
    endif
    do i=1,this%nintegral-1
       call draw_stochastic_event_quota(nevts_to_distribute,res(idx(i)),total,&
            this%integrals(idx(i))%nevts_unw_req,quota_real)
       nevts_to_distribute=nevts_to_distribute-this%integrals(idx(i))%nevts_unw_req
       total=max(0d0,total-res(idx(i)))
    enddo
    this%integrals(idx(this%nintegral))%nevts_unw_req=nevts_to_distribute
  end subroutine channel_update_nevts_unw_req

  subroutine draw_stochastic_event_quota(remaining,rate,total,quota,quota_real)
    implicit none
    integer,intent(in) :: remaining
    real(kind=8),intent(in) :: rate,total
    integer,intent(out) :: quota
    real(kind=8),intent(out) :: quota_real
    real(kind=8) :: fraction
    quota=0
    quota_real=0d0
    if (remaining.le.0 .or. rate.le.0d0 .or. total.le.0d0) return
    fraction=min(1d0,max(0d0,rate/total))
    quota_real=min(dble(remaining),dble(remaining)*fraction)
    quota=min(remaining,max(0,int(quota_real)))
    if (quota.lt.remaining) then
       if (validated_uniform_random('event quota allocation').lt.quota_real-dble(quota)) quota=quota+1
    endif
  end subroutine draw_stochastic_event_quota
  
  ! Combine active channel estimates into the public total result.
  subroutine compute_total_rate(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       if (.not.this%channels(i)%evgen_done) then
          call this%channels(i)%update_result_iter()
          call this%channels(i)%combine_iters()
       endif
    enddo
    do i=1,2
       this%res(i)=sum(this%channels(1:this%nchannel)%res(i))
       this%unc(i)=sqrt(sum(this%channels(1:this%nchannel)%unc(i)**2))
    enddo
  end subroutine compute_total_rate
  
  ! Print per-channel and total integration progress to stdout and unit 99.
  subroutine print_results(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    real(kind=8) :: rel_unc
    do i=1,this%nchannel
       if (.not. this%channels(i)%evgen_done) call this%channels(i)%print_result_iter()
       call this%channels(i)%print_combined_result()
    enddo
    if (this%res(1).gt.0d0) then
       rel_unc=this%unc(1)/this%res(1)*100d0
    else
       rel_unc=0d0
    endif
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',rel_unc,'%)'
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral     (accum):',this%res(2),'+/-',this%unc(2)
    write(99,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',rel_unc,'%)'
    write(99,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral     (accum):',this%res(2),'+/-',this%unc(2)
    call flush()
  end subroutine print_results
  
  ! Choose how many non-zero points to request in the next iteration.
  subroutine update_points_requested(this)
    implicit none
    class(integrator),intent(inout) :: this
    real(kind=8) :: total,total_channel
    real(kind=8) :: rel_unc
    integer :: i,j
    integer(kind=8) :: npoints,npoints_channel
    if (this%integration_only) then
       if (this%force_external_refinement) call double_point_budget(this%npoints_requested)
       if (this%res(1).gt.0d0) then
          rel_unc=this%unc(1)/this%res(1)
          if (.not.this%force_external_refinement .and. &
               rel_unc.gt.2d0*this%requested_accuracy) &
               call double_point_budget(this%npoints_requested)
       else
          if (.not.this%force_external_refinement) call double_point_budget(this%npoints_requested)
       endif
       call update_integration_points_requested(this)
       return
    else
       call double_point_budget(this%npoints_requested)
    endif
    npoints=0_8
    total=sum(this%channels%res(1),mask=.not.this%channels%evgen_done)
    if (total.le.0d0) then
       do i=1,this%nchannel
          if (this%channels(i)%evgen_done) cycle
          do j=1,this%channels(i)%nintegral
             if (this%channels(i)%integrals(j)%evgen_done) cycle
             this%channels(i)%integrals(j)%npoints_requested=min_points_per_integral
             npoints=checked_point_sum(npoints,&
                  this%channels(i)%integrals(j)%npoints_requested,&
                  'zero-rate point allocation')
          enddo
       enddo
       this%npoints_requested=npoints
       return
    endif
    do i=1,this%nchannel
       if (this%channels(i)%evgen_done) cycle
       npoints_channel=max(int(this%channels(i)%res(1)/total*dble(this%npoints_requested),kind=8),&
            min_points_per_channel)
       total_channel=sum(this%channels(i)%integrals%res(1),mask=.not.this%channels(i)%integrals%evgen_done)
       do j=1,this%channels(i)%nintegral
          if (this%channels(i)%integrals(j)%evgen_done) cycle
          if (total_channel.gt.0d0) then
             this%channels(i)%integrals(j)%npoints_requested=&
                  max(int(this%channels(i)%integrals(j)%res(1)/total_channel*dble(npoints_channel),kind=8),&
                  min_points_per_integral)
          else
             this%channels(i)%integrals(j)%npoints_requested=min_points_per_integral
          endif
          npoints=checked_point_sum(npoints,&
               this%channels(i)%integrals(j)%npoints_requested,&
               'event-generation point allocation')
       enddo
    enddo
    this%npoints_requested=npoints
  end subroutine update_points_requested

  subroutine double_point_budget(point_budget)
    implicit none
    integer(kind=8),intent(inout) :: point_budget
    if (point_budget.lt.0_8 .or. point_budget.gt.shiftr(max_exact_point_budget,1)) then
       write (*,*) 'ERROR: integration point budget cannot be doubled safely:',point_budget
       stop 1
    endif
    point_budget=2_8*point_budget
  end subroutine double_point_budget

  ! In integration-only mode, allocate directly to channel/integral leaves.
  ! For a leaf with accumulated error sigma/sqrt(n), q=sigma is estimated as
  ! unc(1)*sqrt(n).  Most of the next budget follows the corresponding
  ! variance-optimal water filling.  A fixed exploration share remains uniform
  ! across all leaves, and a per-iteration growth cap prevents one newly seen
  ! tail point from starving the rest of the integration.
  subroutine update_integration_points_requested(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i,j,k,nleaf,ibest,iter,ios
    integer(kind=8) :: budget,nminimum,nleft,exploration_each,growth_budget,&
         quota_sum
    integer(kind=8),allocatable :: nold(:),quota(:),base(:),cap(:),previous_quota(:)
    real(kind=8),allocatable :: q(:),target(:),remainder(:)
    real(kind=8) :: lo,hi,mid,total_target,best_remainder,qscale
    character(len=256) :: allocation_message

    nleaf=0
    do i=1,this%nchannel
       nleaf=nleaf+this%channels(i)%nintegral
    enddo
    if (nleaf.eq.0) return
    nminimum=min_accuracy_points_per_stratum
    budget=max(this%npoints_requested,checked_count_product(nminimum,int(nleaf,kind=8),&
         'minimum accuracy point budget'))
    if (budget.gt.max_exact_point_budget) then
       write (*,*) 'ERROR: accuracy point budget exceeds exact supported range:',budget
       stop 1
    endif
    allocate(nold(nleaf),quota(nleaf),base(nleaf),cap(nleaf),previous_quota(nleaf),&
         q(nleaf),target(nleaf),remainder(nleaf),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate accuracy-allocation workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    k=0
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          k=k+1
          nold(k)=this%channels(i)%integrals(j)%npoints
          previous_quota(k)=this%channels(i)%integrals(j)%npoints_requested
          q(k)=this%channels(i)%integrals(j)%unc(1)*sqrt(dble(max(nold(k),1_8)))
       enddo
    enddo

    ! If a new tail estimate asks for an abrupt global budget increase, reach it
    ! over several iterations.  This keeps the advertised per-leaf growth cap
    ! feasible without silently relaxing it.
    if (any(max(previous_quota,nminimum).gt.&
         max_exact_point_budget/accuracy_quota_growth_limit)) then
       write (*,*) 'ERROR: per-leaf accuracy quota exceeds exact supported range'
       stop 1
    endif
    growth_budget=0_8
    do k=1,nleaf
       growth_budget=checked_point_sum(growth_budget,&
            max(previous_quota(k),nminimum),'accuracy quota growth budget')
    enddo
    if (growth_budget.gt.max_exact_point_budget/accuracy_quota_growth_limit) then
       growth_budget=max_exact_point_budget
    else
       growth_budget=accuracy_quota_growth_limit*growth_budget
    endif
    budget=min(budget,growth_budget)
    exploration_each=int(floor(accuracy_exploration_fraction*dble(budget)/dble(nleaf)),kind=8)
    cap=accuracy_quota_growth_limit*max(previous_quota,nminimum)
    base=min(cap,max(nminimum,exploration_each))

    qscale=maxval(q)
    if (qscale.le.tiny(1d0)) then
       q=1d0
    else
       q=max(q/qscale,1d-12)
    endif

    lo=0d0
    hi=maxval((dble(nold)+dble(cap))/q)+1d0
    do iter=1,100
       mid=0.5d0*(lo+hi)
       total_target=sum(max(dble(base),min(dble(cap),mid*q-dble(nold))))
       if (total_target.lt.dble(budget)) then
          lo=mid
       else
          hi=mid
       endif
    enddo
    target=max(dble(base),min(dble(cap),lo*q-dble(nold)))
    quota=int(floor(target),kind=8)
    remainder=target-dble(quota)
    quota_sum=0_8
    do k=1,nleaf
       quota_sum=checked_point_sum(quota_sum,quota(k),&
            'initial accuracy quota total')
    enddo
    nleft=budget-quota_sum
    do while (nleft.gt.0_8)
       ibest=0
       best_remainder=-huge(1d0)
       do k=1,nleaf
          if (quota(k).ge.cap(k)) cycle
          if (remainder(k).gt.best_remainder) then
             ibest=k
             best_remainder=remainder(k)
          endif
       enddo
       if (ibest.eq.0) then
          write (*,*) 'ERROR: accuracy allocation cannot satisfy its point budget'
          stop 1
       endif
       quota(ibest)=quota(ibest)+1_8
       remainder(ibest)=-1d0
       nleft=nleft-1_8
    enddo
    quota_sum=0_8
    do k=1,nleaf
       quota_sum=checked_point_sum(quota_sum,quota(k),&
            'final accuracy quota total')
    enddo
    if (quota_sum.ne.budget) then
       write (*,*) 'ERROR: accuracy allocation changed its point budget:',&
            quota_sum,budget
       stop 1
    endif

    k=0
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          k=k+1
          this%channels(i)%integrals(j)%npoints_requested=quota(k)
       enddo
    enddo
    this%npoints_requested=quota_sum
    deallocate(nold,quota,base,cap,previous_quota,q,target,remainder)
  end subroutine update_integration_points_requested
  
  ! Build the next iteration's adaptive grids for one channel.
  subroutine channel_update_grids(this,accuracy_refining)
    implicit none
    class(channel),intent(inout) :: this
    logical,intent(in) :: accuracy_refining
    type(grid) :: new_grid
    integer :: i,iadapt
    integer(kind=8) :: adaptation_points
    logical update_grids
    if (this%integration_only) then
       ! All grids follow the global absolute-envelope convergence gate.
       update_grids=accuracy_refining
    elseif (this%res(1).gt.0d0) then
       update_grids=(((.not.this%evgen_done) .and. &
            this%unc(1)/this%res(1).gt.1d0/(sqrt(dble(max(this%nevts_unw_req,1)))*required_accuracy_factor)) .or. &
            this%current_iter.le.this%warmup_iterations) .and. &
            this%npoints_iter.gt.int(this%npoints*0.2d0)
    else
       update_grids=(this%current_iter.le.this%warmup_iterations)
    endif
    if (.not.update_grids) then
       write (99,*) 'keeping grids fixed for channel',this%number
    endif
    do iadapt=1,this%nadaptation
       adaptation_points=0_8
       do i=1,this%nintegral
          if (this%adaptation_class(i).ne.iadapt) cycle
          adaptation_points=checked_point_sum(adaptation_points,&
               this%integrals(i)%npoints_iter,'adaptation-class point total')
       enddo
       do i=1,this%ndim
          if (update_grids .and. adaptation_points.gt.0_8) then
             call this%grids(i,this%current_iter,iadapt)%update(adaptation_points,new_grid)
             if (this%current_iter .lt. this%max_iters) then
                this%grids(i,this%current_iter+1,iadapt)=new_grid
             endif
          elseif (this%current_iter .lt. this%max_iters) then
             this%grids(i,this%current_iter+1,iadapt)=this%grids(i,this%current_iter,iadapt)
          endif
       enddo
    enddo
  end subroutine channel_update_grids

  ! Recompute f_max values, unweight stored candidates, and decide whether one
  ! channel has produced enough acceptable events.
  subroutine channel_check_gen_evnts(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i,j
    logical :: done
    do i=1,this%nintegral
       call this%integrals(i)%compute_fmax(this)
       if (this%integrals(i)%nevts_unw_req.eq.0) then
          this%integrals(i)%evgen_done=.true.
          this%integrals(i)%nevts_unw_gen=0
          this%integrals(i)%overweight=0d0
          do j=1,this%integrals(i)%nevnt_in_list
             this%integrals(i)%evnt_list(j)%unwgt=.false.
          enddo
          cycle
       endif
       call this%integrals(i)%unwgt()
       if (this%integrals(i)%nevts_unw_gen.ge.this%integrals(i)%nevts_unw_req) then
          call this%integrals(i)%check_overweight(done)
          if (done) then
             this%integrals(i)%evgen_done=.true.
          else
             this%integrals(i)%evgen_done=.false.
             if (this%integrals(i)%nevts_unw_gen.ge.this%integrals(i)%nevts_unw_req) then
                if (.not.ieee_is_finite(this%integrals(i)%overweight)) then
                   write (*,*) 'ERROR: invalid candidate-event overweight estimate'
                   stop 1
                endif
                if (this%integrals(i)%overweight.le.0d0) then
                   write (*,*) 'ERROR: invalid candidate-event overweight estimate'
                   stop 1
                endif
                this%integrals(i)%nevts_unw_gen=min(int(this%integrals(i)%nevts_unw_req*0.8d0),&
                     int(this%integrals(i)%nevts_unw_req*allowed_overweight_factor/&
                     this%integrals(i)%overweight))
             endif
          endif
       else
          this%integrals(i)%evgen_done=.false.
       endif
    enddo
    this%overweight=sum(this%integrals%overweight)
  end subroutine channel_check_gen_evnts

  
  ! Accept/reject the best candidate events and measure the overweight excess.
  subroutine check_overweight(this,done)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    logical,intent(out) :: done
    integer :: j,k,eligible,ios
    real(kind=8),allocatable :: fmax(:),log_score(:),log_score_top(:)
    integer,allocatable :: top_idx(:)
    real(kind=8) :: log_score_required,log_overweight,excess
    logical,allocatable :: to_include(:)
    character(len=256) :: allocation_message
    allocate(fmax(this%current_iter),log_score(this%nevnt_in_list),&
         log_score_top(this%nevts_unw_req),top_idx(this%nevts_unw_req),&
         to_include(this%current_iter),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate final-unweighting workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    ! check which iterations to include (only the final
    ! 'final_n_iters_for_evnt_gen' that generated events for this
    ! integral will be included)
    to_include=.false.
    k=0
    do j=this%nevnt_in_list,1,-1
       if (this%evnt_list(j)%iter.lt.1 .or. this%evnt_list(j)%iter.gt.this%current_iter) then
          write (*,*) 'ERROR: invalid stored candidate-event iteration',j,this%evnt_list(j)%iter
          stop 1
       endif
       if (.not.to_include(this%evnt_list(j)%iter)) then
          k=k+1
          to_include(this%evnt_list(j)%iter)=.true.
       endif
       if (k.eq.final_n_iters_for_evnt_gen) exit
    enddo
    ! rescale all f_abs such that they are equivalent for all iterations
    fmax=0d0
    do j=1,this%nevnt_in_list
       this%evnt_list(j)%unwgt=.false.
       this%evnt_list(j)%overwgt=0d0
       do k=this%warmup_iterations+1,this%current_iter
          fmax(k)=max(fmax(k),this%evnt_list(j)%f_abs(k))
       enddo
    enddo
    if (.not.all(ieee_is_finite(fmax))) then
       write (*,*) 'ERROR: invalid candidate-event maximum during final unweighting'
       stop 1
    endif
    if (any(fmax.lt.0d0)) then
       write (*,*) 'ERROR: invalid candidate-event maximum during final unweighting'
       stop 1
    endif
    ! Rank in log space.  This preserves the original ordering by
    ! (f_abs/rnd)/f_max without overflowing for very small random numbers.
    log_score=-huge(1d0)
    eligible=0
    do j=1,this%nevnt_in_list
       if (to_include(this%evnt_list(j)%iter)) then
          k=this%evnt_list(j)%iter
          if (.not.ieee_is_finite(this%evnt_list(j)%rnd)) then
             write (*,*) 'ERROR: invalid candidate-event random number during final unweighting',j
             stop 1
          endif
          if (this%evnt_list(j)%f_abs(k).gt.0d0 .and. fmax(k).gt.0d0 .and. &
               this%evnt_list(j)%rnd.gt.0d0) then
             log_score(j)=log(this%evnt_list(j)%f_abs(k))-log(this%evnt_list(j)%rnd)-log(fmax(k))
             if (.not.ieee_is_finite(log_score(j))) then
                write (*,*) 'ERROR: invalid candidate-event score during final unweighting',j
                stop 1
             endif
             eligible=eligible+1
          endif
       endif
    enddo
    if (eligible.lt.this%nevts_unw_req) then
       ! integral_unwgt counts every stored iteration, whereas the final
       ! overweight test deliberately retains only the latest event-generating
       ! iterations.  Keep sampling when that restricted pool is still short.
       this%nevts_unw_gen=eligible
       this%overweight=0d0
       done=.false.
       deallocate(fmax,log_score,log_score_top,top_idx,to_include)
       return
    endif
    ! Take the nevts_unw_req largest
    k=this%nevts_unw_req
    call topk_largest(log_score,k,log_score_top,top_idx)
    ! Find the logarithm of the threshold such that all selected events remain.
    log_score_required=log_score_top(k)
    ! check the overweight fraction
    excess=0d0
    do j=1,this%nevts_unw_req
       this%evnt_list(top_idx(j))%unwgt=.true.
       k=this%evnt_list(top_idx(j))%iter
       log_overweight=log(this%evnt_list(top_idx(j))%f_abs(k))-log(fmax(k))-log_score_required
       if (log_overweight.gt.log(integrator_value_limit)) then
          this%evnt_list(top_idx(j))%overwgt=integrator_value_limit
       elseif (log_overweight.lt.log(tiny(1d0))) then
          this%evnt_list(top_idx(j))%overwgt=0d0
       else
          this%evnt_list(top_idx(j))%overwgt=exp(log_overweight)
       endif
       if (this%evnt_list(top_idx(j))%overwgt.gt.1d0) &
            excess=excess+this%evnt_list(top_idx(j))%overwgt-1d0
    enddo
    this%overweight=excess/dble(this%nevts_unw_req)
    if (this%overweight.lt.allowed_overweight_factor) then
       done=.true.
    else
       done=.false.
    endif
    deallocate(fmax,log_score,log_score_top,top_idx,to_include)
  end subroutine check_overweight
  
  ! Count events that pass the current iteration-by-iteration f_max thresholds.
  subroutine integral_unwgt(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: j,iter
    this%nevts_unw_gen=0
    do j=1,this%nevnt_in_list
       iter=this%evnt_list(j)%iter
       if (iter.lt.1 .or. iter.gt.this%current_iter) then
          write (*,*) 'ERROR: invalid stored candidate-event iteration',j,iter
          stop 1
       endif
       if (.not.ieee_is_finite(this%evnt_list(j)%f_abs(iter)) .or. &
            .not.ieee_is_finite(this%f_max(iter)) .or. &
            .not.ieee_is_finite(this%evnt_list(j)%rnd)) then
          write (*,*) 'ERROR: invalid stored candidate event during unweighting',j,iter
          stop 1
       endif
       if (this%evnt_list(j)%f_abs(iter).lt.0d0 .or. this%f_max(iter).lt.0d0 .or. &
            this%evnt_list(j)%rnd.le.0d0) then
          write (*,*) 'ERROR: invalid stored candidate event during unweighting',j,iter
          stop 1
       endif
       if (this%evnt_list(j)%f_abs(iter).le.0d0) then
          this%evnt_list(j)%unwgt=.false.
       elseif (this%f_max(iter).eq.0d0) then
          this%evnt_list(j)%unwgt=.true.
       else
          this%evnt_list(j)%unwgt=log(this%evnt_list(j)%f_abs(iter))-log(this%f_max(iter)).gt.&
               log(this%evnt_list(j)%rnd)
       endif
       if (this%evnt_list(j)%unwgt) then
            this%nevts_unw_gen=this%nevts_unw_gen+1
          this%evnt_list(j)%overwgt=1d0
       else
          this%evnt_list(j)%overwgt=0d0
       endif
    enddo
  end subroutine integral_unwgt
  
  ! Recompute candidate-event envelopes for all active iterations and set f_max.
  subroutine integral_compute_fmax(this,thischan)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevnt,iter,ios
    logical :: rescale_valid
    integer,allocatable,dimension(:) :: index_fmax_top
    real(kind=8),allocatable,dimension(:) :: fmax_top
    real(kind=8),allocatable,dimension(:,:) :: fabs
    character(len=256) :: allocation_message
    nevnt=this%nevnt_in_list
    if (nevnt.eq.0) return
    allocate(fabs(nevnt,this%current_iter),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate candidate-event envelope workspace',&
            nevnt,this%current_iter,trim(allocation_message)
       stop 1
    endif
    fabs=0d0
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       do k=this%warmup_iterations+1,this%current_iter
          if ( (k.eq.iter .and. iter.ne.this%current_iter) .or. &
               (k.ne.iter .and. iter.eq.this%current_iter) ) then
             call thischan%recompute_wgt_from_x(k,this%adaptation_class,x,wgt_new)
             call rescale_event_value(this%evnt_list(j)%f_abs(iter),wgt_new,wgt,&
                  this%evnt_list(j)%f_abs(k),rescale_valid)
             if (.not.rescale_valid) then
                write (*,*) 'ERROR: invalid candidate-event grid reweighting',j,k,iter,wgt_new,wgt
                stop 1
             endif
          endif
          fabs(j,k)=this%evnt_list(j)%f_abs(k)
       enddo
    enddo
    nevnt=max(int(write_evnt_fraction*nevnt),1)
    allocate(fmax_top(nevnt),index_fmax_top(nevnt),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate candidate-event selection workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    do k=this%warmup_iterations+1,this%current_iter
       call topk_largest(fabs(:,k),nevnt,fmax_top,index_fmax_top)
       this%f_max(k)=fmax_top(nevnt)
    enddo
    deallocate(fmax_top)
    deallocate(index_fmax_top)
    deallocate(fabs)
  end subroutine integral_compute_fmax

  ! Predict the next iteration's f_max from stored candidate events.
  subroutine integral_compute_fmax_next_iter(this,thischan)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,nevnt,iter,next_iter,ios
    logical :: rescale_valid
    integer,allocatable,dimension(:) :: index_fmax_top
    real(kind=8),allocatable,dimension(:) :: fmax_top,fabs
    character(len=256) :: allocation_message
    next_iter=this%current_iter
    nevnt=this%nevnt_in_list
    if (nevnt.le.200) then
       this%f_max(next_iter)=this%f_max(next_iter-1)
       return
    endif
    allocate(fabs(nevnt),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate next-iteration envelope workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       call thischan%recompute_wgt_from_x(next_iter,this%adaptation_class,x,wgt_new)
       call rescale_event_value(this%evnt_list(j)%f_abs(iter),wgt_new,wgt,&
            this%evnt_list(j)%f_abs(next_iter),rescale_valid)
       if (.not.rescale_valid) then
          write (*,*) 'ERROR: invalid next-iteration candidate reweighting',j,next_iter,iter,wgt_new,wgt
          stop 1
       endif
       fabs(j)=this%evnt_list(j)%f_abs(next_iter)
    enddo
    nevnt=max(int(write_evnt_fraction*nevnt),1)
    nevnt=max(int(nevnt*dble(thischan%max_iters-this%current_iter)/dble(thischan%max_iters)),1)
    allocate(fmax_top(nevnt),index_fmax_top(nevnt),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate next-iteration selection workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    call topk_largest(fabs,nevnt,fmax_top,index_fmax_top)
    this%f_max(next_iter)=fmax_top(nevnt)
    deallocate(fabs)
    deallocate(fmax_top)
    deallocate(index_fmax_top)
  end subroutine integral_compute_fmax_next_iter
  
  ! Re-evaluate the adaptive-grid weight for an existing point in a given
  ! channel iteration.
  subroutine recompute_wgt_from_x(this,iter,adaptation_class,x,wgt)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: iter,adaptation_class
    real(kind=8),dimension(this%ndim),intent(in) :: x
    real(kind=8),intent(out) :: wgt
    integer :: i
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i,iter,adaptation_class)%get_wgt(x(i),wgt)
    enddo
  end subroutine recompute_wgt_from_x

  ! Multiply wgt by the Jacobian contribution for x in this grid.
  subroutine get_wgt(this,x,wgt)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x
    real(kind=8),intent(inout) :: wgt
    real(kind=8) :: dx,factor
    integer :: cell
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(wgt)) then
       wgt=0d0
       return
    endif
    if (x.lt.0d0 .or. x.gt.1d0 .or. wgt.le.0d0) then
       wgt=0d0
       return
    endif
    call this%find_cell(x,cell)
    dx=this%current(cell)-this%current(cell-1)
    factor=dx*dble(this%size)
    if (.not.ieee_is_finite(factor)) then
       wgt=0d0
       return
    endif
    if (factor.le.0d0 .or. factor.gt.integrator_value_limit) then
       wgt=0d0
       return
    endif
    if (wgt.gt.integrator_value_limit/factor) then
       wgt=0d0
       return
    endif
    wgt=wgt*factor
  end subroutine get_wgt

  ! Locate the interpolation cell containing x in the full grid.
  subroutine find_cell(this,x,cell)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x
    integer,intent(out) :: cell
    integer :: lo,hi,mid
    lo=0
    hi=this%size
    do
       mid=(lo+hi)/2
       if (x.lt.this%current(mid)) then
          hi=mid
       else
          lo=mid+1
       end if
       if (lo.ge.hi) exit
    enddo
    cell=min(max(lo,1),this%size)
  end subroutine find_cell
  
  ! Locate the adaptation fill cell containing x.
  subroutine find_cell_to_fill(this,x,cell)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x
    integer,intent(out) :: cell
    integer :: lo,hi,mid
    lo=0
    hi=this%size_fill
    do
       mid=(lo+hi)/2
       if (x.lt.this%current_for_fillcell(mid)) then
          hi=mid
       else
          lo=mid+1
       end if
       if (lo.ge.hi) exit
    enddo
    cell=min(max(lo,1),this%size_fill)
  end subroutine find_cell_to_fill
  
  ! Print accumulated channel and integral results to the log unit.
  subroutine channel_print_combined_result(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    real(kind=8) :: rel_unc
    if (this%res(1).gt.0d0) then
       rel_unc=this%unc(1)/this%res(1)*100d0
    else
       rel_unc=0d0
    endif
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS (accum):',this%res(1),'+/-',this%unc(1),'(',rel_unc,'%)'
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel     (accum):',this%res(2),'+/-',this%unc(2),'(',rel_unc,'%)'
    do i=1,this%nintegral
       this%integrals(i)%npoints_nonzero_total=checked_point_sum(&
            this%integrals(i)%npoints_nonzero_total,&
            this%integrals(i)%npoints_nonzero,'accumulated non-zero point count')
       if (.not.this%integration_only .and. &
            this%integrals(i)%nevts_unw_gen.ge.this%integrals(i)%nevts_unw_req) then
          write(99,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,f8.6,1x,a,1x,i10,1x,a)') &
               i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
               '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevnt_in_list,&
               '--',this%integrals(i)%overweight,'--',this%integrals(i)%nevts_unw_req,'-- DONE'
       else
          if (this%integrals(i)%nevnt_in_list.lt.this%integrals(i)%nevts_unw_req) then
             write(99,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,a,1x,a,1x,i10)') &
                  i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
                  '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevnt_in_list,&
                  '--','    N/A ','--',this%integrals(i)%nevts_unw_req

          else
             write(99,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,f8.6,1x,a,1x,i10)') &
                  i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
                  '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevnt_in_list,&
                  '--',this%integrals(i)%overweight,'--',this%integrals(i)%nevts_unw_req
          endif
       endif
    enddo
  end subroutine channel_print_combined_result
    
  ! Print the current iteration-only result for one active channel.
  subroutine channel_print_result_iter(this)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8) :: rel_unc
    if (this%res_iter(1).gt.0d0) then
       rel_unc=this%unc_iter(1)/this%res_iter(1)*100d0
    else
       rel_unc=0d0
    endif
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS:',this%res_iter(1),'+/-',this%unc_iter(1),'(',rel_unc,'%)'
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel    :',this%res_iter(2),'+/-',this%unc_iter(2),'(',rel_unc,'%)'
  end subroutine channel_print_result_iter

  ! Combine all integral estimates inside one channel.
  subroutine channel_combine_iters(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    do i=1,this%nintegral
       if (.not. this%integrals(i)%evgen_done) &
            call this%integrals(i)%combine_iters(this%current_iter)
    enddo
    do i=1,2
       this%res(i)=sum(this%integrals(1:this%nintegral)%res(i))
       this%unc(i)=sqrt(sum(this%integrals(1:this%nintegral)%unc(i)**2))
    enddo
    this%aux_res=0d0
    this%aux_unc=0d0
    do i=1,this%nintegral
       this%aux_res=this%aux_res+this%integrals(i)%aux_res
       this%aux_unc=this%aux_unc+this%integrals(i)%aux_unc**2
    enddo
    this%aux_unc=sqrt(this%aux_unc)
    if (this%current_iter.eq.1) then
       this%npoints=this%npoints_iter
    else
       this%npoints=checked_point_sum(this%npoints,this%npoints_iter,&
            'accumulated channel point count')
    endif
  end subroutine channel_combine_iters

  ! Combine a new iteration estimate into one integral's accumulated estimate.
  subroutine integral_combine_iters(this,iter)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: iter
    integer :: i
    if (iter.eq.1) then
       this%res=this%res_iter
       this%unc=this%unc_iter
       this%npoints=this%npoints_iter
    else
       do i=1,2
          call update_res_and_unc(this%res(i),this%unc(i),this%npoints,this%res_iter(i),this%unc_iter(i),this%npoints_iter)
       enddo
    endif
    if (iter.eq.1) then
       this%aux_res=this%aux_res_iter
       this%aux_unc=this%aux_unc_iter
    else
       do i=1,size(this%aux_res)
          call update_res_and_unc(this%aux_res(i),this%aux_unc(i),this%npoints,&
               this%aux_res_iter(i),this%aux_unc_iter(i),this%npoints_iter)
       enddo
    endif
    if (iter.ne.1) this%npoints=checked_point_sum(this%npoints,&
         this%npoints_iter,'accumulated integral point count')
  end subroutine integral_combine_iters
  
  ! Combine two independent sample means and their standard errors.
  subroutine update_res_and_unc(res,unc,npoints,res_iter,unc_iter,npoints_iter)
    implicit none
    real(kind=8),intent(inout) :: res,unc
    real(kind=8),intent(in) :: res_iter,unc_iter
    integer(kind=8),intent(in) :: npoints,npoints_iter
    integer(kind=8) :: np
    real(kind=8) :: old_fraction,new_fraction,variance,mean_difference
    if (npoints.lt.1_8 .or. npoints_iter.lt.1_8) then
       write (*,*) 'ERROR: cannot combine empty integration samples:',npoints,npoints_iter
       stop 1
    endif
    if (npoints.gt.huge(np)-npoints_iter) then
       write (*,*) 'ERROR: accumulated integration point count overflows 64-bit integer'
       stop 1
    endif
    np=npoints+npoints_iter
    old_fraction=dble(npoints)/dble(np)
    new_fraction=dble(npoints_iter)/dble(np)
    mean_difference=res-res_iter
    variance=(old_fraction*unc)**2+(new_fraction*unc_iter)**2+&
         old_fraction*new_fraction*mean_difference**2/dble(np)
    if (.not.ieee_is_finite(variance)) then
       write (*,*) 'ERROR: invalid combined integration variance:',variance
       stop 1
    endif
    if (variance.lt.0d0) then
       write (*,*) 'ERROR: invalid combined integration variance:',variance
       stop 1
    endif
    unc=sqrt(variance)
    res=old_fraction*res+new_fraction*res_iter
    if (.not.ieee_is_finite(res) .or. .not.ieee_is_finite(unc)) then
       write (*,*) 'ERROR: invalid combined integration result:',res,unc
       stop 1
    endif
  end subroutine update_res_and_unc
  
  ! Compute the current iteration estimate for every integral in a channel.
  subroutine channel_update_result_iter(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    do i=1,this%nintegral
       call this%integrals(i)%update_result_iter()
    enddo
    do i=1,2
       this%res_iter(i)=sum(this%integrals(1:this%nintegral)%res_iter(i))
       this%unc_iter(i)=sqrt(sum(this%integrals(1:this%nintegral)%unc_iter(i)**2))
    enddo
    this%aux_res_iter=0d0
    this%aux_unc_iter=0d0
    do i=1,this%nintegral
       this%aux_res_iter=this%aux_res_iter+this%integrals(i)%aux_res_iter
       this%aux_unc_iter=this%aux_unc_iter+this%integrals(i)%aux_unc_iter**2
    enddo
    this%aux_unc_iter=sqrt(this%aux_unc_iter)
  end subroutine channel_update_result_iter

  ! Compute the standard error of the mean from first and second moments.
  subroutine compute_uncertainty(acc,acc2,np,unc)
    implicit none
    real(kind=8),intent(in) :: acc,acc2
    integer(kind=8),intent(in) :: np
    real(kind=8),intent(out) :: unc
    unc=sqrt(abs(acc2-acc**2)/dble(np))
  end subroutine compute_uncertainty
  
  ! Convert one integral's accumulated sums into an iteration estimate.
  subroutine integral_update_result_iter(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: i
    if (this%npoints_iter.ne.0_8) then
       this%res_iter=this%accum/dble(this%npoints_iter)
       this%res2_iter=this%accum2/dble(this%npoints_iter)
       do i=1,2
          call compute_uncertainty(this%res_iter(i),this%res2_iter(i),this%npoints_iter,this%unc_iter(i))
       enddo
       this%aux_res_iter=this%aux_accum/dble(this%npoints_iter)
       this%aux_res2_iter=this%aux_accum2/dble(this%npoints_iter)
       do i=1,size(this%aux_res_iter)
          call compute_uncertainty(this%aux_res_iter(i),this%aux_res2_iter(i),this%npoints_iter,this%aux_unc_iter(i))
       enddo
    endif
  end subroutine integral_update_result_iter
  
  ! Add one evaluated point to a channel grid and to its active integral.
  subroutine channel_add_point(this,x,wgt,f_abs,f,to_write,accepted,event_label,&
       candidate_event_limit,f_aux)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(in) :: x
    real(kind=8),intent(in) :: f_abs,f,wgt
    logical,intent(out) :: to_write
    logical,intent(in) :: accepted
    integer,intent(inout) :: event_label
    integer,intent(in) :: candidate_event_limit
    real(kind=8),dimension(:),intent(in),optional :: f_aux
    integer :: i
    this%npoints_iter=checked_point_sum(this%npoints_iter,1_8,&
         'channel iteration point count')
    do i=1,this%ndim
       if (this%integration_only) then
          ! For a signed integration the variance-optimal density is
          ! proportional to abs(f), after real--dipole cancellation.
          call this%grids(i,this%current_iter,this%adaptation_class(this%current_integral))%add_point(&
               x(i),abs(f),.true.)
       else
          call this%grids(i,this%current_iter,this%adaptation_class(this%current_integral))%add_point(x(i),f_abs)
       endif
    enddo
    if (present(f_aux)) then
       call this%integrals(this%current_integral)%add_point(x,wgt,f_abs,f,to_write,&
            accepted,event_label,candidate_event_limit,f_aux)
    else
       call this%integrals(this%current_integral)%add_point(x,wgt,f_abs,f,to_write,&
            accepted,event_label,candidate_event_limit)
    endif
    if (all(this%integrals%done)) this%done=.true.
  end subroutine channel_add_point

  ! Accumulate one point in an integral and optionally save it as a candidate
  ! event for later final unweighting.
  subroutine integral_add_point(this,x,wgt,f_abs,f,to_write,accepted,event_label,&
       candidate_event_limit,f_aux)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,f,wgt
    real(kind=8),dimension(this%ndim),intent(in) :: x
    logical,intent(out) :: to_write
    logical,intent(in) :: accepted
    integer,intent(inout) :: event_label
    integer,intent(in) :: candidate_event_limit
    real(kind=8),dimension(:),intent(in),optional :: f_aux
    logical :: enough
    this%npoints_iter=checked_point_sum(this%npoints_iter,1_8,&
         'integral iteration point count')
    if (accepted) this%npoints_nonzero=checked_point_sum(&
         this%npoints_nonzero,1_8,'integral non-zero point count')
    this%accum(1)=this%accum(1)+f_abs
    this%accum(2)=this%accum(2)+f
    this%accum2(1)=this%accum2(1)+f_abs**2
    this%accum2(2)=this%accum2(2)+f**2
    if (present(f_aux)) then
       this%aux_accum=this%aux_accum+f_aux
       this%aux_accum2=this%aux_accum2+f_aux**2
    endif
    if (this%current_iter.le.this%warmup_iterations) call this%update_max_value(f_abs)
    call this%check_write_evnt(x,wgt,f_abs,to_write,enough,event_label,candidate_event_limit)
    if (this%integration_only) then
       if (this%npoints_iter.ge.this%npoints_requested) this%done=.true.
    elseif (this%npoints_iter.ge.this%npoints_requested .or. enough) then
       ! Zero-valued and cut points are valid Monte Carlo observations.  Counting
       ! only nonzero points makes an identically zero leaf loop forever.
       this%done=.true.
    endif
  end subroutine integral_add_point

  ! Track the largest absolute integrand value seen before event generation.
  subroutine update_max_value(this,f_abs)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs
    this%max_value=max(this%max_value,f_abs)
  end subroutine update_max_value
  
  ! Decide whether a point should be written as a candidate event.
  subroutine check_write_evnt(this,x,wgt,f_abs,to_write,enough,event_label,candidate_event_limit)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,wgt
    real(kind=8),dimension(this%ndim),intent(in) :: x
    logical,intent(out) :: to_write,enough
    integer,intent(inout) :: event_label
    integer,intent(in) :: candidate_event_limit
    real(kind=8) :: rnd
    integer :: ios
    character(len=256) :: allocation_message
    to_write=.false.
    enough=.false.
    if (this%integration_only) return
    if (this%current_iter.le.this%warmup_iterations) return
    rnd=max(validated_uniform_random('candidate-event selection'),event_random_floor)
    if (f_abs.gt.this%f_max(this%current_iter)*rnd) then
       to_write=.true.
       if (event_label.ge.candidate_event_limit) then
          write (*,*) 'ERROR: retained candidate events exceed the supported workspace:',&
               event_label,candidate_event_limit
          stop 1
       endif
       if (event_label.eq.huge(event_label)) then
          write (*,*) 'ERROR: candidate-event label counter overflow'
          stop 1
       endif
       event_label=event_label+1
       this%evnt=this%evnt+1
       this%nevnt_in_list=this%nevnt_in_list+1
       if (this%nevnt_in_list.gt.size(this%evnt_list)) &
            call this%increase_size_evnt_list(candidate_event_limit)
       allocate(this%evnt_list(this%nevnt_in_list)%f_abs(this%max_iters),&
            this%evnt_list(this%nevnt_in_list)%x(this%ndim),&
            stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'ERROR: could not allocate retained candidate event: ',&
               trim(allocation_message)
          stop 1
       endif
       this%evnt_list(this%nevnt_in_list)%f_abs=0d0
       this%evnt_list(this%nevnt_in_list)%f_abs(this%current_iter)=f_abs
       this%evnt_list(this%nevnt_in_list)%rnd=rnd
       this%evnt_list(this%nevnt_in_list)%wgt=wgt
       this%evnt_list(this%nevnt_in_list)%x=x
       this%evnt_list(this%nevnt_in_list)%iter=this%current_iter
       this%evnt_list(this%nevnt_in_list)%label=event_label
       this%evnt_list(this%nevnt_in_list)%unwgt=.false.
       this%evnt_list(this%nevnt_in_list)%overwgt=0d0
       this%nevts_unw_gen=this%nevts_unw_gen+1
       if (this%nevts_unw_gen.gt.1.5d0*this%nevts_unw_req) enough=.true.
    endif
  end subroutine check_write_evnt

  ! Return final per-event weights for all candidate events written by the
  ! caller during the get_points/fill_points loop.
  subroutine assign_evnt_wgts(this,wgts)
    implicit none
    class(integrator) :: this
    real(kind=8),allocatable,dimension(:,:),intent(out) :: wgts
    integer :: i,j,selected_events,ios
    real(kind=8) :: nominal_wgt
    character(len=256) :: allocation_message
    allocate(wgts(3,this%event_label),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate final event weights: ',trim(allocation_message)
       stop 1
    endif
    wgts=0d0
    selected_events=0
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          if (this%channels(i)%integrals(j)%nevnt_in_list.gt.0) &
               selected_events=selected_events+count(&
               this%channels(i)%integrals(j)%evnt_list(&
               1:this%channels(i)%integrals(j)%nevnt_in_list)%unwgt)
       enddo
    enddo
    if (selected_events.ne.this%nevts_unw_req) then
       write (*,*) 'ERROR: event generation ended without the requested selected sample:',&
            selected_events,this%nevts_unw_req
       stop 1
    endif
    nominal_wgt=this%res(1)
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          call this%channels(i)%integrals(j)%compute_wgts(nominal_wgt,wgts)
       enddo
    enddo
  end subroutine assign_evnt_wgts

  ! Fill the final weight columns for candidate events belonging to one integral.
  subroutine compute_wgts(this,nominal_wgt,wgts)
    implicit none
    class(integral) :: this
    real(kind=8),dimension(:,:),intent(inout) :: wgts
    real(kind=8),intent(in) :: nominal_wgt
    real(kind=8) :: number_of_evnts,number_of_wgts
    integer :: i
    number_of_evnts=0d0
    number_of_wgts=0d0
    do i=1,this%nevnt_in_list
       if (this%evnt_list(i)%unwgt) then
          number_of_evnts=number_of_evnts+1d0
          number_of_wgts=number_of_wgts+max(1d0,this%evnt_list(i)%overwgt)
       endif
    enddo
    do i=1,this%nevnt_in_list
       if (this%evnt_list(i)%unwgt) then
          wgts(1,this%evnt_list(i)%label)=nominal_wgt
          wgts(2,this%evnt_list(i)%label)=nominal_wgt*max(1d0,this%evnt_list(i)%overwgt) &
               *number_of_evnts/number_of_wgts
          wgts(3,this%evnt_list(i)%label)=max(0d0,this%evnt_list(i)%overwgt-1d0)
       else
          wgts(1:3,this%evnt_list(i)%label)=0d0
       endif
    enddo
  end subroutine compute_wgts

  ! Double the candidate-event buffer while preserving existing events.
  subroutine increase_size_evnt_list(this,candidate_event_limit)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: candidate_event_limit
    type(evnt),allocatable,dimension(:) :: tmp_list
    integer :: isize,new_size,ios
    character(len=256) :: allocation_message
    isize=size(this%evnt_list)
    if (isize.ge.candidate_event_limit) then
       write (*,*) 'ERROR: candidate-event buffer reached the supported workspace:',&
            isize,candidate_event_limit
       stop 1
    endif
    if (isize.gt.shiftr(huge(isize),1)) then
       new_size=candidate_event_limit
    else
       new_size=min(2*isize,candidate_event_limit)
    endif
    if (new_size.le.isize) then
       write (*,*) 'ERROR: candidate-event buffer cannot be enlarged safely:',isize,new_size
       stop 1
    endif
    allocate(tmp_list(new_size),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not enlarge candidate-event buffer: ',&
            trim(allocation_message)
       stop 1
    endif
    tmp_list%unwgt=.false.
    tmp_list%overwgt=0d0
    tmp_list(1:isize)=this%evnt_list(1:isize)
    deallocate(this%evnt_list)
    this%evnt_list=tmp_list
  end subroutine increase_size_evnt_list
  
  ! Accumulate one sampled point into the grid-adaptation histogram.
  subroutine grid_add_point(this,x,f_abs,use_mean_abs)
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x,f_abs
    logical,intent(in),optional :: use_mean_abs
    integer :: cell
    logical :: mean_abs_metric
    mean_abs_metric=.false.
    if (present(use_mean_abs)) mean_abs_metric=use_mean_abs
    call this%find_cell_to_fill(x,cell)
    if (mean_abs_metric) then
       ! The accuracy-only metric is a cell mean of abs(f), rather than the
       ! peak-oriented event envelope used by the strategies below.
       this%use_mean_abs=.true.
       this%accum(cell)=this%accum(cell)+f_abs
    elseif (importance_sampling_strategy.eq.1) then
       this%accum(cell)=this%accum(cell)+f_abs
    elseif (importance_sampling_strategy.eq.2) then
       this%accum(cell)=max(this%accum(cell),f_abs)
    elseif (importance_sampling_strategy.eq.3) then
       if (f_abs.gt.this%accum(cell)) then
          this%accum(cell)=this%accum(cell)+(f_abs-this%accum(cell))*0.1d0
       endif
    elseif (importance_sampling_strategy.eq.4) then
       if (this%accum(cell).le.0d0) then
          this%accum(cell)=f_abs*1d-4
       elseif (f_abs.gt.this%accum(cell)) then
          this%accum(cell)=this%accum(cell)*1.1d0
       endif
    endif
    this%nhits(cell)=checked_point_sum(this%nhits(cell),1_8,&
         'adaptive-grid cell hit count')
  end subroutine grid_add_point

  ! Convert accumulated grid information into a new monotone grid.
  subroutine grid_update(this,npoints,new_grid)
    implicit none
    class(grid),intent(inout) :: this
    integer(kind=8),intent(in) :: npoints
    class(grid),intent(out) :: new_grid
    real(kind=8),dimension(0:this%size_fill) :: current
    integer :: i,j
    real(kind=8) :: r
    call this%massage_accum()
    current(0)=0d0
    do i=1,this%size_fill
       r=dble(i)/dble(this%size_fill)
       do j=1,this%size_fill
          if (r.lt.this%accum(j)) then
             current(i)=this%current_for_fillcell(j-1)+(r-this%accum(j-1))/ &
                  (this%accum(j)-this%accum(j-1))*(this%current_for_fillcell(j)-this%current_for_fillcell(j-1))
             exit
          endif
       enddo
    enddo
    deallocate(this%accum)
    deallocate(this%nhits)
    current(this%size_fill)=1d0
    call new_grid%init(npoints,current)
  end subroutine grid_update
  
  ! Resize a monotone grid with shape-preserving interpolation.
  subroutine interpolate_current(this,size_in,size_out,current_in,current_out)
    use pchip_uniform_strict
    implicit none
    class(grid),intent(inout) :: this
    integer,intent(in) :: size_in,size_out
    real(kind=8),dimension(0:size_in),intent(in) :: current_in
    real(kind=8),dimension(0:size_out),intent(out) :: current_out
    call resize_arr_pchip_strict(current_in,size_out,current_out)
  end subroutine interpolate_current
  
  ! Smooth and normalise the adaptation histogram into a cumulative map.
  subroutine massage_accum(this)
    implicit none
    class(grid),intent(inout) :: this
    integer :: i
    real(kind=8) :: total
    real(kind=8), parameter :: tiny=1d-8
    do i=1,this%size_fill
       if (this%nhits(i).eq.0_8) cycle
       if (this%use_mean_abs .or. importance_sampling_strategy.eq.1) then
          this%accum(i)=this%accum(i)/this%nhits(i)
       else
          this%accum(i)=this%accum(i)
       endif
    enddo
    total=sum(this%accum)
    if (.not.all(ieee_is_finite(this%accum))) then
       call set_uniform_accumulation()
       return
    endif
    if (.not.ieee_is_finite(total)) then
       call set_uniform_accumulation()
       return
    endif
    if (any(this%accum.lt.0d0) .or. total.le.0d0) then
       ! A zero-rate channel contains no information with which to deform its
       ! grid.  Preserve a valid uniform cumulative map instead of normalising
       ! 0/0 (or propagating a damaged adaptation histogram).
       call set_uniform_accumulation()
       return
    endif
    do i=1,this%size_fill
       if (this%accum(i).lt.1d-12*total) then
          this%accum(i)=0d0
       elseif (this%accum(i).lt.(1d0-1d-12)*total) then
          this%accum(i)=((this%accum(i)/total-1d0)/log(this%accum(i)/total))**1.5d0
       else
          this%accum(i)=1d0
       endif
       this%accum(i)=this%accum(i-1)+max(this%accum(i),0d0)
    enddo
    this%accum=this%accum/this%accum(this%size_fill)
    ! make sure the elements are at least 'tiny' apart
    do i=1,this%size_fill
       if (this%accum(i).lt.this%accum(i-1)+tiny) then
          this%accum(i)=this%accum(i-1)+tiny
       endif
    enddo
    this%accum=this%accum/this%accum(this%size_fill)
  contains
    subroutine set_uniform_accumulation()
      implicit none
      integer :: j
      this%accum(0)=0d0
      do j=1,this%size_fill
         this%accum(j)=dble(j)/dble(this%size_fill)
      enddo
    end subroutine set_uniform_accumulation
  end subroutine massage_accum

  ! Generate one full point for a channel, including flat extra coordinates.
  subroutine channel_get_point(this,x,wgt)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim+this%ndim_extra),intent(out) :: x
    real(kind=8),intent(out) :: wgt
    integer :: i,iadapt
    real(kind=8) :: uniform
    iadapt=this%adaptation_class(this%current_integral)
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i,this%current_iter,iadapt)%get_x(x(i),wgt)
    enddo
    do i=this%ndim+1,this%ndim+this%ndim_extra
       uniform=validated_uniform_random('flat-coordinate generation')
       x(i)=uniform
    enddo
  end subroutine channel_get_point

  ! Generate one adapted coordinate and multiply by its Jacobian.
  subroutine get_x(this,x,wgt)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(out) :: x
    integer :: cell
    real(kind=8),intent(inout) :: wgt
    real(kind=8) :: rnd,dx,uniform,factor
    uniform=validated_uniform_random('adaptive-coordinate generation')
    if (.not.ieee_is_finite(wgt)) then
       x=0.5d0
       wgt=0d0
       return
    endif
    if (wgt.le.0d0) then
       x=0.5d0
       wgt=0d0
       return
    endif
    rnd=this%size*uniform
    cell=int(rnd)+1
    cell=min(max(cell,1),this%size)
    rnd=rnd-dble(cell-1)
    dx=this%current(cell)-this%current(cell-1)
    x=this%current(cell-1)+rnd*dx
    factor=dx*dble(this%size)
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(factor)) then
       x=0.5d0
       wgt=0d0
       return
    endif
    if (x.lt.0d0 .or. x.gt.1d0 .or. factor.le.0d0 .or. &
         factor.gt.integrator_value_limit) then
       x=0.5d0
       wgt=0d0
       return
    endif
    if (wgt.gt.integrator_value_limit/factor) then
       x=0.5d0
       wgt=0d0
       return
    endif
    wgt=wgt*factor
  end subroutine get_x

  pure subroutine rescale_event_value(value,new_weight,old_weight,result,valid)
    real(kind=8),intent(in) :: value,new_weight,old_weight
    real(kind=8),intent(out) :: result
    logical,intent(out) :: valid
    real(kind=8) :: logarithm
    result=0d0
    valid=.false.
    if (.not.ieee_is_finite(value) .or. .not.ieee_is_finite(new_weight) .or. &
         .not.ieee_is_finite(old_weight)) return
    if (value.lt.0d0 .or. new_weight.le.0d0 .or. old_weight.le.0d0) return
    if (value.eq.0d0) then
       valid=.true.
       return
    endif
    logarithm=log(value)+log(new_weight)-log(old_weight)
    if (.not.ieee_is_finite(logarithm)) return
    if (logarithm.gt.log(integrator_value_limit)) return
    if (logarithm.lt.log(tiny(1d0))) then
       valid=.true.
       return
    endif
    result=exp(logarithm)
    valid=ieee_is_finite(result)
    if (valid) valid=result.le.integrator_value_limit
    if (.not.valid) result=0d0
  end subroutine rescale_event_value

  ! Randomly select an unfinished channel/integral pair for the next batch.
  subroutine get_channel_and_integral(this,ichan,iint,wgt_chan)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(out) :: ichan,iint
    real(kind=8),intent(out) :: wgt_chan
    if (this%integration_only) then
       call select_accuracy_leaf(this,ichan,iint)
    else
       do
          ichan=int(validated_uniform_random('channel selection')*this%nchannel)+1
          if (.not.(this%channels(ichan)%done.or.this%channels(ichan)%evgen_done)) exit
       enddo
    endif
    wgt_chan=1d0!dble(this%nchannel)
    if (.not.this%integration_only) then
       do
          iint=int(validated_uniform_random('integral selection')*this%channels(ichan)%nintegral)+1
          if (.not.(this%channels(ichan)%integrals(iint)%done.or.&
               this%channels(ichan)%integrals(iint)%evgen_done)) exit
       enddo
    endif
    wgt_chan=wgt_chan*1d0!dble(this%channels(ichan)%nintegral)
    this%channels(ichan)%current_integral=iint
  end subroutine get_channel_and_integral

  ! Deterministically choose the least-filled accuracy-only leaf.  This
  ! makes every evaluated point count against its exact leaf quota; ties are
  ! rotated in flattened channel/integral order for reproducibility.
  subroutine select_accuracy_leaf(this,ichan,iint)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(out) :: ichan,iint
    integer :: i,j,k,nleaf,chosen,best_rank,rank
    real(kind=8) :: fraction,best_fraction

    nleaf=0
    do i=1,this%nchannel
       nleaf=nleaf+this%channels(i)%nintegral
    enddo
    ichan=0
    iint=0
    chosen=0
    best_fraction=huge(1d0)
    best_rank=nleaf
    k=0
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          k=k+1
          if (this%channels(i)%integrals(j)%done.or.this%channels(i)%integrals(j)%evgen_done) cycle
          fraction=dble(this%channels(i)%integrals(j)%npoints_iter)/&
               dble(max(this%channels(i)%integrals(j)%npoints_requested,1_8))
          rank=modulo(k-this%allocation_cursor-1,nleaf)
          if (fraction.lt.best_fraction-1d-14 .or. &
               (abs(fraction-best_fraction).le.1d-14 .and. rank.lt.best_rank)) then
             ichan=i
             iint=j
             chosen=k
             best_fraction=fraction
             best_rank=rank
          endif
       enddo
    enddo
    if (chosen.eq.0) then
       write (*,*) 'ERROR: no unfinished accuracy-only integration leaf'
       stop 1
    endif
    this%allocation_cursor=chosen
  end subroutine select_accuracy_leaf
  
  ! Public helper: compute the current adaptive-grid weight for an existing
  ! point in channel ichan.
  subroutine compute_wgt_from_x(this,ichan,x,wgt,iint,adaptation_class)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: ichan
    integer,intent(in),optional :: iint,adaptation_class
    real(kind=8),dimension(:),intent(in) :: x
    real(kind=8),intent(out) :: wgt
    integer :: iadapt
    wgt=0d0
    if (.not.allocated(this%channels)) then
       write (*,*) 'ERROR: compute_wgt_from_x called before integrator initialisation'
       stop 1
    endif
    if (ichan.lt.1 .or. ichan.gt.this%nchannel) then
       write (*,*) 'ERROR: invalid channel in compute_wgt_from_x',ichan,this%nchannel
       stop 1
    endif
    if (size(x).ne.this%channels(ichan)%ndim) then
       write (*,*) 'ERROR: invalid coordinate extent in compute_wgt_from_x',&
            size(x),this%channels(ichan)%ndim
       stop 1
    endif
    if (present(iint) .and. present(adaptation_class)) then
       write (*,*) 'ERROR: compute_wgt_from_x accepts either an integral or an adaptation class'
       stop 1
    endif
    iadapt=1
    if (present(iint)) then
       if (iint.lt.1 .or. iint.gt.this%channels(ichan)%nintegral) then
          write (*,*) 'ERROR: invalid integral in compute_wgt_from_x',iint,&
               this%channels(ichan)%nintegral
          stop 1
       endif
       iadapt=this%channels(ichan)%adaptation_class(iint)
    endif
    if (present(adaptation_class)) iadapt=adaptation_class
    if (iadapt.lt.1 .or. iadapt.gt.this%channels(ichan)%nadaptation) then
       write (*,*) 'ERROR: invalid adaptation class in compute_wgt_from_x',iadapt,&
            this%channels(ichan)%nadaptation
       stop 1
    endif
    if (this%channels(ichan)%current_iter.lt.1 .or. &
         this%channels(ichan)%current_iter.gt.this%channels(ichan)%max_iters) then
       write (*,*) 'ERROR: invalid grid iteration in compute_wgt_from_x',&
            this%channels(ichan)%current_iter
       stop 1
    endif
    if (.not.all(ieee_is_finite(x))) return
    if (any(x.lt.0d0) .or. any(x.gt.1d0)) return
    call this%channels(ichan)%recompute_wgt_from_x(this%channels(ichan)%current_iter,iadapt,x,wgt)
  end subroutine compute_wgt_from_x

  ! Return the signed and absolute combined estimates for every channel.
  subroutine get_channel_results(this,res,unc,success)
    implicit none
    class(integrator),intent(in) :: this
    real(kind=8),allocatable,dimension(:,:),intent(out) :: res,unc
    logical,intent(out),optional :: success
    integer :: i,ios
    character(len=256) :: allocation_message
    if (present(success)) success=.false.
    if (.not.allocated(this%channels) .or. this%nchannel.lt.1) then
       allocate(res(2,0),unc(2,0))
       return
    endif
    if (size(this%channels).ne.this%nchannel) then
       allocate(res(2,0),unc(2,0))
       return
    endif
    allocate(res(2,this%nchannel),unc(2,this%nchannel),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate channel results: ',trim(allocation_message)
       stop 1
    endif
    do i=1,this%nchannel
       res(:,i)=this%channels(i)%res
       unc(:,i)=this%channels(i)%unc
    enddo
    if (present(success)) success=.true.
  end subroutine get_channel_results

  ! Return optional signed auxiliary estimates accumulated alongside the
  ! primary scalar result.  These observables share the sampled points but do
  ! not affect adaptation or event generation.
  subroutine get_channel_aux_results(this,res,unc,success)
    implicit none
    class(integrator),intent(in) :: this
    real(kind=8),allocatable,dimension(:,:),intent(out) :: res,unc
    logical,intent(out),optional :: success
    integer :: i,naux,ios
    character(len=256) :: allocation_message
    if (present(success)) success=.false.
    if (.not.allocated(this%channels) .or. this%nchannel.lt.1) then
       allocate(res(0,0),unc(0,0))
       return
    endif
    if (size(this%channels).ne.this%nchannel) then
       allocate(res(0,0),unc(0,0))
       return
    endif
    naux=this%channels(1)%naux
    if (naux.lt.0) then
       allocate(res(0,0),unc(0,0))
       return
    endif
    do i=1,this%nchannel
       if (this%channels(i)%naux.ne.naux .or. &
            .not.allocated(this%channels(i)%aux_res) .or. &
            .not.allocated(this%channels(i)%aux_unc)) then
          allocate(res(0,0),unc(0,0))
          return
       endif
       if (size(this%channels(i)%aux_res).ne.naux .or. &
            size(this%channels(i)%aux_unc).ne.naux) then
          allocate(res(0,0),unc(0,0))
          return
       endif
    enddo
    allocate(res(naux,this%nchannel),unc(naux,this%nchannel),&
         stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: could not allocate auxiliary channel results: ',&
            trim(allocation_message)
       stop 1
    endif
    do i=1,this%nchannel
       res(:,i)=this%channels(i)%aux_res
       unc(:,i)=this%channels(i)%aux_unc
    enddo
    if (present(success)) success=.true.
  end subroutine get_channel_aux_results
  
  ! Placeholder for future grid restart support.
  subroutine read_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine read_all_grids
  
  ! Placeholder for future grid checkpoint support.
  subroutine write_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine write_all_grids
  
end module simple_integrator_mod
