program simple_integrator_multichannel_regression
  use simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: batch_size=64
  type(integrator) :: family_integrator,independent_integrator,&
       one_iter_integrator
  integer,dimension(1) :: ndim,ndim_extra,nintegral
  integer,dimension(2) :: independent_ndim,independent_ndim_extra,&
       independent_nintegral
  integer,dimension(batch_size) :: selected_map
  integer,allocatable,dimension(:) :: event_map_ids
  real(kind=8),dimension(batch_size) :: f,f_abs
  logical,dimension(batch_size) :: to_write
  real(kind=8),allocatable,dimension(:,:) :: event_wgts
  real(kind=8) :: live_probe,reference_probe,reference_snapshot
  real(kind=8) :: live_at_snapshot,tolerance,max_partition_error
  real(kind=8) :: physical_x,selected_jacobian,combined_uncertainty
  integer(kind=8),dimension(2) :: map_counts
  integer :: ichan,iint,ip,accepted,io_status
  logical :: done,reference_separated,live_changed_after_snapshot
  logical :: extended_production,restored_event_target
  logical :: saw_independent_completion
  logical :: saw_warmup_log,saw_production_log
  character(len=256) :: log_line

  open(unit=99,status='scratch',action='readwrite')

  ! One integration channel is one family. Its adapted coordinate is shared
  ! by both phase maps; a second, flat coordinate selects either map with
  ! prior probability 1/2 and is deliberately absent from the adaptive grid.
  ndim=1
  ndim_extra=1
  nintegral=1
  call family_integrator%init(1,ndim,ndim_extra,nintegral,100,7)
  if (family_integrator%nchannel.ne.1 .or.&
       family_integrator%channels(1)%ndim.ne.1 .or.&
       family_integrator%channels(1)%ndim_extra.ne.1) then
     write (*,*) 'Integration family did not create exactly one grid pair'
     stop 1
  endif

  reference_separated=.false.
  live_changed_after_snapshot=.false.
  extended_production=.false.
  restored_event_target=.false.
  reference_snapshot=0d0
  live_at_snapshot=0d0
  max_partition_error=0d0
  map_counts=0_8
  done=.false.
  do while (.not.done)
     call family_integrator%get_points(batch_size,ichan,iint)
     if (ichan.ne.1) then
        write (*,*) 'A single integration family selected another channel',ichan
        stop 1
     endif
     do ip=1,batch_size
        selected_map(ip)=min(int(2d0*family_integrator%x(2,ip))+1,2)
        map_counts(selected_map(ip))=map_counts(selected_map(ip))+1_8
        call evaluate_family_point(family_integrator,selected_map(ip),&
             family_integrator%x(1,ip),family_integrator%wgt(ip),&
             f(ip),f_abs(ip),max_partition_error)
     enddo
     call family_integrator%fill_points(batch_size,f_abs,f,to_write,done,&
          selected_map)

     ! A very efficient synthetic integrand would normally finish in its first
     ! production iteration.  Temporarily enlarge the internal event target at
     ! the warm-up boundary so that two live-grid updates occur after the
     ! reference snapshot has frozen; then restore the requested sample size.
     if (.not.extended_production .and.&
          family_integrator%channels(1)%nevts_unw_req.eq.100) then
        family_integrator%nevts_unw_req=100000000
        family_integrator%channels(1)%nevts_unw_req=100000000
        family_integrator%channels(1)%integrals(1)%nevts_unw_req=100000000
        extended_production=.true.
     elseif (extended_production .and. .not.restored_event_target .and.&
          family_integrator%channels(1)%current_iter.ge.&
          family_integrator%channels(1)%reference_grid_iter+2) then
        family_integrator%nevts_unw_req=100
        family_integrator%channels(1)%nevts_unw_req=100
        family_integrator%channels(1)%integrals(1)%nevts_unw_req=100
        restored_event_target=.true.
     endif

     call family_integrator%compute_wgt_from_x(1,[0.317d0],live_probe)
     call family_integrator%compute_reference_wgt_from_x(&
          1,[0.317d0],reference_probe)
     tolerance=2d-13*max(1d0,abs(live_probe),abs(reference_probe))
     if (.not.reference_separated) then
        if (abs(live_probe-reference_probe).gt.tolerance) then
           reference_separated=.true.
           reference_snapshot=reference_probe
           live_at_snapshot=live_probe
        endif
     else
        if (abs(reference_probe-reference_snapshot).gt.tolerance) then
           write (*,*) 'Reference grid changed during production',&
                reference_snapshot,reference_probe
           stop 1
        endif
        if (abs(live_probe-live_at_snapshot).gt.tolerance) &
             live_changed_after_snapshot=.true.
     endif
  enddo

  if (max_partition_error.gt.2d-15) then
     write (*,*) 'Family density partitions do not sum pointwise to one',&
          max_partition_error
     stop 1
  endif
  if (abs(dble(map_counts(1))/dble(sum(map_counts))-0.5d0).gt.0.06d0) then
     write (*,*) 'Submap selector is not consistent with an equal prior',map_counts
     stop 1
  endif
  if (.not.extended_production .or. .not.restored_event_target .or.&
       .not.reference_separated .or. .not.live_changed_after_snapshot) then
     write (*,*) 'Live/reference family-grid lifecycle regression failed',&
          extended_production,restored_event_target,reference_separated,&
          live_changed_after_snapshot
     stop 1
  endif
  if (.not.ieee_is_finite(family_integrator%res(2)) .or.&
       .not.ieee_is_finite(family_integrator%unc(2)) .or.&
       abs(family_integrator%res(2)-3d0).gt.&
       3d0*family_integrator%unc(2)+3d-3) then
     write (*,*) 'Synthetic family estimator disagrees with its analytic rate',&
          family_integrator%res(2),family_integrator%unc(2)
     stop 1
  endif

  call family_integrator%assign_evnt_wgts(event_wgts,event_map_ids)
  if (any(.not.ieee_is_finite(event_wgts))) then
     write (*,*) 'Stored-event family weights are non-finite'
     stop 1
  endif
  accepted=count(event_wgts(1,:).gt.0d0)
  if (accepted.ne.100) then
     write (*,*) 'Unexpected number of accepted family events',accepted
     stop 1
  endif
  if (any(event_map_ids.lt.1) .or. any(event_map_ids.gt.2)) then
     write (*,*) 'Stored event lost its selected phase-map ID'
     stop 1
  endif
  if (.not.any(event_map_ids.eq.1) .or. .not.any(event_map_ids.eq.2)) then
     write (*,*) 'Stored family events do not exercise both submaps'
     stop 1
  endif
  deallocate(event_wgts,event_map_ids)

  ! Compare the hierarchical family with two explicit proposal channels.  Each
  ! explicit channel carries half of the target integral, matching the fixed
  ! 1/2 proposal prior used above; their summed estimator must agree with the
  ! family result.  Batched channel selection also exercises independent early
  ! completion before the global iteration can finish.
  independent_ndim=[1,1]
  independent_ndim_extra=[0,0]
  independent_nintegral=[1,1]
  call independent_integrator%init(2,independent_ndim,&
       independent_ndim_extra,independent_nintegral,100,7)
  done=.false.
  saw_independent_completion=.false.
  do while (.not.done)
     call independent_integrator%get_points(batch_size,ichan,iint)
     do ip=1,batch_size
        if (ichan.eq.1) then
           physical_x=independent_integrator%x(1,ip)
           selected_jacobian=1d0
        else
           physical_x=sqrt(independent_integrator%x(1,ip))
           selected_jacobian=1d0/(2d0*physical_x)
        endif
        f(ip)=0.5d0*independent_integrator%wgt(ip)*selected_jacobian*&
             (1d0+4d0*physical_x)
        f_abs(ip)=f(ip)
     enddo
     call independent_integrator%fill_points(&
          batch_size,f_abs,f,to_write,done)
     if (any(independent_integrator%channels%done) .and.&
          .not.all(independent_integrator%channels%done)) &
          saw_independent_completion=.true.
  enddo
  combined_uncertainty=sqrt(family_integrator%unc(2)**2+&
       independent_integrator%unc(2)**2)
  if (.not.saw_independent_completion .or.&
       abs(independent_integrator%res(2)-3d0).gt.&
       3d0*independent_integrator%unc(2)+3d-3 .or.&
       abs(independent_integrator%res(2)-family_integrator%res(2)).gt.&
       3d0*combined_uncertainty+3d-3) then
     write (*,*) 'Hierarchical/independent-channel estimator mismatch',&
          saw_independent_completion,family_integrator%res(2),&
          family_integrator%unc(2),independent_integrator%res(2),&
          independent_integrator%unc(2)
     stop 1
  endif

  ! A one-iteration configuration has no warm-up: its uniform initial grid is
  ! both the reference snapshot and the sole production grid.
  ndim=1
  ndim_extra=0
  nintegral=1
  call one_iter_integrator%init(1,ndim,ndim_extra,nintegral,8,1)
  done=.false.
  do while (.not.done)
     call one_iter_integrator%get_points(batch_size,ichan,iint)
     do ip=1,batch_size
        f(ip)=2d0*one_iter_integrator%wgt(ip)
        f_abs(ip)=f(ip)
     enddo
     call one_iter_integrator%fill_points(batch_size,f_abs,f,to_write,done)
  enddo
  call one_iter_integrator%compute_wgt_from_x(1,[0.417d0],live_probe)
  call one_iter_integrator%compute_reference_wgt_from_x(&
       1,[0.417d0],reference_probe)
  if (abs(live_probe-reference_probe).gt.2d-14 .or.&
       abs(one_iter_integrator%res(2)-2d0).gt.1d-12) then
     write (*,*) 'One-iteration production/reference-grid regression failed',&
          live_probe,reference_probe,one_iter_integrator%res(2)
     stop 1
  endif

  rewind(99)
  saw_warmup_log=.false.
  saw_production_log=.false.
  do
     read(99,'(a)',iostat=io_status) log_line
     if (io_status.ne.0) exit
     if (index(log_line,'warm-up/current').gt.0) saw_warmup_log=.true.
     if (index(log_line,'production/accumulated').gt.0) &
          saw_production_log=.true.
  enddo
  if (.not.saw_warmup_log .or. .not.saw_production_log) then
     write (*,*) 'Integrator log did not label warm-up and production results'
     stop 1
  endif
  close(99)
  write (*,*) 'simple-integrator family regression: PASS'

contains

  subroutine evaluate_family_point(integrator_state,map_id,x,live_weight,&
       value,value_abs,max_error)
    type(integrator),intent(inout) :: integrator_state
    integer,intent(in) :: map_id
    real(kind=8),intent(in) :: x,live_weight
    real(kind=8),intent(out) :: value,value_abs
    real(kind=8),intent(inout) :: max_error
    real(kind=8) :: physical_x,map_one_x,map_two_x,selected_jacobian
    real(kind=8) :: reference_one,reference_two,density_one,density_two
    real(kind=8) :: alpha_one,alpha_two,partition_factor

    if (map_id.eq.1) then
       physical_x=x
       selected_jacobian=1d0
    else
       physical_x=sqrt(x)
       selected_jacobian=1d0/(2d0*physical_x)
    endif
    map_one_x=physical_x
    map_two_x=physical_x**2
    call integrator_state%compute_reference_wgt_from_x(&
         1,[map_one_x],reference_one)
    call integrator_state%compute_reference_wgt_from_x(&
         1,[map_two_x],reference_two)
    density_one=1d0/reference_one
    density_two=2d0*physical_x/reference_two
    alpha_one=density_one/(density_one+density_two)
    alpha_two=density_two/(density_one+density_two)
    max_error=max(max_error,abs(alpha_one+alpha_two-1d0))
    max_error=max(max_error,abs(&
         alpha_one*(1d0+4d0*physical_x)+&
         alpha_two*(1d0+4d0*physical_x)-&
         (1d0+4d0*physical_x)))
    if (map_id.eq.1) then
       partition_factor=2d0*alpha_one
    else
       partition_factor=2d0*alpha_two
    endif
    value=live_weight*selected_jacobian*partition_factor*&
         (1d0+4d0*physical_x)
    value_abs=value
  end subroutine evaluate_family_point

end program simple_integrator_multichannel_regression

real(kind=8) function ran2()
  implicit none
  integer(kind=8),save :: state=104729_8
  state=mod(16807_8*state,2147483647_8)
  ran2=dble(state)/2147483647d0
end function ran2
