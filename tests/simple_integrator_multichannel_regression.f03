program simple_integrator_multichannel_regression
  use simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: batch_size=64,max_observed_iters=12
  type(integrator) :: test_integrator,one_iter_integrator
  integer,dimension(2) :: ndim,ndim_extra,nintegral
  integer,dimension(1) :: one_ndim,one_extra,one_integral
  real(kind=8),dimension(batch_size) :: f,f_abs
  logical,dimension(batch_size) :: to_write
  logical,dimension(max_observed_iters,2) :: selected_seen
  real(kind=8),allocatable,dimension(:,:) :: event_wgts
  real(kind=8) :: live_probe,reference_probe,last_live_probe
  real(kind=8) :: reference_snapshot,live_at_snapshot,tolerance
  integer :: ichan,iint,ip,observed_iter,accepted,i
  logical :: done,reference_separated,live_changed_after_snapshot
  logical :: early_channel_completion,saw_warmup_log,saw_production_log
  character(len=256) :: log_line
  integer :: io_status

  open(unit=99,status='scratch',action='readwrite')
  ndim=1
  ndim_extra=0
  nintegral=1
  call test_integrator%init(2,ndim,ndim_extra,nintegral,100,7)
  selected_seen=.false.
  observed_iter=1
  reference_separated=.false.
  live_changed_after_snapshot=.false.
  reference_snapshot=0d0
  live_at_snapshot=0d0
  call test_integrator%compute_wgt_from_x(2,[0.317d0],last_live_probe)

  done=.false.
  do while (.not.done)
     call test_integrator%get_points(batch_size,ichan,iint)
     if (observed_iter.le.max_observed_iters) &
          selected_seen(observed_iter,ichan)=.true.
     do ip=1,batch_size
        call evaluate_multichannel_point(test_integrator,ichan,&
             test_integrator%x(1,ip),test_integrator%wgt(ip),f(ip),f_abs(ip))
     enddo
     call test_integrator%fill_points(batch_size,f_abs,f,to_write,done)

     call test_integrator%compute_wgt_from_x(2,[0.317d0],live_probe)
     call test_integrator%compute_reference_wgt_from_x(2,[0.317d0],reference_probe)
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
     if (abs(live_probe-last_live_probe).gt.tolerance) then
        observed_iter=observed_iter+1
        last_live_probe=live_probe
     endif
  enddo

  if (.not.reference_separated) then
     write (*,*) 'Reference and live grids never separated in production'
     stop 1
  endif
  if (.not.live_changed_after_snapshot) then
     write (*,*) 'Live grids did not keep adapting after the snapshot'
     stop 1
  endif
  if (.not.ieee_is_finite(test_integrator%res(2)) .or. &
       .not.ieee_is_finite(test_integrator%unc(2))) then
     write (*,*) 'Non-finite accumulated production rate',&
          test_integrator%res(2),test_integrator%unc(2)
     stop 1
  endif
  if (abs(test_integrator%res(2)-3d0).gt.&
       3d0*test_integrator%unc(2)+2d-3) then
     write (*,*) 'Stationary multichannel rate disagrees with analytic result',&
          test_integrator%res(2),test_integrator%unc(2)
     stop 1
  endif

  early_channel_completion=.false.
  do i=2,max_observed_iters
     if (selected_seen(i,2) .and. .not.selected_seen(i,1) .and. &
          any(selected_seen(1:i-1,1))) early_channel_completion=.true.
  enddo
  if (.not.early_channel_completion) then
     write (*,*) 'No production iteration observed one channel completing early'
     stop 1
  endif

  call test_integrator%assign_evnt_wgts(event_wgts)
  if (any(.not.ieee_is_finite(event_wgts))) then
     write (*,*) 'Stored-event grid re-evaluation produced non-finite weights'
     stop 1
  endif
  accepted=count(event_wgts(1,:).gt.0d0)
  if (accepted.ne.100) then
     write (*,*) 'Unexpected number of accepted stored events',accepted
     stop 1
  endif
  do i=1,size(event_wgts,2)
     if (event_wgts(1,i).le.0d0) cycle
     if (abs(event_wgts(1,i)-test_integrator%res(1)).gt.&
          1d-12*max(1d0,abs(test_integrator%res(1)))) then
        write (*,*) 'Accepted stored events have inconsistent nominal weights'
        stop 1
     endif
  enddo
  deallocate(event_wgts)

  ! A one-iteration configuration has no warm-up: its uniform initial grid is
  ! the reference snapshot and the sole iteration contributes to production.
  one_ndim=1
  one_extra=0
  one_integral=1
  call one_iter_integrator%init(1,one_ndim,one_extra,one_integral,8,1)
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
  call one_iter_integrator%compute_reference_wgt_from_x(1,[0.417d0],reference_probe)
  if (abs(live_probe-reference_probe).gt.2d-14 .or. &
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
  write (*,*) 'simple-integrator multichannel regression: PASS'

contains

  subroutine evaluate_multichannel_point(integrator_state,channel_number,x,&
       live_weight,value,value_abs)
    type(integrator),intent(inout) :: integrator_state
    integer,intent(in) :: channel_number
    real(kind=8),intent(in) :: x,live_weight
    real(kind=8),intent(out) :: value,value_abs
    real(kind=8) :: physical_x,partner_x,reference_one,reference_two
    real(kind=8) :: density_one,density_two,partition,envelope

    if (channel_number.eq.1) then
       physical_x=x
    else
       physical_x=1d0-x
    endif
    partner_x=1d0-physical_x
    call integrator_state%compute_reference_wgt_from_x(1,[physical_x],reference_one)
    call integrator_state%compute_reference_wgt_from_x(2,[partner_x],reference_two)
    density_one=1d0/reference_one
    density_two=1d0/reference_two
    if (channel_number.eq.1) then
       partition=density_one/(density_one+density_two)
    else
       partition=density_two/(density_one+density_two)
    endif
    value=live_weight*partition*(1d0+4d0*physical_x)
    if (channel_number.eq.1) then
       envelope=1d0
    elseif (modulo(97d0*physical_x,1d0).lt.0.008d0) then
       envelope=125d0
    else
       envelope=1d-9
    endif
    value_abs=value*envelope
  end subroutine evaluate_multichannel_point

end program simple_integrator_multichannel_regression

real(kind=8) function ran2()
  implicit none
  integer(kind=8),save :: state=104729_8
  state=mod(16807_8*state,2147483647_8)
  ran2=dble(state)/2147483647d0
end function ran2
