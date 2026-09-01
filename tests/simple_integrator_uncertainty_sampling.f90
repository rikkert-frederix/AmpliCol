program simple_integrator_uncertainty_sampling_test
  use simple_integrator_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_finite
  implicit none
  integer,parameter :: nchannel=2
  integer,parameter :: accuracy_pilot_points=4096
  integer,parameter :: event_first_iteration_points=2048
  integer :: accuracy_counts(nchannel,2),event_counts(nchannel,2)
  integer :: zero_counts_a(nchannel,2),zero_counts_b(nchannel,2),signed_counts(nchannel,2)
  integer :: growth_counts(nchannel,2)

  call sample_second_iteration(.true.,accuracy_counts)
  call sample_second_iteration(.false.,event_counts)
  call sample_accuracy_pilot(1,.true.,zero_counts_a)
  call sample_accuracy_pilot(10000,.true.,zero_counts_b)
  call sample_accuracy_pilot(1,.false.,signed_counts)
  call sample_quota_growth(growth_counts)
  call check_external_convergence_gate()
  call check_independent_adaptation_classes()
  call check_invalid_point_rejection()
  call check_zero_rate_event_termination()
  call check_single_event_multichannel_completion()
  call check_integrator_instance_isolation()

  call assert_true(sum(accuracy_counts(1,:)).gt.3*sum(accuracy_counts(2,:)),&
       'accuracy mode favours the channel with the largest uncertainty')
  call assert_true(accuracy_counts(1,1).gt.3*accuracy_counts(1,2),&
       'accuracy mode favours the integral with the largest uncertainty')
  call assert_true(all(accuracy_counts.gt.0),&
       'accuracy mode keeps every channel and integral active')
  call assert_true(minval(accuracy_counts).ge.512,&
       'accuracy mode reserves one quarter of the second iteration for exploration')
  call assert_true(maxval(event_counts).lt.1.2d0*minval(event_counts),&
       'event mode keeps rate-based balanced allocation')
  call assert_true(all(zero_counts_a.eq.1024) .and. all(zero_counts_b.eq.1024),&
       'accuracy pilot counts all evaluated zero points')
  call assert_true(all(zero_counts_a.eq.zero_counts_b),&
       'accuracy allocation is independent of --nevents')
  call assert_true(all(signed_counts.eq.1024),&
       'signed cancellation does not affect absolute-envelope pilot allocation')
  call assert_true(minval(growth_counts).ge.1024,&
       'accuracy exploration remains active after the pilot iteration')
  call assert_true(growth_counts(2,2).le.8192,&
       'a newly noisy leaf cannot grow by more than sixteen times in one iteration')
  call assert_true(growth_counts(2,2).gt.max(growth_counts(1,1),growth_counts(2,1),growth_counts(1,2)),&
       'the growth cap still lets the newly noisy leaf receive the largest quota')

  write(*,'(a)') 'simple integrator uncertainty sampling test: PASS'

contains

  subroutine sample_second_iteration(integration_only,counts)
    logical,intent(in) :: integration_only
    integer,intent(out) :: counts(nchannel,2)
    type(integrator) :: integ
    integer :: ichan,iint,nfirst,nfirst_target
    integer :: ndim(nchannel),ndim_extra(nchannel),nintegral(nchannel)
    real(kind=8) :: f(1),f_abs(1)
    logical :: to_write(1),done
    character(len=512) :: log_line
    integer :: ios

    ndim=1
    ndim_extra=0
    nintegral=2
    open(unit=99,status='scratch',action='readwrite')
    if (integration_only) then
       ! A request for only one event would normally mark most strata as
       ! event-generation DONE.  Accuracy mode must ignore that bookkeeping.
       call integ%init(nchannel,ndim,ndim_extra,nintegral,1,2,&
            accuracy=1d-12)
    else
       call integ%init(nchannel,ndim,ndim_extra,nintegral,100,2)
    endif

    if (integration_only) then
       nfirst_target=accuracy_pilot_points
    else
       nfirst_target=event_first_iteration_points
    endif

    done=.false.
    do nfirst=1,nfirst_target
       call integ%get_points(1,ichan,iint)
       call test_integrand(integ%x(1,1),ichan,iint,f_abs(1))
       f=f_abs
       call integ%fill_points(1,f_abs,f,to_write,done)
       if (done) then
          write(*,*) 'FAIL: integrator stopped before the second iteration'
          stop 1
       endif
    enddo

    counts=0
    do while (.not.done)
       call integ%get_points(1,ichan,iint)
       counts(ichan,iint)=counts(ichan,iint)+1
       call test_integrand(integ%x(1,1),ichan,iint,f_abs(1))
       f=f_abs
       call integ%fill_points(1,f_abs,f,to_write,done)
    enddo

    if (integration_only) then
       flush(99)
       rewind(99)
       do
          read(99,'(a)',iostat=ios) log_line
          if (ios.ne.0) exit
          if (index(log_line,'DONE').ne.0) then
             write(*,*) 'FAIL: accuracy log contains a DONE label'
             stop 1
          endif
       enddo
    endif
    close(99)
  end subroutine sample_second_iteration

  subroutine sample_accuracy_pilot(nevts,all_zero,counts)
    integer,intent(in) :: nevts
    logical,intent(in) :: all_zero
    integer,intent(out) :: counts(nchannel,2)
    type(integrator) :: integ
    integer :: ichan,iint,npoint
    integer :: ndim(nchannel),ndim_extra(nchannel),nintegral(nchannel)
    real(kind=8) :: f(1),f_abs(1)
    logical :: to_write(1),done,accepted(1)

    ndim=1
    ndim_extra=0
    nintegral=2
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(nchannel,ndim,ndim_extra,nintegral,nevts,1,accuracy=1d-12)
    counts=0
    done=.false.
    do npoint=1,accuracy_pilot_points
       call integ%get_points(1,ichan,iint)
       counts(ichan,iint)=counts(ichan,iint)+1
       if (all_zero) then
          f_abs=0d0
          f=0d0
          accepted=.false.
       else
          f_abs=0.5d0+integ%x(1,1)
          if (ichan.eq.1) then
             f=f_abs
          else
             f=-f_abs
          endif
          accepted=.true.
       endif
       call integ%fill_points(1,f_abs,f,to_write,done,accepted)
    enddo
    if (.not.done) then
       write(*,*) 'FAIL: accuracy pilot did not finish at its exact quota'
       stop 1
    endif
    close(99)
  end subroutine sample_accuracy_pilot

  subroutine sample_quota_growth(counts)
    integer,intent(out) :: counts(nchannel,2)
    type(integrator) :: integ
    integer :: ichan,iint,stage
    integer :: ndim(nchannel),ndim_extra(nchannel),nintegral(nchannel)
    real(kind=8) :: f(1),f_abs(1)
    logical :: to_write(1),done,iteration_finished

    ndim=1
    ndim_extra=0
    nintegral=2
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(nchannel,ndim,ndim_extra,nintegral,1,3,accuracy=1d-12)
    counts=0
    stage=1
    done=.false.
    do while (.not.done)
       call integ%get_points(1,ichan,iint)
       if ((stage.eq.1 .and. ichan.eq.1 .and. iint.eq.1) .or. &
            (stage.eq.2 .and. ichan.eq.2 .and. iint.eq.2)) then
          f_abs(1)=0.5d0+integ%x(1,1)
       else
          f_abs(1)=1d0
       endif
       f=f_abs
       if (stage.eq.3) counts(ichan,iint)=counts(ichan,iint)+1
       call integ%fill_points(1,f_abs,f,to_write,done,iteration_finished=iteration_finished)
       if (iteration_finished) stage=stage+1
    enddo
    close(99)
  end subroutine sample_quota_growth

  subroutine check_external_convergence_gate()
    type(integrator) :: integ
    integer :: ichan,iint,npoint
    integer :: ndim(1),ndim_extra(1),nintegral(1)
    real(kind=8) :: f(1),f_abs(1)
    logical :: to_write(1),done,iteration_finished

    ndim=1
    ndim_extra=0
    nintegral=1
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(1,ndim,ndim_extra,nintegral,1,2,accuracy=0.9d0)
    done=.false.
    do npoint=1,1024
       call integ%get_points(1,ichan,iint)
       f=1d0
       f_abs=1d0
       call integ%fill_points(1,f_abs,f,to_write,done,iteration_finished=iteration_finished,&
            external_converged=.false.)
    enddo
    call assert_true(iteration_finished,'external convergence test did not finish its pilot iteration')
    call assert_true(.not.done,'external convergence gate did not defer an otherwise converged integration')
    do while (.not.done)
       call integ%get_points(1,ichan,iint)
       f=1d0
       f_abs=1d0
       call integ%fill_points(1,f_abs,f,to_write,done,external_converged=.true.)
    enddo
    close(99)
  end subroutine check_external_convergence_gate

  subroutine check_independent_adaptation_classes()
    type(integrator) :: integ
    integer :: ichan,iint
    integer :: ndim(1),ndim_extra(1),nintegral(1),classes(2,1),counts(2)
    real(kind=8) :: f(1),f_abs(1),xsum(2),weight_class1,weight_class2,weight_integral1
    logical :: to_write(1),done,iteration_finished,second_iteration

    ndim=1
    ndim_extra=0
    nintegral=2
    classes(:,1)=[1,2]
    counts=0
    xsum=0d0
    second_iteration=.false.
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(1,ndim,ndim_extra,nintegral,1,2,accuracy=0.9d0,adaptation_classes=classes)
    done=.false.
    do while (.not.done)
       call integ%get_points(1,ichan,iint)
       if (second_iteration) then
          counts(iint)=counts(iint)+1
          xsum(iint)=xsum(iint)+integ%x(1,1)
       endif
       if (iint.eq.1) then
          f=exp(12d0*integ%x(1,1))*integ%wgt(1)
       else
          f=exp(12d0*(1d0-integ%x(1,1)))*integ%wgt(1)
       endif
       f_abs=abs(f)
       call integ%fill_points(1,f_abs,f,to_write,done,iteration_finished=iteration_finished,&
            external_converged=second_iteration)
       if (iteration_finished .and. .not.second_iteration) second_iteration=.true.
    enddo
    close(99)
    call assert_true(all(counts.gt.0),'adaptation-class test did not sample both second-iteration leaves')
    call assert_true(xsum(1)/counts(1).gt.0.6d0,'upper-tail adaptation class did not move its grid')
    call assert_true(xsum(2)/counts(2).lt.0.4d0,'lower-tail adaptation class did not move its grid')
    call integ%compute_wgt_from_x(1,[0.25d0],weight_class1,adaptation_class=1)
    call integ%compute_wgt_from_x(1,[0.25d0],weight_class2,adaptation_class=2)
    call integ%compute_wgt_from_x(1,[0.25d0],weight_integral1,iint=1)
    call assert_true(weight_class1.gt.0d0 .and. weight_class2.gt.0d0,&
         'adaptation-class grid lookup returned a non-positive density weight')
    call assert_true(abs(weight_class1-weight_class2).gt.&
         1d-3*max(weight_class1,weight_class2),&
         'adaptation-class grid lookup ignored the requested class')
    call assert_true(abs(weight_class1-weight_integral1).le.&
         1d-13*max(1d0,weight_class1),&
         'integral and explicit adaptation-class grid lookups disagree')
  end subroutine check_independent_adaptation_classes

  subroutine check_invalid_point_rejection()
    type(integrator) :: integ
    integer :: ichan,iint,npoint
    integer :: ndim(1),ndim_extra(1),nintegral(1)
    real(kind=8) :: f(1),f_abs(1),aux(1,1)
    real(kind=8),allocatable :: result(:,:),uncertainty(:,:),aux_result(:,:),aux_uncertainty(:,:)
    real(kind=8) :: expected,recomputed_weight
    logical :: to_write(1),done

    ndim=1
    ndim_extra=0
    nintegral=1
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(1,ndim,ndim_extra,nintegral,1,1,accuracy=0.9d0,naux=1)
    call integ%compute_wgt_from_x(1,[ieee_value(0d0,ieee_quiet_nan)],recomputed_weight)
    call assert_true(recomputed_weight.eq.0d0,'invalid coordinate rejected by grid reweighting')
    done=.false.
    do npoint=1,1024
       call integ%get_points(1,ichan,iint)
       f=2d0
       f_abs=2d0
       aux=3d0
       select case (npoint)
       case (1)
          f=ieee_value(0d0,ieee_quiet_nan)
       case (2)
          f_abs=ieee_value(0d0,ieee_positive_inf)
       case (3)
          aux=ieee_value(0d0,ieee_quiet_nan)
       case (4)
          integ%x(1,1)=ieee_value(0d0,ieee_quiet_nan)
       case (5)
          integ%wgt(1)=ieee_value(0d0,ieee_positive_inf)
       case (6)
          ! Finite is not enough: the integrator also squares every value.
          f=sqrt(huge(1d0))
       case (7)
          integ%wgt(1)=0d0
       case (8)
          integ%wgt(1)=-1d0
       case (9)
          f=0.5d0*sqrt(tiny(1d0))
       case (10)
          f=3d0
          f_abs=2d0
       end select
       call integ%fill_points(1,f_abs,f,to_write,done,f_aux=aux)
    enddo
    call assert_true(done,'invalid-point guard did not complete its exact accuracy quota')
    call integ%get_channel_results(result,uncertainty)
    call integ%get_channel_aux_results(aux_result,aux_uncertainty)
    expected=2d0*1014d0/1024d0
    call assert_true(all(ieee_is_finite(result)) .and. all(ieee_is_finite(uncertainty)),&
         'invalid point contaminated the primary integrator estimate')
    call assert_true(abs(result(1,1)-expected).lt.1d-12 .and. abs(result(2,1)-expected).lt.1d-12,&
         'invalid points were not replaced by rejected zero contributions')
    expected=3d0*1014d0/1024d0
    call assert_true(all(ieee_is_finite(aux_result)) .and. all(ieee_is_finite(aux_uncertainty)),&
         'invalid point contaminated an auxiliary integrator estimate')
    call assert_true(abs(aux_result(1,1)-expected).lt.1d-12,&
         'invalid auxiliary point was not replaced by a rejected zero contribution')
    close(99)
  end subroutine check_invalid_point_rejection

  subroutine check_zero_rate_event_termination()
    type(integrator) :: integ
    integer :: ichan,iint,npoint
    integer :: ndim(1),ndim_extra(1),nintegral(1)
    real(kind=8) :: f(1),f_abs(1)
    real(kind=8),allocatable :: result(:,:),uncertainty(:,:)
    logical :: to_write(1),done

    ndim=1
    ndim_extra=0
    nintegral=1
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(1,ndim,ndim_extra,nintegral,1,1)
    f=0d0
    f_abs=0d0
    done=.false.
    do npoint=1,4096
       call integ%get_points(1,ichan,iint)
       call integ%fill_points(1,f_abs,f,to_write,done)
       if (done) exit
    enddo
    call assert_true(done,'zero-rate event integration did not terminate')
    call assert_true(npoint.lt.4096,'zero-rate event integration exhausted its safety bound')
    call integ%get_channel_results(result,uncertainty)
    call assert_true(all(ieee_is_finite(result)) .and. all(result.eq.0d0),&
         'zero-rate event result is not finite zero')
    call assert_true(all(ieee_is_finite(uncertainty)) .and. all(uncertainty.eq.0d0),&
         'zero-rate event uncertainty is not finite zero')
    close(99)
  end subroutine check_zero_rate_event_termination

  subroutine check_single_event_multichannel_completion()
    integer,parameter :: nch=4
    type(integrator) :: integ
    integer :: ichan,iint,npoint
    integer :: ndim(nch),ndim_extra(nch),nintegral(nch)
    real(kind=8) :: f(1),f_abs(1)
    real(kind=8),allocatable :: weights(:,:)
    logical :: to_write(1),done

    ndim=1
    ndim_extra=0
    nintegral=1
    open(unit=99,status='scratch',action='readwrite')
    call integ%init(nch,ndim,ndim_extra,nintegral,1,8)
    f=1d0
    f_abs=1d0
    done=.false.
    do npoint=1,100000
       call integ%get_points(1,ichan,iint)
       call integ%fill_points(1,f_abs,f,to_write,done)
       if (done) exit
    enddo
    call assert_true(done,'single-event multichannel integration did not terminate')
    call integ%assign_evnt_wgts(weights)
    call assert_true(count(weights(1,:).ne.0d0).eq.1,&
         'single-event quota lost its selected candidate across iterations')
    call assert_true(all(ieee_is_finite(weights)),&
         'single-event multichannel weights are non-finite')
    close(99)
  end subroutine check_single_event_multichannel_completion

  subroutine check_integrator_instance_isolation()
    type(integrator) :: event_integ,accuracy_integ
    integer :: ichan,iint,npoint
    integer :: ndim(1),ndim_extra(1),nintegral(1)
    real(kind=8) :: f(1),f_abs(1)
    real(kind=8),allocatable :: weights(:,:)
    logical :: to_write(1),done

    ndim=1
    ndim_extra=0
    nintegral=1
    open(unit=99,status='scratch',action='readwrite')
    ! Initialising the accuracy driver after the event driver used to mutate
    ! module SAVE variables still consulted by the first live instance.
    call event_integ%init(1,ndim,ndim_extra,nintegral,1,8)
    call accuracy_integ%init(1,ndim,ndim_extra,nintegral,10000,1,accuracy=0.9d0)

    f=1d0
    f_abs=1d0
    done=.false.
    do npoint=1,100000
       call event_integ%get_points(1,ichan,iint)
       call event_integ%fill_points(1,f_abs,f,to_write,done)
       if (done) exit
    enddo
    call assert_true(done,'a second integrator corrupted live event-generation state')
    call event_integ%assign_evnt_wgts(weights)
    call assert_true(count(weights(1,:).ne.0d0).eq.1,&
         'isolated event integrator did not retain its requested sample')

    done=.false.
    do npoint=1,1024
       call accuracy_integ%get_points(1,ichan,iint)
       call accuracy_integ%fill_points(1,f_abs,f,to_write,done)
       call assert_true(.not.any(to_write),&
            'isolated accuracy integrator attempted to write an event')
    enddo
    call assert_true(done,'live event integration corrupted accuracy-only state')
    close(99)
  end subroutine check_integrator_instance_isolation

  subroutine test_integrand(x,ichan,iint,value)
    real(kind=8),intent(in) :: x
    integer,intent(in) :: ichan,iint
    real(kind=8),intent(out) :: value

    if (ichan.eq.1 .and. iint.eq.1) then
       ! Same mean as the other strata, but a much larger uncertainty.
       value=0.5d0+x
    else
       value=1d0
    endif
  end subroutine test_integrand

  subroutine assert_true(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if (.not.condition) then
       write(*,*) 'FAIL: ',trim(label)
       write(*,*) ' accuracy counts:',accuracy_counts
       write(*,*) ' event counts:',event_counts
       stop 1
    endif
  end subroutine assert_true

end program simple_integrator_uncertainty_sampling_test
