program simple_integrator_uncertainty_sampling_test
  use simple_integrator_mod
  implicit none
  integer,parameter :: nchannel=2
  integer,parameter :: accuracy_pilot_points=4096
  integer,parameter :: event_first_iteration_points=2048
  integer :: accuracy_counts(nchannel,2),event_counts(nchannel,2)
  integer :: zero_counts_a(nchannel,2),zero_counts_b(nchannel,2),signed_counts(nchannel,2)

  call sample_second_iteration(.true.,accuracy_counts)
  call sample_second_iteration(.false.,event_counts)
  call sample_accuracy_pilot(1,.true.,zero_counts_a)
  call sample_accuracy_pilot(10000,.true.,zero_counts_b)
  call sample_accuracy_pilot(1,.false.,signed_counts)

  call assert_true(sum(accuracy_counts(1,:)).gt.3*sum(accuracy_counts(2,:)),&
       'accuracy mode favours the channel with the largest uncertainty')
  call assert_true(accuracy_counts(1,1).gt.3*accuracy_counts(1,2),&
       'accuracy mode favours the integral with the largest uncertainty')
  call assert_true(all(accuracy_counts.gt.0),&
       'accuracy mode keeps every channel and integral active')
  call assert_true(maxval(event_counts).lt.1.2d0*minval(event_counts),&
       'event mode keeps rate-based balanced allocation')
  call assert_true(all(zero_counts_a.eq.1024) .and. all(zero_counts_b.eq.1024),&
       'accuracy pilot counts all evaluated zero points')
  call assert_true(all(zero_counts_a.eq.zero_counts_b),&
       'accuracy allocation is independent of --nevents')
  call assert_true(all(signed_counts.eq.1024),&
       'signed cancellation does not affect absolute-envelope pilot allocation')

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
