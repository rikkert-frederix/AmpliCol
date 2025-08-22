! 1. All variables should be allocatables
! 2. Keep track of channels and integrals per channel
! 3. Update grids based on maximum weight found
! 4. First iteration should compute all channels and all integrals
! 5. After iteration three (or based on uncertainty in the integral?), start event generation based on maximum weight found so far. Update number of points accordingly...
! 6. After each iteration, re-unweight the events, and check if we have enough.
! 7. Use a single maximum weight (not a grid) -- one per channel and integral. 
! 8. For which channel and integral to throw the next point????
! 9. Allow for multiple points "in parallel"
! 10. Recycle events between iterations????
! 11. Keep track of overweights (optionally pass them to event file)


module simple_integrator_mod
  implicit none
  private
  type :: channel
     integer :: ndim,nintegral,current_integral,current_iteration&
          &,number,max_iterations
     integer(kind=8) :: npoints,npoints_iter
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,unc_iter
     logical :: done,evgen_done
     type(grid),allocatable,dimension(:,:) :: grids
     type(integral),allocatable,dimension(:) :: integrals
   contains
     procedure,private :: init => channel_init
     procedure,private :: add_point => channel_add_point
     procedure,private :: get_point => channel_get_point
     procedure,private :: finalise_iteration => channel_finalise_iteration
     procedure,private :: update_result_iteration => channel_update_result_iteration
     procedure,private :: combine_iterations => channel_combine_iterations
     procedure,private :: print_result_iteration => channel_print_result_iteration
     procedure,private :: print_combined_result => channel_print_combined_result
     procedure,private :: init_next_iteration => channel_init_next_iteration
     procedure,private :: check_generated_events => channel_check_generated_events
     procedure,private :: recompute_wgt_from_x
  end type channel
  type :: integral
     real(kind=8) :: max_value
     real(kind=8),dimension(:),allocatable :: f_max
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,accum&
          &,accum2,unc_iter
     integer :: ichan,n_unwgt,event,nevent_in_list,ndim &
          &,current_iteration,max_iterations,nevents_unw
     integer(kind=8) :: npoints_iter,npoints,npoints_requested&
          &,npoints_nonzero,npoints_nonzero_total
     logical :: done,evgen_done
     type(event),dimension(:),allocatable :: event_list
   contains
     procedure,private :: init => integral_init
     procedure,private :: add_point => integral_add_point
     procedure,private :: update_result_iteration => integral_update_result_iteration
     procedure,private :: combine_iterations => integral_combine_iterations
     procedure,private :: update_max_value,check_write_event,increase_size_event_list
     procedure,private :: init_next_iteration => integral_init_next_iteration
     procedure,private :: compute_fmax => integral_compute_fmax
     procedure,private :: compute_fmax_next_iteration => integral_compute_fmax_next_iteration
     procedure,private :: unwgt => integral_unwgt
  end type integral
  type :: grid
     integer :: size
     real(kind=8),allocatable,dimension(:) :: current,accum
     integer,allocatable,dimension(:) :: nhits
   contains
     procedure,private :: init => grid_init
     procedure,private :: add_point => grid_add_point
     procedure,private :: get_x,get_wgt,massage_accum
     procedure,private :: update => grid_update
  end type grid
  type :: event
     real(kind=8),allocatable,dimension(:) :: x,f_abs
     real(kind=8) :: wgt,rnd
     integer :: iter
     logical :: unwgt
  end type event
  type,public :: integrator
     integer :: nchannel,current_channel,nevents_unw,npoints_generated
     integer(kind=8) :: npoints_requested
     real(kind=8),dimension(2) :: res,unc
     real(kind=8),allocatable,dimension(:,:),public :: x
     real(kind=8),allocatable,dimension(:),public :: wgt
     integer,allocatable,dimension(:,:) :: cell
     type(channel),allocatable,dimension(:) :: channels
     ! keep track of maximum weights of all events written---integral-by-integral---and approx how many more we should generate in that channel.
   contains
     procedure,public :: init,get_points,fill_points,compute_wgt_from_x
     procedure,private :: read_all_grids,write_all_grids&
          &,get_channel_and_integral,update_points_requested&
          &,print_results,compute_total_rate,check_if_done&
          &,init_next_iteration,get_npoints_nonzero_iter
     ! fill_points should return 'done' when ready; also it should
     ! keep track of number of points thrown and determine if it needs
     ! re-gridding; furthermore, it should tell the main codes which
     ! events should be kept on disk?  get_points should give a set of
     ! points (determined by an argument)--all for the same channel
     ! (and integral); fill_points should return at least a subset of
     ! those --- it may be assumed that these are in the beginning of
     ! the list. It can be assumed that the points that were not used
     ! were 'lost' and should be treated as not-even-generated.
  end type integrator
  double precision, external :: ran2
  integer,parameter :: importance_sampling_strategy=3
  real(kind=8),parameter :: write_event_fraction=1d0
  integer :: iterations_without_events
  integer,parameter :: min_points_per_channel=1024
  integer,parameter :: min_points_per_integral=128
contains

  subroutine init(this,nchannel,ndim,nintegral,nevents_unw,niters)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: nchannel,nevents_unw,niters
    integer,dimension(nchannel),intent(in) :: ndim,nintegral
    integer :: i
    this%nchannel=nchannel
    this%nevents_unw=nevents_unw
    ! if we assume 0.3% unweighting efficiency, we expect ~10% time
    ! spend in iterations that do not produce events:
    iterations_without_events=5
    this%npoints_requested=int(nevents_unw/(0.03*2**iterations_without_events),kind=8)
    do while (this%npoints_requested/this%nchannel.lt.min_points_per_channel &
         .and. iterations_without_events.gt.3)
       iterations_without_events=iterations_without_events-1
       this%npoints_requested=nevents_unw/(0.03*2**iterations_without_events)
    enddo
    this%npoints_requested=max(this%npoints_requested,min_points_per_channel*this%nchannel)
    allocate(this%channels(this%nchannel))
    do i=1,this%nchannel
       call this%channels(i)%init(ndim(i),nintegral(i),this%npoints_requested/nchannel,niters,i)
    enddo
    this%current_channel=0
    this%npoints_generated=0
    this%res=0d0
  end subroutine init

  subroutine channel_init(this,ndim,nintegral,npoints,niters,ichan)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: ndim,nintegral,niters,ichan
    integer(kind=8) :: npoints
    integer :: i
    this%ndim=ndim
    this%max_iterations=niters
    this%nintegral=nintegral
    this%number=ichan
    allocate(this%grids(1:this%ndim,1:this%max_iterations))
    allocate(this%integrals(1:this%nintegral))
    do i=1,this%ndim
       call this%grids(i,1)%init()
    enddo
    do i=1,this%nintegral
       call this%integrals(i)%init(ndim,npoints/this%nintegral,this%number,this%max_iterations)
    enddo
    this%current_integral=0
    this%current_iteration=0
    this%npoints=0_8
    this%res=0d0
    this%unc=0d0
    call this%init_next_iteration()
    this%evgen_done=.false.
  end subroutine channel_init

  subroutine channel_init_next_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    this%res_iter=0d0
    this%res2_iter=0d0
    this%unc_iter=0d0
    this%npoints_iter=0_8
    this%done=.false.
    if (all(this%integrals%evgen_done)) then
       this%evgen_done=.true.
       this%done=.true.
    else
       this%evgen_done=.false.
    endif
    do i=1,this%nintegral
       call this%integrals(i)%init_next_iteration()
    enddo
    this%current_iteration=this%current_iteration+1
  end subroutine channel_init_next_iteration
  
  subroutine integral_init_next_iteration(this)
    implicit none
    class(integral),intent(inout) :: this
    this%current_iteration=this%current_iteration+1
    this%res_iter=0d0
    this%res2_iter=0d0
    this%accum=0d0
    this%accum2=0d0
    this%unc_iter=0d0
    this%npoints_iter=0_8
    this%npoints_nonzero=0_8
    this%event=0
    if (.not. this%evgen_done) this%done=.false.
    if (this%current_iteration.eq.1) then
       this%f_max(this%current_iteration)=-1d0
    elseif (this%current_iteration.le.iterations_without_events+1) then
       this%f_max(this%current_iteration)=this%max_value
    endif
    this%max_value=0
  end subroutine integral_init_next_iteration
  
  subroutine grid_init(this,current)
    implicit none
    integer,parameter :: grid_size=32
    class(grid),intent(inout) :: this
    real(kind=8),dimension(0:grid_size),optional :: current
    integer :: i
    this%size=grid_size
    allocate(this%current(0:this%size))
    if (present(current)) then
       this%current=current
    else
       do i=0,this%size
          this%current(i)=dble(i)/this%size
       enddo
    endif
    allocate(this%accum(0:this%size))
    allocate(this%nhits(this%size))
    this%accum=0d0
    this%nhits=0
  end subroutine grid_init
  
  subroutine integral_init(this,ndim,npoints,ichan,niters)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: ndim,ichan,niters
    integer(kind=8) :: npoints
    this%ndim=ndim
    this%ichan=ichan
    this%npoints=0_8
    this%npoints_requested=npoints
    this%max_iterations=niters
    allocate(this%f_max(this%max_iterations))
    this%f_max=-1d0
    allocate(this%event_list(npoints))
    this%nevent_in_list=0
    this%evgen_done=.false.
    this%current_iteration=0
    this%npoints_nonzero_total=0_8
  end subroutine integral_init
  
  subroutine get_points(this,npoints,ichan,iint)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: npoints
    integer,intent(out) :: ichan,iint
    integer :: i
    real(kind=8) :: wgt_chan
    
    call this%get_channel_and_integral(ichan,iint,wgt_chan)
    this%current_channel=ichan
    
    allocate(this%x(1:this%channels(this%current_channel)%ndim,1:npoints))
    allocate(this%cell(1:this%channels(this%current_channel)%ndim,1:npoints))
    allocate(this%wgt(1:npoints))

    do i=1,npoints
       call this%channels(this%current_channel)%get_point(this%x(1,i),this%cell(1,i),this%wgt(i))
    enddo
    this%wgt=this%wgt*wgt_chan
    this%npoints_generated=npoints

  end subroutine get_points
  
  subroutine fill_points(this,npoints,f_abs,f,to_write,done)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: npoints
    real(kind=8),dimension(npoints),intent(in) :: f,f_abs
    logical,dimension(npoints),intent(out) :: to_write
    logical,intent(out) :: done
    integer :: i
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    character(len=19) :: formatted
    integer(kind=8) :: npoints_nonzero
    done=.false.
    if (npoints.gt.this%npoints_generated) then
       write (*,*) 'ERROR: too many points returned'
       stop 1
    endif
    do i=1,npoints
       call this%channels(this%current_channel)%add_point(this%x(1,i),this%wgt(i),this%cell(1,i),f_abs(i),f(i),to_write(i))
    enddo
    this%npoints_generated=0
    if (all(this%channels%done)) then
       call date_and_time(date, time, zone)
       write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
            date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
       call this%get_npoints_nonzero_iter(npoints_nonzero)
       write (*,'(a,x,i4,x,a,x,i10,x,a)') &
            'iteration',this%channels(1)%current_iteration,'(',npoints_nonzero, &
            'points) '//trim(formatted)//' :'
       do i=1,this%nchannel
          call this%channels(i)%finalise_iteration()
       enddo
       call this%compute_total_rate()
       call this%check_if_done()
       call this%print_results()
       write (*,*) ''
       call this%init_next_iteration()
       if (all(this%channels%evgen_done)) done=.true.
    endif
    if (all(this%channels%done)) done=.true.
    deallocate(this%x)
    deallocate(this%cell)
    deallocate(this%wgt)
  end subroutine fill_points
  
  subroutine get_npoints_nonzero_iter(this,npoints_nonzero)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i,j
    integer(kind=8) :: npoints_nonzero
    npoints_nonzero=0_8
    do i=1,this%nchannel
       npoints_nonzero=npoints_nonzero+sum(this%channels(i)%integrals(1:this%channels(i)%nintegral)%npoints_nonzero)
    enddo
  end subroutine get_npoints_nonzero_iter
  
  subroutine init_next_iteration(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       if (this%channels(i)%current_iteration.lt.this%channels(i)%max_iterations) then
          call this%channels(i)%init_next_iteration()
!!$          this%channels(i)%done=.false.
!!$          this%channels(i)%integrals%done=.false.
       endif
    enddo
    call this%update_points_requested()
  end subroutine init_next_iteration
    
  subroutine check_if_done(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i,nevents_unw
    do i=1,this%nchannel
       nevents_unw=int(this%nevents_unw*this%channels(i)%res(1)/this%res(1))+1
       call this%channels(i)%check_generated_events(nevents_unw)
    enddo
  end subroutine check_if_done

  subroutine compute_total_rate(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,2
       this%res(i)=sum(this%channels(1:this%nchannel)%res(i))
       this%unc(i)=sqrt(sum(this%channels(1:this%nchannel)%unc(i)**2))
    enddo
  end subroutine compute_total_rate
  
  subroutine print_results(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       if (.not. this%channels(i)%evgen_done) call this%channels(i)%print_result_iteration()
       call this%channels(i)%print_combined_result()
    enddo
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    call flush()
  end subroutine print_results
  
  subroutine update_points_requested(this)
    implicit none
    class(integrator),intent(inout) :: this
    real(kind=8) :: total,total_channel
    integer :: i,j
    integer(kind=8) :: npoints,npoints_channel
    this%npoints_requested=this%npoints_requested*2
    npoints=0_8
    total=this%res(1)
    do i=1,this%nchannel
       npoints_channel=max(int(this%channels(i)%res(1)/total*dble(this%npoints_requested)),&
            min_points_per_channel)
       total_channel=this%channels(i)%res(1)
       do j=1,this%channels(i)%nintegral
          this%channels(i)%integrals(j)%npoints_requested=&
               max(int(this%channels(i)%integrals(j)%res(1)/total_channel*dble(npoints_channel)),&
               min_points_per_integral)
          npoints=npoints+this%channels(i)%integrals(j)%npoints_requested
       enddo
    enddo
    this%npoints_requested=npoints
  end subroutine update_points_requested
  
  subroutine channel_finalise_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    type(grid) :: new_grid
    integer :: i
    if (.not. this%evgen_done) call this%update_result_iteration()
    do i=1,this%ndim
       if (.not. this%evgen_done) then
          call this%grids(i,this%current_iteration)%update(new_grid)
          if (this%current_iteration .lt. this%max_iterations) then
             this%grids(i,this%current_iteration+1)=new_grid
          endif
       else
          this%grids(i,this%current_iteration+1)=this%grids(i,this%current_iteration)
       endif
    enddo
    if (.not. this%evgen_done) call this%combine_iterations()
  end subroutine channel_finalise_iteration

  subroutine channel_check_generated_events(this,nevents_unw)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: nevents_unw
    integer :: i
    do i=1,this%nintegral
       call this%integrals(i)%compute_fmax(this)
       call this%integrals(i)%unwgt()
       this%integrals(i)%nevents_unw=int(this%integrals(i)%res(1)/this%res(1)*nevents_unw)+1
       if (this%integrals(i)%n_unwgt.gt.this%integrals(i)%nevents_unw) then
          this%integrals(i)%evgen_done=.true.
       else
          this%integrals(i)%evgen_done=.false.
       endif
       if (this%current_iteration .lt. this%max_iterations) &
            call this%integrals(i)%compute_fmax_next_iteration(this)
    enddo
  end subroutine channel_check_generated_events

  subroutine integral_unwgt(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: j,iter
    this%n_unwgt=0
    do j=1,this%nevent_in_list
       iter=this%event_list(j)%iter
       this%event_list(j)%unwgt=this%event_list(j)%f_abs(iter).gt.this%f_max(iter)*this%event_list(j)%rnd
       if (this%event_list(j)%unwgt) this%n_unwgt=this%n_unwgt+1
    enddo
  end subroutine integral_unwgt
  
  subroutine integral_compute_fmax(this,thischan)
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevent,iter
    nevent=this%nevent_in_list
    do j=1,nevent
       iter=this%event_list(j)%iter
       if (iter.ne.this%current_iteration) cycle
       x=this%event_list(j)%x
       wgt=this%event_list(j)%wgt
       do k=iterations_without_events+1,this%current_iteration
          if (k.ne.iter) then
             call thischan%recompute_wgt_from_x(k,x,wgt_new)
             this%event_list(j)%f_abs(k)=this%event_list(j)%f_abs(iter)*wgt_new/wgt
          endif
          this%f_max(k)=max(this%f_max(k),this%event_list(j)%f_abs(k))
       enddo
    enddo
  end subroutine integral_compute_fmax

  subroutine integral_compute_fmax_next_iteration(this,thischan)
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevent,iter,next_iter
    next_iter=this%current_iteration+1
    nevent=this%nevent_in_list
    do j=1,nevent
       iter=this%event_list(j)%iter
       x=this%event_list(j)%x
       wgt=this%event_list(j)%wgt
       call thischan%recompute_wgt_from_x(next_iter,x,wgt_new)
       this%event_list(j)%f_abs(next_iter)=this%event_list(j)%f_abs(iter)*wgt_new/wgt
       this%f_max(next_iter)=max(this%f_max(next_iter),this%event_list(j)%f_abs(next_iter))
    enddo
  end subroutine integral_compute_fmax_next_iteration
  
  subroutine recompute_wgt_from_x(this,iter,x,wgt)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: iter
    real(kind=8),dimension(this%ndim),intent(in) :: x
    real(kind=8),intent(out) :: wgt
    integer :: i
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i,iter)%get_wgt(x(i),wgt)
    enddo
  end subroutine recompute_wgt_from_x

  subroutine get_wgt(this,x,wgt)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x
    real(kind=8),intent(inout) :: wgt
    real(kind=8) :: dx
    integer :: cell
    cell=1
    do while (this%current(cell).lt.x)
       cell=cell+1
    enddo
    dx=this%current(cell)-this%current(cell-1)
    wgt=wgt*dx*this%size
  end subroutine get_wgt
 
  subroutine channel_print_combined_result(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    do i=1,this%nintegral
       this%integrals(i)%npoints_nonzero_total=this%integrals(i)%npoints_nonzero_total+this%integrals(i)%npoints_nonzero
       if (this%integrals(i)%n_unwgt.gt.this%integrals(i)%nevents_unw) then
          write(*,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a)') &
               i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
               '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevent_in_list,&
               '--',this%integrals(i)%n_unwgt,'(',this%integrals(i)%nevents_unw,') DONE'
       else
          write(*,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a)') &
               i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
               '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevent_in_list,&
               '--',this%integrals(i)%n_unwgt,'(',this%integrals(i)%nevents_unw,')'
       endif
    enddo
  end subroutine channel_print_combined_result
    
  subroutine channel_print_result_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS:',this%res_iter(1),'+/-',this%unc_iter(1),'(',this%unc_iter(1)/this%res_iter(1)*100d0,'%)'
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel    :',this%res_iter(2),'+/-',this%unc_iter(2),'(',this%unc_iter(1)/this%res_iter(1)*100d0,'%)'
  end subroutine channel_print_result_iteration

  subroutine channel_combine_iterations(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    do i=1,this%nintegral
       if (.not. this%integrals(i)%evgen_done) &
            call this%integrals(i)%combine_iterations(this%current_iteration)
    enddo
    do i=1,2
       this%res(i)=sum(this%integrals(1:this%nintegral)%res(i))
       this%unc(i)=sqrt(sum(this%integrals(1:this%nintegral)%unc(i)**2))
    enddo
    if (this%current_iteration.eq.1) then
       this%npoints=this%npoints_iter
    else
       this%npoints=this%npoints+this%npoints_iter
    endif
  end subroutine channel_combine_iterations

  subroutine integral_combine_iterations(this,iter)
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
       this%npoints=this%npoints+this%npoints_iter
    endif
  end subroutine integral_combine_iterations
  
  subroutine update_res_and_unc(res,unc,npoints,res_iter,unc_iter,npoints_iter)
    implicit none
    real(kind=8),intent(inout) :: res,unc
    real(kind=8),intent(in) :: res_iter,unc_iter
    integer(kind=8),intent(inout) :: npoints,npoints_iter
    integer(kind=8) :: np
    np=npoints+npoints_iter
    unc=sqrt((unc**2*dble(npoints)**2+unc_iter**2*dble(npoints_iter)**2)&
         &/dble(np)**2+npoints*(res-res_iter)**2*dble(npoints_iter)&
         &**2/(dble(npoints_iter)*dble(np)**3))
    res=(npoints*res+npoints_iter*res_iter)/dble(np)
  end subroutine update_res_and_unc
  
  subroutine channel_update_result_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    do i=1,this%nintegral
       call this%integrals(i)%update_result_iteration()
    enddo
    do i=1,2
       this%res_iter(i)=sum(this%integrals(1:this%nintegral)%res_iter(i))
       this%unc_iter(i)=sqrt(sum(this%integrals(1:this%nintegral)%unc_iter(i)**2))
    enddo
  end subroutine channel_update_result_iteration

  subroutine compute_uncertainty(acc,acc2,np,unc)
    implicit none
    real(kind=8),intent(in) :: acc,acc2
    integer(kind=8),intent(in) :: np
    real(kind=8),intent(out) :: unc
    unc=sqrt(abs(acc2-acc**2)/dble(np))
  end subroutine compute_uncertainty
  
  subroutine integral_update_result_iteration(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: i
    if (this%npoints_iter.ne.0_8) then
       this%res_iter=this%accum/dble(this%npoints_iter)
       this%res2_iter=this%accum2/dble(this%npoints_iter)
       do i=1,2
          call compute_uncertainty(this%res_iter(i),this%res2_iter(i),this%npoints_iter,this%unc_iter(i))
       enddo
    endif
  end subroutine integral_update_result_iteration
  
  subroutine channel_add_point(this,x,wgt,cell,f_abs,f,to_write)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(in) :: x
    integer,dimension(this%ndim),intent(in) :: cell
    real(kind=8),intent(in) :: f_abs,f,wgt
    logical,intent(out) :: to_write
    integer :: i
    this%npoints_iter=this%npoints_iter+1
    do i=1,this%ndim
       call this%grids(i,this%current_iteration)%add_point(x(i),cell(i),f_abs)
    enddo
    call this%integrals(this%current_integral)%add_point(x,wgt,f_abs,f,to_write)
    if (all(this%integrals%done)) this%done=.true.
  end subroutine channel_add_point

  subroutine integral_add_point(this,x,wgt,f_abs,f,to_write)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,f,wgt
    real(kind=8),dimension(this%ndim),intent(in) :: x
    logical,intent(out) :: to_write
    logical :: enough
    this%npoints_iter=this%npoints_iter+1
    if (f_abs.ne.0d0) this%npoints_nonzero=this%npoints_nonzero+1
    this%accum(1)=this%accum(1)+f_abs
    this%accum(2)=this%accum(2)+f
    this%accum2(1)=this%accum2(1)+f_abs**2
    this%accum2(2)=this%accum2(2)+f**2
    call this%update_max_value(f_abs)
    call this%check_write_event(x,wgt,f_abs,to_write,enough)
    if (this%npoints_nonzero.ge.this%npoints_requested .or. enough) this%done=.true.
  end subroutine integral_add_point

  subroutine update_max_value(this,f_abs)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs
    this%max_value=max(this%max_value,f_abs)
  end subroutine update_max_value
  
  subroutine check_write_event(this,x,wgt,f_abs,to_write,enough)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,wgt
    real(kind=8),dimension(this%ndim),intent(in) :: x
    logical,intent(out) :: to_write,enough
    real(kind=8) :: rnd
    to_write=.false.
    enough=.false.
    if (this%current_iteration.le.iterations_without_events) return
    rnd=ran2()
    if (f_abs.gt.this%f_max(this%current_iteration)*rnd*write_event_fraction) then
       this%event=this%event+1
       this%nevent_in_list=this%nevent_in_list+1
       if (this%nevent_in_list.gt.size(this%event_list)) call this%increase_size_event_list()
       allocate(this%event_list(this%nevent_in_list)%f_abs(this%max_iterations))
       allocate(this%event_list(this%nevent_in_list)%x(this%ndim))
       this%event_list(this%nevent_in_list)%f_abs=0d0
       this%event_list(this%nevent_in_list)%f_abs(this%current_iteration)=f_abs
       this%event_list(this%nevent_in_list)%rnd=rnd
       this%event_list(this%nevent_in_list)%wgt=wgt
       this%event_list(this%nevent_in_list)%x=x
       this%event_list(this%nevent_in_list)%iter=this%current_iteration
       if (f_abs.gt.this%f_max(this%current_iteration)*rnd) this%n_unwgt=this%n_unwgt+1
!!$       write (*,*) this%n_unwgt,this%nevents_unw
       if (this%n_unwgt.gt.1.1d0*this%nevents_unw) enough=.true.
       if (f_abs.gt.this%f_max(this%current_iteration)) this%f_max(this%current_iteration)=f_abs
    endif
  end subroutine check_write_event
  
  subroutine increase_size_event_list(this)
    implicit none
    class(integral),intent(inout) :: this
    type(event),allocatable,dimension(:) :: tmp_list
    integer :: isize
    isize=size(this%event_list)
    allocate(tmp_list(2*isize))
    tmp_list(1:isize)=this%event_list(1:isize)
    deallocate(this%event_list)
    this%event_list=tmp_list
  end subroutine increase_size_event_list
  
  subroutine grid_add_point(this,x,cell,f_abs)
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x,f_abs
    integer,intent(in) :: cell
    if (importance_sampling_strategy.eq.1) then
       this%accum(cell)=this%accum(cell)+f_abs
    elseif (importance_sampling_strategy.eq.2) then
       this%accum(cell)=max(this%accum(cell),f_abs)
    elseif (importance_sampling_strategy.eq.3) then
       if (f_abs.gt.this%accum(cell)) then
          this%accum(cell)=this%accum(cell)+(f_abs-this%accum(cell))*0.1d0
       endif
    elseif (importance_sampling_strategy.eq.4) then
       if (this%accum(cell).eq.0d0) then
          this%accum(cell)=f_abs*1d-4
       elseif (f_abs.gt.this%accum(cell)) then
          this%accum(cell)=this%accum(cell)*1.1d0
       endif
    endif
    this%nhits(cell)=this%nhits(cell)+1
  end subroutine grid_add_point

  subroutine grid_update(this,new_grid)
    implicit none
    class(grid),intent(inout) :: this
    class(grid),intent(out) :: new_grid
    real(kind=8),dimension(0:this%size) :: current
    integer :: i,j
    real(kind=8) :: r
    call this%massage_accum()
    current(0)=0d0
    do i=1,this%size
       r=dble(i)/dble(this%size)
       do j=1,this%size
          if(r.lt.this%accum(j)) then
             current(i)=this%current(j-1)+(r-this%accum(j-1))/ &
                  (this%accum(j)-this%accum(j-1))*(this%current(j)-this%current(j-1))
             exit
          endif
       enddo
    enddo
    deallocate(this%accum)
    deallocate(this%nhits)
    current(this%size)=1d0
    call new_grid%init(current)
  end subroutine grid_update

  subroutine massage_accum(this)
    implicit none
    class(grid),intent(inout) :: this
    integer :: i
    real(kind=8) :: total
    real(kind=8), parameter :: tiny=1d-8
    do i=1,this%size
       if (this%nhits(i).eq.0) cycle
       if (importance_sampling_strategy.eq.1) then
          this%accum(i)=this%accum(i)/this%nhits(i)
       else
          this%accum(i)=this%accum(i)
       endif
    enddo
    total=sum(this%accum)
    do i=1,this%size
       if (this%accum(i).lt.1d-12*total) then
          this%accum(i)=0d0
       elseif (this%accum(i).lt.(1d0-1d-12)*total) then
          this%accum(i)=((this%accum(i)/total-1d0)/log(this%accum(i)/total))**1.5
       else
          this%accum(i)=1d0
       endif
       this%accum(i)=this%accum(i-1)+max(this%accum(i),0d0)
    enddo
    this%accum=this%accum/this%accum(this%size)
    ! make sure the elements are at least 'tiny' apart
    do i=1,this%size
       if (this%accum(i).lt.this%accum(i-1)+tiny) then
          this%accum(i)=this%accum(i-1)+tiny
       endif
    enddo
    this%accum(this%size)=1d0
    do i=this%size-1,1,-1
       if (this%accum(i).gt.this%accum(i+1)-tiny) then
          this%accum(i)=1d0-dble(i)*tiny
       else
          exit
       endif
    enddo
  end subroutine massage_accum

  subroutine channel_get_point(this,x,cell,wgt)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(out) :: x
    integer,dimension(this%ndim),intent(out) :: cell
    real(kind=8),intent(out) :: wgt
    integer :: i
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i,this%current_iteration)%get_x(x(i),cell(i),wgt)
    enddo
  end subroutine channel_get_point

  subroutine get_x(this,x,cell,wgt)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(out) :: x
    integer,intent(out) :: cell
    real(kind=8),intent(inout) :: wgt
    real(kind=8) :: rnd,dx
    rnd=this%size*ran2()
    cell=int(rnd)+1
    rnd=rnd-dble(cell-1)
    dx=this%current(cell)-this%current(cell-1)
    x=this%current(cell-1)+rnd*dx
    wgt=wgt*dx*this%size
  end subroutine get_x

  subroutine get_channel_and_integral(this,ichan,iint,wgt_chan)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(out) :: ichan,iint
    real(kind=8),intent(out) :: wgt_chan
    logical :: done
    do
       ichan=int(ran2()*this%nchannel)+1
       if (.not.(this%channels(ichan)%done.or.this%channels(ichan)%evgen_done)) exit
    enddo
    wgt_chan=1d0!dble(this%nchannel)
    done=.false.
    do
       iint=int(ran2()*this%channels(ichan)%nintegral)+1
       if (.not.(this%channels(ichan)%integrals(iint)%done.or.this%channels(ichan)%integrals(iint)%evgen_done)) exit
    enddo
    wgt_chan=wgt_chan*1d0!dble(this%channels(ichan)%nintegral)
    this%channels(ichan)%current_integral=iint
  end subroutine get_channel_and_integral
  
  subroutine compute_wgt_from_x(this,ichan,x,wgt)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: ichan
    real(kind=8),dimension(this%channels(ichan)%ndim),intent(in) :: x
    real(kind=8),intent(out) :: wgt
    call this%channels(ichan)%recompute_wgt_from_x(this%channels(ichan)%current_iteration,x,wgt)
  end subroutine compute_wgt_from_x
  
  subroutine read_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine read_all_grids
  
  subroutine write_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine write_all_grids
  
end module simple_integrator_mod
