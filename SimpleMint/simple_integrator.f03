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
     integer :: ndim,nintegral,current_integral,current_iter&
          &,number,max_iters,nevts_unw_req
     integer(kind=8) :: npoints,npoints_iter
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,unc_iter
     real(kind=8) :: overweight
     logical :: done,evgen_done
     type(grid),allocatable,dimension(:,:) :: grids
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
  type :: integral
     real(kind=8) :: max_value,overweight
     real(kind=8),dimension(:),allocatable :: f_max
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,accum&
          &,accum2,unc_iter
     integer :: ichan,nevts_unw_gen,evnt,nevnt_in_list,ndim &
          &,current_iter,max_iters,nevts_unw_req
     integer(kind=8) :: npoints_iter,npoints,npoints_requested&
          &,npoints_nonzero,npoints_nonzero_total
     logical :: done,evgen_done
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
  type :: grid
     integer :: size,size_fill
     real(kind=8),allocatable,dimension(:) :: current,accum,current_for_fillcell
     integer,allocatable,dimension(:) :: nhits
   contains
     procedure,private :: init => grid_init
     procedure,private :: add_point => grid_add_point
     procedure,private :: get_x,get_wgt,massage_accum,find_cell&
          &,interpolate_current,find_cell_to_fill
     procedure,private :: update => grid_update
  end type grid
  type :: evnt
     real(kind=8),allocatable,dimension(:) :: x,f_abs
     real(kind=8) :: wgt,rnd,overwgt
     integer :: iter,label
     logical :: unwgt
  end type evnt
  type,public :: integrator
     integer :: nchannel,current_channel,nevts_unw_req,npoints_gen
     integer(kind=8) :: npoints_requested
     real(kind=8),dimension(2) :: res,unc
     real(kind=8),allocatable,dimension(:,:),public :: x
     real(kind=8),allocatable,dimension(:),public :: wgt
     type(channel),allocatable,dimension(:) :: channels
     ! keep track of maximum weights of all evnts written---integral-by-integral---and approx how many more we should generate in that channel.
   contains
     procedure,public :: init,get_points,fill_points,compute_wgt_from_x,assign_evnt_wgts
     procedure,private :: read_all_grids,write_all_grids&
          &,get_channel_and_integral,update_points_requested&
          &,print_results,compute_total_rate,update_nevts_unw_req&
          &,count_unweighted_evnts,init_next_iter&
          &,get_npoints_nonzero_iter,finalise_iter,update_grids
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
  integer,save :: iters_without_evnts,evnt_label=0
  integer,parameter :: importance_sampling_strategy=3
  real(kind=8),parameter :: write_evnt_fraction=0.05d0 ! neglect write_evnt_fraction of largest weights to determine fmax for writing
  integer,parameter :: min_points_per_channel=1024
  integer,parameter :: min_points_per_integral=128
  logical,parameter :: turn_off_evnt_generation=.false.
  real(kind=8),parameter :: required_accuracy_factor=10d0
  integer,parameter :: min_grid_size=8
  integer,parameter :: max_grid_size=2048
  real(kind=8),parameter :: allowed_overweight_factor=0.001d0
  integer,parameter :: final_n_iters_for_evnt_gen=8
contains

  subroutine init(this,nchannel,ndim,nintegral,nevts_unw_req,niters)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: nchannel,nevts_unw_req,niters
    integer,dimension(nchannel),intent(in) :: ndim,nintegral
    integer :: i
    this%nchannel=nchannel
    this%nevts_unw_req=nevts_unw_req
    ! if we assume 1% unweighting efficiency, we expect ~10% time
    ! spend in iterations that do not produce events:
    iters_without_evnts=5
    this%npoints_requested=int(nevts_unw_req/(0.1d0*2**iters_without_evnts),kind=8)
    do while (this%npoints_requested/this%nchannel.lt.max(min_points_per_channel,min_points_per_integral*maxval(nintegral)) &
         .and. iters_without_evnts.gt.3)
       iters_without_evnts=iters_without_evnts-1
       this%npoints_requested=nevts_unw_req/(0.03*2**iters_without_evnts)
    enddo
    this%npoints_requested=max(this%npoints_requested,min_points_per_channel*this%nchannel,&
         min_points_per_integral*maxval(nintegral)*this%nchannel)
    allocate(this%channels(this%nchannel))
    do i=1,this%nchannel
       call this%channels(i)%init(ndim(i),nintegral(i),this%npoints_requested/nchannel,niters,i)
    enddo
    this%current_channel=0
    this%npoints_gen=0
    this%res=0d0
  end subroutine init

  subroutine channel_init(this,ndim,nintegral,npoints,niters,ichan)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: ndim,nintegral,niters,ichan
    integer(kind=8) :: npoints
    integer :: i
    this%ndim=ndim
    this%max_iters=niters
    this%nintegral=nintegral
    this%number=ichan
    allocate(this%grids(1:this%ndim,1:this%max_iters+1))
    allocate(this%integrals(1:this%nintegral))
    do i=1,this%ndim
       call this%grids(i,1)%init(npoints)
    enddo
    do i=1,this%nintegral
       call this%integrals(i)%init(ndim,npoints/this%nintegral,this%number,this%max_iters)
    enddo
    this%current_integral=0
    this%current_iter=0
    this%npoints=0_8
    this%res=0d0
    this%unc=0d0
    call this%init_next_iter()
    this%evgen_done=.false.
  end subroutine channel_init

  subroutine channel_init_next_iter(this)
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
       call this%integrals(i)%init_next_iter(this)
    enddo
    this%current_iter=this%current_iter+1
  end subroutine channel_init_next_iter
  
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
    this%npoints_iter=0_8
    this%npoints_nonzero=0_8
    this%evnt=0
    if (.not. this%evgen_done) this%done=.false.
    if (this%current_iter.eq.1) then
       this%f_max(this%current_iter)=-1d0
    elseif (this%current_iter.le.iters_without_evnts+1) then
       this%f_max(this%current_iter)=this%max_value
    else
       call this%compute_fmax_next_iter(thischan)
    endif
    this%max_value=0d0
  end subroutine integral_init_next_iter
  
  subroutine grid_init(this,npoints,current)
    implicit none
    class(grid),intent(inout) :: this
    integer(kind=8),intent(in) :: npoints
    real(kind=8),dimension(:),intent(in),optional :: current
    integer :: i,isize
    isize=size(current)-1
    this%size=max_grid_size
    call determine_sizefill(npoints,this%size_fill)
    allocate(this%current(0:this%size))
    allocate(this%current_for_fillcell(0:this%size_fill))
    if (present(current)) then
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
    allocate(this%accum(0:this%size_fill))
    allocate(this%nhits(this%size_fill))
    this%accum=0d0
    this%nhits=0
  end subroutine grid_init

  subroutine determine_sizefill(npoints,isize)
    implicit none
    integer(kind=8),intent(in) :: npoints
    integer,intent(out) :: isize
    isize=int(sqrt(dble(npoints))/10)
    isize=max(isize,min_grid_size)
    isize=min(isize,max_grid_size)
  end subroutine determine_sizefill
  
  subroutine integral_init(this,ndim,npoints,ichan,niters)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: ndim,ichan,niters
    integer(kind=8) :: npoints
    this%ndim=ndim
    this%ichan=ichan
    this%npoints=0_8
    this%npoints_requested=npoints
    this%max_iters=niters
    allocate(this%f_max(this%max_iters))
    this%f_max=-1d0
    allocate(this%evnt_list(npoints))
    this%nevnt_in_list=0
    this%evgen_done=.false.
    this%current_iter=0
    this%npoints_nonzero_total=0_8
    this%overweight=0d0
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
    allocate(this%wgt(1:npoints))

    do i=1,npoints
       call this%channels(this%current_channel)%get_point(this%x(1,i),this%wgt(i))
    enddo
    this%wgt=this%wgt*wgt_chan
    this%npoints_gen=npoints

  end subroutine get_points
  
  subroutine fill_points(this,npoints,f_abs,f,to_write,done)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: npoints
    real(kind=8),dimension(npoints),intent(in) :: f,f_abs
    logical,dimension(npoints),intent(out) :: to_write
    logical,intent(out) :: done
    integer :: i
    done=.false.
    if (npoints.gt.this%npoints_gen) then
       write (*,*) 'ERROR: too many points returned'
       stop 1
    endif
    do i=1,npoints
       call this%channels(this%current_channel)%add_point(this%x(1,i),this%wgt(i),f_abs(i),f(i),to_write(i))
    enddo
    this%npoints_gen=0
    if (all(this%channels%done)) then
       call this%finalise_iter(done)
    endif
    deallocate(this%x)
    deallocate(this%wgt)
  end subroutine fill_points

  subroutine finalise_iter(this,done)
    implicit none
    class(integrator),intent(inout) :: this
    logical,intent(out) :: done
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    character(len=19) :: formatted
    integer(kind=8) :: npoints_nonzero
    call date_and_time(date, time, zone)
    write(formatted, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2)') &
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
    call this%get_npoints_nonzero_iter(npoints_nonzero)
    write (*,*) ''
    write (*,'(a,x,i4,x,a,x,i10,x,a)') &
         'iteration',this%channels(1)%current_iter,'(',npoints_nonzero, &
         'points) '//trim(formatted)//' :'
    write (99,*) ''
    write (99,'(a,x,i4,x,a,x,i10,x,a)') &
         'iteration',this%channels(1)%current_iter,'(',npoints_nonzero, &
         'points) '//trim(formatted)//' :'
    call this%compute_total_rate()
    call this%count_unweighted_evnts()
    call this%print_results()
    call this%update_grids()
    call this%init_next_iter()
    if (all(this%channels%evgen_done)) done=.true.
    if (turn_off_evnt_generation .and. &
         this%unc(1)/this%res(1).lt.1d0/(sqrt(dble(this%nevts_unw_req))*required_accuracy_factor)) done=.true.
    if (all(this%channels%done)) done=.true.
    call flush(99)
  end subroutine finalise_iter

  subroutine update_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       call this%channels(i)%update_grids()
    enddo
  end subroutine update_grids
  
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
    
  subroutine count_unweighted_evnts(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    real(kind=8) :: nominal_evt_wgt
    call this%update_nevts_unw_req
    nominal_evt_wgt=this%res(1)/dble(this%nevts_unw_req)
    do i=1,this%nchannel
       call this%channels(i)%check_gen_evnts()
    enddo
  end subroutine count_unweighted_evnts

  subroutine update_nevts_unw_req(this)
    use sort_array_mod
    implicit none
    class(integrator),intent(inout) :: this
    real(kind=8),dimension(this%nchannel) :: res
    integer :: nevts_to_distribute,i
    integer,dimension(this%nchannel) :: idx
    real(kind=8) :: total
    res=this%channels%res(1)
    call sort_indices_by_values(res,idx)
    nevts_to_distribute=this%nevts_unw_req
    total=this%res(1)
    do i=1,this%nchannel-1
       this%channels(idx(i))%nevts_unw_req=int(nevts_to_distribute*this%channels(idx(i))%res(1)/total)
       if (ran2().lt.nevts_to_distribute*this%channels(idx(i))%res(1)/this%res(1)-&
                     this%channels(idx(i))%nevts_unw_req) then
          this%channels(idx(i))%nevts_unw_req=this%channels(idx(i))%nevts_unw_req+1
       endif
       nevts_to_distribute=nevts_to_distribute-this%channels(idx(i))%nevts_unw_req
       total=total-this%channels(idx(i))%res(1)
    enddo
    this%channels(idx(this%nchannel))%nevts_unw_req=nevts_to_distribute
    do i=1,this%nchannel
       call this%channels(i)%update_nevts_unw_req()
    enddo
  end subroutine update_nevts_unw_req

  subroutine channel_update_nevts_unw_req(this)
    use sort_array_mod
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%nintegral) :: res
    integer :: nevts_to_distribute,i
    integer,dimension(this%nintegral) :: idx
    real(kind=8) :: total
    res=this%integrals%res(1)
    call sort_indices_by_values(res,idx)
    nevts_to_distribute=this%nevts_unw_req
    total=this%res(1)
    do i=1,this%nintegral-1
       this%integrals(idx(i))%nevts_unw_req=int(nevts_to_distribute*this%integrals(idx(i))%res(1)/total)
       if (ran2().lt.nevts_to_distribute*this%integrals(idx(i))%res(1)/this%res(1)-&
                     this%integrals(idx(i))%nevts_unw_req) then
          this%integrals(idx(i))%nevts_unw_req=this%integrals(idx(i))%nevts_unw_req+1
       endif
       nevts_to_distribute=nevts_to_distribute-this%integrals(idx(i))%nevts_unw_req
       total=total-this%integrals(idx(i))%res(1)
    enddo
    this%integrals(idx(this%nintegral))%nevts_unw_req=nevts_to_distribute
  end subroutine channel_update_nevts_unw_req
  
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
  
  subroutine print_results(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,this%nchannel
       if (.not. this%channels(i)%evgen_done) call this%channels(i)%print_result_iter()
       call this%channels(i)%print_combined_result()
    enddo
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(99,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(99,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
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
    total=sum(this%channels%res(1),mask=.not.this%channels%evgen_done)
    do i=1,this%nchannel
       if (this%channels(i)%evgen_done) cycle
       npoints_channel=max(int(this%channels(i)%res(1)/total*dble(this%npoints_requested),kind=8),&
            min_points_per_channel)
       total_channel=sum(this%channels(i)%integrals%res(1),mask=.not.this%channels(i)%integrals%evgen_done)
       do j=1,this%channels(i)%nintegral
          if (this%channels(i)%integrals(j)%evgen_done) cycle
          this%channels(i)%integrals(j)%npoints_requested=&
               max(int(this%channels(i)%integrals(j)%res(1)/total_channel*dble(npoints_channel),kind=8),&
               min_points_per_integral)
          npoints=npoints+this%channels(i)%integrals(j)%npoints_requested
       enddo
    enddo
    this%npoints_requested=npoints
  end subroutine update_points_requested
  
  subroutine channel_update_grids(this)
    implicit none
    class(channel),intent(inout) :: this
    type(grid) :: new_grid
    integer :: i
    logical update_grids
    update_grids=(((.not.this%evgen_done) .and. &
         this%unc(1)/this%res(1).gt.1d0/(sqrt(dble(this%nevts_unw_req))*required_accuracy_factor)) .or. &
         this%current_iter.le.iters_without_evnts) .and. &
         this%npoints_iter.gt.int(this%npoints*0.2d0)
    if (.not.update_grids) then
       write (99,*) 'keeping grids fixed for channel',this%number
    endif
    do i=1,this%ndim
       if (update_grids) then
          call this%grids(i,this%current_iter)%update(this%npoints_iter,new_grid)
          if (this%current_iter .lt. this%max_iters) then
             this%grids(i,this%current_iter+1)=new_grid
          endif
       else
          this%grids(i,this%current_iter+1)=this%grids(i,this%current_iter)
       endif
    enddo
  end subroutine channel_update_grids

  subroutine channel_check_gen_evnts(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    logical :: done
    do i=1,this%nintegral
       call this%integrals(i)%compute_fmax(this)
       if (this%integrals(i)%nevts_unw_req.eq.0) then
          this%integrals(i)%evgen_done=.true.
          this%integrals(i)%nevts_unw_gen=0
          this%integrals(i)%overweight=0d0
          cycle
       endif
       call this%integrals(i)%unwgt()
       if (this%integrals(i)%nevts_unw_gen.gt.this%integrals(i)%nevts_unw_req) then
          call this%integrals(i)%check_overweight(done)
          if (done) then
             this%integrals(i)%evgen_done=.true.
          else
             this%integrals(i)%evgen_done=.false.
             this%integrals(i)%nevts_unw_gen=min(int(this%integrals(i)%nevts_unw_req*0.8d0),&
                  int(this%integrals(i)%nevts_unw_req*allowed_overweight_factor/this%integrals(i)%overweight))
          endif
       else
          this%integrals(i)%evgen_done=.false.
       endif
    enddo
    this%overweight=sum(this%integrals%overweight)
  end subroutine channel_check_gen_evnts

  
  subroutine check_overweight(this,done)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    logical,intent(out) :: done
    integer :: j,k
    real(kind=8),dimension(this%current_iter) :: fmax
    real(kind=8),dimension(this%nevnt_in_list) :: fabs
    real(kind=8),dimension(this%nevts_unw_req) :: fabs_top
    integer,dimension(this%nevts_unw_req) :: top_idx
    real(kind=8) :: fmax_req,tmp
    logical,dimension(this%current_iter) :: to_include
    ! check which iterations to include (only the final
    ! 'final_n_iters_for_evnt_gen' that generated events for this
    ! integral will be included)
    to_include=.false.
    k=0
    do j=this%nevnt_in_list,1,-1
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
       do k=iters_without_evnts+1,this%current_iter
          fmax(k)=max(fmax(k),this%evnt_list(j)%f_abs(k))
       enddo
    enddo
    ! rescale
    do j=1,this%nevnt_in_list
       if (to_include(this%evnt_list(j)%iter)) then
          fabs(j)=(this%evnt_list(j)%f_abs(this%evnt_list(j)%iter)/this%evnt_list(j)%rnd)/fmax(this%evnt_list(j)%iter)
       else
          fabs(j)=0d0
       endif
    enddo
    ! Take the nevts_unw_req largest
    k=this%nevts_unw_req
    call topk_largest(fabs,k,fabs_top,top_idx)
    ! find the fmax such that all remain
    fmax_req=fabs_top(k)
    ! check the overweight fraction
    this%overweight=0d0
    do j=1,this%nevts_unw_req
       this%evnt_list(top_idx(j))%unwgt=.true.
       tmp=this%evnt_list(top_idx(j))%f_abs(this%evnt_list(top_idx(j))%iter)/fmax(this%evnt_list(top_idx(j))%iter) 
       this%evnt_list(top_idx(j))%overwgt=tmp/fmax_req
       if (tmp.lt.fmax_req) cycle
       this%overweight=this%overweight+(tmp/fmax_req-1d0)
    enddo
    this%overweight=this%overweight/dble(this%nevts_unw_req)
    if (this%overweight.lt.allowed_overweight_factor) then
       done=.true.
    else
       done=.false.
    endif
  end subroutine check_overweight
  
  subroutine integral_unwgt(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: j,iter
    this%nevts_unw_gen=0
    do j=1,this%nevnt_in_list
       iter=this%evnt_list(j)%iter
       if (this%evnt_list(j)%f_abs(iter).gt.this%f_max(iter)*this%evnt_list(j)%rnd) &
            this%nevts_unw_gen=this%nevts_unw_gen+1
    enddo
  end subroutine integral_unwgt
  
  subroutine integral_compute_fmax(this,thischan)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevnt,iter
    integer,allocatable,dimension(:) :: index_fmax_top
    real(kind=8),allocatable,dimension(:) :: fmax_top
    real(kind=8),dimension(this%nevnt_in_list,this%current_iter) :: fabs
    nevnt=this%nevnt_in_list
    if (nevnt.eq.0) return
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       do k=iters_without_evnts+1,this%current_iter
          if ( (k.eq.iter .and. iter.ne.this%current_iter) .or. &
               (k.ne.iter .and. iter.eq.this%current_iter) ) then
             call thischan%recompute_wgt_from_x(k,x,wgt_new)
             this%evnt_list(j)%f_abs(k)=this%evnt_list(j)%f_abs(iter)*wgt_new/wgt
          endif
          fabs(j,k)=this%evnt_list(j)%f_abs(k)
       enddo
    enddo
    nevnt=max(int(write_evnt_fraction*nevnt),1)
    allocate(fmax_top(nevnt))
    allocate(index_fmax_top(nevnt))
    do k=iters_without_evnts+1,this%current_iter
       call topk_largest(fabs(:,k),nevnt,fmax_top,index_fmax_top)
       this%f_max(k)=fmax_top(nevnt)
    enddo
    deallocate(fmax_top)
    deallocate(index_fmax_top)
  end subroutine integral_compute_fmax

  subroutine integral_compute_fmax_next_iter(this,thischan)
    use topk_heap_mod
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevnt,iter,next_iter
    integer,allocatable,dimension(:) :: index_fmax_top
    real(kind=8),allocatable,dimension(:) :: fmax_top,fabs
    next_iter=this%current_iter
    nevnt=this%nevnt_in_list
    if (nevnt.le.200) then
       this%f_max(next_iter)=this%f_max(next_iter-1)
       return
    endif
    allocate(fabs(nevnt))
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       call thischan%recompute_wgt_from_x(next_iter,x,wgt_new)
       this%evnt_list(j)%f_abs(next_iter)=this%evnt_list(j)%f_abs(iter)*wgt_new/wgt
       fabs(j)=this%evnt_list(j)%f_abs(next_iter)
    enddo
    nevnt=max(int(write_evnt_fraction*nevnt),1)
    nevnt=max(int(nevnt*dble(thischan%max_iters-this%current_iter)/dble(thischan%max_iters)),1)
    allocate(fmax_top(nevnt))
    allocate(index_fmax_top(nevnt))
    call topk_largest(fabs,nevnt,fmax_top,index_fmax_top)
    this%f_max(next_iter)=fmax_top(nevnt)
    deallocate(fabs)
    deallocate(fmax_top)
    deallocate(index_fmax_top)
  end subroutine integral_compute_fmax_next_iter
  
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
    call this%find_cell(x,cell)
    dx=this%current(cell)-this%current(cell-1)
    wgt=wgt*dx*this%size
  end subroutine get_wgt

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
    cell=lo
  end subroutine find_cell
  
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
    cell=lo
  end subroutine find_cell_to_fill
  
  subroutine channel_print_combined_result(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    do i=1,this%nintegral
       this%integrals(i)%npoints_nonzero_total=this%integrals(i)%npoints_nonzero_total+this%integrals(i)%npoints_nonzero
       if (this%integrals(i)%nevts_unw_gen.ge.this%integrals(i)%nevts_unw_req) then
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
    
  subroutine channel_print_result_iter(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS:',this%res_iter(1),'+/-',this%unc_iter(1),'(',this%unc_iter(1)/this%res_iter(1)*100d0,'%)'
    write(99,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel    :',this%res_iter(2),'+/-',this%unc_iter(2),'(',this%unc_iter(1)/this%res_iter(1)*100d0,'%)'
  end subroutine channel_print_result_iter

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
    if (this%current_iter.eq.1) then
       this%npoints=this%npoints_iter
    else
       this%npoints=this%npoints+this%npoints_iter
    endif
  end subroutine channel_combine_iters

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
       this%npoints=this%npoints+this%npoints_iter
    endif
  end subroutine integral_combine_iters
  
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
  end subroutine channel_update_result_iter

  subroutine compute_uncertainty(acc,acc2,np,unc)
    implicit none
    real(kind=8),intent(in) :: acc,acc2
    integer(kind=8),intent(in) :: np
    real(kind=8),intent(out) :: unc
    unc=sqrt(abs(acc2-acc**2)/dble(np))
  end subroutine compute_uncertainty
  
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
    endif
  end subroutine integral_update_result_iter
  
  subroutine channel_add_point(this,x,wgt,f_abs,f,to_write)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(in) :: x
    real(kind=8),intent(in) :: f_abs,f,wgt
    logical,intent(out) :: to_write
    integer :: i
    this%npoints_iter=this%npoints_iter+1
    do i=1,this%ndim
       call this%grids(i,this%current_iter)%add_point(x(i),f_abs)
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
    if (this%current_iter.le.iters_without_evnts) call this%update_max_value(f_abs)
    call this%check_write_evnt(x,wgt,f_abs,to_write,enough)
    if (this%npoints_nonzero.ge.this%npoints_requested .or. enough) this%done=.true.
  end subroutine integral_add_point

  subroutine update_max_value(this,f_abs)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs
    this%max_value=max(this%max_value,f_abs)
  end subroutine update_max_value
  
  subroutine check_write_evnt(this,x,wgt,f_abs,to_write,enough)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,wgt
    real(kind=8),dimension(this%ndim),intent(in) :: x
    logical,intent(out) :: to_write,enough
    real(kind=8) :: rnd
    to_write=.false.
    enough=.false.
    if (turn_off_evnt_generation) return
    if (this%current_iter.le.iters_without_evnts) return
    rnd=ran2()
    if (f_abs.gt.this%f_max(this%current_iter)*rnd) then
       to_write=.true.
       evnt_label=evnt_label+1
       this%evnt=this%evnt+1
       this%nevnt_in_list=this%nevnt_in_list+1
       if (this%nevnt_in_list.gt.size(this%evnt_list)) call this%increase_size_evnt_list()
       allocate(this%evnt_list(this%nevnt_in_list)%f_abs(this%max_iters))
       allocate(this%evnt_list(this%nevnt_in_list)%x(this%ndim))
       this%evnt_list(this%nevnt_in_list)%f_abs=0d0
       this%evnt_list(this%nevnt_in_list)%f_abs(this%current_iter)=f_abs
       this%evnt_list(this%nevnt_in_list)%rnd=rnd
       this%evnt_list(this%nevnt_in_list)%wgt=wgt
       this%evnt_list(this%nevnt_in_list)%x=x
       this%evnt_list(this%nevnt_in_list)%iter=this%current_iter
       this%evnt_list(this%nevnt_in_list)%label=evnt_label
       this%nevts_unw_gen=this%nevts_unw_gen+1
       if (this%nevts_unw_gen.gt.1.5d0*this%nevts_unw_req) enough=.true.
    endif
  end subroutine check_write_evnt

  subroutine assign_evnt_wgts(this,wgts)
    implicit none
    class(integrator) :: this
    real(kind=8),allocatable,dimension(:,:),intent(out) :: wgts
    integer :: i,j
    real(kind=8) :: nominal_wgt
    allocate(wgts(3,evnt_label))
    nominal_wgt=this%res(1)
    do i=1,this%nchannel
       do j=1,this%channels(i)%nintegral
          call this%channels(i)%integrals(j)%compute_wgts(nominal_wgt,wgts)
       enddo
    enddo
  end subroutine assign_evnt_wgts

  subroutine compute_wgts(this,nominal_wgt,wgts)
    implicit none
    class(integral) :: this
    real(kind=8),dimension(3,evnt_label),intent(inout) :: wgts
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

  subroutine increase_size_evnt_list(this)
    implicit none
    class(integral),intent(inout) :: this
    type(evnt),allocatable,dimension(:) :: tmp_list
    integer :: isize
    isize=size(this%evnt_list)
    allocate(tmp_list(2*isize))
    tmp_list(1:isize)=this%evnt_list(1:isize)
    deallocate(this%evnt_list)
    this%evnt_list=tmp_list
  end subroutine increase_size_evnt_list
  
  subroutine grid_add_point(this,x,f_abs)
    class(grid),intent(inout) :: this
    real(kind=8),intent(in) :: x,f_abs
    integer :: cell
    call this%find_cell_to_fill(x,cell)
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
  
  subroutine interpolate_current(this,size_in,size_out,current_in,current_out)
    use pchip_uniform_strict
    implicit none
    class(grid),intent(inout) :: this
    integer,intent(in) :: size_in,size_out
    real(kind=8),dimension(0:size_in),intent(in) :: current_in
    real(kind=8),dimension(0:size_out),intent(out) :: current_out
    call resize_arr_pchip_strict(current_in,size_out,current_out)
  end subroutine interpolate_current
  
  subroutine massage_accum(this)
    implicit none
    class(grid),intent(inout) :: this
    integer :: i
    real(kind=8) :: total
    real(kind=8), parameter :: tiny=1d-8
    do i=1,this%size_fill
       if (this%nhits(i).eq.0) cycle
       if (importance_sampling_strategy.eq.1) then
          this%accum(i)=this%accum(i)/this%nhits(i)
       else
          this%accum(i)=this%accum(i)
       endif
    enddo
    total=sum(this%accum)
    do i=1,this%size_fill
       if (this%accum(i).lt.1d-12*total) then
          this%accum(i)=0d0
       elseif (this%accum(i).lt.(1d0-1d-12)*total) then
          this%accum(i)=((this%accum(i)/total-1d0)/log(this%accum(i)/total))**1.5
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
  end subroutine massage_accum

  subroutine channel_get_point(this,x,wgt)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(out) :: x
    real(kind=8),intent(out) :: wgt
    integer :: i
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i,this%current_iter)%get_x(x(i),wgt)
    enddo
  end subroutine channel_get_point

  subroutine get_x(this,x,wgt)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),intent(out) :: x
    integer :: cell
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
    call this%channels(ichan)%recompute_wgt_from_x(this%channels(ichan)%current_iter,x,wgt)
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



