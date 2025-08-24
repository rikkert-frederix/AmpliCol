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
     procedure,private :: recompute_wgt_from_x
  end type channel
  type :: integral
     real(kind=8) :: max_value
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
     procedure,private :: update_max_value,check_write_evnt,increase_size_evnt_list
     procedure,private :: init_next_iter => integral_init_next_iter
     procedure,private :: compute_fmax => integral_compute_fmax
     procedure,private :: compute_fmax_next_iter => integral_compute_fmax_next_iter
     procedure,private :: unwgt => integral_unwgt
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
     real(kind=8) :: wgt,rnd
     integer :: iter
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
     procedure,public :: init,get_points,fill_points,compute_wgt_from_x
     procedure,private :: read_all_grids,write_all_grids&
          &,get_channel_and_integral,update_points_requested&
          &,print_results,compute_total_rate&
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
  integer :: iters_without_evnts
  integer,parameter :: importance_sampling_strategy=3
  real(kind=8),parameter :: write_evnt_fraction=1d0
  integer,parameter :: min_points_per_channel=1024
  integer,parameter :: min_points_per_integral=128
  logical,parameter :: turn_off_evnt_generation=.false.
  real(kind=8),parameter :: required_accuracy_factor=10d0
  integer,parameter :: min_grid_size=8
  integer,parameter :: max_grid_size=2048
contains

  subroutine init(this,nchannel,ndim,nintegral,nevts_unw_req,niters)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: nchannel,nevts_unw_req,niters
    integer,dimension(nchannel),intent(in) :: ndim,nintegral
    integer :: i
    this%nchannel=nchannel
    this%nevts_unw_req=nevts_unw_req
    ! if we assume 0.3% unweighting efficiency, we expect ~10% time
    ! spend in iterations that do not produce events:
    iters_without_evnts=5
    this%npoints_requested=int(nevts_unw_req/(0.03*2**iters_without_evnts),kind=8)
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
         'iter',this%channels(1)%current_iter,'(',npoints_nonzero, &
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
    do i=1,this%nchannel
       this%channels(i)%nevts_unw_req=int(this%nevts_unw_req*this%channels(i)%res(1)/this%res(1))+1
       call this%channels(i)%check_gen_evnts()
    enddo
  end subroutine count_unweighted_evnts

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
    update_grids=((.not.this%evgen_done) .and. &
         this%unc(1)/this%res(1).gt.1d0/(sqrt(dble(this%nevts_unw_req))*required_accuracy_factor)) .or. &
         this%current_iter.le.iters_without_evnts
    if (.not.update_grids) then
       write (*,*) 'keeping grids fixed for channel',this%number
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
    do i=1,this%nintegral
       call this%integrals(i)%compute_fmax(this)
       call this%integrals(i)%unwgt()
       this%integrals(i)%nevts_unw_req=int(this%integrals(i)%res(1)/this%res(1)*this%nevts_unw_req)+1
       if (this%integrals(i)%nevts_unw_gen.gt.this%integrals(i)%nevts_unw_req) then
          this%integrals(i)%evgen_done=.true.
       else
          this%integrals(i)%evgen_done=.false.
       endif
    enddo
  end subroutine channel_check_gen_evnts

  subroutine integral_unwgt(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: j,iter
    this%nevts_unw_gen=0
    do j=1,this%nevnt_in_list
       iter=this%evnt_list(j)%iter
       this%evnt_list(j)%unwgt=this%evnt_list(j)%f_abs(iter).gt.this%f_max(iter)*this%evnt_list(j)%rnd
       if (this%evnt_list(j)%unwgt) this%nevts_unw_gen=this%nevts_unw_gen+1
    enddo
  end subroutine integral_unwgt
  
  subroutine integral_compute_fmax(this,thischan)
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevnt,iter
    nevnt=this%nevnt_in_list
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       if (iter.ne.this%current_iter) cycle
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       do k=iters_without_evnts+1,this%current_iter
          if (k.ne.iter) then
             call thischan%recompute_wgt_from_x(k,x,wgt_new)
             this%evnt_list(j)%f_abs(k)=this%evnt_list(j)%f_abs(iter)*wgt_new/wgt
          endif
          this%f_max(k)=max(this%f_max(k),this%evnt_list(j)%f_abs(k))
       enddo
    enddo
  end subroutine integral_compute_fmax

  subroutine integral_compute_fmax_next_iter(this,thischan)
    implicit none
    class(integral),intent(inout) :: this
    class(channel),intent(inout) :: thischan
    real(kind=8),dimension(this%ndim) :: x
    real(kind=8) :: wgt,wgt_new
    integer :: j,k,nevnt,iter,next_iter
    next_iter=this%current_iter
    nevnt=this%nevnt_in_list
    do j=1,nevnt
       iter=this%evnt_list(j)%iter
       x=this%evnt_list(j)%x
       wgt=this%evnt_list(j)%wgt
       call thischan%recompute_wgt_from_x(next_iter,x,wgt_new)
       this%evnt_list(j)%f_abs(next_iter)=this%evnt_list(j)%f_abs(iter)*wgt_new/wgt
       this%f_max(next_iter)=max(this%f_max(next_iter),this%evnt_list(j)%f_abs(next_iter))
    enddo
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
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    do i=1,this%nintegral
       this%integrals(i)%npoints_nonzero_total=this%integrals(i)%npoints_nonzero_total+this%integrals(i)%npoints_nonzero
       if (this%integrals(i)%nevts_unw_gen.gt.this%integrals(i)%nevts_unw_req) then
          write(*,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a)') &
               i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
               '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevnt_in_list,&
               '--',this%integrals(i)%nevts_unw_gen,'(',this%integrals(i)%nevts_unw_req,') DONE'
       else
          write(*,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a,1x,i10,1x,a)') &
               i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2),&
               '--',this%integrals(i)%npoints_nonzero_total,'--',this%integrals(i)%nevnt_in_list,&
               '--',this%integrals(i)%nevts_unw_gen,'(',this%integrals(i)%nevts_unw_req,')'
       endif
    enddo
  end subroutine channel_print_combined_result
    
  subroutine channel_print_result_iter(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS:',this%res_iter(1),'+/-',this%unc_iter(1),'(',this%unc_iter(1)/this%res_iter(1)*100d0,'%)'
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
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
    call this%update_max_value(f_abs)
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
    if (f_abs.gt.this%f_max(this%current_iter)*rnd*write_evnt_fraction) then
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
       if (f_abs.gt.this%f_max(this%current_iter)*rnd) this%nevts_unw_gen=this%nevts_unw_gen+1
       if (this%nevts_unw_gen.gt.1.1d0*this%nevts_unw_req) enough=.true.
       if (f_abs.gt.this%f_max(this%current_iter)) this%f_max(this%current_iter)=f_abs
    endif
  end subroutine check_write_evnt
  
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



