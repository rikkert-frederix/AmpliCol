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


module integrator_mod
  implicit none
  private
  type :: channel
     integer :: ndim,nintegral,current_integral,current_iteration&
          &,npoints,npoints_iter,number,max_iterations
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,unc_iter&
          &,chi2
     logical :: done,regrid
     type(grid),allocatable,dimension(:) :: grids
     type(integral),allocatable,dimension(:) :: integrals
   contains
     procedure,private :: init => channel_init
     procedure,private :: add_point => channel_add_point
     procedure,private :: get_jacobian => channel_get_jacobian
     procedure,private :: get_point => channel_get_point
     procedure,private :: finalise_iteration => channel_finalise_iteration
     procedure,private :: update_result_iteration => channel_update_result_iteration
     procedure,private :: combine_iterations => channel_combine_iterations
     procedure,private :: print_result_iteration => channel_print_result_iteration
     procedure,private :: print_combined_result => channel_print_combined_result
     procedure,private :: init_next_iteration => channel_init_next_iteration
  end type channel
  type :: integral
     real(kind=8) :: max_value,f_max
     real(kind=8),dimension(2) :: res,unc,res_iter,res2_iter,accum&
          &,accum2,unc_iter,chi2
     integer :: npoints_iter,npoints,npoints_requested&
          &,npoints_nonzero,over_wgt,event
     logical :: done
   contains
     procedure,private :: init => integral_init
     procedure,private :: add_point => integral_add_point
     procedure,private :: update_result_iteration => integral_update_result_iteration
     procedure,private :: combine_iterations => integral_combine_iterations
     procedure,private :: update_max_value,check_write_event
     procedure,private :: init_next_iteration => integral_init_next_iteration
  end type integral
  type :: grid
     integer :: size
     real(kind=8),allocatable,dimension(:) :: current,accum
     integer,allocatable,dimension(:) :: nhits
   contains
     procedure,private :: init => grid_init
     procedure,private :: add_point => grid_add_point
     procedure,private :: get_jacobian => grid_get_jacobian
     procedure,private :: get_x,massage_accum
     procedure,private :: init_next_iteration => grid_init_next_iteration
     procedure,private :: update => grid_update
  end type grid
  type,public :: integrator
     integer :: nchans,current_channel,npoints_generated,npoints_requested
     real(kind=8),dimension(2) :: res,unc
     real(kind=8),allocatable,dimension(:,:),public :: x
     real(kind=8),allocatable,dimension(:),public :: wgt
     integer,allocatable,dimension(:,:) :: cell
     type(channel),allocatable,dimension(:) :: channels
     ! keep track of maximum weights of all events written---integral-by-integral---and approx how many more we should generate in that channel.
   contains
     procedure,public :: init,get_points,fill_points,get_jacobian
     procedure,private :: read_all_grids,write_all_grids&
          &,get_channel_and_integral,update_points_requested&
          &,print_results
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
  integer,parameter :: importance_sampling_strategy=1
  integer,parameter :: iterations_for_regrid=1000
contains

  subroutine init(this,nchans,ndims,nintegrals,npoints,niters)
    implicit none
    class(integrator),intent(inout) :: this
    integer,intent(in) :: nchans,npoints,niters
    integer,dimension(nchans),intent(in) :: ndims,nintegrals
    integer :: i
    this%nchans=nchans
    allocate(this%channels(this%nchans))
    do i=1,this%nchans
       call this%channels(i)%init(ndims(i),nintegrals(i),npoints/nchans,niters)
       this%channels(i)%number=i
    enddo
    this%current_channel=0
    this%npoints_generated=0
    this%npoints_requested=npoints
    this%res=0d0
  end subroutine init

  subroutine channel_init(this,ndim,nintegral,npoints,niters)
    implicit none
    class(channel),intent(inout) :: this
    integer,intent(in) :: ndim,nintegral,npoints,niters
    integer :: i
    this%ndim=ndim
    this%nintegral=nintegral
    allocate(this%grids(1:this%ndim))
    allocate(this%integrals(1:this%nintegral))
    do i=1,this%ndim
       call this%grids(i)%init()
    enddo
    do i=1,this%nintegral
       call this%integrals(i)%init(npoints/this%nintegral)
    enddo
    this%current_integral=0
    this%current_iteration=1
    this%max_iterations=niters
    this%npoints=0
    this%res=0d0
    this%unc=0d0
    this%regrid=.true.
    call this%init_next_iteration()
  end subroutine channel_init

  subroutine channel_init_next_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    this%res_iter=0d0
    this%res2_iter=0d0
    this%unc_iter=0d0
    this%chi2=0d0
    this%npoints_iter=0
    this%done=.false.
    do i=1,this%nintegral
       call this%integrals(i)%init_next_iteration()
    enddo
    do i=1,this%ndim
       call this%grids(i)%init_next_iteration()
    enddo
  end subroutine channel_init_next_iteration
  
  subroutine integral_init_next_iteration(this)
    implicit none
    class(integral),intent(inout) :: this
    this%res_iter=0d0
    this%res2_iter=0d0
    this%accum=0d0
    this%accum2=0d0
    this%unc_iter=0d0
    this%npoints_iter=0
    this%npoints_nonzero=0
    this%max_value=0
    this%over_wgt=0
    this%event=0
    this%done=.false.
  end subroutine integral_init_next_iteration
  
  subroutine grid_init(this)
    implicit none
    class(grid),intent(inout) :: this
    integer :: i
    this%size=32
    allocate(this%current(0:this%size))
    do i=0,this%size
       this%current(i)=dble(i)/this%size
    enddo
    allocate(this%accum(0:this%size))
    allocate(this%nhits(this%size))
  end subroutine grid_init
  
  subroutine grid_init_next_iteration(this)
    implicit none
    class(grid),intent(inout) :: this
    this%accum=0d0
    this%nhits=0
  end subroutine grid_init_next_iteration
  
  subroutine integral_init(this,npoints)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: npoints
    this%npoints=0
    this%npoints_requested=npoints
    this%f_max=99d99
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
    done=.false.
    if (npoints.gt.this%npoints_generated) then
       write (*,*) 'ERROR: too many points returned'
       stop 1
    endif
    do i=1,npoints
       call this%channels(this%current_channel)%add_point(this%x(1,i),this%cell(1,i),f_abs(i),f(i),to_write(i))
    enddo
    this%npoints_generated=0
    if (all(this%channels%done)) then
       write (*,'(a,x,i4,x,a,x,i10,x,a)') &
            'iteration',this%channels(1)%current_iteration,'(',this%npoints_requested,'points) :'
       do i=1,this%nchans
          call this%channels(i)%finalise_iteration()
       enddo
       call this%print_results()
       write (*,*) ''
       call this%update_points_requested()
    endif
    if (all(this%channels%done)) done=.true.
    deallocate(this%x)
    deallocate(this%cell)
    deallocate(this%wgt)
  end subroutine fill_points

  subroutine print_results(this)
    implicit none
    class(integrator),intent(inout) :: this
    integer :: i
    do i=1,2
       this%res(i)=sum(this%channels(1:this%nchans)%res(i))
       this%unc(i)=sqrt(sum(this%channels(1:this%nchans)%unc(i)**2))
    enddo
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,a,1x,e12.6,1x,a,1x,e10.4,1x,a,f8.4,1x,a)') &
         'Integral     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
  end subroutine print_results
  
  subroutine update_points_requested(this)
    implicit none
    class(integrator),intent(inout) :: this
    real(kind=8) :: total
    integer :: i,j
    this%npoints_requested=this%npoints_requested*2
    total=sum(this%channels(1:this%nchans)%res(1))
    do i=1,this%nchans
       do j=1,this%channels(i)%nintegral
          this%channels(i)%integrals(j)%npoints_requested=&
               int(this%channels(i)%integrals(j)%res(1)/total*dble(this%npoints_requested))
       enddo
    enddo
  end subroutine update_points_requested
  
  subroutine channel_finalise_iteration(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    call this%update_result_iteration()
    if (this%current_iteration.ge.iterations_for_regrid) this%regrid=.false.
    if (this%regrid) then
       do i=1,this%ndim
          call this%grids(i)%update()
       enddo
    endif
    call this%print_result_iteration()
    call this%combine_iterations()
    call this%print_combined_result()
    if (this%current_iteration.lt.this%max_iterations) then
       if(this%regrid) call this%init_next_iteration()
       this%done=.false.
       this%integrals%done=.false.
       this%current_iteration=this%current_iteration+1
    endif
  end subroutine channel_finalise_iteration

  subroutine channel_print_combined_result(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel ABS (accum):',this%res(1),'+/-',this%unc(1),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(4x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,f7.3,1x,a)') &
         this%number,'channel     (accum):',this%res(2),'+/-',this%unc(2),'(',this%unc(1)/this%res(1)*100d0,'%)'
    write(*,'(24x,a,1x,f7.3)') 'chi2:',this%chi2(1)
    do i=1,this%nintegral
       write(*,'(23x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4)') &
            i,':',this%integrals(i)%res(2),'+/-',this%integrals(i)%unc(2)
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
    do i=1,this%nintegral
       write(*,'(15x,i4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,i10)') &
            i,':',this%integrals(i)%res_iter(2),'+/-',this%integrals(i)%unc_iter(2),&
            '--',this%integrals(i)%max_value,'--',this%integrals(i)%event
    enddo
  end subroutine channel_print_result_iteration
  
  subroutine channel_combine_iterations(this)
    implicit none
    class(channel),intent(inout) :: this
    integer :: i
    if (this%current_iteration.eq.1 .or. .not.this%regrid) then
       this%res=this%res_iter
       this%unc=this%unc_iter
       this%npoints=this%npoints_iter
    else
       do i=1,2
          call update_res_and_unc(this%res(i),this%unc(i),this%npoints,this%res_iter(i),this%unc_iter(i),this%npoints_iter)
       enddo
       this%npoints=this%npoints+this%npoints_iter
       this%chi2=this%chi2+(this%res_iter(1)-this%res(1))**2/this%unc_iter(1)**2
    endif
    do i=1,this%nintegral
       call this%integrals(i)%combine_iterations(this%current_iteration,this%regrid)
    enddo
  end subroutine channel_combine_iterations

  subroutine integral_combine_iterations(this,iter,regrid)
    implicit none
    class(integral),intent(inout) :: this
    integer,intent(in) :: iter
    logical,intent(in) :: regrid
    integer :: i
    if (iter.eq.1 .or. .not.regrid) then
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
    integer,intent(inout) :: npoints,npoints_iter
    integer :: np
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
    do i=1,2
       this%res_iter(i)=sum(this%integrals(1:this%nintegral)%accum(i))/dble(this%npoints_iter)
       this%res2_iter(i)=sum(this%integrals(1:this%nintegral)%accum2(i))/dble(this%npoints_iter)
       call compute_uncertainty(this%res_iter(i),this%res2_iter(i),this%npoints_iter,this%unc_iter(i))
    enddo
    do i=1,this%nintegral
       call this%integrals(i)%update_result_iteration()
    enddo
  end subroutine channel_update_result_iteration

  subroutine compute_uncertainty(acc,acc2,np,unc)
    implicit none
    real(kind=8),intent(in) :: acc,acc2
    integer,intent(in) :: np
    real(kind=8),intent(out) :: unc
    unc=sqrt(abs(acc2-acc**2)/dble(np))
  end subroutine compute_uncertainty
  
  subroutine integral_update_result_iteration(this)
    implicit none
    class(integral),intent(inout) :: this
    integer :: i
    this%res_iter=this%accum/dble(this%npoints_iter)
    this%res2_iter=this%accum2/dble(this%npoints_iter)
    do i=1,2
       call compute_uncertainty(this%res_iter(i),this%res2_iter(i),this%npoints_iter,this%unc_iter(i))
    enddo
    this%f_max=this%max_value
  end subroutine integral_update_result_iteration
  
  subroutine channel_add_point(this,x,cell,f_abs,f,to_write)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim) :: x
    integer,dimension(this%ndim) :: cell
    real(kind=8) :: f_abs,f
    logical,intent(out) :: to_write
    integer :: i
    this%npoints_iter=this%npoints_iter+1
    do i=1,this%ndim
       call this%grids(i)%add_point(x(i),cell(i),f_abs)
    enddo
    call this%integrals(this%current_integral)%add_point(f_abs,f,to_write)
    if (all(this%integrals%done)) this%done=.true.
  end subroutine channel_add_point

  subroutine integral_add_point(this,f_abs,f,to_write)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs,f
    logical,intent(out) :: to_write
    this%npoints_iter=this%npoints_iter+1
    if (f_abs.ne.0d0) this%npoints_nonzero=this%npoints_nonzero+1
    this%accum(1)=this%accum(1)+f_abs
    this%accum(2)=this%accum(2)+f
    this%accum2(1)=this%accum2(1)+f_abs**2
    this%accum2(2)=this%accum2(2)+f**2
    call this%update_max_value(f_abs)
    call this%check_write_event(f_abs,to_write)
    if (this%npoints_nonzero.ge.this%npoints_requested) this%done=.true.
  end subroutine integral_add_point

  subroutine update_max_value(this,f_abs)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs
    this%max_value=max(this%max_value,f_abs)
  end subroutine update_max_value
  
  subroutine check_write_event(this,f_abs,to_write)
    implicit none
    class(integral),intent(inout) :: this
    real(kind=8),intent(in) :: f_abs
    logical,intent(out) :: to_write
    if (f_abs.gt.this%f_max) then
       this%over_wgt=this%over_wgt+1
       write (*,*) 'found overweight',this%over_wgt,this%f_max,f_abs
    elseif (f_abs.gt.this%f_max*ran2()) then
       this%event=this%event+1
    endif
    to_write=.false.
  end subroutine check_write_event
  
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

  subroutine grid_get_jacobian(this)
    implicit none
    class(grid),intent(inout) :: this
  end subroutine grid_get_jacobian

  subroutine grid_update(this)
    implicit none
    class(grid),intent(inout) :: this
    real(kind=8),dimension(0:this%size) :: new_grid
    integer :: i,j
    real(kind=8) :: r
    call this%massage_accum()
    new_grid(0)=0d0
    do i=1,this%size
       r=dble(i)/dble(this%size)
       do j=1,this%size
          if(r.lt.this%accum(j)) then
             new_grid(i)=this%current(j-1)+(r-this%accum(j-1))/ &
                  (this%accum(j)-this%accum(j-1))*(this%current(j)-this%current(j-1))
             exit
          endif
       enddo
    enddo
    new_grid(this%size)=1d0
    this%current=new_grid
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

  subroutine channel_get_jacobian(this)
    implicit none
    class(channel),intent(inout) :: this
  end subroutine channel_get_jacobian

  subroutine channel_get_point(this,x,cell,wgt)
    implicit none
    class(channel),intent(inout) :: this
    real(kind=8),dimension(this%ndim),intent(out) :: x
    integer,dimension(this%ndim),intent(out) :: cell
    real(kind=8),intent(out) :: wgt
    integer :: i
    wgt=1d0
    do i=1,this%ndim
       call this%grids(i)%get_x(x(i),cell(i),wgt)
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
       ichan=int(ran2()*this%nchans)+1
       if (.not.this%channels(ichan)%done) exit
    enddo
    wgt_chan=1d0!dble(this%nchans)
    done=.false.
    do
       iint=int(ran2()*this%channels(ichan)%nintegral)+1
       if (.not.this%channels(ichan)%integrals(iint)%done) exit
    enddo
    wgt_chan=wgt_chan*1d0!dble(this%channels(ichan)%nintegral)
    this%channels(ichan)%current_integral=iint
  end subroutine get_channel_and_integral
  
  subroutine get_jacobian(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine get_jacobian
  
  subroutine read_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine read_all_grids
  
  subroutine write_all_grids(this)
    implicit none
    class(integrator),intent(inout) :: this
  end subroutine write_all_grids
  
end module integrator_mod
