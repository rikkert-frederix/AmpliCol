
program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_base
  use phase_space_gen23_mod
  use phase_space_genpt_mod
  use phase_space_haag_mod
  use math_functions
  use particles
  implicit none
  type(physics_model) :: phys_model
  integer :: next,iproc
  real(kind=8) :: weight
  integer :: j,c_o,i
  integer(kind=4),dimension(:),allocatable :: o,part
  real(kind=8),dimension(:),allocatable :: mass,width
  real(kind=8) :: s_cut(2),sqrts,evt_sign
  character(len=80) :: filename
  integer(kind=4) :: PS_choice,nquarks
  integer,parameter :: nevent_hel_filter=10
  integer,dimension(2) :: hel_picked
  integer :: iproc_picked,iproc_iden_picked
  integer :: ngroups,igroup,integration_step
  logical :: read_amps_from_file=.false.,write_amps_to_file=.false.
  logical :: read_proc_from_file
  integer :: nproc_unique
  integer,dimension(:,:),allocatable :: unique_procs
  integer :: nprocs
  integer,dimension(:),allocatable :: iden_iproc,phase_space_orders
  integer,dimension(:,:),allocatable :: processes,color_orders,multi_chans
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: factors
  
  type phase_space_order_group
     type(amplitude_QCD) :: amps
     class(phase_space_type),allocatable :: phase_space
     integer,dimension(:,:),allocatable :: multi_factor,processes,color_orders,multi_chans,multichans_all_chans_proc
     integer,dimension(:),allocatable :: iden_iproc,phase_space_orders,multichans_allchans,multichans_multifactor
     integer :: nproc,max_chans,multichans_nchan,multichans_nchans_proc
     real(kind=8),dimension(:,:),allocatable :: val_procs
     integer,dimension(:,:,:),allocatable :: iden_processes
     integer(kind=4),dimension(:,:),allocatable :: spin
     integer(kind=8),dimension(:),allocatable :: iden
     logical,dimension(-6:7,2) :: ipdgs
     integer(kind=4) :: nhel
     integer,dimension(:),allocatable :: col_fac
     real(kind=8),dimension(:),allocatable :: amp2,amp2_hel
     integer(kind=4),dimension(:),allocatable :: hel,hel_fac
     integer(kind=4) :: passed=0,all_evt=0
     integer,dimension(:),allocatable :: include_hel
  end type phase_space_order_group
  type(phase_space_order_group),dimension(:),allocatable :: pgl
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map

  call cpu_time(tTot_B)

  ! relevant input parameters for integration
  ! Number of events to generate. (If negative, start
  ! from a small number of points and double it each
  ! iteration. If positive, this is the number of
  ! points per iteration as well).
!!$  if (integration_step.eq.0 .or. integration_step.eq.2) then
     ncalls0=-100000
!!$  else
!!$     ncalls0=640000
!!$  endif

  itmax=16         ! Number of iterations. (If ncalls0 < 0, the
                   ! integration is aborted if accuracy (next line)
                   ! has been reached.

  ! setting energy
  sqrts=14000.d0

  ! cuts on invariants (Note: these assume massless particles!)
  s_cut(1)=0d0 ! cut on invariant between initial and final state particle
  s_cut(2)=0d0 ! cut on invariant of two final state particles.
  if (sqrt_s_min.gt.0d0) then
     s_cut(1)=max(s_cut(1),sqrt_s_min**2)
     s_cut(2)=max(s_cut(2),sqrt_s_min**2)
  endif
  if (pt_min.gt.0d0) then
     s_cut(1)=max(s_cut(1),pt_min**2)
     s_cut(2)=max(s_cut(2),2d0*pt_min**2*(1d0-cos(DRjj_min)))
  endif

  if (sqrt_s_min.gt.0d0) then
     s_cut(1:2)=sqrt_s_min**2
  endif

  if (include_pdf) call PDF_initialise

  call phys_model%init_part(173d0,1.491500d0)

  call get_run_arguments()

  if (integration_step.eq.0) then
     accuracy=0.001d0 ! Accuracy of the integration. (Ignored if ncalls0 > 0).
     write_amps_to_file=.true.
  else
     accuracy=max(1d0/sqrt(dble(abs(ncalls0))),0.0005d0)
     read_amps_from_file=.true.
  endif

  
  ! Not so relevant mint-module parameters: only used in special cases.
  call set_mint_module_special_parameters()
  nchans=ngroups ! overwrite the number of mint-channels to the number of
                 ! needed integration channels.

  call create_run_tag()
  
  if (read_amps_from_file .or. write_amps_to_file) then
       open(file='Outputs'//trim(adjustl(add_arg))//'/Res_files/amplitudes.bin',&
            unit=32,access='stream',form='unformatted',status='UNKNOWN')
  endif

  do igroup=1,ngroups
     if (pgl(igroup)%nproc.eq.0) cycle
     ! allocate the amplitudes and the phase-space for each of the integration channels
     if (PS_choice.eq.1) then
        allocate(phase_space_gen23 :: pgl(igroup)%phase_space)
     elseif  (PS_choice.eq.2) then
        allocate(phase_space_haag :: pgl(igroup)%phase_space)
     elseif (PS_choice.eq.3) then
        allocate(phase_space_genpt :: pgl(igroup)%phase_space)
     elseif (PS_choice.eq.4) then
        allocate(phase_space_gen23 :: pgl(igroup)%phase_space)
     endif

     allocate(mass(next))
     allocate(width(next))
     do i=1,next
        mass(i)=phys_model%get_mass(pgl(igroup)%processes(i,1))
        width(i)=phys_model%get_width(pgl(igroup)%processes(i,1))
        do iproc=2,pgl(igroup)%nproc
           if (mass(i).ne.phys_model%get_mass(pgl(igroup)%processes(i,iproc)) .or. &
                width(i).ne.phys_model%get_width(pgl(igroup)%processes(i,iproc))) then
              write (*,*) 'masses and widths not compatible among processes'
              stop 1
           endif
        enddo
     enddo
     call setup_spin(pgl(igroup))

     ! Initialise the phase-space parametrisation
     call cpu_time(tBefore)
     if (PS_choice.ge.1 .and. PS_choice.le.3) then
        call pgl(igroup)%phase_space%init(sqrts,next,mass,pgl(igroup)%phase_space_orders,&
             s_cut,pt_min,eta_max,DRjj_min,sqrt_s_min,.false.,include_pdf,&
             pgl(igroup)%iden_processes)
     elseif (PS_choice.eq.4) then
        call pgl(igroup)%phase_space%init(sqrts,next,mass,pgl(igroup)%phase_space_orders,&
             s_cut,pt_min,eta_max,DRjj_min,sqrt_s_min,.true.,include_pdf,&
             pgl(igroup)%iden_processes)
     endif
     call cpu_time(tAfter)
     t_PS_init=t_PS_init+tAfter-tBefore
     deallocate(mass)
     deallocate(width)

     allocate(pgl(igroup)%iden(pgl(igroup)%nproc))
     pgl(igroup)%iden(1:pgl(igroup)%nproc)=1
     call set_final_state_identical_particle_factor(pgl(igroup)) ! updates 'iden()'
     call set_initial_state_average_factor(pgl(igroup))          ! updates 'iden()'

     if (include_pdf) then
        call set_ipdgs_for_PDF(pgl(igroup))
     endif

     ! initialize the amplitudes. This creates the whole tree-structure from
     ! which the amps%evaluation() can compute the amplitudes for given
     ! phase-space points.
     call cpu_time(tBefore)
     if (read_amps_from_file) then
        call pgl(igroup)%amps%read_init_amps_from_file(next,32)
     else
        call pgl(igroup)%amps%init(1,next,pgl(igroup)%nproc,pgl(igroup)%processes,&
                pgl(igroup)%spin,pgl(igroup)%color_orders,phys_model,read_proc_from_file)
     endif

     !if (.not.read_proc_from_file .and. pgl(ichan)%amps%same_flav(3)) pgl(igroup)%nproc=3
     
     if (write_amps_to_file) then
        call pgl(igroup)%amps%write_init_amps_to_file(next,32)
     endif

     call cpu_time(tAfter)
     t_amp_init=t_amp_init+tAfter-tBefore

     ! Total number of amplitudes is stored in 'nhel'
     pgl(igroup)%nhel=pgl(igroup)%amps%n_amps

!     if (.not.read_proc_from_file.and.pgl(igroup)%amps%same_flav(3)) then
!        allocate(pgl(igroup)%col_fac(pgl(igroup)%amps%nprocs))
!     else
        allocate(pgl(igroup)%col_fac(pgl(igroup)%nproc))
!     endif

     call compute_LC_colour_factor(pgl(igroup))  ! updates 'col_fac()'

     allocate(pgl(igroup)%amp2(pgl(igroup)%nproc))

     ! number of helicities to sum over
     allocate(pgl(igroup)%amp2_hel(1:pgl(igroup)%nhel))
     allocate(pgl(igroup)%hel(1:next))
     allocate(pgl(igroup)%hel_fac(1:pgl(igroup)%nhel))
     pgl(igroup)%hel_fac(1:pgl(igroup)%nhel)=1

  enddo ! loop over phase-space-order groups
  
  if (read_amps_from_file .or. write_amps_to_file) then
     close(32)
  endif

  if (integration_step.le.1) then
     ! grid setup, or computation of upper bounding envelope
     call mint(integrand)
  else
     ! actual (unweighted) event generation
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise counters
     filename='Outputs'//trim(adjustl(add_arg))//'/events'//trim(adjustl(tag))//'.lhe'
     open(unit=11,file=filename,status='unknown')
     if (COMMAND_ARGUMENT_COUNT().le.10) call write_unique_in_file()
     do j=1,abs(ncalls0)
        call gen(integrand,1,2) ! generate an unweighted event
        call unwgt_process      ! pick a random process
        call unwgt_helicity     ! pick a random helicity for the process picked
        call write_event(11,sign(ans(1,0),evt_sign))
     enddo
     close(11)
     call gen(integrand,3,-1) ! print counters
  endif
     
  call cpu_time(tTot_a)
  t_all=tTot_a-tTot_b
  write(*,*) 'Time spent in phase-space initialisation:',t_PS_init 
  write(*,*) 'Time spent in amplitude initialisation',t_Amp_init
  write(*,*) 'Time spent in phase-space generation:',t_PS
  write(*,*) 'Time spent in amplitude evaluation',t_Amp
  write(*,*) 'Time spent in squaring amplitudes',t_mat
  write(*,*) 'Total time:',t_all
  write(*,*) 'Number of events:',pgl(1:ngroups)%all_evt
  write(*,*) 'Number passing cuts:',pgl(1:ngroups)%passed
  write(*,*) 'Fraction passing:',float(pgl(1:ngroups)%passed)/float(pgl(1:ngroups)%all_evt)
 
contains
  subroutine write_unique_in_file()
    implicit none
    integer :: iproc
    write(11,*) next,pgl_unique%nproc
    do iproc=1,pgl_unique%nproc
       write(11,*) unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:next,iproc)
    enddo
    deallocate(pgl_unique)
  end subroutine write_unique_in_file
  
  subroutine check_unique_processes()
    implicit none
    integer :: i,iproc,ievent,ih,nqq
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    real(kind=8),dimension(:),allocatable :: mass,width
    real(kind=8),dimension(ndim) :: x
    real(kind=8),external :: ran2
    allocate(phase_space_genpt :: pgl_unique%phase_space)
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(mass(next))
    allocate(width(next))
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes=unique_procs
    do i=1,next
        mass(i)=phys_model%get_mass(pgl_unique%processes(i,1))
        width(i)=phys_model%get_width(pgl_unique%processes(i,1))
        do iproc=2,pgl_unique%nproc
           if ( mass(i).ne.phys_model%get_mass(pgl_unique%processes(i,iproc)) .or. &
                width(i).ne.phys_model%get_width(pgl_unique%processes(i,iproc))) then
              write (*,*) 'masses and widths not compatible among processes'
              stop 1
           endif
        enddo
     enddo
     call setup_spin(pgl_unique)
     call setup_color_order(pgl_unique)

     do iproc=1,pgl_unique%nproc
        do i=1,2
              pgl_unique%processes(i,iproc)=phys_model%get_antipart(pgl_unique%processes(i,iproc))
        enddo
     enddo

     ! No multi-channel needed to check: simply use the color_orders for the phase-space order
     pgl_unique%phase_space_orders(1:next)=pgl_unique%color_orders(1:next,1)

     call pgl_unique%phase_space%init(sqrts,next,mass,pgl_unique%phase_space_orders,&
          s_cut,pt_min,eta_max,DRjj_min,sqrt_s_min,.false.,include_pdf, &
          pgl_unique%iden_processes)

     call pgl_unique%amps%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
             pgl_unique%spin,pgl_unique%color_orders,phys_model,read_proc_from_file)
     
     allocate(amp2(nevent,pgl_unique%nproc))

     ievent=0
     do while (ievent.lt.nevent)
        do i=1,ndim
           x(i)=ran2()
        enddo
        call pgl_unique%phase_space%generate_momenta(x)
        if (pgl_unique%phase_space%jac.lt.0d0) cycle
        ievent=ievent+1
        call pgl_unique%amps%evaluate(next,pgl_unique%phase_space%p,pgl_unique%hel,read_proc_from_file)
        iproc=0
        amp2(ievent,:)=0d0
        if (use_real_gluons .and. all(pgl_unique%amps%n_qqbar(1:pgl_unique%amps%nprocs).eq.0)) then
           do ih=1,pgl_unique%amps%n_amps
              do while (pgl_unique%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
              amp2(ievent,iproc)=amp2(ievent,iproc)+pgl_unique%amps%amps_r(ih)*pgl_unique%amps%amps_r(ih)
           enddo
        else
           do ih=1,pgl_unique%amps%n_amps
              do while (pgl_unique%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
              amp2(ievent,iproc)=amp2(ievent,iproc)+dble(pgl_unique%amps%amps(ih)*dconjg(pgl_unique%amps%amps(ih)))
           enddo
        endif
     enddo
     allocate(unique_map(1:pgl_unique%nproc))
     allocate(unique_map_value(1:pgl_unique%nproc))
     call find_unique(pgl_unique,nevent,amp2,unique_map,unique_map_value)

     do iproc=1,pgl_unique%nproc
        write (*,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:next,iproc),&
                ':',pgl_unique%color_orders(1:next,iproc)
     enddo

     deallocate(pgl_unique%spin)
     deallocate(pgl_unique%phase_space)
     deallocate(amp2)
   end subroutine check_unique_processes

   subroutine find_unique(pgl,nevent,amp2,unique_map,unique_map_value)
     implicit none
     type(phase_space_order_group) :: pgl
     integer :: nevent
     real(kind=8),dimension(nevent,pgl%nproc) :: amp2
     real(kind=8),dimension(pgl%nproc) :: unique_map_value
     integer,dimension(pgl%nproc) :: unique_map
     integer :: i,j,n
     real(kind=8),dimension(nevent) :: ratio
     real(kind=8) :: ave
     real(kind=8),parameter :: tiny=1d-6
     unique_map=-1d0
     do i=1,pgl%nproc
        do j=1,i-1
           if (all(amp2(1:nevent,j).eq.0d0)) cycle
           ratio(1:nevent)=amp2(1:nevent,i)/amp2(1:nevent,j)
           ave=sum(ratio(1:nevent))/nevent
           if (all(abs(ratio(1:nevent)/ave-1d0).lt.tiny)) then
              unique_map_value(i)=ave
              unique_map(i)=j
              exit
           endif
        enddo
        if (j.eq.i) then
           unique_map(i)=-1
           unique_map_value(i)=1d0
        endif
     enddo
   end subroutine find_unique
  
  real(kind=8) function integrand(x,vol,ifirst,f1)
    implicit none
    integer :: ifirst
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(nintegrals) :: f1
    real(kind=8), dimension(:),allocatable,save :: val,val_abs,vol_ichan
    real(kind=8),dimension(pgl(ichan)%nproc) :: colour_singlet_multichannel_weight
    integer :: ih,iproc,i
    real(kind=8) :: vol,cuts_wgt
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real(kind=4) :: tBefore,tAfter

    if (.not.allocated(val)) then
       allocate(val(1:maxval(pgl(1:ngroups)%nproc)))
       allocate(val_abs(1:maxval(pgl(1:ngroups)%nproc)))
       allocate(vol_ichan(1:ngroups))
    endif
    ! some point-by-point initialisation
    f1(1:nintegrals)=0d0
    if (pgl(ichan)%nproc.eq.0) return
    if (ifirst.eq.2) then
       ! use previously computed integrand
       f1(1)=sum(val_abs(1:pgl(ichan)%nproc))
       f1(2)=sum(val(1:pgl(ichan)%nproc))
       f1(3:pgl(ichan)%nproc+2)=val(1:pgl(ichan)%nproc)
       return
    endif
    new_point=.true.
    pass_cuts_check=.true.
    val_abs(1:pgl(ichan)%nproc)=0d0

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    call cpu_time(tBefore)
    call pgl(ichan)%phase_space%generate_momenta(x)
    call cpu_time(tAfter)
    t_PS= t_PS +tAfter-tBefore
    
    if (debug ) then
       write (*,*) pgl(ichan)%phase_space%jac
       stop 1
    endif
    
    pgl(ichan)%all_evt=pgl(ichan)%all_evt+1

    if (pgl(ichan)%phase_space%jac.lt.0d0) then
       val(1:pgl(ichan)%nproc)=0d0
       return
    endif

    cuts_wgt=pass_cuts(next,pgl(ichan)%phase_space%p)
    if ((pgl(ichan)%phase_space%jac.lt.0d0) .or. &
         (smooth_cuts .and. cuts_wgt.lt.0d0) .or. &
         (.not.smooth_cuts .and. cuts_wgt.lt.1d0)) then
       pass_cuts_check=.false.
       val(1:pgl(ichan)%nproc)=0d0
       return
    endif

    pgl(ichan)%passed = pgl(ichan)%passed + 1

    ! compute amplitudes
    call cpu_time(tBefore)

    call pgl(ichan)%amps%evaluate(next,pgl(ichan)%phase_space%p,pgl(ichan)%hel,read_proc_from_file)
    call cpu_time(tAfter)
    t_amp=t_amp+tAfter-tBefore
    
    call compute_multichannel_weight(ichan,pgl(ichan)%nproc,pgl(ichan)%max_chans,pgl(ichan)%multi_chans, &
                                     pgl(ichan)%phase_space%x,pgl(ichan)%phase_space%p,pgl(ichan)%phase_space%jac, &
                                     colour_singlet_multichannel_weight)
    
    call cpu_time(tBefore)
    iproc=0
    pgl(ichan)%amp2(:)=0d0
    if (use_real_gluons .and. all(pgl(ichan)%amps%n_qqbar(1:pgl(ichan)%amps%nprocs).eq.0)) then
       do ih=1,pgl(ichan)%amps%n_amps
          do while (pgl(ichan)%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          pgl(ichan)%amp2_hel(ih)=pgl(ichan)%amps%amps_r(ih)*pgl(ichan)%col_fac(iproc)*pgl(ichan)%amps%amps_r(ih)*&
                  pgl(ichan)%hel_fac(ih)
          pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
       enddo
    else
       do ih=1,pgl(ichan)%amps%n_amps
          do while (pgl(ichan)%amps%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          pgl(ichan)%amp2_hel(ih)=dble(pgl(ichan)%amps%amps(ih)*pgl(ichan)%col_fac(iproc)*dconjg(pgl(ichan)%amps%amps(ih)))*&
               pgl(ichan)%hel_fac(ih)
!          if (.not.read_proc_from_file .and. pgl(ichan)%amps%nprocs.eq.3 .and. &
!                  .not. pgl(ichan)%amps%same_flav(iproc)) then
!               cycle
!          elseif (.not.read_proc_from_file .and. pgl(ichan)%amps%nprocs.eq.3 .and. &
!                  pgl(ichan)%amps%same_flav(iproc)) then 
!               iproc=1
!          endif
          pgl(ichan)%amp2(iproc)=pgl(ichan)%amp2(iproc)+pgl(ichan)%amp2_hel(ih)
       enddo
    endif
    
    if (pgl(ichan)%passed.le.nevent_hel_filter) then
       call setup_helicity_filter(pgl(ichan))
       if (integration_step.eq.2 .and. pgl(ichan)%passed.eq.nevent_hel_filter) then
          ! since we update the helicities we need to compute when
          ! passed==nevent_hel_filter, the unweighting of the helicities goes
          ! wrong for this phase-space point. Hence, we need to skip it.
          pgl(ichan)%amp2(1:pgl(ichan)%nproc)=0d0
       endif
    endif

    ! MINT weight, phase-space jacobian and GeV -> pb conversion factor
    weight=vol*pgl(ichan)%phase_space%jac*conv

    ! multiply by the strong coupling
    if (pgl(ichan)%amps%n_sing(1).lt.next-2) then
       weight=weight*(4*pi*alphas)**(next-2-pgl(ichan)%amps%n_sing(1))
    endif
    
    ! multiply by the EW coupling
    if (pgl(ichan)%amps%n_sing(1).ge.1) then
       weight=weight*(2d0*4d0*pi*alphaEW)**pgl(ichan)%amps%n_sing(1)
    endif

    val(1:pgl(ichan)%nproc)=pgl(ichan)%amp2(1:pgl(ichan)%nproc)*weight/dble(pgl(ichan)%iden(1:pgl(ichan)%nproc))
    val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*colour_singlet_multichannel_weight(1:pgl(ichan)%nproc)
    
    ! Apply the weight from the cuts
    if (smooth_cuts) val(1:pgl(ichan)%nproc)=val(1:pgl(ichan)%nproc)*cuts_wgt

    call include_PDF_and_identical_procs(val,val_abs,pgl(ichan))

    ! pass the result to the mint module
    f1(1)=sum(val_abs(1:pgl(ichan)%nproc))
    f1(2)=sum(val(1:pgl(ichan)%nproc))
    f1(3:pgl(ichan)%nproc+2)=val(1:pgl(ichan)%nproc)
    
    integrand=f1(1)
    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore
  end function integrand

  double precision function pass_cuts(n,p)
    ! Cuts on the phase-space point. Note that these cuts need to be symmetric
    ! under pz -> -pz.
    implicit none
    integer :: i,j,n
    real(kind=8),dimension(0:3,n) :: p
    double precision :: frac,y,steep

    pass_cuts=1d0
    if (sqrt_s_min.gt.0d0) then
       do i=1,n-1
          do j=i+1,n
             if (abs(2d0*dot(p(0,i),p(0,j))).lt.sqrt_s_min**2) then
                pass_cuts=-1d0
                return
             endif
          enddo
       enddo
    endif

    do i=3,n
!!$       if (abs(part(i)).ge.0.and.abs(part(i)).le.6) then ! for quarks
         frac  = 0.9d0
         steep = 0.1d0
!!$       elseif (part(i).eq.21.or.part(i).eq.22) then ! for gluons and photons
!!$         frac  = 0.8d0
!!$         steep = 0.1d0
!!$       endif
         
       if (pt_min.gt.0d0) then
          if (pt(p(0,i)).lt.frac*pt_min) then
             pass_cuts=-1d0
             return
          endif
          if (pt(p(0,i)).gt.frac*pt_min.and.pt(p(0,i)).lt.pt_min) then
             y=(pt(p(0,i))-frac*pt_min)/(pt_min*(1d0-frac))
             if (integration_step.le.0) then
               !if (abs(part(i)).ge.0.and.abs(part(i)).le.6) then ! for quarks
                  pass_cuts=pass_cuts*((steep)*y/(steep+1d0-y)) ! 1/x damping function
               !elseif (part(i).eq.21.or.part(i).eq.22) then ! for gluons and photons
               !   pass_cuts=pass_cuts*((steep)*y/((steep+1d0-y)**2)) ! 1/x2 damping function
               !endif
             else
               pass_cuts=-1d0
             endif
             return
          endif
       endif

       if (eta_max.gt.0d0) then
          if (abs(eta(p(0,i))).gt.eta_max) then
             pass_cuts=-1d0
             return
          endif
       endif
       if (drjj_min.gt.0d0) then
          if (i.ne.n) then
             do j=i+1,n
                if (DeltaR(p(0,i),p(0,j)).lt.drjj_min) then
                   pass_cuts=-1d0
                   return
                endif
             enddo
          endif
       endif
    enddo
  end function pass_cuts

  subroutine compute_multichannel_weight(ichan,nproc,max_chans,chans,x,p,jac,weight)
    ! Computes the multichannel weight 'weight' when there are
    ! 'chans(0)' channels (that are listed in the array 'chans(1:)') and
    ! the current channel is 'ichan'. The momenta 'p' have been
    ! generate with phase-space jacobian 'jac' using the random
    ! variables 'x' within the channel 'ichan'.
    !
    ! weight = 1/( J_{ichan}*[ sum_{i=1}^{chans(0)} 1/J_i ] )
    !
    ! with J_i the combined Jacobian coming from MINT and the
    ! phase-space.
    implicit none
    integer,intent(in) :: ichan,nproc,max_chans
    integer,dimension(0:max_chans,nproc),intent(in) :: chans
    real(kind=8),dimension(0:3,next),intent(in) :: p
    real(kind=8),dimension(ndim),intent(in) :: x
    real(kind=8),intent(in) :: jac
    real(kind=8),dimension(pgl(ichan)%multichans_nchan) :: factors
    real(kind=8),dimension(pgl(ichan)%multichans_nchans_proc) :: weight_factors
    real(kind=8),dimension(nproc),intent(out) :: weight
    integer :: i,j,iproc,ii
    real(kind=8) :: vol_ichan,vol
    if (.not. use_colour_singlet_multichannel) then
       weight(1:nproc)=1d0/dble(chans(0,1:nproc))
       return
    endif
    call mint_get_jacobian_from_x(ichan,x,vol_ichan)
    do j=1,pgl(ichan)%multichans_nchan
       i=pgl(ichan)%multichans_allchans(j)
       if (i.eq.ichan) then
          ii=j
          cycle
       endif
       call pgl(i)%phase_space%compute_x_from_momenta(p)
       if (pgl(i)%phase_space%jac.lt.0d0) then
          ! The x's could not be correctly computed from the momenta
          write (*,*) 'WARNING: multi-channel weight not included'
          weight(1:nproc)=1d0/dble(chans(0,1:nproc))
          return
       endif
       call mint_get_jacobian_from_x(i,pgl(i)%phase_space%x,vol)
       factors(j)=pgl(i)%phase_space%jac*vol
    enddo
    do i=1,pgl(ichan)%multichans_nchans_proc
       weight_factors(i)=1d0
       do j=1,pgl(ichan)%multichans_all_chans_proc(0,i)
          if (pgl(ichan)%multichans_all_chans_proc(j,i).eq.ii) cycle
          weight_factors(i)=weight_factors(i)+jac*vol_ichan/factors(pgl(ichan)%multichans_all_chans_proc(j,i))
       enddo
    enddo
    weight(1:nproc)=1d0/weight_factors(pgl(ichan)%multichans_multifactor(1:nproc))
    
    
!!$    if (.not. use_colour_singlet_multichannel) then
!!$       weight(1:nproc)=1d0/dble(chans(0,1:nproc))
!!$       return
!!$    endif
!!$    do iproc=1,nproc
!!$       if (all(ichan.ne.chans(1:chans(0,iproc),iproc))) then
!!$          write (*,*) 'Current channel not among the multi-channel channels',ichan,iproc
!!$          write (*,*) chans(:,iproc)
!!$          stop 1
!!$       endif
!!$       if (chans(0,iproc).eq.1) then
!!$          weight(iproc)=1d0
!!$          cycle
!!$       endif
!!$       call mint_get_jacobian_from_x(ichan,x,vol_ichan)
!!$       weight(iproc)=1d0
!!$       do i=1,chans(0,iproc)
!!$          if (chans(i,iproc).eq.ichan) cycle
!!$          call pgl(chans(i,iproc))%phase_space%compute_x_from_momenta(p)
!!$          if (pgl(chans(i,iproc))%phase_space%jac.lt.0d0) then
!!$             ! The x's could not be correctly computed from the momenta
!!$             write (*,*) 'WARNING: multi-channel weight set to 1'
!!$             weight(iproc)=1d0
!!$             return
!!$          endif
!!$          call mint_get_jacobian_from_x(chans(i,iproc),pgl(chans(i,iproc))%phase_space%x,vol)
!!$          weight(iproc)=weight(iproc)+jac/pgl(chans(i,iproc))%phase_space%jac * vol_ichan/vol
!!$       enddo
!!$       weight(iproc)=1d0/weight(iproc)
!!$    enddo
  end subroutine compute_multichannel_weight

  subroutine setup_helicity_filter(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8) :: max_value
    integer :: ih1,ih2,iproc1,iproc2
    if (.not.allocated(pgl%include_hel)) then
       allocate(pgl%include_hel(pgl%nhel))
       pgl%include_hel(1:pgl%nhel)=0
    endif
    ! filter zero helicities and helicities that are identical
    max_value=maxval(pgl%amp2_hel(1:pgl%nhel))
    do ih1=1,pgl%nhel
       if (pgl%include_hel(ih1).ne.0) cycle
       if (pgl%amp2_hel(ih1)/max_value.gt.1d-10) then
          ! non-zero
          pgl%include_hel(ih1)=1
       else
          cycle
       endif
       do ih2=ih1+1,pgl%nhel
          if (abs(pgl%amp2_hel(ih1)-pgl%amp2_hel(ih2))/abs(pgl%amp2_hel(ih1)+pgl%amp2_hel(ih2)).lt.1d-10) then
             ! identical value. Now check that they belong to the same process
             iproc1=1; do while (iproc1.lt.pgl%nproc .and. (pgl%amps%iproc_start(iproc1+1)-ih1).le.0) ; iproc1=iproc1+1 ; enddo
             iproc2=1; do while (iproc2.lt.pgl%nproc .and. (pgl%amps%iproc_start(iproc2+1)-ih2).le.0) ; iproc2=iproc2+1 ; enddo
             if (iproc1.ne.iproc2) cycle
             ! identical process
             pgl%include_hel(ih2)=-ih1
             pgl%include_hel(ih1)=pgl%include_hel(ih1)+1
          endif
       enddo
    enddo

    if (pgl%passed.lt.nevent_hel_filter) return

    ih2=0
    do ih1=1,pgl%nhel
       if (pgl%include_hel(ih1).gt.0) ih2=ih2+1
    enddo

    call pgl%amps%filter_helicity(next,pgl%nhel,pgl%include_hel) ! this updates 'nhel' and 'include_hel'
    deallocate(pgl%hel_fac)
    allocate(pgl%hel_fac(pgl%nhel))
    pgl%hel_fac(1:pgl%nhel)=pgl%include_hel(1:pgl%nhel)
    deallocate(pgl%include_hel)
  end subroutine setup_helicity_filter

  real(kind=8) function pt(p)
    ! transverse momentum of 'p'
    implicit none
    real(kind=8), dimension(0:3) :: p
    pt=sqrt(p(1)**2+p(2)**2)
  end function pt
  
  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot

  real(kind=8) function eta(p)
    ! pseudo-rapidity of 'p'
    implicit none
    real(kind=8), dimension(0:3) :: p
    real(kind=8) :: theta
    theta=acos(p(3)/sqrt(p(1)**2+p(2)**2+p(3)**2))
    eta=-log(dtan(theta/2d0))
  end function eta

  real(kind=8) function delta_phi(p1,p2)
    ! azimuthal difference of 'p1' and 'p2'
    implicit none
    real(kind=8), dimension(0:3) :: p1,p2
    real(kind=8) :: denom
    denom=pt(p1)*pt(p2)
    delta_phi=acos((p1(1)*p2(1)+p1(2)*p2(2))/denom)
  end function delta_phi

  real(kind=8) function deltaR(p1,p2)
    ! Distance (Delta-R) between 'p1' and 'p2'
    implicit none
    real(kind=8), dimension(0:3) :: p1,p2
    deltaR=sqrt(delta_phi(p1,p2)**2+(eta(p1)-eta(p2))**2)
  end function deltaR

  subroutine write_event(iunit,wgt)
    implicit none
    integer :: i,iunit
    real(kind=8) :: wgt
    real(kind=8),external :: ran2
    write (iunit,*) '<event>'
    write (iunit,*) next,wgt!,pgl(ichan)%amp2(iproc_picked)*weight,pgl(ichan)%amp2(iproc_picked),weight
    write (iunit,'(100i3)') pgl(ichan)%amps%spins(1:next,hel_picked(1),hel_picked(2))
    if (.not.read_proc_from_file) iproc_picked=1
    write (iunit,'(100i3)') pgl(ichan)%color_orders(1:next,iproc_picked)
    ! Since some of the symmetry factors (in particular for gg->qqbar+ng)
    ! compensate for reducing the number of integration channels assuming
    ! symmetric initial states, we need to randomly flip all z-components in
    ! those cases. Easiest to always do this if the two incoming particles are
    ! identical.
    if (pgl(ichan)%processes(1,iproc_picked).ne.pgl(ichan)%processes(2,iproc_picked) .or. ran2().lt.0.5d0) then
       ! do not flip
       do i=1,next
          write (iunit,*) pgl(ichan)%iden_processes(i,iproc_iden_picked,iproc_picked),&
               pgl(ichan)%phase_space%p(1:3,i),pgl(ichan)%phase_space%p(0,i)
       enddo
    else
       ! do flip
       do i=1,next
          if (i.le.2) then
             write (iunit,*) pgl(ichan)%iden_processes(i,iproc_iden_picked,iproc_picked),&
                  pgl(ichan)%phase_space%p(1:2,3-i),-pgl(ichan)%phase_space%p(3,3-i),pgl(ichan)%phase_space%p(0,3-i)
          else
             write (iunit,*) pgl(ichan)%iden_processes(i,iproc_iden_picked,iproc_picked),&
                  pgl(ichan)%phase_space%p(1:2,i),-pgl(ichan)%phase_space%p(3,i),pgl(ichan)%phase_space%p(0,i)
          endif
       enddo
    endif
    write (iunit,*) '</event>'
  end subroutine write_event

  subroutine unwgt_process
    implicit none
    integer :: i,iproc
    real(kind=8) :: random,accum,target
    real(kind=8),external :: ran2
    target=0d0
    do iproc=1,pgl(ichan)%nproc
       do i=1,pgl(ichan)%iden_iproc(iproc)
          target=target+abs(pgl(ichan)%val_procs(i,iproc))
       enddo
    enddo
    random=ran2()*target
    iproc=1
    i=1
    accum=abs(pgl(ichan)%val_procs(i,iproc))
    do
       if (accum.gt.random) then
          exit
       else
          i=i+1
          if (i.gt.pgl(ichan)%iden_iproc(iproc)) then
             i=1
             iproc=iproc+1
          endif
          accum=accum+abs(pgl(ichan)%val_procs(i,iproc))
       endif
    enddo
    iproc_picked=iproc
    iproc_iden_picked=i
    if (iproc_picked.gt.pgl(ichan)%amps%nprocs) then
       write (*,*) "Could not unweight process",iproc_picked,pgl(ichan)%amps%nprocs
       stop 1
    endif
    if (iproc_iden_picked.gt.pgl(ichan)%iden_iproc(iproc)) then
       write (*,*) "Could not unweight process",iproc,iproc_iden_picked,pgl(ichan)%iden_iproc(iproc)
       stop 1
    endif
    if (pgl(ichan)%val_procs(iproc_iden_picked,iproc_picked).lt.0d0) then
       evt_sign=-1d0
    else
       evt_sign=+1d0
    endif
  end subroutine unwgt_process
  

  subroutine unwgt_helicity
    implicit none
    integer :: i
    real(kind=8) :: random
    real(kind=8),external :: ran2

    if (.not.read_proc_from_file .and. pgl(ichan)%amps%nprocs.eq.3) then
            do iproc=1,pgl(ichan)%amps%nprocs
               if (pgl(ichan)%amps%same_flav(iproc)) then
                       iproc_picked=iproc
                       exit
                endif
            enddo
    endif
    random=ran2()*pgl(ichan)%amp2(iproc_picked)
    i=pgl(ichan)%amps%iproc_start(iproc_picked)
    do
       if (pgl(ichan)%amp2_hel(i).gt.random) then
          exit
       else
          i=i+1
          pgl(ichan)%amp2_hel(i)=pgl(ichan)%amp2_hel(i)+pgl(ichan)%amp2_hel(i-1)
       endif
    enddo
    hel_picked(2)=i
    if ( hel_picked(2).lt.pgl(ichan)%amps%iproc_start(iproc_picked) .or. &
         hel_picked(2).ge.pgl(ichan)%amps%iproc_start(iproc_picked+1)) then
       write (*,*) 'Could not unweight helicity',hel_picked,iproc_picked,pgl(ichan)%amps%iproc_start(iproc_picked),&
            pgl(ichan)%amps%iproc_start(iproc_picked+1)
       stop 1
    endif
    if (pgl(ichan)%hel_fac(hel_picked(2)).gt.1) then
       hel_picked(1)=1+int(ran2()*pgl(ichan)%hel_fac(hel_picked(2)))
    else
       hel_picked(1)=1
    endif
  end subroutine unwgt_helicity
  
  subroutine get_run_arguments()
    implicit none
    integer :: argc,n_ps
    integer :: i,k,iproc
    character(len=256) :: argv
    logical :: found_1
    integer(kind=8) :: sym_fac
    integer(kind=8) iseed
    integer,dimension(:),allocatable :: process,order,ichans
    integer,dimension(:,:),allocatable :: ps_o
    integer :: factor,nproc_in_group,icheck,max_chans
    logical :: same_flavour
    character(len=1024) :: buff
    common /to_seed/iseed
    iseed=0
    ! integration steps:
    ! integration_step=0  (Setting up grids)
    ! integration_step=-1 (same as integration_step=0, but starting from existing grids)
    ! integration_step=1  (computing bounding envelope)
    ! integration_step=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    if (argc.eq.3) then
       read_proc_from_file=.true.
       do i=1,argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,'(a)') filename
          if (i.eq.2) read(argv,*) PS_choice
          if (i.eq.3) read(argv,*) integration_step
       enddo
       open(unit=10,file=filename,status='old')
       read (10,*) next,nproc_unique
       ndim=3*(next-2)-4
       if (include_pdf) ndim=ndim+2
       allocate(unique_procs(1:next,1:nproc_unique))
       do iproc=1,nproc_unique
          read(10,*) unique_procs(1:next,iproc)
       enddo
       allocate(pgl_unique)
       call check_unique_processes()
       read(10,*)
       read(10,*)
       read (10,*) ngroups
       allocate(pgl(ngroups))

       allocate(process(1:next))
       allocate(order(1:next))
       
       read (10,*) 
       do igroup=1,ngroups
          nprocs=0
          allocate(phase_space_orders(1:next))
          read(10,*) icheck,nproc_in_group,max_chans,phase_space_orders(1:next)
          if (icheck.ne.igroup) then
             write (*,*) 'ERROR in processes file',icheck,igroup
             stop 1
          endif
          allocate(iden_iproc(nproc_in_group))
          allocate(processes(1:next,nproc_in_group))
          allocate(color_orders(1:next,nproc_in_group))
          allocate(iden_processes(1:next,nproc_in_group,nproc_in_group))
          allocate(factors(nproc_in_group,nproc_in_group))
          allocate(multi_chans(0:max_chans,nproc_in_group))
          allocate(ichans(0:max_chans))
          do iproc=1,nproc_in_group
             read(10,'(a)') buff
             read(buff,*) ichans(0)
             read(buff,*) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),factor
             call add_to_process_list(process,order,factor,max_chans,ichans)
          enddo
          pgl(igroup)%nproc=nprocs
          pgl(igroup)%max_chans=max_chans
          allocate(pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%phase_space_orders(1:next))
          allocate(pgl(igroup)%multi_factor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%val_procs(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%multi_chans(0:max_chans,1:pgl(igroup)%nproc))
          pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc)=processes(1:next,1:pgl(igroup)%nproc)
          pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc)=color_orders(1:next,1:pgl(igroup)%nproc)
          pgl(igroup)%phase_space_orders(1:next)=phase_space_orders(1:next)
          pgl(igroup)%multi_factor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
               factors(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
          pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc)=iden_iproc(1:pgl(igroup)%nproc)
          pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
               iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
          pgl(igroup)%multi_chans(0:max_chans,1:pgl(igroup)%nproc)=multi_chans(0:max_chans,1:pgl(igroup)%nproc)
          deallocate(iden_iproc)
          deallocate(processes)
          deallocate(color_orders)
          deallocate(phase_space_orders)
          deallocate(iden_processes)
          deallocate(factors)
          deallocate(multi_chans)
          deallocate(ichans)
          write (*,*) '****************************************************'
          do iproc=1,pgl(igroup)%nproc
             write(*,*) iproc,':',pgl(igroup)%processes(1:next,iproc),' ; ',&
                  pgl(igroup)%color_orders(1:next,iproc),' ; ',pgl(igroup)%iden_iproc(iproc),' ; ',&
                  pgl(igroup)%multi_chans(1:pgl(igroup)%multi_chans(0,iproc),iproc)
          enddo
          write (*,*) '****************************************************'
          read(10,*)
          read(10,*)
          read(10,*)
       enddo
       if (integration_step.le.1) deallocate(pgl_unique)
999    continue
       close(10)
    elseif (argc.le.10) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'PS_choice, mode, next, *process*, *order*'
       stop 2
    else
       read_proc_from_file=.false.
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,*) PS_choice
          if (i.eq.2) read(argv,*) integration_step
          if (i.eq.3) then
             read(argv,*) next
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             ndim=3*(next-2)-4
             if (include_pdf) ndim=ndim+2
             allocate(part(1:next))
             allocate(o(1:next))
          endif
          do k=0,next-1
             if (i.eq.4+k) then
                read(argv,*) part(k+1)
             endif
          enddo
          do k=0,next-1
             if (i.eq.4+next+k) then
                read(argv,*) o(k+1)
             endif
          enddo
          if (argc.eq.3+2*next +1 .and. i.eq.argc) then
             ! Special case: we have an additional argument. Use it as a special tag
             read(argv,*) add_arg
             read (add_arg((index(add_arg,'S'))+1:(index(add_arg,'I'))-1),*,err=99) iseed
99           continue             
          endif
       enddo
       write (*,*) '******************************************'
       write (*,*) 'Process is     ',part
       write (*,*) 'Colour order is',o
       write (*,*) '******************************************'
       nquarks = 0
       do i=1,next
          if ((abs(part(i)).ge.1).and.(abs(part(i)).le.6)) then
             nquarks = nquarks + 1
          endif
       enddo
       if (nquarks.eq.0) then
          ! count the number of gluons between '1' and '2' in the colour order
          c_o=0
          i=0
          found_1=.false.
          do
             i=i+1
             if (i.gt.next) i=1
             if (found_1) c_o=c_o+1
             if (o(i).eq.1) then
                found_1=.true.
             endif
             if (o(i).eq.2 .and. found_1) then
                c_o=c_o-1
                exit
             endif
          enddo
       endif

       do i=2,next-1
          if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
              (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
             if (abs(part(o(i))).eq.abs(part(o(1))) .and. abs(part(o(i))).eq.abs(part(o(next)))) then
                same_flavour=.true.
             else
                same_flavour=.false.
             endif
             exit
          endif
       enddo

       if ((nquarks.ne.0 .and. nquarks.ne.2 .and. nquarks.ne.4) .or. (nquarks.gt.next)) then
          write (*,*) 'Not consistent number of external quarks (up to 2)',nquarks
          stop
       endif
       call compute_multichannel_symmetry_factor(sym_fac)

       call determine_multi_channel_size(part,n_ps)
       
       allocate(ps_o(1:next,1:n_ps))
       call determine_phase_space_orders(part,o,n_ps,ps_o)
       
       ngroups=n_ps
       allocate(pgl(ngroups))
       do i=1,ngroups
          pgl(i)%nproc=1
          allocate(pgl(i)%multi_chans(0:n_ps,1))
          pgl(i)%multi_chans(0,1)=n_ps
          do k=1,n_ps
             pgl(i)%multi_chans(k,1)=k
          enddo
          allocate(pgl(i)%processes(1:next,pgl(i)%nproc))
          pgl(i)%processes(1:next,1)=part(1:next)
          allocate(pgl(i)%color_orders(1:next,pgl(i)%nproc))
          pgl(i)%color_orders(1:next,1)=o(1:next)
          allocate(pgl(i)%phase_space_orders(1:next))
          pgl(i)%phase_space_orders(1:next)=ps_o(1:next,i)
          allocate(pgl(i)%multi_factor(1,1))
          pgl(i)%multi_factor(1,1)=sym_fac
          allocate(pgl(i)%iden_iproc(1))
          pgl(i)%iden_iproc(1)=1
          allocate(pgl(i)%val_procs(1,1))
          allocate(pgl(i)%iden_processes(1:next,1,1))
          pgl(i)%iden_processes(1:next,1,1)=pgl(i)%processes(1:next,1)
       enddo
    endif

    do i=1,ngroups
       call setup_optimised_multichannel_weight_computation(pgl(i))
    enddo
    
    
    ! basic checks:
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (integration_step.ne.0 .and. integration_step.ne.1 .and. integration_step.ne.2) then
       write (*,*) 'Incorrect integration_step',integration_step
       stop
    else
       imode=integration_step ! imode is used in MINT
    endif
    if (PS_choice.ne.1 .and. PS_choice.ne.2 .and. PS_choice.ne.3 .and. PS_choice.ne.4) then
       write (*,*) 'PS_Choice modes only 1, 2, 3 or 4',PS_choice
       stop
    endif
  end subroutine get_run_arguments

  subroutine setup_optimised_multichannel_weight_computation(pgl)
    implicit none
    type(phase_space_order_group) :: pgl
    integer,dimension(pgl%nproc*pgl%max_chans) :: all_chans
    integer,dimension(pgl%nproc*ngroups) :: all_chans_inv
    integer,dimension(0:pgl%max_chans,pgl%nproc) :: all_chans_proc
    integer,dimension(pgl%nproc) :: multifactor
    integer :: iproc,ichan,nchans,nchans_proc
    logical :: found
    nchans=0
    nchans_proc=0
    do iproc=1,pgl%nproc
       do ichan=1,pgl%multi_chans(0,iproc)
          if (all(all_chans(1:nchans).ne.pgl%multi_chans(ichan,iproc))) then
             nchans=nchans+1
             all_chans(nchans)=pgl%multi_chans(ichan,iproc)
             all_chans_inv(pgl%multi_chans(ichan,iproc))=nchans
          endif
       enddo
    enddo
    do iproc=1,pgl%nproc
       found=.false.
       do i=1,nchans_proc
          if (all_chans_proc(0,i).ne.pgl%multi_chans(0,iproc))cycle
          if (all(all_chans_inv(all_chans_proc(1:pgl%multi_chans(0,iproc),i)).eq. &
               pgl%multi_chans(1:pgl%multi_chans(0,iproc),iproc))) then
             multifactor(iproc)=i
             found=.true.
             exit
          endif
       enddo
       if (.not.found) then
          nchans_proc=nchans_proc+1
          all_chans_proc(0,nchans_proc)= &
               pgl%multi_chans(0,iproc)
          all_chans_proc(1:pgl%multi_chans(0,iproc),nchans_proc)= &
               all_chans_inv(pgl%multi_chans(1:pgl%multi_chans(0,iproc),iproc))
          multifactor(iproc)=nchans_proc
       endif
    enddo

    pgl%multichans_nchan=nchans
    pgl%multichans_nchans_proc=nchans_proc
    allocate(pgl%multichans_allchans(nchans))
    pgl%multichans_allchans=all_chans(1:nchans)
    allocate(pgl%multichans_all_chans_proc(0:pgl%max_chans,nchans_proc))
    pgl%multichans_all_chans_proc(0:pgl%max_chans,1:nchans_proc)=all_chans_proc(0:pgl%max_chans,1:nchans_proc)
    allocate(pgl%multichans_multifactor(pgl%nproc))
    pgl%multichans_multifactor(1:pgl%nproc)=multifactor(1:pgl%nproc)

  end subroutine setup_optimised_multichannel_weight_computation
    

  
  subroutine determine_multi_channel_size(part,n_ps)
    implicit none
    integer,dimension(1:next),intent(in) :: part
    integer,intent(out) :: n_ps
    integer :: i
    if (.not. use_colour_singlet_multichannel) then
       n_ps=1
       return
    endif
    n_ps=0
    do i=1,next
       if (is_singlet(part(i))) n_ps=n_ps+1
    enddo
    n_ps=factorial(n_ps)
  end subroutine determine_multi_channel_size
  
  subroutine determine_phase_space_orders(part,col_o,n_ps,PS_o)
    implicit none
    integer,dimension(1:next),intent(in) :: part,col_o
    integer,intent(in) :: n_ps
    integer,dimension(1:next,1:n_ps),intent(inout) :: ps_o
    integer :: i,n_sing,i_sing
    integer,dimension(:),allocatable :: jmap,ips,ips_out
    if (n_ps.eq.1) then
       PS_o(1:next,1)=col_o(1:next)
       return
    endif
    n_sing=0
    do i=1,next
       if (is_singlet(part(i))) n_sing=n_sing+1
    enddo
    allocate(jmap(n_sing))
    j=0
    do i=1,next
       if (is_singlet(part(o(i)))) then
          j=j+1
          jmap(j)=i
       endif
    enddo
    allocate(ips(n_sing))
    allocate(ips_out(n_sing))
    do i=1,n_sing
       ips(i)=i
    enddo
    do j=1,n_ps
       i_sing=0
       do i=1,next
          if (is_singlet(part(o(i)))) then
             i_sing=i_sing+1
             ps_o(i,j)=col_o(jmap(ips(i_sing)))
          else
             ps_o(i,j)=col_o(i)
          endif
       enddo
       if (j.eq.n_ps) exit
       call get_next_iperm(n_sing,ips,ips_out,n_sing)
       ips(1:n_sing)=ips_out(1:n_sing)
    enddo
  end subroutine determine_phase_space_orders

  subroutine add_to_process_list(process,order,factor,max_chans,ichans)
    implicit none
    integer :: max_chans
    integer,dimension(0:max_chans) :: ichans
    integer,dimension(next) :: process,order,process_mapped,process_unique,mapping
    integer :: factor
    real(kind=8) :: multi_factor
    call map_to_canonical_form(process,process_mapped,mapping)
    call get_unique_process(process,process_mapped,process_unique,multi_factor,mapping)
    multi_factor=multi_factor*factor
    call add_to_unique_process_list(process,process_unique,order,multi_factor,max_chans,ichans)
  end subroutine add_to_process_list

  subroutine move_colour_singlet_in_order(process,order)
    ! Move the colour singlet(s) to go *before* the final anti-quark in the order
    implicit none
    integer,dimension(next),intent(in) :: process
    integer,dimension(next),intent(inout) :: order
    integer :: i,iord,aq,iaq,itmp,ipart
    ! find the final anti-quark
    do i=next,1,-1
       iord=order(i)
       ipart=process(iord)
       if (((iord.le.2 .and. is_quark(ipart)).or.(iord.gt.2 .and. is_antiquark(ipart)))) then
          aq=i
          iaq=iord
          exit
       endif
    enddo
    ! move the colour_singlet to after the anti-quark
    i=1
    do while (i.le.next)
       iord=order(i)
       ipart=process(iord)
       if(is_singlet(ipart)) then
          if (i.gt.aq) then
             order(aq:i)=[order(i),order(aq)]
             aq=aq+1
             i=i+1
          else
             order(i:aq)=[order(i+1:aq),order(i)]
             aq=aq-1
          endif
       else
          i=i+1
       endif
    enddo
  end subroutine move_colour_singlet_in_order
          
  subroutine add_to_unique_process_list(process,process_unique,order,multi_factor,max_chans,ichans)
    implicit none
    integer,intent(in) :: max_chans
    integer,dimension(0:max_chans),intent(in) :: ichans
    integer,dimension(next) :: process,process_unique,order
    real(kind=8) :: multi_factor
    integer :: iproc
    call move_colour_singlet_in_order(process,order)
    do iproc=1,nprocs
       if (all(process_unique(1:next).eq.processes(1:next,iproc)) &
            .and. all(order(1:next).eq.color_orders(1:next,iproc))) exit
    enddo
    if (iproc.gt.nprocs) then
       ! new matrix element to generate
       nprocs=nprocs+1
       processes(1:next,iproc)=process_unique(1:next)
       color_orders(1:next,iproc)=order(1:next)
       iden_iproc(iproc)=1
       iden_processes(1:next,iden_iproc(iproc),iproc)=process(1:next)
       factors(iden_iproc(iproc),iproc)=multi_factor
       multi_chans(0:ichans(0),iproc)=ichans(0:ichans(0))
    else
       ! identical to another matrix element
       iden_iproc(iproc)=iden_iproc(iproc)+1
       iden_processes(1:next,iden_iproc(iproc),iproc)=process(1:next)
       factors(iden_iproc(iproc),iproc)=multi_factor
       if (ichans(0).ne.multi_chans(0,iproc)) then
          write (*,*) 'Number of multi-channels not the same among identical processes',&
               ichans(0),multi_chans(0,iproc)
          stop 1
       endif
       if (any(ichans(1:ichans(0)).ne.multi_chans(1:ichans(0),iproc))) then
          write (*,*) 'Multi-channels not the same among identical processes'
          write (*,*) ichans(1:ichans(0))
          write (*,*) multi_chans(1:ichans(0),iproc)
          stop 1
       endif
    endif
  end subroutine add_to_unique_process_list

  
  subroutine get_unique_process(process,process_mapped,process_unique,multi_factor,mapping)
    implicit none
    integer,dimension(next) :: process,process_mapped,process_unique,mapping
    integer :: iproc,map_from,map_to,i,j
    real(kind=8) :: multi_factor
    do iproc=1,pgl_unique%nproc
       if (all(process_mapped(1:next).eq.pgl_unique%processes(1:next,iproc))) exit
    enddo
    if (iproc.gt.pgl_unique%nproc) then
       write (*,*) 'Process not found'
       write (*,*) process(1:next)
       write (*,*) process_mapped(1:next)
       stop 1
    endif
    multi_factor=unique_map_value(iproc)
    if (unique_map(iproc).gt.0) then
       do i=1,next
          if ((i.le.2 .and. mapping(i).gt.2) .or. (i.gt.2 .and. mapping(i).le.2)) then
             process_unique(mapping(i))=phys_model%get_antipart(pgl_unique%processes(i,unique_map(iproc)))
          else
             process_unique(mapping(i))=pgl_unique%processes(i,unique_map(iproc))
          endif
       enddo
    else
       do i=1,next
          if ((i.le.2 .and. mapping(i).gt.2) .or. (i.gt.2 .and. mapping(i).le.2)) then
             process_unique(mapping(i))=phys_model%get_antipart(pgl_unique%processes(i,iproc))
          else
             process_unique(mapping(i))=pgl_unique%processes(i,iproc)
          endif
       enddo
    endif
  end subroutine get_unique_process


  subroutine sort_with_mapping(n,array,mapping)
    !
    ! EXAMPLE:
    !
    ! input:
    ! n=5
    ! array=[4, 1, 8, 2, 3]
    !
    ! output:
    ! array=[1, 2, 3, 4, 8]
    ! mapping=[2, 4, 5, 1, 3]
    !
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(inout) :: array
    integer,dimension(n),intent(out) :: mapping
    integer :: i, j, temp
    ! Initialize mapping
    mapping = [(i,i=1,n)]
    ! Sort the array and mapping using a simple bubble sort
    do i=1,n-1
       do j=1,n-i
          if (array(j) .gt. array(j+1)) then
             ! Swap array elements
             temp = array(j)
             array(j) = array(j+1)
             array(j+1) = temp
             ! Swap mapping
             temp = mapping(j)
             mapping(j) = mapping(j+1)
             mapping(j+1) = temp
          endif
       enddo
    enddo
  end subroutine sort_with_mapping

  
  subroutine map_to_canonical_form(process,part,mapping)
    ! cross the two initial state particle PDGs, order according to
    ! the PDG value, (and reflip the two initial states again)
    implicit none
    integer,dimension(next) :: process,part,mapping
    part(1:next)=process(1:next)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
    call sort_with_mapping(next,part,mapping)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
  end subroutine map_to_canonical_form
  
  
  subroutine setup_spin(pgl)
    ! Use the first process in the processes() array to setup all the possible
    ! spin states. Note that this assumes that all the processes() have the
    ! same number of spin states
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc
    if (.not. allocated(pgl%spin)) allocate(pgl%spin(0:3,1:next))
    do i=1,next
       pgl%spin(0,i)=phys_model%get_spin(pgl%processes(i,1))
       if (pgl%spin(0,i).eq.2) then
          pgl%spin(1,i)=-1
          pgl%spin(2,i)=1
       else
          write (*,*) 'spin state not known',i,pgl%processes(i,1),pgl%spin(0,i)
          stop 1
       endif
    enddo
    do iproc=2,pgl%nproc
       do i=1,next
          if (pgl%spin(0,i).ne.phys_model%get_spin(pgl%processes(i,iproc))) then
             write (*,*) 'Spin states of particles in different processes not compatible',iproc
             stop 1
          endif
       enddo
    enddo
  end subroutine setup_spin

  subroutine setup_color_order(pgl_unique)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl_unique
    integer :: i,iproc,nq,ng,nsing,iq,iaq,is,ig
    do iproc=1,pgl_unique%nproc
       nq=0
       ng=0
       nsing=0
       do i=1,next
          if (is_quark(abs(pgl_unique%processes(i,iproc)))) then
             nq=nq+1
          elseif(is_gluon(pgl_unique%processes(i,iproc))) then
             ng=ng+1
          elseif(is_singlet(pgl_unique%processes(i,iproc))) then
             nsing=nsing+1
          else
             write (*,*) 'unknown particle type:',pgl_unique%processes(i,iproc)
             stop 1
          endif
       enddo
       if (nq.eq.0 .and. nsing.ne.0) then
          ig=nsing+1
          is=1
          do i=1,next
             if (is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             else
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             endif
          enddo
!!$          write (*,*) 'when there are colour singlets, there should be quarks'
!!$          stop 1
       elseif (nq.eq.0) then
          do i=1,next
             pgl_unique%color_orders(i,iproc)=i
          enddo
       elseif (nq.eq.2) then
          ig=2
          is=ng+2
          do i=1,next
             if (is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(1,iproc)=i
             elseif (is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(next,iproc)=i
             elseif (is_gluon(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             elseif (is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             endif
          enddo
       elseif (nq.eq.4) then
          iq=1
          iaq=2
          ig=4
          is=ng+4
          do i=1,next
             if (is_quark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iq,iproc)=i
                iq=iq+2
             elseif (is_antiquark(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(iaq,iproc)=i
                iaq=next
             elseif (is_gluon(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(ig,iproc)=i
                ig=ig+1
             elseif (is_singlet(pgl_unique%processes(i,iproc))) then
                pgl_unique%color_orders(is,iproc)=i
                is=is+1
             endif
          enddo
       endif
    enddo
  end subroutine setup_color_order

  
  subroutine create_run_tag()
    implicit none
    integer :: i1,i2,i
    tag='_'       ! tag of current run
    tag_read='_'  ! same as 'tag', but with previous integration_step (i.e., defines the file to read the integration grids from)
!    call add_to_string(tag,PS_choice,.true.)
!    call add_to_string(tag_read,PS_choice,.true.)
    call add_to_string(tag,next,.true.)
    call add_to_string(tag_read,next,.true.)
    call add_to_string(tag,integration_step,.true.)
    if(integration_step.gt.0) then
       call add_to_string(tag_read,integration_step-1,.true.)
    else
       call add_to_string(tag_read,integration_step,.true.)
    endif
    if (allocated(part)) then
       do i=1,next
          call add_to_string(tag,part(i),.true.)
          call add_to_string(tag_read,part(i),.true.)
       enddo
       do i=1,next-1
          call add_to_string(tag,o(i),.true.)
          call add_to_string(tag_read,o(i),.true.)
       enddo
       call add_to_string(tag,o(next),.false.)
       call add_to_string(tag_read,o(next),.false.)
    else
!!$        ! just to the first process; the should all give the same value for 'i'
!!$       i1=ifindloc(color_orders(1:next,1),next,1)
!!$       i2=ifindloc(color_orders(1:next,1),next,2)
!!$       i=i2-i1 -1
!!$       if (i.lt.0) i=i+next
!!$       call add_to_string(tag,i,.false.)
!!$       call add_to_string(tag_read,i,.false.)
    endif
    write (*,*) 'File tag is: ',tag,'   ',add_arg
  end subroutine create_run_tag

  subroutine add_to_string(string,inter,add_underscore)
    ! Adds an integer 'inter' to the end of the string 'string' (followed by
    ! an underscore if 'add_underscore=.true.')
    implicit none
    character(len=string_len) :: string
    integer :: inter
    logical :: add_underscore
    character(len=1) :: s1
    character(len=2) :: s2
    character(len=3) :: s3
    if (inter.ge.0 .and. inter.le.9) then
       write(s1,'(i1)') inter
       string=trim(adjustl(string))//trim(adjustl(s1))
       if (add_underscore) string=trim(adjustl(string))//'_'
    elseif(inter.ge.-9 .and. inter.le.99) then
       write(s2,'(i2)') inter
       string=trim(adjustl(string))//trim(adjustl(s2))
       if (add_underscore) string=trim(adjustl(string))//'_'
    elseif(inter.ge.-99 .and. inter.le.999) then
       write(s3,'(i3)') inter
       string=trim(adjustl(string))//trim(adjustl(s3))
       if (add_underscore) string=trim(adjustl(string))//'_'
    else
       write (*,*) 'value too large to add to the run tag',inter
    endif
  end subroutine add_to_string

  subroutine set_initial_state_average_factor(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iproc
    do iproc=1,pgl%nproc
       do i=1,2
          if (pgl%processes(i,iproc).eq.21) then
             ! gluon: two polarisations and 8 colours
             pgl%iden(iproc)=pgl%iden(iproc)*2*8
          elseif (abs(pgl%processes(i,iproc)).ge.1 .and. abs(pgl%processes(i,iproc)).le.6) then
             ! (anti-)quark: two helicities and 3 colours
             pgl%iden(iproc)=pgl%iden(iproc)*2*3
          else
             ! assume two spin states and no colour:
             pgl%iden(iproc)=pgl%iden(iproc)*2
          endif
       enddo
    enddo
  end subroutine set_initial_state_average_factor
  
  subroutine set_final_state_identical_particle_factor(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,j,ni,iproc
    integer,dimension(:,:),allocatable :: iden_part
    allocate(iden_part(1:next,2))
    do iproc=1,pgl%nproc
       ni=0
       do i=3,next
          do j=1,ni
             if (iden_part(j,1).eq.pgl%processes(i,iproc)) then
                iden_part(j,2)=iden_part(j,2)+1
                exit
             endif
          enddo
          if (j.eq.ni+1) then
             ni=ni+1
             iden_part(j,1)=pgl%processes(i,iproc)
             iden_part(j,2)=1
          endif
       enddo
       do i=1,ni
          pgl%iden(iproc)=pgl%iden(iproc)*factorial8(iden_part(i,2))
       enddo
    enddo
    deallocate(iden_part)
  end subroutine set_final_state_identical_particle_factor

  subroutine compute_LC_colour_factor(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,ifac,iproc
    real(kind=8) :: fac
    integer :: it,lim

    lim=pgl%nproc
!    if (.not.read_proc_from_file.and.pgl%amps%same_flav(3)) lim=pgl%amps%nprocs
    do iproc=1,lim
       fac=0d0
       do i=1,next
          if (pgl%processes(i,iproc).eq.21) then
             fac=fac+1d0
          elseif (abs(pgl%processes(i,iproc)).ge.1 .and. abs(pgl%processes(i,iproc)).le.6) then
             fac=fac+0.5d0
          endif
       enddo
       ifac=nint(fac)
       if (dble(ifac).ne.fac) then
          write (*,*) 'There is some issue with the LC colour factor computation: '// &
               'colour factor is not an integer',ifac,fac
          stop 1
       endif
       if (abs(pgl%processes(pgl%color_orders(1,iproc),iproc)).ne.abs(pgl%processes(pgl%color_orders(next,iproc),iproc)) &
            .and. .not.pgl%amps%same_flav(iproc)) then
          ifac=ifac-2
       endif
       pgl%col_fac(iproc)=3**ifac
    enddo
  end subroutine compute_LC_colour_factor
  
  subroutine set_mint_module_special_parameters()
    ! these parameters need to be set for the mint-module to work correctly,
    ! but are irrelevant for any LO process; except for 'nchans'!)
    implicit none
    fixed_order=.false.
    nlo_ps=.true.
    n_ord_virt=1
    nchans=1
    iconfig=1
    ichan=1
    ifold_energy=1
    ifold_yij=1
    ifold_phi=1
    ifold(1:ndimmax)=1
    iconfigs(1:maxchannels)=1
    min_virt_fraction_mint=1d0
    virt_fraction=1d0
    wgt_mult=1d0
    average_virtual(0:n_ave_virt,maxchannels)=0d0
    virt_wgt_mint(0:n_ave_virt)=0d0
    born_wgt_mint(0:n_ave_virt)=0d0
    virtual_fraction(1:maxchannels)=1d0
    ans(1:nintegrals,0:maxchannels)=0d0
    unc(1:nintegrals,0:maxchannels)=0d0
    only_virt=.false.
  end subroutine set_mint_module_special_parameters
  
  subroutine set_ipdgs_for_PDF(pgl)
    ! determines for which flavours the PDFs should be evolved
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: iflav,iproc
    pgl%ipdgs(-6:7,1:2)=.false.
    do iflav=-6,7
       do iproc=1,pgl%nproc
          if (iflav.eq.0) then    ! gluon
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.21)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.21)) pgl%ipdgs(iflav,2)=.true.
          elseif(iflav.eq.7) then ! photon
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.22)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.22)) pgl%ipdgs(iflav,2)=.true.
          else                    ! quarks and anti-quarks
             if (any(pgl%iden_processes(1,1:pgl%iden_iproc(iproc),iproc).eq.iflav)) pgl%ipdgs(iflav,1)=.true.
             if (any(pgl%iden_processes(2,1:pgl%iden_iproc(iproc),iproc).eq.iflav)) pgl%ipdgs(iflav,2)=.true.
          endif
       enddo
    enddo
  end subroutine set_ipdgs_for_PDF

  subroutine include_PDF_and_identical_procs(val,val_abs,pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    real(kind=8),intent(inout),dimension(*) :: val,val_abs
    integer :: iproc,ip
    real(kind=8) :: xmu_fac
    real(kind=8), dimension(-6:7,2) :: PDF
    if (include_pdf) then
       ! Include the PDFs
       xmu_fac=91.188d0 ! factorisation scale
       call PDF_eval(1,pgl%ipdgs(-6,1),pgl%phase_space%xbjrk(1),xmu_fac,PDF(-6,1))
       call PDF_eval(1,pgl%ipdgs(-6,2),pgl%phase_space%xbjrk(2),xmu_fac,PDF(-6,2))
    endif
    do iproc=1,pgl%nproc
       pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)=val(iproc)*pgl%multi_factor(1:pgl%iden_iproc(iproc),iproc)
       if (include_pdf) then
          do ip=1,pgl%iden_iproc(iproc)
             ! first incoming particle
             if (pgl%iden_processes(1,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(0,1)
             elseif(pgl%iden_processes(1,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(7,1)
             else
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(pgl%iden_processes(1,ip,iproc),1)
             endif
             ! second incoming particle
             if (pgl%iden_processes(2,ip,iproc).eq.21) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(0,2)
             elseif(pgl%iden_processes(2,ip,iproc).eq.22) then
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(7,2)
             else
                pgl%val_procs(ip,iproc)=pgl%val_procs(ip,iproc)*PDF(pgl%iden_processes(2,ip,iproc),2)
             endif
          enddo
       endif
       val(iproc)=sum(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc))
       val_abs(iproc)=sum(abs(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)))
    enddo
  end subroutine include_PDF_and_identical_procs

  subroutine define_identical_procs(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: iproc,ip,n
    ! first fill the number of identical processes per iproc (so that we can
    ! allocate the array with the right size)
    allocate(pgl%iden_iproc(1:pgl%nproc))
    do iproc=1,pgl%nproc
       pgl%iden_iproc(iproc)=1
       if (any(abs(pgl%processes(1:next,iproc)).eq.1)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*5
       endif
       if (any(abs(pgl%processes(1:next,iproc)).eq.2)) then
          pgl%iden_iproc(iproc)=pgl%iden_iproc(iproc)*4
       endif
    enddo
    allocate(pgl%val_procs(1:maxval(pgl%iden_iproc(1:pgl%nproc)),1:pgl%nproc))
    allocate(pgl%iden_processes(1:next,1:maxval(pgl%iden_iproc(1:pgl%nproc)),1:pgl%nproc))
    ! Loop again and actually fill the iden_processes()
    do iproc=1,pgl%nproc
       do ip=0,pgl%iden_iproc(iproc)-1
          do n=1,next
             if (abs(pgl%processes(n,iproc)).eq.1) then
                pgl%iden_processes(n,ip+1,iproc)=sign(mod(ip,5)+1,pgl%processes(n,iproc))
             elseif (abs(pgl%processes(n,iproc)).eq.2) then
                if (mod(ip,5)+1.eq.ip/5+2) then
                   pgl%iden_processes(n,ip+1,iproc)=sign(1,pgl%processes(n,iproc))
                else
                   pgl%iden_processes(n,ip+1,iproc)=sign(ip/5+2,pgl%processes(n,iproc))
                endif
             else
                pgl%iden_processes(n,ip+1,iproc)=pgl%processes(n,iproc)
             endif
          enddo
       enddo
    enddo
  end subroutine define_identical_procs
  
  subroutine compute_multichannel_symmetry_factor(sym_fac)
    implicit none
    integer(kind=8),intent(out) :: sym_fac
    integer :: ngl=0,ngl_tot=0,n_sing=0
    integer,dimension(6) :: nq,naq
    integer :: i,j
    integer(kind=8) :: tot_ord
    integer :: ic,i_qq,iaq,i_ini,i_inv,i_swap,k,l,ip
    integer,dimension(next) :: fgluons,ips,ips_out
    logical :: same_flavour
    logical,dimension(next) :: fgluon
    integer,dimension(:,:),allocatable :: io_list,io

    nq=0
    naq=0
    ! count the number of final state gluons and quarks
    do i=3,next
       if (part(i).eq.21) then
          ngl=ngl+1
       endif
       do j=1,6
         if (part(i).eq.j) nq(j)=nq(j)+1
         if (part(i).eq.-j) naq(j)=naq(j)+1
       enddo
    enddo
    do i=1,next
       if (part(i).eq.21) then
          ngl_tot=ngl_tot+1
       endif
    enddo

    ! Since we only need to include a subset of all the colour-orderings, we
    ! need to compensate with a symmetry factor
    if (nquarks.eq.0) then
       tot_ord=factorial8(ngl_tot-1)
       ! All gluon process. This assumes that the only channels we are
       ! including are strictly different. We distinguish them by considering
       ! how many (final state) gluons are attached to the two colour lines
       ! that link the two incoming gluons. Hence, we only include
       ! floor(next/2) channels, e.g., for next=6 we only consider:
       ! i   --> 1,2,3,4,5,6   (0 and 4 gluons on the two lines)
       ! ii  --> 1,3,2,4,5,6   (1 and 3 gluons on the two lines)
       ! iii --> 1,3,4,2,5,6   (2 and 2 gluons on the two lines)
       ! And, e.g., for next=9, we only consider:
       ! i   --> 1,2,3,4,5,6,7,8,9   (0 and 7 gluons on the two lines)
       ! ii  --> 1,3,2,4,5,6,7,8,9   (1 and 6 gluons on the two lines)
       ! iii --> 1,3,4,2,5,6,7,8,9   (2 and 5 gluons on the two lines)
       ! iv  --> 1,3,4,5,2,6,7,8,9   (3 and 4 gluons on the two lines)
       ! This means that the sym_fac should be equal to the number of final
       ! state gluon permutations, multiplied by 2 (except if we have an equal
       ! number of gluons on both colour lines that attached the two incoming
       ! gluons).
       if (c_o*2.eq.(ngl)) then
          sym_fac=factorial8(ngl)
       else
          sym_fac=2*factorial8(ngl)
       endif
    elseif (nquarks.eq.2) then
       tot_ord=factorial8(ngl_tot)
       if ((abs(part(1)).ge.1 .and. abs(part(1)).le.6) .and. &
           (abs(part(2)).ge.1 .and. abs(part(2)).le.6) )then
          ! quark and anti-quark are incoming. Only 1 channel needed,
          ! which would result in the following symmetry factor:
          sym_fac=factorial8(ngl)
       elseif ((abs(part(1)).ge.1 .and. abs(part(1)).le.6) .or. &
               (abs(part(2)).ge.1 .and. abs(part(2)).le.6) )then
          ! one incoming quark (or anti-quark). There are ngluons
          ! channels needed: they correspond to having the incoming
          ! gluon at all possible positions between the quark and
          ! anti-quark in the colour order. Hence, each channel comes
          ! with an (ngluons-1)! symmetry factor:
          sym_fac=factorial8(ngl)
       else
          ! both quark and anti-quark are final state. This is similar
          ! to the all-gluon case above, treating the q-qbar pair as
          ! another gluon. This special gluon is identifiable! So, for
          ! next=6 (and assuming that the qqbar pair are particles 5
          ! and 6) one has the following possibilities:
          !
          ! ia   --> 1,2,3,4,(5,6) = 5,4,3,2,1,6 ---- : both gluons on the same
          ! ib   --> 1,2,3,(5,6),4 = 5,3,2,1,4,6 --/         line as the qqbar pair
          ! ic   --> 1,2,(5,6),3,4 = 5,2,1,4,3,6 -/
          ! iia  --> 1,3,2,4,(5,6) = 5,4,2,3,1,6 ---- : one gluon on the same 
          ! iib  --> 1,3,2,(5,6),4 = 5,2,3,1,4,6 -/          line as the qqbar pair
          ! iii  --> 1,3,4,2,(5,6) = 5,2,4,3,1,6 ---- : both gluons on the other quark line
          !
          ! Note, however, that in the above situation ia and ic results in
          ! the same cross section. Also iia and iib result in the same
          ! rate. Furthermore, there is an additional factor 2, since one can
          ! interchange the two incoming particles. Hence, there are only
          ! trully 4 independent colour orders to consider:
          ! ia, with symmetry factor (ngluon-2)!*2*2
          ! ib, with symmetry factor (ngluon-2)!*2
          ! iaa, with symmetry factor (ngluon-2)!*2*2
          ! iii, with symmetry factor (ngluon-2)!*2
          ! Hence
          sym_fac=factorial8(ngl)*2
          if (ngl.gt.0) then
               sym_fac=factorial8(ngl)*2
          endif
          if (ifindloc(o,next,1).ne.next-ifindloc(o,next,2)+1) then
             sym_fac=sym_fac*2
          endif
       endif
    elseif (nquarks.eq.4) then
       ! total number of potentially different orders: two ways of connecting
       ! quarks, (n-4)! orderings for the gluons, n-3 ways for an order to
       ! distribute the gluons among the two quark lines
       tot_ord=2*factorial8(ngl_tot)*(ngl_tot+1)
       n_sing=next-4-ngl_tot

       do i=2,next-1
          if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
              (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
             if (abs(part(o(i))).eq.abs(part(o(1))) .and. abs(part(o(i))).eq.abs(part(o(next)))) then
                same_flavour=.true.
             else
                same_flavour=.false.
             endif
             exit
          endif
       enddo
       
       allocate(io_list(1:next,tot_ord))
       allocate(io(1:next,5))
       ic=0
       ! 1. two ways of connecting quarks with anti-quarks
       do_i_qq: do i_qq=1,2
          io(1:next,1)=o(1:next)
          if (i_qq.eq.2) then
             if (.not. same_flavour) cycle
             do i=2,next-1
                if ((o(i).gt.2 .and. part(o(i)).le.-1 .and. part(o(i)).ge.-6) .or. &
                    (o(i).le.2 .and. part(o(i)).ge. 1 .and. part(o(i)).le. 6) ) then
                   ! check if it would give the same result:
                   if ( ((o(1   ).le.2 .and. o(i+1).gt.2) .or. (o(1   ).gt.2 .and. o(i+1).le.2)) .and. &
                        ((o(next).le.2 .and. o(i  ).gt.2) .or. (o(next).gt.2 .and. o(i  ).le.2)) ) then
                      ! not both initial or both final state
                      cycle do_i_qq
                   endif
                   iaq=o(i)
                   io(i,1)=o(next)
                   io(next,1)=iaq
                   exit
                endif
             enddo
          endif

          ! 2. invert order of two initial states
          do i_ini=1,2
             io(1:next,2)=io(1:next,1)
             if (part(1).ne.part(2) .and. i_ini.eq.2) cycle
             if (i_ini.eq.2) then
                do i=1,next
                   if (io(i,1).eq.1) then
                      io(i,2)=2
                   elseif (io(i,1).eq.2) then
                      io(i,2)=1
                   endif
                enddo
             endif
             ! 3. invert order of both quark lines
             do_i_inv: do i_inv=1,2
                io(1:next,3)=io(1:next,2)
                if (i_inv.eq.2) then
                   ! invert order
                   do i=2,next-n_sing-1
                      if ((io(i,2).gt.2 .and. part(io(i,2)).le.-1 .and. part(io(i,2)).ge.-6) .or. &
                           (io(i,2).le.2 .and. part(io(i,2)).ge. 1 .and. part(io(i,2)).le. 6) ) then
                         if ( ((io(1   ,2).le.2 .and. io(i  ,2).gt.2) .or. (io(1   ,2).gt.2 .and. io(i  ,2).le.2)) .or. &
                              ((io(next,2).le.2 .and. io(i+1,2).gt.2) .or. (io(next,2).gt.2 .and. io(i+1,2).le.2)) ) then
                            ! not both initial or both final state. Cannot invert order.
                            cycle do_i_inv
                         endif
                         io(2:i-1,3)=io(i-1:2:-1,2)
                         io(i+2:next-n_sing-1,3)=io(next-n_sing-1:i+2:-1,2)
                         exit
                      endif
                   enddo
                endif
                ! 4. If both quark lines are similar (identical quarks and
                ! FF+FF or IF+FI or FI+IF or FI+IF or IF+FI), we can swap the gluons from one line to the
                ! other
                do i_swap=1,2
                   io(1:next,4)=io(1:next,3)
                   if (i_swap.eq.2) then
                      if (.not.same_flavour) cycle
                      do i=2,next-1
                         if ((io(i,3).gt.2 .and. part(io(i,3)).le.-1 .and. part(io(i,3)).ge.-6) .or. &
                              (io(i,3).le.2 .and. part(io(i,3)).ge. 1 .and. part(io(i,3)).le. 6) ) then
                            if ((io(1,3).gt.2 .and. io(i,3).gt.2 .and. io(i+1,3).gt.2 .and. io(next,3).gt.2) .or. &
                                (io(1,3).gt.2 .and. io(i,3).le.2 .and. io(i+1,3).le.2 .and. io(next,3).gt.2) .or. &
                                (io(1,3).le.2 .and. io(i,3).gt.2 .and. io(i+1,3).gt.2 .and. io(next,3).le.2) .or. &
                                (io(1,3).gt.2 .and. io(i,3).le.2 .and. io(i+1,3).gt.2 .and. io(next,3).le.2) .or. &
                                (io(1,3).le.2 .and. io(i,3).gt.2 .and. io(i+1,3).le.2 .and. io(next,3).gt.2) ) then
                               io(2:next-i-1,4)=io(i+2:next-1,3) ! gluons
                               io(next-i:next-i+1,4)=io(i:i+1,3) ! qbarq
                               io(2+next-i:next-1,4)=io(2:i-1,3) ! gluons
                               exit
                            endif
                         endif
                      enddo
                      if (i.eq.next) cycle
                   endif
                   ! 5. permute all final state gluons
                   k=0
                   l=0
                   do i=1,next
                      if (part(io(i,4)).eq.21 .and. io(i,4).gt.2) then
                         k=k+1
                         fgluons(k)=io(i,4)
                         fgluon(i)=.true.
                      else
                         fgluon(i)=.false.
                      endif
                   enddo
                   io(1:next,5)=io(1:next,4)
                   do ip=1,factorial(k)
                      if (ip.eq.1) then
                         do i=1,k
                            ips(i)=i
                         enddo
                      else
                         call get_next_iperm(k,ips,ips_out,k)
                         ips(1:k)=ips_out(1:k)
                      endif
                      i=0
                      l=1
                      do
                         i=i+1
                         if (i.eq.next) exit
                         if (fgluon(i)) then
                            j=0
                            do
                               if (fgluon(i+j+1)) then
                                  j=j+1
                               else
                                  exit
                               endif
                            enddo
                            io(i:i+j,5)=fgluons(ips(l:l+j))
                            l=l+j+1
                            i=i+j
                         endif
                      enddo
                      ! ---> if not yet in list of identical contributions, add it!
                      do i=1,ic
                         if (all(io(1:next,5).eq.io_list(1:next,i))) exit
                      enddo
                      if (i.eq.ic+1) then
                         ! new identical contribution
                         io_list(1:next,i)=io(1:next,5)
                         ic=i
                      endif
                   enddo
                enddo
             enddo do_i_inv
          enddo
       enddo do_i_qq
       sym_fac=ic
    else        
       write (*,*) 'WARNING: symmetry factor missing',nquarks
    endif
    write (*,*) 'total number of orders is',tot_ord,' and multi-channel symmetry factor is',sym_fac,&
         '. They should be the same when including all channels.'
  end subroutine compute_multichannel_symmetry_factor

  integer(kind=4) function ifindloc(a,n,i)
    ! returns the location of the value 'i' in the array 'a' (of size 'n'). If
    ! the array 'a' does not contain 'i', the value 'n+1' is returned.
    implicit none
    integer,intent(in) :: i,n
    integer,dimension(n),intent(in) :: a
    do ifindloc=1,n
       if (a(ifindloc).eq.i) return
    enddo
  end function ifindloc

end program matrix_integrate_QCD
