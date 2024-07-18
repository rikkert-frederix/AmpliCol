
program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_gen23
  use phase_space_genpt
  use haag
  use math_functions
  implicit none
  type(amplitude_QCD) :: amps
  type(amplitude_QCD) :: amps_sf
  integer :: next
  real(kind=8) :: amp2,weight
  real(kind=8),dimension(:),allocatable :: amp2_hel
  integer :: j,c_o,i
  integer(kind=8) :: sym_fac,iden
  integer(kind=4),dimension(:),allocatable :: o,part,orig_part,part_sf,hel,hel_fac
  integer(kind=4),dimension(:,:),allocatable :: spin
  real(kind=8),dimension(:),allocatable :: mass,width
  real(kind=8) :: s_cut(2),sqrts
  logical :: t_chan
  character(len=80) :: filename
  integer(kind=4) :: integration, nquarks
  logical,dimension(-6:7,2) :: ipdgs
  integer :: col_fac,nhel
  integer :: it ! quark order type
  integer,parameter :: nevent_hel_filter=5
  integer,dimension(2) :: hel_picked
  
  call get_run_arguments()
  call compute_multichannel_symmetry_factor()
  call create_run_tag()

  allocate(mass(next))
  allocate(width(next))

  call cpu_time(tTot_B)

! relevant input parameters for integration
  ! Number of events to generate. (If negative, start
  ! from a small number of points and double it each
  ! iteration. If positive, this is the number of
  ! points per iteration as well).
  if (imode.eq.0 .or. imode.eq.2) then
     ncalls0=-100
  else
     ncalls0=640000
  endif

  ndim=3*(next-2)-4   ! Number of dimensions of the integration.

  itmax=10         ! Number of iterations. (If ncalls0 < 0, the
                   ! integration is aborted if accuracy (next line)
                   ! has been reached.

  accuracy=0.00001d0 ! Accuracy of the integration. (Ignored if ncalls0 > 0).


! relevant physics input parameters and initialisation of amplitudes

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

  mass(1:next)=0d0
  width(1:next)=0d0
!  mass(1:2) = 0d0
!  mass(3:4) = 173d0
!  mass(5) = 0d0
!  width(1:2) = 0d0
!  width(3:4) = 1.491500d0
!  width(5) = 0d0

  call cpu_time(tBefore)
  t_chan=.false.
  if (integration.eq.1) then
     call gen23_init(sqrts,next,mass,o,part,s_cut,pt_min,DRjj_min,.false.,include_pdf)
  elseif  (integration.eq.2) then
     call haag_init(sqrts,next,mass,o,part,s_cut,t_chan,include_pdf)
  elseif (integration.eq.3) then
     call genpt_init(sqrts,next,mass,pt_min,eta_max,DRjj_min,include_pdf)
  elseif (integration.eq.4) then
     call gen23_init(sqrts,next,mass,o,part,s_cut,pt_min,DRjj_min,.true.,include_pdf)
  endif
  call cpu_time(tAfter)
  t_PS_init=t_PS_init+tAfter-tBefore

  iden=1
  call set_final_state_identical_particle_factor(iden)
  call set_initial_state_average_factor(iden)

  ! initialize the amplitudes (sets up the imaps(), helicity maps,
  ! colour factors, etc.)
  call cpu_time(tBefore)
  it = 0 ! dummy
  orig_part(:)=part(:)

  if (include_pdf) then
     ndim=ndim+2
     call PDF_initialise
     call set_ipdgs_for_PDF(ipdgs)
  endif

  ! counting of quark flavours in process
  call fill_quark_info()

  if (amps%n_qqbar.eq.2) then
    call define_symm_2qq(next,part,1)
  endif
  call amps%init(1,next,part,spin,mass,width,o,it)

  if (amps%n_qqbar.eq.2.and.amps%same_flav) then
     part_sf(:) = orig_part(:)
     amps_sf%n_qqbar=amps%n_qqbar
     amps_sf%same_flav=amps%same_flav
     call define_symm_2qq(next,part_sf,2)
     call amps_sf%init(1,next,part_sf,spin,mass,width,o,it)
  endif

  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  ! Compute the leading colour factor
  if (amps%n_qqbar.eq.2) then
      if (abs(part(o(1))).ne.abs(part(o(next)))) it = 2
  endif

  call compute_LC_colour_factor(col_fac,it)
  
  ! number of helicities to sum over
  nhel=(amps%n_cur_end(next-1)-amps%n_cur_start(next-1)+1)*(amps%n_cur_end(next)-amps%n_cur_start(next)+1)
  allocate(amp2_hel(1:nhel))
  allocate(hel(1:next))
  allocate(hel_fac(1:nhel))
  hel_fac(1:nhel)=1
  

  ! Not so relevant mint-module parameters: only used in special cases.
  call set_mint_module_special_parameters()

  if (imode.le.1) then
     ! grid setup, or computation of upper bounding envelope
     call mint(integrand)
  else
     ! actual (unweighted) event generation
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise counters
     filename='Outputs'//trim(adjustl(add_arg))//'/events'//trim(adjustl(tag))//'.lhe'
     open(unit=11,file=filename,status='unknown')
     do j=1,abs(ncalls0)
        call gen(integrand,1,2) ! generate an unweighted event
        call unwgt_helicity     ! pick a random helicity
        call write_event(11,ans(1,0))
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
  write(*,*) 'Number of events:',all_evt
  write(*,*) 'Number passing cuts:',passed
  write(*,*) 'Fraction passing:',float(passed)/float(all_evt)
  write(*,*) 'Number of numerical errors:',num_error
 
contains
  real(kind=8) function integrand(x,vol,ifirst,f1)
    implicit none
    integer :: ifirst
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(nintegrals) :: f1
    real(kind=8), save :: val
    integer :: ih
    real(kind=8) :: vol,cuts_wgt
    real(kind=8), parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real*4 :: tBefore,tAfter
    real(kind=8) :: Q
   
    ! some point-by-point initialisation
    f1(1:nintegrals)=0d0
    if (ifirst.eq.2) then
       ! use previously computed integrand
       f1(1)=abs(val)
       f1(2)=val
       return
    endif
    new_point=.true.
    pass_cuts_check=.true.

    ! Generate phase-space point based on the random numbers 'x(1:ndim)'
    call cpu_time(tBefore)
    if (integration.eq.1 .or. integration.eq.4)then
       call gen23_phase_space(x)
    elseif (integration.eq.2) then
        call PS_haag(x)
    elseif (integration.eq.3) then
        call genpt_phase_space(x)
    endif
    call cpu_time(tAfter)
    t_PS= t_PS +tAfter-tBefore


    if (debug ) then
       write (*,*) jac
       stop 1
    endif
    
    all_evt=all_evt+1

    cuts_wgt=pass_cuts(next,p)
    if ((jac.lt.0d0) .or. (smooth_cuts .and. cuts_wgt.lt.0d0) .or. (.not.smooth_cuts .and. cuts_wgt.lt.1d0)) then
       pass_cuts_check=.false.
       val=0d0
       return
    endif

    passed = passed + 1

    ! compute amplitudes
    call cpu_time(tBefore)

    call amps%evaluate(next,p,mass,width,hel,part)

    if (amps%n_qqbar.eq.2.and.amps%same_flav) then
      call amps_sf%evaluate(next,p,mass,width,hel,part_sf)
      do ih=1,nhel
        if (it.eq.2) then
           amps%amps(ih)=(1d0/3d0)*amps%amps(ih)+amps_sf%amps(ih)
        else
           amps%amps(ih)=amps%amps(ih)+(1d0/3d0)*amps_sf%amps(ih)
        endif
      enddo
    endif
    call cpu_time(tAfter)
    t_amp=t_amp+tAfter-tBefore

    call cpu_time(tBefore)
    amp2_hel(1:nhel)=0d0
    do ih=1,nhel
       if (use_real_gluons .and. amps%n_qqbar.eq.0) then
          amp2_hel(ih)=amp2_hel(ih)+amps%amps_r(ih)*col_fac*amps%amps_r(ih)
       else
          amp2_hel(ih)=amp2_hel(ih)+dble(amps%amps(ih)*col_fac*dconjg(amps%amps(ih)))
       endif
       amp2_hel(ih)=amp2_hel(ih)*hel_fac(ih)
    enddo

    amp2=sum(amp2_hel(1:nhel))
    
    if (passed.le.nevent_hel_filter) then
       call setup_helicity_filter(passed)
       if (imode.eq.2 .and. passed.eq.nevent_hel_filter) then
          ! since we update the helicities we need to compute when
          ! passed==nevent_hel_filter, the unweighting of the helicities goes
          ! wrong for this phase-space point. Hence, we need to skip it.
          amp2=0d0
       endif
    endif
    
    weight=vol*jac*(4*pi*alphas)**(next-2-amps%n_sing)/dble(iden)*conv
    
    if (amps%n_sing.ge.1) then
       do i=1,next
          if (abs(part(i)).le.6) then
             if (mod(abs(part(i)),2).eq.0) Q=2d0/3d0
             if (mod(abs(part(i)),2).eq.1) Q=-1d0/3d0
          endif
       enddo
       weight=weight*(Q**2*2d0*4d0*pi*alphaEW)**amps%n_sing
    endif

    val=amp2*weight

    ! Apply the weight from the cuts
    if (smooth_cuts) val=val*cuts_wgt

    ! Since we only need to include a subset of all the colour-orderings, we
    ! need to compensate with a symmetry factor
    val=val*sym_fac

    if (include_PDF) then
       call multiply_by_PDF_value(val)
    endif
    ! pass the result to the mint module
    f1(1)=abs(val)
    f1(2)=val
    integrand=f1(2)
    
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
       if (abs(part(i)).ge.0.and.abs(part(i)).le.6) then ! for quarks
         frac  = 0.9d0
         steep = 0.1d0
       elseif (part(i).eq.21.or.part(i).eq.22) then ! for gluons and photons
         frac  = 0.8d0
         steep = 0.1d0
       endif
       if (pt_min.gt.0d0) then
          if (pt(p(0,i)).lt.frac*pt_min) then
             pass_cuts=-1d0
             return
          endif
          if (pt(p(0,i)).gt.frac*pt_min.and.pt(p(0,i)).lt.pt_min) then
             y=(pt(p(0,i))-frac*pt_min)/(pt_min*(1d0-frac))
             if (imode.le.0) then
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

  subroutine setup_helicity_filter(nevent)
    implicit none
    real(kind=8) :: max_value
    integer :: ih1,ih2,nevent,nhel_sf
    integer,dimension(:,:),allocatable,save :: include_hel
    integer,dimension(:),allocatable :: include_hel_sf
    if (.not.allocated(include_hel)) allocate(include_hel(nhel,nevent_hel_filter))
    ! filter zero helicities and helicities that are identical
    include_hel(1:nhel,nevent)=1
    max_value=maxval(amp2_hel(1:nhel))
    do ih1=1,nhel
       if (include_hel(ih1,nevent).ne.1) cycle
       if (amp2_hel(ih1)/max_value.lt.1d-28) then
          ! zero
          include_hel(ih1,nevent)=0
       else
          do ih2=ih1+1,nhel
             if (abs(amp2_hel(ih1)-amp2_hel(ih2))/abs(amp2_hel(ih1)+amp2_hel(ih2)).lt.1d-10) then
                ! identical
                include_hel(ih2,nevent)=-ih1
                include_hel(ih1,nevent)=include_hel(ih1,nevent)+1
             endif
          enddo
       endif
    enddo

    if (nevent.lt.nevent_hel_filter) return
    
    do ih1=1,nhel
       if (any(include_hel(ih1,2:nevent_hel_filter).ne.include_hel(ih1,1))) then
          write (*,*) 'inconsistent helicity. Cannot setup helicity filter.'
          write (*,*) ih1,nevent_hel_filter,':',include_hel(ih1,1:nevent_hel_filter)
          stop 1
       endif
    enddo
    nhel_sf=nhel
    allocate(include_hel_sf(1:nhel))
    include_hel_sf(1:nhel)=include_hel(1:nhel,1)
    call amps%filter_helicity(next,nhel,include_hel(1,1)) ! this updates 'nhel' and 'include_hel'
    if (amps%n_qqbar.eq.2.and.amps%same_flav) then
       call amps_sf%filter_helicity(next,nhel_sf,include_hel_sf)
       if (nhel.ne.nhel_sf) then
          write (*,*) 'number of helicity not consistent',nhel,nhel_sf
          stop 1
       endif
    endif
    deallocate(hel_fac)
    allocate(hel_fac(nhel))
    hel_fac(1:nhel)=include_hel(1:nhel,1)
    deallocate(include_hel)
    deallocate(include_hel_sf)
  end subroutine setup_helicity_filter

  subroutine define_symm_2qq(next,part,chan)
    implicit none
    integer :: next,chan
    integer, dimension(next) :: part
    integer :: i,j,sgn
    logical :: first
    if (amps%same_flav) then
       if (chan.eq.2) then
          do i=1,next
             if (abs(part(i)).gt.0.and.abs(part(i)).lt.6) then
                first=.true.
                do j=i+1,next
                   if (i.le.2.and.j.le.2) sgn=-1
                   if (i.le.2.and.j.gt.2) sgn=+1
                   if (i.gt.2.and.j.gt.2) sgn=-1
                   if (part(j).eq.sgn*part(i).and..not.first) then
                      part(i) = sign(abs(orig_part(i))+1,orig_part(i))
                      part(j) = sgn*(part(i))
                      exit
                   endif
                   if (part(j).eq.sgn*part(i).and.first) then
                      first = .false.
                   endif
                enddo
             endif
          enddo
       elseif (chan.eq.1) then
          do i=1,next
             if (abs(orig_part(i)).gt.0.and.abs(orig_part(i)).lt.6) then
                do j=i+1,next
                   if (i.le.2.and.j.le.2) sgn=-1
                   if (i.le.2.and.j.gt.2) sgn=+1
                   if (i.gt.2.and.j.gt.2) sgn=-1
                   if (orig_part(j).eq.sgn*orig_part(i)) then
                      part(i) = sign(abs(orig_part(i))+1,orig_part(i))
                      part(j) = sgn*(part(i))
                      exit
                   endif
                enddo
                exit
             endif
          enddo
       endif
    endif
  end subroutine define_symm_2qq

  subroutine fill_quark_info()
    implicit none
    integer,dimension(8) :: flav
    integer :: k
    flav = 0
    k = 1
    amps%n_qqbar= 0
    amps%same_flav=.true.
    do i=1,next
       if (i.le.2) then
          if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
             flav(k) = abs(orig_part(i))
             k= k+1
             if (orig_part(i).lt.0) amps%n_qqbar=amps%n_qqbar+1
          endif
       else
          if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
             flav(k) = abs(orig_part(i))
             k= k+1
             if (orig_part(i).gt.0) amps%n_qqbar=amps%n_qqbar+1
          endif
       endif
    enddo
    if (any(flav(1:2*amps%n_qqbar).ne.flav(1))) amps%same_flav = .false.
  end subroutine fill_quark_info
  
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
    write (iunit,*) next,wgt,amp2*weight,amp2,weight
    write (iunit,'(100i3)') amps%spins(1:next,hel_picked(1),hel_picked(2))
    write (iunit,'(100i3)') o(1:next)
    ! Since some of the symmetry factors (in particular for gg->qqbar+ng)
    ! compensate for reducing the number of integration channels (see
    ! sym_fac()) assuming symmetric initial states, we need to randomly flip
    ! all z-components in those cases. Easiest to always do this if the two
    ! incoming particles are identical.
    if (orig_part(1).ne.orig_part(2) .or. ran2().lt.0.5d0) then
       ! do not flip
       do i=1,next
          write (iunit,*) orig_part(i),p(1:3,i),p(0,i)
       enddo
    else
       ! do flip
       do i=1,next
          if (i.le.2) then
             write (iunit,*) orig_part(i),p(1:2,3-i),-p(3,3-i),p(0,3-i)
          else
             write (iunit,*) orig_part(i),p(1:2,i),-p(3,i),p(0,i)
          endif
       enddo
    endif
    write (iunit,*) '</event>'
  end subroutine write_event

  subroutine unwgt_helicity
    implicit none
    integer :: i
    real(kind=8) :: random
    real(kind=8),external :: ran2
    random=ran2()*amp2
    i=1
    do
       if (amp2_hel(i).gt.random) then
          exit
       else
          i=i+1
          amp2_hel(i)=amp2_hel(i)+amp2_hel(i-1)
       endif
    enddo
    hel_picked(2)=i
    if (hel_picked(2).gt.nhel) then
       write (*,*) 'Could not unweight helicity',hel_picked,nhel
       stop 1
    endif
    if (hel_fac(hel_picked(2)).gt.1) then
       hel_picked(1)=1+int(ran2()*hel_fac(hel_picked(2)))
    endif
  end subroutine unwgt_helicity
  
  subroutine get_run_arguments()
    implicit none
    integer :: argc
    integer :: i,k
    character(len=256) :: argv
    logical :: found_1
    integer(kind=8) iseed
    common /to_seed/iseed
    iseed=0
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    if (argc.le.10) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'integration, mode, next, *process*, *order*'
       stop 2
    else
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,*) integration
          if (i.eq.2) read(argv,*) imode
          if (i.eq.3) then
             read(argv,*) next
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             allocate(part(1:next))
             allocate(orig_part(1:next))
             allocate(part_sf(1:next))
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
    endif

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

    call setup_spin()
    
    ! basic checks:
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (imode.ne.0 .and. imode.ne.1 .and. imode.ne.2) then
       write (*,*) 'Incorrect imode',imode
       stop
    endif

    if (integration.ne.1 .and. integration.ne.2 .and. integration.ne.3 .and. integration.ne.4) then
       write (*,*) 'Integration modes only 1, 2, 3 or 4',integration
       stop
    endif
    if ((nquarks.ne.0 .and. nquarks.ne.2 .and. nquarks.ne.4) .or. (nquarks.gt.next)) then
       write (*,*) 'Not consistent number of external quarks (up to 2)',nquarks
       stop
    endif
  end subroutine get_run_arguments

  subroutine setup_spin()
    implicit none
    if (.not. allocated(spin)) allocate(spin(0:3,1:next))
    do i=1,next
       if (abs(part(i)).le.6 .or. part(i).eq.21 .or. part(i).eq.22) then
          spin(0,i)=2   ! two spin states: '-1' and '1'
          spin(1,i)=-1
          spin(2,i)=1
       else
          write (*,*) 'spin state not known',i,part(i)
          stop 1
       endif
    enddo
  end subroutine setup_spin
  
  subroutine create_run_tag()
    implicit none
    tag='_'       ! tag of current run
    tag_read='_'  ! same as 'tag', but with previous imode (i.e., defines the file to read the integration grids from)
!    call add_to_string(tag,integration,.true.)
!    call add_to_string(tag_read,integration,.true.)
    call add_to_string(tag,next,.true.)
    call add_to_string(tag_read,next,.true.)
    call add_to_string(tag,imode,.true.)
    if(imode.gt.0) then
       call add_to_string(tag_read,imode-1,.true.)
    else
       call add_to_string(tag_read,imode,.true.)
    endif
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

  subroutine set_initial_state_average_factor(iden)
    implicit none
    integer(kind=8),intent(inout) :: iden
    integer :: i
    do i=1,2
       if (part(i).eq.21) then
          ! gluon: two polarisations and 8 colours
          iden=iden*2*8
       elseif (abs(part(i)).ge.1 .and. abs(part(i)).le.6) then
          ! (anti-)quark: two helicities and 3 colours
          iden=iden*2*3
       else
          ! assume two helicities:
          iden=iden*2
       endif
    enddo
  end subroutine set_initial_state_average_factor
  
  subroutine set_final_state_identical_particle_factor(iden)
    implicit none
    integer(kind=8),intent(inout) :: iden
    integer :: i,j,ni=0
    integer,dimension(:,:),allocatable :: iden_part
    allocate(iden_part(1:next,2))
    do i=3,next
       do j=1,ni
          if (iden_part(j,1).eq.part(i)) then
             iden_part(j,2)=iden_part(j,2)+1
             exit
          endif
       enddo
       if (j.eq.ni+1) then
          ni=ni+1
          iden_part(j,1)=part(i)
          iden_part(j,2)=1
       endif
    enddo
    do i=1,ni
       iden=iden*factorial8(iden_part(i,2))
    enddo
    deallocate(iden_part)
  end subroutine set_final_state_identical_particle_factor

  subroutine compute_LC_colour_factor(col_fac,it)
    implicit none
    integer,intent(inout) :: col_fac
    integer :: i,ifac
    real(kind=8) :: fac=0d0
    integer :: it
    do i=1,next
       if (part(i).eq.21) then
          fac=fac+1d0
       elseif (abs(part(i)).ge.1 .and. abs(part(i)).le.6) then
          fac=fac+0.5d0
       endif
    enddo
    ifac=nint(fac)
    if (dble(ifac).ne.fac) then
       write (*,*) 'There is some issue with the LC colour factor computation: '// &
            'colour factor is not an integer',ifac,fac
       stop 1
    endif
    if (it.eq.2.and..not.amps%same_flav) then
        ifac=(ifac-2) 
    endif
    col_fac=3**ifac
  end subroutine compute_LC_colour_factor
  
  subroutine set_mint_module_special_parameters()
    ! these parameters need to be set for the mint-module to work correctly,
    ! but are irrelevant for any LO process
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
  
  subroutine set_ipdgs_for_PDF(ipdgs)
    implicit none
    logical,dimension(-6:7,2) :: ipdgs
    ipdgs(-6:7,1:2)=.false.
    if (orig_part(1).eq.21) then
       ipdgs(0,1)=.true.    ! gluon is '0'
    elseif (orig_part(1).eq.22) then
       ipdgs(7,1)=.true.    ! photon is '7'
    elseif (abs(orig_part(1)).ge.1 .and. abs(orig_part(1)).le.6) then
       ipdgs(orig_part(1),1)=.true.
    else
       write (*,*) 'unknown PDF 1',orig_part(1)
       stop 1
    endif
    if (orig_part(2).eq.21) then
       ipdgs(0,2)=.true.    ! gluon is '0'
    elseif (orig_part(2).eq.22) then
       ipdgs(7,2)=.true.    ! photon is '7'
    elseif (abs(orig_part(2)).ge.1 .and. abs(orig_part(2)).le.6) then
       ipdgs(orig_part(2),2)=.true.
    else
       write (*,*) 'unknown PDF 2',orig_part(2)
       stop 1
    endif
  end subroutine set_ipdgs_for_PDF
  
  subroutine multiply_by_PDF_value(val)
    implicit none
    real(kind=8),intent(inout) :: val
    real(kind=8) :: xmu_fac
    real(kind=8), dimension(-6:7,2) :: PDF
    ! Include the PDFs
    xmu_fac=91.188d0 ! factorisation scale

    call PDF_eval(1,ipdgs(-6,1),xbjrk(1),xmu_fac,PDF(-6,1))
    call PDF_eval(1,ipdgs(-6,2),xbjrk(2),xmu_fac,PDF(-6,2))
    if (orig_part(1).eq.21) then
       val=val*PDF(0,1)
    elseif (orig_part(1).eq.22) then
       val=val*PDF(7,1)
    else
       !write(*,*) 'orig pdf',orig_part(1)
       val=val*PDF(orig_part(1),1)
    endif
    if (orig_part(2).eq.21) then
       val=val*PDF(0,2)
    elseif (orig_part(2).eq.22) then
       val=val*PDF(7,2)
    else
       val=val*PDF(orig_part(2),2)
    endif
  end subroutine multiply_by_PDF_value
  
  subroutine compute_multichannel_symmetry_factor()
    implicit none
    integer :: ngl=0
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

    ! Since we only need to include a subset of all the colour-orderings, we
    ! need to compensate with a symmetry factor
    if (nquarks.eq.0) then
       tot_ord=factorial8(next-1)
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
       tot_ord=factorial8(next-2)
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
          if (ifindloc(o,next,1).ne.next-ifindloc(o,next,2)+1) then
             sym_fac=sym_fac*2
          endif
       endif
    elseif (nquarks.eq.4) then
       ! total number of potentially different orders: two ways of connecting
       ! quarks, (n-4)! orderings for the gluons, n-3 ways for an order to
       ! distribute the gluons among the two quark lines
       tot_ord=2*factorial8(next-4)*(next-3)

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
                   do i=2,next-1
                      if ((io(i,2).gt.2 .and. part(io(i,2)).le.-1 .and. part(io(i,2)).ge.-6) .or. &
                           (io(i,2).le.2 .and. part(io(i,2)).ge. 1 .and. part(io(i,2)).le. 6) ) then
                         if ( ((io(1   ,2).le.2 .and. io(i  ,2).gt.2) .or. (io(1   ,2).gt.2 .and. io(i  ,2).le.2)) .or. &
                              ((io(next,2).le.2 .and. io(i+1,2).gt.2) .or. (io(next,2).gt.2 .and. io(i+1,2).le.2)) ) then
                            ! not both initial or both final state. Cannot invert order.
                            cycle do_i_inv
                         endif
                         io(2:i-1,3)=io(i-1:2:-1,2)
                         io(i+2:next-1,3)=io(next-1:i+2:-1,2)
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
!!$                      if (.not.same_flavour) cycle
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
