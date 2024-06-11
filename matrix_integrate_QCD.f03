
program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_gen23
  use phase_space_genpt
  use haag
  use math_functions
  implicit none
  integer :: j,c_o,i,c_o_t,c_o_i,c_o_j,c_o_k
  integer(kind=8) :: sym_fac,iden
  real*4 :: tBefore,tAfter,tTot_A,tTot_B
  integer(kind=4),dimension(:),allocatable :: o,part,orig_part,part_sf
  real(kind=8),dimension(:),allocatable :: mass,width
  real(kind=8) :: s_cut(2),sqrts
  logical :: t_chan
  character(len=80) :: filename
  integer(kind=4) :: integration, nquarks
  logical,dimension(-6:7,2) :: ipdgs
  integer :: col_fac,nhel
  integer :: it ! quark order type
 
  call get_run_arguments()
  call compute_multichannel_symmetry_factor()
  call create_run_tag()

  allocate(mass(next))
  allocate(width(next))

  call cpu_time(tTot_B)

! relevant input parameters for integration
!!$  ncalls0=-10000   ! Number of events to generate. (If negative, start
  ncalls0=1000000   ! Number of events to generate. (If negative, start
                   ! from a small number of points and double it each
                   ! iteration. If positive, this is the number of
                   ! points per iteration as well).

  ndim=3*(next-2)-4   ! Number of dimensions of the integration.

  itmax=20         ! Number of iterations. (If ncalls0 < 0, the
                   ! integration is aborted if accuracy (next line)
                   ! has been reached.

  accuracy=0.003d0 ! Accuracy of the integration. (Ignored if ncalls0 > 0).


! relevant physics input parameters and initialisation of amplitudes
  sqrts=14000.d0
  ! setting energy

  s_cut(1)=max(sqrt_s_min,pt_min)**2
  s_cut(2)=max(sqrt_s_min**2,2d0*pt_min**2*(1d0-cos(DRjj_min)))

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
     call gen23_init(sqrts,next,mass,o,part,s_cut,pt_min,DRjj_min,t_chan,include_pdf)
  elseif  (integration.eq.2) then
     call haag_init(sqrts,next,mass,o,part,s_cut,t_chan,include_pdf)
  elseif (integration.eq.3) then
     call genpt_init(sqrts,next,mass,pt_min,eta_max,DRjj_min,include_pdf)
  endif
  call cpu_time(tAfter)
  t_PS_init=t_PS_init+tAfter-tBefore

  ! colour, polarisation incoming gluons: 8, 2
  ! colour, polarisation incoming quarks: 3, 2
  ! identical final state particle factor (gluons): nfin_glu!
  nfin_glu=0
  do i=3,next
     if (part(i).eq.21) then
        nfin_glu=nfin_glu+1
     endif
  enddo

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
  call amps%init(1,next,orig_part,part,mass,width,o,it)

  if (amps%n_qqbar.eq.2.and.amps%same_flav) then
    part_sf(:) = orig_part(:)
    call define_symm_2qq(next,part_sf,2)
    call amps_sf%init(1,next,orig_part,part_sf,mass,width,o,it)
  endif

  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  ! Compute the leading colour factor
  if (amps%n_qqbar.eq.2) then
      if (abs(part(o(1))).ne.abs(part(o(next)))) it = 2
  endif

  call compute_LC_colour_factor(col_fac,it)
  
  ! number of helicities to sum over
  nhel=amps%current_list(amps%n_cur)%nhel*amps%current_list(next)%nhel
  allocate(amp2_hel(1:nhel))

  ! Not so relevant mint-module parameters: only used in special cases.
  call set_mint_module_special_parameters()

  if (imode.le.1) then
     ! grid setup, or computation of upper bounding envelope
     call mint(integrand)
  else
     ! actual (unweighted) event generation
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise counters
     filename='Outputs/events'//trim(tag)//'.lhe'
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
  function integrand(x,vol,ifirst,nit,f1)
    implicit none
    real*8 :: integrand
    integer :: ifirst
    real*8, dimension(ndim) :: x
    real*8, dimension(nintegrals) :: f1
    real*8, save :: val
    integer :: iperm,ih
    real*8 :: vol,cuts_wgt
    real*8, parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real*4 :: tBefore,tAfter
    integer :: nit ! iteration number
    double precision :: y,frac,steep,cuts_wgt_1,Q
   
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
    if (integration.eq.1)then
       call gen23_phase_space(x)
    elseif (integration.eq.2) then
        call PS_haag(x)
    elseif (integration.eq.3) then
        call genpt_phase_space(x)
    endif
    
    call cpu_time(tAfter)
    t_PS= t_PS +tAfter-tBefore

    all_evt=all_evt+1

    cuts_wgt=pass_cuts(next,p,nit)
    if ((jac.lt.0d0) .or. (smooth_cuts .and. cuts_wgt.lt.0d0) .or. (.not.smooth_cuts .and. cuts_wgt.lt.1d0)) then
       pass_cuts_check=.false.
       val=0d0
       return
    endif
    
    passed = passed + 1

    ! compute amplitudes
    call cpu_time(tBefore)

    call amps%evaluate(next,p,mass,width,0,part)

    if (amps%n_qqbar.eq.2.and.amps%same_flav) then
      call amps_sf%evaluate(next,p,mass,width,0,part_sf)
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
    enddo
    amp2=sum(amp2_hel(1:nhel))

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

    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore
    
    if (include_PDF) then
       call multiply_by_PDF_value(val)
    endif

    ! pass the result to the mint module

    f1(1)=abs(val)
    f1(2)=val

  end function integrand

  double precision function pass_cuts(n,p,nit)
    ! Cuts on the phase-space point.
    implicit none
    integer :: i,j,n
    real*8,dimension(0:3,n) :: p
    double precision :: frac,y,steep
    integer :: nit

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
  
  real*8 function pt(p)
    ! transverse momentum of 'p'
    implicit none
    real*8, dimension(0:3) :: p
    pt=sqrt(p(1)**2+p(2)**2)
  end function pt
  
  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot

  real*8 function eta(p)
    ! pseudo-rapidity of 'p'
    implicit none
    real*8, dimension(0:3) :: p
    real*8 :: theta
    theta=acos(p(3)/sqrt(p(1)**2+p(2)**2+p(3)**2))
    eta=-log(dtan(theta/2d0))
  end function eta

  real*8 function delta_phi(p1,p2)
    ! azimuthal difference of 'p1' and 'p2'
    implicit none
    real*8, dimension(0:3) :: p1,p2
    real*8 :: denom
    denom=pt(p1)*pt(p2)
    delta_phi=acos((p1(1)*p2(1)+p1(2)*p2(2))/denom)
  end function delta_phi

  real*8 function deltaR(p1,p2)
    ! Distance (Delta-R) between 'p1' and 'p2'
    implicit none
    real*8, dimension(0:3) :: p1,p2
    deltaR=sqrt(delta_phi(p1,p2)**2+(eta(p1)-eta(p2))**2)
  end function deltaR

  subroutine write_event(iunit,wgt)
    implicit none
    integer :: i,iunit
    real(kind=8) :: wgt
    write (iunit,*) '<event>'
    write (iunit,*) next,hel_picked,wgt,amp2*weight,amp2,weight
    write (iunit,'(100i3)') o(1:next)
    do i=1,next
       if (i.le.2) then
          write (iunit,*) orig_part(i) ,p(1:3,i),p(0,i)
       else
          write (iunit,*) orig_part(i) ,p(1:3,i),p(0,i)
       endif
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event

  subroutine unwgt_helicity
    implicit none
    integer :: i
    real*8 :: random
    real*8,external :: ran2
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
    hel_picked=i
    if (hel_picked.gt.nhel) then
       write (*,*) 'Could not unweight helicity',hel_picked,nhel
       stop 1
    endif
  end subroutine unwgt_helicity
  
  subroutine get_run_arguments()
    implicit none
    integer :: argc
    integer :: i,k
    character(len=256) :: argv
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

    if (read_from_file) then
       call read_process_from_file
    !else
    !   call get_process_from_arguments
    endif

    ! basic checks:
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (imode.ne.0 .and. imode.ne.1 .and. imode.ne.2) then
       write (*,*) 'Incorrect imode',imode
       stop
    endif
    if (integration.ne.1 .and. integration.ne.2) then
       write (*,*) 'Integration modes only 1 or 2',integration
       stop
    endif
    if ((nquarks.ne.0 .and. nquarks.ne.2 .and. nquarks.ne.4) .or. (nquarks.gt.next)) then
       write (*,*) 'Not consistent number of external quarks (up to 2)',nquarks
       stop
    endif
  end subroutine get_run_arguments

  subroutine create_run_tag()
    implicit none
    tag=''       ! tag of current run
    tag_read=''  ! same as 'tag', but with previous imode (i.e., defines the file to read the integration grids from)
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
    do i=1,next
       call add_to_string(tag,o(i),.true.)
       call add_to_string(tag_read,o(i),.true.)
    enddo
    call fill_string(tag,len(trim(tag)))
    call fill_string(tag_read,len(trim(tag)))
    write (*,*) 'File tag is: ',tag
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
  subroutine fill_string(string,size)
    ! Fills the string 'string' with leading underscores until the string has
    ! size 'size'. The declaration of the string must be at least size 'size'.
    implicit none
    character(len=string_len) :: string
    integer :: size,n_to_add
    if (size.gt.len(string)) then
       write (*,*) 'Size greater than string',size,string
       stop 1
    endif
    n_to_add=len(trim(string))+2-len(trim(string))
    do i=1,n_to_add
       string='_'//trim(adjustl(string))
    enddo
  end subroutine fill_string

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
    real*8, dimension(-6:7,2) :: PDF
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

    nq=0
    naq=0
    ! count the number of final state gluons
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
          sym_fac=factorial8(ngl) !2*factorial8(ngl)  ! TO CHANGE: put back factor 2 here
       endif
    elseif (nquarks.eq.2) then
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
          ! Furthermore all these can have the quark and anti-quark order
          ! reversed (or, equivalently, the two incoming particles
          ! interchanged in the colour order), so there are in total 12 truly
          ! different colour orders to consider.
          !
          ! All these come with a symmetry factor of (ngluon-2)! =
          ! 2!. Hence we have:
          sym_fac=factorial8(ngl)
       endif
    elseif (nquarks.eq.4) then
       sym_fac=factorial8(ngl)
    else        
       write (*,*) 'WARNING: symmetry factor missing',nquarks
    endif
  end subroutine compute_multichannel_symmetry_factor

  subroutine read_process_from_file
    implicit none
    integer :: i,end,start,glu
    integer,dimension(:),allocatable :: ord
    open (unit=99, file='process.txt', status='old', action='read')
    read(99, *) next
    allocate(part(next))
    allocate(o(next))
    read(99, *) part
    read(99, *) o
    nquarks = 0
    do i=1,next
       if (abs(part(i)).ge.1 .and. abs(part(i)).le.6) then
          nquarks=nquarks+1
          if (i.le.2) part(i)=-part(i)
       endif
    enddo
    if (nquarks.eq.0) then
       do i=1,next
          if (o(i).eq.1) start=i
          if (o(i).eq.2) end=i
       enddo
       c_o_t=0
       c_o_k=0
       c_o_i=abs(end-start)-1
       c_o_j=next-2-c_o_i
       c_o=min(c_o_i,c_o_j)
    elseif (nquarks.eq.2) then
       do i=1,next
          if (o(i).eq.1) start=i
          if (o(i).eq.2) end=i
       enddo
       if (start.lt.end) then
          c_o_t=1
          if (start.ne.1) c_o_i=abs(start-1)-1
          if (start.eq.1) c_o_i=next+1
          if (end.ne.next) c_o_k=abs(next-end)-1
          if (end.eq.next) c_o_k=next+1
       endif
       if (start.gt.end) then
          c_o_t=2
          if (end.ne.1) c_o_k=abs(end-1)-1
          if (end.eq.1) c_o_k=next+1
          if (start.ne.next) c_o_i=abs(next-start)-1
          if (start.eq.next) c_o_i=next+1
       endif
       c_o_j=abs(start-end)-1
    endif
  end subroutine read_process_from_file

  subroutine check_input_colour_order_consistency
    ! Checks that the c_o_t, c_o_i, c_o_j, and c_o_k are in their correct
    ! ranges.
    implicit none
    if (c_o_t.eq.0 .and. (c_o_j.gt.int((next-2)/2) .or. c_o_j.lt.0)) then
       write (*,*) 'Inconsistent colour order #1', c_o_t,c_o_i,c_o_j,c_o_k
       stop 1
    elseif (c_o_i.eq.next+1 .and. c_o_k.eq.next+1) then
       if (c_o_t.ne.1 .and. c_o_t.ne.2) then
          write (*,*) 'Inconsistent colour order #2a', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_j.lt.0 .or. c_o_j.gt.next-2) then
          write (*,*) 'Inconsistent colour order #2b', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
    elseif (c_o_i.eq.next+1) then
       if (c_o_t.lt.1 .or. c_o_t.gt.4) then
          write (*,*) 'Inconsistent colour order #3a', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_k.gt.next-3 .or. c_o_k.lt.0) then
          write (*,*) 'Inconsistent colour order #3b', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_j.gt.next-3 .or. c_o_j.lt.0) then
          write (*,*) 'Inconsistent colour order #3c', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_k+c_o_j.gt.next-3) then
          write (*,*) 'Inconsistent colour order #3d', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
    elseif (c_o_k.eq.next+1) then
       if (c_o_t.lt.1 .or. c_o_t.gt.4) then
          write (*,*) 'Inconsistent colour order #4a', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_i.gt.next-3 .or. c_o_i.lt.0) then
          write (*,*) 'Inconsistent colour order #4b', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_j.gt.next-3 .or. c_o_j.lt.0) then
          write (*,*) 'Inconsistent colour order #4c', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_i+c_o_j.gt.next-3) then
          write (*,*) 'Inconsistent colour order #4d', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
    elseif (c_o_t.ge.1) then
       if (c_o_t.lt.1 .or. c_o_t.gt.8) then
          write (*,*) 'Inconsistent colour order #5a', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_i.gt.next-4 .or. c_o_i.lt.0) then
          write (*,*) 'Inconsistent colour order #5b', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_j.gt.next-4 .or. c_o_j.lt.0) then
          write (*,*) 'Inconsistent colour order #5c', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_k.gt.next-4 .or. c_o_k.lt.0) then
          write (*,*) 'Inconsistent colour order #5d', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
       if (c_o_i+c_o_j+c_o_k.gt.next-4) then
          write (*,*) 'Inconsistent colour order #5e', c_o_t,c_o_i,c_o_j,c_o_k
          stop 1
       endif
    endif
  end subroutine check_input_colour_order_consistency
  
  subroutine get_process_from_arguments
    implicit none
    ! c_o_t == 0 --> all gluon process; c_o_j is the minimum between the
    !                numbers of gluons on each of the two colour lines that
    !                connect the two incoming gluons
    !
    ! c_o_t > 0 --> there is 1 qqbar particle:
    !      c_o_i == n+1 && c_o_k == n+1 --> the two quarks are in the initial
    !                state. (c_o_t==1 or 2 defines if it is qqbar->ngluons or
    !                qbarq->ngluons). c_o_j defines the number of gluons that
    !                are colour-ordered; hence if this is less than the number
    !                of gluons in the process, the remaining ones are photons.
    !
    !      c_o_i == n+1 && c_o_k < n-2 --> left incoming particle is a quark
    !                (if c_o_t==2) or anti-quark (if c_o_t==1) and the other
    !                incoming particle a gluon (add 2 to c_o_t to make this a
    !                photon). c_o_j defines the number of gluons between the
    !                incoming (anti-)quark and the incoming gluon in the
    !                colour order; c_o_k is the number of gluons between the
    !                incoming gluon and the final state anti-quark.
    !
    !      c_o_i < n-2 && c_o_k == n+1 --> Same as previous, but with the
    !                right-incoming particle the (anti-quark), and the
    !                left-incoming the gluon.
    !
    !      c_o_i < n-2 && c_o_k < n-2 --> Both the quark and anti-quark are
    !                final state. As before, c_o_i is the number of gluons in
    !                the colour order between the left-incoming gluon and the
    !                quark; c_o_k between the right-incoming gluon and the
    !                anti-quark; and c_o_j the number of gluons between the
    !                two incoming gluons. If c_o_t==2, the role of
    !                left-incoming and right-incoming is interchanged. If
    !                c_o_i+c_o_j+c_o_k < n-4, the remaining particles are
    !                photons. Add 2 to c_o_t to make the left incoming initial
    !                state particle a photon; Add 4 to c_o_t to make the right
    !                incoming initial state particle a photon; add 6 to c_o_t
    !                to make both incoming particles photons.
    !
    ! NOTE: for qqbar process, the colour order is such that the first
    ! particle should be a quark (if final state) or anti-quark (if initial
    ! state), while the last particle is that anti-quark (if final state) or
    ! quark (if initial state). This is even true when there are photons
    ! around: these photons are always put just before the final particle in
    ! the colour (i.e., just before the anti-quark (if final state) or quark
    ! (if initial state)). THE LATTER IS SOMEWHAT COUNTER-INTUITIVE, since the
    ! photons shouldn't be part of the colour order whatsoever...
    integer :: i,k
    integer,dimension(:),allocatable :: ord
    integer :: nsing
    !call check_input_colour_order_consistency
    if (c_o_t.eq.0) then
       nquarks=0       
    else
       nquarks=2
       ! count the number of colour-singlets:
       if (c_o_i.eq.next+1 .and. c_o_k.eq.next+1) then
          ! We gave next-2 final state gluons, of which c_o_j are colour-ordered
          nsing=(next-2)-c_o_j
       elseif(c_o_i.eq.next+1 .and. c_o_k.ne.next+1) then
          ! We have next-3 final state gluons, of which (c_o_j+c_o_k) are colour-ordered
          nsing=(next-3)-(c_o_j+c_o_k)
       elseif(c_o_i.ne.next+1 .and. c_o_k.eq.next+1) then
          ! We have next-3 final state gluons, of which (c_o_i+c_o_j) are colour-ordered
          nsing=(next-3)-(c_o_i+c_o_j)
       else
          ! We have next-4 final state gluons, of which (c_o_i+c_o_j+c_o_k) are colour ordered
          nsing=(next-4)-(c_o_i+c_o_j+c_o_k)
       endif
    endif
       
    allocate(part(next))
    allocate(o(next))
    allocate(ord(next))
    if (nquarks.eq.2) then
       if (c_o_i.eq.next+1.and.c_o_k.eq.next+1) then
          if (c_o_t.eq.1) then
             part(1)=-1
             part(2)=1
             ord(1)=1
             ord(next)=2
          elseif (c_o_t.eq.2) then
             part(1)=1
             part(2)=-1
             ord(1)=2
             ord(next)=1
          endif
          do i=3,next
             if (i.le.c_o_j+2) then
                part(i)=21
             else
                part(i)=22
             endif
          enddo
          do i=2,next-1
             ord(i)=i+1
          enddo
       elseif (c_o_i.eq.next+1.and.c_o_k.ne.next+1) then
          if (mod(c_o_t,2).eq.1) then
             part(1)=-1
             if (c_o_t.eq.1) then
                part(2)=21
             else
                part(2)=22
             endif
             part(3)=-1
             ord(1)=1
             ord(next)=3
             ord(2+c_o_j)=2
             k=4
             do i=2,2+c_o_j-1
                ord(i)=k
                k=k+1
             enddo
             do i=2+c_o_j+1,next-1
                ord(i)=k
                k=k+1
             enddo
          elseif (mod(c_o_t,2).eq.0) then
             part(1)=1
             if (c_o_t.eq.2) then
                part(2)=21
             else
                part(2)=22
             endif
             part(3)=1
             ord(1)=3
             ord(next)=1
             ord(2+c_o_k)=2
             k=4
             do i=2,2+c_o_k-1
                ord(i)=k
                k=k+1
             enddo
             do i=2+c_o_k+1,next-1
                ord(i)=k
                k=k+1
             enddo
          endif
          do i=4,next
             if (i.le.c_o_j+c_o_k+3) then
                part(i)=21
             else
                part(i)=22
             endif
          enddo
       elseif (c_o_i.ne.next+1.and.c_o_k.eq.next+1) then
          if (mod(c_o_t,2).eq.1) then
             if (c_o_t.eq.1) then
                part(1)=21
             else
                part(1)=22
             endif
             part(2)=1
             part(3)=1
             ord(1)=3
             ord(next)=2
             ord(2+c_o_i)=1
             k=4
             do i=2,2+c_o_i-1
                ord(i)=k
                k=k+1
             enddo
             do i=2+c_o_i+1,next-1
                ord(i)=k
                k=k+1
             enddo
          elseif (mod(c_o_t,2).eq.0) then
             if (c_o_t.eq.2) then
                part(1)=21
             else
                part(1)=22
             endif
             part(2)=-1
             part(3)=-1
             ord(1)=2
             ord(next)=3
             ord(2+c_o_j)=1
             k=4
             do i=2,2+c_o_j-1
                ord(i)=k
                k=k+1
             enddo
             do i=2+c_o_j+1,next-1
                ord(i)=k
                k=k+1
             enddo
          endif
          do i=4,next
             if (i.le.c_o_j+c_o_k+3) then
                part(i)=21
             else
                part(i)=22
             endif
          enddo
       else
          if (c_o_t.le.2) then
             part(1)=21
             part(2)=21
          elseif (c_o_t.le.4) then
             part(1)=22
             part(2)=21
          elseif (c_o_t.le.6) then
             part(1)=21
             part(2)=22
          elseif (c_o_t.le.8) then
             part(1)=22
             part(2)=22
          endif
          part(3)=1
          part(4)=-1
          do i=5,next
             if (i.le.c_o_i+c_o_j+c_o_k+4) then
                part(i)=21
             else
                part(i)=22
             endif
          enddo
          ord(1)=3
          ord(next)=4
          if (mod(c_o_t,2).eq.1) then
             ord(2+c_o_i)=1
             ord(3+c_o_i+c_o_j)=2
          else
             ord(2+c_o_i)=2
             ord(3+c_o_i+c_o_j)=1
          endif
          k=5
          do i=2,2+c_o_i-1
             ord(i)=k
             k=k+1
          enddo
          do i=2+c_o_i+1,3+c_o_i+c_o_j-1
             ord(i)=k
             k=k+1
          enddo
          do i=3+c_o_i+c_o_j+1,next-1
             ord(i)=k
             k=k+1
          enddo
       endif
    elseif (nquarks.eq.0) then
       do i=1,next
          part(i)=21
       enddo
       ord(1)=1
       ord(2+c_o_j)=2
       k=3
       do i=2,2+c_o_j-1
          ord(i)=k
          k=k+1
       enddo
       do i=2+c_o_j+1,next
          ord(i)=k
          k=k+1
       enddo
       c_o=min(next-2-c_o_j,c_o_j)
       if (c_o_j.lt.0 .or. c_o_j.gt.next-2) then
          write(*,*) 'Incorrect colour order for all gluons: ',c_o_j
          stop
       elseif (c_o_i.gt.0 .or. c_o_k.gt.0) then
          write(*,*) 'Incorrect colour order for all gluons (c_i and c_k must be 0): ',c_o_i,c_o_k
          stop
       endif
    endif
    o=ord
    write (*,*) '******************************************'
    write (*,*) 'Process is     ',part
    write (*,*) 'Colour order is',o
    write (*,*) '******************************************'
  end subroutine get_process_from_arguments
  
end program matrix_integrate_QCD
