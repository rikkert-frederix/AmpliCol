
program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_gen23
  use phase_space_genpt
  use haag
  use math_functions
  implicit none
  integer :: j,c_o,i,c_o_t,c_o_i,c_o_j,c_o_k
  integer(kind=8) :: sym_fac
  real*4 :: tBefore,tAfter,tTot_A,tTot_B
  integer(kind=4),dimension(:),allocatable :: o,part
  real(kind=8),dimension(:),allocatable :: mass
  real(kind=8) :: s_cut(2),sqrts
  logical :: t_chan
  character(len=30) :: filename
  integer(kind=4) :: integration, nquarks
  logical,dimension(-6:7,2) :: ipdgs=.false.
  integer :: col_fac,nhel
  integer*8 :: iden
  integer(kind=4) :: n_events
 

  call get_run_arguments()
  call create_run_tag()

  allocate(mass(next))

  call cpu_time(tTot_B)

! relevant input parameters for integration
  ncalls0=-10000   ! Number of events to generate. (If negative, start
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

  mass(1:next)=0d0

  ! include pdfs?
  if (include_pdf) then
     call PDF_initialise
     ndim=ndim+2
     if (part(1).eq.21) then
        ipdgs(0,1)=.true.
     else
        ipdgs(part(1),1)=.true.
     endif
     if (part(2).eq.21) then
        ipdgs(0,2)=.true.
     else
        ipdgs(part(2),2)=.true.
     endif
  endif
  
  call cpu_time(tBefore)
  t_chan=.false.
  if (integration.eq.1) then
     call gen23_init(sqrts,next,mass,o,part,s_cut,t_chan,include_pdf)
  elseif  (integration.eq.2) then
     call  haag_init(sqrts,next,mass,o,part,s_cut,t_chan,include_pdf)
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
  do i=1,2
     if (part(i).eq.21) then
        iden=iden*8*2
     elseif (abs(part(i)).ge.1 .and. abs(part(i)).le.6) then
        iden=iden*3*2
     endif
  enddo
  iden=iden*factorial8(nfin_glu)

  ! initialize the amplitudes (sets up the imaps(), helicity maps,
  ! colour factors, etc.)
  call cpu_time(tBefore)
  call amps%init(1,next,part,o)
  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  ! Compute the leading colour factor
  if (nquarks.eq.0) then
     col_fac=3**next
  elseif (nquarks.eq.2) then
     col_fac=3**(next-1)
  else
     write (*,*) 'Leading colour factor not implemented'
  endif

  ! number of helicities to sum over
  nhel=amps%current_list(amps%n_cur)%nhel*amps%current_list(next)%nhel
  allocate(amp2_hel(1:nhel))

! Not so relevant mint-module parameters: only used in special cases.
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

  open(unit=14,file='pt_list.txt',status='replace')
  if (imode.le.1) then
     call mint(integrand)
  else
     call read_grids_from_file
     call gen(integrand,0,-1) ! initialise countersi
     filename='Outputs/events'//tag//'.lhe'
     open(unit=11,file=filename,status='unknown')
     do j=1,abs(ncalls0)
        call gen(integrand,1,2) ! generate an unweighted event
        call unwgt_helicity
        call write_event(11,ans(1,0))
     enddo
     close(11)
     call gen(integrand,3,-1) ! print counters
  endif
     
  close(14)
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
  function integrand(x,vol,ifirst,f1)
    implicit none
    real*8 :: integrand
    integer :: ifirst
    real*8, dimension(ndim) :: x
    real*8, dimension(nintegrals) :: f1
    real*8, save :: val
    integer :: icol,iperm,jperm,ih
    real*8 :: vol,xmu_fac,cuts_wgt
    real*8, dimension(-6:7,2) :: PDF
    real*8, parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real*4 :: tBefore,tAfter

    double precision :: y,frac,steep,cuts_wgt_1
   
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

    call cpu_time(tBefore)
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

    cuts_wgt=pass_cuts(next,p)
    if ((jac.lt.0d0) .or. (smooth_cuts .and. cuts_wgt.lt.0d0) .or. (.not.smooth_cuts .and. cuts_wgt.lt.1d0)) then
       pass_cuts_check=.false.
       val=0d0
       return
    endif
    
    passed = passed + 1



    ! compute amplitudes
    call cpu_time(tBefore)
    call amps%evaluate(next,p,0)
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

    ! include the jacobian from vegas ('vol') and the wgt from the phase-space ('jac')
    weight=vol*jac*(4*pi*alphas)**(next-2)/dble(iden)*conv
    val=amp2*weight

    frac=0.8d0
    steep=0.01d0
    i=5
    y=(pt(p(0,i))-frac*pt_min)/(pt_min*(1d0-frac))
    if (pt(p(0,i)).gt.frac*pt_min.and.pt(p(0,i)).lt.pt_min) then
      cuts_wgt_1=((steep)*y/(steep+1d0-y))
    elseif (pt(p(0,i)).gt.pt_min) then
      cuts_wgt_1 = 1d0
    elseif (pt(p(0,i)).lt.frac*pt_min) then
      cuts_wgt_1 = 0d0
    endif
    !if (pt(p(0,3)).lt.pt_min) then
    !   if (cuts_wgt_1.gt.0d0) write(*,*) 'STOP'
    !endif
    write(14,*) pt(p(0,i)),cuts_wgt_1

    ! Apply the weight from the cuts
    if (smooth_cuts) val=val*cuts_wgt

    ! Since we only need to include a subset of all the colour-orderings, we
    ! need to compensate with a symmetry factor
    val=val*sym_fac


    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore

    
    if (include_PDF) then
       ! Include the PDFs
       xmu_fac=91.188d0 ! factorisation scale
       call PDF_eval(1,ipdgs(-6,1),xbjrk(1),xmu_fac,PDF(-6,1))
       call PDF_eval(1,ipdgs(-6,2),xbjrk(2),xmu_fac,PDF(-6,2))
       if (part(1).eq.21) then
          val=val*PDF(0,1)
       else
          val=val*PDF(part(1),1)
       endif
       if (part(2).eq.21) then
          val=val*PDF(0,2)
       else
          val=val*PDF(part(2),2)
       endif
    endif

    ! pass the result to the mint module
    f1(1)=abs(val)
    f1(2)=val

  end function integrand

  double precision function pass_cuts(n,p)
    ! Cuts on the phase-space point.
    implicit none
    integer :: i,j,n
    real*8,dimension(0:3,n) :: p
    double precision :: frac,y,steep

    frac=0.8d0
    steep=0.1d0
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
       if (pt_min.gt.0d0) then
          if (pt(p(0,i)).lt.frac*pt_min) then
             pass_cuts=-1d0
             return
          endif
          if (pt(p(0,i)).gt.frac*pt_min.and.pt(p(0,i)).lt.pt_min) then
             y=(pt(p(0,i))-frac*pt_min)/(pt_min*(1d0-frac))
             if (imode.le.0) then
               pass_cuts=pass_cuts*((steep)*y/(steep+1d0-y)) ! 1/x damping function
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
          write (iunit,*) part(i) ,p(1:3,i),p(0,i)
       else
          write (iunit,*) part(i) ,p(1:3,i),p(0,i)
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
    integer :: argc,start,end,glu
    character(len=256) :: argv
    integer, dimension(:), allocatable :: process,ord
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    if (argc.ne.2) then
       write(*,*)  'imode'
       write(*,*) 'integration mode (1, 2 or 3):'
       read (*,*)  imode,integration
    else
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,*) imode
          if (i.eq.2) read(argv,*) integration
       enddo
    endif

    open (unit=99, file='process.txt', status='old', action='read')
    read(99, *) next
    allocate(process(next))
    allocate(o(next))
    allocate(part(next))
    allocate(ord(next))
    read(99, *) process
    part=process
    read(99, *) ord
    nquarks = 0
    do i=1,next
       if (abs(process(i)).ge.1 .and. abs(process(i)).le.6) then
           nquarks=nquarks+1
       endif
       if ((i.le.2) .and. (abs(process(i)).ge.1 .and. abs(process(i)).le.6))  then
          process(i)=-process(i)
       endif
    enddo

    o=ord
    if (nquarks.eq.0) then
      do i=1,next
        if (ord(i).eq.1) start=i
        if (ord(i).eq.2) end=i
      enddo
      c_o=abs(end-start)-1
      c_o_t=0
      c_o_k=0
      c_o_i=abs(end-start)-1
      c_o_j=next-2-c_o_i
    elseif (nquarks.eq.2) then
      c_o=0 ! dummy value
      glu=1
      do i=1,next
        if (process(i).lt.0) then
          o(next)=i
          end=i
        endif
        if (process(i).gt.0 .and. process(i).ne.21) then
        o(1)=i
        start=i
        endif
        if (process(i).eq.21) then
            o(1+glu)=i
            glu=glu+1
        endif
      enddo
      if ((ord(next).eq.end) .and. (ord(1).eq.start)) then
        write(*,*) 'VALID ORDER!!!'
        o=ord ! the input order was a valid one, use that instead
      endif

      write(*,*) o
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
       if (c_o*2.eq.(next-2)) then
          sym_fac=factorial8(next-2)
       else
          sym_fac=2*factorial8(next-2)
       endif
    elseif (nquarks.eq.2) then
       if ((abs(process(1)).ge.1 .and. abs(process(1)).le.6) .and. &
           (abs(process(2)).ge.1 .and. abs(process(2)).le.6) )then
          ! quark and anti-quark are incoming. Only 1 channel needed,
          ! which would result in the following symmetry factor:
          sym_fac=factorial8(next-2)
       elseif ((abs(process(1)).ge.1 .and. abs(process(1)).le.6) .or. &
               (abs(process(2)).ge.1 .and. abs(process(2)).le.6) )then
          ! one incoming quark (or anti-quark). There are ngluons
          ! channels needed: they correspond to having the incoming
          ! gluon at all possible positions between the quark and
          ! anti-quark in the colour order. Hence, each channel comes
          ! with an (ngluons-1)! symmetry factor:
          sym_fac=factorial8(next-3)
       else
          ! both quark and anti-quark are final state. This is similar
          ! to the all-gluon case above, treating the q-qbar pair as
          ! another gluon. This special gluon is identifiable! So, for
          ! next=6 (and assuming that the qqbar pair are particles 5
          ! and 6) one has the following possibilities:
          !
          ! ia   --> 1,2,3,4,(5,6)  ---- : both gluons on the same
          ! ib   --> 1,2,3,(5,6),4  --/         line as the qqbar pair
          ! ic   --> 1,2,(5,6),3,4  -/
          ! iia  --> 1,3,2,4,(5,6)  ---- : one gluon on the same 
          ! iib  --> 1,3,2,(5,6),4  -/          line as the qqbar pair
          ! iii  --> 1,3,4,2,(5,6)  ---- : both gluons on the other quark line
          !
          ! Furthermore all these can have the quark and anti-quark
          ! order reversed, so there are in total 12 truly different
          ! colour orders to consider.
          !
          ! All these come with a symmetry factor of (ngluon-2)! =
          ! 2!. Hence we have:
          sym_fac=factorial8(next-4)
       endif
    else
       write (*,*) 'WARNING: symmetry factor missing',nquarks
    endif
    
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (imode.ne.0 .and. imode.ne.1 .and. imode.ne.2) then
       write (*,*) 'Incorrect imode',imode
       stop
    endif
    if (c_o.lt.0 .or. c_o .gt. next-2) then
       write (*,*) 'inconsistent color-ordering',c_o
       stop
    endif
    if (integration.ne.1 .and. integration.ne.2 .and. integration.ne.3) then
       write (*,*) 'Integration modes only 1, 2 or 3',integration
       stop
    endif
    if ((nquarks.ne.0 .and. nquarks.ne.2) .or. (nquarks.gt.next)) then
       write (*,*) 'Not consistent number of external quarks (up to 2)',nquarks
       stop
    endif

  end subroutine get_run_arguments

  subroutine create_run_tag()
    implicit none
    character(len=1) :: s1
    character(len=2) :: s2
    if (next.le.9) then
       write(s1,'(i1)') next
       tag=trim(adjustl(s1))//'_'
       tag_read=trim(adjustl(s1))//'_'
    else
       write(s2,'(i2)') next
       tag=trim(adjustl(s2))//'_'
       tag_read=trim(adjustl(s2))//'_'
    endif
    write(s1,'(i1)') imode
    tag=trim(adjustl(tag))//trim(adjustl(s1))//'_'
    if (imode.gt.0) write(s1,'(i1)') imode-1
    tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))//'_'
    write(s1,'(i1)') c_o_t
    tag=trim(adjustl(tag))//trim(adjustl(s1))//'_'
    tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))//'_'
    if (c_o_i.le.9) then
       write(s1,'(i1)') c_o_i
       tag=trim(adjustl(tag))//trim(adjustl(s1))//'_'
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))//'_'
    else
       write(s2,'(i2)') c_o_i
       tag=trim(adjustl(tag))//trim(adjustl(s2))//'_'
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s2))//'_'
    endif
    if (c_o_j.le.9) then
       write(s1,'(i1)') c_o_j
       tag=trim(adjustl(tag))//trim(adjustl(s1))//'_'
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))//'_'
    else
       write(s2,'(i2)') c_o_j
       tag=trim(adjustl(tag))//trim(adjustl(s2))//'_'
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s2))//'_'
    endif
    if (c_o_k.le.9) then
       write(s1,'(i1)') c_o_k
       tag=trim(adjustl(tag))//trim(adjustl(s1))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))
    else
       write(s2,'(i2)') c_o_k
       tag=trim(adjustl(tag))//trim(adjustl(s2))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s2))
    endif
    write(*,*) len(trim(tag_read))
    if (len(trim(tag_read)).lt.13) then
       if (13-len(trim(tag)).eq.1) then
          tag='_'//trim(adjustl(tag))
          tag_read='_'//trim(adjustl(tag_read))
       elseif(13-len(trim(tag)).eq.2) then
          tag='__'//trim(adjustl(tag))
          tag_read='__'//trim(adjustl(tag_read))
       elseif(13-len(trim(tag)).eq.3) then
          tag='___'//trim(adjustl(tag))
          tag_read='___'//trim(adjustl(tag_read))
       endif
    endif
    write (*,*) tag
  end subroutine create_run_tag

end program matrix_integrate_QCD
