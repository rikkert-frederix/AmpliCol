
! gfortran -ffast-math -O3 -o matrix_integrate_QCD simple_mint/mint_module.f90 simple_mint/MC_integer.f simple_mint/ranmar.f simple_mint/HwU.f PhaseSpace_BycklingKajantie/LUPdecompose.f90 PhaseSpace_BycklingKajantie/phase_space_gen23.f90 PhaseSpace_haag/haag.f90 color_algebra.f95 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 amplitude_real.f03 matrix_integrate_QCD.f03


module common
  use amplitude_mod
  use amplitude_QCD_mod
  implicit none
  real*8,parameter  :: alphaS=0.12d0
  integer :: next,nfin,hel_picked

  type(amplitude) :: amplitudes
  type(amplitude_QCD)  :: amps

  ! timing
  real*4 :: t_PS_init=0.,t_Amp_init=0.,t_PS=0.,t_Amp=0.,t_all=0.,t_mat=0.
  real*8 :: amp2,weight
  real*8,dimension(:),allocatable :: amp2_hel
  real(kind=8),dimension(:,:),allocatable,public :: p
  real(kind=8),public :: jac
  

  ! counting events
  integer(kind=4) :: passed=0
  
end module common


program matrix_integrate_QCD
  use common
  use mint_module
  use phase_space_gen23
  use haag
  implicit none
  integer :: col_acc,j,c_o,i
  integer(kind=8) :: sym_fac
  real*4 :: tBefore,tAfter,tTot_A,tTot_B
  integer(kind=4),dimension(:),allocatable :: o,part
  real(kind=8),dimension(:),allocatable :: mass
  real(kind=8) :: s_cut(2),sqrtshat
  logical :: t_chan
  character(len=30) :: filename
  real(kind=8) :: sqrt_s_min,pt_min,drjj_min,eta_max
  integer(kind=4) :: integration, nquarks

  call get_run_arguments()
  call create_run_tag()

  allocate(mass(next))

!  nperm=3*2*1

!  allocate(helmap(2**n,nperm))
!  allocate(amps(nperm))
!  allocate(col_fac(nperm,nperm))

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
  sqrtshat=1000.d0

  pt_min=-1d0
  DRjj_min=-1d0
  eta_max=-1d0
  sqrt_s_min=30d0

  s_cut(1)=max(sqrt_s_min,pt_min)**2
  s_cut(2)=max(sqrt_s_min,pt_min*DRjj_min)**2

  mass(1:next)=0d0
  !do i=1,next
  !   if (i.eq.1) then
  !      o(i)=1
  !   elseif (i.lt.2+c_o) then
  !      o(i)=i+1
  !   elseif (i.eq.2+c_o) then
  !      o(i)=2
  !   else
  !      o(i)=i
  !   endif
  !enddo
  t_chan=.false.

  call cpu_time(tBefore)
  if (integration.eq.1) then
        call gen23_init(sqrtshat,next,mass,o,s_cut,t_chan)
  elseif  (integration.eq.2) then
        call  haag_init(sqrtshat,next,mass,o,s_cut,t_chan)
  endif

  call cpu_time(tAfter)

  t_PS_init=t_PS_init+tAfter-tBefore
  
  col_acc=0 ! Leading colour

  ! initialize the amplitudes (sets up the imaps(), helicity maps,
  ! colour factors, etc.)
  call cpu_time(tBefore)

  !if (nquarks.gt.0) then
  !  part(1)=-1
  !  part(next)=-1
  !  do i=2,next-1
  !    part(i) = 21
  !  enddo
  !else
  ! do i=1,next
  !    part(i)=21
  !  enddo
  !endif

  nfin=0
  do i=3,next
     if (part(i).eq.21) then
        nfin=nfin+1
      endif
  enddo

  call amps%init(1,next,part,o)

!  call amps%init_col(next,0) ! LC colour factors only
  allocate(amp2_hel(0:2**next))

  if (c_o*2.eq.(next-2)) then
     sym_fac=factorial8(next-2)
  else
     sym_fac=2*factorial8(next-2)
  endif
  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore


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
        call write_event(11,ans(1,0)*sym_fac)
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
  write(*,*) 'Number passing cuts:',passed
 
contains
  function integrand(x,vol,ifirst,f1)
    implicit none
    real*8 :: integrand
    integer :: ifirst
    real*8, dimension(ndim) :: x
    real*8, dimension(nintegrals) :: f1
    real*8, save :: val
    integer :: icol,iperm,jperm,ih
    integer*8 :: iden
    real*8 :: vol
    real*8, parameter :: pi=3.14159265358979323846d0,conv=389379660d0
    real*4 :: tBefore,tAfter
    real(kind=8),dimension(2) :: ztemp
    integer :: ih1, ih2
    integer :: col_fac


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
    endif
    
    call cpu_time(tAfter)
    t_PS= t_PS +tAfter-tBefore

    if ((jac.lt.0d0) .or. (.not.pass_cuts(next,p))) then
       pass_cuts_check=.false.
       val=0d0
       return
    endif
    passed = passed + 1

    ! colour, polarisation incoming gluons: 8*8, 2*2
    ! identical final state particle factor: nfin!

    iden=1
    do i=1,2
      if (part(i).eq.21) then
         iden=iden*8*2
      elseif (part(i).eq.-1) then
         iden=iden*3*2
      endif
    enddo
    iden=iden* factorial8(nfin)



    ! compute amplitudes
    call cpu_time(tBefore)

    !do iperm=1,nperm
    !p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
    !p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03/)
    !p(0:3,3)=(/0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03/)
    !p(0:3,4)=(/0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03/)

    call amps%evaluate(next,p,0)

!    write(*,*) 'amp eval',amps%amps(1)
    !enddo

    call cpu_time(tAfter)
    t_amp=t_amp+tAfter-tBefore
    call cpu_time(tBefore)
    amp2_hel=0d0

!    do icol=1,amplitudes%colmap(0,0)
!       iperm=amplitudes%colmap(1,icol)
!       jperm=amplitudes%colmap(2,icol)
!       do ih=0,amplitudes%nhel(amplitudes%isize+1)-1
!          amp2_hel(ih)=amp2_hel(ih)+amplitudes%amps(amplitudes%helmap(iperm,ih),iperm)* &
!               amplitudes%colmap(0,icol)* &
!               amplitudes%amps(amplitudes%helmap(jperm,ih),jperm)
!       enddo
!    enddo
 
    if (nquarks .gt.0) then
      col_fac = 3**(next-1)
    else
      col_fac=3**next
    endif

    !write(*,*) 'nqruarks',nquarks

    do ih1=1,amps%current_list(amps%n_cur)%nhel
      do ih2=1,amps%current_list(next)%nhel
        ih=(ih2-1)*amps%current_list(amps%n_cur)%nhel+ih1
        amp2_hel(ih)=amp2_hel(ih)+dble(amps%amps(amps%helmap(ih))*col_fac*dconjg(amps%amps(amps%helmap(ih))))
        !write(*,*) col_fac
      enddo
    enddo
    amp2=sum(amp2_hel(1:2**next))

    ! include the jacobian from vegas ('vol') and the wgt from the phase-space ('jac')
   
    weight=vol*jac*(4*pi*alphas)**nfin/dble(iden)*conv
    val=amp2*weight

    call cpu_time(tAfter)
    t_mat=t_mat+tAfter-tBefore

    ! pass the result to the mint module
    f1(1)=abs(val)
    f1(2)=val

  end function integrand

  logical function pass_cuts(n,p)
    ! Cuts on the phase-space point.
    implicit none
    integer :: i,j,n
    real*8,dimension(0:3,n) :: p
    pass_cuts=.true.
    if (sqrt_s_min.gt.0d0) then
       do i=1,n-1
          do j=i+1,n
             if (abs(2d0*dot(p(0,i),p(0,j))).lt.sqrt_s_min**2) then
                pass_cuts=.false.
                return
             endif
          enddo
       enddo
    endif
    do i=3,n
       if (pt_min.gt.0d0) then
          if (pt(p(0,i)).lt.pt_min) then
             pass_cuts=.false.
             return
          endif
       endif
       if (eta_max.gt.0d0) then
          if (abs(eta(p(0,i))).gt.eta_max) then
             pass_cuts=.false.
             return
          endif
       endif
       if (drjj_min.gt.0d0) then
          if (i.ne.n) then
             do j=i+1,n
                if (DeltaR(p(0,i),p(0,j)).lt.drjj_min) then
                   pass_cuts=.false.
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
    i=0
    do
       if (amp2_hel(i).gt.random) then
          exit
       else
          i=i+1
          amp2_hel(i)=amp2_hel(i)+amp2_hel(i-1)
       endif
    enddo
    hel_picked=i
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
       write(*,*) 'integration mode (1 or 2):'
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
       if ((abs(process(i)).ge.1) .and. abs(process(i)).le.6) then
           nquarks=nquarks+1
       endif
       if ((i.le.2) .and. ((abs(process(i)).ge.1) .and. abs(process(i)).le.6))  then
          process(i)=-process(i)
       endif
    enddo

    if (nquarks.eq.0) then
      do i=1,next
        if (ord(i).eq.1) start=i
        if (ord(i).eq.2) end=i
      enddo
      c_o=abs(end-start)-1
    else
      c_o=0 ! dummy value
    endif

    o=ord
    glu=1
    if (nquarks.eq.2) then
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
        o=ord ! the input order was a valid one, use that instead
      endif
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
    if (integration.ne.1 .and. integration.ne.2) then
       write (*,*) 'Integration modes only 1 or 2',integration
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
    if (c_o.le.9) then
       write(s1,'(i1)') c_o
       tag=trim(adjustl(tag))//trim(adjustl(s1))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))
    else
       write(s2,'(i2)') c_o
       tag=trim(adjustl(tag))//trim(adjustl(s2))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s2))
    endif
    if (len(trim(tag_read)).lt.8) then
       if (8-len(trim(tag)).eq.1) then
          tag='_'//trim(adjustl(tag))
          tag_read='_'//trim(adjustl(tag_read))
       elseif(8-len(trim(tag)).eq.2) then
          tag='__'//trim(adjustl(tag))
          tag_read='__'//trim(adjustl(tag_read))
       elseif(8-len(trim(tag)).eq.3) then
          tag='___'//trim(adjustl(tag))
          tag_read='___'//trim(adjustl(tag_read))
       endif
    endif
    write (*,*) tag
  end subroutine create_run_tag

end program matrix_integrate_QCD
