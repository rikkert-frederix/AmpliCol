! gfortran -fbounds-check -o test_QCD math_functions.f03 color_algebra.f95 feynmanrules.f03 amplitude_QCD.f03 amplitude_real.f03 test_QCD.f03 

program test_QCD
  use amplitude_mod
  use amplitude_QCD_mod
  use color_algebra
  implicit none
  type(amplitude_QCD),dimension(:),allocatable :: amps
  type(amplitude_QCD) :: amps_col
  type(amplitude_cache) :: ampplitudes_cache
  integer :: n
  integer,dimension(:),allocatable :: part
  integer,dimension(:,:),allocatable :: helmap,order
  integer :: i,ih,iden,iperm,jperm,nperm
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: amp2,t
  real(kind=8),parameter :: pi=3.14159265358979323846d0,alphas=0.118d0,alphaEW=7.547E-003
  real(kind=8),dimension(:,:),allocatable :: col_fac
  complex(kind=8) :: ztemp
  logical :: quarks, single_perm,photon
  real(kind=8) :: tBefore,tAfter,t_eval,t_init
  real(kind=8) :: Q

  n=4
  quarks=.true.
  photon=.true.
  if (.not.quarks) photon=.false.
  single_perm = .false.
  
  allocate(part(n))

  if (n.eq.4) then
     if (quarks) then
        nperm = 2
        if (photon) nperm = 1 
     else
        nperm = 3*2
     endif
     if (single_perm) nperm=1
  elseif (n.eq.5) then
     if (quarks) then
        nperm = 3*2   
        if (photon) nperm = 2
     else
        nperm = 4*3*2
     endif
     if (single_perm) nperm=1
  elseif (n.eq.6) then
     if (quarks) then
        nperm = 4*3*2   
     else
        nperm = 5*4*3*2
     endif
     if (single_perm) nperm=1
  endif

  allocate(order(n,nperm))
  allocate(helmap(2**n,nperm))
  allocate(amps(nperm))
  allocate(col_fac(nperm,nperm))


!!$  part(1:n)=21
!!$  part(1)=-1
!!$  part(n)=-1
!!$  call amps(1)%init_AllColOrder(n,part)
!!$
!!$  stop
  
  
  allocate(p(0:3,n))
  
  if (n.eq.4) then
     p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
     p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03/)
     p(0:3,3)=(/0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03/)
     p(0:3,4)=(/0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03/)
  elseif (n.eq.5) then
     p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
     p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  -0.5000000E+03/)
     p(0:3,3)=(/0.4959179E+03, -0.1999048E+02,  0.7981352E+02, -0.4890448E+03/)
     p(0:3,4)=(/0.1328543E+03, -0.2648250E+02, -0.4416981E+02,  0.1224662E+03/)
     p(0:3,5)=(/0.3712277E+03,  0.4647298E+02, -0.3564372E+02,  0.3665785E+03/)
  elseif (n.eq.6) then
     p(0:3,1)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,   500.00000000000000d0/)
     p(0:3,2)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,  -500.00000000000000d0/)
     p(0:3,3)=(/88.551333054502976d0,  -22.100690287689982d0,   40.080353191685326d0,  -75.805430956936632d0/)
     p(0:3,4)=(/328.32941922709853d0,  -103.84961188345629d0,  -301.93375538954007d0,   76.494921387165888d0/)
     p(0:3,5)=(/152.35810946743064d0,  -105.88095966659219d0,  -97.709638326975707d0,   49.548385226792817d0/)
     p(0:3,6)=(/430.76113825096763d0,   231.83126183773845d0,   359.56304052483051d0,  -50.237875657022109d0/)
  endif

  if (n.eq.4) then
     if (quarks) then
       !part(1:n)=[-1,21,21,-1]
       !iden=3*8*2*2
       !order(1:n,1)= [1,2,3,4]
       !if (.not. single_perm) then
       !  order(1:n,2)= [1,3,2,4]
       !  iden=iden*1
       !endif

       if (.not.photon) then
       part(1:n)=[-1,1,21,21]
       iden=3*3*2*2
       order(1:n,1)= [1,3,4,2]
       if (.not. single_perm) then
         order(1:n,2)= [1,4,3,2]
         iden=iden*2
       endif
       endif

       !part(1:n)=[21,21,-1,1]
       !iden=8*8*2*2
       !order(1:n,1)= [4,1,2,3]
       !if (.not. single_perm) then
       !  order(1:n,2)= [4,2,1,3]
       !  iden=iden*1
       !endif

       if (photon) then
       part(1:n)=[-1,21,-1,22]
       iden=3*8*2*2
       order(1:n,1)= [1,2,4,3]
       if (.not. single_perm) then
         iden=iden*1
       endif

       !part(1:n)=[-1,21,-1,22]
       !iden=3*7*2*2
       !order(1:n,1)= [1,2,4,3]
       !if (.not. single_perm) then
       !  iden=iden*1
       !endif
       endif

     else
       part(1:n)=[21,21,21,21]
       iden=8*8*2*2
       order(1:n,1)= [1,2,3,4]
       if (.not.single_perm) then
         order(1:n,2)= [1,3,2,4]
         order(1:n,3)= [3,1,2,4]
         order(1:n,4)= [3,2,1,4]
         order(1:n,5)= [2,1,3,4]
         order(1:n,6)= [2,3,1,4]
         iden=iden*2
       endif
     endif

  elseif (n.eq.5) then

    if (quarks) then
     if (photon) then
       part(1:n)=[-1,1,21,21,22]
       iden=3*3*2*2
       order(1:n,1)= [1,3,4,5,2]
       if (.not. single_perm) then
         order(1:n,2)= [1,4,3,5,2]
         iden=iden*2
       endif

     else
       part(1:n)=[-1,1,21,21,21]
       iden=3*8*2*2 ! initial status colours and helicities/polarisation
       order(1:n,1)=  [1,3,4,5,2]
       if (.not.single_perm) then
         order(1:n,2)=[1,3,5,4,2]
         order(1:n,3)=[1,4,3,5,2]
         order(1:n,4)=[1,4,5,3,2]
         order(1:n,5)=[1,5,3,4,2]
         order(1:n,6)=[1,5,4,3,2]
         iden=iden*2  ! final state identical particles
       endif

       !part(1:n)=[-1,1,21,21,21]
       !iden=3*3*2*2 ! initial status colours and helicities/polarisation
       !order(1:n,1)=[1,3,4,5,2]
       !if (.not.single_perm) then
       !  order(1:n,2)=[1,3,5,4,2]
       !  order(1:n,3)=[1,4,3,5,2]
       !  order(1:n,4)=[1,4,5,3,2]
       !  order(1:n,5)=[1,5,3,4,2]
       !  order(1:n,6)=[1,5,4,3,2]
       !  iden=iden*3*2   ! final state identical particles
       !endif

       !part(1:n)=[21,21,1,-1,21]
       !iden=8*8*2*2
       !  order(1:n,1)=[3,1,2,5,4]
       !if (.not.single_perm) then
       !  order(1:n,2)=[3,1,5,2,4]
       !  order(1:n,3)=[3,2,1,5,4]
       !  order(1:n,4)=[3,2,5,1,4]
       !  order(1:n,5)=[3,5,1,2,4]
       !  order(1:n,6)=[3,5,2,1,4]
       !  iden=iden*1 
       !endif

      endif

     else
       part(1:n)=[21,21,21,21,21]   
       iden=8*8*2*2  
       order(1:n,1)= [1,2,3,4,5]
       if (.not.single_perm) then
         order(1:n,2)= [1,2,3,5,4]
         order(1:n,3)= [1,2,4,3,5]
         order(1:n,4)= [1,2,4,5,3]
         order(1:n,5)= [1,2,5,3,4]
         order(1:n,6)= [1,2,5,4,3]
         order(1:n,7)= [1,3,2,4,5]
         order(1:n,8)= [1,3,2,5,4]
         order(1:n,9)= [1,4,2,3,5]
         order(1:n,10)=[1,4,2,5,3]
         order(1:n,11)=[1,5,2,3,4]
         order(1:n,12)=[1,5,2,4,3]
         order(1:n,13)=[1,3,4,2,5]
         order(1:n,14)=[1,3,5,2,4]
         order(1:n,15)=[1,4,3,2,5]
         order(1:n,16)=[1,4,5,2,3]
         order(1:n,17)=[1,5,3,2,4]
         order(1:n,18)=[1,5,4,2,3]
         order(1:n,19)=[1,3,4,5,2]
         order(1:n,20)=[1,3,5,4,2]
         order(1:n,21)=[1,4,3,5,2]
         order(1:n,22)=[1,4,5,3,2]
         order(1:n,23)=[1,5,3,4,2]
         order(1:n,24)=[1,5,4,3,2]
         iden=iden*3*2
       endif
     endif

     
  elseif (n.eq.6) then
     if (quarks) then
       part(1:n)=[1,21,21,21,21,1]
       iden=3*8*2*2
       iden=iden*3*2
     else
       part(1:n)=[21,21,21,21,21,21]
       order(1:n,1)=[1,2,3,4,5,6]
       order(1:n,2)=[1,2,3,5,4,6]
       order(1:n,3)=[1,2,4,3,5,6]
       order(1:n,4)=[1,2,4,5,3,6]
       order(1:n,5)=[1,2,5,3,4,6]
       order(1:n,6)=[1,2,5,4,3,6]
       order(1:n,7)=[1,3,2,4,5,6]
       order(1:n,8)=[1,3,2,5,4,6]
       order(1:n,9)=[1,3,4,2,5,6]
       order(1:n,10)=[1,3,4,5,2,6]
       order(1:n,11)=[1,3,5,2,4,6]
       order(1:n,12)=[1,3,5,4,2,6]
       order(1:n,13)=[1,4,2,3,5,6]
       order(1:n,14)=[1,4,2,5,3,6]
       order(1:n,15)=[1,4,3,2,5,6]
       order(1:n,16)=[1,4,3,5,2,6]
       order(1:n,17)=[1,4,5,2,3,6]
       order(1:n,18)=[1,4,5,3,2,6]
       order(1:n,19)=[1,5,2,3,4,6]
       order(1:n,20)=[1,5,2,4,3,6]
       order(1:n,21)=[1,5,3,2,4,6]
       order(1:n,22)=[1,5,3,4,2,6]
       order(1:n,23)=[1,5,4,2,3,6]
       order(1:n,24)=[1,5,4,3,2,6]

       i=24
       order(1:n,i+1)=[2,1,3,4,5,6]
       order(1:n,i+2)=[2,1,3,5,4,6]
       order(1:n,i+3)=[2,1,4,3,5,6]
       order(1:n,i+4)=[2,1,4,5,3,6]
       order(1:n,i+5)=[2,1,5,3,4,6]
       order(1:n,i+6)=[2,1,5,4,3,6]
       order(1:n,i+7)=[2,3,1,4,5,6]
       order(1:n,i+8)=[2,3,1,5,4,6]
       order(1:n,i+9)=[2,3,4,1,5,6]
       order(1:n,i+10)=[2,3,4,5,1,6]
       order(1:n,i+11)=[2,3,5,1,4,6]
       order(1:n,i+12)=[2,3,5,4,1,6]
       order(1:n,i+13)=[2,4,1,3,5,6]
       order(1:n,i+14)=[2,4,1,5,3,6]
       order(1:n,i+15)=[2,4,3,1,5,6]
       order(1:n,i+16)=[2,4,3,5,1,6]
       order(1:n,i+17)=[2,4,5,1,3,6]
       order(1:n,i+18)=[2,4,5,3,1,6]
       order(1:n,i+19)=[2,5,1,3,4,6]
       order(1:n,i+20)=[2,5,1,4,3,6]
       order(1:n,i+21)=[2,5,3,1,4,6]
       order(1:n,i+22)=[2,5,3,4,1,6]
       order(1:n,i+23)=[2,5,4,1,3,6]
       order(1:n,i+24)=[2,5,4,3,1,6]

       i=48
       order(1:n,i+1)=[3,1,2,4,5,6]
       order(1:n,i+2)=[3,1,2,5,4,6]
       order(1:n,i+3)=[3,1,4,2,5,6]
       order(1:n,i+4)=[3,1,4,5,2,6]
       order(1:n,i+5)=[3,1,5,2,4,6]
       order(1:n,i+6)=[3,1,5,4,2,6]
       order(1:n,i+7)=[3,2,1,4,5,6]
       order(1:n,i+8)=[3,2,1,5,4,6]
       order(1:n,i+9)=[3,2,4,1,5,6]
       order(1:n,i+10)=[3,2,4,5,1,6]
       order(1:n,i+11)=[3,2,5,1,4,6]
       order(1:n,i+12)=[3,2,5,4,1,6]
       order(1:n,i+13)=[3,4,1,2,5,6]
       order(1:n,i+14)=[3,4,1,5,2,6]
       order(1:n,i+15)=[3,4,2,1,5,6]
       order(1:n,i+16)=[3,4,2,5,1,6]
       order(1:n,i+17)=[3,4,5,1,2,6]
       order(1:n,i+18)=[3,4,5,2,1,6]
       order(1:n,i+19)=[3,5,1,2,4,6]
       order(1:n,i+20)=[3,5,1,4,2,6]
       order(1:n,i+21)=[3,5,2,1,4,6]
       order(1:n,i+22)=[3,5,2,4,1,6]
       order(1:n,i+23)=[3,5,4,1,2,6]
       order(1:n,i+24)=[3,5,4,2,1,6]

       i=72
       order(1:n,i+1)=[4,1,3,2,5,6]
       order(1:n,i+2)=[4,1,3,5,2,6]
       order(1:n,i+3)=[4,1,2,3,5,6]
       order(1:n,i+4)=[4,1,2,5,3,6]
       order(1:n,i+5)=[4,1,5,3,2,6]
       order(1:n,i+6)=[4,1,5,2,3,6]
       order(1:n,i+7)=[4,3,1,2,5,6]
       order(1:n,i+8)=[4,3,1,5,2,6]
       order(1:n,i+9)=[4,3,2,1,5,6]
       order(1:n,i+10)=[4,3,2,5,1,6]
       order(1:n,i+11)=[4,3,5,1,2,6]
       order(1:n,i+12)=[4,3,5,2,1,6]
       order(1:n,i+13)=[4,2,1,3,5,6]
       order(1:n,i+14)=[4,2,1,5,3,6]
       order(1:n,i+15)=[4,2,3,1,5,6]
       order(1:n,i+16)=[4,2,3,5,1,6]
       order(1:n,i+17)=[4,2,5,1,3,6]
       order(1:n,i+18)=[4,2,5,3,1,6]
       order(1:n,i+19)=[4,5,1,3,2,6]
       order(1:n,i+20)=[4,5,1,2,3,6]
       order(1:n,i+21)=[4,5,3,1,2,6]
       order(1:n,i+22)=[4,5,3,2,1,6]
       order(1:n,i+23)=[4,5,2,1,3,6]
       order(1:n,i+24)=[4,5,2,3,1,6]

       i=96
       order(1:n,i+1)=[5,1,3,2,4,6]
       order(1:n,i+2)=[5,1,3,4,2,6]
       order(1:n,i+3)=[5,1,2,3,4,6]
       order(1:n,i+4)=[5,1,2,4,3,6]
       order(1:n,i+5)=[5,1,4,3,2,6]
       order(1:n,i+6)=[5,1,4,2,3,6]
       order(1:n,i+7)=[5,3,1,2,4,6]
       order(1:n,i+8)=[5,3,1,4,2,6]
       order(1:n,i+9)=[5,3,2,1,4,6]
       order(1:n,i+10)=[5,3,2,4,1,6]
       order(1:n,i+11)=[5,3,4,1,2,6]
       order(1:n,i+12)=[5,3,4,2,1,6]
       order(1:n,i+13)=[5,2,1,3,4,6]
       order(1:n,i+14)=[5,2,1,4,3,6]
       order(1:n,i+15)=[5,2,3,1,4,6]
       order(1:n,i+16)=[5,2,3,4,1,6]
       order(1:n,i+17)=[5,2,4,1,3,6]
       order(1:n,i+18)=[5,2,4,3,1,6]
       order(1:n,i+19)=[5,4,1,3,2,6]
       order(1:n,i+20)=[5,4,1,2,3,6]
       order(1:n,i+21)=[5,4,3,1,2,6]
       order(1:n,i+22)=[5,4,3,2,1,6]
       order(1:n,i+23)=[5,4,2,1,3,6]
       order(1:n,i+24)=[5,4,2,3,1,6]

       iden=8*8*2*2
       iden=iden*4*3*2

     endif

  elseif (n.eq.7) then
   if (quarks) then
       part(1:n)=[1,21,21,21,21,1]
       iden=3*8*2*2
       iden=iden*3*2
   else
       part(1:n)=[21,21,21,21,21,21,21]
       order(1:n,1)=[1,2,3,4,5,6,7]
       order(1:n,2)=[1,2,3,5,4,6,7]
       order(1:n,3)=[1,2,4,3,5,6,7]
       order(1:n,4)=[1,2,4,5,3,6,7]
       order(1:n,5)=[1,2,5,3,4,6,7]
       order(1:n,6)=[1,2,5,4,3,6,7]
       order(1:n,7)=[1,3,2,4,5,6,7]
       order(1:n,8)=[1,3,2,5,4,6,7]
       order(1:n,9)=[1,3,4,2,5,6,7]

       iden=8*8*2*2
       iden=iden*5*4*3*2
  endif

  endif

  if (single_perm) then
    call cpu_time(tBefore) 
    call amps_col%init(2,n,part,order(1:n,1))
    call cpu_time(tAfter)
    t_init = tAfter-tBefore
    call cpu_time(tBefore)
    call amps_col%evaluate(n,p,5)
    call cpu_time(tAfter)
    t_eval = tAfter-tBefore
  else
    do iperm=1,nperm
     call cpu_time(tBefore)
     call amps(iperm)%init(1,n,part,order(1:n,iperm))
     call cpu_time(tAfter)
     t_init=tAfter-tBefore
     call cpu_time(tBefore)
     call amps(iperm)%evaluate(n,p,0)
     call cpu_time(tAfter)
     t_eval=tAfter-tBefore
    enddo
  endif

!  do i = 1, amps_col%nColOrd
!     write (*,*) amps_col%amps(i),amps_col%perm(1:n,i)
!  enddo
  
!  stop
  
  call Tr_allocate(n)
  do jperm=1,nperm
     do iperm=1,nperm
        if (quarks) then
          col_fac(iperm,jperm)=color_factor(iperm,jperm) ! colour factor for qqbar+gluons
          if (photon) col_fac(iperm,jperm)=color_factor_photon(iperm,jperm)
        else
          col_fac(iperm,jperm)=color_factor_gluons(iperm,jperm) ! colour factor for all-gluons
        endif
     enddo
  enddo

  amp2=0d0
  if (.not.single_perm)then
    do ih=1,2**n
     t=0d0
     do jperm=1,nperm    ! loop over permutations of conjugated amplitude
        ztemp=(0d0,0d0)
        do iperm=1,nperm ! loop over permutations of amplitude
           ztemp=ztemp+amps(iperm)%amps(ih)*col_fac(iperm,jperm)
        enddo
        t=t+dble(ztemp*dconjg(amps(jperm)%amps(ih)))
     enddo
     amp2=amp2+t
    enddo
  else
    t=0d0
    do jperm=1,nperm    ! loop over permutations of conjugated amplitude
       ztemp=(0d0,0d0)
       do iperm=1,nperm ! loop over permutations of amplitude
          ztemp=ztemp+amps_col%amps(iperm)*col_fac(iperm,jperm)
       enddo
       t=t+dble(ztemp*dconjg(amps_col%amps(jperm)))
    enddo
    !write (*,'(e20.12,x,a,x,b6.6)') t*(4*pi*alphas)**(n-2)
   amp2=t
  endif

  if (photon) then
  do i=1,n
     if (abs(part(i)).le.6) then
        if (mod(abs(part(i)),2).eq.0) Q=2d0/3d0
        if (mod(abs(part(i)),2).eq.1) Q=-1d0/3d0
     endif
  enddo
  endif


  if (photon) then 
       amp2=amp2*(4*pi*alphas)**(n-3)/dble(iden)
       amp2=amp2*(Q*dsqrt(4*pi*alphaEW))**2 
       amp2=amp2*dsqrt(2d0)*dsqrt(2d0) ! do to normalization of fund. matrices
  else
       amp2=amp2*(4*pi*alphas)**(n-2)/dble(iden)
  endif

  write(*,*) 'IDEN',iden
  write (*,*) 'Matrix element =',amp2,'GeV^',-(2*n-8)
  write(*,*) 'Init time: ',t_init
  write(*,*) 'Evaluation time: ',t_eval




contains
  double precision function color_factor(iperm,jperm)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm', for qqbar+gluons
    use color_algebra
    implicit none
    integer :: iperm,jperm
    integer, dimension(n) :: iper,jper
    real*16 :: col_factor
    iper(1:n-2)=order(2:n-1,iperm)
    jper(1:n-2)=order(2:n-1,jperm)
    Tr(0,0,0)=1 ! one term
    Tr(0,0,1)=1 ! that term is single string of matrices
    Tr(0,1,1)=2*(n-2)
    Tr(1:n-2,1,1)=iper(1:n-2) ! the order of the matrices in each term
    Tr(n-1:2*(n-2),1,1)=jper(n-2:1:-1)
    coef(1)=1
    call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
    color_factor=dble(col_factor)
    if (.not.photon) write(*,*) color_factor
  end function color_factor

  double precision function color_factor_gluons(iperm,jperm)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm', for qqbar+gluons
    use color_algebra
    implicit none
    integer :: iperm,jperm
    integer, dimension(n) :: iper,jper
    real*16 :: col_factor
    Tr(0,0,0)=1 ! one term
    Tr(0,0,1)=2 ! that term is a product of two terms
    Tr(0,1,1)=n ! both terms in the product are a trace with n matrices
    Tr(0,2,1)=n
    iper(1:n) = order(1:n,iperm)
    jper(1:n) = order(1:n,jperm)
    Tr(1:n,1,1)=iper(1:n) ! the order of the matrices in each term
    Tr(1:n,2,1)=jper(n:1:-1)
    coef(1)=1
    call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the colour string
    color_factor_gluons=dble(col_factor)
  end function color_factor_gluons

    double precision function color_factor_photon(iperm,jperm)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm', for qqbar+gluons
    use color_algebra
    implicit none
    integer :: iperm,jperm
    integer, dimension(n) :: iper,jper
    real*16 :: col_factor

    iper(1:n-3)=order(2:n-2,iperm)
    jper(1:n-3)=order(2:n-2,jperm)
    Tr(0,0,0)=1 ! one term
    Tr(0,0,1)=1 ! that term is single string of matrices
    Tr(0,1,1)=2*(n-3)
    Tr(1:n-3,1,1)=iper(1:n-3) ! the order of the matrices in each term
    Tr(n-2:2*(n-3),1,1)=jper(n-3:1:-1)
    coef(1)=1
    call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
    color_factor_photon=dble(col_factor)
  end function color_factor_photon
  
end program test_QCD
