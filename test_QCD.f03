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
  real(kind=8),dimension(:),allocatable :: mass,width
  integer,dimension(:,:),allocatable :: helmap,order
  integer :: i,ih,iden,iperm,jperm,nperm
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: amp2,t
  real(kind=8),parameter :: pi=3.14159265358979323846d0,alphas=0.118d0,alphaEW=7.547E-003
  real(kind=8),dimension(:,:),allocatable :: col_fac
  complex(kind=8) :: ztemp
  logical :: one_qq, single_perm,photon, two_qq,heavy,same_flav
  real(kind=8) :: tBefore,tAfter,t_eval,t_init
  real(kind=8) :: Q
  integer :: gi,gj,ui,uj
  integer :: it ! quark order
  real(kind=8) :: top_mass

  n=6
  one_qq=.false.
  photon=.false.

  two_qq=.true.
  heavy=.false.
  
  same_flav=.false.

  single_perm = .false.
  
  allocate(part(n))
  allocate(mass(n))
  allocate(width(n))

  mass(1:n)=0d0
  width(1:n)=0d0

  if (n.eq.4) then
     if (one_qq) then
        nperm = 2
        if (photon) nperm = 1 
     else
        nperm = 3*2
     endif
     if (single_perm) nperm=1
  elseif (n.eq.5) then
     if (one_qq) then
        nperm = 3*2   
        if (photon) nperm = 2
     else
        nperm = 4*3*2
     endif
     if (single_perm) nperm=1
  elseif (n.eq.6) then
     if (one_qq) then
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

  allocate(p(0:3,n))
  
  if (n.eq.4) then
     p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
     p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03/)
     p(0:3,3)=(/0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03/)
     p(0:3,4)=(/0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03/)
     if (heavy) then
       p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
       p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03/)
       p(0:3,3)=(/0.5000000E+03,  0.2252977E+03,  0.3101302E+03, -0.2704278E+03/)
       p(0:3,4)=(/0.5000000E+03, -0.2252977E+03, -0.3101302E+03,  0.2704278E+03/)
     endif
  elseif (n.eq.5) then
     p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
     p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  -0.5000000E+03/)
     p(0:3,3)=(/0.4959179E+03, -0.1999048E+02,  0.7981352E+02, -0.4890448E+03/)
     p(0:3,4)=(/0.1328543E+03, -0.2648250E+02, -0.4416981E+02,  0.1224662E+03/)
     p(0:3,5)=(/0.3712277E+03,  0.4647298E+02, -0.3564372E+02,  0.3665785E+03/)
     if (heavy) then
       p(0:3,1)=(/500.0000d0,  0.0000000d0,  0.0000000d0,  500.0000d0/)
       p(0:3,2)=(/500.0000d0,  0.0000000d0,  0.0000000d0, -500.0000d0/)
       p(0:3,3)=(/457.7589d0,  156.6052d0,   350.8682d0, -178.8310d0/)
       p(0:3,4)=(/378.3336d0, -16.94009d0, -321.3412d0,  98.28614d0/)
       p(0:3,5)=(/163.9075d0, -139.6651d0, -29.52695d0,  80.54490d0/)
     endif
  elseif (n.eq.6) then
     p(0:3,1)=(/500.0000000d0,   0.00000000d0,  0.0000000000d0,   500.0000000d0/)
     p(0:3,2)=(/500.0000000d0,   0.00000000d0,  0.0000000000d0,  -500.0000000d0/)
     p(0:3,3)=(/88.55133305d0,  -22.1006902d0,  40.080353191d0,  -75.80543095d0/)
     p(0:3,4)=(/328.3294192d0,  -103.849611d0, -301.93375538d0,   76.49492138d0/)
     p(0:3,5)=(/152.3581094d0,  -105.880959d0, -97.709638326d0,   49.54838522d0/)
     p(0:3,6)=(/430.7611382d0,   231.831261d0,  359.56304052d0,  -50.23787565d0/)
     if (heavy) then
       p(0:3,1)=(/500.0000000d0,   0.00000000d0,  0.0000000000d0,   500.0000000d0/)
       p(0:3,2)=(/500.0000000d0,   0.00000000d0,  0.0000000000d0,  -500.0000000d0/)
       p(0:3,3)=(/0.1881691E+03, -0.1847334E+02,  0.3350204E+02, -0.6336362E+02/)
       p(0:3,4)=(/0.3244180E+03, -0.8680496E+02, -0.2523779E+03,  0.6393995E+02/)
       p(0:3,5)=(/0.1273518E+03, -0.8850291E+02, -0.8167273E+02,  0.4141610E+02/)
       p(0:3,6)=(/0.3600611E+03,  0.1937812E+03,  0.3005486E+03, -0.4199242E+02/)
     endif
  endif

  if (n.eq.4) then

     if (one_qq) then
         part(1:n)=[-1,1,21,21]
         iden=3*3*2*2
         order(1:n,1)= [1,3,4,2]
         order(1:n,2)= [1,4,3,2]
         iden=iden*2
         if (heavy) then
            nperm=2
            part(1:n)=[21,21,6,-6]
            order(1:n,1)= [3,1,2,4]
            order(1:n,2)= [3,2,1,4]
            iden=8*8*2*2
         endif
     endif

     if (two_qq) then
         nperm=2
         part(1:n)=[-1,1,-2,2] !  1 -1 -2 2
         iden=3*3*2*2
         order(1:n,1)=[1,3,4,2]
         order(1:n,2)=[1,2,4,3]
         if (same_flav) then
            nperm=1
            part(1:n)=[-1,1,-1,1] !  1 -1 1 -1
            iden=3*3*2*2
            order(1:n,1)=[1,3,4,2]
            order(1:n,2)=[1,2,4,3]
         endif
         if (heavy) then
            nperm=2
            part(1:n)=[-1,1,6,-6] !  1 -1 1 -1
            iden=3*3*2*2
            order(1:n,1)=[1,4,3,2]
            order(1:n,2)=[1,2,3,4]

            top_mass=173d0

            mass(1:n)=[0d0,0d0,top_mass,top_mass]
            width(1:n)=[0d0,0d0,1.4915d0,1.4915d0]
         endif

     endif

     

     if (photon) then
       part(1:n)=[-1,21,-1,22]
       iden=3*8*2*2
       order(1:n,1)= [1,2,4,3]
       if (.not. single_perm) then
         iden=iden*1
       endif
     endif

  elseif (n.eq.5) then

     if (two_qq) then
         nperm=4
         part(1:n)=[-1,1,2,-2,21] !    -1 1 2 -2 21
         iden=3*3*2*2
         order(1:n,1)=[1,5,2,3,4]
         order(1:n,2)=[1,2,3,5,4]
         order(1:n,3)=[1,5,4,3,2]
         order(1:n,4)=[1,4,3,5,2]
         if (same_flav) then
            nperm=4
            part(1:n)=[-1,1,1,-1,21] !    1 -1 1 -1 21
            iden=3*3*2*2
            order(1:n,1)=[1,5,2,3,4]
            order(1:n,2)=[1,2,3,5,4]
            order(1:n,3)=[1,5,4,3,2]
            order(1:n,4)=[1,4,3,5,2]
         endif
         if (heavy) then
           nperm=4
           part(1:n)=[-1,1,6,-6,21] !    -1 1 2 -2 21
           iden=3*3*2*2
           order(1:n,4)=[1,5,2,3,4]
           order(1:n,3)=[1,2,3,5,4]
           order(1:n,2)=[1,5,4,3,2]
           order(1:n,1)=[1,4,3,5,2]

           top_mass=173.d0

           mass(1:n)=[0d0,0d0,top_mass,top_mass,0d0]
           width(1:n)=[0d0,0d0,1.4915d0,1.4915d0,0d0]
         endif
     endif

     if (photon) then
       part(1:n)=[-1,1,21,21,22]
       iden=3*3*2*2
       order(1:n,1)= [5,4,3,1,2]
       if (.not. single_perm) then
         order(1:n,2)= [1,4,3,5,2]
         iden=iden*2
       endif
     endif

    if (one_qq) then
       !nperm=1
       part(1:n)=[-1,1,21,21,21]
       iden=3*3*2*2 ! initial status colours and helicities/polarisation
       ! 1
       !order(1:n,1)=  [1,3,4,5,2]
       !if (.not.single_perm) then
       !  order(1:n,2)=[1,3,5,4,2]
       !  order(1:n,3)=[1,4,3,5,2]
       !  order(1:n,4)=[1,4,5,3,2]
       !  order(1:n,5)=[1,5,3,4,2]
       !  order(1:n,6)=[1,5,4,3,2]
       ! 2
       !order(1:n,1)=  [5,4,3,1,2]
       !if (.not.single_perm) then
       !  order(1:n,2)=[5,3,4,1,2]
       !  order(1:n,3)=[4,3,5,1,2]
       !  order(1:n,4)=[3,4,5,1,2]
       !  order(1:n,5)=[3,5,4,1,2]
       !  order(1:n,6)=[4,5,3,1,2]
       ! 3
       !order(1:n,1)=  [4,3,2,1,5]
       !if (.not.single_perm) then
       !  order(1:n,2)=[3,4,2,1,5]
       !  order(1:n,3)=[3,5,2,1,4]
       !  order(1:n,4)=[4,5,2,1,3]
       !  order(1:n,5)=[5,4,2,1,3]
       !  order(1:n,6)=[5,3,2,1,4]
       ! 4
       !order(1:n,1)=  [3,1,2,5,4]
       !if (.not.single_perm) then
       !  order(1:n,2)=[3,1,2,4,5]
       !  order(1:n,3)=[4,1,2,5,3]
       !  order(1:n,4)=[4,1,2,3,5]
       !  order(1:n,5)=[5,1,2,4,3]
       !  order(1:n,6)=[5,1,2,3,4]
       ! 5
       order(1:n,1)=  [2,1,3,4,5]
       if (.not.single_perm) then
         order(1:n,2)=[2,1,3,5,4]
         order(1:n,3)=[2,1,4,3,5]
         order(1:n,4)=[2,1,4,5,3]
         order(1:n,5)=[2,1,5,4,3]
         order(1:n,6)=[2,1,5,3,4]
         iden=iden*2*3  ! final state identical particles
       endif
       if (heavy) then
         nperm=6
         part(1:n)=[21,21,6,-6,21]
         order(1:n,1)= [3,1,2,5,4]
         order(1:n,2)= [3,1,5,2,4]
         order(1:n,3)= [3,2,1,5,4]
         order(1:n,4)= [3,2,5,1,4]
         order(1:n,5)= [3,5,1,2,4]
         order(1:n,6)= [3,5,2,1,4]
         iden=8*8*2*2

         mass(1:n)=[0d0,0d0,173d0,173d0,0d0]
         width(1:n)=[0d0,0d0,1.4915d0,1.4915d0,0d0]

       endif
     endif

     if ((.not.photon).and.(.not.one_qq).and.(.not.two_qq)) then
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
     if (one_qq) then
       part(1:n)=[-1,1,21,21,21,21]
       iden=3*3*2*2 ! initial status colours and helicities/polarisation
         order(1:n,1)=  [1,3,4,5,6,2]
       if (.not.single_perm) then
         order(1:n,2)=  [1,3,4,6,5,2]
         order(1:n,3)=  [1,3,5,4,6,2]
         order(1:n,4)=  [1,3,5,6,4,2]
         order(1:n,5)=  [1,3,6,4,5,2]
         order(1:n,6)=  [1,3,6,5,4,2]
         order(1:n,7)=  [1,4,3,5,6,2]
         order(1:n,8)=  [1,4,3,6,5,2]
         order(1:n,9)=  [1,4,5,3,6,2]
         order(1:n,10)= [1,4,5,6,3,2]
         order(1:n,11)= [1,4,6,3,5,2]
         order(1:n,12)= [1,4,6,5,3,2]
         order(1:n,13)= [1,5,3,4,6,2]
         order(1:n,14)= [1,5,3,6,4,2]
         order(1:n,15)= [1,5,4,3,6,2]
         order(1:n,16)= [1,5,4,6,3,2]
         order(1:n,17)= [1,5,6,3,4,2]
         order(1:n,18)= [1,5,6,4,3,2]
         order(1:n,19)= [1,6,3,5,4,2]
         order(1:n,20)= [1,6,3,4,5,2]
         order(1:n,21)= [1,6,5,3,4,2]
         order(1:n,22)= [1,6,5,4,3,2]
         order(1:n,23)= [1,6,4,3,5,2]
         order(1:n,24)= [1,6,4,5,3,2]

         iden=iden*2*3*4  ! final state identical particles
       endif


       if (heavy) then
         nperm=24
         part(1:n)=[21,21,6,-6,21,21]
         order(1:n,1)= [3,1,2,5,6,4]
         order(1:n,2)= [3,1,5,2,6,4]
         order(1:n,3)= [3,2,1,5,6,4]
         order(1:n,4)= [3,2,5,1,6,4]
         order(1:n,5)= [3,5,1,2,6,4]
         order(1:n,6)= [3,5,2,1,6,4]

         order(1:n,7)= [3,1,2,6,5,4]
         order(1:n,8)= [3,1,5,6,2,4]
         order(1:n,9)= [3,2,1,6,5,4]
         order(1:n,10)= [3,2,5,6,1,4]
         order(1:n,11)= [3,5,1,6,2,4]
         order(1:n,12)= [3,5,2,6,1,4]

         order(1:n,13)= [3,1,6,2,5,4]
         order(1:n,14)= [3,1,6,5,2,4]
         order(1:n,15)= [3,2,6,1,5,4]
         order(1:n,16)= [3,2,6,5,1,4]
         order(1:n,17)= [3,5,6,1,2,4]
         order(1:n,18)= [3,5,6,2,1,4]

         order(1:n,19)= [3,6,1,2,5,4]
         order(1:n,20)= [3,6,1,5,2,4]
         order(1:n,21)= [3,6,2,1,5,4]
         order(1:n,22)= [3,6,2,5,1,4]
         order(1:n,23)= [3,6,5,1,2,4]
         order(1:n,24)= [3,6,5,2,1,4]
         iden=8*8*2*2
         iden=iden*2

         mass(1:n)=[0d0,0d0,173d0,173d0,0d0,0d0]
         width(1:n)=[0d0,0d0,1.4915d0,1.4915d0,0d0,0d0]

       endif

     endif

     if (two_qq) then
         nperm=12
         part(1:n)=[-1,1,2,-2,21,21] !    -1 1 2 -2 21
         iden=3*3*2*2
         order(1:n,1)=[1,5,6,2,3,4]
         order(1:n,2)=[1,6,5,2,3,4]
         order(1:n,3)=[1,2,3,5,6,4]
         order(1:n,4)=[1,2,3,6,5,4]
         order(1:n,5)=[1,5,2,3,6,4]
         order(1:n,6)=[1,6,2,3,5,4]

         order(1:n,7)=[1,5,6,4,3,2]
         order(1:n,8)=[1,6,5,4,3,2]
         order(1:n,9)=[1,4,3,5,6,2]
         order(1:n,10)=[1,4,3,6,5,2]
         order(1:n,11)=[1,5,4,3,6,2]
         order(1:n,12)=[1,6,4,3,5,2]
         iden=iden*2

         if (heavy) then
         nperm=12
         part(1:n)=[-1,1,6,-6,21,21] !    -1 1 2 -2 21
         iden=3*3*2*2
         order(1:n,1)=[1,5,6,2,3,4]
         order(1:n,2)=[1,6,5,2,3,4]
         order(1:n,3)=[1,2,3,5,6,4]
         order(1:n,4)=[1,2,3,6,5,4]
         order(1:n,5)=[1,5,2,3,6,4]
         order(1:n,6)=[1,6,2,3,5,4]

         order(1:n,7)=[1,5,6,4,3,2]
         order(1:n,8)=[1,6,5,4,3,2]
         order(1:n,9)=[1,4,3,5,6,2]
         order(1:n,10)=[1,4,3,6,5,2]
         order(1:n,11)=[1,5,4,3,6,2]
         order(1:n,12)=[1,6,4,3,5,2]
         iden=iden*2

         top_mass=173.d0

         mass(1:n)=[0d0,0d0,top_mass,top_mass,0d0,0d0]
         width(1:n)=[0d0,0d0,1.4915d0,1.4915d0,0d0,0d0]
         endif
     endif



  endif

  it = 0 ! dummy
  if (single_perm) then
    call cpu_time(tBefore) 
    call amps_col%init(2,n,part,mass,width,order(1:n,1),it)
    call cpu_time(tAfter)
    t_init = tAfter-tBefore
    call cpu_time(tBefore)
    call amps_col%evaluate(n,p,mass,width,5,part)
    call cpu_time(tAfter)
    t_eval = tAfter-tBefore
  else
    do iperm=1,nperm
     call cpu_time(tBefore)
     call amps(iperm)%init(1,n,part,mass,width,order(1:n,iperm),it)
     call cpu_time(tAfter)
     t_init=tAfter-tBefore
     call cpu_time(tBefore)
     call amps(iperm)%evaluate(n,p,mass,width,0,part)
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
       if (one_qq) col_fac(iperm,jperm)=color_factor_one_qq(iperm,jperm) ! colour factor for qqbar+gluons
       if (photon) col_fac(iperm,jperm)=color_factor_photon(iperm,jperm)
       if (two_qq) then
             call get_perm_params(iperm,jperm,gi,gj,ui,uj)
             col_fac(iperm,jperm)=color_factor_two_qq(iperm,gi,ui,jperm,gj,uj)
       endif
       if ((.not.photon).and.(.not.one_qq).and.(.not.two_qq)) then
          col_fac(iperm,jperm)=color_factor_gluons(iperm,jperm) ! colour factor for all-gluons
       endif
     enddo
  enddo


  !do iperm=1,nperm
  !    write(*,*) 'iperm',iperm
  !    write(*,*) amps(iperm)%amps(10)
  !    do jperm=1,nperm
  !       write(*,*) 'jperm',jperm
  !       write(*,*) 'color ',col_fac(iperm,jperm)
  !    enddo
  !enddo

  amp2=0d0
  if (.not.single_perm)then
    do ih=1,2**n
     t=0d0
     !write(*,*) '**********************'
     !write(*,*) 'ih',ih
     !do i=1,n
     !   write(*,*) 'particle number',i
     !   if (btest(ih-1,i-1)) then
     !      write(*,*) 'hel +1'
     !   else
     !      write(*,*) 'hel -1'
     !   endif
     !enddo
     do jperm=1,nperm    ! loop over permutations of conjugated amplitude
        ztemp=(0d0,0d0)
        do iperm=1,nperm ! loop over permutations of amplitude
           !write(*,*) 'added',iperm,jperm
           !write(*,*) 'helmap: ',amps(iperm)%helmap(ih)
           ztemp=ztemp+amps(iperm)%amps(ih) *col_fac(iperm,jperm)*(4*pi*alphas)**(n-2)
      !     write(*,*) amps(1)%amps(ih)
      !     write(*,*) amps(2)%amps(ih)
      !     write(*,*) '****'
           !write(*,*) 'color',col_fac(iperm,jperm)
           !write(*,*) 'amp for iperm',iperm
        !write(*,*) amps(iperm)%amps(ih)*(4*pi*alphas)**(n-2)*col_fac(iperm,jperm)*dconjg(amps(jperm)%amps(ih))
        enddo
        !write(*,*) 'total ztemp',ztemp
        t=t+dble(ztemp*dconjg(amps(jperm)%amps(ih)))
        !write(*,*) 'closing with amp',dconjg(amps(jperm)%amps(ih))
     enddo
     amp2=amp2+t
     !write(*,*) 'T',T
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
       amp2=amp2/dble(iden) !*(4*pi*alphas)**(n-2)/dble(iden)
  endif

  write(*,*) 'IDEN',iden
  write (*,*) 'Matrix element =',amp2,'GeV^',-(2*n-8)
  write(*,*) 'Init time: ',t_init
  write(*,*) 'Evaluation time: ',t_eval



contains
  double precision function color_factor_one_qq(iperm,jperm)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm', for qqbar+gluons
    use color_algebra
    implicit none
    integer :: iperm,jperm
    integer, dimension(n) :: iper,jper
    real*16 :: col_factor
    ! 1
    iper(1:n-2)=order(2:n-1,iperm)
    jper(1:n-2)=order(2:n-1,jperm)

    ! 2
    !iper(1:n-2)=order(1:n-2,iperm)
    !jper(1:n-2)=order(1:n-2,jperm)

    ! 3
    !iper(1:1)=order(n:n,iperm)
    !jper(1:1)=order(n:n,jperm)
    !iper(2:n-2)=order(1:n-3,iperm)
    !jper(2:n-2)=order(1:n-3,jperm)

    ! 4
    !iper(1:2)=order(n-1:n,iperm)
    !jper(1:2)=order(n-1:n,jperm)
    !iper(3:n-2)=order(1:n-4,iperm)
    !jper(3:n-2)=order(1:n-4,jperm)

    !5
    !iper(1:n-2)=order(3:n,iperm)
    !jper(1:n-2)=order(3:n,jperm)


    Tr(0,0,0)=1 ! one term
    Tr(0,0,1)=1 ! that term is single string of matrices
    Tr(0,1,1)=2*(n-2)
    Tr(1:n-2,1,1)=iper(1:n-2) ! the order of the matrices in each term
    Tr(n-1:2*(n-2),1,1)=jper(n-2:1:-1)
    coef(1)=1
    call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
    color_factor_one_qq=dble(col_factor)

    !if (iperm.ne.jperm) color_factor_one_qq=0d0
  end function color_factor_one_qq

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

  double precision function color_factor_two_qq(iperm,gi,ui,jperm,gj,uj)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm', for qqbar+gluons
    use color_algebra
    implicit none
    integer :: iperm,jperm,gi,gj,ui,uj
    integer, dimension(n-4) :: iper,jper
    real*16 :: col_factor
    integer :: i,q_i,q_j
    integer,dimension(n) :: temp_part
    integer,dimension(2*n) :: dorder_i,dorder_j

    dorder_i(1:n)=order(1:n,iperm)
    dorder_i(n+1:2*n)=order(1:n,iperm)
    dorder_j(1:n)=order(1:n,jperm)
    dorder_j(n+1:2*n)=order(1:n,jperm)

    temp_part = part
    do i=1,n
       if ((i.le.2).and.(abs(temp_part(i)).le.6)) then
               temp_part(i) = -temp_part(i)
        endif
    enddo

    do i=1,n
       if ((temp_part(dorder_i(i)).le.6.and.temp_part(dorder_i(i)).ge.0)) then ! found first quark
          q_i = i
          exit
       endif
    enddo
    do i=1,n
       if ((temp_part(dorder_j(i)).le.6.and.temp_part(dorder_j(i)).ge.0)) then ! found first quark
          q_j = i
          exit
       endif
    enddo

    iper(1:gi)    =dorder_i(q_i+1:q_i+gi)
    iper(gi+1:n-4)=dorder_i(q_i+gi+3:q_i+gi+2+(n-4-gi))
    jper(1:gj)    =dorder_j(q_j+1:q_j+gj)
    jper(gj+1:n-4)=dorder_j(q_j+gj+3:q_j+gj+2+(n-4-gj))

    if (ui.eq.uj.and.ui.eq.1) then
       Tr(0,0,0)=1 ! one term
       Tr(0,0,1)=2 
       Tr(0,1,1) = gi+gj  ! number of generators in first trace
       Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
       Tr(1:gi,1,1) = iper(1:gi)
       Tr(gi+1:gi+gj,1,1) = jper(gj:1:-1)
       Tr(1:n-4-gi,2,1) = iper(gi+1:n-4)
       Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-4:gj+1:-1)
       coef(1)=(1d0,0d0)
       call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
       color_factor_two_qq=dble(col_factor)
    elseif (ui.eq.uj.and.ui.eq.2) then
       Tr(0,0,0) = 1 ! one term
       Tr(0,0,1) = 2 ! product of two traces
       Tr(0,1,1) = gi+gj  ! number of generators in first trace
       Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
       Tr(1:gi,1,1) = iper(1:gi)
       Tr(gi+1:gi+gj,1,1) = jper(gj:1:-1)
       Tr(1:n-4-gi,2,1) = iper(gi+1:n-4)
       Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-4:gj+1:-1)
       coef(1)=(1d0,0d0) 
       if (.not.same_flav) coef(1)=coef(1)*((1d0/3d0)**2)*(1d0,0d0)
       call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
       color_factor_two_qq=dble(col_factor)

    elseif (ui.eq.2.and.uj.eq.1) then
       Tr(0,0,0)=1 ! one term
       Tr(0,0,1)=1 ! a single trace
       Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
       Tr(1:gi,1,1) = iper(1:gi)
       Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-4:gj+1:-1)
       Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(gi+1:n-4)
       Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(gj:1:-1)
       coef(1)=(1d0,0d0)
       if (.not.same_flav) coef(1)=((1d0/3d0))*(1d0,0d0)
       coef(1)=coef(1)*(-1d0)
       call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
       color_factor_two_qq=dble(col_factor)
    elseif (ui.eq.1.and.uj.eq.2) then
       Tr(0,0,0)=1 ! one term
       Tr(0,0,1)=1 ! a single trace
       Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
       Tr(1:gi,1,1) = iper(1:gi)
       Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-4:gj+1:-1)
       Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(gi+1:n-4)
       Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(gj:1:-1)
       coef(1)=(1d0,0d0)
       if (.not.same_flav) coef(1)=((1d0/3d0))*(1d0,0d0)
       coef(1)=coef(1)*(-1d0)
       call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
       color_factor_two_qq=dble(col_factor)
    endif

    if (same_flav) color_factor_two_qq = ((1d0+1d0/3d0)**2)*color_factor_two_qq

  end function color_factor_two_qq
  
  subroutine get_perm_params(iperm,jperm,gi,gj,ui,uj)
    implicit none
    integer :: iperm,jperm,ui,uj,gi,gj
    integer :: i,q_i,q_j
    integer,dimension(n) :: temp_part

    temp_part = part
    do i=1,n
       if ((i.le.2).and.(abs(temp_part(i)).le.6)) then
               temp_part(i) = -temp_part(i)
        endif
    enddo

    do i=1,n
       if ((temp_part(order(i,iperm)).le.6.and.temp_part(order(i,iperm)).ge.0)) then ! found first quark
          q_i = i
          exit
       endif
    enddo

    do i=q_i+1,n
       if (temp_part(order(i,iperm)).le.-1) then ! found first a-quark after q_i
          if (abs(temp_part(order(i,iperm))).eq. abs(temp_part(order(q_i,iperm)))) then
                  ui = 2
          elseif (abs(temp_part(order(i,iperm))).ne. abs(temp_part(order(q_i,iperm)))) then
                  ui = 1
          endif
          gi=i-q_i-1
          exit
       endif
    enddo

    if (amps(1)%same_flav) then
      if (order(n,iperm).eq.order(n,1)) ui = 1
      if (order(n,iperm).ne.order(n,1)) ui = 2
    endif

    do i=1,n
       if ((temp_part(order(i,jperm)).le.6.and.temp_part(order(i,jperm)).ge.0)) then ! found first quark
          q_j = i
          exit
       endif
    enddo

    do i=q_j+1,n
       if (temp_part(order(i,jperm)).le.-1) then ! found first a-quark after q_i
          if (abs(temp_part(order(i,jperm))).eq. abs(temp_part(order(q_j,jperm)))) then
                  uj = 2
          elseif (abs(temp_part(order(i,jperm))).ne. abs(temp_part(order(q_j,jperm)))) then
                  uj = 1
          endif
          gj=i-q_j-1
          exit
       endif
    enddo

    if (amps(1)%same_flav) then
      if (order(n,jperm).eq.order(n,1)) uj = 1
      if (order(n,jperm).ne.order(n,1)) uj = 2
    endif

  end subroutine get_perm_params

end program test_QCD
