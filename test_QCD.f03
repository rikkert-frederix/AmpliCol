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
  logical :: one_qq, single_perm,photon, two_qq
  real(kind=8) :: tBefore,tAfter,t_eval,t_init
  real(kind=8) :: Q
  integer :: gi,gj,ui,uj

  n=5
  one_qq=.false.
  photon=.false.
  two_qq=.true.

  single_perm = .false.
  
  allocate(part(n))

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
  elseif (n.eq.5) then
     p(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
     p(0:3,2)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  -0.5000000E+03/)
     p(0:3,4)=(/0.4959179E+03, -0.1999048E+02,  0.7981352E+02, -0.4890448E+03/)
     p(0:3,3)=(/0.1328543E+03, -0.2648250E+02, -0.4416981E+02,  0.1224662E+03/)
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

     if (one_qq) then
         part(1:n)=[-1,1,21,21]
         iden=3*3*2*2
         order(1:n,1)= [1,3,4,2]
         if (.not. single_perm) then
           order(1:n,2)= [1,4,3,2]
           iden=iden*2
         endif
     endif

     if (two_qq) then
         nperm=2
         part(1:n)=[-1,1,-2,2]
         iden=3*3*2*2
         order(1:n,1)=[1,3,4,2]
         order(1:n,2)=[1,2,4,3]
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
         part(1:n)=[1,-1,-2,2,21]
         iden=3*3*2*2
         order(1:n,1)=[2,5,3,4,1]
         order(1:n,2)=[2,3,4,5,1]
         order(1:n,3)=[2,1,4,5,3]
         order(1:n,4)=[2,5,1,4,3]
     endif

     if (photon) then
       part(1:n)=[-1,1,21,21,22]
       iden=3*3*2*2
       order(1:n,1)= [1,3,4,5,2]
       if (.not. single_perm) then
         order(1:n,2)= [1,4,3,5,2]
         iden=iden*2
       endif
     endif

    if (one_qq) then
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

  amp2=0d0
  if (.not.single_perm)then
    do ih=1,2**n
     t=0d0
     do jperm=1,nperm    ! loop over permutations of conjugated amplitude
        ztemp=(0d0,0d0)
        do iperm=1,nperm ! loop over permutations of amplitude
           write(*,*) 'added',iperm,jperm
           ztemp=ztemp+amps(iperm)%amps(ih)*col_fac(iperm,jperm)
           write(*,*) 'amp',amps(iperm)%amps(ih)
           write(*,*) 'ih',ih
           write(*,*) col_fac(iperm,jperm)
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
  double precision function color_factor_one_qq(iperm,jperm)
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
    color_factor_one_qq=dble(col_factor)
    if (.not.photon) write(*,*) color_factor_one_qq
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

    iper(1:gi)    =order(2:1+gi,iperm)
    iper(gi+1:n-4)=order(4+gi:n-1,iperm)
    jper(1:gj)   =order(2:1+gj,jperm)
    jper(gj+1:n-4)=order(4+gj:n-1,jperm)
    if (ui.eq.uj.and.ui.eq.1) then
       Tr(0,0,0)=1 ! one term
       Tr(0,0,1)=2 
       Tr(0,1,1) = gi+gj  ! number of generators in first trace
       Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
       Tr(1:gi,1,1) = iper(1:gi)
       Tr(gi+1:gi+gj,1,1) = jper(gj:1:-1)
       Tr(1:n-4-gi,2,1) = iper(gi+1:n-4)
       Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-4:gj+1:-1)
       coef(1)=(1d0)
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
       coef(1)=((-1/3d0)**2)*(1d0)
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
       coef(1)=((-1/3d0))*(1d0,0d0)
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
       coef(1)=((-1/3d0))*(1d0,0d0)
       call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
       color_factor_two_qq=dble(col_factor)
    endif
  end function color_factor_two_qq
  
  subroutine get_perm_params(iperm,jperm,gi,gj,ui,uj)
    implicit none
    integer :: iperm,jperm,ui,uj,gi,gj
    integer :: i
    integer,dimension(n) :: temp_part

    if (abs(part(order(1,iperm))).eq.abs(part(order(n,iperm)))) then 
       ui = 1
    endif
    if (abs(part(order(1,iperm))).ne.abs(part(order(n,iperm)))) then
            ui = 2
    endif
    if (abs(part(order(1,jperm))).eq.abs(part(order(n,jperm)))) uj = 1
    if (abs(part(order(1,jperm))).ne.abs(part(order(n,jperm)))) uj = 2

    temp_part = part
    do i=1,n
       if ((i.le.2).and.(abs(temp_part(i)).le.6)) then
               temp_part(i) = -temp_part(i)
        endif
    enddo
    do i=1,n
    if (temp_part(order(i,iperm)).lt.0) then
        gi = i-2
        exit
    endif
    enddo
    do i=1,n
    if (temp_part(order(i,jperm)).lt.0) then
        gj = i-2
        exit
    endif
    enddo
  end subroutine get_perm_params

end program test_QCD
