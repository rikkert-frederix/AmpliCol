! gfortran -fbounds-check -o test_QCD math_functions.f03 color_algebra.f95 feynmanrules.f03 amplitude_QCD.f03 test_QCD.f03 

program test_QCD
  use amplitude_QCD_mod
  use color_algebra
  implicit none
  type(amplitude_QCD),dimension(:),allocatable :: amps
  type(amplitude_QCD) :: amps_col
  integer :: n
  integer,dimension(:),allocatable :: part
  integer,dimension(:,:),allocatable :: helmap,order
  integer :: i,ih,iden,iperm,jperm,nperm
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: amp2,t
  real(kind=8),parameter :: pi=3.14159265358979323846d0,alphas=0.118d0
  real(kind=8),dimension(:,:),allocatable :: col_fac
  complex(kind=8) :: ztemp
  n=6
  allocate(part(n))

  if (n.eq.4) then
     nperm=1
  elseif (n.eq.5) then
     nperm=6
!!$     nperm=1
  elseif (n.eq.6) then
     nperm=24
!!$     nperm=1
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
  
  if (n.eq.5) then
     p(0:3,1)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,   500.00000000000000d0/)
     p(0:3,2)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,  -500.00000000000000d0/)
     p(0:3,3)=(/362.35913008614813d0,   211.67753703547245d0,  -165.00336466515242d0,   243.45564097092699d0/)
     p(0:3,4)=(/248.22961274839184d0,   68.500483807498668d0,   237.98757697004567d0,   16.956932838275925d0/)
     p(0:3,5)=(/389.41125716545980d0,  -280.17802084297108d0,  -72.984212304893376d0,  -260.41257380920291d0/)
  elseif (n.eq.6) then
     p(0:3,1)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,   500.00000000000000d0/)
     p(0:3,2)=(/500.00000000000000d0,   0.0000000000000000d0,   0.0000000000000000d0,  -500.00000000000000d0/)
     p(0:3,3)=(/88.551333054502976d0,  -22.100690287689982d0,   40.080353191685326d0,  -75.805430956936632d0/)
     p(0:3,4)=(/328.32941922709853d0,  -103.84961188345629d0,  -301.93375538954007d0,   76.494921387165888d0/)
     p(0:3,5)=(/152.35810946743064d0,  -105.88095966659219d0,  -97.709638326975707d0,   49.548385226792817d0/)
     p(0:3,6)=(/430.76113825096763d0,   231.83126183773845d0,   359.56304052483051d0,  -50.237875657022109d0/)
  endif

  if (n.eq.5) then
!!$  part(1:n)=[-1,21,21,21,-1]
!!$  order(1:n,1)=[1,2,3,4,5]
!!$  order(1:n,2)=[1,2,4,3,5]
!!$  order(1:n,3)=[1,3,2,4,5]
!!$  order(1:n,4)=[1,3,4,2,5]
!!$  order(1:n,5)=[1,4,2,3,5]
!!$  order(1:n,6)=[1,4,3,2,5]
!!$  iden=3*8*2*2 ! initial status colours and helicities/polarisations
!!$  iden=iden*2  ! final state identical particles

!!$  part(1:n)=[21,21,1,-1,21]
!!$  order(1:n,1)=[3,1,2,5,4]
!!$  order(1:n,2)=[3,1,5,2,4]
!!$  order(1:n,3)=[3,2,1,5,4]
!!$  order(1:n,4)=[3,2,5,1,4]
!!$  order(1:n,5)=[3,5,1,2,4]
!!$  order(1:n,6)=[3,5,2,1,4]
!!$  iden=8*8*2*2 ! initial status colours and helicities/polarisations

     part(1:n)=[1,-1,21,21,21]
     order(1:n,1)=[2,3,4,5,1]
     order(1:n,2)=[2,3,5,4,1]
     order(1:n,3)=[2,4,3,5,1]
     order(1:n,4)=[2,4,5,3,1]
     order(1:n,5)=[2,5,3,4,1]
     order(1:n,6)=[2,5,4,3,1]
     iden=3*3*2*2
     iden=iden*6

!!$     part(1:n)=[21,21,21,21,21]
!!$     order(1:n,1)=[1,2,3,4,5]
!!$     iden=1
     
  elseif (n.eq.6) then
     part(1:n)=[-1,21,21,21,21,-1]
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
     iden=3*8*2*2
     iden=iden*6

!!$     part(1:n)=[21,21,-1,21,1,21]
!!$     order(1:n,1)=[5,2,6,4,1,3]
!!$     order(1:n,2)=[5,2,6,1,4,3]
!!$     order(1:n,3)=[5,2,4,6,1,3]
!!$     order(1:n,4)=[5,2,4,1,6,3]
!!$     order(1:n,5)=[5,2,1,6,4,3]
!!$     order(1:n,6)=[5,2,1,4,6,3]
!!$     order(1:n,7)=[5,6,2,4,1,3]
!!$     order(1:n,8)=[5,6,2,1,4,3]
!!$     order(1:n,9)=[5,6,4,2,1,3]
!!$     order(1:n,10)=[5,6,4,1,2,3]
!!$     order(1:n,11)=[5,6,1,2,4,3]
!!$     order(1:n,12)=[5,6,1,4,2,3]
!!$     order(1:n,13)=[5,4,2,6,1,3]
!!$     order(1:n,14)=[5,4,2,1,6,3]
!!$     order(1:n,15)=[5,4,6,2,1,3]
!!$     order(1:n,16)=[5,4,6,1,2,3]
!!$     order(1:n,17)=[5,4,1,2,6,3]
!!$     order(1:n,18)=[5,4,1,6,2,3]
!!$     order(1:n,19)=[5,1,2,6,4,3]
!!$     order(1:n,20)=[5,1,2,4,6,3]
!!$     order(1:n,21)=[5,1,6,2,4,3]
!!$     order(1:n,22)=[5,1,6,4,2,3]
!!$     order(1:n,23)=[5,1,4,2,6,3]
!!$     order(1:n,24)=[5,1,4,6,2,3]
!!$     iden=8*8*2*2
!!$     iden=iden*2


!!$     part(1:n)=[21,21,21,21,21,21]
!!$     order(1:n,1)=[1,5,3,4,2,6]
!!$     iden=1

  endif
  
!!$  do iperm=1,nperm
!!$     call amps(iperm)%init(1,n,part,order(1:n,iperm))
!!$     call amps(iperm)%evaluate(n,p,0)
!!$  enddo

  call amps_col%init(2,n,part,order(1:n,1))
  call amps_col%evaluate(n,p,5)

  do i = 1, amps_col%nColOrd
     write (*,*) amps_col%amps(i),amps_col%perm(1:n-1,i)
  enddo
  
  stop
  
  call Tr_allocate(n)
  do jperm=1,nperm
     do iperm=1,nperm
        col_fac(iperm,jperm)=color_factor(iperm,jperm) ! colour factor for qqbar+gluons
     enddo
  enddo

  amp2=0d0
  do ih=1,2**n
     t=0d0
!!$  do ih=0,2**n-1
     do jperm=1,nperm    ! loop over permutations of conjugated amplitude
        ztemp=(0d0,0d0)
        do iperm=1,nperm ! loop over permutations of amplitude
           ztemp=ztemp+amps(iperm)%amps(ih)*col_fac(iperm,jperm)
        enddo
        t=t+dble(ztemp*dconjg(amps(jperm)%amps(ih)))
     enddo
     write (*,'(e20.12,x,a,x,b6.6)') t*(4*pi*alphas)**(n-2),'for helicity',(ih-1)
     amp2=amp2+t
  enddo
  

  t=0d0
  do jperm=1,nperm    ! loop over permutations of conjugated amplitude
     ztemp=(0d0,0d0)
     do iperm=1,nperm ! loop over permutations of amplitude
        ztemp=ztemp+amps_col%amps(iperm)*col_fac(iperm,jperm)
     enddo
     t=t+dble(ztemp*dconjg(amps_col%amps(jperm)))
  enddo
  write (*,'(e20.12,x,a,x,b6.6)') t*(4*pi*alphas)**(n-2)

  
  
  amp2=amp2*(4*pi*alphas)**(n-2)/dble(iden)

  write (*,*) 'Matrix element =',amp2,'GeV^',-(2*n-8)





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
!!$    if (jperm.ne.iperm) col_factor=col_factor*2d0 ! add factor 2 for off-diagonal terms
    color_factor=dble(col_factor)
  end function color_factor

  
end program test_QCD
