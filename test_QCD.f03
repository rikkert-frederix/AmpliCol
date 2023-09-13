! gfortran -fbounds-check -o test_QCD color_algebra.f95 feynmanrules.f03 amplitude_QCD.f03 test_QCD.f03 

program test_QCD
  use amplitude_mod
  use color_algebra
  implicit none
  type(amplitude) :: amplitudes
  type(amplitude),dimension(1:6) :: amps
  integer :: n
  integer,dimension(:),allocatable :: part
  integer,dimension(:,:),allocatable :: helmap,order
  integer :: i,ih,iden,iperm,jperm
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: amp2
  real(kind=8),parameter :: pi=3.14159265358979323846d0,alphas=0.118d0
  real(kind=8),dimension(6,6) :: col_fac
  complex(kind=8) :: ztemp
  n=5
  allocate(part(n))
  allocate(order(n,6))
  allocate(helmap(2**n,6))

  part(1:n)=[-1,21,21,21,-1]
  order(1:n,1)=[1,2,3,4,5]


!!$  call amplitudes%init_OneOrder(n,part,order(1,1))
!!$
!!$  do i=1,amplitudes%n_vert
!!$     write (*,*) i,':',amplitudes%interaction_list(i)%currents(1:2),':',amplitudes%interaction_list(i)%type
!!$  enddo
!!$  
!!$  do i=1,amplitudes%n_cur
!!$     write (*,*) i,':',amplitudes%current_list(i)%order(1:popcnt(amplitudes%current_list(i)%bin)),':',&
!!$          amplitudes%current_list(i)%type,':',amplitudes%current_list(i)%bin
!!$     write (*,*) amplitudes%current_list(i)%vertices(1:amplitudes%current_list(i)%n_vert)
!!$  enddo

  
  allocate(p(0:3,n))
  
  if (n.eq.5) then
     p(0:3,1)=(/   500.00000000000000      ,   0.0000000000000000      ,   0.0000000000000000      ,   500.00000000000000      /)
     p(0:3,2)=(/   500.00000000000000      ,   0.0000000000000000      ,   0.0000000000000000      ,  -500.00000000000000      /)
     p(0:3,3)=(/   362.35913008614813      ,   211.67753703547245      ,  -165.00336466515242      ,   243.45564097092699      /)
     p(0:3,4)=(/   248.22961274839184      ,   68.500483807498668      ,   237.98757697004567      ,   16.956932838275925      /)
     p(0:3,5)=(/   389.41125716545980      ,  -280.17802084297108      ,  -72.984212304893376      ,  -260.41257380920291      /)
  endif
  
!!$  call amplitudes%evaluate_OneOrder(n,p)



  

  part(1:n)=[-1,21,21,21,-1]

  order(1:n,1)=[1,2,3,4,5]
  order(1:n,2)=[1,2,4,3,5]
  order(1:n,3)=[1,3,2,4,5]
  order(1:n,4)=[1,3,4,2,5]
  order(1:n,5)=[1,4,2,3,5]
  order(1:n,6)=[1,4,3,2,5]

  do iperm=1,6
     call amps(iperm)%init_OneOrder(n,part,order(1:n,iperm))
     call amps(iperm)%evaluate_OneOrder(n,p)


     do ih=1,2**n
        helmap(ih,iperm)=0
        do i=1,n
           if (btest(ih-1,i-1)) helmap(ih,iperm)=ibset(helmap(ih,iperm),order(i,iperm)-1)
        enddo
        helmap(ih,iperm)=helmap(ih,iperm)+1
     enddo

     write (*,*) helmap(:,iperm)

     do ih=1,2**n
        helmap(ih,iperm)=0
        do i=1,n
           call mvbits(ih-1,order(i,iperm)-1,1,helmap(ih,iperm),i-1)
        enddo
        helmap(ih,iperm)=helmap(ih,iperm)+1
     enddo
     
     write (*,*) helmap(:,iperm)
  enddo

  
  call Tr_allocate(n)
  do jperm=1,6
     do iperm=1,6
        col_fac(iperm,jperm)=color_factor(iperm,jperm)
     enddo
  enddo


  
  do ih=1,32
     write (*,*) (amps(iperm)%amps(helmap(ih,iperm)),iperm=1,6)
  enddo

  
  amp2=0d0
  do ih=1,2**n
     do jperm=1,6    ! loop over permutations of conjugated amplitude
        ztemp=(0d0,0d0)
        do iperm=1,6 ! loop over permutations of amplitude
           ztemp=ztemp+amps(iperm)%amps(helmap(ih,iperm))*col_fac(iperm,jperm)
        enddo
        amp2=amp2+dble(ztemp*dconjg(amps(jperm)%amps(helmap(ih,jperm))))
     enddo
  enddo
  

  
  iden=3*8*2*2 ! initial status colours and helicities/polarisations
  iden=iden*2  ! final state identical particles
  amp2=amp2*(4*pi*alphas)**(n-2)/dble(iden)

  write (*,*) 'Matrix element =',amp2,'GeV^',-(2*n-8)





contains
  double precision function color_factor(iperm,jperm)
    ! Compute colour factor for permutation numbers 'iperm' and
    ! 'jperm'
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
