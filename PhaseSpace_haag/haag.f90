module haag
  use common
  private
  real(kind=8),parameter :: pi=3.1415926535897932d0
  real(kind=8) :: ksi_m

  real(kind=8),public :: s0
  logical,public :: debug=.false.,  flat=.false.,  open=.false.
  logical,public :: flat_split=.false., a1_split=.true.
  logical,public :: haag_style=.true.
  logical,public :: one_pt=.false.
  real(kind=8),dimension(:),allocatable :: masses
  real(kind=8),public :: tot_mass
  real(kind=8),public :: mass_sum
  integer(kind=4) :: mm
  integer(kind=4) :: ix, ndim

  integer(kind=4),dimension(:),allocatable :: order
  real(kind=8),dimension(:),allocatable :: invm,invm_min,invm_max,x
  real(kind=8),dimension(:,:),allocatable :: pp
  integer(kind=4),dimension(:,:),allocatable :: sets
  logical :: t_channel,includePDF
  real(kind=8) :: sqrtshat,sqrts,tau,ycm
  integer :: nquarks

  logical,parameter :: verbose=.true.
  logical,parameter :: exper=.false.
  real(kind=8),parameter :: ip=-1d0,ip_shat=-2d0
  real(kind=8),parameter :: vtiny=1d-12

  integer(kind=4), public :: n
  real(kind=8),public :: jaco,soft
  public :: haag_init,PS_haag

contains

  subroutine haag_init(sqrtsh,nn,m,o,part,s_cut,s_chan,include_pdf)
    implicit none
    real(kind=8),intent(in) :: sqrtsh
    integer(kind=4),intent(in) :: nn
    integer(kind=4),dimension(nn),intent(in) :: o,part
    integer(kind=4),dimension(nn) :: process,ord,temp_order
    integer(kind=4) :: glu,end,start
    real(kind=8),intent(in) :: s_cut(2)
    real(kind=8),dimension(nn),intent(in) :: m
    logical,intent(in) :: s_chan
    ! Should we include a PDF set? Currently, only the NNPDF2.3 NLO QED is available.
    logical,intent(in) :: include_pdf
    integer(kind=4) :: i,j
    sqrtshat=sqrtsh
    sqrts=sqrtsh
    if (verbose) then
       write (*,*) 'Setting up',nn,'particle phase-space'
       write (*,*) 'Total available energy, sqrt(s-hat) =',sqrtshat
       write (*,*) 'Cut on invariants used in the phase-space generation: abs((p_i+p_j)^2) >=',s_cut
    endif
    includePDF=include_pdf
    call gen23_deallocate
    next=nn
    n=nn-2
    ndim=3*(next-2)-4
    if (includePDF) ndim=ndim+2 ! the two Bjorken x's
    allocate(order(next))
    allocate(invm(maskr(next)))
    allocate(invm_min(maskr(next)))
    allocate(invm_max(maskr(next)))
    allocate(pp(0:3,0:maskr(next)))
    allocate(masses(next))
    pp(0:3,0:maskr(next))=0d0
    allocate(p(0:3,next))
    allocate(x(ndim))
    allocate(sets(0:next-2,2))
    ! masses of external particles
    do i=1,n
       if ((i.eq.1 .or. i.eq.2) .and. m(i).ne.0d0) then
          write (*,*) 'ERROR in gen23_init() -- ', &
               & 'incoming particles should be massless'
          write (*,*) m
          stop 1
       endif
       invm(ibset(0,i-1))=m(i)**2
    enddo
    masses=m
    tot_mass=sum(masses)
    call get_approx_s0()

    ! Bring the colour order to a canonical order (first in the list
    ! should be particle 1, i.e., the first incoming particle).
    do i=1,next
       if (o(i).eq.1) then
          do j=0,next-1
             order(j+1)=o(1+mod(i+j-1,next))
          enddo
          exit
       endif
    enddo
    ord=order
    process=part
    nquarks=0
    do i=1,next
       if ((abs(process(i)).ge.1) .and. abs(process(i)).le.6) then
           nquarks=nquarks+1
       endif
       if ((i.le.2) .and. ((abs(process(i)).ge.1) .and. abs(process(i)).le.6))  then
          process(i)=-process(i)
       endif
    enddo
    if (verbose) write (*,*) 'Canonical order',order
    ! Define the sets from the colour order. Set 1 contains all the
    ! particles between the first and second incoming particles. Set 2
    ! contains the particles between the second and first incoming
    ! particles.
    call define_gen_order(next,process,o,order)
    do i=1,next
       if (order(i).eq.1) then
          do j=0,next-1
             temp_order(j+1)=order(1+mod(i+j-1,next))
          enddo
          exit
       endif
    enddo

    temp_order=ord ! REMOVE comment for OLD
    sets=0
    i=0
    do i=2,next
       if (temp_order(i).eq.2) then
          do j=i+1,next
             sets(0,2)=ibset(sets(0,2),temp_order(j)-1)
          enddo
          sets(1:i-2,1)=temp_order(2:i-1)
          sets(1:next-i,2)=temp_order(i+1:next)
          exit
       endif
       sets(0,1)=ibset(sets(0,1),temp_order(i)-1)
    enddo

    if (verbose) then
       write (*,*) "set 1:",sets(:,1)
       write (*,*) "set 2:",sets(:,2)
    endif
    !stop 1
    if (verbose) then
       write (*,*) "Power in importance sampling:",ip
    endif
  end subroutine haag_init

  subroutine get_approx_s0
   implicit none
   double precision :: eta_max,deltan,theta_1,theta_2,Phi
   eta_max = log(tan(cos(asin(pT_min/sqrtshat))))
   deltan=0.00001d0
   Phi=DRjj_min
   theta_1 = 2d0*atan(exp(-eta_max))
   theta_2 = 2d0*atan(exp(-eta_max-deltan))
   s0=2d0*(pT_min*pT_min/(sin(theta_1)*sin(theta_2))*(1d0-cos(theta_1)*cos(theta_2))-pT_min*pT_min*cos(Phi))
   if (sqrt_s_min.gt.0d0) then
      s0=sqrt_s_min**2
   endif
  end subroutine get_approx_s0

  subroutine define_gen_order(next,process,o,order)
    implicit none
    integer :: next,i,j,glu,quark,aquark
    integer(kind=4),dimension(next) :: process,o,order
    integer :: in1,in2
    order=0
    if (all(process.eq.21)) then
       do i=1,next
        if (o(i).eq.1) then
          do j=0,next-1
             order(j+1)=o(1+mod(i+j-1,next))
          enddo
          exit
        endif
       enddo
    else
      do i=1,next
        if (process(i).ge.1.and.process(i).le.6)then
            quark=i
        elseif (-process(i).ge.1.and.-process(i).le.6)then
            aquark=i
        endif
      enddo
        if ((quark.le.2) .and. (aquark.le.2))then
            order(1)=aquark
            order(2)=quark
            glu=0
            do i=1,next
               if (o(i).gt.2) then
                  order(3+glu)=o(i)
                  glu=glu+1
                endif
            enddo
         elseif ((quark.le.2).and.(aquark.gt.2))then
            order(1)=quark
            order(2)=mod(quark,2)+1
            order(next)=aquark
            glu=0
            do i=1,next
               if (o(i).ne.quark .and. o(i).ne.aquark.and.o(i).gt.2) then
                  order(3+glu)=o(i)
                  glu=glu+1
                endif
            enddo
         elseif ((quark.gt.2).and.(aquark.le.2))then
            order(1)=aquark
            order(2)=mod(aquark,2)+1
            order(next)=quark
            glu=0
            do i=1,next
               if (o(i).ne.quark .and. o(i).ne.aquark.and.o(i).gt.2) then
                  order(3+glu)=o(i)
                  glu=glu+1
                endif
            enddo
         elseif ((quark.gt.2).and.(aquark.gt.2))then
            if (quark.lt.aquark) then
              order(1)=2
              order(2)=1
              order(next-1)=aquark
              order(next)=quark
            else
              order(1)=1
              order(2)=2
              order(next-1)=aquark
              order(next)=quark
            endif
            glu=0
            do i=1,next
               if (o(i).ne.quark .and. o(i).ne.aquark.and.o(i).gt.2) then
                  order(3+glu)=o(i)
                  glu=glu+1
                endif
            enddo
         endif
         do i=1,next
            if (o(i).eq.1) in1=i
            if (o(i).eq.2) in2=i
         enddo
         if (in1.lt.in2) then
            order(1)=2
            order(2)=1
         elseif (in1 .gt. in2) then
            order(1)=1
            order(2)=2
         endif
         glu=0
         do i=1,next
            if (o(i).eq.quark.and.o(i).gt.2) then
               order(next-1)=quark
            elseif (o(i).eq.aquark.and.o(i).gt.2) then
               order(next)=aquark
            elseif (o(i).gt.2) then
               order(3+glu)=o(i)
               glu=glu+1
            endif
         enddo
    endif
    order=o
  end subroutine define_gen_order

  subroutine setup_PS_cuts(s_cut)
    ! Given s_cut = abs((p_i+p_j)^2), fills the minimum (s-channel)
    ! and/or maximum (t-channel) values the invariants can be in the
    ! phase-space generation. Does not apply these cuts on invariants
    ! not used in the phase-space generation.
    implicit none
    real(kind=8),intent(in) :: s_cut(2)
    real(kind=8) :: mass
    integer(kind=4) :: i,j,npart
    invm_min=0d0
    invm_max=0d0
    do i=1,maskr(next)
       npart=popcnt(i)
       if (btest(i,0).and.btest(i,1)) then
          invm_min(i)=0d0
       elseif (btest(i,0).or.btest(i,1)) then
          if (npart.eq.2 .or. npart.eq.next-2) then
             invm_max(i)=-s_cut(1)
          endif
       else
          mass=0d0
          do j=0,next-1
             if (btest(i,j)) then
                mass=mass+sqrt(invm(ibset(0,j)))
             endif
          enddo
          invm_min(i)=max(s_cut(2)*(npart)*(npart-1)/2d0,mass**2)
       endif
    enddo
  end subroutine setup_PS_cuts

  subroutine gen23_deallocate
    implicit none
    if (allocated(order)) deallocate(order)
    if (allocated(invm)) deallocate(invm)
    if (allocated(invm_min)) deallocate(invm_min)
    if (allocated(invm_max)) deallocate(invm_max)
    if (allocated(pp)) deallocate(pp)
    if (allocated(p)) deallocate(p)
    if (allocated(x)) deallocate(x)
    if (allocated(sets)) deallocate(sets)
  end subroutine gen23_deallocate

  subroutine PS_haag(xx)
    ! Wrapper for the routine that generates the momenta.
    implicit none
    real(kind=8),dimension(99),intent(in) :: xx
    integer(kind=4) :: i
    x(1:ndim)=xx(1:ndim)
    jac=1d0
    ix=0
    if (includePDF) call generate_initial_state
    call generate_momenta
    do i=1,next
       p(0:3,i)=pp(0:3,i)
    enddo
  end subroutine PS_haag

  subroutine generate_momenta
    implicit none
    real(kind=8),dimension(0:3,next) :: q, qk, qlab
    real(kind=8),dimension(0:3) :: qtot,qtotm, tot, bst, bst_back, qin1,qin2
    real(kind=8) :: costheta,phi,sintheta,dum,scale,a_sum
    integer(kind=4):: i, t,j,k,m,first,second
    real(kind=8) :: sk, E, pz, pT, xy,kappa,Atilde,wsq,v,esum,R,min
    real(kind=8),dimension(0:3) :: Qm,Qnm
    real(kind=8) :: Qz,E1,E2,s1,s2,s,Qt,antenna,test,mass1,mass2,mass_in
    logical :: m1
    integer(kind=8),dimension(:),allocatable :: subperm1(:), subperm2(:),subperm,subperm_rest
    integer(kind=8),dimension(next-2) :: perm_final
    real(kind=8),dimension(next-2) :: schan_ran,schan_ran_sorted
    real(kind=8),dimension(0:3) :: q1_ref,q2_ref
    
    if(.not.allocated(subperm1)) allocate(subperm1(1:next-2))
    if(.not.allocated(subperm2)) allocate(subperm2(1:next-2))
    if(.not.allocated(subperm)) allocate(subperm(1:next-2))
    if(.not.allocated(subperm_rest)) allocate(subperm_rest(1:next-2))
    mm = 0
    do i=1,next-2
       if (sets(i,1) .ne. 0) then
          mm = mm+1
       endif
    enddo
    subperm1 = sets(1:,1)
    subperm2 = sets(1:,2)

    !! Initial momenta in lab frame
    qin1(0) = sqrtshat/2d0
    qin1(1:2) = 0d0
    qin1(3) = qin1(0)
    qin2(0) = sqrtshat/2d0
    qin2(1:2) = 0d0
    qin2(3) = -qin2(0)
    q(:,1) = qin1
    q(:,2) = qin2

    ! Total momentum P=p_1+p_2
    qk(:,n) = qin1 + qin2
    jaco = 1d0 ! jacobian only for the antenna part
    soft = 1d0

! First split the antenna to two subantennae: Q_m and Q_{n-m}
if ((mm .gt. 1).and.(n-mm .gt. 1)) then ! Do m>1 type splitting
     !write(*,*) 'doing m>1 type splitting'
     mass1=0d0
     do i=1,mm
       mass1=mass1+masses(subperm1(i))
     enddo
       mass2=0d0
     do i=1,n-mm
       mass2=mass2+masses(subperm2(i))
     enddo
     call generate_split_Qm_Qnm(sqrtshat**2,mass1,mass2,mm,Qm,Qnm)

     mass_sum=0d0

! Generate Q_{m} antenna
  if (mm .gt. 2) then
      call basic_antenna(q(0:3,subperm1(mm)),masses(subperm1(mm)),qk(0:3,mm-1),-1d0,&
         Qm,qin2,qin1,0,.false.,mm)
      mass_sum=mass_sum+masses(subperm1(mm))
  else
      call basic_antenna(q(0:3,subperm1(2)),masses(subperm1(2)),qk(0:3,1),&
           masses(subperm1(1)),Qm,qin2,qin1,0,.false.,mm)
      mass_sum=mass_sum+masses(subperm1(2))
      mass_sum=mass_sum+masses(subperm1(1))
  endif
  do i=1,mm-3
      call basic_antenna(q(0:3,subperm1(mm-i)),masses(subperm1(mm-i)),qk(0:3,mm-i-1),-1d0,&
           qk(0:3,mm-i),q(0:3,subperm1(mm-i+1)),qin1,i,.false.,mm)
      mass_sum=mass_sum+masses(subperm1(mm-i))
  enddo
  if (mm .gt. 2) then
      call basic_antenna(q(0:3,subperm1(2)),masses(subperm1(2)),qk(0:3,1),masses(subperm1(1)),&
               qk(0:3,2),q(0:3,subperm1(3)),qin1,mm-2,.false.,mm)
      mass_sum=mass_sum+masses(subperm1(2))
      mass_sum=mass_sum+masses(subperm1(1))
  endif
      q(0:3,subperm1(1)) = qk(0:3,1)

! Generate Q_{n-m} antenna
  if (n-mm .gt. 2) then
    call basic_antenna(q(0:3,subperm2(n-mm)),masses(subperm2(n-mm)),qk(0:3,n-1),-1d0,&
           Qnm,qin1,qin2,0,.false.,n-mm)
    mass_sum=mass_sum+masses(subperm2(n-mm))
  else
    call basic_antenna(q(0:3,subperm2(2)),masses(subperm2(2)),qk(0:3,n-1),&
        masses(subperm2(1)),Qnm,qin1,qin2,0,.false.,n-mm)
    mass_sum=mass_sum+masses(subperm2(1))
    mass_sum=mass_sum+masses(subperm2(2))
  endif
  do i=1,(n-mm)-3
    call basic_antenna(q(0:3,subperm2(n-mm-i)),masses(subperm2(n-mm-i)),qk(0:3,n-1-i),-1d0,&
               qk(0:3,n-i),q(0:3,subperm2(n-mm-i+1)),qin2,i,.false.,n-mm)
    mass_sum=mass_sum+masses(subperm2(n-mm-i))
  enddo
  if (n-mm .gt. 2) then
    call basic_antenna(q(0:3,subperm2(2)),masses(subperm2(2)),qk(0:3,mm+1),masses(subperm2(1)),&
               qk(0:3,mm+2),q(0:3,subperm2(3)),qin2,n-mm-2,.false.,n-mm)
    mass_sum=mass_sum+masses(subperm2(2))
    mass_sum=mass_sum+masses(subperm2(1))
  endif
    q(0:3,subperm2(1)) = qk(0:3,mm+1)


! Do m=1 type splitting
elseif (((mm .eq. 1).or.(n-mm .eq. 1))) then 
  !write(*,*) 'doing m=1 type splitting'
  m1 = .true.
  if (mm .eq. 1) then
     subperm=subperm1
     subperm_rest=subperm2
     q1_ref = qin2
     q2_ref = qin1
  elseif (next-2-mm .eq. 1) then
     subperm = subperm2
     subperm_rest=subperm1
     q1_ref = qin1
     q2_ref = qin2
  endif

  mass_sum=0d0
  if (n .gt. 2) then
     call basic_antenna(q(0:3,subperm(1)),masses(subperm(1)),qk(0:3,n-1),-1d0,&
                  qk(0:3,n),q1_ref,q2_ref,0,m1,n)
     mass_sum=mass_sum+masses(subperm(1))
  else
     call basic_antenna(q(0:3,subperm(1)),masses(subperm(1)),qk(0:3,n-1),&
                 masses(subperm_rest(1)),qk(0:3,n),q1_ref,q2_ref,0,m1,n)
     mass_sum=mass_sum+masses(subperm(1))
     mass_sum=mass_sum+masses(subperm_rest(1))
  endif
  if (n .gt. 3) then
       call basic_antenna(q(0:3,subperm_rest(n-1)),masses(subperm_rest(n-1)),qk(0:3,n-1-1),&
                 -1d0,qk(0:3,n-1),q2_ref,q1_ref,1,m1,n)
       mass_sum=mass_sum+masses(subperm_rest(n-1))
  elseif (n .eq. 3) then 
       call basic_antenna(q(0:3,subperm_rest(2)),masses(subperm_rest(2)),qk(0:3,n-2),&
         masses(subperm_rest(1)),qk(0:3,n-1),q2_ref,q1_ref,1,.false.,n)
       mass_sum=mass_sum+masses(subperm_rest(2))
       mass_sum=mass_sum+masses(subperm_rest(1))
  endif
  do i=2,n-3
   call basic_antenna(q(0:3,subperm_rest(n-i)),masses(subperm_rest(n-i)),qk(0:3,n-i-1),-1d0,&
                  qk(0:3,n-i),q(0:3,subperm_rest(n-i+1)),q1_ref,i,.false.,n)
   mass_sum=mass_sum+masses(subperm_rest(n-i))
  enddo
  if (n .gt. 3) then
          call basic_antenna(q(0:3,subperm_rest(2)),masses(subperm_rest(2)),qk(0:3,1),&
               masses(subperm_rest(1)),qk(0:3,2),&
               q(0:3,subperm_rest(3)),q1_ref,n-2,.false.,n)
     mass_sum=mass_sum+masses(subperm_rest(2))
     mass_sum=mass_sum+masses(subperm_rest(1))
  endif
     q(0:3,subperm_rest(1)) = qk(0:3,1)

! Do m=0 type splitting
else  
  !write(*,*) 'doing m=0 type splitting'
  m1 = .false.
  mass_sum=0d0

  if (mm .eq. 0) then
     perm_final=subperm2
     q1_ref = qin1
     q2_ref = qin2
  elseif (n-mm .eq. 0) then
     perm_final=subperm1
     q1_ref = qin2
     q2_ref = qin1
  endif
  if (n .gt. 2) then
     call basic_antenna(q(0:3,perm_final(n)),masses(perm_final(n)),qk(0:3,n-1),mass_in,&
                qk(0:3,n),q1_ref,q2_ref,0,m1,n)
     mass_sum=mass_sum+masses(perm_final(n))
  else
     call basic_antenna(q(0:3,perm_final(n)),masses(perm_final(n)),qk(0:3,n-1),&
               masses(perm_final(1)),qk(0:3,n),q1_ref,q2_ref,0,m1,n)
     mass_sum=mass_sum+masses(perm_final(n))
  endif
  do i=1,n-3
     call basic_antenna(q(0:3,perm_final(n-i)),masses(perm_final(n-i)),qk(0:3,n-i-1),&
               mass_in,qk(0:3,n-i),q(0:3,perm_final(n-i+1)),q2_ref,i,.false.,n)
     mass_sum=mass_sum+masses(perm_final(n-i))
  enddo
  if (n .gt. 2) then
        call basic_antenna(q(0:3,perm_final(2)),masses(perm_final(2)),qk(0:3,1),&
        masses(perm_final(1)),qk(0:3,2),q(0:3,perm_final(3)),q2_ref,n-2,.false.,n)
  endif
  mass_sum=mass_sum+masses(perm_final(2))
  mass_sum=mass_sum+masses(perm_final(1))
  q(0:3,perm_final(1)) = qk(0:3,1)
endif

    do i=1,n+2
      pp(0:3,i) = q(0:3,i)
    enddo

    !call check_momenta(p,masses)

    ! Compute the weight (i.e. jacobian)
    ! The usual 2*pi factors for the phase-space
    jac=jac*jaco*soft/((2d0*pi)**(3*n-4)) !wgt
    jac=jac/(2d0*sqrtshat**2)

    !write(*,*) ' '
    !do i=1,next
    !   write(*,*) pp(:,i)
    !enddo


  end subroutine generate_momenta

  subroutine generate_initial_state
    implicit none
    call generate_tau
    call generate_y
    sqrtshat=sqrt(tau)*sqrts
    xbjrk(1)=sqrt(tau)*exp(ycm)
    xbjrk(2)=sqrt(tau)*exp(-ycm)
  end subroutine generate_initial_state

  subroutine generate_tau
    implicit none
    integer :: i
    real(kind=8) :: smin,smax,shat
    smin=(next-2)*(next-3)*s0
    smax=sqrts**2
    ix=ix+1
    call random_to_var(x(ix),ip_shat,smin,smax,shat,jac)
    tau=shat/smax
    jac=jac/smax
  end subroutine generate_tau

  subroutine generate_y
    implicit none
    real(kind=8) ::  ymin,ymax
    ymin= log(tau)/2d0
    ymax=-log(tau)/2d0
    ix=ix+1
    call random_to_var(x(ix),0d0,ymin,ymax,ycm,jac)
  end subroutine generate_y

  subroutine basic_antenna(p1,mass1,p2,mass2,P,q1,q2,i,m1,maxn)
    ! Incoming momentum: P
    ! Reference momenta: q1,q2
    ! Outgoing momenta: p1, p2
    implicit none
    real(kind=8),dimension(0:3),intent(in) :: P, q1,q2
    real(kind=8),dimension(0:3),intent(out) :: p1, p2
    real(kind=8),intent(in) :: mass1,mass2
    real(kind=8),dimension(0:3) :: qtot,qtotm,q1_cmf,q2_cmf,p1_cmf, p2_cmf,Pm,P_cmf
    real(kind=8),dimension(0:3) :: qmir_cmf,qmir
    real(kind=8),dimension(0:3) :: test1,test2,test_out
    real(kind=8) :: esum,dum,costheta,sintheta,phi
    real(kind=8) :: s1, s2, a1, a2, s, gs, a1max,a1min
    integer :: i,k
    real(kind=8):: h1,h
    real(kind=8),dimension(3) :: solution
    real(kind=8) ::  plot,beta,a1cut,a2cut
    integer :: maxn
    logical :: m1
    real(kind=8) :: schan_ran
    real(kind=8) z_sign
    real(kind=8),dimension(2) :: buff,min_point1,min_point2
    double precision :: inv,w1,w2,w,R
    integer :: term 
    double precision,dimension(2) :: a2pm,a1pm
    real(kind=8),external :: ran2

    real(kind=8) :: tmin,tmax,t1,t2,yr,pt2
    real(kind=8),dimension(0:3) :: pass1,pass2
    real(kind=8) :: pt2min,pt2max,pt,y,ymin,ymax

    k = maxn - i   !  k is the number of particles remaining to generate
    s1 = mass1     ! mass of the final state particle to be generated
    z_sign=sign(1d0,q2(3))
    s = (dot(P,P)) ! Ingoing inv mass

    ! boost qi to CMF (P rest frame)
    if (dot(P,P).le.0d0) return
    esum=dsqrt(dot(P,P))
    Pm(0)=P(0)
    Pm(1:3)=-P(1:3)
    P_cmf(0)=esum
    P_cmf(1:3)=(/0d0,0d0,0d0/)

    call boostm(q1,Pm,esum,q1_cmf)
    call boostm(q2,Pm,esum,q2_cmf)

    ! Split into the long (L) decomposition of the q1 (relevant for massive)
    beta = dsqrt(threedot(q1_cmf(1:3),q1_cmf(1:3)))/q1_cmf(0) ! needed only when massive
    q1_cmf(1:3) = q1_cmf(1:3)/beta ! if massless, beta=1
    if (mass1 .eq. 0d0) then
       beta=1.d0
    endif

    if (threedot(q1_cmf(1:3),q1_cmf(1:3))*threedot(q2_cmf(1:3),q2_cmf(1:3)).eq.0d0) return
    ! angles between q1,q2 in CMF_k frame
    costheta = threedot(q1_cmf(1:3),q2_cmf(1:3))/ &
            (sqrt(threedot(q1_cmf(1:3),q1_cmf(1:3))*threedot(q2_cmf(1:3),q2_cmf(1:3))))
    sintheta = dsqrt(1D0-costheta**2)

    !Generate s2
    if (m1 .and. (i .eq. 0)) then 
        ! If m=1 type first splitting -> do also a1,phi sampling here!
         if (haag_style.and..not.one_pt) then
            call generate_first_single(k,s,s1,q1_cmf,P_cmf,a1,s2)
            if (k.eq.2) s2=0d0
            a2 = 300d0
            goto 20
         elseif(.not.haag_style.and. .not.one_pt) then
            yr = dsqrt(kallen(sqrtshat**2,0d0,(next-3)*(next-4)/2d0*s0))
            tmin=(-sqrtshat**2+(next-3)*(next-4)/2d0*s0-yr)/2d0
            tmax=(-sqrtshat**2+(next-3)*(next-4)/2d0*s0+yr)/2d0
            if (-s0.ne.0d0) tmax=min(-s0,tmax)
            if (0d0.ne.0d0) tmin=max(tmin,0d0)
            ix=ix+1
            call random_to_var(x(ix),-1d0,tmin,tmax,t1,jac)
            if (tmin.ge.tmax) then
                jac=-2d0
                return
            endif
            tmin=-sqrtshat**2-t1+(next-3)*(next-4)/2d0*s0
            tmax=0d0*(0d0-sqrtshat-t1)/(0d0-t1)

            if (-s0.ne.0d0) tmax=min(-s0,tmax)
            if (0d0.ne.0d0) tmin=max(tmin,0d0)

            if (tmin.ge.tmax) then
                jac=-2d0
                return
            endif
            ix=ix+1
            call random_to_var(x(ix),-1d0,tmin,tmax,t2,jac)
            ix=ix+1
            call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,jac)
            pt2=t1*t2/sqrtshat**2+ &
               & 0d0**2/sqrtshat**2-(t1+t2)*0d0/sqrtshat**2-0d0
            pass1(0) =(-t1-t2+2d0*0d0)/(2d0*sqrtshat)
            pass1(1)=sqrt(pt2)*cos(phi)
            pass1(2)=sqrt(pt2)*sin(phi)
            pass1(3)=(t1-t2)/(2d0*sqrtshat)

            pass2(0)=sqrtshat-pass1(0)
            pass2(1:3)=-pass1(1:3)

            
            jac = jac/(4d0*sqrt(kallen(sqrtshat**2,0d0,0d0)))
            p1 = pass1
            p2 = pass2
            
            return
         elseif (one_pt) then
            pt2min=pt_min**2
            pt2max=sqrtshat**2/4d0
            ix=ix+1
            call random_to_var(x(ix),-1d0,pt2min,pt2max,pt2,jac)
            pt=sqrt(pt2)
            ! generate phi
            ix=ix+1
            call random_to_var(x(ix),0d0,-pi,pi,phi,jac)
            ! generate rapidity
            ix=ix+1
            ymin=-min(eta_max,acosh(sqrtshat/pt))
            ymax=min(eta_max,acosh(sqrtshat/pt))
            call random_to_var(x(ix),0d0,ymin,ymax,y,jac)
            ! fill momentum
            pass1(1)=pt*cos(phi)
            pass1(2)=pt*sin(phi)
            pass1(3)=pt*sinh(y)
            pass1(0)=pt*cosh(y)

            pass2(0)=sqrtshat-pass1(0)
            pass2(1:3)=-pass1(1:3)

            if (dot(pass2,pass2).lt.0d0) return

            p1 = pass1
            p2 = pass2
            jac=jac/(32d0*pi*pi*pi)
            return

         endif
    else
       if (k .ge. 3) then
         call generate_s2(k,s,s1,s2,q1_cmf,P_cmf)
       else
         s2 = mass2
         gs = 1d0
         jaco = jaco/gs
       endif
    endif
     



    ! cuts on a1 and a2
    a2cut = (k-1)*(s0/2d0)/(dot(q2_cmf,P_cmf))  ! should be (k+1) dont change!
    a1cut = (s0/2d0)/(dot(q1_cmf,P_cmf))

    h = 1d-6
    if (h .gt. get_min_a2_bound(s,s1,s2,costheta)) then
        h = get_min_a2_bound(s,s1,s2,costheta)
    elseif (h .gt. a2cut) then
        h = a2cut
    endif
    h= 1d-8 ! for using the "h"-technique
    !h = 0d0  ! for using the partial decomposition
    !h = -1d0  ! for using actually h=0

    R = ran2()
    call get_partial_weights(w1,w2,s,s1,s2,a1cut,a2cut,h,costheta)
    w=w1+w2
    if (R.lt.w1/w) then
       term=1
    elseif ((R.gt.w1/w).and.(R.lt.(w1+w2)/w)) then
       term=2
    endif

    if (h.eq.0d0) then
    if (((i .eq. 0) .or. (m1 .and. (i .le. 1)))) then
      call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
      call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
    else
      if (term.eq.1) then
        call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
        call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
      elseif (term.eq.2) then
        call generate_a2_term2(i,m1,maxn,s,s1,s2,costheta,a2cut,beta,a2)
        call generate_a1_term2(i,m1,maxn,a2,s,s1,s2,costheta,a1cut,a1)
      endif
    endif

    elseif (abs(h).gt.0d0) then
        call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
        call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
    endif


    ! Mapping back to the momenta p1,p2 in CMF
    20    p1_cmf(0) = (s+s1-s2)/(2D0*sqrt(s))
    solution = solver(s,s1,s2,q1_cmf,q2_cmf,z_sign,a1,a2,p1_cmf(0),P_cmf,costheta)
    p1_cmf(1) = solution(1)
    p1_cmf(2) = solution(2)
    p1_cmf(3) = solution(3)

    p2_cmf(0) = sqrt(s) - p1_cmf(0)
    p2_cmf(1:3) = -p1_cmf(1:3)

    qmir_cmf(0)=q2_cmf(0)
    qmir_cmf(1:3)=-q2_cmf(0)

    ! Boost back to lab frame
    call boostm(p1_cmf,P,esum,p1)
    call boostm(p2_cmf,P,esum,p2)

    qmir(0)=q2(0)
    qmir(1:3)=-q2(0)

  end subroutine basic_antenna

  subroutine get_partial_weights(w1,w2,s,s1,s2,a1cut,a2cut,h,cos)
    implicit none
    integer :: i
    real(kind=8) :: w1,w2
    real(kind=8) :: s,s1,s2,cos
    real(kind=8) :: h1,a1min,a1max,a2min,a2max,f_h1,beta,dum,h
    real(kind=8) :: a1cut,a2cut
    real(kind=8),dimension(3) :: buff
    real(kind=8) :: amin,amax,bmin,bmax
    
    h1=0d0
    a1max = 0.5d0*(1d0-(s2-s1)/(s)+dsqrt(kallen(1d0,s1/s,s2/s)))
    a1min = 0.5d0*(1d0-(s2-s1)/(s)-dsqrt(kallen(1d0,s1/s,s2/s)))
    if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) then
       a1min=a1cut
    endif

    buff = f_func_term1(-h1,cos,s,s1,s2,h)
    f_h1 = buff(2)
    buff = f_func_term1(a1min,cos,s,s1,s2,h)
    Amin= (a1min+ h1 + buff(2)-f_h1) &
         /(a1min+ h1 + buff(2) +f_h1)
    buff = f_func_term1(a1max,cos,s,s1,s2,h)
    Amax= (a1max + h1 + buff(2)-f_h1) &
          /(a1max + h1 + buff(2) +f_h1)
    w1 = (1d0/f_h1)*(log(Amax)-log(Amin))

    a2max = 0.5d0*((1d0+(s2-s1)/s)+dsqrt(kallen(1d0,s1/s,s2/s)))
    a2min = 0.5d0*((1d0+(s2-s1)/s)-dsqrt(kallen(1d0,s1/s,s2/s)))
    if ((a2cut.gt.a2min).and.(a2cut.lt.a2max)) then
       a2min=a2cut
    endif

    buff = f_func_term2(-h1,cos,s,s1,s2)
    f_h1 = buff(2)
    buff = f_func_term2(a2min,cos,s,s1,s2)
    bmin= (a2min+ h1 + buff(2)-f_h1) &
         /(a2min+ h1 + buff(2) +f_h1)
    buff = f_func_term2(a2max,cos,s,s1,s2)
    bmax= (a2max + h1 + buff(2)-f_h1) &
          /(a2max + h1 + buff(2) +f_h1)
    w2=(1d0/f_h1)*(log(bmax)-log(bmin))

  end subroutine get_partial_weights

  real(kind=8) function kallen(a,b,c)
    implicit none
    real(kind=8) :: a,b,c
    kallen=a**2d0+b**2d0+c**2d0-2d0*a*b-2d0*a*c-2d0*b*c
    if (c.eq.0d0) then
       kallen=(a-b)**2
    endif
  end function kallen

  real(kind=8) function get_min_a2_bound(s,s1,s2,cos)
    implicit none
    real(kind=8),dimension(2) :: min_point1,min_point2
    real(kind=8) :: a1,s,s1,s2,cos,sin,r,A
    real(kind=8) :: p,q,a1plus,a1minus

    sin = dsqrt(1d0-cos**2)
    r = (s2-s1)/s
    A = 1d0 - r
    p = -A
    q = (A**2*sin**2 + 4d0*cos**2*s1/s)/(4d0)
    a1minus = -0.5d0*p + dsqrt((0.5d0*p)**2-q)
    a1plus =  -0.5d0*p - dsqrt((0.5d0*p)**2-q)
    min_point1 = a2_pm(a1minus,s,s1,s2,cos)
    min_point2 = a2_pm(a1plus,s,s1,s2,cos)
    get_min_a2_bound = max(min_point1(2),min_point2(2))
  end function get_min_a2_bound

  subroutine generate_split_Qm_Qnm(s,mass1,mass2,mn,Qm,Qnm)
    implicit none
    real(kind=8),dimension(0:3),intent(out) :: Qm,Qnm
    real(kind=8),dimension(0:3) :: dummy
    integer :: mn
    real(kind=8) :: mass1,mass2
    real(kind=8) :: c1,c2,m,dum,s,s1,s2,E1,E2,phi,Qt,Qz
    real(kind=8) :: RHS,max,min,sum_w,R
    real(kind=8),dimension(4) :: g1,g2,g3,d,e,w    
    integer :: i,pick
    real(kind=8),external :: ran2
    real(kind=8) :: dm,a1,a1min,a1max,a2,comm,root,a1cut
    real(kind=8),dimension(3) :: solution

    ! Generate s1,s2 invariants
    m=dsqrt(s)
    c1 = dsqrt(mass1 + mn*(mn-1)*s0/2d0)
    c2 = dsqrt(mass2 + (n-mn)*(n-mn-1)*s0/2d0)

    if (.not. flat_split) then
      ix = ix +1
      call random_to_var(x(ix),-1d0,c1**2,(m-c2)**2,s1,dum)
      ix = ix +1
      call random_to_var(x(ix),-1d0,c2**2,(m-dsqrt(s1))**2,s2,dum)
      soft = soft*(log((m-sqrt(s1))**2)-log(c2**2))
      soft = soft*log((m-c2)**2/(c1**2))
      jaco = jaco*s1*s2
    endif

    if (flat_split) then
      ix = ix + 1
      call random_to_var(x(ix),0d0,c1**2,(m-c2)**2,s1,dum)
      soft = soft*((m-c2)**2-c1**2)
      ix = ix + 1
      call random_to_var(x(ix),0d0,c2**2,(m-dsqrt(s1))**2,s2,dum)
      soft = soft*((m-dsqrt(s1))**2-c2**2)
    endif

    if (.not.a1_split) then
      E1 = (s+s1-s2)/(2d0*dsqrt(s))
      E2 = dsqrt(s) - E1

      ix = ix + 1
      call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,dum)
      soft = soft*2d0*pi

      ! multichanneling for Qz sampling
      ! Not the same input as from paper!!! 
      g1 = (/-1d0,+1d0,+1d0,-1d0/)
      g2 = (/1d0,1d0,-1d0,-1d0/)
      g3 = (/1d0,-1d0,1d0,-1d0/)
      max = dsqrt(E1**2-s1)
      min = -dsqrt(E1**2-s1)

      do i=1,4
         w(i) = (1d0/(4d0*E1*E2))*1d0/(E2+g1(i)*E1)&
                  *((g2(i)*log(E1+g2(i)*max)-g2(i)*log(E2+g3(i)*max))&
                   -(g2(i)*log(E1+g2(i)*min)-g2(i)*log(E2+g3(i)*min)))
      enddo
      !w(1) = (1d0/(4d0*E1*E2))*1d0/(E1-E2)*((log(E1+max)-log(E2+max))-(log(E1+min)-log(E2+min)))
      !w(2) = (1d0/(4d0*E1*E2))*1d0/(E1+E2)*((log(E1+max)-log(E2-max))-(log(E1+min)-log(E2-min)))
      !w(3) = (1d0/(4d0*E1*E2))*1d0/(E1+E2)*((-log(E1-max)+log(E2+max))-(-log(E1-min)+log(E2+min)))
      !w(4) = (1d0/(4d0*E1*E2))*1d0/(E1-E2)*((-log(E1-max)+log(E2-max))-(-log(E1-min)+log(E2-min)))
      sum_w = sum(w)
      R = ran2()
     
      if (R .lt. w(1)/sum_w) then
        pick = 1
      elseif ((R .lt. (w(2)+w(1))/sum_w) .and. (R .gt. w(1)/sum_w)) then
        pick = 2
      elseif ((R .lt. (w(3)+w(2)+w(1))/sum_w) .and. (R .gt. (w(2)+w(1))/sum_w)) then 
        pick = 3
      elseif (R .gt. (w(3)+w(2)+w(1))/sum_w) then
        pick = 4
      endif

      if (.not. flat_split) then
        ix = ix + 1
        call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
        RHS = ((E1+g2(pick)*min)/(E2+g3(pick)*min))&
           *exp(R*w(pick)*4d0*E1*E2*(E2+g1(pick)*E1)/g2(pick))
        Qz = (E1-E2*RHS)/(g3(pick)*RHS-g2(pick))
        soft = soft*sum_w
        jaco = jaco*(E1**2-Qz**2)*(E2**2-Qz**2)
      endif

      if (flat_split) then
        ix = ix + 1
        call random_to_var(x(ix),0d0,min,max,Qz,dum)
        soft = soft*(max-min)
      endif

      soft = soft/(4d0*dsqrt(s)) !! Double check where this comes from???
      Qt = dsqrt(E1**2-s1-Qz**2)
      Qm = (/E1,Qt*cos(phi),Qt*sin(phi),Qz/)
      Qnm = (/E2,-Qt*cos(phi),-Qt*sin(phi),-Qz/)
    else 
      ! do a1 sampling
      comm = 0.5d0*(s+s1-s2)/s
      root = dsqrt(comm**2-s1/s)
      a1min = comm - root
      a1cut = 0.5d0*s0*mn/(s/2d0)
      if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) a1min=a1cut
      a1max = comm + root
      a1cut = 1d0-0.5d0*s0*(n-mn)/(s/2d0)
      if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) a1max=a1cut
      ix = ix + 1
      call random_to_var(x(ix),-1d0,a1min,a1max,a1,dum)
      jaco = jaco*a1
      soft = soft*log(a1max/a1min)

      dummy=(/0d0,0d0,0d0,0d0/)
      dm=1d0
      a2=300d0
      E1 = (s+s1-s2)/(2d0*dsqrt(s))
      E2 = dsqrt(s) - E1
      solution = solver(s,s1,s2,dummy,dummy,dm,a1,a2,E1,dummy,dm) ! use same solver 

      Qm = (/E1,solution(1),solution(2),solution(3)/)
      Qnm = (/E2,-solution(1),-solution(2),-solution(3)/)
    endif
  end subroutine generate_split_Qm_Qnm


  subroutine generate_first_single(k,s,s1,q1_cmf,P_cmf,a1,s2)
    implicit none
    real(kind=8) :: s,s1
    integer :: k
    real(kind=8),intent(out) :: a1,s2
    real(kind=8) :: Lambda,Sigma,Delta,sigmak,smin,smax,smax_force
    real(kind=8) :: A,B,C,R,gs,a2,mu,a1min,a1max,dum
    real(kind=8),dimension(0:3) :: q1_cmf,P_cmf

    Lambda = (tot_mass-mass_sum-s1)+(k-1)*(k-2)*s0/2D0
    Sigma = (tot_mass-mass_sum-s1) 
    Delta = s1+2D0*(k-1)*s0/2D0
    sigmak =  s1  
    ! NOTE: added extra upper limit for massive case!
    if (Delta .lt. (2d0*dsqrt(s1*s)-s1)) then
            Delta = (2d0*dsqrt(s1*s)-s1)
    endif
    smax_force = s*(1-(k)*(s0/2d0)/dot(q1_cmf,P_cmf))
    smin = Lambda
    smax = s - Delta

    A = Sigma
    if (k.gt.2) then
      if (.not. flat) then
        ix = ix +1
        call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
        C = ((smax-A)/(smin-A))**R
        s2 = (smin-A)*C + A
        soft = soft*(log(smax-A)-log(smin-A))
        gs = 1d0/(s2-Sigma)
        jaco = jaco/gs
        a1 = a1_m1(s,s1,s2,q1_cmf,P_cmf,k)
        write(14,*) a1
        mu = (s2-s1)/s
        a2 = a1 + mu
        jaco = jaco*a1*(1d0-a1)*(1d0-a2)*a2
      elseif (flat) then
        ix = ix +1
        call random_to_var(x(ix),0d0,smax,smin,s2,dum)
        soft = soft*(smax-smin)
        gs = 1d0
        jaco = jaco/gs
        a1max = 0.5d0*(1d0+(s1-s2)/(s)+dsqrt(kallen(1d0,s1/s,s2/s)))
         a1max = a1max-0.00001d0
         a1min = 0.5d0*(1d0+(s1-s2)/(s)-dsqrt(kallen(1d0,s1/s,s2/s)))
         if ((a1min .lt. (s0/2d0)/(dot(q1_cmf,P_cmf)) ) .and.&
              a1max .gt. (s0/2d0)/(dot(q1_cmf,P_cmf))) then
              a1min = (s0/2d0)/(dot(q1_cmf,P_cmf))
         endif
         ix = ix +1
         call random_to_var(x(ix),0d0,a1min,a1max,a1,dum)
         soft = soft*(a1max-a1min)
      endif

    elseif (k.eq.2) then
         a1 = a1_m1(s,s1,s2,q1_cmf,P_cmf,k)
         jaco = jaco*a1 *(1d0-a1)
    endif
   
  end subroutine generate_first_single

  subroutine generate_s2(k,s,s1,s2,q1_cmf,P_cmf)
    implicit none
    integer :: k
    real(kind=8),intent(out) :: s2
    real(kind=8),dimension(0:3) :: q1_cmf,P_cmf
    real(kind=8) :: s,s1,A,B,C,dum,R,gs
    real(kind=8) :: Lambda,Delta,Sigma,sigmak,Sigmaold,smin,smax,smax_force
    double precision :: scut,inv

    double precision :: p1,p2,p3,bb

    ! Include also initial momenta in the limits!
   
    scut = s0
    Lambda = mass_sum+(k-1)*(k-2)*scut/2D0
    Sigma = mass_sum    !always just a sum of particle masses
    Sigmaold = mass_sum-s1
    Delta = s1+2d0*(k-1)*scut/2D0
    sigmak =  s1   !always just the previous particle mass
    ! NOTE: added extra upper limit for massive case!
    if (Delta .lt. (2d0*dsqrt(s1*s)-s1)) then
            Delta = (2d0*dsqrt(s1*s)-s1)
    endif
    smin = Lambda
    smax = s - Delta
    smax_force = s*(1-(s0/2d0)/dot(q1_cmf,P_cmf))
    if (smax .gt. smax_force) then
        smax = smax_force
    endif
    A = Sigma
    B = s - sigmak

    ! S limits exactly same as in COMIX! 

    if (smin.gt.smax) return
    if ((.not.open) .and. (.not. flat)) then
       ix = ix +1
       call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
       C = ((smax - A)*(B-smin)/((B-smax) * (smin - A)))**R
       s2 = (A*(B - smin) + B*(smin - A)*C)/(B - smin + (smin-A)*C)
       soft = soft*(log((smax-Sigma)/(s-sigmak-smax))-&
                    log((smin-Sigma)/(s-sigmak-smin)))
       gs = (s-Sigmaold)/((s-sigmak-s2)*(s2-Sigma))
       jaco = jaco/gs
    endif

    if (flat) then
       ix = ix +1
       call random_to_var(x(ix),0d0,smin,smax,s2,dum)
       soft = soft*(smax-smin)
       gs = 1d0
       jaco = jaco/gs
    endif
    if (open) then
       ix = ix +1
       call random_to_var(x(ix),-1d0,smin,smax,s2,dum)
       soft = soft*log(smax/smin)
       jaco = jaco*s2
    endif
  end subroutine generate_s2


  subroutine generate_a1_term1(i,m1,maxn,s,s1,s2,cos,a1cut,beta,h_in,a1)
    implicit none
    integer :: i,maxn
    real(kind=8),intent(out) :: a1
    real(kind=8) :: s,s1,s2,cos
    logical :: m1
    real(kind=8) :: a1cut
    real(kind=8) :: h1,h,h_in,a1min,a1max,a1max_force,f_h1,beta,dum
    real(kind=8),dimension(3) :: buff
    real(kind=8) :: Amin,Amax,Atilde,R,v,wsq,kappa
    double precision :: low,upp,E,dum1,dum2

    h=h_in
    h1 = ((1d0-beta)/(2d0*beta))*(1d0+(s1-s2)/s)
    a1max = 0.5d0*(1d0-(s2-s1)/(s)+dsqrt(kallen(1d0,s1/s,s2/s)))
    a1min = 0.5d0*(1d0-(s2-s1)/(s)-dsqrt(kallen(1d0,s1/s,s2/s)))
    if (a1min .lt. 0d0) then ! just for numerical stability!
        a1min = 0d0
    endif
    if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) then
        a1min = a1cut
    endif
    if (abs(h).gt.0) then
      a1max_force = 1d0 - a1cut*(maxn-i-1)
      if ((a1max_force.lt.a1max).and.(a1max_force.gt.a1min)) then
        a1max=a1max_force
      endif
    endif

    buff = f_func_term1(-h1,cos,s,s1,s2,h)
    f_h1 = buff(2)
    buff = f_func_term1(a1min,cos,s,s1,s2,h)
    Amin= (a1min+ h1 + buff(2)-f_h1) &
         /(a1min+ h1 + buff(2) +f_h1)
    buff = f_func_term1(a1max,cos,s,s1,s2,h)
    Amax= (a1max + h1 + buff(2)-f_h1) &
          /(a1max + h1 + buff(2) +f_h1)

    ! Now generate a1, depending on which distribution
    if ((.not. open) .and. (.not. flat)) then
      if (( ((i .ne. 0) .and. (.not. m1)) .or. (maxn .ne. n))) then
        ! Sample with Pi^{-1/2} factor
        ix = ix + 1
        call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
        buff = f_func_term1(a1,cos,s,s1,s2,h)
        v = buff(1)/2d0
        buff = f_func_term1(0d0,cos,s,s1,s2,h)
        wsq = buff(2)**2  ! w^2
        if (Amin .le. 0d0) then
          Amin = Amax*1d-8 ! for numerical stability
        endif
        Atilde = (Amin**(1d0-R))*(Amax)**R
        kappa = - h1 + f_h1*(1d0+Atilde)/(1d0-Atilde)
        buff = f_func_term1(a1,cos,s,s1,s2,h)
        if (v+kappa.eq.0d0) then
           a1=a1max
        else
           a1 = ((kappa**2) - wsq)/(2d0*(v+kappa))
        endif
        if ((a1min-a1)/a1min.le.1d-6.and.a1min-a1.gt.0d0) a1=a1min 
        if ((a1-a1max)/a1max.le.1d-6.and.a1-a1max.gt.0d0) a1=a1max
        soft = soft*(1d0/f_h1)*(log(Amax)-log(Amin))
        jaco = jaco*a1
        jaco = jaco*beta
      elseif (((i .eq. 0) .or. (m1 .and. (i .le. 1)))) then
        ! Sample with 1/x 
        ix = ix + 1
        call random_to_var(x(ix),-1d0,a1min,a1max,a1,dum)
        soft = soft*log(a1max/a1min)
        jaco = jaco*a1
      endif
      if ((a1 .gt. a1max) .or. (a1 .lt. a1min)) then
              num_error = num_error+1
      endif
    endif

    if ((flat)) then
      ix = ix +1
      call random_to_var(x(ix),0d0,a1min,a1max,a1,dum)
      soft = soft*(a1max-a1min)
    endif
    if (open) then
        ix = ix + 1
        call random_to_var(x(ix),-1d0,a1min,a1max,a1,dum)
        soft = soft*log(a1max/a1min)
        jaco = jaco*a1
    endif
  end subroutine generate_a1_term1

  subroutine generate_a1_term2(i,m1,maxn,a2,s,s1,s2,cos,a1cut,a1)
    implicit none
    integer :: i,maxn
    real(kind=8) :: a2,s,s1,s2,cos,h,a1cut
    logical :: m1
    real(kind=8),intent(out) :: a1
    real(kind=8),dimension(2) :: a1pm
    real(kind=8) :: a1maxbar,a1minbar,R,xy,a1max,a1min,dum

    a1pm = a1_pm(a2,s,s1,s2,cos)
    a1min = a1pm(2)
    a1max = a1pm(1)
    if ((a1cut.lt.a1pm(1)).and.(a1cut.gt.a1pm(2)))  then
       a1min = a1cut
    endif
    h=a2
    a1maxbar = a1max + h
    a1minbar = a1min + h

    ! Now generate a1
    ix= ix + 1
    call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
    xy = tan(-pi/2d0 * R)**2
    a1 = a1maxbar*a1minbar*(1d0+xy)/(a1minbar + xy*a1maxbar) - h
    if ((a1min-a1)/a1min.le.1d-8.and.a1min-a1.gt.0d0) a1=a1min
    if ((a1-a1max)/a1max.le.1d-8.and.a1-a1max.gt.0d0) a1=a1max
    jaco = jaco*a1
  end subroutine generate_a1_term2  

  subroutine generate_a2_term1(i,m1,maxn,a1,s,s1,s2,cos,a2cut,h_in,a2)
    implicit none
    integer :: i,maxn
    real(kind=8) :: a1,s,s1,s2,cos,h,h_in,a2cut
    logical :: m1
    real(kind=8),intent(out) :: a2
    real(kind=8),dimension(2) :: a2pm
    real(kind=8) :: a2maxbar,a2minbar,R,xy,a2max,a2min,dum,a2max_force

    a2pm = a2_pm(a1,s,s1,s2,cos)
    a2min = a2pm(2)
    a2max = a2pm(1)

    h=h_in
    if (h_in.eq.0d0) then
       h = a1
    elseif (h_in.lt.0d0) then
       h = 0d0
    endif

    if (h_in.le.0d0) then
     if (maxn-i-1.le.1) then
      if ((a2cut.lt.a2pm(1)).and.(a2cut.gt.a2pm(2)))  then
       a2minbar = a2cut + h
       a2min = a2cut
      endif
     endif
     a2max_force = 1d0-a2cut/(maxn-i-1)
     if ((a2max_force.lt.a2max).and.(a2max_force.gt.a2min)) then
       a2maxbar = a2max_force + h
       a2max = a2max_force
     endif
    endif

    a2maxbar = a2max + h
    a2minbar = a2min + h

    ! Pretty much same as COMIX, except h stuff

    ! Now generate a2
    if ((.not. open) .and. (.not. flat)) then
    if ((m1 .and. (i .ge. 2)) .or.&
            ((.not. m1) .and. (maxn .eq. n) .and. (i .ge. 1))&
            .or. (maxn .ne. n))  then
        ix = ix + 1
        call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
        xy = tan(-pi/2d0 * R)**2
        a2 = a2maxbar*a2minbar*(1d0+xy)/(a2minbar + xy*a2maxbar) - h
        if ((a2min-a2)/a2min.le.1d-8.and.a2min-a2.gt.0d0) a2=a2min
        if ((a2-a2max)/a2max.le.1d-8.and.a2-a2max.gt.0d0) a1=a2max
        jaco = jaco*a2
        soft = soft*(pi/2d0)
    elseif( ((i .eq. 0) .and. (maxn .eq. n)) .or. ((m1 .and. (i .le. 1)))) then
     ! Do phi-integration instead of a2
      a2 = 300d0          ! dummy value to do phi-integration
    endif
    endif

    if (flat) then
      if ((m1 .and. (i .ge. 2)) .or.&
            ((.not. m1) .and. (maxn .eq. n) .and. (i .ge. 1))&
            .or. (maxn .ne. n))  then
        ix = ix +1
        call random_to_var(x(ix),0d0,a2min,a2max,a2,dum)
        soft = soft*(a2max-a2min)
        jaco = jaco/dsqrt(4d0*(a2max-a2)*(a2-a2min))
      elseif( ((i .eq. 0) .and. (maxn .eq. n)) .or. ((m1 .and. (i .le. 1)))) then
        ! Do phi-integration instead of a2
        a2 = 300d0          ! dummy value to do phi-integration
      endif
    endif
    if (open) then
      a2 = 300d0          ! dummy value to do phi-integration
    endif

  end subroutine generate_a2_term1

  subroutine generate_a2_term2(i,m1,maxn,s,s1,s2,cos,a2cut,beta,a2)
    implicit none
    integer :: i,maxn
    real(kind=8),intent(out) :: a2
    real(kind=8) :: a1,s,s1,s2,cos
    logical :: m1
    real(kind=8) :: a2cut
    real(kind=8) :: h1,a2min,a2max,f_h1,beta,dum
    real(kind=8),dimension(3) :: buff
    real(kind=8) :: Amin,Amax,Atilde,R,v,wsq,kappa
    double precision :: low,upp,E,dum1,dum2

    double precision t1,t2,aa,sin

    h1=0d0
    a2max = 0.5d0*((1d0+(s2-s1)/s)+dsqrt(kallen(1d0,s1/s,s2/s)))
    a2min = 0.5d0*((1d0+(s2-s1)/s)-dsqrt(kallen(1d0,s1/s,s2/s)))
    if ((a2cut.gt.a2min).and.(a2cut.lt.a2max)) then
        a2min=a2cut
    endif

    buff = f_func_term2(-h1,cos,s,s1,s2)
    f_h1 = buff(2)
    buff = f_func_term2(a2min,cos,s,s1,s2)
    Amin= (a2min+ h1 + buff(2)-f_h1) &
         /(a2min+ h1 + buff(2) +f_h1)
    buff = f_func_term2(a2max,cos,s,s1,s2)
    Amax= (a2max + h1 + buff(2)-f_h1) &
          /(a2max + h1 + buff(2) +f_h1)

    ! Sample with Pi^{-1/2} factor
    ix = ix + 1
    call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
    buff = f_func_term2(0d0,cos,s,s1,s2)
    wsq = buff(3)  ! w^2
    if (Amin .le. 0d0) then
       Amin = Amax*1d-8
    endif
    Atilde = (Amin**(1d0-R))*(Amax)**R
    kappa = - h1 + f_h1*(1d0+Atilde)/(1d0-Atilde)
    buff = f_func_term2(a2,cos,s,s1,s2)
    a2 = ((kappa**2) - wsq)/(1d0*buff(1)+2d0*kappa)
    if ((a2min-a2)/a2min.le.1d-8.and.a2min-a2.gt.0d0) a2=a2min
    if ((a2-a2max)/a2max.le.1d-8.and.a2-a2max.gt.0d0) a2=a2max
    soft = soft*(1d0/f_h1)*(log(Amax)-log(Amin))
    jaco = jaco*a2
    jaco = jaco*beta
    if ((a2 .gt. a2max) .or. (a2 .lt. a2min)) then
              num_error = num_error+1
    endif
  end subroutine generate_a2_term2

  function f_func_term1(a1,cos,s,s1,s2,h)
    implicit none
    real(kind=8),dimension(3) :: f_func_term1
    real(kind=8) :: beta, a1,cos,sin,s,s1,s2,b,lin_coeff,con_term,h
    sin = dsqrt(1d0-cos**2)
    beta = (1d0+(s2-s1)/s) + (1d0-(s2-s1)/s)*cos
    if (h.gt.0d0) then ! Old approach
      lin_coeff = -cos*beta-sin**2+sin**2*(s2-s1)/s-2d0*h*cos
      con_term = 0.25d0*(beta**2) + h*beta + h**2 + (sin**2)*s1/s
    elseif (h.eq.0d0) then ! New approach
      lin_coeff = (-cos*beta-sin**2+sin**2*(s2-s1)/s + beta)/(2d0-2d0*cos)
      con_term = (0.25d0*(beta**2) + (sin**2)*s1/s)/(2d0-2d0*cos)
    elseif (h.lt.0d0) then
      lin_coeff = -cos*beta-sin**2+sin**2*(s2-s1)/s
      con_term = 0.25d0*(beta**2) + (sin**2)*s1/s
    endif
    f_func_term1 = (/lin_coeff,dsqrt(a1**2 + lin_coeff*a1 + con_term),con_term/)
  end function f_func_term1

  function f_func_term2(a2,cos,s,s1,s2)
    implicit none
    real(kind=8),dimension(3) :: f_func_term2
    real(kind=8) :: A,B, a2,cos,sin,s,s1,s2,lin_coeff,con_term,beta
    sin = dsqrt(1d0-cos**2)
    beta = (1d0-(s2-s1)/s)+(1d0+(s2-s1)/s)*cos
    lin_coeff = (-cos*beta - sin**2*(1d0 + (s2-s1)/s) + beta)/(2d0-2d0*cos)
    con_term = (0.25d0*(beta**2) + (sin**2)*s2/s)/(2d0-2d0*cos)
    f_func_term2 = (/lin_coeff,dsqrt(a2**2 + lin_coeff*a2 + con_term),con_term/)
  end function f_func_term2

  function a2_pm(a1,s,s1,s2,c)
    ! cross-checked with pi function!
    implicit none
    real(kind=8) :: comm,s1,s2,s,a1,a2plus,a2minus,c
    real(kind=8), dimension(2) :: a2_pm
    comm = 0.5*(1.+((s2-s1)/s)+c*(1.-2.*a1-(s2-s1)/s))
    a2plus  = comm + dsqrt(1.-c**2)*dsqrt(a1*(1.-a1-(s2-s1)/s)-s1/s)
    a2minus = comm - dsqrt(1.-c**2)*dsqrt(a1*(1.-a1-(s2-s1)/s)-s1/s)
    a2_pm = (/a2plus,a2minus/)
  end function a2_pm

  function a1_pm(a2,s,s1,s2,c)
    ! cross-checked with pi function!
    implicit none
    real(kind=8) :: comm,s1,s2,s,a2,a1plus,a1minus,c
    real(kind=8), dimension(2) :: a1_pm
    comm = 0.5*(1.-(s2-s1)/s+c*(1.-2.*a2+(s2-s1)/s))
    a1plus  = comm + dsqrt(1.-c**2)*dsqrt(a2*(1.-a2+(s2-s1)/s)-s2/s)
    a1minus = comm - dsqrt(1.-c**2)*dsqrt(a2*(1.-a2+(s2-s1)/s)-s2/s)
    a1_pm = (/a1plus,a1minus/)
  end function a1_pm

  real(kind=8) function pi_func(a1,a2,s1,s2,s,cos)
      ! note: slightly different than in haag paper (same as in COMIX)
      implicit none
      real(kind=8) :: a1,a2,s1,s2,s,cos,sin
      sin=dsqrt(1.-cos**2)
      pi_func = 4.*sin**2*((1.-a2+s2/s-s1/s)*a2 - s2/s)- &
      (1.-2.*a1+s1/s-s2/s + cos*(1.-2.*a2+s2/s-s1/s))**2
  end function pi_func

  function solver(s,s1,s2,q1_cmf,q2_cmf,z_sign,a1,a2,E,P_cmf,co)
    implicit none
    real(kind=8),dimension(3) :: solver
    real(kind=8) z_sign
    real(kind=8) :: s,s1,s2,a1,a2,E,co
    real(kind=8), dimension(0:3) :: q1_cmf,q2_cmf,P_cmf,q2_rot,q1_rot,q2_rot_x0,q2_rot_x0_new
    real(kind=8), dimension(0:3) :: p1,p2,p1_rot,p2_rot,p1_xy_rot,p2_xy_rot
    real(kind=8), dimension(0:3) :: test1,test2,test3
    real(kind=8) :: a_1,a_2,b_1,b_2,c_1,c_2,e_cmf,r1,r2,pmag
    real(kind=8) :: A,B,P,R,T,Q,t1,t2,t3,xxx,y,z,dum,phi,sgn
    logical,parameter :: old_mapping=.false.
    real(kind=8),external :: ran2
          
    if (a2 .gt. 100d0) then
        z = E - dsqrt(s)*a1
        ix = ix + 1
        call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,dum)
        soft = soft*(pi/2d0)
        xxx =  dsqrt(E**2 - s1 - z**2)*cos(phi)
        y = dsqrt(E**2 - s1 - z**2)*sin(phi)
        p1=(/E,xxx,y,z/)
        p2=(/dsqrt(s)-E,-xxx,-y,-z/)
        call rotxxx(p1,q1_cmf,p1_rot)
        call rotxxx(p2,q1_cmf,p2_rot)
    else
        z = E - dsqrt(s)*a1
        y = (-dsqrt(s)+E+dsqrt(s)*a2-co*z)/dsqrt(1d0-co**2)
        R = RAN2()
        if (R .lt. 0.5d0) then
           sgn = 1d0
        else
           sgn=-1d0
        endif

        if (E**2-s1-y**2-z**2 .le. 0d0) then
           xxx = 0d0
        else
           xxx = sgn*dsqrt(E**2-s1-y**2-z**2)
        endif
        p1 = (/E,xxx,y,z/)

        p2 = (/sqrt(s)-E,-xxx,-y,-z/)
        e_cmf = dsqrt(s)/2d0
        r1 = dot(P_cmf,P_cmf)/(2d0*dot(P_cmf,q1_cmf))
        r2 = dot(P_cmf,P_cmf)/(2d0*dot(P_cmf,q2_cmf))
        ! q1,q2 in the frame in which p1,p2 given above
        q1_rot = (/e_cmf/r1,0d0,0d0,e_cmf/r1/)
        q2_rot_x0=(/e_cmf/r2,0d0,e_cmf/r2*dsqrt(1d0-co**2),e_cmf/r2*co/) ! temporary
        call rot_to_z(q1_cmf,q2_cmf,q2_rot,1d0)
        call rot_xy_plane(q2_rot,q2_rot,q2_rot_x0,q2_rot_x0_new,-1d0)
        call rot_xy_plane(p1,q2_rot,q2_rot_x0_new,p1_xy_rot,1d0)
        call rot_xy_plane(p2,q2_rot,q2_rot_x0_new,p2_xy_rot,1d0)
        call rot_to_z(q1_cmf,p1_xy_rot,p1_rot,-1d0)
        call rot_to_z(q1_cmf,p2_xy_rot,p2_rot,-1d0)

    endif

    solver = (/p1_rot(1),p1_rot(2),p1_rot(3)/)
  end function solver

  subroutine rot_to_z(q1,q2,q2_rot,sgnin)
         implicit none
         real(kind=8),dimension(0:3) :: q1,q2
         real(kind=8),dimension(0:3),intent(out) :: q2_rot
         real(kind=8),dimension(1:3) :: q1_s,q2_s,z,n,p,t,tp
         real(kind=8),dimension(1:3) :: t_pz,t_pz_p
         real(kind=8) :: costheta,sintheta,sgn,sgnin
        
         z = (/0d0,0d0,1d0/)
         q1_s(1:3) = q1(1:3)
         q2_s(1:3) = q2(1:3)

         n = threecross(q1_s,z) 
         n = n/dsqrt(threedot(n,n)) ! Normal to q1-z plane 
         t = q2_s - threedot(n,q2_s)*n ! Component of q2 in q1-z plane

         p = q1_s - threedot(q1_s,z)*z
         p = p/dsqrt(threedot(p,p)) ! Orthogonal vector to z in q1-z plane

         ! Decompone t to p- and z-components
         t_pz(1) = threedot(t,p)
         t_pz(2) = threedot(t,z)
         t_pz(3) = 0d0

         ! Angle of rotation in q1-z plane
         costheta=threedot(q1_s,z)/dsqrt(threedot(q1_s,q1_s))
         sintheta=dsqrt(1d0-costheta**2)

         if ((q1(2) .ge. 0d0) .and. (q2(2) .ge. 0d0)) then
            sgn = sgnin*1d0
         elseif ((q1(2) .lt. 0d0) .and. (q2(2) .lt. 0d0)) then
            sgn = sgnin*1d0
         else
            sgn = sgnin*1d0
         endif

         t_pz_p(1) = costheta*t_pz(1) + sgn*(-sintheta)*t_pz(2)
         t_pz_p(2) = sgn*sintheta*t_pz(1) + costheta*t_pz(2)
         t_pz_p(3) = t_pz(3) ! Dummy variable

         q2_rot(1:3) = threedot(n,q2_s)*n + t_pz_p(1)*p + t_pz_p(2)*z ! Rotated q2
         q2_rot(0) = q2(0)
  end subroutine rot_to_z

  function a1_m1(s,s1,s2,q1,P,k)
    implicit none
    real(kind=8) :: s,s1,s2,RHS
    real(kind=8) :: mu,a1max,a1min,sum_w,R,dum,sgn,a1cut
    real(kind=8),dimension(4) :: g1,g2,d,e,w
    real(kind=8),dimension(0:3) :: q1,P
    real(kind=8) :: a1_m1
    integer :: i,pick,k
    real(kind=8),external :: ran2
    real(kind=8) :: w1,w2,w_tot,a1,A

    mu=(s2-s1)/s
    a1max = 0.5d0*(1d0+(s1-s2)/(s)+dsqrt(kallen(1d0,s1/s,s2/s)))
    a1min = 0.5d0*(1d0+(s1-s2)/(s)-dsqrt(kallen(1d0,s1/s,s2/s)))

    a1cut = (s0/2d0)/(dot(q1,P))

    if ((a1min.lt.a1cut ).and.(a1max.gt.a1cut)) then
        a1min = a1cut
    endif

    !if ((1d0-(k-1)*a1cut-mu.gt.a1min).and.(1d0-(k-1)*a1cut-mu.lt.a1max)) then
    !     a1max = 1d0-(k-1)*a1cut-mu
    !endif

    a1max = a1max-1d-8  !! This is a hard cut, but it works

    if (k.gt.2) then
      g1 = (/1d0,1d0,1d0,-1d0/)
      g2 = (/1d0,-1d0,-1d0,-1d0/)
      e = (/mu, 1d0-mu, 1d0, 1d0-mu/)
      d = (/0d0,0d0,mu,1d0/)

      do i=1,4
        w(i) = 1d0/(g2(i)*e(i)-g1(i)*d(i))*&
          (log((g1(i)*a1max+d(i))/(g2(i)*a1max+e(i)))-&
           log((g1(i)*a1min+d(i))/(g2(i)*a1min+e(i))))
      enddo
      ! Correction 
      w(2) = -w(2)
      w(3) = -w(3)
      sum_w = sum(w)
  
      R = ran2()
      if (R .lt. w(1)/sum_w) then
         pick = 1
      elseif ((R .lt. (w(1)+w(2))/sum_w) .and. (R .gt. w(1)/sum_w)) then
         pick = 2
      elseif ((R .lt. (w(1)+w(2)+w(3))/sum_w) .and. (R .gt. (w(1)+w(2))/sum_w)) then
         pick = 3
      elseif (R .gt. (w(1)+w(2)+w(3))/sum_w) then
         pick = 4
      endif
      ix = ix +1
      call random_to_var(x(ix),0d0,0d0,1d0,R,dum)
      sgn = 1d0
      if ((pick .eq. 2) .or. (pick .eq. 3)) then
            sgn = -1d0
      endif
      RHS = exp(w(pick)*R*sgn*(g2(pick)*e(pick)-g1(pick)*d(pick)))*&
            (g1(pick)*a1min+d(pick))/(g2(pick)*a1min+e(pick))
      a1_m1 = (e(pick)*RHS-d(pick))/(g1(pick)- RHS*g2(pick))
      if ((a1min-a1_m1)/a1min.le.1d-8.and.a1min-a1_m1.gt.0d0) a1_m1=a1min
      if ((a1_m1-a1max)/a1max.le.1d-8.and.a1_m1-a1max.gt.0d0) a1_m1=a1max
      soft = soft*sum_w

    elseif (k.eq.2) then
      a1min = (s0/2d0)/(dot(q1,P))
      a1max = 1d0
      if ((1d0-(s0/2d0)/(dot(q1,P)) .lt.a1max).and.(1d0-(s0/2d0)/(dot(q1,P)).gt.a1min)) then
         a1max = 1d0-(s0/2d0)/(dot(q1,P))
      endif
      R=ran2()
      w1= log(a1max/a1min)
      w2 = log((1d0-a1min)/(1d0-a1max))
      w_tot=w1+w2

      if (R.lt.(w1/w_tot)) pick=1
      if (R.gt.(w1/w_tot)) pick=2  

      if (pick.eq.1) then
         ix = ix +1
         call random_to_var(x(ix),-1d0,a1min,a1max,a1,dum)
         a1_m1 = a1
         if ((a1min-a1_m1)/a1min.le.1d-8.and.a1min-a1_m1.gt.0d0) a1_m1=a1min
         if ((a1_m1-a1max)/a1max.le.1d-8.and.a1_m1-a1max.gt.0d0) a1_m1=a1max
      elseif (pick.eq.2) then
         ix = ix +1
         call random_to_var(x(ix),-1d0,1d0-a1max,1d0-a1min,a1,dum)
         a1_m1 = a1
         if ((a1min-a1_m1)/a1min.le.1d-8.and.a1min-a1_m1.gt.0d0) a1_m1=a1min
         if ((a1_m1-a1max)/a1max.le.1d-8.and.a1_m1-a1max.gt.0d0) a1_m1=a1max
      endif
      soft = soft*w_tot

    endif

  end function a1_m1

  function polylog(x,n)
    implicit none
    real(kind=8) :: x,polylog
    integer :: i,n
    polylog = 0d0
    do i=1,n
      polylog = polylog + ((1d0-x)**i)/(i**2)
    enddo
  end function polylog

  function threecross(a,b)
       implicit none
       real(kind=8),dimension(1:3) :: a,b,c
       real(kind=8),dimension(1:3) :: threecross
       c(1) = a(2)*b(3)-a(3)*b(2)
       c(2) = a(3)*b(1)-a(1)*b(3)
       c(3) = a(1)*b(2)-a(2)*b(1)
       threecross=c
  end function

  FUNCTION DDILOG(X)
!*
!* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
!*
!* $Log: imp64.inc,v $
!* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
!* Mathlib gen
!*
!*
!* imp64.inc
!*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:19)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159265358979324D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)
      DATA C( 0) / 0.42996693560813697D0/
      DATA C( 1) / 0.40975987533077105D0/
      DATA C( 2) /-0.01858843665014592D0/
      DATA C( 3) / 0.00145751084062268D0/
      DATA C( 4) /-0.00014304184442340D0/
      DATA C( 5) / 0.00001588415541880D0/
      DATA C( 6) /-0.00000190784959387D0/
      DATA C( 7) / 0.00000024195180854D0/
      DATA C( 8) /-0.00000003193341274D0/
      DATA C( 9) / 0.00000000434545063D0/
      DATA C(10) /-0.00000000060578480D0/
      DATA C(11) / 0.00000000008612098D0/
      DATA C(12) /-0.00000000001244332D0/
      DATA C(13) / 0.00000000000182256D0/
      DATA C(14) /-0.00000000000027007D0/
      DATA C(15) / 0.00000000000004042D0/
      DATA C(16) /-0.00000000000000610D0/
      DATA C(17) / 0.00000000000000093D0/
      DATA C(18) /-0.00000000000000014D0/
      DATA C(19) /+0.00000000000000002D0/
      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG=H
      RETURN
  END

  subroutine  rot_xy_plane(p,q,q_x0,p_rot,sgnin)
!     input:
!          p(0:3) : four-momentum in the frame q_x0(1) = 0
!          q(0:3) : four-momentum to which the q_x0 is rotated to
!          q_x0(0:3) : four-momentum which is rotated to q
!          sgnin : determines which direction to rotate
!      output:
!          p_rot(0:3) : four-momentum p rotated from p according to 
!                       q_x0 -> q rotation
          implicit none
          real(kind=8),dimension(0:3) :: p,q,q_x0,p_rot
          real(kind=8),dimension(1:3) :: q_xy,q_x0_xy
          real(kind=8) :: sgn,costheta,sintheta,sgnin

          q_xy(1) = q(1)
          q_xy(2) = q(2)
          q_xy(3) = 0d0
          
          q_x0_xy(1) = q_x0(1)
          q_x0_xy(2) = q_x0(2)
          q_x0_xy(3) = 0d0

          if (q(1) .ge. 0d0) then
             sgn = sgnin*(-1d0)
          else
             sgn = sgnin*1d0
          endif

          costheta = threedot(q_xy,q_x0_xy)/&
                  (dsqrt(threedot(q_xy,q_xy))*dsqrt(threedot(q_x0_xy,q_x0_xy)))
          if (costheta .gt. 1d0) then
              if ((abs(costheta)-1d0) .gt. 0.000001d0) then
                      write(*,*) 'ERROR in rot_xy_plane',costheta
                      stop 2
              endif
              costheta = 1d0
          elseif (costheta .lt. -1d0)then
              if ((abs(costheta)-1d0) .gt. 0.000001d0) then
                      write(*,*) 'ERROR in rot_xy_plane',abs(costheta-1d0)
                      stop 2
              endif
              costheta = -1d0
          endif
          sintheta=dsqrt(1d0-costheta**2)

          p_rot(1) = costheta*p(1) + sgn*(-sintheta)*p(2)
          p_rot(2) = sgn*sintheta*p(1) + costheta*p(2)
          p_rot(3) = p(3)
          p_rot(0) = p(0)
  end subroutine

  subroutine check_momenta(p,mass)
     implicit none
     real(kind=8), dimension(0:3,next) :: p
     real(kind=8), dimension(0:3) :: tot_mom
     real(kind=8), dimension(next) :: mass
     integer i 
     real(kind=8) :: curr_mass

     tot_mom = (/0d0,0d0,0d0,0d0/)
     do i=1,next
       curr_mass = dot(p(0:3,i),p(0:3,i))
       if (abs(curr_mass - masses(i)**2) .gt. 10d0**(-9d0)) then
          write(*,*) 'ERROR in mass!',abs(curr_mass - masses(i)**2)
          write(*,*) abs(curr_mass - masses(i)**2) .gt. 10d0**(-9d0)
       endif
       tot_mom = tot_mom+p(:,i)
     enddo

     if (abs(tot_mom(0)-2d0*sqrtshat)/2d0*sqrtshat .gt.0.0001d0) then
        write(*,*) 'ERROR in energy conservation!'
     endif
     do i=1,3
        if (abs(tot_mom(i)).gt. 0.0000010) then
                write(*,*) 'ERROR in momentum conservation!'
        endif
     enddo
  end subroutine check_momenta

  subroutine rotxxx(p,q,prot)  ! from HELAS library
! This subroutine performs the spacial rotation of a four-momentum.
! the momentum p is assumed to be given in the frame where the spacial
! component of q points the positive z-axis.  prot is the momentum p
! rotated to the frame where q is given.
! input:
!       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame
!       real    q(0:3)         : four-momentum q in the rotated frame
! output:
!       real    prot(0:3)      : four-momentum p in the rotated frame
    implicit none
    real(kind=8),dimension(0:3) :: p,q
    real(kind=8),dimension(0:3) :: prot
    real(kind=8) :: qt2,qt,psgn,qq
    real(kind=8),parameter :: rZero=0d0,rOne=1d0
    prot(0) = p(0)
    qt2 = q(1)**2 + q(2)**2
    if ( qt2.eq.rZero ) then
       if ( q(3).eq.rZero ) then
          prot(1:3) = p(1:3)
       else
          psgn = sign(rOne,q(3))
          prot(1:3) = p(1:3)*psgn
       endif
    else
       qq = sqrt(qt2+q(3)**2)
       qt = sqrt(qt2)
       prot(1) = q(1)*q(3)/qq/qt*p(1) -q(2)/qt*p(2) +q(1)/qq*p(3)
       prot(2) = q(2)*q(3)/qq/qt*p(1) +q(1)/qt*p(2) +q(2)/qq*p(3)
       prot(3) =          -qt/qq*p(1)               +q(3)/qq*p(3)
    endif
  end subroutine rotxxx

  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot

  real(kind=8) function threedot(p1,p2)
    ! Inner product between two 3-vectors
    implicit none
    real(kind=8),intent(in),dimension(1:3) :: p1,p2
    threedot=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
  end function threedot

  subroutine boostm(p,q,m,pboost)  ! from HELAS library
! This subroutine performs the Lorentz boost of a four-momentum.  The
! momentum p is assumed to be given in the rest frame of q.  pboost is
! the momentum p boosted to the frame in which q is given.  q must be a
! timelike momentum.
! input:
!       real    p(0:3)         : four-momentum p in the q rest  frame
!       real    q(0:3)         : four-momentum q in the boosted frame
!       real    m              : mass of q (for numerical stability)
! output:
!       real    pboost(0:3)    : four-momentum p in the boosted frame
    implicit none
    real(kind=8) ,dimension(0:3) :: p,q,pboost
    real(kind=8) :: pq,qq,m,lf
    real(kind=8),parameter :: rZero=0d0
    qq = q(1)**2+q(2)**2+q(3)**2
    if ( qq.ne.rZero ) then
       pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
       lf = ((q(0)-m)*pq/qq+p(0))/m
       pboost(0) = (p(0)*q(0)+pq)/m
       pboost(1:3) =  p(1:3)+q(1:3)*lf
    else
       pboost(0:3) = p(0:3)
    endif
  end subroutine boostm

  subroutine random_to_var(x,power,var_min,var_max,var,jac)
    ! Given a random number x, it generates var in the range var_min
    ! <= var <= var_max according to var^(power). 'jac' is the
    ! corresponding Jacobian.
    implicit none
    real(kind=8),intent(in) :: x,power,var_min,var_max
    real(kind=8),intent(out) :: var
    real(kind=8),intent(inout) :: jac
    integer(kind=4) :: ip
    real(kind=8) :: varmin,varmax
    if (var_min.lt.0d0 .and. var_max.le.0d0) then
       varmin=-var_max
       varmax=-var_min
    elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (power.ne.0d0)) then
       write (*,*) 'ERROR: in random_to_var one of the two limits '/&
            &/'is negative',var_min,var_max
       stop 8
    else
       varmin=var_min
       varmax=var_max
    endif
    ip=nint(power)
    if (dble(ip).eq.power) then
       ! integer
       if (ip.eq.-1) then
          var=varmax**x * varmin**(1d0-x)
          jac=jac*var*log(varmax/varmin)
       elseif (ip.eq.-2) then
          var=1d0/( (1d0-x)/varmin + x/varmax )
          jac=jac*var**2*(varmax-varmin)/(varmax*varmin)
       elseif (ip.eq.0) then
          var=varmin+x*(varmax-varmin)
          jac=jac*(varmax-varmin)
       else
          var=(varmin**(1+ip)*(1d0-x)+varmax**(1+ip)*x)**(1d0/(1d0+power))
          jac=jac*(varmax**(1+ip)-varmin**(1+ip))* &
               (varmin**(1+ip)*(1d0-x)+varmax**(1+ip)*x)**(-power/(1d0+power))/ &
               (1d0+power)
       endif
    else
       var=(varmin**(1d0+power)*(1d0-x)+varmax**(1d0+power)*x)**(1d0/(1d0+power))
       jac=jac*(varmax**(1d0+power)-varmin**(1d0+power))* &
            (varmin**(1d0+power)*(1d0-x)+varmax**(1d0+power)*x)**(-power/(1d0+power))/&
            (1d0+power)
    endif
    if (var_min.le.0d0 .and. var_max.le.0d0) then
       var=-var
    endif
  end subroutine random_to_var


!!$  FUNCTION RN(IDUMMY)
!!$      REAL*8 RN,RAN
!!$      SAVE INIT
!!$      DATA INIT /1/
!!$      IF (INIT.EQ.1) THEN
!!$        INIT=0
!!$        CALL RMARIN(1802,9373)
!!$      END IF
!!$!*
!!$  10  CALL RANMAR(RAN)
!!$      IF (RAN.LT.1D-16) GOTO 10
!!$      RN=RAN
!!$!*
!!$   END
!!$
!!$   SUBROUTINE RANMAR(RVEC)
!!$!*     -----------------
!!$!* Universal random number generator proposed by Marsaglia and Zaman
!!$!* in report FSU-SCRI-87-50
!!$!* In this version RVEC is a double precision variable.
!!$      IMPLICIT REAL*8(A-H,O-Z)
!!$      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
!!$      COMMON/ RASET2 / IRANMR,JRANMR
!!$      SAVE /RASET1/,/RASET2/
!!$      UNI = RANU(IRANMR) - RANU(JRANMR)
!!$      IF(UNI .LT. 0D0) UNI = UNI + 1D0
!!$      RANU(IRANMR) = UNI
!!$      IRANMR = IRANMR - 1
!!$      JRANMR = JRANMR - 1
!!$      IF(IRANMR .EQ. 0) IRANMR = 97
!!$      IF(JRANMR .EQ. 0) JRANMR = 97
!!$      RANC = RANC - RANCD
!!$      IF(RANC .LT. 0D0) RANC = RANC + RANCM
!!$      UNI = UNI - RANC
!!$      IF(UNI .LT. 0D0) UNI = UNI + 1D0
!!$      RVEC = UNI
!!$   END
!!$
!!$   SUBROUTINE RMARIN(IJ,KL)
!!$!*     -----------------
!!$!* Initializing routine for RANMAR, must be called before generating
!!$!* any pseudorandom numbers with RANMAR. The input values should be in
!!$!* the ranges 0<=ij<=31328 ; 0<=kl<=30081
!!$      IMPLICIT REAL*8(A-H,O-Z)
!!$      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
!!$      COMMON/ RASET2 / IRANMR,JRANMR
!!$      SAVE /RASET1/,/RASET2/
!!$!* This shows correspondence between the simplified input seeds IJ, KL
!!$!* and the original Marsaglia-Zaman seeds I,J,K,L.
!!$!* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
!!$!* k=56,l=78) put ij=1802, kl=9373
!!$      I = MOD( IJ/177 , 177 ) + 2
!!$      J = MOD( IJ     , 177 ) + 2
!!$      K = MOD( KL/169 , 178 ) + 1
!!$      L = MOD( KL     , 169 )
!!$      DO 300 II = 1 , 97
!!$        S =  0D0
!!$        T = .5D0
!!$        DO 200 JJ = 1 , 24
!!$          M = MOD( MOD(I*J,179)*K , 179 )
!!$          I = J
!!$          J = K
!!$          K = M
!!$          L = MOD( 53*L+1 , 169 )
!!$          IF(MOD(L*M,64) .GE. 32) S = S + T
!!$          T = .5D0*T
!!$  200   CONTINUE
!!$        RANU(II) = S
!!$  300 CONTINUE
!!$      RANC  =   362436D0 / 16777216D0
!!$      RANCD =  7654321D0 / 16777216D0
!!$      RANCM = 16777213D0 / 16777216D0
!!$      IRANMR = 97
!!$      JRANMR = 33
!!$      END




end module haag
