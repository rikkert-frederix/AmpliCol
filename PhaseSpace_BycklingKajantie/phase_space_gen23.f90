module phase_space_gen23
  use common
  private
  integer(kind=4) :: ix,ndim
  integer(kind=4),dimension(:),allocatable :: order
  real(kind=8),dimension(:),allocatable :: invm,invm_min,invm_max,x
  real(kind=8),dimension(:,:),allocatable :: pp
  integer(kind=4),dimension(:,:),allocatable :: sets
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical :: t_channel,includePDF
  real(kind=8) :: sqrtshat,sqrts,tau,ycm
  integer :: nquarks

  ! TECHNIAL PARAMETERS
  ! vebose:
  logical,parameter :: verbose=.true., debug=.false.
  ! importance sampling (0d0=flat transformation; -1d0=1/x transformation):
  real(kind=8),parameter :: ip=-1d0,ip_shat=-2d0
  ! tiny parameter cutoff to prevent/reduce numerical instabilities:
  real(kind=8),parameter :: vtiny=1d-12

  ! OUTPUT
  ! phase-space point
!  real(kind=8),dimension(:,:),allocatable,public :: p
  ! phase-space weight for the phase-space point (includes all factors
  ! of 2*pi and flux factor)
!  real(kind=8),public :: jac
  
  public :: gen23_init,gen23_phase_space
contains
  subroutine gen23_init(sqrtsh,n,m,o,part,s_cut,t_chan,include_pdf)
    ! Phase-space initialisation routines.
    implicit none
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrtsh
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! the colour order:
    integer(kind=4),dimension(n),intent(in) :: o
    ! cut on the minimum invariant mass of all pairs of particles,
    ! abs((p_i+p_j)^2)>s_cut, (initial and final state). s_cut(1) is between
    ! an initial and a final state particle; s_cut(2) is between two final
    ! state particles.
    real(kind=8),intent(in) :: s_cut(2)
    ! masses of all the particles. The two incoming particles must be
    ! massless.
    real(kind=8),dimension(n),intent(in) :: m
    ! ** With t_chan=.false., use the Byckling Kayantie 2->3 phase-space
    ! generation, i.e., eq (8) of E.~Byckling and K.~Kajantie,
    ! ``Reductions of the phase-space integral in terms of simpler
    ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    ! ** With t_chan=.true., use equation (13) instead.
    ! ** Note, for the first splitting, there is an additional logical
    ! parameter in the generate_momenta() subtroutine that allows one
    ! to choose between s-channel and t-channel.
    logical,intent(in) :: t_chan
    ! Should we include a PDF set? Currently, only the NNPDF2.3 NLO QED is available.
    logical,intent(in) :: include_pdf
    integer(kind=4) :: i,j
    integer(kind=4),dimension(n) :: part,process,temp_order,ord
    sqrtshat=sqrtsh
    sqrts=sqrtsh
    t_channel=t_chan
    if (verbose) then
       write (*,*) 'Setting up',n,'particle phase-space'
       write (*,*) 'Total available energy, sqrt(s-hat) =',sqrtshat
       write (*,*) 'Cut on invariants used in the phase-space generation: abs((p_i+p_j)^2) >=',s_cut
       write (*,*) 'Use the simple t-channel?',t_channel
    endif
    includePDF=include_pdf
    call gen23_deallocate
    next=n
    ndim=3*(next-2)-4
    if (includePDF) ndim=ndim+2 ! the two Bjorken x's
    allocate(order(next))
    allocate(invm(maskr(next)))
    allocate(invm_min(maskr(next)))
    allocate(invm_max(maskr(next)))
    allocate(pp(0:3,0:maskr(next)))
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
    if (verbose) write (*,*) 'masses:',m(1:n)
    call setup_PS_cuts(s_cut)
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
    write(*,*) 'TEMP ORDER',temp_order
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
    if (verbose) then
       write (*,*) "Power in importance sampling:",ip
    endif
    !stop 1
  end subroutine gen23_init

  subroutine define_gen_order(next,process,o,order)
    implicit none
    integer :: next,i,j,glu,quark,aquark
    integer(kind=4),dimension(next) :: process,o,order

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
            write(*,*) 'one'
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
            write(*,*) 'two',quark,aquark
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
            write(*,*) 'both'
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
            order(next-1)=quark
            order(next)=aquark
         endif
    endif

    order=o
    write(*,*) 'new order',order
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

  subroutine gen23_phase_space(xx)
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
       p(0:3,i)=pp(0:3,ibset(0,i-1))
    enddo
  end subroutine gen23_phase_space

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
    real(kind=8) :: smin,smax,shat
    smin=invm_min(maskr(next)-3)
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
  
  subroutine generate_momenta
    implicit none
    integer(kind=4) :: i,j,inext,im1
    integer(kind=4),dimension(2) :: set
    logical,parameter :: use_t_channel_at_start=.true.

    ! incoming momenta
    pp(0,ibset(0,0))=sqrtshat/2d0
    pp(0,ibset(0,1))=sqrtshat/2d0
    pp(1:2,ibset(0,0))=0d0
    pp(1:2,ibset(0,1))=0d0
    pp(3,ibset(0,0))= pp(0,ibset(0,0))
    pp(3,ibset(0,1))=-pp(0,ibset(0,1))
    ! incoming momenta when labeled from the other side
    pp(0:3,maskr(next)-ibset(0,0))=pp(0:3,1)
    pp(0:3,maskr(next)-ibset(0,1))=pp(0:3,2)
    ! invariant mass of all final state particles combined
    invm(ibset(0,0)+ibset(0,1))=sqrtshat**2
    invm(maskr(next)-ibset(0,0)-ibset(0,1))=sqrtshat**2

    set(1)=sets(0,1)
    set(2)=sets(0,2)

    ! Generate the central 2->2 process in case both set(1) and set(2) are not empty
    if (popcnt(set(1)).gt.1 .and. popcnt(set(2)).gt.1) then
       if (debug) write (*,*) 'two sets with at least two ',&
            & 'particles',popcnt(sets(0,1)),popcnt(sets(0,2))
       if (use_t_channel_at_start) then
          call gent_one_step(set(2),set(1),1)
       else
          call gens_one_step(set(2),set(1))
       endif
       if (jac.le.0d0) return
       pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
       invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
    elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).gt.1) then
       if (debug) write (*,*) 'special double t-channel (1)'&
            &,popcnt(sets(0,1)),popcnt(sets(0,2))
       call double_t(set(1),set(2),1,2)
       if (jac.le.0d0) return
       pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
       invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
    elseif (popcnt(set(1)).gt.1 .and. popcnt(set(2)).eq.1) then
       if (debug) write (*,*) 'special double t-channel (2)'&
            &,popcnt(sets(0,1)),popcnt(sets(0,2))
       call double_t(set(2),set(1),1,2)
       if (jac.le.0d0) return
       pp(0:3,set(1)+1)=pp(0:3,2)-pp(0:3,set(2))
       invm(set(1)+1)=dot(pp(0:3,set(1)+1),pp(0:3,set(1)+1))
    elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).eq.1) then
       if (debug) write (*,*) '2->2 scattering with one particle in each set'&
            &,popcnt(sets(0,1)),popcnt(sets(0,2))
       call gens_one_step(set(2),set(1))
       if (jac.le.0d0) return
       pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
       invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
    endif
    do i=1,2
       if (popcnt(set(i)).le.1) cycle ! at least 2 particles in a set
       inext=ibset(0,sets(1,i)-1)
       set(i)=set(i)-inext
       if (popcnt(set(i)).ge.2) then
          ! at least 3 particles in a set
          if (debug) write (*,*) 'At least 3 particles in a set',&
               & popcnt(sets(0,i)),popcnt(sets(0,3-i))
          call gent_one_step(set(i),inext,i)
          if (jac.le.0d0) return
          pp(0:3,(3-i)+set(i)+inext)=pp(0:3,i)-pp(0:3,sets(0,(3-i)))
          invm((3-i)+set(i)+inext)=dot(pp(0:3,(3-i)+set(i)+inext),pp(0:3,(3-i)+set(i)+inext))
          pp(0:3,(3-i)+set(i))=pp(0:3,i)-pp(0:3,sets(0,(3-i)))-pp(0:3,inext)
          invm((3-i)+set(i))=dot(pp(0:3,(3-i)+set(i)),pp(0:3,(3-i)+set(i)))
          do j=2,popcnt(set(i))
             ! loop over the remaining particles in the set
             inext=ibset(0,sets(j,i)-1)
             im1=ibset(0,sets(j-1,i)-1)
             set(i)=set(i)-inext
             if (t_channel) then
                call gent_one_step(inext,set(i),3-i)
             else
                call gen23_one_step(inext,set(i),3-i,im1)
             endif
             if (jac.le.0d0) return
          enddo
       elseif (popcnt(set(i)).eq.1 .and. popcnt(sets(0,3-i)).ne.0) then
          ! Exactly 2 particles in a set (and the other set contains at least one)
          if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
               & 'the other set contains at least one)', &
               & popcnt(sets(0,i)),popcnt(sets(0,3-i))
          im1=3-i
          pp(0:3,set(i)+inext+im1)=pp(0:3,set(i)+inext)-pp(0:3,im1)
          pp(0:3,set(i)+inext+i+im1)=pp(0:3,set(i)+inext+im1)-pp(0:3,i)
          invm(set(i)+inext+im1)=dot(pp(0:3,set(i)+inext+im1),pp(0:3,set(i)+inext+im1))
          invm(set(i)+inext+i+im1)=dot(pp(0:3,set(i)+inext+im1+i),pp(0:3,set(i)+inext+im1+i))
          if (t_channel) then
             call gent_one_step(set(i),inext,i)
          else
             call gen23_one_step(set(i),inext,i,im1)
          endif
          if (jac.le.0d0) return
       elseif (popcnt(set(i)).eq.1 .and. popcnt(sets(0,3-i)).eq.0) then
          ! Exactly 2 particles in a set (and the other set contains none)
          if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
               & 'the other set contains none)', &
               & popcnt(sets(0,i)),popcnt(sets(0,3-i))
          call gent_one_step(set(i),inext,i)
          if (jac.le.0d0) return
       else
          write (*,*) 'Inconsistent sets'
          write (*,*) i,':',sets(:,i)
          write (*,*) 3-i,':',sets(:,3-i)
          stop 
       endif
       ! We need to get the momentum of the final particle of the set.
       pp(0:3,set(i))=pp(0:3,set(i)+inext+(3-i))+pp(0:3,(3-i))-pp(0:3,inext)
    enddo
    if (debug) call test_momenta

! Add factors of 2*pi
    jac=jac/((2d0*pi)**(3*(next-2)-4))
! Add flux factor
    jac=jac/(2d0*sqrtshat**2)
    
  end subroutine generate_momenta

  subroutine test_momenta
    ! Writes the momenta to the screen.
    implicit none
    integer(kind=4) :: i,i1,i2,i1b,i2b
    real(kind=8),dimension(0:3) :: ptot
    write (*,*) 'Momenta check:'
    if (ix.ne.ndim) then
       write (*,*) 'ERROR: number of random numbers used not consistent',ix,ndim
       stop 1
    endif
    ptot(0:3)=0d0
    do i=0,next-1
       ptot(0:3)=ptot(0:3)+pp(0:3,ibset(0,i))
       write (*,*) i+1,pp(0:3,ibset(0,i)),dot(pp(0,ibset(0,i)),pp(0,ibset(0,i)))
    enddo
    write (*,*) 'ptot',ptot(0:3)
    do i=1,next
       i1=order(i)
       i2=order(mod(i,next)+1)
       i1b=ibset(0,i1-1)
       i2b=ibset(0,i2-1)
       write (*,*) '(pi+pj)^2',i1,i2,invm(i1b+i2b),invm(maskr(next)-(i1b+i2b))&
            &,dot(pp(0:3,i1b)+pp(0:3,i2b),pp(0:3,i1b)+pp(0:3,i2b))
    enddo
    write (*,*) 'number of dimensions',ndim
    write (*,*) ''
    write (*,*) ''
  end subroutine test_momenta


  subroutine double_t(i,ir,ia,ib)
! Generates the p_a+p_b->p_i+p_ir process, with (p_a+p_i)^2 and
! (p_b+p_i)^2 (and phi) as integration variables. Note that this means
! that the mass of the system ir is derived from the integration
! variables.
!
! p_a, p_b and p_i are assumed to be massless, and p_a and p_b
! back-to-back incoming particles.
    implicit none
    integer(kind=4),intent(in) :: i,ir,ia,ib
    real(kind=8) :: tmin,tmax,phi,pt2,esum,yr
    esum=sqrt(invm(ia+ib))
    yr=sqrt(lambda(invm(ia+ib),invm(i),invm_min(ir)))
    tmin=(-invm(ia+ib)+invm(i)+invm_min(ir)-yr)/2d0
    tmax=(-invm(ia+ib)+invm(i)+invm_min(ir)+yr)/2d0
    if (invm_max(i+ia).ne.0d0) tmax=min(invm_max(i+ia),tmax)
    if (invm_min(i+ia).ne.0d0) tmin=max(tmin,invm_min(i+ia))
    if (tmin.ge.tmax) then
       jac=-1d0
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(i+ia),jac)
    if (debug) then
       write (*,*) 'dt- i+ia',i+ia,invm(i+ia),invm_min(i+ia),invm_max(i+ia)
    endif
    tmin=-invm(ia+ib)-invm(i+ia)+invm(i)+invm_min(ir)
    tmax=invm(i)*(invm(i)-invm(ia+ib)-invm(i+ia))/(invm(i)-invm(i+ia))
    if (invm_max(i+ib).ne.0d0) tmax=min(invm_max(i+ib),tmax)
    if (invm_min(i+ib).ne.0d0) tmin=max(tmin,invm_min(i+ib))
    if (tmin.ge.tmax) then
       jac=-2d0
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(i+ib),jac)
    if (debug) then
       write (*,*) 'dt- i+ib',i+ib,invm(i+ib),invm_min(i+ib),invm_max(i+ib)
    endif
    ix=ix+1
    call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,jac)
    if (debug) then
       write (*,*) 'dt- phi',phi
    endif
    pt2=invm(i+ia)*invm(i+ib)/invm(ia+ib)+ &
         & invm(i)**2/invm(ia+ib)-(invm(i+ia)+invm(i+ib))*invm(i)/invm(ia+ib)-invm(i)
    if (pt2/invm(ia+ib).lt.-vtiny) then
       write (*,*) "ERROR in double_t: pt2 should be larger than zero", &
            & pt2,invm(i+ia),invm(i+ib),invm(ia+ib),invm(i)
       write (*,*) invm(i+ia)*invm(i+ib)/invm(ia+ib), &
            & invm(i)**2/invm(ia+ib),(invm(i+ia)+invm(i+ib))*invm(i)/invm(ia+ib),invm(i)
       stop 1
    elseif (pt2.lt.0d0) then
       pt2=0d0
    endif
    pp(0,i)=(-invm(i+ia)-invm(i+ib)+2d0*invm(i))/(2d0*esum)
    pp(1,i)=sqrt(pt2)*cos(phi)
    pp(2,i)=sqrt(pt2)*sin(phi)
    pp(3,i)=(invm(i+ia)-invm(i+ib))/(2d0*esum)
    pp(0,ir)=esum-pp(0,i)
    pp(1:3,ir)=-pp(1:3,i)
    invm(ir)=dot(pp(0,ir),pp(0,ir))
    if (invm(ir).le.0d0) then
       write (*,*) "ERROR in double_t: invariant mass of system", &
            & " must be larger than zero",ir,invm(ir),i
       write (*,*) invm(ir),invm(ia+ib)+invm(i+ia)+invm(i+ib)-invm(i)&
            &,invm(ia+ib),invm(i +ia),invm(i+ib),invm(i)
       stop
    endif
    jac = jac/(4d0*sqrt(lambda(invm(ir+i),0d0,0d0)))
  end subroutine double_t
  
  subroutine gen23_one_step(i,ir,ib,im1)
! Generates one step using the 2->3 setup, using the invariants
! shat(i), s(i) and t(i) as defined in E.~Byckling and K.~Kajantie,
! ``Reductions of the phase-space integral in terms of simpler
! processes,'' Phys. Rev. 187 (1969), 2008-2016,
! doi:10.1103/PhysRev.187.2008.  Assumes massless incoming particles.
    implicit none
    integer(kind=4),intent(in) :: im1,i,ir,ib
    real(kind=8) :: tmin,tmax,smin,smax,phi,gram4,V,sqrtGG
    call generate_masses(i,ir)
    if (debug) then
       write (*,*) '23- i    ',i,invm(i)
       write (*,*) '23- ir   ',ir,invm(ir)
    endif
    if (jac.le.0d0) return
    call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin,tmax)
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    if (tmin.ge.tmax) then
       jac=-3d0
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ir+ib),jac)
    if (debug) then
       write (*,*) '23- ir+ib',ir+ib,invm(ir+ib),invm_min(ir+ib),invm_max(ir+ib)
    endif
    call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
         &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin,smax,V,sqrtGG)
    if (invm_min(i+im1).ne.0d0) smin=max(smin,invm_min(i+im1))
    if (invm_max(i+im1).ne.0d0) smax=min(smax,invm_max(i+im1))
    if (smin.ge.smax) then
       jac=-4d0
       if (debug) write (*,*) 'smin.ge.smax',smin,smax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,smin,smax,invm(i+im1),jac)
    if (debug) then
       write (*,*) '23- i+im1',i+im1,invm(i+im1),invm_min(i+im1),invm_max(i+im1)
    endif
    phi=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi &
            &,sqrt(invm(i)),sqrt(invm(ir)),pp(0,i),pp(0,ib+ir))
    if (im1.le.2) then
       pp(0:3,ir)=pp(0:3,ir+i)-pp(0:3,i)
    endif
    gram4=gram_determinant4(invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib)&
         &,invm(ir+i),invm(i+im1),invm(ir+ib+i+im1),invm(ir),invm(i)&
         &,invm(im1))
    if (gram4.ge.0d0) then 
       write (*,*) 'error, gram4 greater than or equal to zero',gram4,i,ir
       jac=-5d0
       return
    endif
    jac=jac/(8d0*sqrt(-gram4))
  end subroutine gen23_one_step


  subroutine gent_one_step(i,ir,ib)
    ! One step in the usual MadGraph t-channel phase-space generation.
    implicit none
    integer(kind=4),intent(in) :: i,ir,ib
    real(kind=8) :: tmin,tmax,phi
    call generate_masses(i,ir)
    if (debug) then
       write (*,*) 't - i    ',i,invm(i)
       write (*,*) 't - ir   ',ir,invm(ir)
    endif
    if (jac.le.0d0) return
    call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin,tmax)
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    if (tmin.ge.tmax) then
       jac=-6d0
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ir+ib),jac)
    if (debug) then
       write (*,*) 't - ir+ib',ir+ib,invm(ir+ib),invm_min(ir+ib),invm_max(ir+ib)
    endif
    ix=ix+1
    call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,jac)
    call gentcms(pp(0,ib+ir+i),pp(0,ib),invm(ib+ir),phi,sqrt(invm(i)) &
         &,sqrt(invm(ir)),pp(0,i),pp(0,ib+ir))
    pp(0:3,ir)=pp(0:3,ib+ir+i)+pp(0:3,ib)-pp(0:3,i)
    jac = jac/(4d0*sqrt(lambda(invm(ir+i),0d0,invm(ir+i+ib))))
  end subroutine gent_one_step

  subroutine gens_one_step(i,ir)
    implicit none
    integer(kind=4),intent(in) :: i,ir
    real(kind=8) :: esum,costh,phi
    real(kind=8),dimension(0:3) :: p_i,p_ir
    call generate_masses(i,ir)
    esum=sqrt(invm(i+ir))
    ix=ix+1
    call random_to_var(x(ix),0d0,-1d0,1d0,costh,jac)
    ix=ix+1
    call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,jac)
    call mom2cx(esum,sqrt(invm(i)),sqrt(invm(ir)),costh,phi,p_i,p_ir)
    call boostm(p_i,p(0:3,i+ir),esum,pp(0:3,i))
    call boostm(p_ir,p(0:3,i+ir),esum,pp(0:3,ir))
    jac=jac*sqrt(lambda(invm(i+ir),invm(i),invm(ir)))/(8d0*invm(i+ir))
    ! fill t-channel stuff to be safe.
    pp(0:3,i+1)=pp(0:3,i)-pp(0:3,1)
    pp(0:3,i+2)=pp(0:3,i)-pp(0:3,2)
    pp(0:3,ir+1)=pp(0:3,ir)-pp(0:3,1)
    pp(0:3,ir+2)=pp(0:3,ir)-pp(0:3,2)
    invm(i+1)=dot(pp(0:3,i+1),pp(0:3,i+1))
    invm(i+2)=dot(pp(0:3,i+2),pp(0:3,i+2))
    invm(ir+1)=dot(pp(0:3,ir+1),pp(0:3,ir+1))
    invm(ir+2)=dot(pp(0:3,ir+2),pp(0:3,ir+2))
  end subroutine gens_one_step

  subroutine generate_masses(i,ir)
    ! Generates the two invariant masses for (combined) particles i and ir. 
    implicit none
    integer(kind=4),intent(in) :: i,ir
    integer(kind=4) :: j,j1,j2
    real(kind=8) :: shatmin,shatmax
    if (popcnt(i).gt.1) invm(i)=0d0
    if (popcnt(ir).gt.1) invm(ir)=0d0
    do j=1,2
       if (j.eq.1) then
          j1=ir
          j2=i
       else
          j1=i
          j2=ir
       endif
       if (popcnt(j1).ge.2) then
          call shatminmax(j1,j2,shatmin,shatmax)
          if (invm_min(j1).ne.0d0) shatmin=max(shatmin,invm_min(j1))
          if (invm_max(j1).ne.0d0) shatmax=max(shatmax,invm_max(j1))
          if (shatmin.ge.shatmax) then
             jac=-7d0
             if (debug) write (*,*) 'shatmin.ge.shatmax',j,i,ir,shatmin,shatmax,invm(j2)
             return
          endif
          ix=ix+1
          call random_to_var(x(ix),ip,shatmin,shatmax,invm(j1),jac)
       endif
    enddo
  end subroutine generate_masses
  
  subroutine shatminmax(j1,j2,shatmin,shatmax)
    ! Determines minimum and maximum allowed s-channel invariant
    ! masses based on previously generated masses and the masses of
    ! the final-state particles.
    implicit none
    integer(kind=4),intent(in) :: j1,j2
    real(kind=8),intent(out) :: shatmin,shatmax
    integer(kind=4) :: j
    shatmin=0d0
    do j=0,next-1
       if (btest(j1,j)) then
          shatmin=shatmin+sqrt(invm(ibset(0,j)))
       endif
    enddo
    shatmin=shatmin**2
    shatmax=(sqrt(invm(j1+j2))-sqrt(max(invm(j2),invm_min(j2))))**2
  end subroutine shatminmax

    


  






  real(kind=8) function lambda(s,xa2,xb2)
! The usual two dimensional phase-space volume factor. See, e.g.,
! Eq.(A2) of E.~Byckling and K.~Kajantie, ``Reductions of the
! phase-space integral in terms of simpler processes,'' Phys. Rev. 187
! (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8) :: xa2,xb2,S,rat
    real(kind=8),parameter :: tiny=1.d-8
    lambda=s**2-2d0*(xa2+xb2)*s+(xa2-xb2)**2
    if(lambda.le.0.d0)then
       if(xa2.lt.0.d0.or.xb2.lt.0.d0)then
          write(6,*)'Error #1 in function Lambda:',s,xa2,xb2
          stop
       endif
       rat=1-(sqrt(xa2)+sqrt(xb2))/s
       if(rat.gt.-tiny)then
          lambda=0.d0
       else
          write(6,*)'Error #2 in function Lambda:',s,xa2,xb2,rat
       endif
    endif
  end function lambda
  
  subroutine boostm(p,q,m,pboost)
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

  subroutine mom2cx(esum,mass1,mass2,costh1,phi1,p1,p2)
! This subroutine sets up two four-momenta in the two particle rest
! frame.
! input:
!       real    esum           : energy sum of particle 1 and 2
!       real    mass1          : mass            of particle 1
!       real    mass2          : mass            of particle 2
!       real    costh1         : cos(theta)      of particle 1
!       real    phi1           : azimuthal angle of particle 1
! output:
!       real    p1(0:3)        : four-momentum of particle 1
!       real    p2(0:3)        : four-momentum of particle 2
    implicit none
    real(kind=8),dimension(0:3) :: p1,p2
    real(kind=8) :: esum,mass1,mass2,costh1,phi1,md2,ed,pp,sinth1
    real(kind=8),parameter :: rZero = 0d0,rHalf = 0.5d0,rOne = 1d0,rTwo = 2d0
    md2 = (mass1-mass2)*(mass1+mass2)
    ed = md2/esum
    if ( mass1*mass2.eq.rZero ) then
       pp = (esum-abs(ed))*rHalf
    else
       pp = sqrt((md2/esum)**2-rTwo*(mass1**2+mass2**2)+esum**2)*rHalf
    endif
    sinth1 = sqrt((rOne-costh1)*(rOne+costh1))
    p1(0) = max((esum+ed)*rHalf,rZero)
    p1(1) = pp*sinth1*cos(phi1)
    p1(2) = pp*sinth1*sin(phi1)
    p1(3) = pp*costh1
    p2(0) = max((esum-ed)*rHalf,rZero)
    p2(1:3) = -p1(1:3)
  end subroutine mom2cx

  subroutine gentcms(pa,pb,t,phi,m1,m2,p1,pr)
! Generates 4 momentum for particle p1, and remainder pr given the
! values t, and phi in the process pb+pa -> pr+p1.  Assuming incoming
! particles with momenta pa, pb and outgoing particles with mass
! m1,m2; t=(pa-p1)^2 ; phi is the azimuthal angle between p1 and pa in
! the pa+pb rest frame, with pa aligned with the positive z-axis.
    implicit none
    real(kind=8),intent(in) :: t,phi,m1,m2
    real(kind=8),intent(in),dimension(0:3) :: pa,pb
    real(kind=8),intent(out),dimension(0:3) :: p1,pr
    real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,ma2,pt,pt2
    real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,p1_rot
    real(kind=8),parameter :: tiny=1d-8
    ptot(0:3)=pa(0:3)+pb(0:3)
    ptotm(0)=ptot(0)
    ptotm(1:3)=-ptot(1:3)
    ma2=dot(pa,pa)
    ! determine magnitude of p1 in cms frame (from dhelas routine mom2cx)
    ESUM2 = dot(ptot,ptot)
    if (esum2 .le. 0d0) then
       write (*,*) "error :: must be time-like momentum in gentcms",esum2
       stop 1
    endif
    esum=sqrt(esum2)
    MD2=(M1-M2)*(M1+M2)
    ED=MD2/ESUM
    IF (M1*M2.EQ.0.d0) THEN
       PP2=0.25d0*(ESUM-ABS(ED))**2
    ELSE
       PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
       if(pp2.lt.0d0) then
          write(*,*) 'Error #12 in genps_fks.f: magnitude^2 of '/&
               &/'3-momentum smaller than 0',pp2
          stop 1
       endif
    ENDIF

    call boostm(pa,ptotm,esum,pa_cms)
    E_acms = pa_cms(0)
    p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)
    
    ! define p1 in the frame where pa_cms is aligned with the positive z axis.
    p1(0) = MAX((ESUM+ED)*0.5d0,0.d0)
    if (esum+ed.le.0d0) then
       write (*,*) 'Error #14 in genps_fks.f: negative energy',esum,ed
       jac=-8d0
       return
       write (*,*) pa(0:3)
       write (*,*) pb(0:3)
       write (*,*) m1,m2,t,phi
       write (*,*) pa_cms(0:3)
       stop 1
    endif
    p1(3) = -(m1**2+ma2-t-2d0*p1(0)*E_acms)/(2d0*p_acms)
    pt2=pp2-p1(3)**2
    if (pt2/esum2.lt.-tiny) then
       write (*,*) 'Error #13 in genps_fks.f: relative pt^2 smaller than 0',pt2,esum2
       stop 1
    elseif (pt2.lt.0d0) then
       pt2=0d0
    endif
    pt = sqrt(pt2)
    p1(1) = pt*cos(phi)
    p1(2) = pt*sin(phi)
    call rotxxx(p1,pa_cms,p1_rot)       !Rotate p1 to the pa_cms frame
    call boostm(p1_rot,ptot,esum,p1)    !boost back to lab fram
    pr(0:3)=pa(0:3)-p1(0:3)         !Return remainder of momentum
  end subroutine gentcms


  subroutine gentcms2(pa,pb,pc,t,phi,m1,m2,p1,pr)
! Generates 4 momentum for particle p1, and remainder pr given the
! values t, and phi in the process pa+pb -> pr+p1.  Assumes incoming
! particles with momenta pa, pb and outgoing particles with mass
! m1,m2; t=(pb-p1)^2. Assumes that pa is a massless momentum; phi is
! the azimuthal angle between pr and pc in the pa+pb rest frame, with
! pa aligned with the positive z-axis.
    implicit none
    real(kind=8),intent(in) :: t,phi,m1,m2
    real(kind=8),intent(in),dimension(0:3) :: pa,pb,pc
    real(kind=8),intent(out),dimension(0:3) :: p1,pr
    real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,pt,pt2,phi_off
    real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,Pii,pc_cms,pc_rot,pii_rot
    real(kind=8),parameter :: tiny=1d-8
    ptot(0:3)=pa(0:3)+pb(0:3)
    ptotm(0)=ptot(0)
    ptotm(1:3)=-ptot(1:3)
    ! determine magnitude of Pii in cms frame (from dhelas routine mom2cx)
    ESUM2 = dot(ptot,ptot)
    if (esum2 .le. 0d0) then
       jac=-9d0
       write (*,*) "error :: must be time-like momentum in gentcms2",esum2
       return
       stop 1
    endif
    esum=sqrt(esum2)
    MD2=(M2-M1)*(M1+M2)
    ED=MD2/ESUM
    IF (M1*M2.EQ.0.d0) THEN
       PP2=0.25d0*(ESUM-ABS(ED))**2
    ELSE
       PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
       if(pp2.lt.0d0) then
          write(*,*) 'Error #12 in genps_fks.f: magnitude^2 of '/&
               &/'3-momentum smaller than 0',pp2
          stop 1
       endif
    ENDIF
    call boostm(pa,ptotm,esum,pa_cms)
    E_acms = pa_cms(0)
    p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)

    ! determine the offset in phi; the frame in which phi is defined
    ! is in the ptot rest-frame, with pa_cms aligned with the z-axis,
    ! and pc having zero phi angle.
    call boostm(pc,ptotm,esum,pc_cms)
    call rotxxx_inv(pc_cms,pa_cms,pc_rot)
    if (pc_rot(1).ne.0d0) then
       phi_off=atan(pc_rot(2)/pc_rot(1))
    else
       phi_off=0d0
    endif
    if (pc_rot(1).lt.0d0) then
       phi_off=phi_off+pi
    endif
    
    ! define Pii in the frame where pa_cms is aligned with the positive z axis
    Pii(0) = MAX((ESUM+ED)*0.5d0,0.d0)
    if (esum+ed.le.0d0) then
       write (*,*) 'Error #15 in genps_fks.f: negative energy',esum,ed
       jac=-10d0
       return
       stop 1
    endif
    Pii(3) = -(m2**2-t-2d0*Pii(0)*E_acms)/(2d0*p_acms)
    pt2=pp2-Pii(3)**2
    if (pt2/esum2.lt.-tiny) then
       write (*,*) 'Error #16 in genps_fks.f: relative pt^2 smaller than 0',pt2
       jac=-11d0
       return
       stop 1
    elseif (pt2.lt.0d0) then
       pt2=0d0
    endif
    pt = sqrt(pt2)
    Pii(1) = pt*cos(phi+phi_off)
    Pii(2) = pt*sin(phi+phi_off)
    call rotxxx(Pii,pa_cms,Pii_rot)       !Rotate Pii to the pa_cms frame
    call boostm(Pii_rot,ptot,esum,Pii)    !boost back to lab fram
    p1(0:3)=ptot(0:3)-pii(0:3)
    pr(0:3)=pb(0:3)-p1(0:3)         !Return remainder of momentum
  end subroutine gentcms2






!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c Functions not to be touched (basic functions)
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine rotxxx(p,q,prot)
! This subroutine performs the spacial rotation of a four-momentum.
! the momentum p is assumed to be given in the frame where the spacial
! component of q points the positive z-axis.  prot is the momentum p
! rotated to the frame where q is given.
! input:
!       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame
!!       real    q(0:3)         : four-momentum q in the rotated frame
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


 subroutine rotxxx_inv(p,q,prot)
    ! Same as rotxxx, but inverse. That is, first doing
    ! rotxxx(p,q,prot) and then rotxxx_inv(prot,q,p) should give you
    ! back the original p.
    implicit none
    real(kind=8),dimension(0:3),intent(in) :: p,q
    real(kind=8),dimension(0:3),intent(out) :: prot
    real(kind=8) :: qt2,qt,psgn,qq
    prot(0) = p(0)
    qt2 = q(1)**2 + q(2)**2
    if ( qt2.eq.0d0 ) then
       if ( q(3).eq.0d0 ) then
          prot(1:3)=p(1:3)
       else
          psgn = sign(1d0,q(3))
          prot(1:3)=p(1:3)*psgn
       endif
    else
       qq = sqrt(qt2+q(3)**2)
       qt = sqrt(qt2)
       prot(1) = q(1)*q(3)/qq/qt*p(1) +q(2)*q(3)/qq/qt*p(2) -  qt/qq*p(3)
       prot(2) =        -q(2)/qt*p(1) +        q(1)/qt*p(2)
       prot(3) =   qt*q(1)/qq/qt*p(1) +        q(2)/qq*p(2) +q(3)/qq*p(3)
    endif
  end subroutine rotxxx_inv

  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot

  subroutine random_to_var(x,power,var_min,var_max,var,jac)
    ! Given a random number x, it generates var in the range var_min
    ! <= var <= var_max according to var^(power)
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
       stop 1
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
  
  real(kind=8) function getphifroms(si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG)
    ! Given s_i (invariant mass of p_i and p_i+1, it transforms it
    ! into phi_i. Note that there are two possibilities for phi: need
    ! to pick one at random.
    ! Based on eq.(11) of E.~Byckling and K.~Kajantie, ``Reductions of
    ! the phase-space integral in terms of simpler processes,''
    ! Phys. Rev. 187 (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG
    real(kind=8) :: cosphi,x
    real(kind=8),external :: ran2
    cosphi=((si-shat_im1-shat_ip1)*0.5d0*lambda(shat_i,t_i,0d0)-4d0*V)/sqrtGG
    x=ran2()
    if (x.gt.0.5d0) then
       getphifroms=acos(cosphi)
    else
       getphifroms=-acos(cosphi)+2d0*pi
    endif
  end function getphifroms

  real(kind=8) function computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    ! Computes the determinant of V (with m_a=0) based on eq.(11) of
    ! E.~Byckling and K.~Kajantie, ``Reductions of the phase-space
    ! integral in terms of simpler processes,'' Phys. Rev. 187 (1969),
    ! 2008-2016, doi:10.1103/PhysRev.187.2008
    use LUPdecomposition
    real(kind=8),intent(in) :: shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2
    real(kind=8) :: deter
    integer(kind=4),parameter :: n=3
    real(kind=8),dimension(n,n) :: a
    integer(kind=4),dimension(0:n) :: p
    real(kind=8),parameter :: tol=1d-10
    logical :: success
    a(1:3,1)=(/2d0*shat_i             , shat_i-t_i    , shat_i+shat_im1-m_i_2/)
    a(1:3,2)=(/shat_i-t_i             , 0d0           , shat_im1-t_im1       /)
    a(1:3,3)=(/shat_ip1+shat_i-m_ip1_2, shat_ip1-t_ip1, 0d0                  /)
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       computeV=-deter/8d0
    else
       computeV=-1d0
    endif
  end function computeV
  
  real(kind=8) function G(x,y,z,u,v,w)
    ! The G-function, eq.(A5) of E.~Byckling and K.~Kajantie,
    ! ``Reductions of the phase-space integral in terms of simpler
    ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: x,y,z,u,v,w
    G=x**2*y+x*y**2+z**2*u+z*u**2+v**2*w+v*w**2&
         & +x*z*w+x*u*v+y*z*v+y*u*w &
         & -x*y*(z+u+v+w)-z*u*(x+y+v+w)-v*w*(x+y+z+u)
  end function G
  
  subroutine sminmax(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2&
       &,m_ip1_2,smin,smax,V,sqrtGG)
    ! Determines the integration limits for s_i, as determined by
    ! setting cos(phi) to +/- 1 in eq.(11) of E.~Byckling and
    ! K.~Kajantie, ``Reductions of the phase-space integral in terms
    ! of simpler processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    implicit none
    real(kind=8),intent(in) :: shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2
    real(kind=8),intent(out) :: smin,smax,V,sqrtGG
    real(kind=8) :: GG,s1,s2
    V=computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    GG = G(t_i  , shat_ip1, shat_i  , t_ip1, m_ip1_2, 0d0) &
      & *G(t_im1, shat_i  , shat_im1, t_i  , m_i_2  , 0d0)
    if (GG.le.0d0) then
       write (*,*) 'No allowed range for s: smin=smax',GG
!!$       stop 1
       GG=0d0
    endif
    sqrtGG=sqrt(GG)
    s1=shat_im1+shat_ip1+2d0/lambda(shat_i,t_i,0d0) * (4d0*V + sqrtGG)
    s2=shat_im1+shat_ip1+2d0/lambda(shat_i,t_i,0d0) * (4d0*V - sqrtGG)
    smin=min(s1,s2)
    smax=max(s1,s2)
  end subroutine sminmax

  

  subroutine tminmax(X,Z,U,V,W,tmin,tmax)
    ! Determines the integration limits for t_i, as determined by
    ! setting cos(theta) to +/- 1 in eq.(10) of E.~Byckling and
    ! K.~Kajantie, ``Reductions of the phase-space integral in terms
    ! of simpler processes,'' Phys. Rev. 187 (1969), 2008-2016,
    ! doi:10.1103/PhysRev.187.2008.
    ! x = shat_i
    ! z = t_i
    ! u = shat_{i-1}
    ! v = m_i^2
    ! w = m_a^2
    implicit none
    real(kind=8),intent(in) :: x,z,u,v,w 
    real(kind=8),intent(out):: tmin,tmax
    real(kind=8) ::  t1,t2,yr
    yr = lambda(x,u,v)*lambda(x,w,z)
    if (yr.le.0d0) then
       write (*,*) 'No allowed range for t: tmin=tmax',yr
!!$       stop 1
       yr=0d0
    endif
    yr=sqrt(yr)
    t1 = u+w - ((x+u-v)*(x+w-z) - yr)/(2d0*x)
    t2 = u+w - ((x+u-v)*(x+w-z) + yr)/(2d0*x)
    tmin = min(t1,t2)
    tmax = max(t1,t2)
  end subroutine tminmax
    
  real(kind=8) function gram_determinant4(shat_ip1,t_im1,t_i,shat_i,s_i,t_ip1&
       &,shat_im1,m_i_2,m_ip1_2)
    ! Computes the 4x4 Gram determinant (with m_a=0) as in eq.(B6), of
    ! E.~Byckling and K.~Kajantie, ``Reductions of the phase-space
    ! integral in terms of simpler processes,'' Phys. Rev. 187 (1969),
    ! 2008-2016, doi:10.1103/PhysRev.187.2008.
    use LUPdecomposition
    implicit none
    real(kind=8) :: shat_ip1,t_im1,t_i,shat_i,s_i,t_ip1,shat_im1,m_i_2,m_ip1_2
    integer(kind=4),parameter :: n=4
    real(kind=8),dimension(n,n) :: a
    integer(kind=4),dimension(0:n) :: p
    real(kind=8),parameter :: tol=1d-10
    real(kind=8) :: deter
    logical :: success
    a(1:4,1)=(/ 0d0           , t_im1-shat_im1 , t_i-shat_i       , t_ip1-shat_ip1    /)
    a(1:4,2)=(/ t_im1-shat_im1, 2d0*t_im1      , t_i+t_im1-m_i_2  , t_im1+t_ip1-s_i   /)
    a(1:4,3)=(/ t_i-shat_i    , t_i+t_im1-m_i_2, 2d0*t_i          , t_i+t_ip1-m_ip1_2 /)
    a(1:4,4)=(/ t_ip1-shat_ip1, t_im1+t_ip1-s_i, t_i+t_ip1-m_ip1_2, 2d0*t_ip1         /)
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       gram_determinant4=deter/16d0
    else
       gram_determinant4=1d0
    endif
  end function gram_determinant4
  
end module phase_space_gen23
