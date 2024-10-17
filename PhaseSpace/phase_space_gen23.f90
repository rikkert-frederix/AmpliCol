module phase_space_gen23
  use common
  private
  integer(kind=4) :: ix,ndim,next
  integer(kind=4),dimension(:),allocatable :: order
  real(kind=8),dimension(:),allocatable :: invm,invm_min,invm_max,x,ETmin
  real(kind=8),dimension(:,:),allocatable :: pp
  integer(kind=4),dimension(:,:),allocatable :: sets
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical :: t_channel,includePDF
  real(kind=8) :: sqrtshat,sqrts,ycm,ptcut,drcut

  ! TECHNIAL PARAMETERS
  ! vebose:
  logical,parameter :: verbose=.false.
  logical,parameter,public :: debug=.false.
  ! importance sampling (0d0=flat transformation; -1d0=1/x transformation):
  real(kind=8),parameter :: ip=-1d0,ip_shat=-1.2d0
  ! tiny parameter cutoff to prevent/reduce numerical instabilities:
  real(kind=8),parameter :: vtiny=1d-12,tiny=1d-8

  ! OUTPUT
  ! phase-space point
!  real(kind=8),dimension(:,:),allocatable,public :: p
  ! phase-space weight for the phase-space point (includes all factors
  ! of 2*pi and flux factor)
!  real(kind=8),public :: jac
  
  public :: gen23_init,gen23_phase_space
contains
  subroutine gen23_init(sqrtsh,n,m,o,s_cut,pt_cut,dr_cut,t_chan,include_pdf)
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
    real(kind=8),intent(in) :: s_cut(2),pt_cut,dr_cut
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
    allocate(ETmin(maskr(next)))
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
       invm(ibclr(maskr(next),i-1))=m(i)**2
    enddo
    if (verbose) write (*,*) 'masses:',m(1:n)
    if (pt_cut.gt.0d0) then
       drcut=dr_cut
       ptcut=pt_cut
    else
       drcut=0d0
       ptcut=0d0
    endif
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
    sets=0
    i=0
    do i=2,next
       if (order(i).eq.2) then
          do j=i+1,next
             sets(0,2)=ibset(sets(0,2),order(j)-1)
          enddo
          sets(1:i-2,1)=order(2:i-1)
          sets(1:next-i,2)=order(i+1:next)
          exit
       endif
       sets(0,1)=ibset(sets(0,1),order(i)-1)
    enddo
    if (verbose) then
       write (*,*) "set 1:",sets(:,1)
       write (*,*) "set 2:",sets(:,2)
    endif
    if (verbose) then
       write (*,*) "Power in importance sampling:",ip
    endif
  end subroutine gen23_init

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
          if (npart.eq.next-2) then
             invm_min(i)=max(s_cut(2)*npart**2,mass**2)
          else
             invm_min(i)=max(s_cut(2)*(npart)*(npart-1)/2d0,mass**2)
          endif
       endif
    enddo
    call setup_ETmin
  end subroutine setup_PS_cuts

  subroutine setup_ETmin
    ! Setup the minimum required (transverse) energy for each (combination of)
    ! final state particle(s) in the collision c.o.m. frame. Based on the
    ! masses and the ptcut (i.e., assumes that all pz=0 and pT=pTcut)
    integer :: i,j
    ETmin(1:maskr(next))=0d0
    do i=1,maskr(next)
       if (btest(i,0).or.btest(i,1)) cycle ! skip the ones that include incoming particles
       do j=0,next-1
          if (btest(i,j)) ETmin(i)=ETmin(i)+sqrt(invm(ibset(0,j))+ptcut**2)
       enddo
       ETmin(i)=max(ETmin(i),sqrt(invm_min(i)))
    enddo
  end subroutine setup_ETmin
  
  subroutine gen23_deallocate
    implicit none
    if (allocated(order)) deallocate(order)
    if (allocated(invm)) deallocate(invm)
    if (allocated(invm_min)) deallocate(invm_min)
    if (allocated(invm_max)) deallocate(invm_max)
    if (allocated(ETmin)) deallocate(ETmin)
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
       if (includePDF) then
          ! Note: 'ycm' is the rapidity needed to go from lab to CM
          ! frame. Hence, here we boost from CM to lab frame with '-ycm'
          call boostz(pp(0:3,ibset(0,i-1)),-ycm,p(0:3,i))
       else
          p(0:3,i)=pp(0:3,ibset(0,i-1))
       endif
    enddo
  end subroutine gen23_phase_space

  subroutine generate_initial_state
    implicit none
    real(kind=8) :: tau
    call generate_tau(tau)
    call generate_y(tau)
    sqrtshat=sqrt(tau)*sqrts
    xbjrk(1)=sqrt(tau)*exp(ycm)
    xbjrk(2)=sqrt(tau)*exp(-ycm)
    if (debug) write (*,*) 'sqrtshat :',sqrtshat,xbjrk(1:2),sqrtshat**2
  end subroutine generate_initial_state

  subroutine generate_tau(tau)
    implicit none
    real(kind=8),intent(out) :: tau
    real(kind=8) :: smin,smax,shat
    smin=max(invm_min(maskr(next)-3),ETmin(maskr(next)-3)**2)
    smax=sqrts**2
    ix=ix+1
    call random_to_var(x(ix),ip_shat,smin,smax,shat,jac)
    tau=shat/smax
    jac=jac/smax
  end subroutine generate_tau
  
  subroutine generate_y(tau)
    implicit none
    real(kind=8),intent(in) :: tau
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
!!$       call gens_one_step(set(2),set(1))
       call gent_one_step(set(2),set(1),1)
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
          do j=2,popcnt(set(i))-1
             ! loop over the remaining particles in the set
             inext=ibset(0,sets(j,i)-1)
             im1=ibset(0,sets(j-1,i)-1)
             set(i)=set(i)-inext
             if (t_channel) then
                call gent_one_step(inext,set(i),3-i)
             else
!!$                call gen23_one_step_v2(inext,set(i),3-i,im1)
                call gen23_one_step(inext,set(i),3-i,im1)
!!$                call genpt_one_step(inext,set(i),3-i,im1)
!!$                call gent_one_step_v2(inext,set(i),3-i,im1)
             endif
             if (jac.le.0d0) return
          enddo
          inext=ibset(0,sets(j,i)-1)
          im1=ibset(0,sets(j-1,i)-1)
          set(i)=set(i)-inext
          if (t_channel) then
             call gent_one_step(inext,set(i),3-i)
          else
!!$             call gen23_one_step_v2(inext,set(i),3-i,im1)
             call gen23_one_step(inext,set(i),3-i,im1)
!!$             call genpt_one_step(inext,set(i),3-i,im1)
          endif
          if (jac.le.0d0) return
!!$          inext=ibset(0,sets(j+1,i)-1)
!!$          im1=ibset(0,sets(j,i)-1)
!!$          set(i)=set(i)-inext
!!$          if (t_channel) then
!!$             call gent_one_step(inext,set(i),3-i)
!!$          else
!!$!!$             call gen23_one_step_v2(inext,set(i),3-i,im1)
!!$             call gen23_one_step(inext,set(i),3-i,im1)
!!$          endif
!!$          if (jac.le.0d0) return
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
! p_a and p_b are assumed to be massless, and p_a and p_b
! back-to-back incoming particles. The mass of i is fixed.
    implicit none
    integer(kind=4),intent(in) :: i,ir,ia,ib
    real(kind=8) :: tmin,tmax,phi,pt2,yr,Eimax,pzmax
    if (popcnt(i).ne.1 .or. popcnt(ir).le.1) then
       write (*,*) 'Subroutine only for i is a single particle '&
            //'and ir is more than 1',i,ir,popcnt(i),popcnt(ir)
       stop 1
    endif
    yr=sqrt(lambda(invm(ia+ib),invm(i),invm_min(ir)))
    tmin=(-invm(ia+ib)+invm(i)+invm_min(ir)-yr)/2d0
    tmax=(-invm(ia+ib)+invm(i)+invm_min(ir)+yr)/2d0
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))

    ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
    ! being larger than ETmin(i) and ETmin(ir), respectively:
    pzmax=sqrt(lambda(sqrtshat**2,ETmin(i)**2,ETmin(ir)**2))/(2d0*sqrtshat)
    Eimax=sqrtshat-sqrt(ETmin(ir)**2+pzmax**2)
    
    tmin=max(tmin,invm(i)-sqrtshat*(Eimax+pzmax))
    tmax=min(tmax,invm(i)-sqrtshat*(Eimax-pzmax))
    
    if (tmin.ge.tmax) then
       jac=-1d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(i+ia),jac)

    if (debug) then
       write (*,*) 'dt- i+ia',i+ia,invm(i+ia),tmin,tmax
    endif
    
    tmin=-invm(ia+ib)-invm(i+ia)+invm(i)+invm_min(ir)
    tmax=invm(i)*(invm(i)-invm(ia+ib)-invm(i+ia))/(invm(i)-invm(i+ia))
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    
    ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
    ! being larger than ETmin(i) and ETmin(ir), respectively:
    tmin=max(tmin,invm(i)-sqrtshat**2*(1-ETmin(ir)**2/(sqrtshat**2+invm(i+ia)-invm(i))))
    tmax=min(tmax,invm(i)-sqrtshat**2*(ETmin(i)**2/(invm(i)-invm(i+ia))))
    
    if (tmin.ge.tmax) then
       jac=-2d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(i+ib),jac)

    if (debug) then
       write (*,*) 'dt- i+ib',i+ib,invm(i+ib),tmin,tmax
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
    pp(0,i)=(-invm(i+ia)-invm(i+ib)+2d0*invm(i))/(2d0*sqrtshat)
    pp(1,i)=sqrt(pt2)*cos(phi)
    pp(2,i)=sqrt(pt2)*sin(phi)
    pp(3,i)=(invm(i+ia)-invm(i+ib))/(2d0*sqrtshat)
    pp(0,ir)=sqrtshat-pp(0,i)
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
    real(kind=8) :: tmin,tmax,smin,smax,phi1,phi2,gram4,V,sqrtGG,shatmin,shatmax,y,base,root,phi_rot,&
         etminir,etmini
    real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
    real(kind=8),external :: ran2
    if (popcnt(i).gt.1) then
       if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
       call shatminmax(i,ir,shatmin,shatmax)
       call generate_mass(i,shatmin,shatmax)
    endif
    if (popcnt(ir).gt.1) then
       call shatminmax(ir,i,shatmin,shatmax)
       call generate_mass(ir,shatmin,shatmax)
    endif
    if (jac.le.0d0) return
    if (debug) then
       write (*,*) '23- i    ',i,invm(i)
       write (*,*) '23- ir   ',ir,invm(ir)
    endif

    call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin,tmax)
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    ! Make sure that the t-range is compatible with the pT cut. Since t is an
    ! invariant we can compute it in any frame. Let's use the frame in which
    ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
    ! ETmin() is boost invariant in the z-direction)
    pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
    y=log((pp(0,i+ir)+pp(3,i+ir))/(pp(0,i+ir)-pp(3,i+ir)))/2d0
    call boostz(pp(0,i+ir),y,piir)
    call boostz(pp(0,ib),y,pib)
    if ( piir(1)**2+piir(2)**2.lt.etmin(i)**2-invm(i) .and. popcnt(i).eq.1 ) then
       etminir=max(etmin(ir),sqrt(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(i)-invm(i)))**2) )
    else
       etminir=max(etmin(ir),sqrt(invm(ir)))
    endif
    if ( piir(1)**2+piir(2)**2.lt.etmin(ir)**2-invm(ir) .and. popcnt(ir).eq.1 ) then
       etmini=max(etmin(i),sqrt(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(ir)-invm(ir)))**2) )
    else
       etmini=max(etmin(i),sqrt(invm(i)))
    endif
    base=piir(0)**2-ETmini**2+ETminir**2
    ! Note, root=lambda(piir(0)**2,Etmin(i)**2,Etmin(ir)**2), but the
    ! following is more stable:
    root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
         (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
    if (root.lt.0d0) then
       jac=-33d0
       num_error=num_error+1
       if (debug) write (*,*) 'root.lt.0d0',root
       return
    endif
    tmin=max(tmin,invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
    tmax=min(tmax,invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
    if (tmin.ge.tmax) then
       jac=-3d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ir+ib),jac)
    if (debug) then
       write (*,*) '23- ir+ib',ir+ib,invm(ir+ib),tmin,tmax
    endif

    call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
         &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin,smax,V,sqrtGG)
    if (invm_min(i+im1).ne.0d0) smin=max(smin,invm_min(i+im1))
    if (invm_max(i+im1).ne.0d0) smax=min(smax,invm_max(i+im1))

    if (im1.gt.2) then
! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
    y=log((pp(0,im1)+pp(3,im1))/(pp(0,im1)-pp(3,im1)))/2d0
    call boostz(pp(0,i+ir),y,piirr)
    call boostz(pp(0,im1),y,pim1r)
    call boostz(pp(0,ib),y,pib)
    phi_rot=atan(pp(2,im1)/pp(1,im1))
    if(pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
    call rotz(piirr,-phi_rot,piir)
    call rotz(pim1r,-phi_rot,pim1)
    ! Eir > Etmin(ir) + constraint coming from t
    etminir=max(pib(0)*etmin(ir)**2/(invm(ir)-invm(ir+ib))+(invm(ir)-invm(ir+ib))/(4d0*pib(0)),Etmin(ir))
    smax=min(smax,invm(i)+invm(im1)+2d0*(piir(0)-etminir)*pim1(0)+2d0*sqrt((piir(0)-etminir)**2-invm(i))*pim1(1))

    if(invm(i).eq.0d0) then
       smin=max(smin,2d0*ETmin(i)*(pim1(0)-pim1(1)*cos(drcut)))
    endif

    endif
    
    if (smin.ge.smax) then
       jac=-4d0
       num_error=num_error+1
       if (debug) write (*,*) 'smin.ge.smax',smin,smax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,smin,smax,invm(i+im1),jac)

    if (debug) then
       write (*,*) '23- i+im1',i+im1,invm(i+im1),smin,smax
    endif

    ! Generate the momenta from the integration variables. Since there is an
    ! ambiguity in phi, get both of them and pick the one that passes the cuts
    ! (if it's only one). If both pass, simply pick one of the two at random
    ! with a flat prior.
    phi1=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,1d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi1 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi1,ppibir1)
    pr1(0:3)=pp(0:3,ir+i)-pi1(0:3)
    phi2=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,0d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi2 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi2,ppibir2)
    pr2(0:3)=pp(0:3,ir+i)-pi2(0:3)
    if ( pi1(0)**2-pi1(3)**2.ge.ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.ETmin(ir)**2 .and. &
         pi2(0)**2-pi2(3)**2.ge.ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.ETmin(ir)**2 ) then
       if(ran2().gt.0.5d0) then
          pp(0:3,i)=pi1(0:3)
          pp(0:3,ir)=pr1(0:3)
          pp(0:3,ib+ir)=ppibir1(0:3)
       else
          pp(0:3,i)=pi2(0:3)
          pp(0:3,ir)=pr2(0:3)
          pp(0:3,ib+ir)=ppibir2(0:3)
       endif
    elseif (pi1(0)**2-pi1(3)**2.ge.ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.ETmin(ir)**2) then
       pp(0:3,i)=pi1(0:3)
       pp(0:3,ir)=pr1(0:3)
       pp(0:3,ib+ir)=ppibir1(0:3)
       jac=jac/2d0
    elseif (pi2(0)**2-pi2(3)**2.ge.ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.ETmin(ir)**2) then
       pp(0:3,i)=pi2(0:3)
       pp(0:3,ir)=pr2(0:3)
       pp(0:3,ib+ir)=ppibir2(0:3)
       jac=jac/2d0
    else
       jac=-19d0
       if (debug) then
          write (*,*) 'piir',pp(0:3,i+ir)
          write (*,*) 'pim1',pp(0:3,im1)
          write (*,*) '1:',phi1,(phi1+phi2)/(2d0*pi)
          write (*,*) 'i',i,ETmin(i),':',pi1(0:3)
          write (*,*) 'ir',ir,ETmin(ir),':',pr1(0:3)
          write (*,*) '2:',phi2
          write (*,*) 'i',i,ETmin(i),':',pi2(0:3)
          write (*,*) 'ir',ir,ETmin(ir),':',pr2(0:3)
          write (*,*) ''
       endif
       return
    endif
    
    ! Compute the Jacobian
    gram4=gram_determinant4(invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib)&
         &,invm(ir+i),invm(i+im1),invm(ir+ib+i+im1),invm(ir),invm(i)&
         &,invm(im1))
    if (gram4.ge.0d0) then 
       write (*,*) 'error, gram4 greater than or equal to zero',gram4,i,ir
       jac=-5d0
       num_error=num_error+1
       return
    endif
    jac=jac/(8d0*sqrt(-gram4))
  end subroutine gen23_one_step


    
  subroutine gen23_one_step_v2(i,ir,ib,im1)
! Generates one step using the 2->3 setup, using the invariants
! shat(i), s(i) and t(i) as defined in E.~Byckling and K.~Kajantie,
! ``Reductions of the phase-space integral in terms of simpler
! processes,'' Phys. Rev. 187 (1969), 2008-2016,
! doi:10.1103/PhysRev.187.2008.  Assumes massless incoming particles.
    implicit none
    integer(kind=4),intent(in) :: im1,i,ir,ib
    real(kind=8) :: tmin,tmax,smin,smax,phi1,phi2,gram4,V,sqrtGG,shatmin,shatmax,y,base,root,phi_rot,&
         etminir,etmini
    real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
    real(kind=8),external :: ran2
    integer :: ic,irc,ibc,im1c
    common /current_step/ ic,irc,ibc,im1c

    integer,parameter :: n_try=1000
    integer :: nb,icode
    real(kind=8),parameter :: xacc=1d-8
    real(kind=8),dimension(n_try) :: xbb1,xbb2

    ic=i ; irc=ir; ibc=ib; im1c=im1
    if (popcnt(i).gt.1) then
       if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
       call shatminmax(i,ir,shatmin,shatmax)
       call generate_mass(i,shatmin,shatmax)
    endif
    if (debug) then
       write (*,*) '23- i    ',i,invm(i),shatmin,shatmax
    endif
    if (popcnt(ir).gt.1) then
       call shatminmax(ir,i,shatmin,shatmax)
       call generate_mass(ir,shatmin,shatmax)
    endif
    if (jac.le.0d0) return
    if (debug) then
       write (*,*) '23- ir   ',ir,invm(ir),shatmin,shatmax,sqrt(invm(ir))
    endif

    call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin,tmax)
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    ! Make sure that the t-range is compatible with the pT cut. Since t is an
    ! invariant we can compute it in any frame. Let's use the frame in which
    ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
    ! ETmin() is boost invariant in the z-direction)
    pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
    y=log((pp(0,i+ir)+pp(3,i+ir))/(pp(0,i+ir)-pp(3,i+ir)))/2d0
    call boostz(pp(0,i+ir),y,piir)
    call boostz(pp(0,ib),y,pib)
    
    if ( piir(1)**2+piir(2)**2.lt.etmin(i)**2-invm(i) .and. popcnt(i).eq.1 ) then
       etminir=max(etmin(ir),sqrt(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(i)-invm(i)))**2) )
    else
       etminir=max(etmin(ir),sqrt(invm(ir)))
    endif
    if ( piir(1)**2+piir(2)**2.lt.etmin(ir)**2-invm(ir) .and. popcnt(ir).eq.1 ) then
       etmini=max(etmin(i),sqrt(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(ir)-invm(ir)))**2) )
    else
       etmini=max(etmin(i),sqrt(invm(i)))
    endif
    base=piir(0)**2-ETmini**2+ETminir**2
    ! Note, root=lambda(piir(0)**2,Etmin(i)**2,Etmin(ir)**2), but the
    ! following is more stable:
    root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
         (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
    if (root.lt.0d0) then
       jac=-33d0
       num_error=num_error+1
       if (debug) write (*,*) 'root.lt.0d0',root
       return
    endif
    tmin=max(tmin,invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
    tmax=min(tmax,invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
    if (tmin.ge.tmax) then
       jac=-3d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ir+ib),jac)
    if (debug) then
       write (*,*) '23- ir+ib',ir+ib,invm(ir+ib),tmin,tmax
    endif

    call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
         &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin,smax,V,sqrtGG)
    if (invm_min(i+im1).ne.0d0) smin=max(smin,invm_min(i+im1))
    if (invm_max(i+im1).ne.0d0) smax=min(smax,invm_max(i+im1))

! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
    y=log((pp(0,im1)+pp(3,im1))/(pp(0,im1)-pp(3,im1)))/2d0
    call boostz(pp(0,i+ir),y,piirr)
    call boostz(pp(0,im1),y,pim1r)
    call boostz(pp(0,ib),y,pib)
    phi_rot=atan(pp(2,im1)/pp(1,im1))
    if(pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
    call rotz(piirr,-phi_rot,piir)
    call rotz(pim1r,-phi_rot,pim1)
    ! Eir > Etmin(ir) + constraint coming from t
    etminir=max(pib(0)*etmin(ir)**2/(invm(ir)-invm(ir+ib))+(invm(ir)-invm(ir+ib))/(4d0*pib(0)),Etmin(ir))
    smax=min(smax,invm(i)+invm(im1)+2d0*(piir(0)-etminir)*pim1(0)+2d0*sqrt((piir(0)-etminir)**2-invm(i))*pim1(1))

    if(invm(i).eq.0d0) then
       smin=max(smin,2d0*ETmin(i)*(pim1(0)-pim1(1)*cos(drcut)))
    endif
    
    if (smin.ge.smax) then
       jac=-4d0
       num_error=num_error+1
       if (debug) write (*,*) 'smin.ge.smax',smin,smax
       return
    endif

    smin=smin*1.00000001d0
    smax=smax*0.99999999d0

    nb=n_try
    call zbrak(smin_constraint,smin,smax,n_try,xbb1,xbb2,nb)
    
    if (nb.ge.1) then
       if (smin_constraint(smin).lt.0d0) then
          root=rtbis(smin_constraint,xbb1(1),xbb2(1),xacc,icode)
          if (icode.ge.0) smin=max(smin,root)
       endif
       if (smin_constraint(smax).lt.0d0) then
          root=rtbis(smin_constraint,xbb1(nb),xbb2(nb),xacc,icode)
          if (icode.ge.0) smax=min(smax,root)
       endif
    endif
    
    ix=ix+1
    call random_to_var(x(ix),ip,smin,smax,invm(i+im1),jac)
    
    if (debug) then
       write (*,*) '23- i+im1',i+im1,invm(i+im1),smin,smax
    endif

    ! Generate the momenta from the integration variables. Since there is an
    ! ambiguity in phi, get both of them and pick the one that passes the cuts
    ! (if it's only one). If both pass, simply pick one of the two at random
    ! with a flat prior.
    phi1=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,1d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi1 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi1,ppibir1)
    pr1(0:3)=pp(0:3,ir+i)-pi1(0:3)
    phi2=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,0d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi2 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi2,ppibir2)
    pr2(0:3)=pp(0:3,ir+i)-pi2(0:3)
    if ( pi1(0)**2-pi1(3)**2.ge.ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.ETmin(ir)**2 .and. &
         pi2(0)**2-pi2(3)**2.ge.ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.ETmin(ir)**2 ) then
       if(ran2().gt.0.5d0) then
          pp(0:3,i)=pi1(0:3)
          pp(0:3,ir)=pr1(0:3)
          pp(0:3,ib+ir)=ppibir1(0:3)
       else
          pp(0:3,i)=pi2(0:3)
          pp(0:3,ir)=pr2(0:3)
          pp(0:3,ib+ir)=ppibir2(0:3)
       endif
    elseif (pi1(0)**2-pi1(3)**2.ge.ETmin(i)**2 .and. pr1(0)**2-pr1(3)**2.ge.ETmin(ir)**2) then
       pp(0:3,i)=pi1(0:3)
       pp(0:3,ir)=pr1(0:3)
       pp(0:3,ib+ir)=ppibir1(0:3)
       jac=jac/2d0
    elseif (pi2(0)**2-pi2(3)**2.ge.ETmin(i)**2 .and. pr2(0)**2-pr2(3)**2.ge.ETmin(ir)**2) then
       pp(0:3,i)=pi2(0:3)
       pp(0:3,ir)=pr2(0:3)
       pp(0:3,ib+ir)=ppibir2(0:3)
       jac=jac/2d0
    else
       jac=-19d0
       if (debug) then
          write (*,*) 'nb',nb
          write (*,*) 'piir',pp(0:3,i+ir)
          write (*,*) 'pim1',pp(0:3,im1)
          write (*,*) '1:',phi1,(phi1+phi2)/(2d0*pi)
          write (*,*) 'i',i,ETmin(i),':',pi1(0:3)
          write (*,*) 'ir',ir,ETmin(ir),':',pr1(0:3)
          write (*,*) '2:',phi2
          write (*,*) 'i',i,ETmin(i),':',pi2(0:3)
          write (*,*) 'ir',ir,ETmin(ir),':',pr2(0:3)
          write (*,*) ''
       endif
       return
    endif
    
    ! Compute the Jacobian
    gram4=gram_determinant4(invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib)&
         &,invm(ir+i),invm(i+im1),invm(ir+ib+i+im1),invm(ir),invm(i)&
         &,invm(im1))
    if (gram4.ge.0d0) then 
       write (*,*) 'error, gram4 greater than or equal to zero',gram4,i,ir
       jac=-5d0
       num_error=num_error+1
       return
    endif
    jac=jac/(8d0*sqrt(-gram4))


    
  end subroutine gen23_one_step_v2


  real(kind=8) function smin_constraint(siim1)
    implicit none
    real(kind=8),intent(in) :: siim1
    real(kind=8) :: smin,smax,V,sqrtGG,phi1,phi2
    real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2
    integer :: i,ir,ib,im1
    common /current_step/ i,ir,ib,im1
    
    call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
         &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin,smax,V,sqrtGG)
    ! Generate the momenta from the integration variables. Since there is an
    ! ambiguity in phi, get both of them and pick the one that passes the cuts
    ! (if it's only one). If both pass, simply pick one of the two at random
    ! with a flat prior.
    phi1=getphifroms(siim1,invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,1d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi1 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi1,ppibir1)
    pr1(0:3)=pp(0:3,ir+i)-pi1(0:3)
    phi2=getphifroms(siim1,invm(ir+i),invm(ir),invm(ir+i+im1)&
         &,invm(ir+i+ib),V,sqrtGG,0d0)
    call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi2 &
         &,sqrt(invm(i)),sqrt(invm(ir)),pi2,ppibir2)
    pr2(0:3)=pp(0:3,ir+i)-pi2(0:3)
    smin_constraint= &
         max(min(pi1(0)**2-pi1(3)**2-ETmin(i)**2,pr1(0)**2-pr1(3)**2-ETmin(ir)**2), &
             min(pi2(0)**2-pi2(3)**2-ETmin(i)**2,pr2(0)**2-pr2(3)**2-ETmin(ir)**2))
  end function smin_constraint
  

subroutine genpt_one_step(i,ir,ib,im1)
    ! This subroutines assumes that all particles in 'ir' are massless. 
    implicit none
    integer(kind=4),intent(in) :: i,ir,ib,im1
    real(kind=8) :: pt2min,pt2max,phimin,phimax,y,shatmin,shatmax,pt2,phi,phi_rot,&
         xjac,cosphi,pt,root,denom,base,pre,ptiir
    real(kind=8),dimension(0:3) :: pim1,piir,pip,pim,prp,prm,pipr
    real(kind=8),external :: ran2
    logical :: use_plus
    if (invm(i).ne.0d0) then
       write (*,*) 'genpt_one_step only for massless particles',i,invm(i)
       stop 1
    endif
    ! get the energy in the frame where p(:,i+ir) has p_z=0.
    y=log((pp(0,i+ir)+pp(3,i+ir))/(pp(0,i+ir)-pp(3,i+ir)))/2d0
    call boostz(pp(0,i+ir),y,piir)
    call boostz(pp(0,im1),y,pim1)
    ptiir=sqrt(piir(1)**2+piir(2)**2)
    ! generate pT^2
    pt2min=ETmin(i)**2-invm(i)
    pt2max=min((piir(0)-ETmin(ir))**2-invm(i),0.25d0*(piir(0)+ptiir)**2)
    if (pt2min.gt.pt2max) then
       if (debug) write (*,*) 'pt2min,pt2max',pt2min,pt2max
       jac=-14d0
       return
    endif

    ix=ix+1
    call random_to_var(x(ix),ip,pt2min,pt2max,pt2,jac)
    pt=sqrt(pt2)

    if (debug) then
       write (*,*) 'pt2 - i  ',i,pt2,pt2min,pt2max
    endif

    ! generate phi
    phimax=(piir(1)**2+piir(2)**2+pT2-(piir(0)-sqrt(pT2+invm(i)))**2) &
         /(2d0*ptiir*pT)
    if (phimax.gt.1d0) then
       write (*,*) 'ERROR,phimax',phimax,sqrt(piir(1)**2+piir(2)**2),pT,piir(0)
       stop 1
    elseif (phimax.gt.-1d0) then
       phimax=acos(phimax)
    else
       phimax=pi
    endif
    phimin=0d0

    ix=ix+1
    call random_to_var(x(ix),0d0,phimin,phimax,phi,jac)
    if (ran2().lt.0.5d0) phi=-phi
    jac=jac*2d0
    
    if (debug) then
       write (*,*) 'phi - i  ',i,phi,phimin,phimax
    endif

    ! phi is the angle between pT(i+ir) and pT(i), but it needs to be the
    ! angle between pT(im1) and pT(i). Hence, compensate:
    cosphi=(piir(1)*pim1(1)+piir(2)*pim1(2))/&
         (ptiir*sqrt(pim1(1)**2+pim1(2)**2))
    if (cosphi.gt.1d0 .and. cosphi.lt.1d0+tiny) cosphi=1d0
    if (cosphi.lt.-1d0 .and. cosphi.gt.-1d0-tiny) cosphi=-1d0
    phi_rot=acos(cosphi)
    if(piir(1)*pim1(2).lt.piir(2)*pim1(1)) phi_rot=-phi_rot
    phi=phi-phi_rot
    if (phi.lt.-pi) phi=phi+2d0*pi
    if (phi.gt.+pi) phi=phi-2d0*pi

    ! boost to the frame where p(:,im1) has p_z=0.
    y=log((pp(0,im1)+pp(3,im1))/(pp(0,im1)-pp(3,im1)))/2d0
    call boostz(pp(0,im1),y,pim1)
    call boostz(pp(0,i+ir),y,piir)

    shatmin=2d0*pt*(pim1(0)-sqrt(pim1(1)**2+pim1(2)**2)*cos(max(drcut,abs(phi))))

    pre=pim1(0)*piir(0)-2d0*cos(phi)*pim1(0)*pt
    base=-pim1(0)*piir(0)*(ETmin(ir)-pt)*(ETmin(ir)+pt)
    root=abs(pim1(0)*piir(3))*sqrt((piir(0)**2-(etmin(ir)-pt)**2-piir(3)**2)*&
         (piir(0)**2-(etmin(ir)+pt)**2-piir(3)**2))
    denom=(piir(0)-piir(3))*(piir(0)+piir(3))
    
    if ( pt.gt.piir(0)-abs(piir(3)) .or. &
         ( pt.lt.piir(0)-abs(piir(3)) .and. &
           ETmin(ir)**2.gt.(piir(0)-pt)**2-piir(3)**2 ) ) then
       shatmin=max(shatmin,pre+(base-root)/denom)
    elseif (pt2.gt.denom) then
       write (*,*) 'pT2 too large',pt2,piir(0)**2-piir(3)**2
       stop 1
    elseif (Etmin(ir).gt.sqrt(piir(0)**2-piir(3)**2)-pt) then
       write (*,*) 'Not enough energy',Etmin(ir),sqrt(piir(0)**2-piir(3)**2)-pt
       stop 1
    endif

    shatmax=pre+(base+root)/denom
    
    if (shatmin.gt.shatmax) then
       if (debug) write (*,*) shatmin,shatmax,2d0*pt*pim1(0)*(1d0-cos(drcut))
       jac=-13d0
       num_error=num_error+1
       return
    endif
      
    ix=ix+1
    call random_to_var(x(ix),ip,shatmin,shatmax,invm(i+im1),jac)
    
    if (debug) then
       write (*,*) 'shat - i+im1',i+im1,invm(i+im1),shatmin,shatmax
    endif
    
    ! fill momentum, assuming that previous particle is along the x-axis.
    call fill_momentum_ptinvmphi(pt,invm(i+im1),phi,pim1(0),pipr(0),xjac)
    if (xjac.lt.0d0) then
       jac=xjac
       return
    endif
    jac=jac*xjac

    ! rotate about the z-axis
    phi_rot=atan(pp(2,im1)/pp(1,im1))
    if(pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
    call rotz(pipr,phi_rot,pip)

    ! check both +pz and -pz
    pim(0:2)=pip(0:2)
    pim(3)=-pip(3)
    prp(0:3)=piir(0:3)-pip(0:3)
    prm(0:3)=piir(0:3)-pim(0:3)
    if ( prp(0)**2-prp(3)**2 .ge. Etmin(ir)**2 .and. &
         prp(0)**2-prp(3)**2 .ge. Etmin(ir)**2 .and. &
         dot(prp(0:3),prp(0:3)) .ge. 0d0 .and. &
         dot(prm(0:3),prm(0:3)) .ge. 0d0 &
         ) then
       if (ran2().gt.0.5d0) then
          use_plus=.true.
       else
          use_plus=.false.
       endif
       jac=jac*2d0
    elseif ( prp(0)**2-prp(3)**2 .ge. Etmin(ir)**2 .and. &
         dot(prp(0:3),prp(0:3)) .ge. 0d0 &
         ) then
       use_plus=.true.
    elseif ( prm(0)**2-prm(3)**2 .ge. Etmin(ir)**2 .and. &
         dot(prm(0:3),prm(0:3)) .ge. 0d0 &
         ) then
       use_plus=.false.
    else
       ! The constraints shatmin/shatmax are such that
       ! E(ir)>ETmin(ir). However, sometimes, the invariant mass of ir
       ! is smaller than zero for both configurations (i.e., for both
       ! +pz and -pz). This boundary has not been implemented
       ! consistently. In that case we simply return.
       jac=-21d0
       num_error=num_error+1
       return
    endif
    
    ! boost along the z-axis
    if (use_plus) then
       call boostz(pip,-y,pp(0,i))
    else
       call boostz(pim,-y,pp(0,i))
    endif
    
    pp(0:3,ir)=pp(0:3,ir+i)-pp(0:3,i)
    jac=jac/dble(4)
    invm(ir)=dot(pp(0:3,ir),pp(0:3,ir))
    if (debug) then
       write (*,*) 'pp(ir+i+im1)',pp(0:3,ir+i+im1),invm(ir+i+im1)
       write (*,*) 'pp(im1 )    ',pp(0:3,im1),invm(im1),im1
       write (*,*) 'pp(ir+i)    ',pp(0:3,ir+i),invm(ir+i),ir+i
       write (*,*) 'pp(i)       ',pp(0:3,i),invm(i),i
       write (*,*) 'pp(ir)      ',pp(0:3,ir),invm(ir),ir
       write (*,*) ''
    endif
    if (invm(ir).lt.0d0) then
       jac=-12d0
       num_error=num_error+1
       return
    endif
    ! fill t-channel stuff to be safe.
    pp(0:3,i+ib)=pp(0:3,i)-pp(0:3,ib)
    pp(0:3,ir+ib)=pp(0:3,ir)-pp(0:3,ib)
    invm(i+ib)=dot(pp(0:3,i+ib),pp(0:3,i+ib))
    invm(ir+ib)=dot(pp(0:3,ir+ib),pp(0:3,ir+ib))
  end subroutine genpt_one_step
  subroutine fill_momentum_ptinvmphi(pt,invm,phi,Eref,p,xjac)
    implicit none
    real(kind=8) :: invm,phi,xjac,pt,Eref
    real(kind=8),dimension(0:3) :: p
    real(kind=8),external :: ran2
    p(1)=pt*cos(phi)
    p(2)=pt*sin(phi)
    p(0)=invm/(2d0*Eref)+p(1)
    ! There are two values of the pz that correspond to a single
    ! invm. Take one of the two at random.
    p(3)=p(0)**2-pt**2
    if (p(3).lt.0d0 .and. p(3).ge.-tiny) then
       p(3)=0d0
    elseif(p(3).lt.-tiny) then
       xjac=-20d0
       return
    else
       p(3)=sqrt(p(3))
    endif
    xjac=abs(1d0/(2d0*Eref*(p(3)+vtiny)))
  end subroutine fill_momentum_ptinvmphi
  subroutine boostz(p,yb,pb)
    ! boost in the z-direction with rapidity yb
    implicit none
    real(kind=8),dimension(0:3) :: p,pb
    real(kind=8) :: yb
    pb(0)=p(0)*cosh(yb)-p(3)*sinh(yb)
    pb(1:2)=p(1:2)
    pb(3)=p(3)*cosh(yb)-p(0)*sinh(yb)
  end subroutine boostz

  subroutine rotz(p,phi,prot)
    implicit none
    real(kind=8),dimension(0:3) :: p,prot
    real(kind=8) :: phi
    prot(0)=p(0)
    prot(1)=p(1)*cos(phi)-p(2)*sin(phi)
    prot(2)=p(1)*sin(phi)+p(2)*cos(phi)
    prot(3)=p(3)
  end subroutine rotz

  subroutine roty(p,phi,prot)
    implicit none
    real(kind=8),dimension(0:3) :: p,prot
    real(kind=8) :: phi
    prot(0)=p(0)
    prot(1)=p(1)*cos(phi)-p(3)*sin(phi)
    prot(2)=p(2)
    prot(3)=p(1)*sin(phi)+p(3)*cos(phi)
  end subroutine roty

  subroutine gent_one_step(i,ir,ib)
    ! One step in the usual MadGraph t-channel phase-space generation.
    implicit none
    integer(kind=4),intent(in) :: i,ir,ib
    real(kind=8) :: tmin,tmax,phi,Eimax,shatmin,shatmax,base,etminir,root,y,etmini
    real(kind=8),dimension(0:3) :: piir,pib
    if (popcnt(i).gt.1) then
       if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
       call shatminmax(i,ir,shatmin,shatmax)
       if (popcnt(i+ir).eq.next-2) then
          ! The energy of i will be
          ! Ei=(sqrtshat+(invm(i)-invm(ir))/sqrtshat)/2d0. This gives a
          ! constraint on the allowed value of invm(i), since Ei>ETmin(i)
          shatmin=max(shatmin,max(invm(ir),invm_min(ir))+sqrtshat*(2d0*ETmin(i)-sqrtshat))
          Eimax=sqrtshat-ETmin(ir) ! maximum energy for i
          shatmax=min(shatmax,max(invm(ir),invm_min(ir))+sqrtshat*(2d0*Eimax-sqrtshat))
       endif
       call generate_mass(i,shatmin,shatmax)
    endif
    if (popcnt(ir).gt.1) then
       call shatminmax(ir,i,shatmin,shatmax)
       if (popcnt(i+ir).eq.next-2) then
          ! The energy of ir will be
          ! Eir=(sqrtshat+(invm(ir)-invm(i))/sqrtshat)/2d0. This gives a
          ! constraint on the allowed value of invm(ir), since Eir>ETmin(ir)
          shatmin=max(shatmin,invm(i)+sqrtshat*(2d0*ETmin(ir)-sqrtshat))
          shatmax=min(shatmax,invm(i)+sqrtshat*(sqrtshat-2d0*max(sqrt(invm(i)),ETmin(i))))
       endif
       call generate_mass(ir,shatmin,shatmax)
    endif
    if (jac.le.0d0) return
    if (debug) then
       write (*,*) 't - i    ',i,invm(i),invm_min(i),invm_max(i)
       write (*,*) 't - ir   ',ir,invm(ir)
    endif

    call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin,tmax)
    if (invm_max(ir+ib).ne.0d0) tmax=min(tmax,invm_max(ir+ib))
    if (invm_min(ir+ib).ne.0d0) tmin=max(tmin,invm_min(ir+ib))
    ! Make sure that the t-range is compatible with the pT cut. Since t is an
    ! invariant we can compute it in any frame. Let's use the frame in which
    ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
    ! ETmin() is boost invariant in the z-direction)
    pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
    y=log((pp(0,i+ir)+pp(3,i+ir))/(pp(0,i+ir)-pp(3,i+ir)))/2d0
    call boostz(pp(0,i+ir),y,piir)
    call boostz(pp(0,ib),y,pib)
    if ( piir(1)**2+piir(2)**2.lt.etmin(i)**2-invm(i) .and. popcnt(i).eq.1 ) then
       etminir=max(etmin(ir),sqrt(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(i)-invm(i)))**2) )
    else
       etminir=max(etmin(ir),sqrt(invm(ir)))
    endif
    if ( piir(1)**2+piir(2)**2.lt.etmin(ir)**2-invm(ir) .and. popcnt(ir).eq.1 ) then
       etmini=max(etmin(i),sqrt(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-sqrt(etmin(ir)-invm(ir)))**2) )
    else
       etmini=max(etmin(i),sqrt(invm(i)))
    endif
    base=piir(0)**2-ETmini**2+ETminir**2
    ! Note, root=lambda(piir(0)**2,Etmin(i)**2,Etmin(ir)**2), but the
    ! following is more stable:
    root=(piir(0)-ETmini-ETminir)*(piir(0)+ETmini-ETminir)*&
         (piir(0)-ETmini+ETminir)*(piir(0)+ETmini+ETminir)
    if (root.lt.0d0) then
       jac=-33d0
       num_error=num_error+1
       if (debug) write (*,*) 'root.lt.0d0',root
       return
    endif
    tmin=max(tmin,invm(ir)-pib(0)/piir(0)*(base+sqrt(root)))
    tmax=min(tmax,invm(ir)-pib(0)/piir(0)*(base-sqrt(root)))
    if (tmin.ge.tmax) then
       jac=-3d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ir+ib),jac)
    
    if (debug) then
       write (*,*) 't- ir+ib',ir+ib,invm(ir+ib),tmin,tmax
    endif
    ix=ix+1
    call random_to_var(x(ix),0d0,0d0,2d0*pi,phi,jac)
    call gentcms(pp(0,ib+ir+i),pp(0,ib),invm(ib+ir),phi,sqrt(invm(i)) &
         &,sqrt(invm(ir)),pp(0,i),pp(0,ib+ir))
    pp(0:3,ir)=pp(0:3,ib+ir+i)+pp(0:3,ib)-pp(0:3,i)
    jac = jac/(4d0*sqrt(lambda(invm(ir+i),0d0,invm(ir+i+ib))))
  end subroutine gent_one_step

  subroutine gent_one_step_v2(i,ir,ib,im1)
    implicit none
    integer(kind=4),intent(in) :: i,ir,ib,im1
    real(kind=8) :: pt2,y,phi_rot,pt2min,pt2max,tmin,tmax,smin,smax,&
         ea,mi,tr,siim1,mim1,eim1,pxim1,base,root,yb
    real(kind=8),dimension(0:3) :: piirr,pim1r,pib,piir,pim1,pii,pir
    real(kind=8),external :: ran2
    
! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
    y=log((pp(0,im1)+pp(3,im1))/(pp(0,im1)-pp(3,im1)))/2d0
    call boostz(pp(0,i+ir),y,piirr)
    call boostz(pp(0,im1),y,pim1r)
    call boostz(pp(0,ib),y,pib)
    phi_rot=atan(pp(2,im1)/pp(1,im1))
    if (pp(1,im1).lt.0d0) phi_rot=phi_rot+pi
    call rotz(piirr,-phi_rot,piir)
    call rotz(pim1r,-phi_rot,pim1)

    ! generate pT^2
    pt2min=ETmin(i)**2-invm(i)
    pt2max=min((piir(0)-ETmin(ir))**2-invm(i),&
               0.25d0*(piir(0)+sqrt(piir(1)**2+piir(2)**2))**2)
    if (pt2min.gt.pt2max) then
       if (debug) write (*,*) 'pt2min,pt2max',pt2min,pt2max
       jac=-14d0
       return
    endif

    ix=ix+1
    call random_to_var(x(ix),ip,pt2min,pt2max,pt2,jac)
    
    if (debug) then
       write (*,*) 'pt2 - i  ',i,pt2,pt2min,pt2max
    endif

    Ea=pib(0)
    mi=sqrt(invm(i))
    tr=invm(ib+i)
    mim1=sqrt(invm(im1))
    eim1=pim1(0)
    pxim1=pim1(1)
    
    ! Make sure that the t-range is compatible with the pT cut. Since t is an
    ! invariant we can compute it in any frame. Let's use the frame in which
    ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
    ! ETmin() is boost invariant in the z-direction)
    pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
    yb=log((pp(0,i+ir)+pp(3,i+ir))/(pp(0,i+ir)-pp(3,i+ir)))/2d0
    call boostz(pp(0,i+ir),yb,piir)
    call boostz(pp(0,ib),yb,pib)
    base=piir(0)**2+pt2-Etmin(ir)**2+invm(i)*(1d0-piir(0)/pib(0))
    ! Note, root=lambda(piir(0)**2,pt2+invm(i),Etmin(ir)**2), but the
    ! following is more stable:
    root=(piir(0)-sqrt(pt2+invm(i))-Etmin(ir))*(piir(0)+sqrt(pt2+invm(i))-Etmin(ir))*&
         (piir(0)-sqrt(pt2+invm(i))+Etmin(ir))*(piir(0)+sqrt(pt2+invm(i))+Etmin(ir))
    if (root.lt.0d0) then
       jac=-33d0
       num_error=num_error+1
       if (debug) write (*,*) 'root.lt.0d0',root
       return
    endif
    tmin=invm(i)-pib(0)/piir(0)*(base+sqrt(root))
    tmax=invm(i)-pib(0)/piir(0)*(base-sqrt(root))
    if (tmin.ge.tmax) then
       jac=-3d0
       num_error=num_error+1
       if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
       return
    endif

    
    ! generate t
    ix=ix+1
    call random_to_var(x(ix),ip,tmin,tmax,invm(ib+i),jac)
    
    if (debug) then
       write (*,*) 'pt2 - ib+i',ib+i,invm(ib+i),tmin,tmax
    endif


    pii(0)=(Ea*(mi**2 + pt2))/(mi**2 - tr) + (mi**2 - tr)/(4.*Ea)
    pii(3)=(Ea*(mi**2 + pt2))/(mi**2 - tr) + (-mi**2 + tr)/(4.*Ea)

    

    smin=invm_min(i+im1)
    if(invm(i).eq.0d0) then
       smin=max(smin,2d0*sqrt(pt2)*(pim1(0)-pim1(1)*cos(drcut)))
    endif
    smin=max(smin,2*pii(0)*eim1+mi**2+mim1**2-2d0*sqrt(pt2)*pxim1)
    smax=2*pii(0)*eim1+mi**2+mim1**2+2d0*sqrt(pt2)*pxim1
    
    if (smin.ge.smax) then
       jac=-4d0
       num_error=num_error+1
       if (debug) write (*,*) 'smin.ge.smax',smin,smax
       return
    endif
    ix=ix+1
    call random_to_var(x(ix),ip,smin,smax,invm(i+im1),jac)
    
    if (debug) then
       write (*,*) '23- i+im1',i+im1,invm(i+im1),smin,smax
    endif


    siim1=invm(i+im1)

    pii(1)=(2*(mi**2 + mim1**2 - siim1) + &
         (4*Ea*eim1*(mi**2 + pt2))/(mi**2 - tr) + (eim1*(mi**2 - tr))/Ea)/(4.*pxim1)

    if (pii(1)**2.gt.pt2) then
       jac=-101d0
       return
    endif
    
    pii(2)=sqrt(pt2-pii(1)**2)
    if (ran2().gt.0.5d0) pii(2)=-pii(2)
    jac=jac*2d0



    
    call rotz(pii,phi_rot,pir)
    call boostz(pir,-y,pp(0,i))

    pp(0:3,ir)=pp(0:3,i+ir)-pp(0:3,i)
    invm(ir)=dot(pp(0:3,ir),pp(0:3,ir))
    pp(0:3,ir+ib)=pp(0:3,ir)-pp(0:3,ib)
    invm(ir+ib)=dot(pp(0:3,ir+ib),pp(0:3,ir+ib))

    jac=jac*(1d0/abs(8*Ea*pxim1*pii(2)*(-1d0 + pii(3)/pii(0))))*2d0/pii(0)
    jac=jac/4d0

    if (invm(ir).lt.0d0 .or. pp(0,ir).lt.ETmin(ir)) then
       jac=-102d0
       return
    endif

  end subroutine gent_one_step_v2

  double precision function rtbis(func,x1,x2,xacc,icode)
    implicit none
    integer,parameter :: jmax=1000
    real*8,intent(in) :: x1,x2,xacc
    integer,intent(out) :: icode
    real*8,external :: func
    ! Using bisection, find the root of a function func known to lie
    ! between x1 and x2. The root, returned as rtbis, will be refined
    ! until its accuracy is |func(rtbis)|<xacc.
    integer :: j
    real*8 :: dx,f,fmid,xmid
    fmid=func(x2)
    f=func(x1)
    icode=0
    if(f*fmid.ge.0.) then
!!$       write (*,*) "root must be bracketed in rtbis",x1,x2,f,fmid
       icode=-1
!!$       stop 1
    endif
    if(f.lt.0.)then ! Orient the search so that f>0 lies at x+dx.
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    endif
    do j=1,jmax ! Bisection loop.
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmid)
       if(fmid.le.0.) rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) then
          return
       endif
    enddo
    write (*,*)  "WARNING: too many bisections in rtbis"
  end function rtbis

  
  subroutine zbrak(fx,x1,x2,n,xb1,xb2,nb)
! Given a function fx defined on the interval from x1-x2 subdivide the
! interval into n equally spaced segments, and search for zero crossings of
! the function. nb is input as the maximum number of roots sought, and is
! reset to the number of bracketing pairs xb1(1:nb), xb2(1:nb) that are found.
    implicit none
    integer n,nb
    real(kind=8) :: x1,x2,xb1(nb),xb2(nb)
    real(kind=8),external :: fx
    integer :: i,nbb
    real(kind=8) dx,fc,fp,x
    nbb=0
    x=x1
    dx=(x2-x1)/n
    fp=fx(x)
    ! Determine the spacing appropriate to the mesh. Loop over all intervals
    ! If a sign change occurs then record values for the bounds.
    do i=1,n
       x=x+dx
       fc=fx(x)
       if(fc*fp.le.0d0) then
          nbb=nbb+1
          xb1(nbb)=x-dx
          xb2(nbb)=x
          if(nbb.eq.nb) exit
       endif
       fp=fc
    enddo
    nb=nbb
  end subroutine zbrak


  subroutine gentcms_v2(pa,pb,t,phi,m1,m2,p1,pr)
! Generates 4 momentum for particle p1, and remainder pr=pa-p1=p2-pb given the
! values t, and phi in the process pa+pb -> p1+p2.  Assuming incoming
! particles with momenta pa, pb and outgoing particles with mass m1,m2;
! t=(pb-p2)^2 ; phi is the azimuthal angle between p1 and pa in the
! pa+pb rest frame, with pa aligned with the positive z-axis.
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
       num_error=num_error+1
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
    call boostm(p1_rot,ptot,esum,p1)    !boost back to lab frame
    pr(0:3)=pa(0:3)-p1(0:3)         !Return remainder of momentum
  end subroutine gentcms_v2


  
  subroutine gens_one_step(i,ir)
    implicit none
    integer(kind=4),intent(in) :: i,ir
    real(kind=8) :: esum,costh,phi,shatmin,shatmax
    real(kind=8),dimension(0:3) :: p_i,p_ir
    if (popcnt(i).gt.1) then
       if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
       call shatminmax(i,ir,shatmin,shatmax)
       call generate_mass(i,shatmin,shatmax)
    endif
    if (popcnt(ir).gt.1) then
       call shatminmax(ir,i,shatmin,shatmax)
       call generate_mass(ir,shatmin,shatmax)
    endif
    if (jac.le.0d0) return
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

  subroutine generate_mass(i,shatmin,shatmax)
    implicit none
    integer :: i
    real(kind=8) :: shatmin,shatmax
    if (invm_min(i).ne.0d0) shatmin=max(shatmin,invm_min(i))
    if (invm_max(i).ne.0d0) shatmax=min(shatmax,invm_max(i))
    if (shatmin.ge.shatmax) then
       jac=-7d0
       num_error=num_error+1
       if (debug) write (*,*) 'shatmin.ge.shatmax',i,shatmin,shatmax
       return
    endif
    ix=ix+1
!!$          call random_to_var(x(ix),ip,shatmin,shatmax,invm(j1),jac)
    call random_to_var(x(ix),-0.5d0,shatmin,shatmax,invm(i),jac)
  end subroutine generate_mass
  
  
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
    real(kind=8) :: xa2,xb2,S
    lambda=s**2-2d0*(xa2+xb2)*s+(xa2-xb2)**2
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
! Generates 4 momentum for particle p1, and remainder pr=pa-p1=p2-pb given the
! values t, and phi in the process pa+pb -> p1+p2.  Assuming incoming
! particles with momenta pa, pb and outgoing particles with mass m1,m2;
! t=(pb-p2)^2 ; phi is the azimuthal angle between p1 and pa in the
! pa+pb rest frame, with pa aligned with the positive z-axis.
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
       num_error=num_error+1
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
!!$    p1_rot=p1
    call boostm(p1_rot,ptot,esum,p1)    !boost back to lab frame
    pr(0:3)=pa(0:3)-p1(0:3)         !Return remainder of momentum
  end subroutine gentcms


  subroutine gentcms2(pa,pb,pc,t,phi,m1,m2,p1,pr)
! Generates 4 momentum for particle p1, and remainder pr given the
! values t, and phi in the process pa+pb -> pr+p1.  Assumes incoming
! particles with momenta pa, pb and outgoing particles with mass
! m1,m2; t=(pb-p1)^2=(pa-pr)^2. Assumes that pa is a massless
! momentum; phi is the azimuthal angle between pr and pc in the pa+pb
! rest frame, with pa aligned with the positive z-axis.
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
       num_error=num_error+1
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
       num_error=num_error+1
       return
       stop 1
    endif
    Pii(3) = -(m2**2-t-2d0*Pii(0)*E_acms)/(2d0*p_acms)
    pt2=pp2-Pii(3)**2
    if (pt2/esum2.lt.-tiny) then
       write (*,*) 'Error #16 in genps_fks.f: relative pt^2 smaller than 0',pt2
       jac=-11d0
       num_error=num_error+1
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

  subroutine random_to_var(x,power_in,var_min,var_max,var,jac)
    ! Given a random number x, it generates var in the range var_min
    ! <= var <= var_max according to var^(power)
    implicit none
    real(kind=8),intent(in) :: x,power_in,var_min,var_max
    real(kind=8),intent(out) :: var
    real(kind=8),intent(inout) :: jac
    integer(kind=4) :: ip
    real(kind=8) :: varmin,varmax,power
    if (var_min.lt.0d0 .and. var_max.le.0d0) then
       power=power_in
       varmin=-var_max
       varmax=-var_min
    elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (abs(power_in).gt.vtiny)) then
       write (*,*) 'ERROR: in random_to_var one of the two limits '/&
            &/'is negative',var_min,var_max,power_in,jac,x
       write (*,*) 'using flat transformation'
       power=0d0
       varmin=var_min
       varmax=var_max
    else
       power=power_in
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
  
  real(kind=8) function getphifroms(si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran)
    ! Given s_i (invariant mass of p_i and p_i+1, it transforms it
    ! into phi_i. Note that there are two possibilities for phi: need
    ! to pick one at random.
    ! Based on eq.(11) of E.~Byckling and K.~Kajantie, ``Reductions of
    ! the phase-space integral in terms of simpler processes,''
    ! Phys. Rev. 187 (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran
    real(kind=8) :: cosphi,x
    cosphi=((si-shat_im1-shat_ip1)*0.5d0*lambda(shat_i,t_i,0d0)-4d0*V)/sqrtGG
    x=ran
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
       computeV=-99d99
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
    if (GG.le.0d0 .or. V.eq.-99d99) then
!!$       write (*,*) 'No allowed range for s: smin=smax',GG,V
!!$       stop 1
       GG=0d0
       V=0d0
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
    real(kind=8),parameter :: tol=1d-8
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
