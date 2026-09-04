module phase_space_genpt_mod
!  use common
  use phase_space_base
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type,extends(phase_space_type),public :: phase_space_genpt
   contains
     procedure :: init => genpt_init
     procedure :: generate_momenta => genpt_generate_momenta
     procedure :: compute_x_from_momenta => genpt_compute_x_from_momenta
     procedure :: cleanup => genpt_cleanup
  end type phase_space_genpt
  private
  real(kind=8),parameter :: pi=3.1415926535897932d0
  integer,parameter :: use_mode=1 ! all modes use pT^2 and phi, but
                                  ! 1 = uses rapidity (original chili)
                                  ! 2 = uses DeltaR with previous particle
                                  ! 3 = uses invariant mass with previous particle
                                  ! 4 = uses cos(theta) with previous particle
contains
  subroutine genpt_compute_x_from_momenta(this,ps)
    implicit none
    class(phase_space_genpt),intent(inout) :: this
    type(psv),intent(inout) :: ps
    if (allocated(ps%x)) ps%x=0d0
    ps%jac=-14d0
  end subroutine genpt_compute_x_from_momenta
  subroutine genpt_cleanup(this)
    implicit none
    class(phase_space_genpt),intent(inout) :: this
    if (allocated(this%order)) deallocate(this%order)
    if (allocated(this%masses)) deallocate(this%masses)
    if (allocated(this%invm)) deallocate(this%invm)
    if (allocated(this%invm_min)) deallocate(this%invm_min)
    if (allocated(this%invm_max)) deallocate(this%invm_max)
    if (allocated(this%ETmin)) deallocate(this%ETmin)
    if (allocated(this%pp)) deallocate(this%pp)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%sets)) deallocate(this%sets)
    if (allocated(this%ptcut)) deallocate(this%ptcut)
    if (allocated(this%ycut)) deallocate(this%ycut)
    if (allocated(this%drcut)) deallocate(this%drcut)
    if (allocated(this%sqrt_s_min)) deallocate(this%sqrt_s_min)
    this%jac=0d0
    this%xbjrk=0d0
    this%s0=0d0
    this%tot_mass=0d0
    this%sqrtshat=0d0
    this%sqrts=0d0
    this%ndim=0
    this%next=0
    this%ndim_extra=0
    this%t_channel=.false.
    this%can_invert_momenta=.false.
  end subroutine genpt_cleanup
  subroutine genpt_init(this,sqrts,n,m,o,pt_cut,rap_cut,DR_cut,sqrt_s_min,t_chan,include_pdf,flat)
    implicit none
    class(phase_space_genpt),intent(inout) :: this
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrts
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! rapidity and pT cut (and DR and sqrt_s_min) on all the particles
    real(kind=8),dimension(n),intent(in) :: rap_cut,pt_cut
    real(kind=8),dimension(n,n),intent(in) :: DR_cut,sqrt_s_min
    ! masses of all the particles. The two incoming particles must be
    ! massless.
    real(kind=8),dimension(n),intent(in) :: m
    integer,dimension(n),intent(in) :: o
    logical,intent(in) :: include_pdf,t_chan
    logical,intent(in),optional :: flat
    integer :: i,j,allocation_status
    real(kind=8) :: maximum_rapidity
    character(len=256) :: allocation_message
    call this%cleanup()
    if (n.lt.4) then
       write (*,*) 'genpt phase-space requires at least two final-state particles'
       stop 1
    endif
    if (n.gt.max_bitmask_particles) then
       write (*,*) 'genpt phase-space particle multiplicity exceeds the supported limit:',&
            n,max_bitmask_particles
       stop 1
    endif
    if (.not.ieee_is_finite(sqrts)) then
       write (*,*) 'genpt phase-space received an invalid collision energy:',sqrts
       stop 1
    endif
    if (sqrts.le.sqrt(tiny(1d0)) .or. sqrts.gt.phase_space_momentum_limit) then
       write (*,*) 'genpt phase-space received an invalid collision energy:',sqrts
       stop 1
    endif
    if (.not.all(ieee_is_finite(m)) .or. .not.all(ieee_is_finite(pt_cut)) .or. &
         .not.all(ieee_is_finite(rap_cut)) .or. .not.all(ieee_is_finite(DR_cut)) .or. &
         .not.all(ieee_is_finite(sqrt_s_min))) then
       write (*,*) 'genpt phase-space received non-finite masses or cuts'
       stop 1
    endif
    do i=1,n
       if (count(o.eq.i).ne.1) then
          write (*,*) 'genpt phase-space colour order is not a permutation:',o
          stop 1
       endif
    enddo
    if (max(maxval(abs(m)),maxval(abs(pt_cut)),maxval(abs(sqrt_s_min))).gt. &
         phase_space_momentum_limit) then
       write (*,*) 'genpt phase-space mass or dimensionful cut exceeds numerical range'
       stop 1
    endif
    allocate(this%ptcut(1:n),this%ycut(1:n),this%DRcut(n**2),&
         this%p(0:3,n),this%masses(n),this%order(n),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       call this%cleanup()
       write (*,*) 'genpt phase-space could not allocate its workspace:',&
            trim(allocation_message)
       stop 1
    endif
    this%ptcut=pt_cut
    this%ycut=rap_cut
    do i=1,n
       do j=1,n
          this%DRcut(n*(i-1)+j)=DR_cut(i,j)
       enddo
    enddo
    this%sqrts=sqrts
    this%next=n
    this%ndim=3*n-8
    this%can_invert_momenta=.false.
    if (.not.include_pdf) then
       write (*,*) 'genpt phase-space only with include_pdf=.true.'
       stop 1
    endif
    if (any(m(1:n).ne.0d0)) then
       write (*,*) 'genpt phase-space only for all massless particles'
       stop 1
    endif
    if (use_mode.eq.1 .and. any(rap_cut(3:n).le.0d0)) then
       write (*,*) 'genpt phase-space must have rapidity cut'
       stop 1
    endif
    do i=3,n
       do j=3,n
          if (i.eq.j) cycle
          if (use_mode.ne.1 .and. this%DRcut(n*(i-1)+j).le.0d0) then
             write (*,*) 'genpt phase-space must have DeltaR cut'
             stop 1
          endif
       enddo
    enddo
    if (any(this%ptcut(3:n).le.sqrt(tiny(1d0)))) then
       write (*,*) 'genpt phase-space must have pT cut'
       stop 1
    endif
    if (any(this%ptcut(3:n).ge.0.5d0*sqrts)) then
       write (*,*) 'genpt phase-space pT cuts must be below half the collision energy'
       stop 1
    endif
    ! A configured rapidity interval may contain both safe and unsafe points,
    ! especially under an overall rescaling of sqrt(s).  Keep the input range
    ! representable here and reject individual samples below using their actual
    ! transverse momentum; tying this bound directly to sqrt(s) breaks scale
    ! covariance even for perfectly safe generated points.
    maximum_rapidity=log(2d0*phase_space_momentum_limit/sqrt(tiny(1d0)))
    if (any(abs(this%ycut(3:n)).ge.maximum_rapidity)) then
       write (*,*) 'genpt phase-space rapidity range exceeds its safe numerical range'
       stop 1
    endif
    this%p=0d0
    this%masses=m
    this%order=o
    select case (use_mode)
    case (1)
       this%ndim_extra=0
    case (2,3)
       this%ndim_extra=max(0,n-3)
    case (4)
       this%ndim_extra=1
    case default
       write (*,*) 'genpt phase-space has an invalid internal mode:',use_mode
       stop 1
    end select
  end subroutine genpt_init

  subroutine genpt_generate_momenta(this,ps)
    implicit none
    class(phase_space_genpt),intent(inout) :: this
    type(psv),intent(inout) :: ps
    real(kind=8) :: pt2min,pt2max,pt2,ymin,ymax,y,phimin,phimax,phi&
         &,drmin,drmax,dr,phi2,invmmin,invmmax,invm,tau
    real(kind=8) :: costheta,costhetamin,costhetamax,theta,pzmax
    integer(kind=4) :: i,ix,extra_index
    real(kind=8),dimension(0:3) :: ptot,pb
    if (allocated(ps%p)) ps%p=0d0
    ps%xbjrk=0d0
    ps%jac=-10d0
    if (.not.allocated(ps%p) .or. .not.allocated(ps%x)) return
    if (size(ps%p,1).ne.4 .or. size(ps%p,2).ne.this%next .or. &
         size(ps%x).lt.this%ndim+this%ndim_extra) return
    if (.not.all(ieee_is_finite(ps%x(1:this%ndim+this%ndim_extra)))) return
    if (any(ps%x(1:this%ndim+this%ndim_extra).lt.0d0) .or. &
         any(ps%x(1:this%ndim+this%ndim_extra).gt.1d0)) return
    ix=0
    extra_index=0
    ps%jac=1d0
    do i=3,this%next-1
       ! generate pT^2
       pt2min=this%ptcut(i)**2
       pt2max=this%sqrts**2/4d0
       ix=ix+1
       call random_to_var(ps%x(ix),-1.5d0,pt2min,pt2max,pt2,ps%jac)
       if (.not.ieee_is_finite(pt2) .or. .not.ieee_is_finite(ps%jac)) then
          call reject_point(-20d0)
          return
       endif
       if (pt2.le.0d0 .or. ps%jac.le.0d0) then
          call reject_point(-20d0)
          return
       endif
       ! generate phi
       phimin=-pi
       phimax=pi
       ix=ix+1
       call random_to_var(ps%x(ix),0d0,phimin,phimax,phi,ps%jac)
       if (.not.ieee_is_finite(phi) .or. .not.ieee_is_finite(ps%jac)) then
          call reject_point(-20d0)
          return
       endif
       if (ps%jac.le.0d0) then
          call reject_point(-20d0)
          return
       endif
       if (use_mode.eq.1 .or.i.eq.3) then
          ! generate rapidity
          ymin=-this%ycut(i)
          ymax=+this%ycut(i)
          ix=ix+1
          call random_to_var(ps%x(ix),0d0,ymin,ymax,y,ps%jac)
          if (.not.ieee_is_finite(y) .or. .not.ieee_is_finite(ps%jac)) then
             call reject_point(-20d0)
             return
          endif
          if (ps%jac.le.0d0) then
             call reject_point(-20d0)
             return
          endif
          if (abs(y).ge.log(2d0*phase_space_momentum_limit/sqrt(pt2))) then
             call reject_point(-20d0)
             return
          endif
          ! fill momentum
          call fill_momentum_pt2yphi(pt2,y,phi,ps%p(0,i))

       elseif (use_mode.eq.2) then
          ! generate deltaR w.r.t. previously generated particle
          y=log((ps%p(0,i-1)+ps%p(3,i-1))/(ps%p(0,i-1)-ps%p(3,i-1)))/2d0
          drmin=max(this%DRcut(this%next*(i-1)+i-1),abs(phi))
          drmax=sqrt((this%ycut(i)+abs(y))**2+phi**2)
          ix=ix+1
          call random_to_var(ps%x(ix),0d0,drmin,drmax,dr,ps%jac)
          ! fill momentum, assuming that previous particle is along the x-axis.
          call fill_momentum_pt2drphi(pt2,dr,phi,ps%p(0,i),ps%jac)
          ! boost along the z-axis
          call boostz(ps%p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(ps%p(2,i-1)/ps%p(1,i-1))
          if(ps%p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,ps%p(0,i))

       elseif (use_mode.eq.3) then
          ! get the energy in the frame where p(:,i-1) has p_z=0.
          y=log((ps%p(0,i-1)+ps%p(3,i-1))/(ps%p(0,i-1)-ps%p(3,i-1)))/2d0
          call boostz(ps%p(0,i-1),y,pb)
          invmmin=2d0*sqrt(pt2)*pb(0)*(1d0-cos(max(this%DRcut(this%next*(i-1)+i-1),abs(phi))))
          invmmax=this%sqrts**2
          ix=ix+1
          call random_to_var(ps%x(ix),-1d0,invmmin,invmmax,invm,ps%jac)
          ! fill momentum, assuming that previous particle is along the x-axis.
          call fill_momentum_pt2invmphi(pt2,invm,phi,pb(0),ps%p(0,i),ps%jac)
          if (ps%jac.lt.0d0) return
          ! boost along the z-axis
          call boostz(ps%p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(ps%p(2,i-1)/ps%p(1,i-1))
          if(ps%p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,ps%p(0,i))

       elseif (use_mode.eq.4) then
          y=log((ps%p(0,i-1)+ps%p(3,i-1))/(ps%p(0,i-1)-ps%p(3,i-1)))/2d0
          costhetamin=-1d0 
          costhetamax=cos(this%DRcut(this%next*(i-1)+i-1))
          ix=ix+1
          call random_to_var(ps%x(ix),0d0,costhetamin,costhetamax,costheta,ps%jac)
          call fill_momentum_pt2cosphi(pt2,costheta,phi,ps%p(0,i),ps%jac)
          if (ps%jac.lt.0d0) return
          call boostz(ps%p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(ps%p(2,i-1)/ps%p(1,i-1))
          if(ps%p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,ps%p(0,i))

       endif
    enddo


!!$    if (.true.) then
    if (use_mode.eq.1.or.use_mode.eq.4) then
       ! final particle: generate rapidity
       ymin=-this%ycut(this%next)
       ymax=+this%ycut(this%next)
       ix=ix+1
       call random_to_var(ps%x(ix),0d0,ymin,ymax,y,ps%jac)
       if (.not.ieee_is_finite(y) .or. .not.ieee_is_finite(ps%jac)) then
          call reject_point(-20d0)
          return
       endif
       if (ps%jac.le.0d0) then
          call reject_point(-20d0)
          return
       endif
       ! final particle: fill momentum
       ps%p(1,this%next)=-sum(ps%p(1,3:this%next-1))
       ps%p(2,this%next)=-sum(ps%p(2,3:this%next-1))
       pt2=hypot(ps%p(1,this%next),ps%p(2,this%next))
       if (.not.ieee_is_finite(pt2)) then
          call reject_point(-20d0)
          return
       endif
       if (pt2.le.sqrt(tiny(1d0))) then
          call reject_point(-20d0)
          return
       endif
       if (abs(y).ge.log(2d0*phase_space_momentum_limit/pt2)) then
          call reject_point(-20d0)
          return
       endif
       ps%p(3,this%next)=pt2*sinh(y)
       ps%p(0,this%next)=pt2*cosh(y)
       pt2=pt2*pt2

    elseif(use_mode.eq.2) then
       ps%p(1,this%next)=-sum(ps%p(1,3:this%next-1))
       ps%p(2,this%next)=-sum(ps%p(2,3:this%next-1))
       pt2=ps%p(1,this%next)**2+ps%p(2,this%next)**2
       phi=atan(ps%p(2,this%next-1)/ps%p(1,this%next-1))
       if(ps%p(1,this%next-1).lt.0d0) phi=phi+pi
       if (phi.gt.pi) phi=phi-2d0*pi
       phi2=atan(ps%p(2,this%next)/ps%p(1,this%next))
       if(ps%p(1,this%next).lt.0d0) phi2=phi2+pi
       if (phi2.gt.pi) phi2=phi2-2d0*pi
       phi=phi2-phi ! aximuthal separation particle 'next' and 'next-1'
       ! generate deltaR w.r.t. previously generated particle
       y=log((ps%p(0,this%next-1)+ps%p(3,this%next-1))/(ps%p(0,this%next-1)-ps%p(3,this%next-1)))/2d0
       drmin=max(this%DRcut(this%next*(this%next-1)+this%next-1),abs(phi))
       drmax=sqrt((this%ycut(this%next)+abs(y))**2+phi**2)
       ix=ix+1
       call random_to_var(ps%x(ix),0d0,drmin,drmax,dr,ps%jac)
       y=sqrt(dr**2-phi**2)
       if (next_extra_random().gt.0.5d0) y=-y
       ps%jac=ps%jac*2d0
       ps%p(3,this%next)=sqrt(pt2)*sinh(y)
       ps%p(0,this%next)=sqrt(pt2)*cosh(y)
       ps%jac=ps%jac*abs(dr/y)
       ! boost along the z-axis
       y=log((ps%p(0,this%next-1)+ps%p(3,this%next-1))/(ps%p(0,this%next-1)-ps%p(3,this%next-1)))/2d0
       call boostz(ps%p(0,this%next),-y,pb)
       ps%p(0:3,this%next)=pb(0:3)

    elseif(use_mode.eq.3) then
       ps%p(1,this%next)=-sum(ps%p(1,3:this%next-1))
       ps%p(2,this%next)=-sum(ps%p(2,3:this%next-1))
       pt2=ps%p(1,this%next)**2+ps%p(2,this%next)**2
       y=log((ps%p(0,this%next-1)+ps%p(3,this%next-1))/(ps%p(0,this%next-1)-ps%p(3,this%next-1)))/2d0
       call boostz(ps%p(0,this%next-1),y,pb)
       phi=delta_phi(ps%p(0,this%next),ps%p(0,this%next-1))
       invmmin=2d0*sqrt(pt2)*pb(0)*(1d0-cos(max(this%DRcut(this%next*(this%next-1)+this%next-1),phi)))
       invmmax=this%sqrts**2
       if(invmmin.ge.invmmax) then
          ps%jac=-1d0
          return
       endif
       ix=ix+1
       call random_to_var(ps%x(ix),-1d0,invmmin,invmmax,invm,ps%jac)
       ps%p(0,this%next)=(invm/2d0+pb(1)*ps%p(1,this%next)+pb(2)*ps%p(2,this%next))/pb(0)
       ! There are two values of the pz that correspond to a single
       ! invm. Take one of the two at random.
!!$       if (ps%p(0,this%next)**2.lt.pt2) then
!!$          ps%jac=-1d0
!!$          return
!!$       endif
       ps%p(3,this%next)=sqrt(ps%p(0,this%next)**2-pt2)
       if (next_extra_random().gt.0.5d0) ps%p(3,this%next)=-ps%p(3,this%next)
       ps%jac=ps%jac*2d0
       ps%jac=ps%jac/abs(2d0*pb(0)*ps%p(3,this%next))
       ! boost along the z-axis
       call boostz(ps%p(0,this%next),-y,pb)
       ps%p(0:3,this%next)=pb(0:3)

    elseif (use_mode.eq.4) then
       ps%p(1,this%next)=-sum(ps%p(1,3:this%next-1))
       ps%p(2,this%next)=-sum(ps%p(2,3:this%next-1))
       ps%p(0,this%next)=1d0 ! dummy value
       ps%p(3,this%next)=1d0 ! dummy value

       phi =atan(ps%p(2,this%next-1)/ps%p(1,this%next-1))
       if(ps%p(1,this%next-1).lt.0d0) phi=phi+pi
       phi2 =atan(ps%p(2,this%next)/ps%p(1,this%next))
       if(ps%p(1,this%next).lt.0d0) phi2=phi2+pi

       call rotz(ps%p(0,this%next),-phi,pb)
       ps%p(:,this%next)=pb

       phi =atan(ps%p(2,this%next-1)/ps%p(1,this%next-1))
       if(ps%p(1,this%next-1).lt.0d0) phi=phi+pi
       phi2 =atan(ps%p(2,this%next)/ps%p(1,this%next))
       if(ps%p(1,this%next).lt.0d0) phi2=phi2+pi

       pt2=ps%p(1,this%next)**2+ps%p(2,this%next)**2
       y=log((ps%p(0,this%next-1)+ps%p(3,this%next-1))/(ps%p(0,this%next-1)-ps%p(3,this%next-1)))/2d0
       ymax=this%ycut(this%next)+abs(y)
       pzmax = sqrt(pt2)/tan(2d0*atan(exp(-ymax)))

       costhetamin= abs(ps%p(1,this%next))/(dsqrt(pt2+pzmax**2))
       costhetamax= cos(max(abs(phi2),this%DRcut(this%next*(this%next-1)+this%next-1)))

       if (abs(phi2).ge.pi/2d0) then
          costhetamin= cos(max(abs(phi2),this%DRcut(this%next*(this%next-1)+this%next-1)))
          costhetamax = cos(pi-acos(abs(ps%p(1,this%next))/(dsqrt(pt2+pzmax**2))))
       endif
       if (costhetamin.gt.costhetamax) return

       ix=ix+1
       call random_to_var(ps%x(ix),0d0,costhetamin,costhetamax,costheta,ps%jac)

       ps%p(0,this%next) = abs(ps%p(1,this%next)/costheta)
       if (costheta.lt.0d0) then
          theta=acos(costheta)
          theta = pi-theta
          phi = acos(ps%p(2,this%next)/(ps%p(0,this%next)*dsqrt(1d0-cos(theta)**2)))
          ps%p(3,this%next) = dsqrt((ps%p(1,this%next)**2/costheta**2)-pt2)
       else
          theta=acos(costheta)
          phi = acos(ps%p(2,this%next)/(ps%p(0,this%next)*dsqrt(1d0-costheta**2)))
          ps%p(3,this%next) = dsqrt((ps%p(1,this%next)**2/costheta**2)-pt2)
       endif
       if (next_extra_random().gt.0.5d0) ps%p(3,this%next)=-ps%p(3,this%next)
       ps%jac=ps%jac*2d0
       ps%jac = ps%jac/(1+costheta**2-(costheta**2-1d0)*cos(2d0*phi))
       !ps%jac = ps%jac*ps%p(1,this%next)/sqrt(pt2)*0.5d0/(sin(2d0*atan(exp(-ymax)))**2)/((1d0+tan)**1.5d0)
       !ps%jac = ps%jac*cosh(y)
       call boostz(ps%p(0,this%next),-y,pb)
       phi=atan(ps%p(2,this%next-1)/ps%p(1,this%next-1))
       if(ps%p(1,this%next-1).lt.0d0) phi=phi+pi
       call  rotz(pb,phi,ps%p(0,this%next))

    endif

    ! initial states
    ptot(0:3)=sum(ps%p(0:3,3:this%next),dim=2)
    if (.not.all(ieee_is_finite(ptot))) then
       call reject_point(-20d0)
       return
    endif
    if (ptot(0).le.abs(ptot(3))) then
       call reject_point(-20d0)
       return
    endif
    ps%xbjrk(1)=(ptot(0)+ptot(3))/this%sqrts
    ps%xbjrk(2)=(ptot(0)-ptot(3))/this%sqrts
    if (.not.all(ieee_is_finite(ps%xbjrk))) then
       call reject_point(-1d0)
       return
    endif
    if (any(ps%xbjrk.le.0d0) .or. any(ps%xbjrk.ge.1d0)) then
       call reject_point(-1d0)
       return
    endif
    tau=ps%xbjrk(1)*ps%xbjrk(2)
    if (.not.ieee_is_finite(tau)) then
       call reject_point(-20d0)
       return
    endif
    if (tau.le.tiny(1d0)) then
       call reject_point(-20d0)
       return
    endif
    ps%p(0,1)=ps%xbjrk(1)*this%sqrts/2d0
    ps%p(1,1)=0d0
    ps%p(2,1)=0d0
    ps%p(3,1)=+ps%xbjrk(1)*this%sqrts/2d0
    ps%p(0,2)=ps%xbjrk(2)*this%sqrts/2d0
    ps%p(1,2)=0d0
    ps%p(2,2)=0d0
    ps%p(3,2)=-ps%xbjrk(2)*this%sqrts/2d0

    ! Ps%Jacobian factor (corresponds to the full ps%jacobian for
    ! use_mode=1. The other use_modes already have a partially computed
    ! ps%jacobian above)
    ps%jac=ps%jac/this%sqrts
    ps%jac=ps%jac/this%sqrts
    do i=3,this%next-1
       ps%jac=ps%jac/4d0
    enddo
    ! Add factors of 2*pi
    do i=1,3*(this%next-2)-4
       ps%jac=ps%jac/(2d0*pi)
    enddo
    ! Add flux factor
    ps%jac=ps%jac/(2d0*tau)
    ps%jac=ps%jac/this%sqrts
    ps%jac=ps%jac/this%sqrts
    if (.not.ieee_is_finite(ps%jac)) then
       call reject_point(-20d0)
       return
    endif
    if (ps%jac.le.0d0 .or. &
         .not.generated_momenta_are_valid(ps%p,this%masses,ps%xbjrk,.true.)) then
       call reject_point(-20d0)
       return
    endif
  contains
    subroutine reject_point(code)
      real(kind=8),intent(in) :: code
      if (allocated(ps%p)) ps%p=0d0
      ps%xbjrk=0d0
      ps%jac=code
    end subroutine reject_point

    real(kind=8) function next_extra_random()
      extra_index=extra_index+1
      if (extra_index.gt.this%ndim_extra) then
         next_extra_random=0.5d0
         ps%jac=-10d0
      else
         next_extra_random=ps%x(this%ndim+extra_index)
      endif
    end function next_extra_random

    subroutine fill_momentum_pt2cosphi(pt2,costheta,phi,p,jac)
      implicit none
      real(kind=8) :: pt2,costheta,phi,jac,pt
      real(kind=8),dimension(0:3) :: p

      pt=dsqrt(pt2)
      p(0) = pt/(dsqrt(costheta**2+(1d0-costheta**2)*cos(phi)**2))
      p(1) = p(0)*costheta
      p(2) = p(0)*dsqrt(1d0-costheta**2)*cos(phi)
      p(3) = p(0)*dsqrt(1d0-costheta**2)*sin(phi)
      jac=jac*2d0 ! left-over factor from the overall jacobian of pi's and 2's
      jac = jac/(1+costheta**2-(costheta**2-1d0)*cos(2d0*phi))
    end subroutine fill_momentum_pt2cosphi

    subroutine fill_momentum_pt2invmphi(pt2,invm,phi,Eref,p,jac)
      implicit none
      real(kind=8) :: pt2,invm,phi,jac,pt,Eref
      real(kind=8),dimension(0:3) :: p
      pt=sqrt(pt2)
      p(1)=pt*cos(phi)
      p(2)=pt*sin(phi)
      p(0)=invm/(2d0*Eref)+p(1)
      ! There are two values of the pz that correspond to a single
      ! invm. Take one of the two at random.
      p(3)=sqrt(p(0)**2-pt2)
      if (next_extra_random().gt.0.5d0) p(3)=-p(3)
      jac=jac*2d0
      jac=jac*abs(1d0/(2d0*Eref*p(3)))
    end subroutine fill_momentum_pt2invmphi
    subroutine fill_momentum_pt2drphi(pt2,dr,phi,p,jac)
      implicit none
      real(kind=8) :: pt2,dr,phi,jac,pt,y
      real(kind=8),dimension(0:3) :: p
      pt=sqrt(pt2)
      p(1)=pt*cos(phi)
      p(2)=pt*sin(phi)
      ! There are two values of the rapidity that correspond to a single
      ! DeltaR. Take one of the two at random.
      y=sqrt(dr**2-phi**2)
      if (next_extra_random().gt.0.5d0) y=-y
      jac=jac*2d0
      p(3)=pt*sinh(y)
      p(0)=pt*cosh(y)
      jac=jac*abs(dr/y)
    end subroutine fill_momentum_pt2drphi
    subroutine fill_momentum_pt2yphi(pt2,y,phi,p)
      implicit none
      real(kind=8) :: pt2,y,phi,pt
      real(kind=8),dimension(0:3) :: p
      pt=sqrt(pt2)
      p(1)=pt*cos(phi)
      p(2)=pt*sin(phi)
      p(3)=pt*sinh(y)
      p(0)=pt*cosh(y)
    end subroutine fill_momentum_pt2yphi

    subroutine random_to_var(x,power,var_min,var_max,var,jac)
      ! Given a random number x, it generates var in the range var_min
      ! <= var <= var_max according to var^(power)
      implicit none
      real(kind=8),intent(in) :: x,power,var_min,var_max
      real(kind=8),intent(out) :: var
      real(kind=8),intent(inout) :: jac
      real(kind=8) :: varmin,varmax,xloc,jac_factor,interval_tolerance
      real(kind=8) :: scale,updated_jac
      logical :: valid_power_map
      var=0d0
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(power) .or. &
           .not.ieee_is_finite(var_min) .or. .not.ieee_is_finite(var_max) .or. &
           .not.ieee_is_finite(jac)) then
         jac=-20d0
         return
      endif
      if (x.lt.0d0 .or. x.gt.1d0 .or. var_min.ge.var_max .or. jac.le.0d0) then
         jac=-20d0
         return
      endif
      if (var_min.lt.0d0 .and. var_max.le.0d0) then
         varmin=-var_max
         varmax=-var_min
      elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (power.ne.0d0)) then
         jac=-20d0
         return
      else
         varmin=var_min
         varmax=var_max
      endif
      scale=max(abs(varmin),abs(varmax))
      interval_tolerance=128d0*epsilon(1d0)*max(spacing(0d0),scale)
      if (varmax-varmin.le.interval_tolerance .or. &
           (power.lt.0d0 .and. varmin.le.0d0)) then
         jac=-20d0
         return
      endif
      xloc=max(0d0,min(1d0,x))
      if (power.eq.0d0) then
         var=varmin+xloc*(varmax-varmin)
         jac_factor=varmax-varmin
      elseif (power.eq.-1d0) then
         jac_factor=log(varmax)-log(varmin)
         if (.not.ieee_is_finite(jac_factor) .or. jac_factor.le.0d0) then
            jac=-20d0
            return
         endif
         var=exp((1d0-xloc)*log(varmin)+xloc*log(varmax))
         jac_factor=var*jac_factor
      else
         call stable_phase_space_power_map(xloc,power,varmin,varmax,var,&
              jac_factor,valid_power_map)
         if (.not.valid_power_map) then
            jac=-20d0
            return
         endif
      endif
      if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac_factor) .or. &
           jac_factor.le.0d0) then
         var=0d0
         jac=-20d0
         return
      endif
      if (.not.safe_phase_space_product(jac,jac_factor,updated_jac)) then
         var=0d0
         jac=-20d0
         return
      endif
      jac=updated_jac
      if (var_min.le.0d0 .and. var_max.le.0d0) then
         var=-var
      endif
      if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac)) then
         var=0d0
         jac=-20d0
         return
      endif
      if (jac.le.0d0) then
         var=0d0
         jac=-20d0
      endif
    end subroutine random_to_var

    real(kind=8) function dot(p1,p2)
      ! Inner product between two 4-vectors
      implicit none
      real(kind=8),intent(in),dimension(0:3) :: p1,p2
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
    end function dot

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

    real*8 function deltaR(p1,p2)
      ! Distance (Delta-R) between 'p1' and 'p2'
      implicit none
      real*8, dimension(0:3) :: p1,p2
      deltaR=sqrt(delta_phi(p1,p2)**2+(eta(p1)-eta(p2))**2)
    end function deltaR

    real*8 function eta(p)
      ! pseudo-rapidity of 'p'
      implicit none
      real*8, dimension(0:3) :: p
      real*8 :: theta
      theta=acos(p(3)/sqrt(p(1)**2+p(2)**2+p(3)**2))
      eta=-log(dtan(theta/2d0))
    end function eta

    real*8 function pt(p1)
      ! transverse momentum of 'p1'
      implicit none
      real*8, dimension(0:3) :: p1
      pt=sqrt(p1(1)**2+p1(2)**2)
    end function pt

    real*8 function delta_phi(p1,p2)
      ! azimuthal difference of 'p1' and 'p2'
      implicit none
      real*8, dimension(0:3) :: p1,p2
      real*8 :: denom
      denom=pt(p1)*pt(p2)
      delta_phi=acos((p1(1)*p2(1)+p1(2)*p2(2))/denom)
    end function delta_phi

  end subroutine genpt_generate_momenta
end module phase_space_genpt_mod
