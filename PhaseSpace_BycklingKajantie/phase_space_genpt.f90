module phase_space_genpt
  use common
  private
  integer(kind=4) :: ix
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical :: includePDF
  real(kind=8) :: sqrtshat,sqrts,tau,ycm,ptcut,ycut,DRcut
  integer,parameter :: use_mode=2 ! all modes use pT^2 and phi, but
                                  ! 1 = uses rapidity (original chili)
                                  ! 2 = uses DeltaR with previous particle
                                  ! 3 = uses invariant mass with previous particle
                                  ! 4 = uses cos(theta) with previous particle
  public :: genpt_init,genpt_phase_space
contains
  subroutine genpt_phase_space(xx)
    implicit none
    real(kind=8),dimension(99),intent(in) :: xx
    real(kind=8) :: pt2min,pt2max,pt2,ymin,ymax,y,phimin,phimax,phi&
         &,ycm,pt,drmin,drmax,dr,phi2,invmmin,invmmax,invm
    real(kind=8) :: costheta,costhetamin,costhetamax,theta,pzmax
    integer(kind=4) :: i,ix
    real(kind=8),dimension(0:3) :: ptot,pb
    real(kind=8),external :: ran2
    ix=0
    jac=1d0
    do i=3,next-1
       ! generate pT^2
       pt2min=ptcut**2
       pt2max=sqrts**2/4d0
       ix=ix+1
       call random_to_var(xx(ix),-1d0,pt2min,pt2max,pt2,jac)
       ! generate phi
       phimin=-pi
       phimax=pi
       ix=ix+1
       call random_to_var(xx(ix),0d0,phimin,phimax,phi,jac)
       if (use_mode.eq.1 .or.i.eq.3) then
          ! generate rapidity
          ymin=-ycut
          ymax=+ycut
          ix=ix+1
          call random_to_var(xx(ix),0d0,ymin,ymax,y,jac)
          ! fill momentum
          call fill_momentum_pt2yphi(pt2,y,phi,p(0,i))

       elseif (use_mode.eq.2) then
          ! generate deltaR w.r.t. previously generated particle
          y=log((p(0,i-1)+p(3,i-1))/(p(0,i-1)-p(3,i-1)))/2d0
          drmin=max(drcut,abs(phi))
          drmax=sqrt((ycut+abs(y))**2+phi**2)
          ix=ix+1
          call random_to_var(xx(ix),0d0,drmin,drmax,dr,jac)
          ! fill momentum, assuming that previous particle is along the x-axis.
          call fill_momentum_pt2drphi(pt2,dr,phi,p(0,i),jac)
          ! boost along the z-axis
          call boostz(p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(p(2,i-1)/p(1,i-1))
          if(p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,p(0,i))

       elseif (use_mode.eq.3) then
          ! get the energy in the frame where p(:,i-1) has p_z=0.
          y=log((p(0,i-1)+p(3,i-1))/(p(0,i-1)-p(3,i-1)))/2d0
          call boostz(p(0,i-1),y,pb)
          invmmin=2d0*sqrt(pt2)*pb(0)*(1d0-cos(max(drcut,abs(phi))))
          invmmax=sqrts**2
          ix=ix+1
          call random_to_var(xx(ix),-1d0,invmmin,invmmax,invm,jac)
          ! fill momentum, assuming that previous particle is along the x-axis.
          call fill_momentum_pt2invmphi(pt2,invm,phi,pb(0),p(0,i),jac)
          if (jac.lt.0d0) return
          ! boost along the z-axis
          call boostz(p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(p(2,i-1)/p(1,i-1))
          if(p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,p(0,i))

       elseif (use_mode.eq.4) then
          y=log((p(0,i-1)+p(3,i-1))/(p(0,i-1)-p(3,i-1)))/2d0
          costhetamin=-1d0 
          costhetamax=cos(drcut)
          ix=ix+1
          call random_to_var(xx(ix),0d0,costhetamin,costhetamax,costheta,jac)
          call fill_momentum_pt2cosphi(pt2,costheta,phi,p(0,i),jac)
          if (jac.lt.0d0) return
          call boostz(p(0,i),-y,pb)
          ! rotate about the z-axis
          phi=atan(p(2,i-1)/p(1,i-1))
          if(p(1,i-1).lt.0d0) phi=phi+pi
          call rotz(pb,phi,p(0,i))

       endif
    enddo


!!$    if (.true.) then
    if (use_mode.eq.1.or.use_mode.eq.4) then
! final particle: generate rapidity
       ymin=-ycut
       ymax=+ycut
       ix=ix+1
       call random_to_var(xx(ix),0d0,ymin,ymax,y,jac)
       ! final particle: fill momentum
       p(1,next)=-sum(p(1,3:next-1))
       p(2,next)=-sum(p(2,3:next-1))
       pt=sqrt(p(1,next)**2+p(2,next)**2)
       p(3,next)=pt*sinh(y)
       p(0,next)=pt*cosh(y)

    elseif(use_mode.eq.2) then
       p(1,next)=-sum(p(1,3:next-1))
       p(2,next)=-sum(p(2,3:next-1))
       pt=sqrt(p(1,next)**2+p(2,next)**2)
       phi=atan(p(2,next-1)/p(1,next-1))
       if(p(1,next-1).lt.0d0) phi=phi+pi
       if (phi.gt.pi) phi=phi-2d0*pi
       phi2=atan(p(2,next)/p(1,next))
       if(p(1,next).lt.0d0) phi2=phi2+pi
       if (phi2.gt.pi) phi2=phi2-2d0*pi
       phi=phi2-phi ! aximuthal separation particle 'next' and 'next-1'
       ! generate deltaR w.r.t. previously generated particle
       y=log((p(0,next-1)+p(3,next-1))/(p(0,next-1)-p(3,next-1)))/2d0
       drmin=max(drcut,abs(phi))
       drmax=sqrt((ycut+abs(y))**2+phi**2)
       ix=ix+1
       call random_to_var(xx(ix),0d0,drmin,drmax,dr,jac)
       y=sqrt(dr**2-phi**2)
       if (ran2().gt.0.5d0) y=-y
       jac=jac*2d0
       p(3,next)=pt*sinh(y)
       p(0,next)=pt*cosh(y)
       jac=jac*abs(dr/y)
       ! boost along the z-axis
       y=log((p(0,next-1)+p(3,next-1))/(p(0,next-1)-p(3,next-1)))/2d0
       call boostz(p(0,next),-y,pb)
       p(0:3,next)=pb(0:3)

    elseif(use_mode.eq.3) then
       p(1,next)=-sum(p(1,3:next-1))
       p(2,next)=-sum(p(2,3:next-1))
       pt=sqrt(p(1,next)**2+p(2,next)**2)
       y=log((p(0,next-1)+p(3,next-1))/(p(0,next-1)-p(3,next-1)))/2d0
       call boostz(p(0,next-1),y,pb)
       phi=delta_phi(p(0,next),p(0,next-1))
       invmmin=2d0*pt*pb(0)*(1d0-cos(max(drcut,phi)))
       invmmax=sqrts**2
       if(invmmin.ge.invmmax) then
          jac=-1d0
          return
       endif
       ix=ix+1
       call random_to_var(xx(ix),-1d0,invmmin,invmmax,invm,jac)
       p(0,next)=(invm/2d0+pb(1)*p(1,next)+pb(2)*p(2,next))/pb(0)
       ! There are two values of the pz that correspond to a single
       ! invm. Take one of the two at random.
!!$       if (p(0,next).lt.pt) then
!!$          jac=-1d0
!!$          return
!!$       endif
       p(3,next)=sqrt(p(0,next)**2-pt**2)
       if (ran2().gt.0.5d0) p(3,next)=-p(3,next)
       jac=jac*2d0
       jac=jac/abs(2d0*pb(0)*p(3,next))
       ! boost along the z-axis
       call boostz(p(0,next),-y,pb)
       p(0:3,next)=pb(0:3)

   elseif (use_mode.eq.4) then
       p(1,next)=-sum(p(1,3:next-1))
       p(2,next)=-sum(p(2,3:next-1))
       p(0,next)=1d0 ! dummy value
       p(3,next)=1d0 ! dummy value

       phi =atan(p(2,next-1)/p(1,next-1))
       if(p(1,next-1).lt.0d0) phi=phi+pi
       phi2 =atan(p(2,next)/p(1,next))
       if(p(1,next).lt.0d0) phi2=phi2+pi

       call rotz(p(0,next),-phi,pb)
       p(:,next)=pb

       phi =atan(p(2,next-1)/p(1,next-1))
       if(p(1,next-1).lt.0d0) phi=phi+pi
       phi2 =atan(p(2,next)/p(1,next))
       if(p(1,next).lt.0d0) phi2=phi2+pi

       pt=sqrt(p(1,next)**2+p(2,next)**2)
       y=log((p(0,next-1)+p(3,next-1))/(p(0,next-1)-p(3,next-1)))/2d0
       ymax=ycut+abs(y)
       pzmax = pt/tan(2d0*atan(exp(-ymax)))

       costhetamin= abs(p(1,next))/(dsqrt(pt**2+pzmax**2))
       costhetamax= cos(max(abs(phi2),drcut))
       
       if (abs(phi2).ge.pi/2d0) then
           costhetamin= cos(max(abs(phi2),drcut))
           costhetamax = cos(pi-acos(abs(p(1,next))/(dsqrt(pt**2+pzmax**2))))
       endif
       if (costhetamin.gt.costhetamax) return

       ix=ix+1
       call random_to_var(xx(ix),0d0,costhetamin,costhetamax,costheta,jac)
       
       p(0,next) = abs(p(1,next)/costheta)
       if (costheta.lt.0d0) then
           theta=acos(costheta)
           theta = pi-theta
           phi = acos(p(2,next)/(p(0,next)*dsqrt(1d0-cos(theta)**2)))
           p(3,next) = dsqrt((p(1,next)**2/costheta**2)-pt**2)
       else
           theta=acos(costheta)
           phi = acos(p(2,next)/(p(0,next)*dsqrt(1d0-costheta**2)))
           p(3,next) = dsqrt((p(1,next)**2/costheta**2)-pt**2)
       endif
       if (ran2().gt.0.5d0) p(3,next)=-p(3,next)
       jac=jac*2d0
       jac = jac/(1+costheta**2-(costheta**2-1d0)*cos(2d0*phi))
       !jac = jac*p(1,next)/pt*0.5d0/(sin(2d0*atan(exp(-ymax)))**2)/((1d0+tan)**1.5d0)
       !jac = jac*cosh(y)
       call boostz(p(0,next),-y,pb)
       phi=atan(p(2,next-1)/p(1,next-1))
       if(p(1,next-1).lt.0d0) phi=phi+pi
       call  rotz(pb,phi,p(0,next))

    endif

    
! initial states
    ptot(0:3)=sum(p(0:3,3:next),dim=2)
    tau=dot(ptot,ptot)/sqrts**2
    ycm=log((ptot(0)+ptot(3))/(ptot(0)-ptot(3)))/2d0
    xbjrk(1)=sqrt(tau)*exp(ycm)
    xbjrk(2)=sqrt(tau)*exp(-ycm)
    if (xbjrk(1).ge.1d0 .or. xbjrk(2).ge.1d0) then
       jac=-1d0
       return
    endif
    p(0,1)=xbjrk(1)*sqrts/2d0
    p(1,1)=0d0
    p(2,1)=0d0
    p(3,1)=+xbjrk(1)*sqrts/2d0
    p(0,2)=xbjrk(2)*sqrts/2d0
    p(1,2)=0d0
    p(2,2)=0d0
    p(3,2)=-xbjrk(2)*sqrts/2d0

! Jacobian factor (corresponds to the full jacobian for
! use_mode=1. The other use_modes already have a partially computed
! jacobian above)
    jac=jac/(sqrts**2*dble(4**(next-3)))
! Add factors of 2*pi
    jac=jac/((2d0*pi)**(3*(next-2)-4))
! Add flux factor
    jac=jac/(2d0*tau*sqrts**2)
  end subroutine genpt_phase_space

  subroutine fill_momentum_pt2cosphi(pt2,costheta,phi,p,jac)
   implicit none
    real(kind=8) :: pt2,costheta,phi,jac,pt
    real(kind=8),dimension(0:3) :: p
    real(kind=8),external :: ran2
  
    pt=dsqrt(pt2)
    p(0) = pt/(dsqrt(costheta**2+(1d0-costheta**2)*cos(phi)**2))
    p(1) = p(0)*costheta
    p(2)= p(0)*dsqrt(1d0-costheta**2)*cos(phi)
    p(3) = p(0)*dsqrt(1d0-costheta**2)*sin(phi)
    jac=jac*2d0 ! left-over factor from the overall jacobian of pi's and 2's
    jac = jac/(1+costheta**2-(costheta**2-1d0)*cos(2d0*phi))
  end subroutine fill_momentum_pt2cosphi

  subroutine fill_momentum_pt2invmphi(pt2,invm,phi,Eref,p,jac)
    implicit none
    real(kind=8) :: pt2,invm,phi,jac,pt,Eref
    real(kind=8),dimension(0:3) :: p
    real(kind=8),external :: ran2
    pt=sqrt(pt2)
    p(1)=pt*cos(phi)
    p(2)=pt*sin(phi)
    p(0)=invm/(2d0*Eref)+p(1)
    ! There are two values of the pz that correspond to a single
    ! invm. Take one of the two at random.
    p(3)=sqrt(p(0)**2-pt2)
    if (ran2().gt.0.5d0) p(3)=-p(3)
    jac=jac*2d0
    jac=jac*abs(1d0/(2d0*Eref*p(3)))
  end subroutine fill_momentum_pt2invmphi
  subroutine fill_momentum_pt2drphi(pt2,dr,phi,p,jac)
    implicit none
    real(kind=8) :: pt2,dr,phi,jac,pt,y
    real(kind=8),dimension(0:3) :: p
    real(kind=8),external :: ran2
    pt=sqrt(pt2)
    p(1)=pt*cos(phi)
    p(2)=pt*sin(phi)
    ! There are two values of the rapidity that correspond to a single
    ! DeltaR. Take one of the two at random.
    y=sqrt(dr**2-phi**2)
    if (ran2().gt.0.5d0) y=-y
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
  subroutine genpt_init(sqrtsh,n,m,ptmin,rapcut,DRmin,include_pdf)
    implicit none
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrtsh
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! rapidity and pT cut on all final state particles
    real(kind=8),intent(in) :: rapcut,ptmin,DRmin
    ! masses of all the particles. The two incoming particles must be
    ! massless.
    real(kind=8),dimension(n),intent(in) :: m
    logical,intent(in) :: include_pdf
    ptcut=ptmin
    ycut=rapcut
    DRcut=DRmin
    sqrts=sqrtsh
    next=n
    if (.not.include_pdf) then
       write (*,*) 'genpt phase-space only with include_pdf=.true.'
       stop 1
    endif
    if (any(m(1:n).ne.0d0)) then
       write (*,*) 'genpt phase-space only for all massless particles'
       stop 1
    endif
    if (use_mode.eq.1 .and. rapcut.le.0d0) then
       write (*,*) 'genpt phase-space must have rapidity cut'
       stop 1
    endif
    if (use_mode.ne.1 .and. DRcut.le.0d0) then
       write (*,*) 'genpt phase-space must have DeltaR cut'
       stop 1
    endif
    if (ptcut.le.0d0) then
       write (*,*) 'genpt phase-space must have pT cut'
       stop 1
    endif
    allocate(p(0:3,next))
  end subroutine genpt_init

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

  real*8 function delta_phi(p1,p2)
    ! azimuthal difference of 'p1' and 'p2'
    implicit none
    real*8, dimension(0:3) :: p1,p2
    real*8 :: denom
    denom=pt(p1)*pt(p2)
    delta_phi=acos((p1(1)*p2(1)+p1(2)*p2(2))/denom)
  end function delta_phi

  real*8 function pt(p)
    ! transverse momentum of 'p'
    implicit none
    real*8, dimension(0:3) :: p
    pt=sqrt(p(1)**2+p(2)**2)
  end function pt
end module phase_space_genpt
