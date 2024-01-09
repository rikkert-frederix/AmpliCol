module phase_space_genpt
  use common
  private
  integer(kind=4) :: ix
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical :: t_channel,includePDF
  real(kind=8) :: sqrtshat,sqrts,tau,ycm,ptcut,ycut
  integer :: nquarks

  public :: genpt_init,genpt_phase_space
contains
  subroutine genpt_phase_space(xx)
    implicit none
    real(kind=8),dimension(99),intent(in) :: xx
    real(kind=8) :: pt2min,pt2max,pt2,ymin,ymax,y,phimin,phimax,phi,ycm,pt
    integer(kind=4) :: i,ix
    real(kind=8),dimension(0:3) :: ptot
    ix=0
    jac=1d0
    do i=3,next-1
       ! generate pT^2
       pt2min=ptcut**2
       pt2max=sqrts**2/4d0
       ix=ix+1
       call random_to_var(xx(ix),-1d0,pt2min,pt2max,pt2,jac)
       ! generate rapidity
       ymin=-ycut
       ymax=+ycut
       ix=ix+1
       call random_to_var(xx(ix),0d0,ymin,ymax,y,jac)
       ! generate phi
       phimin=0d0
       phimax=2d0*pi
       ix=ix+1
       call random_to_var(xx(ix),0d0,phimin,phimax,phi,jac)
       ! fill momentum
       call fill_momentum_pt2yphi(pt2,y,phi,p(0,i))
    enddo
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
! Jacobian factor
    jac=jac/(sqrts**2*dble(4**(next-3)))
! Add factors of 2*pi
    jac=jac/((2d0*pi)**(3*(next-2)-4))
! Add flux factor
    jac=jac/(2d0*tau*sqrts**2)
  end subroutine genpt_phase_space
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
  subroutine genpt_init(sqrtsh,n,m,ptmin,rapcut,include_pdf)
    implicit none
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrtsh
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! rapidity and pT cut on all final state particles
    real(kind=8),intent(in) :: rapcut,ptmin
    ! masses of all the particles. The two incoming particles must be
    ! massless.
    real(kind=8),dimension(n),intent(in) :: m
    logical,intent(in) :: include_pdf
    ptcut=ptmin
    ycut=rapcut
    sqrts=sqrtsh
    if (.not.include_pdf) then
       write (*,*) 'genpt phase-space only with include_pdf=.true.'
       stop 1
    endif
    if (any(m(1:n).ne.0d0)) then
       write (*,*) 'genpt phase-space only for all massless particles'
       stop 1
    endif
    if (rapcut.le.0d0) then
       write (*,*) 'genpt phase-space must have rapidity cut'
       stop 1
    endif
    if (ptcut.le.0d0) then
       write (*,*) 'genpt phase-space must have pT cut'
       stop 1
    endif
    next=n
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
end module phase_space_genpt
