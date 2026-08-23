module phase_space_onebody_mod
  use phase_space_base
  implicit none

  type,extends(phase_space_type),public :: phase_space_onebody
   contains
     procedure :: init => onebody_init
     procedure :: generate_momenta => onebody_generate_momenta
     procedure :: compute_x_from_momenta => onebody_compute_x_from_momenta
     procedure :: cleanup => onebody_cleanup
  end type phase_space_onebody

  private
  real(kind=8),parameter :: pi=3.1415926535897932384626433832795d0

contains

  subroutine onebody_init(this,sqrts,n,m,o,pt_cut,rap_cut,dr_cut,&
       sqrt_s_min,t_chan,include_pdf,flat)
    implicit none
    class(phase_space_onebody),intent(inout) :: this
    real(kind=8),intent(in) :: sqrts
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(n),intent(in) :: o
    real(kind=8),dimension(n),intent(in) :: m,pt_cut,rap_cut
    real(kind=8),dimension(n,n),intent(in) :: dr_cut,sqrt_s_min
    logical,intent(in) :: t_chan,include_pdf
    logical,intent(in),optional :: flat

    if (n.ne.3) then
       write (*,*) 'One-body phase space requires exactly two incoming and one final particle',n
       stop 1
    endif
    if (.not.include_pdf) then
       write (*,*) 'One-body cross sections require PDFs to integrate the partonic delta function'
       stop 1
    endif
    if (m(1).ne.0d0 .or. m(2).ne.0d0) then
       write (*,*) 'One-body phase space requires massless incoming partons',m(1:2)
       stop 1
    endif
    if (m(3).le.0d0 .or. m(3).ge.sqrts) then
       write (*,*) 'One-body phase space requires 0 < final mass < collider energy',m(3),sqrts
       stop 1
    endif

    this%next=n
    this%ndim=1
    this%ndim_extra=0
    this%sqrts=sqrts
    this%sqrtshat=m(3)
    this%s0=m(3)**2
    this%tot_mass=m(3)
    this%t_channel=t_chan
    this%can_invert_momenta=.true.
    allocate(this%masses(n))
    allocate(this%order(n))
    allocate(this%p(0:3,n))
    allocate(this%x(1))
    this%masses=m
    this%order=o
    this%p=0d0
    this%x=0d0
  end subroutine onebody_init

  subroutine onebody_generate_momenta(this,ps)
    implicit none
    class(phase_space_onebody),intent(inout) :: this
    type(psv),intent(inout) :: ps
    real(kind=8) :: mass,tau,ymax,y,root_tau

    mass=this%masses(3)
    tau=(mass/this%sqrts)**2
    root_tau=sqrt(tau)
    ymax=-0.5d0*log(tau)
    y=(2d0*ps%x(1)-1d0)*ymax
    ps%xbjrk(1)=root_tau*exp(y)
    ps%xbjrk(2)=root_tau*exp(-y)
    if (ps%xbjrk(1).le.0d0 .or. ps%xbjrk(1).ge.1d0 .or.&
         ps%xbjrk(2).le.0d0 .or. ps%xbjrk(2).ge.1d0) then
       ps%jac=-1d0
       return
    endif

    ps%p=0d0
    ps%p(0,1)=ps%xbjrk(1)*this%sqrts/2d0
    ps%p(3,1)=ps%p(0,1)
    ps%p(0,2)=ps%xbjrk(2)*this%sqrts/2d0
    ps%p(3,2)=-ps%p(0,2)
    ps%p(0:3,3)=ps%p(0:3,1)+ps%p(0:3,2)

    ! dPhi_1 = 2*pi*delta(shat-mass**2).  Integrating the delta
    ! function over tau=x1*x2 and mapping the remaining rapidity range
    ! to x in [0,1] gives this phase-space, flux and Bjorken Jacobian.
    ps%jac=(2d0*ymax)*pi/(mass**2*this%sqrts**2)
  end subroutine onebody_generate_momenta

  subroutine onebody_compute_x_from_momenta(this,ps)
    implicit none
    class(phase_space_onebody),intent(inout) :: this
    type(psv),intent(inout) :: ps
    real(kind=8),parameter :: tolerance=1d-10
    real(kind=8) :: mass,tau,ymax,y,x1,x2

    mass=this%masses(3)
    tau=(mass/this%sqrts)**2
    ymax=-0.5d0*log(tau)
    x1=2d0*ps%p(0,1)/this%sqrts
    x2=2d0*ps%p(0,2)/this%sqrts
    if (x1.le.0d0 .or. x1.ge.1d0 .or. x2.le.0d0 .or. x2.ge.1d0 .or.&
         abs(x1*x2-tau).gt.tolerance*tau) then
       ps%jac=-1d0
       return
    endif
    y=0.5d0*log(x1/x2)
    ps%x(1)=(y+ymax)/(2d0*ymax)
    if (ps%x(1).lt.-tolerance .or. ps%x(1).gt.1d0+tolerance) then
       ps%jac=-1d0
       return
    endif
    ps%x(1)=max(0d0,min(1d0,ps%x(1)))
    ps%xbjrk=[x1,x2]
    ps%jac=(2d0*ymax)*pi/(mass**2*this%sqrts**2)
  end subroutine onebody_compute_x_from_momenta

  subroutine onebody_cleanup(this)
    implicit none
    class(phase_space_onebody),intent(inout) :: this
    if (allocated(this%masses)) deallocate(this%masses)
    if (allocated(this%order)) deallocate(this%order)
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
  end subroutine onebody_cleanup

end module phase_space_onebody_mod
