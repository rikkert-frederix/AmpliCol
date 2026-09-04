module phase_space_base
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  ! GEN23 forms a product of two cubic invariant polynomials.  Keep every
  ! momentum far enough below huge() that this twelfth-power scale remains
  ! representable even before the final output validator can run.
  real(kind=8),parameter,public :: phase_space_momentum_limit= &
       0.125d0*huge(1d0)**(1d0/12d0)
  ! GEN23 and HAAG store every subset in default-integer bit masks.  This
  ! ceiling keeps MASKR/IBSET defined and caps each subset workspace at about
  ! one million entries instead of allowing corrupt metadata to request an
  ! unbounded exponential allocation.
  integer,parameter,public :: max_bitmask_particles=20
  private :: stable_phase_space_rotation
  type :: psv
     real(kind=8),dimension(:,:),allocatable,public :: p
     real(kind=8),dimension(:),allocatable,public :: x
     real(kind=8),public :: jac=0d0,xbjrk(2)=0d0
   contains
     final :: finalize_psv
  end type psv
  type,abstract :: phase_space_type
     ! If adding variables here, make sure to also update the
     ! 'cleanup' subroutines for all phase-space parametrisation
     ! module.
     real(kind=8),dimension(:,:),allocatable,public :: p
     real(kind=8),public :: jac=0d0,xbjrk(2)=0d0
     real(kind=8) :: s0=0d0,tot_mass=0d0,sqrtshat=0d0,sqrts=0d0
     real(kind=8),dimension(:),allocatable :: masses,ptcut,ycut,drcut
     integer(kind=4) :: ndim=0,next=0,ndim_extra=0
     integer(kind=4),dimension(:,:),allocatable :: sets
     real(kind=8),dimension(:,:),allocatable :: pp,sqrt_s_min
     real(kind=8),dimension(:),allocatable :: x,invm
     real(kind=8),dimension(:,:),allocatable :: invm_min,invm_max,ETmin
     integer(kind=4),dimension(:),allocatable :: order
     logical :: t_channel=.false.
     logical :: can_invert_momenta=.false.
   contains
     procedure(phase_space_interface_init),deferred :: init
     procedure(phase_space_interface_generate_momenta),deferred :: generate_momenta
     procedure(phase_space_interface_compute_x_from_momenta),deferred :: compute_x_from_momenta
     procedure(phase_space_interface_cleanup),deferred :: cleanup
  end type phase_space_type
  
  ! Declare the abstract interface for the procedures
  abstract interface
     subroutine phase_space_interface_init(this,sqrts,n,m,o,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,flat)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
       real(kind=8),intent(in) :: sqrts
       integer(kind=4),intent(in) :: n
       integer(kind=4),dimension(n),intent(in) :: o
       real(kind=8),dimension(n,n),intent(in) :: dr_cut,sqrt_s_min
       real(kind=8),dimension(n),intent(in) :: m,pt_cut,rap_cut
       logical,intent(in) :: t_chan
       logical,intent(in) :: include_pdf
       logical,intent(in),optional :: flat
     end subroutine phase_space_interface_init
     subroutine phase_space_interface_generate_momenta(this,ps)
       import :: phase_space_type,psv
       class(phase_space_type),intent(inout) :: this
       type(psv),intent(inout) :: ps
     end subroutine phase_space_interface_generate_momenta
     subroutine phase_space_interface_compute_x_from_momenta(this,ps)
       import :: phase_space_type,psv
       class(phase_space_type),intent(inout) :: this
       type(psv),intent(inout) :: ps
     end subroutine phase_space_interface_compute_x_from_momenta
     subroutine phase_space_interface_cleanup(this)
       import :: phase_space_type
       class(phase_space_type),intent(inout) :: this
     end subroutine phase_space_interface_cleanup
  end interface
contains
  subroutine stable_rotate_from_z_axis(p,q,rotated,valid)
    real(kind=8),intent(in) :: p(0:3),q(0:3)
    real(kind=8),intent(out) :: rotated(0:3)
    logical,intent(out) :: valid
    call stable_phase_space_rotation(p,q,.false.,rotated,valid)
  end subroutine stable_rotate_from_z_axis

  subroutine stable_rotate_to_z_axis(p,q,rotated,valid)
    real(kind=8),intent(in) :: p(0:3),q(0:3)
    real(kind=8),intent(out) :: rotated(0:3)
    logical,intent(out) :: valid
    call stable_phase_space_rotation(p,q,.true.,rotated,valid)
  end subroutine stable_rotate_to_z_axis

  subroutine stable_phase_space_rotation(p,q,inverse,rotated,valid)
    real(kind=8),intent(in) :: p(0:3),q(0:3)
    logical,intent(in) :: inverse
    real(kind=8),intent(out) :: rotated(0:3)
    logical,intent(out) :: valid
    real(kind=8) :: spatial_scale,qx,qy,qz,qt,qnorm
    real(kind=8) :: cos_phi,sin_phi,cos_theta,sin_theta,direction_sign

    rotated=0d0
    valid=.false.
    if (.not.all(ieee_is_finite(p)) .or. .not.all(ieee_is_finite(q))) return
    rotated(0)=p(0)
    spatial_scale=maxval(abs(q(1:3)))
    if (spatial_scale.eq.0d0) then
       rotated(1:3)=p(1:3)
       valid=.true.
       return
    endif
    qx=q(1)/spatial_scale
    qy=q(2)/spatial_scale
    qz=q(3)/spatial_scale
    qt=sqrt(qx*qx+qy*qy)
    if (qt.eq.0d0) then
       direction_sign=sign(1d0,qz)
       rotated(1:3)=direction_sign*p(1:3)
       valid=.true.
       return
    endif
    qnorm=sqrt(qt*qt+qz*qz)
    if (.not.ieee_is_finite(qnorm) .or. qnorm.le.0d0) return
    cos_phi=qx/qt
    sin_phi=qy/qt
    cos_theta=qz/qnorm
    sin_theta=qt/qnorm
    if (inverse) then
       rotated(1)=cos_phi*cos_theta*p(1)+sin_phi*cos_theta*p(2)-sin_theta*p(3)
       rotated(2)=-sin_phi*p(1)+cos_phi*p(2)
       rotated(3)=cos_phi*sin_theta*p(1)+sin_phi*sin_theta*p(2)+cos_theta*p(3)
    else
       rotated(1)=cos_phi*cos_theta*p(1)-sin_phi*p(2)+cos_phi*sin_theta*p(3)
       rotated(2)=sin_phi*cos_theta*p(1)+cos_phi*p(2)+sin_phi*sin_theta*p(3)
       rotated(3)=-sin_theta*p(1)+cos_theta*p(3)
    endif
    valid=all(ieee_is_finite(rotated))
    if (.not.valid) rotated=0d0
  end subroutine stable_phase_space_rotation

  logical function safe_phase_space_product(first,second,value) result(valid)
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    real(kind=8) :: abs_first,abs_second

    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       valid=.true.
       return
    endif
    abs_first=abs(first)
    abs_second=abs(second)
    if (abs_second.gt.1d0) then
       if (abs_first.gt.huge(1d0)/abs_second) return
    elseif (abs_first.gt.1d0) then
       if (abs_second.gt.huge(1d0)/abs_first) return
    endif
    value=first*second
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function safe_phase_space_product

  logical function safe_phase_space_ratio(numerator,denominator,value) result(valid)
    real(kind=8),intent(in) :: numerator,denominator
    real(kind=8),intent(out) :: value

    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator)) return
    if (denominator.eq.0d0) return
    if (abs(denominator).lt.1d0) then
       if (abs(numerator).gt.huge(1d0)*abs(denominator)) return
    endif
    value=numerator/denominator
    valid=ieee_is_finite(value)
    if (.not.valid) value=0d0
  end function safe_phase_space_ratio

  logical function safe_phase_space_scaled_ratio(first,second,denominator,value) result(valid)
    ! Evaluate first*second/denominator without committing to an order that
    ! can overflow or underflow even though the final value is representable.
    real(kind=8),intent(in) :: first,second,denominator
    real(kind=8),intent(out) :: value
    real(kind=8) :: intermediate

    valid=.false.
    value=0d0
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second) .or. &
         .not.ieee_is_finite(denominator)) return
    if (denominator.eq.0d0) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       valid=.true.
       return
    endif

    if (safe_phase_space_ratio(second,denominator,intermediate)) then
       if (intermediate.ne.0d0) then
          if (safe_phase_space_product(first,intermediate,value)) then
             if (value.ne.0d0) then
                valid=.true.
                return
             endif
          endif
       endif
    endif
    if (safe_phase_space_product(first,second,intermediate)) then
       if (intermediate.ne.0d0) then
          if (safe_phase_space_ratio(intermediate,denominator,value)) then
             if (value.ne.0d0) then
                valid=.true.
                return
             endif
          endif
       endif
    endif
    value=0d0
  end function safe_phase_space_scaled_ratio

  subroutine stable_phase_space_power_map(x,power,varmin,varmax,var,jac_factor,valid)
    ! Invert a v**power cumulative distribution without ever materialising
    ! dimensionful powers.  The scaled logarithmic form keeps both the
    ! mapped variable and its derivative representable near the ends of the
    ! double-precision range.
    real(kind=8),intent(in) :: x,power,varmin,varmax
    real(kind=8),intent(out) :: var,jac_factor
    logical,intent(out) :: valid
    real(kind=8) :: exponent,scale,log_ratio,ratio_power,base
    real(kind=8) :: log_scale,log_var,log_jac,delta

    valid=.false.
    var=0d0
    jac_factor=0d0
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(power) .or. &
         .not.ieee_is_finite(varmin) .or. .not.ieee_is_finite(varmax)) return
    if (x.lt.0d0 .or. x.gt.1d0 .or. varmin.lt.0d0 .or. &
         varmax.le.varmin) return
    exponent=1d0+power
    if (abs(exponent).le.epsilon(1d0)) return
    if (exponent.lt.0d0 .and. varmin.le.0d0) return
    if (exponent.gt.0d0) then
       scale=varmax
       if (varmin.eq.0d0) then
          ratio_power=0d0
       else
          log_ratio=log(varmin)-log(varmax)
          ratio_power=exp(exponent*log_ratio)
       endif
       base=(1d0-x)*ratio_power+x
    else
       scale=varmin
       log_ratio=log(varmax)-log(varmin)
       ratio_power=exp(exponent*log_ratio)
       base=(1d0-x)+x*ratio_power
    endif
    if (.not.ieee_is_finite(ratio_power)) return
    delta=1d0-ratio_power
    if (delta.le.0d0) return
    log_scale=log(scale)
    if (x.eq.0d0) then
       if (varmin.le.0d0) return
       log_var=log(varmin)
    elseif (x.eq.1d0) then
       log_var=log(varmax)
    else
       if (.not.ieee_is_finite(base) .or. base.le.0d0) return
       log_var=log_scale+log(base)/exponent
    endif
    if (.not.ieee_is_finite(log_var)) return
    if (log_var.gt.log(huge(1d0)) .or. &
         log_var.lt.log(spacing(0d0))) return
    var=exp(log_var)
    if (.not.ieee_is_finite(var) .or. var.le.0d0) then
       var=0d0
       return
    endif
    log_jac=log_scale+log(delta)-log(abs(exponent))-&
         power*(log_var-log_scale)
    if (.not.ieee_is_finite(log_jac)) then
       var=0d0
       return
    endif
    if (log_jac.gt.log(huge(1d0)) .or. &
         log_jac.lt.log(spacing(0d0))) then
       var=0d0
       return
    endif
    jac_factor=exp(log_jac)
    if (.not.ieee_is_finite(jac_factor) .or. jac_factor.le.0d0) then
       var=0d0
       jac_factor=0d0
       return
    endif
    valid=.true.
  end subroutine stable_phase_space_power_map

  logical function generated_momenta_are_valid(p,masses,xbjrk,include_pdf) result(valid)
    ! Validate a generated physical point before it reaches cuts or matrix
    ! elements.  The tolerance follows the arithmetic scale of each
    ! cancellation, so it is independent of the chosen momentum unit.
    real(kind=8),intent(in) :: p(0:,:),masses(:),xbjrk(2)
    logical,intent(in) :: include_pdf
    real(kind=8),parameter :: safety=4096d0
    real(kind=8) :: shell,shell_scale,tolerance
    real(kind=8) :: incoming(0:3),outgoing(0:3),residual(0:3),component_scale
    integer :: i,mu,n

    valid=.false.
    n=size(p,2)
    if (size(p,1).ne.4 .or. size(masses).ne.n .or. n.lt.3) return
    if (.not.all(ieee_is_finite(p)) .or. .not.all(ieee_is_finite(masses))) return
    if (max(maxval(abs(p)),maxval(abs(masses))).gt.phase_space_momentum_limit) return
    if (any(masses.lt.0d0)) return
    if (any(p(0,1:n).le.0d0)) return
    if (include_pdf) then
       if (.not.all(ieee_is_finite(xbjrk))) return
       tolerance=safety*epsilon(1d0)
       if (any(xbjrk.le.0d0) .or. any(xbjrk.gt.1d0+tolerance)) return
    endif

    do i=1,n
       shell=p(0,i)*p(0,i)-p(1,i)*p(1,i)-p(2,i)*p(2,i)-p(3,i)*p(3,i)
       shell_scale=abs(p(0,i)*p(0,i))+abs(p(1,i)*p(1,i))+&
            abs(p(2,i)*p(2,i))+abs(p(3,i)*p(3,i))+masses(i)*masses(i)
       if (.not.ieee_is_finite(shell) .or. .not.ieee_is_finite(shell_scale)) return
       tolerance=safety*epsilon(1d0)*max(tiny(1d0),shell_scale)
       if (abs(shell-masses(i)*masses(i)).gt.tolerance) return
    enddo

    incoming=p(:,1)+p(:,2)
    outgoing=0d0
    do i=3,n
       outgoing=outgoing+p(:,i)
    enddo
    residual=incoming-outgoing
    do mu=0,3
       component_scale=abs(p(mu,1))+abs(p(mu,2))
       do i=3,n
          component_scale=component_scale+abs(p(mu,i))
       enddo
       tolerance=safety*epsilon(1d0)*max(tiny(1d0),component_scale)
       if (abs(residual(mu)).gt.tolerance) return
    enddo
    valid=.true.
  end function generated_momenta_are_valid

  subroutine finalize_psv(this)
    type(psv),intent(inout) :: this
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%x)) deallocate(this%x)
  end subroutine finalize_psv
end module phase_space_base
