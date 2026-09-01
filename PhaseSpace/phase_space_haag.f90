module phase_space_haag_mod
  !  use common
  use phase_space_base
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type,extends(phase_space_type),public :: phase_space_haag
     logical :: flat_mode=.false.
     logical :: include_pdf=.false.
   contains
     procedure :: init => haag_init
     procedure :: generate_momenta => haag_generate_momenta
     procedure :: compute_x_from_momenta => haag_compute_x_from_momenta
     procedure :: cleanup => haag_cleanup
  end type phase_space_haag
  private
  real(kind=8),parameter :: pi=3.1415926535897932d0

  logical,parameter :: open=.false.
  logical,parameter :: flat_split=.false.,a1_split=.false.

  logical,parameter :: verbose=.true.
  logical,parameter :: exper=.false.
  real(kind=8),parameter :: ip=-1d0,ip_shat=-2d0
  real(kind=8),parameter :: vtiny=1d-12

contains
  logical function clamp_haag_interval(value,lower,upper,relative_tolerance) result(ok)
    ! Analytic inversions occasionally overshoot a closed interval by a few
    ! ulps. Clamp only such roundoff-sized excursions; reject a genuinely
    ! out-of-range result without dividing by a possibly zero endpoint.
    implicit none
    real(kind=8),intent(inout) :: value
    real(kind=8),intent(in) :: lower,upper,relative_tolerance
    real(kind=8) :: scale

    ok=.false.
    if (.not.ieee_is_finite(value) .or. .not.ieee_is_finite(lower) .or. &
         .not.ieee_is_finite(upper) .or. .not.ieee_is_finite(relative_tolerance)) return
    if (relative_tolerance.lt.0d0) return
    if (lower.gt.upper) return
    scale=max(1d0,abs(lower),abs(upper))
    if (value.lt.lower) then
       if (lower-value.gt.relative_tolerance*scale) return
       value=lower
    elseif (value.gt.upper) then
       if (value-upper.gt.relative_tolerance*scale) return
       value=upper
    endif
    ok=.true.
  end function clamp_haag_interval

  subroutine haag_compute_x_from_momenta(this,ps)
    implicit none
    class(phase_space_haag),intent(inout) :: this
    type(psv),intent(inout) :: ps
    if (allocated(ps%x)) ps%x=0d0
    ps%jac=-14d0
  end subroutine haag_compute_x_from_momenta
  subroutine haag_cleanup(this)
    implicit none
    class(phase_space_haag),intent(inout) :: this
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
    this%include_pdf=.false.
    this%flat_mode=.false.
  end subroutine haag_cleanup
  subroutine haag_init(this,sqrts,n,m,o,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,flat)
    implicit none
    class(phase_space_haag),intent(inout) :: this
    real(kind=8),intent(in) :: sqrts
    integer(kind=4),intent(in) :: n
    integer(kind=4),dimension(n),intent(in) :: o
    ! rapidity and pT cut (and DR and sqrt_s_min) on all the particles
    real(kind=8),dimension(n),intent(in) :: rap_cut,pt_cut
    real(kind=8),dimension(n,n),intent(in) :: DR_cut,sqrt_s_min
    ! masses
    real(kind=8),dimension(n),intent(in) :: m
    real(kind=8),dimension(2) :: s_cut
    logical,intent(in) :: t_chan
    logical,intent(in) :: include_pdf
    logical,intent(in),optional :: flat
    integer(kind=4) :: i,j,allocation_status
    real(kind=8) :: drjj_min,pt_min,sqrt_smin
    character(len=256) :: allocation_message
    integer(kind=4),dimension(n) :: order
    real(kind=8),dimension(:),allocatable :: invm,invm_min,invm_max,sigma_ij
    if (n.lt.4) then
       write (*,*) 'ERROR in haag_init() -- HAAG requires at least two final-state particles'
       stop 1
    endif
    if (n.gt.max_bitmask_particles) then
       write (*,*) 'ERROR in haag_init() -- particle multiplicity exceeds bit-mask workspace limit:',&
            n,max_bitmask_particles
       stop 1
    endif
    if (.not.ieee_is_finite(sqrts)) then
       write (*,*) 'ERROR in haag_init() -- invalid collider energy:',sqrts
       stop 1
    endif
    if (sqrts.le.0d0) then
       write (*,*) 'ERROR in haag_init() -- invalid collider energy:',sqrts
       stop 1
    endif
    if (sqrts.gt.phase_space_momentum_limit) then
       write (*,*) 'ERROR in haag_init() -- collider energy exceeds numerical range:',sqrts
       stop 1
    endif
    do i=1,n
       if (count(o.eq.i).ne.1) then
          write (*,*) 'ERROR in haag_init() -- colour order is not a permutation:',o
          stop 1
       endif
       if (.not.ieee_is_finite(m(i))) then
          write (*,*) 'ERROR in haag_init() -- invalid external mass:',i,m(i)
          stop 1
       endif
       if (m(i).lt.0d0) then
          write (*,*) 'ERROR in haag_init() -- invalid external mass:',i,m(i)
          stop 1
       endif
    enddo
    if (any(.not.ieee_is_finite(pt_cut)) .or. any(.not.ieee_is_finite(rap_cut)) .or. &
         any(.not.ieee_is_finite(dr_cut)) .or. any(.not.ieee_is_finite(sqrt_s_min))) then
       write (*,*) 'ERROR in haag_init() -- non-finite phase-space cut'
       stop 1
    endif
    if (max(maxval(abs(m)),maxval(abs(pt_cut)),maxval(abs(sqrt_s_min))).gt. &
         phase_space_momentum_limit) then
       write (*,*) 'ERROR in haag_init() -- mass or dimensionful cut exceeds numerical range'
       stop 1
    endif
    this%flat_mode=.false.
    if (present(flat)) this%flat_mode=flat
    this%sqrtshat=sqrts
    this%sqrts=sqrts
    drjj_min=99d99
    do i=3,n-1
       drjj_min=min(drjj_min,minval(dr_cut(i,i+1:n)))
    enddo
    pt_min=minval(pt_cut(3:n))
    sqrt_smin=99d99
    do i=3,n-1
       sqrt_smin=min(sqrt_smin,minval(sqrt_s_min(i,i+1:n)))
    enddo
    if (verbose) then
       write (*,*) 'Setting up',n,'particle phase-space'
       write (*,*) 'Total available energy, sqrt(s-hat) =',this%sqrtshat
    endif
    this%include_pdf=include_pdf
    this%t_channel=t_chan
    this%can_invert_momenta=.false.
    call haag_deallocate
    this%next=n
    this%ndim=3*(this%next-2)-4
    if (this%include_pdf) this%ndim=this%ndim+2 ! the two Bjorken x's
    ! HAAG uses discrete multichannel, orientation, and quadratic-root
    ! choices in addition to its adapted coordinates.  Supply those choices
    ! as flat coordinates so the map is deterministic for a fixed point.
    ! Four per final-state particle plus two is a conservative upper bound;
    ! unused flat coordinates intentionally carry unit weight.
    this%ndim_extra=4*(this%next-2)+2
    allocate(invm(maskr(this%next)),invm_min(maskr(this%next)),&
         invm_max(maskr(this%next)),sigma_ij(maskr(this%next)),&
         this%pp(0:3,0:maskr(this%next)),this%masses(this%next),&
         this%p(0:3,this%next),this%x(this%ndim+this%ndim_extra),&
         this%sets(0:this%next-2,2),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       call haag_deallocate
       write (*,*) 'ERROR in haag_init() -- unable to allocate phase-space workspace:',&
            trim(allocation_message)
       stop 1
    endif
    invm=0d0
    invm_min=0d0
    invm_max=0d0
    sigma_ij=0d0
    this%pp(0:3,0:maskr(this%next))=0d0
    this%p=0d0
    this%x=0d0
    this%sets=0
    ! masses of external particles
    do i=1,n
       if ((i.eq.1 .or. i.eq.2) .and. m(i).ne.0d0) then
          write (*,*) 'ERROR in haag_init() -- ', &
               & 'incoming particles should be massless'
          write (*,*) m
          stop 1
       endif
       invm(ibset(0,i-1))=m(i)**2
    enddo

    s_cut(1)=0d0 ! cut on invariant between initial and final state particle
    s_cut(2)=0d0 ! cut on invariant of two final state particles.
    if (sqrt_smin.gt.0d0) then
       s_cut(1)=max(s_cut(1),sqrt_smin**2)
       s_cut(2)=max(s_cut(2),sqrt_smin**2)
    endif
    if (pt_min.gt.0d0) then
       s_cut(1)=max(s_cut(1),pt_min**2)
       if (drjj_min.gt.0d0) &
            s_cut(2)=max(s_cut(2),2d0*pt_min**2*(1d0-cos(DRjj_min)))
    endif
    call setup_PS_cuts(s_cut)
    this%masses=m
    ! In the HAAG notation Sigma_k is the sum of the diagonal invariants,
    ! i.e. the sum of squared particle masses (hep-ph/0204055, eq. (9)).
    this%tot_mass=sum(this%masses**2)

    call get_approx_s0()

    ! Bring the colour order to a canonical order (first in the list
    ! should be particle 1, i.e., the first incoming particle).
    do i=1,this%next
       if (o(i).eq.1) then
          do j=0,this%next-1
             order(j+1)=o(1+mod(i+j-1,this%next))
          enddo
          exit
       endif
    enddo
    if (verbose) write (*,*) 'Canonical order',order
    ! Define the sets from the colour order. Set 1 contains all the
    ! particles between the first and second incoming particles. Set 2
    ! contains the particles between the second and first incoming
    ! particles.
    this%sets=0
    i=0
    do i=2,this%next
       if (order(i).eq.2) then
          do j=i+1,this%next
             this%sets(0,2)=ibset(this%sets(0,2),order(j)-1)
          enddo
          this%sets(1:i-2,1)=order(2:i-1)
          this%sets(1:this%next-i,2)=order(i+1:this%next)
          exit
       endif
       this%sets(0,1)=ibset(this%sets(0,1),order(i)-1)
    enddo

    if (verbose) then
       write (*,*) "Power in importance sampling:",ip
    endif
  contains

    subroutine get_approx_s0
      implicit none
      ! setup_PS_cuts uses s_cut(2) as the minimum pair invariant.  Reuse
      ! that same bound here; the former trigonometric approximation also
      ! produced a spurious positive cutoff when pT cuts were disabled.
      this%s0=max(s_cut(2),0d0)
    end subroutine get_approx_s0

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
      sigma_ij=0d0
      do i=1,maskr(this%next)
         npart=popcnt(i)
         if (btest(i,0).and.btest(i,1)) then ! both initial particles
            invm_min(i)=0d0
         elseif (btest(i,0).or.btest(i,1)) then  ! one initial particle
            if (npart.eq.2) then
               do j=0,this%next-1
                  if (btest(i,j)) then              
                     invm_max(i)=-s_cut(1) + invm(ibset(0,j))
                  endif
               enddo
            elseif (npart.eq.this%next-2) then
               do j=0,this%next-1
                  if (.not.btest(i,j)) then
                     invm_max(i)=-s_cut(1) + invm(ibset(0,j))
                  endif
               enddo
            endif
         else      ! all final particles
            mass=0d0
            do j=0,this%next-1
               if (btest(i,j)) then
                  mass=mass+sqrt(invm(ibset(0,j)))
                  if (npart.eq.2) then
                     sigma_ij(i) = sigma_ij(i) - invm(ibset(0,j))
                  endif
               endif
            enddo
            invm_min(i)=max(s_cut(2)*(npart)*(npart-1)/2d0,mass**2)
            if (npart.eq.2) then
               sigma_ij(i) = (sigma_ij(i) + invm_min(i))/2d0
            endif
         endif
      enddo
    end subroutine setup_PS_cuts

    subroutine haag_deallocate
      implicit none
      if (allocated(invm)) deallocate(invm)
      if (allocated(invm_min)) deallocate(invm_min)
      if (allocated(invm_max)) deallocate(invm_max)
      if (allocated(sigma_ij)) deallocate(sigma_ij)
      if (allocated(this%masses)) deallocate(this%masses)
      if (allocated(this%pp)) deallocate(this%pp)
      if (allocated(this%p)) deallocate(this%p)
      if (allocated(this%x)) deallocate(this%x)
      if (allocated(this%sets)) deallocate(this%sets)
    end subroutine haag_deallocate
  end subroutine haag_init

  subroutine haag_generate_momenta(this,ps)
    ! Wrapper for the routine that generates the momenta.
    implicit none
    class(phase_space_haag),intent(inout) :: this
    type(psv),intent(inout) :: ps
    real(kind=8) :: mass_sum,mass_sum_linear,total_mass_linear,soft,tau,ycm
    real(kind=8) :: boost_source(0:3),boosted_momentum(0:3)
    real(kind=8),parameter :: random_tolerance=4096d0*epsilon(1d0)
    integer(kind=4) :: i,ix,ix_extra,mm
    ps%jac=-13d0
    if (.not.allocated(ps%p) .or. .not.allocated(ps%x)) return
    if (size(ps%p,1).ne.4 .or. size(ps%p,2).ne.this%next .or. &
         lbound(ps%p,1).ne.0 .or. &
         size(ps%x).lt.this%ndim+this%ndim_extra) return
    ps%p=0d0
    ps%xbjrk=1d0
    if (.not.allocated(this%x) .or. .not.allocated(this%pp) .or. &
         .not.allocated(this%masses) .or. .not.allocated(this%sets)) return
    if (size(this%x).lt.this%ndim+this%ndim_extra .or. &
         size(this%pp,1).ne.4 .or. lbound(this%pp,1).ne.0 .or. &
         lbound(this%pp,2).gt.0 .or. ubound(this%pp,2).lt.this%next .or. &
         size(this%masses).ne.this%next .or. &
         size(this%sets,1).lt.this%next-1 .or. size(this%sets,2).ne.2) return
    if (any(.not.ieee_is_finite(ps%x(1:this%ndim+this%ndim_extra)))) return
    if (any(ps%x(1:this%ndim+this%ndim_extra).lt.-random_tolerance) .or. &
         any(ps%x(1:this%ndim+this%ndim_extra).gt.1d0+random_tolerance)) return
    this%x(1:this%ndim+this%ndim_extra)=max(0d0, &
         min(1d0,ps%x(1:this%ndim+this%ndim_extra)))
    ps%jac=1d0
    ix=0
    ix_extra=this%ndim
    total_mass_linear=sum(this%masses)
    if (this%include_pdf) then
       call generate_initial_state
       if (bad_forward_jac()) return
    endif
    call generate_momenta
    if (bad_forward_jac()) return
    if (ix.ne.this%ndim .or. .not.ieee_is_finite(ps%jac)) then
       ps%p=0d0
       ps%jac=-13d0
       return
    endif
    do i=1,this%next
       if (this%include_pdf) then
          ! Note: 'ycm' is the rapidity needed to go from lab to CM
          ! frame. Hence, here we boost from CM to lab frame with '-ycm'
          ! Passing both allocatable-component sections directly to this
          ! internal routine triggers invalid bounds-check code in gfortran
          ! 16.1 on macOS.  A local result also makes the two workspaces
          ! independent if a compiler elects to use an array temporary.
          boost_source=this%pp(0:3,i)
          call boostz(boost_source,-ycm,boosted_momentum)
          ps%p(0:3,i)=boosted_momentum
       else
          ps%p(0:3,i)=this%pp(0:3,i)
       endif
    enddo
    if (.not.generated_momenta_are_valid(ps%p,this%masses,ps%xbjrk,this%include_pdf)) then
       ps%p=0d0
       ps%jac=-13d0
    endif
  contains

    logical function bad_forward_jac(allow_zero)
      implicit none
      logical,intent(in),optional :: allow_zero
      logical :: zero_is_valid
      zero_is_valid=.false.
      if (present(allow_zero)) zero_is_valid=allow_zero
      if (.not.ieee_is_finite(ps%jac)) then
         ps%jac=-13d0
         bad_forward_jac=.true.
         return
      endif
      if (zero_is_valid) then
         bad_forward_jac=ps%jac.lt.0d0
      else
         bad_forward_jac=ps%jac.le.0d0
      endif
    end function bad_forward_jac

    logical function scale_soft_weight(factor)
      implicit none
      real(kind=8),intent(in) :: factor
      real(kind=8) :: updated_weight
      scale_soft_weight=.false.
      if (.not.ieee_is_finite(factor) .or. factor.le.0d0) then
         ps%jac=-12d0
         return
      endif
      if (.not.safe_phase_space_product(soft,factor,updated_weight)) then
         ps%jac=-12d0
         return
      endif
      soft=updated_weight
      scale_soft_weight=.true.
    end function scale_soft_weight

    logical function divide_soft_weight(denominator)
      implicit none
      real(kind=8),intent(in) :: denominator
      real(kind=8) :: updated_weight
      divide_soft_weight=.false.
      if (denominator.le.0d0 .or. &
           .not.safe_phase_space_ratio(soft,denominator,updated_weight)) then
         ps%jac=-12d0
         return
      endif
      soft=updated_weight
      divide_soft_weight=.true.
    end function divide_soft_weight

    logical function scale_forward_jac(factor)
      implicit none
      real(kind=8),intent(in) :: factor
      real(kind=8) :: updated_jac
      scale_forward_jac=.false.
      if (.not.ieee_is_finite(factor) .or. factor.le.0d0 .or. &
           .not.safe_phase_space_product(ps%jac,factor,updated_jac)) then
         ps%jac=-12d0
         return
      endif
      ps%jac=updated_jac
      scale_forward_jac=.true.
    end function scale_forward_jac

    logical function divide_forward_jac(denominator)
      implicit none
      real(kind=8),intent(in) :: denominator
      real(kind=8) :: updated_jac
      divide_forward_jac=.false.
      if (denominator.le.0d0 .or. &
           .not.safe_phase_space_ratio(ps%jac,denominator,updated_jac)) then
         ps%jac=-12d0
         return
      endif
      ps%jac=updated_jac
      divide_forward_jac=.true.
    end function divide_forward_jac

    real(kind=8) function next_extra_random()
      implicit none
      next_extra_random=0.5d0
      ix_extra=ix_extra+1
      if (ix_extra.gt.size(this%x)) then
         ps%jac=-13d0
         return
      endif
      next_extra_random=this%x(ix_extra)
    end function next_extra_random

    subroutine generate_momenta
      implicit none
      real(kind=8),dimension(0:3,this%next) :: q, qk
      integer(kind=4):: i
      real(kind=8),dimension(0:3) :: Qm,Qnm
      real(kind=8) :: mass1,mass2,mass_linear1,mass_linear2,mass_in,R
      logical :: m1
      integer(kind=8),dimension(:),allocatable :: subperm1(:), subperm2(:),subperm,subperm_rest
      integer(kind=8),dimension(this%next-2) :: perm_final
      real(kind=8),dimension(0:3) :: q1_ref,q2_ref,incoming_first,incoming_second
      integer :: parts

      if(.not.allocated(subperm1)) allocate(subperm1(1:this%next-2))
      if(.not.allocated(subperm2)) allocate(subperm2(1:this%next-2))
      if(.not.allocated(subperm)) allocate(subperm(1:this%next-2))
      if(.not.allocated(subperm_rest)) allocate(subperm_rest(1:this%next-2))
      mass_in=-1d0
      mm = 0
      do i=1,this%next-2
         if (this%sets(i,1) .ne. 0) then
            mm = mm+1
         endif
      enddo
      subperm1 = this%sets(1:,1)
      subperm2 = this%sets(1:,2)

      !! Initial momenta in lab frame
      this%pp(0,ibset(0,0))=this%sqrtshat/2d0
      this%pp(0,ibset(0,1))=this%sqrtshat/2d0
      this%pp(1:2,ibset(0,0))=0d0
      this%pp(1:2,ibset(0,1))=0d0
      this%pp(3,ibset(0,0))= this%pp(0,ibset(0,0))
      this%pp(3,ibset(0,1))=-this%pp(0,ibset(0,1))
      ! Keep fixed local copies for internal procedure calls.  Passing a
      ! section of an allocatable derived-type component directly can trigger
      ! invalid bounds-check code in gfortran 16.1 on macOS.
      incoming_first=this%pp(:,ibset(0,0))
      incoming_second=this%pp(:,ibset(0,1))

      parts = 0

      ! Total momentum P=p_1+p_2
      qk(:,this%next-2) = incoming_first+incoming_second

      ! initiate soft factors
      soft = 1d0

      ! Do m>1 type splitting
      ! First split the antenna to two subantennae: Q_m and Q_{n-m}
      if ((mm .gt. 1).and.(this%next-2-mm .gt. 1)) then 
         mass1=0d0
         mass_linear1=0d0
         do i=1,mm
            mass1=mass1+this%masses(subperm1(i))**2
            mass_linear1=mass_linear1+this%masses(subperm1(i))
         enddo
         mass2=0d0
         mass_linear2=0d0
         do i=1,this%next-2-mm
            mass2=mass2+this%masses(subperm2(i))**2
            mass_linear2=mass_linear2+this%masses(subperm2(i))
         enddo

         call generate_split_Qm_Qnm(this%sqrtshat**2,mass1,mass2,mass_linear1,mass_linear2,mm,Qm,Qnm)
         if (bad_forward_jac()) return
         ! While generating one open sub-antenna, count the squared masses
         ! in the other sub-antenna as already outside the remainder.
         mass_sum=mass2
         mass_sum_linear=mass_linear2

         ! Generate Q_{m} antenna
         if (mm .gt. 2) then
            call basic_antenna(q(0:3,subperm1(mm)),this%masses(subperm1(mm)),qk(0:3,mm-1),-1d0,&
                 Qm,incoming_second,incoming_first,0,.false.,mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm1(mm))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm1(mm))
         else
            call basic_antenna(q(0:3,subperm1(2)),this%masses(subperm1(2)),qk(0:3,1),&
                 this%masses(subperm1(1)),Qm,incoming_second,incoming_first,0,.false.,mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm1(2))**2
            mass_sum=mass_sum+this%masses(subperm1(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm1(2))+this%masses(subperm1(1))
         endif

         do i=1,mm-3
            call basic_antenna(q(0:3,subperm1(mm-i)),this%masses(subperm1(mm-i)),qk(0:3,mm-i-1),-1d0,&
                 qk(0:3,mm-i),q(0:3,subperm1(mm-i+1)),incoming_first,i,.false.,mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm1(mm-i))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm1(mm-i))
         enddo

         if (mm .gt. 2) then
            call basic_antenna(q(0:3,subperm1(2)),this%masses(subperm1(2)),qk(0:3,1),this%masses(subperm1(1)),&
                 qk(0:3,2),q(0:3,subperm1(3)),incoming_first,mm-2,.false.,mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm1(2))**2
            mass_sum=mass_sum+this%masses(subperm1(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm1(2))+this%masses(subperm1(1))
         endif
         q(0:3,subperm1(1)) = qk(0:3,1)

         ! Generate Q_{n-m} antenna
         mass_sum=mass1
         mass_sum_linear=mass_linear1
         if (this%next-2-mm .gt. 2) then
            call basic_antenna(q(0:3,subperm2(this%next-2-mm)),this%masses(subperm2(this%next-2-mm)),qk(0:3,this%next-2-1),&
                 -1d0,Qnm,incoming_first,incoming_second,0,.false.,this%next-2-mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm2(this%next-2-mm))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm2(this%next-2-mm))
         else
            call basic_antenna(q(0:3,subperm2(2)),this%masses(subperm2(2)),qk(0:3,this%next-2-1),&
                 this%masses(subperm2(1)),Qnm,incoming_first,incoming_second,0,.false.,this%next-2-mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm2(1))**2
            mass_sum=mass_sum+this%masses(subperm2(2))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm2(1))+this%masses(subperm2(2))
         endif

         do i=1,(this%next-2-mm)-3
            call basic_antenna(q(0:3,subperm2(this%next-2-mm-i)),this%masses(subperm2(this%next-2-mm-i)),qk(0:3,this%next-2-1-i),&
                 -1d0,qk(0:3,this%next-2-i),q(0:3,subperm2(this%next-2-mm-i+1)),incoming_second,i,.false.,&
                 this%next-2-mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm2(this%next-2-mm-i))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm2(this%next-2-mm-i))
         enddo

         if (this%next-2-mm .gt. 2) then
            call basic_antenna(q(0:3,subperm2(2)),this%masses(subperm2(2)),qk(0:3,mm+1),this%masses(subperm2(1)),&
                 qk(0:3,mm+2),q(0:3,subperm2(3)),incoming_second,this%next-2-mm-2,.false.,this%next-2-mm,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm2(2))**2
            mass_sum=mass_sum+this%masses(subperm2(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm2(2))+this%masses(subperm2(1))
         endif
         q(0:3,subperm2(1)) = qk(0:3,mm+1)




         ! Do m=1 type splitting
      elseif (((mm .eq. 1).or.(this%next-2-mm .eq. 1))) then 
         !write(*,*) 'doing m=1 type splitting'
         m1 = .true.

         R=next_extra_random()
         if (bad_forward_jac()) return
         if (mm .eq. 1) then
            subperm=subperm1
            subperm_rest=subperm2
            q1_ref = incoming_second
            q2_ref = incoming_first
            if (R.lt.0.5d0) then
               q1_ref = incoming_first
               q2_ref = incoming_second
            endif
         elseif (this%next-2-mm .eq. 1) then
            subperm = subperm2
            subperm_rest=subperm1
            q1_ref = incoming_first
            q2_ref = incoming_second
            if (R.lt.0.5d0) then
               q1_ref = incoming_second
               q2_ref = incoming_first
            endif
         endif

         mass_sum=0d0
         mass_sum_linear=0d0
         if (this%next-2 .gt. 2) then
            call basic_antenna(q(0:3,subperm(1)),this%masses(subperm(1)),qk(0:3,this%next-2-1),-1d0,&
                 qk(0:3,this%next-2),q1_ref,q2_ref,0,m1,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm(1))
         else
            call basic_antenna(q(0:3,subperm(1)),this%masses(subperm(1)),qk(0:3,this%next-2-1),&
                 this%masses(subperm_rest(1)),qk(0:3,this%next-2),q1_ref,q2_ref,0,m1,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm(1))**2
            mass_sum=mass_sum+this%masses(subperm_rest(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm(1))+this%masses(subperm_rest(1))
         endif
         if (this%next-2 .gt. 3) then
            call basic_antenna(q(0:3,subperm_rest(this%next-2-1)),this%masses(subperm_rest(this%next-2-1)),&
                 qk(0:3,this%next-2-1-1),-1d0,qk(0:3,this%next-2-1),q2_ref,q1_ref,1,m1,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm_rest(this%next-2-1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm_rest(this%next-2-1))
         elseif (this%next-2 .eq. 3) then 
            call basic_antenna(q(0:3,subperm_rest(2)),this%masses(subperm_rest(2)),qk(0:3,this%next-2-2),&
                 this%masses(subperm_rest(1)),qk(0:3,this%next-2-1),q2_ref,q1_ref,1,.false.,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm_rest(2))**2
            mass_sum=mass_sum+this%masses(subperm_rest(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm_rest(2))+this%masses(subperm_rest(1))
         endif
         do i=2,this%next-2-3
            call basic_antenna(q(0:3,subperm_rest(this%next-2-i)),this%masses(subperm_rest(this%next-2-i)),qk(0:3,this%next-2-i-1),&
                 -1d0,qk(0:3,this%next-2-i),q(0:3,subperm_rest(this%next-2-i+1)),q1_ref,i,.false.,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm_rest(this%next-2-i))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm_rest(this%next-2-i))
         enddo
         if (this%next-2 .gt. 3) then
            call basic_antenna(q(0:3,subperm_rest(2)),this%masses(subperm_rest(2)),qk(0:3,1),&
                 this%masses(subperm_rest(1)),qk(0:3,2),&
                 q(0:3,subperm_rest(3)),q1_ref,this%next-2-2,.false.,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(subperm_rest(2))**2
            mass_sum=mass_sum+this%masses(subperm_rest(1))**2
            mass_sum_linear=mass_sum_linear+this%masses(subperm_rest(2))+this%masses(subperm_rest(1))
         endif
         q(0:3,subperm_rest(1)) = qk(0:3,1)

         ! Do m=0 type splitting
      else  
         !write(*,*) 'doing m=0 type splitting'
         m1 = .false.
         mass_sum=0d0
         mass_sum_linear=0d0

         if (mm .eq. 0) then
            perm_final=subperm2
            q1_ref = incoming_first
            q2_ref = incoming_second
         elseif (this%next-2-mm .eq. 0) then
            perm_final=subperm1
            q1_ref = incoming_second
            q2_ref = incoming_first
         endif

         parts = 0
         do i=1,this%next-2
            parts = parts + ibset(0,perm_final(i)-1)
         enddo

         if (this%next-2 .gt. 2) then
            call basic_antenna(q(0:3,perm_final(this%next-2)),this%masses(perm_final(this%next-2)),&
                 qk(0:3,this%next-2-1),mass_in,qk(0:3,this%next-2),q1_ref,q2_ref,0,m1,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(perm_final(this%next-2))**2
            mass_sum_linear=mass_sum_linear+this%masses(perm_final(this%next-2))
            parts = parts - ibset(0,perm_final(this%next-2)-1)
         else
            call basic_antenna(q(0:3,perm_final(this%next-2)),this%masses(perm_final(this%next-2)),qk(0:3,this%next-2-1),&
                 this%masses(perm_final(1)),qk(0:3,this%next-2),q1_ref,q2_ref,0,m1,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(perm_final(this%next-2))**2
            mass_sum_linear=mass_sum_linear+this%masses(perm_final(this%next-2))
            parts = parts - ibset(0,perm_final(this%next-2)-1)
         endif
         do i=1,this%next-2-3
            call basic_antenna(q(0:3,perm_final(this%next-2-i)),this%masses(perm_final(this%next-2-i)),qk(0:3,this%next-2-i-1),&
                 mass_in,qk(0:3,this%next-2-i),q(0:3,perm_final(this%next-2-i+1)),q2_ref,i,.false.,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(perm_final(this%next-2-i))**2
            mass_sum_linear=mass_sum_linear+this%masses(perm_final(this%next-2-i))
            parts = parts - ibset(0,perm_final(this%next-2-i)-1)
         enddo
         if (this%next-2 .gt. 2) then
            call basic_antenna(q(0:3,perm_final(2)),this%masses(perm_final(2)),qk(0:3,1),&
                 this%masses(perm_final(1)),qk(0:3,2),q(0:3,perm_final(3)),q2_ref,this%next-2-2,.false.,this%next-2,parts)
            if (bad_forward_jac(.true.)) return
            mass_sum=mass_sum+this%masses(perm_final(2))**2
            mass_sum_linear=mass_sum_linear+this%masses(perm_final(2))
            parts = parts - ibset(0,perm_final(2)-1)
         endif
         mass_sum=mass_sum+this%masses(perm_final(1))**2
         mass_sum_linear=mass_sum_linear+this%masses(perm_final(1))
         q(0:3,perm_final(1)) = qk(0:3,1)
      endif

      do i=3,this%next
         this%pp(0:3,i) = q(0:3,i)
      enddo

      ! Compute the weight (i.e. jacobian)
      if (ix.ne.this%ndim .or. .not.ieee_is_finite(soft)) then
         ps%jac=-12d0
         return
      endif
      if (soft.le.0d0) then
         ps%jac=-12d0
         return
      endif
      ! The usual 2*pi factors for the phase-space
      if (.not.scale_forward_jac(soft)) return
      if (.not.divide_forward_jac((2d0*pi)**(3*(this%next-2)-4))) return
      if (.not.divide_forward_jac(2d0*this%sqrtshat**2)) return
      if (.not.ieee_is_finite(ps%jac)) then
         ps%jac=-12d0
      elseif (ps%jac.le.0d0) then
         ps%jac=-12d0
      endif

      !write(*,*) 'jac',ps%jac

    end subroutine generate_momenta

    subroutine generate_initial_state
      implicit none
      call generate_tau
      if (bad_forward_jac()) return
      call generate_y
      if (bad_forward_jac()) return
      this%sqrtshat=sqrt(tau)*this%sqrts
      ps%xbjrk(1)=sqrt(tau)*exp(ycm)
      ps%xbjrk(2)=sqrt(tau)*exp(-ycm)
    end subroutine generate_initial_state

    subroutine generate_tau
      implicit none
      real(kind=8) :: smin,smax,shat
      smin=max(sum(this%masses)**2, &
           (this%next-2)*(this%next-3)*this%s0/2d0)
      smax=this%sqrts**2
      if (smin.lt.0d0 .or. smin.ge.smax) then
         ps%jac=-10d0
         return
      endif
      ix=ix+1
      if (smin.gt.0d0) then
         call random_to_var(this%x(ix),ip_shat,smin,smax,shat,ps%jac)
      else
         call random_to_var(this%x(ix),0d0,smin,smax,shat,ps%jac)
      endif
      if (bad_forward_jac()) return
      if (.not.ieee_is_finite(shat)) then
         ps%jac=-10d0
         return
      endif
      if (shat.le.0d0) then
         ps%jac=-10d0
         return
      endif
      tau=shat/smax
      if (.not.divide_forward_jac(smax)) return
    end subroutine generate_tau

    subroutine generate_y
      implicit none
      real(kind=8) ::  ymin,ymax
      if (tau.le.0d0 .or. tau.gt.1d0) then
         ps%jac=-10d0
         return
      endif
      ymin= log(tau)/2d0
      ymax=-log(tau)/2d0
      ix=ix+1
      call random_to_var(this%x(ix),0d0,ymin,ymax,ycm,ps%jac)
    end subroutine generate_y

    subroutine basic_antenna(p1,mass1,p2,mass2,P,q1,q2,i,m1,maxn,parts)
      ! Incoming momentum: P
      ! Reference momenta: q1,q2
      ! Outgoing momenta: p1, p2
      implicit none
      real(kind=8),dimension(0:3),intent(in) :: P, q1,q2
      real(kind=8),dimension(0:3),intent(out) :: p1, p2
      real(kind=8),intent(in) :: mass1,mass2
      real(kind=8),dimension(0:3) :: q1_cmf,q2_cmf,p1_cmf, p2_cmf,Pm,P_cmf
      real(kind=8) :: esum,costheta,sintheta
      real(kind=8) :: s1, s2, a1, a2, s, gs
      integer :: i,k
      real(kind=8):: h
      real(kind=8),dimension(3) :: solution
      real(kind=8) ::  beta,a1cut,a2cut,q1norm2,q2norm2,angle_den,cut_den1,cut_den2
      integer :: maxn
      logical :: m1
      real(kind=8) z_sign
      double precision :: w1,w2,w,R
      integer :: term 
      integer :: parts

      k = maxn - i   !  k is the number of particles remaining to generate
      if (k.lt.2 .or. mass1.lt.0d0 .or. (k.lt.3 .and. mass2.lt.0d0)) then
         ps%jac=-1d0
         return
      endif
      s1 = mass1**2     ! mass of the final state particle to be generated
      z_sign=sign(1d0,q2(3))
      s = (dot(P,P)) ! Incoming inv mass

      ! boost qi to CMF (P rest frame)
      if (dot(P,P).le.0d0) then
         ps%jac=-1d0
         return
      endif
      esum=dsqrt(dot(P,P))
      Pm(0)=P(0)
      Pm(1:3)=-P(1:3)
      P_cmf(0)=esum
      P_cmf(1:3)=(/0d0,0d0,0d0/)
      call boostm(q1,Pm,esum,q1_cmf)
      call boostm(q2,Pm,esum,q2_cmf)

      ! Split into the long (L) decomposition of the q1 (relevant for massive)
      q1norm2=threedot(q1_cmf(1:3),q1_cmf(1:3))
      q2norm2=threedot(q2_cmf(1:3),q2_cmf(1:3))
      if (q1norm2.le.0d0 .or. q2norm2.le.0d0) then
         ps%jac=-2d0
         return
      endif
      if (mass1.eq.0d0) then
         beta=1d0
      else
         if (q1_cmf(0).le.0d0) then
            ps%jac=-2d0
            return
         endif
         beta=dsqrt(q1norm2)/q1_cmf(0)
         if (beta.le.vtiny) then
            ps%jac=-2d0
            return
         endif
         q1_cmf(1:3)=q1_cmf(1:3)/beta
         q1norm2=threedot(q1_cmf(1:3),q1_cmf(1:3))
      endif

      angle_den=dsqrt(q1norm2)*dsqrt(q2norm2)
      if (angle_den.le.vtiny*max(spacing(0d0),abs(q1_cmf(0)*q2_cmf(0)))) then
         ps%jac=-2d0
         return
      endif
      ! angles between q1,q2 in CMF_k frame
      costheta=threedot(q1_cmf(1:3),q2_cmf(1:3))/angle_den
      if (abs(costheta).gt.1d0+1d-10) then
         ps%jac=-2d0
         return
      endif
      costheta=max(-1d0,min(1d0,costheta))

      !Generate s2
      if (m1 .and. (i .eq. 0)) then 
         ! If m=1 type first splitting -> do also a1,phi sampling here!
         if (k.eq.2) s2=mass2**2
         call generate_first_single(k,s,s1,q1_cmf,P_cmf,a1,s2)
         if (bad_forward_jac()) return
         a2 = 300d0
         goto 20
      else
         if (k .ge. 3) then
            call generate_s2(k,s,s1,s2,q1_cmf,P_cmf)
            if (bad_forward_jac()) return
         else
            s2 = mass2**2
            gs = 1d0
            if (.not.divide_forward_jac(gs)) return
         endif
      endif

      ! cuts on a1 and a2
      cut_den1=dot(q1_cmf,P_cmf)
      cut_den2=dot(q2_cmf,P_cmf)
      if (cut_den1.le.0d0 .or. cut_den2.le.0d0) then
         ps%jac=-2d0
         return
      endif
      a2cut=(k-1)*(this%s0/2d0)/cut_den2 ! should be (k+1) dont change!
      a1cut=this%s0/(2d0*cut_den1)

      h= 1d-8 ! for using the "h"-technique
      !h = 0d0  ! for using the partial decomposition
      !h = -1d0  ! for using actually h=0

      ! The two partial-decomposition weights select term 1 versus term 2
      ! only for h=0.  The default h-technique always uses term 1, so avoid
      ! evaluating the unused (and more singular) term-2 weight there.
      if (h.eq.0d0) then
         R=next_extra_random()
         if (bad_forward_jac()) return
         call get_partial_weights(w1,w2,s,s1,s2,a1cut,a2cut,h,costheta)
         w=w1+w2
         if (.not.ieee_is_finite(w1) .or. .not.ieee_is_finite(w2) .or. &
              .not.ieee_is_finite(w)) then
            ps%jac=-9d0
            return
         endif
         if (w1.lt.0d0 .or. w2.lt.0d0 .or. w.le.0d0) then
            ps%jac=-9d0
            return
         endif
         if (R.lt.w1/w) then
            term=1
         else
            term=2
         endif
      endif

      if (h.eq.0d0) then
         if (((i .eq. 0) .or. (m1 .and. (i .le. 1)))) then
            call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
            if (bad_forward_jac()) return
            call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
            if (bad_forward_jac()) return
         else
            if (term.eq.1) then
               call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
               if (bad_forward_jac()) return
               call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
               if (bad_forward_jac()) return
            elseif (term.eq.2) then
               call generate_a2_term2(i,m1,maxn,s,s1,s2,costheta,a2cut,beta,a2)
               if (bad_forward_jac()) return
               call generate_a1_term2(i,m1,maxn,a2,s,s1,s2,costheta,a1cut,a1)
               if (bad_forward_jac()) return
            endif
         endif

      elseif (abs(h).gt.0d0) then
         call generate_a1_term1(i,m1,maxn,s,s1,s2,costheta,a1cut,beta,h,a1)
         if (bad_forward_jac()) return
         call generate_a2_term1(i,m1,maxn,a1,s,s1,s2,costheta,a2cut,h,a2)
         if (bad_forward_jac()) return
      endif

      ! Mapping back to the momenta p1,p2 in CMF
20    p1_cmf(0) = (s+s1-s2)/(2D0*sqrt(s))
      solution = solver(s,s1,s2,q1_cmf,q2_cmf,z_sign,a1,a2,p1_cmf(0),P_cmf,costheta)
      if (bad_forward_jac()) return

      p1_cmf(1) = solution(1)
      p1_cmf(2) = solution(2)
      p1_cmf(3) = solution(3)

      p2_cmf(0) = sqrt(s) - p1_cmf(0)
      p2_cmf(1:3) = -p1_cmf(1:3)

      ! Boost back to lab frame
      call boostm(p1_cmf,P,esum,p1)
      call boostm(p2_cmf,P,esum,p2)

    end subroutine basic_antenna

    subroutine get_partial_weights(w1,w2,s,s1,s2,a1cut,a2cut,h,cos)
      implicit none
      real(kind=8) :: w1,w2
      real(kind=8) :: s,s1,s2,cos
      real(kind=8) :: h1,a1min,a1max,a2min,a2max,f_h1,h
      real(kind=8) :: a1cut,a2cut
      real(kind=8),dimension(3) :: buff
      real(kind=8) :: amin,amax,bmin,bmax,root,den

      w1=-1d0
      w2=-1d0
      if (s.le.0d0) return
      h1=0d0
      root=sqrt_kallen(1d0,s1/s,s2/s)
      if (bad_forward_jac()) return
      a1max=0.5d0*(1d0-(s2-s1)/s+root)
      a1min=0.5d0*(1d0-(s2-s1)/s-root)
      if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) then
         a1min=a1cut
      endif

      buff = f_func_term1(-h1,cos,s,s1,s2,h)
      if (bad_forward_jac()) return
      f_h1 = buff(2)
      if (f_h1.le.0d0) return
      buff = f_func_term1(a1min,cos,s,s1,s2,h)
      if (bad_forward_jac()) return
      den=a1min+h1+buff(2)+f_h1
      if (den.eq.0d0) return
      Amin=(a1min+h1+buff(2)-f_h1)/den
      buff = f_func_term1(a1max,cos,s,s1,s2,h)
      if (bad_forward_jac()) return
      den=a1max+h1+buff(2)+f_h1
      if (den.eq.0d0) return
      Amax=(a1max+h1+buff(2)-f_h1)/den
      if (Amin.le.0d0 .or. Amax.le.Amin) return
      w1 = (1d0/f_h1)*(log(Amax)-log(Amin))

      a2max=0.5d0*(1d0+(s2-s1)/s+root)
      a2min=0.5d0*(1d0+(s2-s1)/s-root)
      if ((a2cut.gt.a2min).and.(a2cut.lt.a2max)) then
         a2min=a2cut
      endif

      buff = f_func_term2(-h1,cos,s,s1,s2)
      if (bad_forward_jac()) return
      f_h1 = buff(2)
      if (f_h1.le.0d0) return
      buff = f_func_term2(a2min,cos,s,s1,s2)
      if (bad_forward_jac()) return
      den=a2min+h1+buff(2)+f_h1
      if (den.eq.0d0) return
      bmin=(a2min+h1+buff(2)-f_h1)/den
      buff = f_func_term2(a2max,cos,s,s1,s2)
      if (bad_forward_jac()) return
      den=a2max+h1+buff(2)+f_h1
      if (den.eq.0d0) return
      bmax=(a2max+h1+buff(2)-f_h1)/den
      if (bmin.le.0d0 .or. bmax.le.bmin) return
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

    real(kind=8) function sqrt_kallen(a,b,c)
      implicit none
      real(kind=8),intent(in) :: a,b,c
      real(kind=8) :: value,scale
      sqrt_kallen=0d0
      value=kallen(a,b,c)
      scale=max(1d0,a**2,b**2,c**2,abs(a*b),abs(a*c),abs(b*c))
      if (value.lt.-1d-12*scale) then
         ps%jac=-9d0
         return
      endif
      sqrt_kallen=dsqrt(max(value,0d0))
    end function sqrt_kallen

    subroutine generate_split_Qm_Qnm(s,mass1,mass2,mass_linear1,mass_linear2,mn,Qm,Qnm)
      implicit none
      real(kind=8),dimension(0:3),intent(out) :: Qm,Qnm
      real(kind=8),dimension(0:3) :: dummy
      integer :: mn
      real(kind=8) :: mass1,mass2,mass_linear1,mass_linear2
      real(kind=8) :: c1,c2,m,dum,s,s1,s2,E1,E2,phi,Qt,Qz
      real(kind=8) :: RHS,qzmax,qzmin,sum_w,R,weight_den,weight_log,weight_scale
      real(kind=8) :: log_arg1,log_arg2,log_arg3,log_arg4,map_den,exponent,rad
      real(kind=8),dimension(4) :: g1,g2,g3,w    
      integer :: i,pick
      real(kind=8) :: dm,a1,a1min,a1max,a2,comm,root,a1cut
      real(kind=8),dimension(3) :: solution

      ! Generate s1,s2 invariants
      Qm=0d0
      Qnm=0d0
      dum=1d0
      if (s.le.0d0 .or. mass1.lt.0d0 .or. mass2.lt.0d0 .or. &
           mass_linear1.lt.0d0 .or. mass_linear2.lt.0d0) then
         ps%jac=-4d0
         return
      endif
      m=dsqrt(s)
      c1=dsqrt(max(mass_linear1**2,mass1+mn*(mn-1)*this%s0/2d0))
      c2=dsqrt(max(mass_linear2**2, &
           mass2+(this%next-2-mn)*(this%next-2-mn-1)*this%s0/2d0))
      if (c1+c2.ge.m) then
         ps%jac=-4d0
         return
      endif

      if (.not. flat_split) then
         if (c1.le.0d0 .or. c2.le.0d0) then
            ps%jac=-4d0
            return
         endif
         ix = ix +1
         call random_to_var(this%x(ix),-1d0,c1**2,(m-c2)**2,s1,dum)
         if (dum.le.0d0) then
            ps%jac=-4d0
            return
         endif
         ix = ix +1
         call random_to_var(this%x(ix),-1d0,c2**2,(m-dsqrt(s1))**2,s2,dum)
         if (dum.le.0d0 .or. (m-dsqrt(s1))**2.le.c2**2) then
            ps%jac=-4d0
            return
         endif
         if (.not.scale_soft_weight(log((m-sqrt(s1))**2)-log(c2**2))) return
         if (.not.scale_soft_weight(log((m-c2)**2/(c1**2)))) return
         if (.not.scale_forward_jac(s1*s2)) return
      endif

      if (flat_split) then
         ix = ix + 1
         call random_to_var(this%x(ix),0d0,c1**2,(m-c2)**2,s1,dum)
         if (dum.le.0d0) then
            ps%jac=-4d0
            return
         endif
         if (.not.scale_soft_weight((m-c2)**2-c1**2)) return
         ix = ix + 1
         call random_to_var(this%x(ix),0d0,c2**2,(m-dsqrt(s1))**2,s2,dum)
         if (dum.le.0d0) then
            ps%jac=-4d0
            return
         endif
         if (.not.scale_soft_weight((m-dsqrt(s1))**2-c2**2)) return
      endif

      if (.not.a1_split) then
         E1 = (s+s1-s2)/(2d0*dsqrt(s))
         E2 = dsqrt(s) - E1

         ix = ix + 1
         call random_to_var(this%x(ix),0d0,0d0,2d0*pi,phi,ps%jac)
         if (bad_forward_jac()) return
         !soft = soft*2d0*pi

         ! multichanneling for Qz sampling
         ! Not the same input as from paper!!! 
         g1 = (/-1d0,+1d0,+1d0,-1d0/)
         g2 = (/1d0,1d0,-1d0,-1d0/)
         g3 = (/1d0,-1d0,1d0,-1d0/)
         rad=E1**2-s1
         if (E1.le.0d0 .or. E2.le.0d0 .or. &
              rad.lt.-1d-12*max(spacing(0d0),E1**2,s1)) then
            ps%jac=-4d0
            return
         endif
         qzmax=dsqrt(max(rad,0d0))
         qzmin=-qzmax

         do i=1,4
            log_arg1=E1+g2(i)*qzmax
            log_arg2=E2+g3(i)*qzmax
            log_arg3=E1+g2(i)*qzmin
            log_arg4=E2+g3(i)*qzmin
            if (min(log_arg1,log_arg2,log_arg3,log_arg4).le.0d0) then
               ps%jac=-4d0
               return
            endif
            weight_den=E2+g1(i)*E1
            weight_scale=max(E1,E2)
            if (abs(weight_den).le.1d-10*weight_scale) then
               if (g1(i).ne.-1d0 .or. g3(i).ne.g2(i)) then
                  ps%jac=-4d0
                  return
               endif
               w(i)=g2(i)/(4d0*E1*E2)*(1d0/(E1+g2(i)*qzmin)-1d0/(E1+g2(i)*qzmax))
            else
               weight_log=g2(i)*(log(log_arg1)-log(log_arg2)-log(log_arg3)+log(log_arg4))
               w(i)=weight_log/(4d0*E1*E2*weight_den)
            endif
         enddo
         !w(1) = (1d0/(4d0*E1*E2))*1d0/(E1-E2)*((log(E1+max)-log(E2+max))-(log(E1+min)-log(E2+min)))
         !w(2) = (1d0/(4d0*E1*E2))*1d0/(E1+E2)*((log(E1+max)-log(E2-max))-(log(E1+min)-log(E2-min)))
         !w(3) = (1d0/(4d0*E1*E2))*1d0/(E1+E2)*((-log(E1-max)+log(E2+max))-(-log(E1-min)+log(E2+min)))
         !w(4) = (1d0/(4d0*E1*E2))*1d0/(E1-E2)*((-log(E1-max)+log(E2-max))-(-log(E1-min)+log(E2-min)))
         sum_w = sum(w)
         if (any(.not.ieee_is_finite(w)) .or. .not.ieee_is_finite(sum_w)) then
            ps%jac=-4d0
            return
         endif
         if (any(w.lt.0d0) .or. sum_w.le.0d0) then
            ps%jac=-4d0
            return
         endif
         R=next_extra_random()
         if (bad_forward_jac()) return

         if (R .lt. w(1)/sum_w) then
            pick = 1
         elseif (R .lt. (w(2)+w(1))/sum_w) then
            pick = 2
         elseif (R .lt. (w(3)+w(2)+w(1))/sum_w) then
            pick = 3
         else
            pick = 4
         endif

         if (.not. flat_split) then
            ix = ix + 1
            call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
            if (bad_forward_jac()) return
            exponent=R*w(pick)*4d0*E1*E2*(E2+g1(pick)*E1)/g2(pick)
            if (exponent.gt.log(huge(1d0)) .or. exponent.lt.log(tiny(1d0))) then
               ps%jac=-4d0
               return
            endif
            RHS=((E1+g2(pick)*qzmin)/(E2+g3(pick)*qzmin))*exp(exponent)
            map_den=g3(pick)*RHS-g2(pick)
            if (.not.ieee_is_finite(RHS) .or. .not.ieee_is_finite(map_den)) then
               ps%jac=-4d0
               return
            endif
            if (abs(map_den).le.vtiny) then
               ps%jac=-4d0
               return
            endif
            Qz=(E1-E2*RHS)/map_den
            if (.not.ieee_is_finite(Qz)) then
               ps%jac=-4d0
               return
            endif
            if (E1**2-Qz**2.le.0d0 .or. E2**2-Qz**2.le.0d0) then
               ps%jac=-4d0
               return
            endif
            if (.not.scale_soft_weight(sum_w)) return
            if (.not.scale_forward_jac((E1**2-Qz**2)*(E2**2-Qz**2))) return
         endif

         if (flat_split) then
            ix = ix + 1
            call random_to_var(this%x(ix),0d0,qzmin,qzmax,Qz,ps%jac)
            if (bad_forward_jac()) return
            !soft = soft*(max-min)
         endif

         if (.not.divide_soft_weight(4d0*dsqrt(s))) return
         if (E1**2-s1-Qz**2.lt.-1d-10*max(spacing(0d0),E1**2,s1,Qz**2)) then
            ps%jac=-4d0
            return
         endif
         Qt=dsqrt(max(E1**2-s1-Qz**2,0d0))
         Qm = (/E1,Qt*cos(phi),Qt*sin(phi),Qz/)
         Qnm = (/E2,-Qt*cos(phi),-Qt*sin(phi),-Qz/)
      else 
         ! do a1 sampling
         comm = 0.5d0*(s+s1-s2)/s
         rad=comm**2-s1/s
         if (rad.lt.-1d-12*max(1d0,comm**2,abs(s1/s))) then
            ps%jac=-4d0
            return
         endif
         root=dsqrt(max(rad,0d0))
         a1min = comm - root
         a1max = comm + root
         a1cut = 0.5d0*this%s0*mn/(s/2d0)
         if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) a1min=a1cut
         a1cut = 1d0-0.5d0*this%s0*(this%next-2-mn)/(s/2d0)
         if ((a1cut.gt.a1min).and.(a1cut.lt.a1max)) a1max=a1cut
         ix = ix + 1
         call random_to_var(this%x(ix),-1d0,a1min,a1max,a1,ps%jac)
         if (bad_forward_jac()) return
         if (.not.scale_forward_jac(a1)) return
         !soft = soft*log(a1max/a1min)

         dummy=(/0d0,0d0,0d0,0d0/)
         dm=1d0
         a2=300d0
         E1 = (s+s1-s2)/(2d0*dsqrt(s))
         E2 = dsqrt(s) - E1
         solution = solver(s,s1,s2,dummy,dummy,dm,a1,a2,E1,dummy,dm) ! use same solver 
         if (bad_forward_jac()) return

         Qm =  (/E1, solution(1), solution(2), solution(3)/)
         Qnm = (/E2,-solution(1),-solution(2),-solution(3)/)
      endif
    end subroutine generate_split_Qm_Qnm


    subroutine generate_first_single(k,s,s1,q1_cmf,P_cmf,a1,s2)
      implicit none
      real(kind=8) :: s,s1
      integer :: k
      real(kind=8),intent(out) :: a1
      real(kind=8),intent(inout) :: s2
      real(kind=8) :: Lambda,Sigma,Delta,sigmak,smin,smax,smax_force,remaining_mass
      real(kind=8) :: A,C,R,gs,a2,mu,a1min,a1max
      real(kind=8),dimension(0:3) :: q1_cmf,P_cmf

      a1=0d0
      if (s.le.0d0 .or. s1.lt.0d0 .or. dot(q1_cmf,P_cmf).le.0d0) then
         ps%jac=-5d0
         return
      endif
      Sigma=this%tot_mass-mass_sum-s1
      if (Sigma.lt.-1d-12*max(spacing(0d0),this%tot_mass,mass_sum,s1)) then
         ps%jac=-5d0
         return
      endif
      Sigma=max(Sigma,0d0)
      remaining_mass=total_mass_linear-mass_sum_linear-dsqrt(s1)
      if (remaining_mass.lt.-1d-12*max(spacing(0d0),total_mass_linear, &
           mass_sum_linear,dsqrt(s1))) then
         ps%jac=-5d0
         return
      endif
      remaining_mass=max(remaining_mass,0d0)
      Lambda=max(remaining_mass**2,Sigma+(k-1)*(k-2)*this%s0/2d0)
      Delta = s1+2D0*(k-1)*this%s0/2D0
      sigmak =  s1  
      ! NOTE: added extra upper limit for massive case!
      if (Delta .lt. (2d0*dsqrt(s1*s)-s1)) then
         Delta = (2d0*dsqrt(s1*s)-s1)
      endif
      smax_force = s*(1-(k)*(this%s0/2d0)/dot(q1_cmf,P_cmf))
      smin = Lambda
      smax = s - Delta
      if (smax_force.lt.smax) smax=smax_force

      A = Sigma
      if (k.gt.2) then
         if (.not. this%flat_mode) then
            if (smin.ge.smax .or. smax-A.le.0d0 .or. smin-A.le.0d0) then
               ps%jac=-5d0
               return
            endif
            ix = ix +1
            call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
            if (bad_forward_jac()) return
            C = ((smax-A)/(smin-A))**R
            s2 = (smin-A)*C + A
            if (.not.scale_soft_weight(log(smax-A)-log(smin-A))) return
            gs = 1d0/(s2-Sigma)
            if (.not.divide_forward_jac(gs)) return
            a1 = a1_m1(s,s1,s2,q1_cmf,P_cmf,k)
            if (bad_forward_jac()) return
            mu = (s2-s1)/s
            a2 = a1 + mu
            if (min(a1,1d0-a1,a2,1d0-a2).le.0d0) then
               ps%jac=-5d0
               return
            endif
            if (.not.scale_forward_jac(a1*(1d0-a1)*(1d0-a2)*a2)) return
         elseif (this%flat_mode) then
            ix = ix +1
            call random_to_var(this%x(ix),0d0,smin,smax,s2,ps%jac)
            if (bad_forward_jac()) return
            !soft = soft*(smax-smin)
            gs = 1d0
            if (.not.divide_forward_jac(gs)) return
            a1max = 0.5d0*(1d0+(s1-s2)/(s)+sqrt_kallen(1d0,s1/s,s2/s))
            a1max = a1max-0.00001d0
            a1min = 0.5d0*(1d0+(s1-s2)/(s)-sqrt_kallen(1d0,s1/s,s2/s))
            if ((a1min .lt. (this%s0/2d0)/(dot(q1_cmf,P_cmf)) ) .and.&
                 a1max .gt. (this%s0/2d0)/(dot(q1_cmf,P_cmf))) then
               a1min = (this%s0/2d0)/(dot(q1_cmf,P_cmf))
            endif
            ix = ix +1
            call random_to_var(this%x(ix),0d0,a1min,a1max,a1,ps%jac)
            if (bad_forward_jac()) return
            !soft = soft*(a1max-a1min)
         endif

      elseif (k.eq.2) then
         a1 = a1_m1(s,s1,s2,q1_cmf,P_cmf,k)
         if (bad_forward_jac()) return
         if (.not.scale_forward_jac(a1*(1d0-a1))) return
      endif

    end subroutine generate_first_single

    subroutine generate_s2(k,s,s1,s2,q1_cmf,P_cmf)
      implicit none
      integer :: k
      real(kind=8),intent(out) :: s2
      real(kind=8),dimension(0:3) :: q1_cmf,P_cmf
      real(kind=8) :: s,s1,A,B,C,R,gs,ratio_min,ratio_max,den_min,den_max,map_den
      real(kind=8) :: Lambda,Delta,Sigma,sigmak,Sigmaold,smin,smax,smax_force,remaining_mass
      double precision :: scut

      ! Include also initial momenta in the limits!

      scut = this%s0

      ! Sigma_{k-1} is the diagonal-mass sum of the ungenerated
      ! remainder; Sigma_k additionally contains the emitted particle.
      if (s.le.0d0 .or. s1.lt.0d0 .or. dot(q1_cmf,P_cmf).le.0d0) then
         ps%jac=-8d0
         return
      endif
      Sigma=this%tot_mass-mass_sum-s1
      if (Sigma.lt.-1d-12*max(spacing(0d0),this%tot_mass,mass_sum,s1)) then
         ps%jac=-8d0
         return
      endif
      Sigma=max(Sigma,0d0)
      remaining_mass=total_mass_linear-mass_sum_linear-dsqrt(s1)
      if (remaining_mass.lt.-1d-12*max(spacing(0d0),total_mass_linear, &
           mass_sum_linear,dsqrt(s1))) then
         ps%jac=-8d0
         return
      endif
      remaining_mass=max(remaining_mass,0d0)
      Lambda=max(remaining_mass**2,Sigma+(k-1)*(k-2)/2d0*scut)
      Sigmaold=Sigma+s1
      Delta = s1 + (k-1)*scut
      sigmak =  s1   !always just the previous particle mass
      ! NOTE: added extra upper limit for massive case!
      if (Delta .lt. (2d0*dsqrt(s1*s)-s1)) then
         Delta = (2d0*dsqrt(s1*s)-s1)
      endif
      smin = Lambda
      smax = s - Delta
      smax_force = s*(1-(this%s0/2d0)/dot(q1_cmf,P_cmf))
      if (smax .gt. smax_force) then
         smax = smax_force
      endif
      A = Sigma
      B = s - sigmak

      ! S limits exactly same as in COMIX! 

      if (smin.gt.smax) then
         ps%jac=-3d0
         return
      endif

      if ((.not.open) .and. (.not. this%flat_mode)) then
         den_max=s-sigmak-smax
         den_min=s-sigmak-smin
         if (den_max.le.0d0 .or. den_min.le.0d0 .or. &
              smax-Sigma.le.0d0 .or. smin-Sigma.le.0d0) then
            ps%jac=-8d0
            return
         endif
         ratio_max=(smax-Sigma)/den_max
         ratio_min=(smin-Sigma)/den_min
         ix = ix +1
         call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
         if (bad_forward_jac()) return
         C = ((smax - A)*(B-smin)/((B-smax) * (smin - A)))**R
         map_den=B-smin+(smin-A)*C
         if (abs(map_den).le.vtiny*max(spacing(0d0),abs(B),abs(smin),abs(A))) then
            ps%jac=-8d0
            return
         endif
         s2=(A*(B-smin)+B*(smin-A)*C)/map_den
         if (.not.scale_soft_weight(log(ratio_max)-log(ratio_min))) return
         gs = (s-Sigmaold)/((s-sigmak-s2)*(s2-Sigma))
         if (.not.ieee_is_finite(gs)) then
            ps%jac=-8d0
            return
         endif
         if (gs.le.0d0) then
            ps%jac=-8d0
            return
         endif
         if (.not.divide_forward_jac(gs)) return
      endif

      if (this%flat_mode) then
         ix = ix +1
         call random_to_var(this%x(ix),0d0,smin,smax,s2,ps%jac)
         if (bad_forward_jac()) return
         !soft = soft*(smax-smin)
         gs = 1d0
         if (.not.divide_forward_jac(gs)) return
      endif
      if (open) then
         ix = ix +1
         call random_to_var(this%x(ix),-1d0,smin,smax,s2,ps%jac)
         if (bad_forward_jac()) return
         !soft = soft*log(smax/smin)
         if (.not.scale_forward_jac(s2)) return
      endif
    end subroutine generate_s2


    subroutine generate_a1_term1(i,m1,maxn,s,s1,s2,cos,a1cut,beta,h_in,a1)
      implicit none
      integer :: i,maxn
      real(kind=8),intent(out) :: a1
      real(kind=8) :: s,s1,s2,cos
      logical :: m1
      real(kind=8) :: a1cut
      real(kind=8) :: h1,h,h_in,a1min,a1max,a1max_force,f_h1,beta
      real(kind=8),dimension(3) :: buff
      real(kind=8) :: Amin,Amax,Atilde,R,v,wsq,kappa,root,den

      a1=0d0
      if (s.le.0d0 .or. beta.le.0d0) then
         ps%jac=-9d0
         return
      endif
      h=h_in
      h1 = ((1d0-beta)/(2d0*beta))*(1d0+(s1-s2)/s)
      root=sqrt_kallen(1d0,s1/s,s2/s)
      if (bad_forward_jac()) return
      a1max=0.5d0*(1d0-(s2-s1)/s+root)
      a1min=0.5d0*(1d0-(s2-s1)/s-root)

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
      if (bad_forward_jac()) return
      f_h1 = buff(2)
      if (f_h1.le.0d0) then
         ps%jac=-9d0
         return
      endif
      buff = f_func_term1(a1min,cos,s,s1,s2,h)
      if (bad_forward_jac()) return
      den=a1min+h1+buff(2)+f_h1
      if (abs(den).le.vtiny) then
         ps%jac=-9d0
         return
      endif
      Amin=(a1min+h1+buff(2)-f_h1)/den
      buff = f_func_term1(a1max,cos,s,s1,s2,h)
      if (bad_forward_jac()) return
      den=a1max+h1+buff(2)+f_h1
      if (abs(den).le.vtiny) then
         ps%jac=-9d0
         return
      endif
      Amax=(a1max+h1+buff(2)-f_h1)/den
      if (Amax.le.0d0) then
         ps%jac=-9d0
         return
      endif
      if (Amin.le.0d0) Amin=Amax*1d-8
      if (Amin.ge.Amax) then
         ps%jac=-9d0
         return
      endif

      ! Now generate a1, depending on which distribution
      if ((.not. open) .and. (.not. this%flat_mode)) then
         if (( ((i .ne. 0) .and. (.not. m1)) .or. (maxn .ne. this%next-2))) then
            ! Sample with Pi^{-1/2} factor
            ix = ix + 1
            call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
            if (bad_forward_jac()) return
            ! Only the linear coefficient is needed here. It is independent
            ! of a1, which has not been generated yet.
            buff = f_func_term1(0d0,cos,s,s1,s2,h)
            if (bad_forward_jac()) return
            v = buff(1)/2d0
            buff = f_func_term1(0d0,cos,s,s1,s2,h)
            if (bad_forward_jac()) return
            wsq = buff(2)**2  ! w^2
            Atilde = (Amin**(1d0-R))*(Amax)**R
            if (.not.ieee_is_finite(Atilde)) then
               ps%jac=-9d0
               return
            endif
            if (Atilde.le.0d0 .or. Atilde.ge.1d0-vtiny) then
               ps%jac=-9d0
               return
            endif
            kappa = - h1 + f_h1*(1d0+Atilde)/(1d0-Atilde)
            if (v+kappa.lt.h*1d-5) then
               a1=a1max
            else
               a1 = ((kappa**2) - wsq)/(2d0*(v+kappa))
            endif
            if (.not.clamp_haag_interval(a1,a1min,a1max,1d-6)) then
               ps%jac=-9d0
               return
            endif
            if (.not.scale_soft_weight(log(Amax)-log(Amin))) return
            if (.not.divide_soft_weight(f_h1)) return
            if (.not.scale_forward_jac(a1)) return
            if (.not.scale_forward_jac(beta)) return
         elseif (((i .eq. 0) .or. (m1 .and. (i .le. 1)))) then
            ! Sample with 1/x 
            ix = ix + 1
            call random_to_var(this%x(ix),-1d0,a1min,a1max,a1,ps%jac)
            if (bad_forward_jac()) return
            !soft = soft*log(a1max/a1min)
            !ps%jac = ps%jac*a1 

         endif
      endif

      if (this%flat_mode) then
         ix = ix +1
         call random_to_var(this%x(ix),0d0,a1min,a1max,a1,ps%jac)
         if (bad_forward_jac()) return
         !soft = soft*(a1max-a1min)
      endif
      if (open) then
         ix = ix + 1
         call random_to_var(this%x(ix),-1d0,a1min,a1max,a1,ps%jac)
         if (bad_forward_jac()) return
         !soft = soft*log(a1max/a1min)
         !ps%jac = ps%jac*a1
      endif
    end subroutine generate_a1_term1

    subroutine generate_a1_term2(i,m1,maxn,a2,s,s1,s2,cos,a1cut,a1)
      implicit none
      integer :: i,maxn
      real(kind=8) :: a2,s,s1,s2,cos,h,a1cut
      logical :: m1
      real(kind=8),intent(out) :: a1
      real(kind=8),dimension(2) :: a1pm
      real(kind=8) :: a1maxbar,a1minbar,R,xy,a1max,a1min,map_den

      a1=0d0
      a1pm = a1_pm(a2,s,s1,s2,cos)
      if (bad_forward_jac()) return
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
      call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
      if (bad_forward_jac()) return
      xy = tan(-pi/2d0 * R)**2
      map_den=a1minbar+xy*a1maxbar
      if (abs(map_den).le.vtiny*max(1d0,abs(a1minbar),abs(xy*a1maxbar))) then
         ps%jac=-9d0
         return
      endif
      a1=a1maxbar*a1minbar*(1d0+xy)/map_den-h
      if (.not.clamp_haag_interval(a1,a1min,a1max,1d-8)) then
         ps%jac=-9d0
         return
      endif
      if (.not.scale_forward_jac(a1)) return
    end subroutine generate_a1_term2

    subroutine generate_a2_term1(i,m1,maxn,a1,s,s1,s2,cos,a2cut,h_in,a2)
      implicit none
      integer :: i,maxn
      real(kind=8) :: a1,s,s1,s2,cos,h,h_in,a2cut
      logical :: m1
      real(kind=8),intent(out) :: a2
      real(kind=8),dimension(2) :: a2pm
      real(kind=8) :: a2maxbar,a2minbar,R,xy,a2max,a2min,a2max_force,map_den,rad

      a2=0d0
      a2pm = a2_pm(a1,s,s1,s2,cos)
      if (bad_forward_jac()) return
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
         if (maxn-i-1.gt.0) then
            a2max_force = 1d0-a2cut/(maxn-i-1)
            if ((a2max_force.lt.a2max).and.(a2max_force.gt.a2min)) then
               a2maxbar = a2max_force + h
               a2max = a2max_force
            endif
         endif
      endif

      a2maxbar = a2max + h
      a2minbar = a2min + h

      ! Pretty much same as COMIX, except h stuff


      ! Now generate a2
      if ((.not. open) .and. (.not. this%flat_mode)) then
         if ((m1 .and. (i .ge. 2)) .or.&
              ((.not. m1) .and. (maxn .eq. this%next-2) .and. (i .ge. 1))&
              .or. (maxn .ne. this%next-2))  then
            ix = ix + 1
            call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
            if (bad_forward_jac()) return
            xy = tan(-pi/2d0 * R)**2
            map_den=a2minbar+xy*a2maxbar
            if (abs(map_den).le.vtiny*max(1d0,abs(a2minbar),abs(xy*a2maxbar))) then
               ps%jac=-9d0
               return
            endif
            a2=a2maxbar*a2minbar*(1d0+xy)/map_den-h
            if (.not.clamp_haag_interval(a2,a2min,a2max,1d-8)) then
               ps%jac=-9d0
               return
            endif
            if (.not.scale_forward_jac(a2)) return
            if (.not.scale_soft_weight(pi/2d0)) return
         elseif( ((i .eq. 0) .and. (maxn .eq. this%next-2)) .or. ((m1 .and. (i .le. 1)))) then
            ! Do phi-integration instead of a2
            a2 = 300d0          ! dummy value to do phi-integration
         endif
      endif

      if (this%flat_mode) then
         if ((m1 .and. (i .ge. 2)) .or.&
              ((.not. m1) .and. (maxn .eq. this%next-2) .and. (i .ge. 1))&
              .or. (maxn .ne. this%next-2))  then
            ix = ix +1
            call random_to_var(this%x(ix),0d0,a2min,a2max,a2,ps%jac)
            if (bad_forward_jac()) return
            !soft = soft*(a2max-a2min)
            rad=4d0*(a2max-a2)*(a2-a2min)
            if (rad.le.0d0) then
               ps%jac=-9d0
               return
            endif
            if (.not.divide_forward_jac(dsqrt(rad))) return
         elseif( ((i .eq. 0) .and. (maxn .eq. this%next-2)) .or. ((m1 .and. (i .le. 1)))) then
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
      real(kind=8) :: s,s1,s2,cos
      logical :: m1
      real(kind=8) :: a2cut
      real(kind=8) :: h1,a2min,a2max,f_h1,beta
      real(kind=8),dimension(3) :: buff
      real(kind=8) :: Amin,Amax,Atilde,R,wsq,kappa,root,den,map_den

      a2=0d0
      if (s.le.0d0 .or. beta.le.0d0) then
         ps%jac=-9d0
         return
      endif
      h1=0d0
      root=sqrt_kallen(1d0,s1/s,s2/s)
      if (bad_forward_jac()) return
      a2max=0.5d0*(1d0+(s2-s1)/s+root)
      a2min=0.5d0*(1d0+(s2-s1)/s-root)
      if ((a2cut.gt.a2min).and.(a2cut.lt.a2max)) then
         a2min=a2cut
      endif

      buff = f_func_term2(-h1,cos,s,s1,s2)
      if (bad_forward_jac()) return
      f_h1 = buff(2)
      if (f_h1.le.0d0) then
         ps%jac=-9d0
         return
      endif
      buff = f_func_term2(a2min,cos,s,s1,s2)
      if (bad_forward_jac()) return
      den=a2min+h1+buff(2)+f_h1
      if (abs(den).le.vtiny) then
         ps%jac=-9d0
         return
      endif
      Amin=(a2min+h1+buff(2)-f_h1)/den
      buff = f_func_term2(a2max,cos,s,s1,s2)
      if (bad_forward_jac()) return
      den=a2max+h1+buff(2)+f_h1
      if (abs(den).le.vtiny) then
         ps%jac=-9d0
         return
      endif
      Amax=(a2max+h1+buff(2)-f_h1)/den
      if (Amax.le.0d0) then
         ps%jac=-9d0
         return
      endif
      if (Amin.le.0d0) Amin=Amax*1d-8
      if (Amin.ge.Amax) then
         ps%jac=-9d0
         return
      endif

      ! Sample with Pi^{-1/2} factor
      ix = ix + 1
      call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
      if (bad_forward_jac()) return
      buff = f_func_term2(0d0,cos,s,s1,s2)
      if (bad_forward_jac()) return
      wsq = buff(3)  ! w^2
      Atilde = (Amin**(1d0-R))*(Amax)**R
      if (.not.ieee_is_finite(Atilde)) then
         ps%jac=-9d0
         return
      endif
      if (Atilde.le.0d0 .or. Atilde.ge.1d0-vtiny) then
         ps%jac=-9d0
         return
      endif
      kappa = - h1 + f_h1*(1d0+Atilde)/(1d0-Atilde)
      ! The linear coefficient is independent of a2, which is not yet set.
      buff = f_func_term2(0d0,cos,s,s1,s2)
      if (bad_forward_jac()) return
      map_den=buff(1)+2d0*kappa
      if (abs(map_den).le.vtiny) then
         ps%jac=-9d0
         return
      endif
      a2=(kappa**2-wsq)/map_den
      if (.not.clamp_haag_interval(a2,a2min,a2max,1d-8)) then
         ps%jac=-9d0
         return
      endif
      if (.not.scale_soft_weight(log(Amax)-log(Amin))) return
      if (.not.divide_soft_weight(f_h1)) return
      if (.not.scale_forward_jac(a2)) return
      if (.not.scale_forward_jac(beta)) return
    end subroutine generate_a2_term2

    function f_func_term1(a1,cos,s,s1,s2,h)
      implicit none
      real(kind=8),dimension(3) :: f_func_term1
      real(kind=8) :: beta, a1,cos,sin,s,s1,s2,lin_coeff,con_term,h,rad,den
      f_func_term1=0d0
      if (s.le.0d0 .or. 1d0-cos**2.lt.-1d-12) then
         ps%jac=-9d0
         return
      endif
      sin=dsqrt(max(1d0-cos**2,0d0))
      beta = (1d0+(s2-s1)/s) + (1d0-(s2-s1)/s)*cos
      if (h.gt.0d0) then ! Old approach
         if (cos.lt.0.999d0) then
            lin_coeff = -cos*beta-sin**2+sin**2*(s2-s1)/s-2d0*h*cos
         else
            ! use Taylor expansion around cos=1
            lin_coeff = -beta-2d0*h + 2d0*(1-beta-h-(s2-s1)/s)*(cos-1d0)
         endif
         con_term = 0.25d0*(beta**2) + h*beta + h**2 + (sin**2)*s1/s
      elseif (h.eq.0d0) then ! New approach
         den=2d0-2d0*cos
         if (den.le.vtiny) then
            ps%jac=-9d0
            return
         endif
         lin_coeff=(-cos*beta-sin**2+sin**2*(s2-s1)/s+beta)/den
         con_term=(0.25d0*beta**2+sin**2*s1/s)/den
      else
         lin_coeff = -cos*beta-sin**2+sin**2*(s2-s1)/s
         con_term = 0.25d0*(beta**2) + (sin**2)*s1/s
      endif
      rad=a1**2+lin_coeff*a1+con_term
      if (rad.lt.-1d-12*max(1d0,a1**2,abs(lin_coeff*a1),abs(con_term))) then
         ps%jac=-9d0
         return
      endif
      f_func_term1=(/lin_coeff,dsqrt(max(rad,0d0)),con_term/)
    end function f_func_term1

    function f_func_term2(a2,cos,s,s1,s2)
      implicit none
      real(kind=8),dimension(3) :: f_func_term2
      real(kind=8) :: a2,cos,sin,s,s1,s2,lin_coeff,con_term,beta,rad,den
      f_func_term2=0d0
      if (s.le.0d0 .or. 1d0-cos**2.lt.-1d-12) then
         ps%jac=-9d0
         return
      endif
      sin=dsqrt(max(1d0-cos**2,0d0))
      beta = (1d0-(s2-s1)/s)+(1d0+(s2-s1)/s)*cos
      den=2d0-2d0*cos
      if (den.le.vtiny) then
         ps%jac=-9d0
         return
      endif
      lin_coeff=(-cos*beta-sin**2*(1d0+(s2-s1)/s)+beta)/den
      con_term=(0.25d0*beta**2+sin**2*s2/s)/den
      rad=a2**2+lin_coeff*a2+con_term
      if (rad.lt.-1d-12*max(1d0,a2**2,abs(lin_coeff*a2),abs(con_term))) then
         ps%jac=-9d0
         return
      endif
      f_func_term2=(/lin_coeff,dsqrt(max(rad,0d0)),con_term/)
    end function f_func_term2

    function a2_pm(a1,s,s1,s2,c)
      ! cross-checked with pi function!
      implicit none
      real(kind=8) :: comm,s1,s2,s,a1,a2plus,a2minus,c,rad1,rad2
      real(kind=8), dimension(2) :: a2_pm
      a2_pm=0d0
      if (s.le.0d0) then
         ps%jac=-9d0
         return
      endif
      comm = 0.5d0*(1d0+((s2-s1)/s)+c*(1d0-2d0*a1-(s2-s1)/s))
      rad1=1d0-c**2
      rad2=a1*(1d0-a1-(s2-s1)/s)-s1/s
      if (rad1.lt.-1d-12 .or. rad2.lt.-1d-12) then
         ps%jac=-9d0
         return
      endif
      a2plus=comm+dsqrt(max(rad1,0d0))*dsqrt(max(rad2,0d0))
      a2minus=comm-dsqrt(max(rad1,0d0))*dsqrt(max(rad2,0d0))
      a2_pm = (/a2plus,a2minus/)
    end function a2_pm

    function a1_pm(a2,s,s1,s2,c)
      ! cross-checked with pi function!
      implicit none
      real(kind=8) :: comm,s1,s2,s,a2,a1plus,a1minus,c,rad1,rad2
      real(kind=8), dimension(2) :: a1_pm
      a1_pm=0d0
      if (s.le.0d0) then
         ps%jac=-9d0
         return
      endif
      comm = 0.5d0*(1d0-(s2-s1)/s+c*(1d0-2d0*a2+(s2-s1)/s))
      rad1=1d0-c**2
      rad2=a2*(1d0-a2+(s2-s1)/s)-s2/s
      if (rad1.lt.-1d-12 .or. rad2.lt.-1d-12) then
         ps%jac=-9d0
         return
      endif
      a1plus=comm+dsqrt(max(rad1,0d0))*dsqrt(max(rad2,0d0))
      a1minus=comm-dsqrt(max(rad1,0d0))*dsqrt(max(rad2,0d0))
      a1_pm = (/a1plus,a1minus/)
    end function a1_pm

    real(kind=8) function pi_func(a1,a2,s1,s2,s,cos)
      ! note: slightly different than in haag paper (same as in COMIX)
      implicit none
      real(kind=8) :: a1,a2,s1,s2,s,cos,sin
      sin=dsqrt(1d0-cos**2)
      pi_func = 4d0*sin**2*((1d0-a2+s2/s-s1/s)*a2 - s2/s)- &
           (1d0-2d0*a1+s1/s-s2/s + cos*(1d0-2d0*a2+s2/s-s1/s))**2
    end function pi_func

    function solver(s,s1,s2,q1_cmf,q2_cmf,z_sign,a1,a2,E,P_cmf,co)
      implicit none
      real(kind=8),dimension(3) :: solver
      real(kind=8) z_sign
      real(kind=8) :: s,s1,s2,a1,a2,E,co
      real(kind=8), dimension(0:3) :: q1_cmf,q2_cmf,P_cmf,q2_rot
      real(kind=8), dimension(0:3) :: p1,p1_rot,p1_xy_rot
      real(kind=8) :: R,xxx,y,z,phi,sgn,rad,sin2,phi_rotation
      logical,parameter :: old_mapping=.false.

      solver=0d0
      p1=0d0
      p1_rot=0d0
      if (a2 .gt. 100d0) then
         z = E - dsqrt(s)*a1
         ix = ix + 1
         call random_to_var(this%x(ix),0d0,0d0,2d0*pi,phi,ps%jac)
         if (bad_forward_jac()) return
         if (.not.scale_soft_weight(1d0/4d0)) return
         rad=E**2-s1-z**2
         if (rad.lt.-1d-10*max(spacing(0d0),E**2,s1,z**2)) then
            ps%jac=-11d0
            return
         endif
         rad=max(rad,0d0)
         xxx=dsqrt(rad)*cos(phi)
         y=dsqrt(rad)*sin(phi)
         p1=(/E,xxx,y,z/)
         call rotxxx(p1,q1_cmf,p1_rot)
      else
         z = E - dsqrt(s)*a1
         sin2=1d0-co**2
         if (sin2.le.vtiny) then
            ps%jac=-11d0
            return
         endif
         y=(-dsqrt(s)+E+dsqrt(s)*a2-co*z)/dsqrt(sin2)
         R=next_extra_random()
         if (bad_forward_jac()) return
         if (R .lt. 0.5d0) then
            sgn = 1d0
         else
            sgn=-1d0
         endif


         rad=E**2-s1-y**2-z**2
         if (rad.lt.-1d-10*max(spacing(0d0),E**2,s1,y**2,z**2)) then
            ps%jac=-11d0
            return
         endif
         xxx=sgn*dsqrt(max(rad,0d0))
         p1 = (/E,xxx,y,z/)

         ! p1 is expressed in the frame where q1 is +z and the transverse
         ! component of q2 is +y.  Obtain q2's actual azimuth in the
         ! q1-aligned frame, rotate around z, then rotate back to the CMF.
         call rotxxx_inv(q2_cmf,q1_cmf,q2_rot)
         if (q2_rot(1)**2+q2_rot(2)**2.le. &
              vtiny*max(spacing(0d0),q2_rot(0)**2)) then
            ps%jac=-11d0
            return
         endif
         phi_rotation=atan2(q2_rot(2),q2_rot(1))-pi/2d0
         p1_xy_rot(0)=p1(0)
         p1_xy_rot(1)=cos(phi_rotation)*p1(1)-sin(phi_rotation)*p1(2)
         p1_xy_rot(2)=sin(phi_rotation)*p1(1)+cos(phi_rotation)*p1(2)
         p1_xy_rot(3)=p1(3)
         call rotxxx(p1_xy_rot,q1_cmf,p1_rot)

      endif

      solver = (/p1_rot(1),p1_rot(2),p1_rot(3)/)
    end function solver

    function a1_m1(s,s1,s2,q1,P,k)
      implicit none
      real(kind=8) :: s,s1,s2,RHS
      real(kind=8) :: mu,a1max,a1min,sum_w,R,dum,sgn,a1cut,qdot,weight_den,map_den
      real(kind=8) :: arg_num_max,arg_den_max,arg_num_min,arg_den_min
      real(kind=8),dimension(4) :: g1,g2,d,e,w
      real(kind=8),dimension(0:3) :: q1,P
      real(kind=8) :: a1_m1
      integer :: i,pick,k
      real(kind=8) :: w1,w2,w_tot,a1

      a1_m1=0d0
      pick=0
      if (k.lt.2) then
         write(*,*) 'ERROR: invalid HAAG a1 channel:',k
         ps%jac=-9d0
         return
      endif
      dum=1d0
      qdot=dot(q1,P)
      if (s.le.0d0 .or. qdot.le.0d0) then
         ps%jac=-9d0
         return
      endif
      mu=(s2-s1)/s
      a1max = 0.5d0*(1d0+(s1-s2)/(s)+sqrt_kallen(1d0,s1/s,s2/s))
      a1min = 0.5d0*(1d0+(s1-s2)/(s)-sqrt_kallen(1d0,s1/s,s2/s))

      a1cut=(this%s0/2d0)/qdot

      if ((a1min.lt.a1cut ).and.(a1max.gt.a1cut)) then
         a1min = a1cut
      endif

      !if ((1d0-(k-1)*a1cut-mu.gt.a1min).and.(1d0-(k-1)*a1cut-mu.lt.a1max)) then
      !     a1max = 1d0-(k-1)*a1cut-mu
      !endif

      a1max = a1max-1d-8  !! This is a hard cut, but it works
      if (a1min.le.0d0 .or. a1min.ge.a1max) then
         ps%jac=-9d0
         return
      endif

      if (k.gt.2) then
         g1 = (/1d0,1d0,1d0,-1d0/)
         g2 = (/1d0,-1d0,-1d0,-1d0/)
         e = (/mu, 1d0-mu, 1d0, 1d0-mu/)
         d = (/0d0,0d0,mu,1d0/)

         do i=1,4
            weight_den=g2(i)*e(i)-g1(i)*d(i)
            arg_num_max=g1(i)*a1max+d(i)
            arg_den_max=g2(i)*a1max+e(i)
            arg_num_min=g1(i)*a1min+d(i)
            arg_den_min=g2(i)*a1min+e(i)
            if (abs(weight_den).le.vtiny .or. arg_num_max.le.0d0 .or. &
                 arg_den_max.le.0d0 .or. arg_num_min.le.0d0 .or. arg_den_min.le.0d0) then
               ps%jac=-9d0
               return
            endif
            w(i)=(log(arg_num_max/arg_den_max)-log(arg_num_min/arg_den_min))/weight_den
         enddo
         ! Correction 
         w(2) = -w(2)
         w(3) = -w(3)
         sum_w = sum(w)
         if (any(.not.ieee_is_finite(w)) .or. .not.ieee_is_finite(sum_w)) then
            ps%jac=-9d0
            return
         endif
         if (any(w.lt.0d0) .or. sum_w.le.0d0) then
            ps%jac=-9d0
            return
         endif

         R=next_extra_random()
         if (bad_forward_jac()) return
         if (R .lt. w(1)/sum_w) then
            pick = 1
         elseif (R .lt. (w(1)+w(2))/sum_w) then
            pick = 2
         elseif (R .lt. (w(1)+w(2)+w(3))/sum_w) then
            pick = 3
         else
            pick = 4
         endif
         ix = ix +1
         call random_to_var(this%x(ix),0d0,0d0,1d0,R,ps%jac)
         if (bad_forward_jac()) return
         sgn = 1d0
         if ((pick .eq. 2) .or. (pick .eq. 3)) then
            sgn = -1d0
         endif
         RHS = exp(w(pick)*R*sgn*(g2(pick)*e(pick)-g1(pick)*d(pick)))*&
              (g1(pick)*a1min+d(pick))/(g2(pick)*a1min+e(pick))
         map_den=g1(pick)-RHS*g2(pick)
         if (.not.ieee_is_finite(RHS) .or. .not.ieee_is_finite(map_den)) then
            ps%jac=-9d0
            return
         endif
         if (abs(map_den).le.vtiny) then
            ps%jac=-9d0
            return
         endif
         a1_m1=(e(pick)*RHS-d(pick))/map_den
         if (.not.clamp_haag_interval(a1_m1,a1min,a1max,1d-8)) then
            ps%jac=-9d0
            a1_m1=0d0
            return
         endif
         if (.not.scale_soft_weight(sum_w)) return

      elseif (k.eq.2) then
         a1min=(this%s0/2d0)/qdot
         a1max = 1d0
         if ((1d0-(this%s0/2d0)/qdot.lt.a1max).and.(1d0-(this%s0/2d0)/qdot.gt.a1min)) then
            a1max=1d0-(this%s0/2d0)/qdot
         endif
         if (a1min.le.0d0 .or. a1min.ge.a1max .or. a1max.ge.1d0) then
            ps%jac=-9d0
            return
         endif
         R=next_extra_random()
         if (bad_forward_jac()) return
         w1= log(a1max/a1min)
         w2 = log((1d0-a1min)/(1d0-a1max))
         w_tot=w1+w2
         if (.not.ieee_is_finite(w1) .or. .not.ieee_is_finite(w2) .or. &
              .not.ieee_is_finite(w_tot)) then
            ps%jac=-9d0
            return
         endif
         if (w1.lt.0d0 .or. w2.lt.0d0 .or. w_tot.le.0d0) then
            ps%jac=-9d0
            return
         endif

         if (R.lt.(w1/w_tot)) then
            pick=1
         else
            pick=2
         endif

         if (pick.eq.1) then
            ix = ix +1
            call random_to_var(this%x(ix),-1d0,a1min,a1max,a1,dum)
            if (dum.le.0d0) then
               ps%jac=-9d0
               return
            endif
            a1_m1 = a1
         elseif (pick.eq.2) then
            ix = ix +1
            call random_to_var(this%x(ix),-1d0,1d0-a1max,1d0-a1min,a1,dum)
            if (dum.le.0d0) then
               ps%jac=-9d0
               return
            endif
            ! This channel samples b=1-a with density 1/b.
            a1_m1 = 1d0-a1
         endif
         if (.not.clamp_haag_interval(a1_m1,a1min,a1max,1d-8)) then
            ps%jac=-9d0
            a1_m1=0d0
            return
         endif
         if (.not.scale_soft_weight(w_tot)) return

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

    real(kind=8) FUNCTION DDILOG(X)
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
         DO I = 19,0,-1
            B0=C(I)+ALFA*B1-B2
            B2=B1
            B1=B0
         enddo
         H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG=H
      RETURN
    end FUNCTION DDILOG
      subroutine check_momenta(p)
        implicit none
        real(kind=8), dimension(0:3,this%next) :: p
        real(kind=8), dimension(0:3) :: tot_mom
        integer i 
        real(kind=8) :: curr_mass

        tot_mom = (/0d0,0d0,0d0,0d0/)
        do i=1,this%next
           curr_mass = dot(p(0:3,i),p(0:3,i))
           if (abs(curr_mass - this%masses(i)**2) .gt. 1d-9) then
              write(*,*) 'ERROR in mass!',abs(curr_mass - this%masses(i)**2)
              write(*,*) abs(curr_mass - this%masses(i)**2) .gt. 1d-9
           endif
           tot_mom = tot_mom+p(:,i)
        enddo

        if (abs(tot_mom(0)-2d0*this%sqrtshat)/(2d0*this%sqrtshat) .gt.1d-4) then
           write(*,*) 'ERROR in energy conservation!'
        endif
        do i=1,3
           if (abs(tot_mom(i)).gt.1d-6) then
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
        real(kind=8),intent(in),dimension(0:3) :: p,q
        real(kind=8),intent(out),dimension(0:3) :: prot
        logical :: valid
        call stable_rotate_from_z_axis(p,q,prot,valid)
        if (.not.valid) prot=0d0
      end subroutine rotxxx

      subroutine rotxxx_inv(p,q,prot)
        ! Inverse of rotxxx: rotate a momentum from the frame containing q
        ! into the frame where q points along the positive z-axis.
        implicit none
        real(kind=8),dimension(0:3),intent(in) :: p,q
        real(kind=8),dimension(0:3),intent(out) :: prot
        logical :: valid
        call stable_rotate_to_z_axis(p,q,prot,valid)
        if (.not.valid) prot=0d0
      end subroutine rotxxx_inv

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
        real(kind=8),intent(in),dimension(0:3) :: p,q
        real(kind=8),intent(in) :: m
        real(kind=8),intent(out),dimension(0:3) :: pboost
        real(kind=8) :: pq,lf,boost_scale
        real(kind=8),parameter :: boost_tol=128d0*epsilon(1d0)
        pboost=0d0
        if (.not.all(ieee_is_finite(p)) .or. .not.all(ieee_is_finite(q)) .or. &
             .not.ieee_is_finite(m)) return
        boost_scale=max(spacing(0d0),abs(m),maxval(abs(q)))
        if (m.le.boost_tol*boost_scale) return
        if (q(0)+m.le.boost_tol*boost_scale) return
        pq=p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
        if (.not.ieee_is_finite(pq)) return
        lf=(pq/(q(0)+m)+p(0))/m
        pboost(0)=(p(0)*q(0)+pq)/m
        pboost(1:3)=p(1:3)+q(1:3)*lf
        if (.not.all(ieee_is_finite(pboost))) pboost=0d0
      end subroutine boostm

      subroutine random_to_var(x,power,var_min,var_max,var,jac)
        ! Given a random number x, it generates var in the range var_min
        ! <= var <= var_max according to var^(power). 'jac' is the
        ! corresponding Jacobian.
        implicit none
        real(kind=8),intent(in) :: x,power,var_min,var_max
      real(kind=8),intent(out) :: var
      real(kind=8),intent(inout) :: jac
      real(kind=8) :: varmin,varmax,xloc,jac_factor,interval_tolerance,updated_jac,scale
      logical :: valid_power_map
      var=var_min
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(power) .or. &
           .not.ieee_is_finite(var_min) .or. .not.ieee_is_finite(var_max) .or. &
           .not.ieee_is_finite(jac)) then
         jac=-10d0
         return
      endif
      if (jac.le.0d0) return
      if (x.lt.-vtiny .or. x.gt.1d0+vtiny .or. var_min.ge.var_max) then
         jac=-10d0
         return
      endif
      xloc=max(0d0,min(1d0,x))
      if (var_min.lt.0d0 .and. var_max.le.0d0) then
         varmin=-var_max
         varmax=-var_min
       elseif (var_min.lt.0d0 .and. var_max.gt.0d0 .and. (power.ne.0d0)) then
          jac=-10d0
          return
       else
          varmin=var_min
          varmax=var_max
       endif
       scale=max(abs(varmin),abs(varmax))
       interval_tolerance=128d0*epsilon(1d0)*max(spacing(0d0),scale)
       if (varmax-varmin.le.interval_tolerance .or. &
            (power.lt.0d0 .and. varmin.le.0d0)) then
          jac=-10d0
          return
       endif
       if (power.eq.0d0) then
          var=varmin+xloc*(varmax-varmin)
          jac_factor=varmax-varmin
       elseif (power.eq.-1d0) then
          jac_factor=log(varmax)-log(varmin)
          if (.not.ieee_is_finite(jac_factor)) then
             jac=-10d0
             return
          endif
          if (jac_factor.le.0d0) then
             jac=-10d0
             return
          endif
          var=exp((1d0-xloc)*log(varmin)+xloc*log(varmax))
          jac_factor=var*jac_factor
       else
          call stable_phase_space_power_map(xloc,power,varmin,varmax,var,&
               jac_factor,valid_power_map)
          if (.not.valid_power_map) then
             jac=-10d0
             return
          endif
       endif
       if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac_factor)) then
          var=var_min
          jac=-10d0
          return
       endif
       if (jac_factor.le.0d0) then
          var=var_min
          jac=-10d0
          return
       endif
       if (.not.safe_phase_space_product(jac,jac_factor,updated_jac)) then
          var=var_min
          jac=-10d0
          return
       endif
       jac=updated_jac
       if (var_min.le.0d0 .and. var_max.le.0d0) then
          var=-var
       endif
       if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac)) then
          var=var_min
          jac=-10d0
          return
       endif
       if (jac.le.0d0) then
          var=var_min
          jac=-10d0
       endif
      end subroutine random_to_var

      subroutine boostz(p,yb,pb)
        ! boost in the z-direction with rapidity yb
        implicit none
        real(kind=8),dimension(0:3),intent(in) :: p
        real(kind=8),intent(in) :: yb
        real(kind=8),dimension(0:3),intent(out) :: pb
        pb(0)=p(0)*cosh(yb)-p(3)*sinh(yb)
        pb(1:2)=p(1:2)
        pb(3)=p(3)*cosh(yb)-p(0)*sinh(yb)
      end subroutine boostz

    end subroutine haag_generate_momenta

  end module phase_space_haag_mod
