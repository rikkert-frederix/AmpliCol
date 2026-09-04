module phase_space_gen23_mod
  !  use common
  use phase_space_base
  use run_parameters, only: z_mass,z_width
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type,extends(phase_space_type),public :: phase_space_gen23
     ! Whether this phase-space instance uses the cut-aware (soft) bounds
     ! as its actual integration limits.  The module-level setter below is
     ! copied into this component during init so different channels can use
     ! different bound modes in one multichannel integration.
     logical :: use_soft_bounds_as_actual_limits=.false.
     logical :: include_pdf=.false.
     real(kind=8),dimension(-1:1) :: ip=[2d0,-1d0,-2d0]
     real(kind=8),dimension(-1:1) :: ip_shat=[2d0,-1.2d0,-2d0]
     real(kind=8),dimension(-1:1) :: ip_dt=[2d0,-1d0,-2d0]
     real(kind=8),dimension(-1:1) :: ip_mass=[2d0,-0.5d0,-2d0]
   contains
     procedure :: init => gen23_init
     procedure :: generate_momenta => gen23_generate_momenta
     procedure :: compute_x_from_momenta => gen23_compute_x_from_momenta
     procedure :: cleanup => gen23_cleanup
     final :: gen23_finalize
  end type phase_space_gen23
  private
  ! TECHNIAL PARAMETERS
  ! vebose:
  logical,parameter :: verbose=.true.
  logical,parameter,public :: debug=.false.
  ! importance sampling (0d0=flat transformation; -1d0=1/x transformation):
  real(kind=8),dimension(-1:1),parameter :: ip_flat=[0d0,0d0,0d0]
  ! tiny parameter cutoff to prevent/reduce numerical instabilities:
  real(kind=8),parameter :: vtiny=1d-12,tiny=1d-8
  real(kind=8),parameter :: tiny_kin=128d0*epsilon(1d0)
  real(kind=8),parameter :: pi=3.1415926535897932d0
  logical,parameter :: use_t_channel_at_start=.true.
  ! If true, the cut-aware bounds are used as the actual integration limits.
  logical :: use_soft_bounds_as_actual_limits=.false.

  public :: set_use_soft_bounds_as_actual_limits

contains
  subroutine set_use_soft_bounds_as_actual_limits(flag)
    implicit none
    logical,intent(in) :: flag
    use_soft_bounds_as_actual_limits=flag
  end subroutine set_use_soft_bounds_as_actual_limits

  pure elemental real(kind=8) function sqrt0(x)
    implicit none
    real(kind=8),intent(in) :: x
    sqrt0=sqrt(max(x,0d0))
  end function sqrt0

  subroutine longitudinal_rapidity(p,y,valid)
    implicit none
    real(kind=8),intent(in) :: p(0:3)
    real(kind=8),intent(out) :: y
    logical,intent(out) :: valid
    real(kind=8) :: lightcone_plus,lightcone_minus
    y=0d0
    lightcone_plus=p(0)+p(3)
    lightcone_minus=p(0)-p(3)
    valid=.false.
    if (.not.ieee_is_finite(lightcone_plus) .or. .not.ieee_is_finite(lightcone_minus)) return
    if (lightcone_plus.le.0d0 .or. lightcone_minus.le.0d0) return
    y=(log(lightcone_plus)-log(lightcone_minus))/2d0
    valid=ieee_is_finite(y)
  end subroutine longitudinal_rapidity

  subroutine gen23_init(this,sqrts,n,m,o,pt_cut,rap_cut,dr_cut,sqrt_s_min,t_chan,include_pdf,flat)
    ! Phase-space initialisation routines.
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    ! INPUT
    ! Sqrt(s-hat), i.e, the collision energy
    real(kind=8),intent(in) :: sqrts
    ! number of particles (initial state + final state)
    integer(kind=4),intent(in) :: n
    ! the colour order:
    integer(kind=4),dimension(n),intent(in) :: o
    ! rapidity (not used) and pT cut (and DR and sqrt_s_min) on all the particles
    real(kind=8),dimension(n),intent(in) :: rap_cut,pt_cut
    real(kind=8),dimension(n,n),intent(in) :: DR_cut,sqrt_s_min
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
    logical,intent(in),optional :: flat
    integer(kind=4) :: i,j,ndim_extra,cnt1,cnt2,allocation_status
    integer(kind=4),dimension(2) :: iset
    character(len=256) :: allocation_message
    if (n.lt.4) then
       write (*,*) 'ERROR in gen23_init() -- at least two final-state particles are required'
       stop 1
    endif
    if (n.gt.max_bitmask_particles) then
       write (*,*) 'ERROR in gen23_init() -- particle multiplicity exceeds bit-mask workspace limit:',&
            n,max_bitmask_particles
       stop 1
    endif
    if (.not.ieee_is_finite(sqrts)) then
       write (*,*) 'ERROR in gen23_init() -- invalid collider energy:',sqrts
       stop 1
    endif
    if (sqrts.le.0d0) then
       write (*,*) 'ERROR in gen23_init() -- invalid collider energy:',sqrts
       stop 1
    endif
    if (sqrts.gt.phase_space_momentum_limit) then
       write (*,*) 'ERROR in gen23_init() -- collider energy exceeds numerical range:',sqrts
       stop 1
    endif
    do i=1,n
       if (count(o.eq.i).ne.1) then
          write (*,*) 'ERROR in gen23_init() -- colour order is not a permutation:',o
          stop 1
       endif
       if (.not.ieee_is_finite(m(i))) then
          write (*,*) 'ERROR in gen23_init() -- invalid external mass:',i,m(i)
          stop 1
       endif
       if (m(i).lt.0d0) then
          write (*,*) 'ERROR in gen23_init() -- invalid external mass:',i,m(i)
          stop 1
       endif
    enddo
    if (any(.not.ieee_is_finite(pt_cut)) .or. any(.not.ieee_is_finite(rap_cut)) .or. &
         any(.not.ieee_is_finite(dr_cut)) .or. any(.not.ieee_is_finite(sqrt_s_min))) then
       write (*,*) 'ERROR in gen23_init() -- non-finite phase-space cut'
       stop 1
    endif
    if (max(maxval(abs(m)),maxval(abs(pt_cut)),maxval(abs(sqrt_s_min))).gt. &
         phase_space_momentum_limit) then
       write (*,*) 'ERROR in gen23_init() -- mass or dimensionful cut exceeds numerical range'
       stop 1
    endif
    call gen23_cleanup(this)
    this%use_soft_bounds_as_actual_limits=use_soft_bounds_as_actual_limits
    this%can_invert_momenta=.true.
    this%sqrtshat=sqrts
    this%sqrts=sqrts
    this%t_channel=t_chan
    if (verbose) then
       write (99,*) 'Setting up',n,'particle phase-space'
       write (99,*) 'Total available energy, sqrt(s-hat) =',this%sqrtshat
       write (99,*) 'Use the simple t-channel?',this%t_channel
    endif
    this%include_pdf=include_pdf
    this%next=n
    this%ndim=3*(this%next-2)-4
    if (this%include_pdf) this%ndim=this%ndim+2 ! the two Bjorken x's
    allocate(this%order(this%next),this%invm(maskr(this%next)),&
         this%invm_min(maskr(this%next),1:2),this%ETmin(maskr(this%next),1:2),&
         this%invm_max(maskr(this%next),1:2),this%pp(0:3,0:maskr(this%next)),&
         this%p(0:3,this%next),this%masses(this%next),this%x(this%ndim),&
         this%sets(0:this%next-2,2),this%ptcut(1:this%next),&
         this%drcut(maskr(this%next)),&
         this%sqrt_s_min(1:this%next,1:this%next),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       call gen23_cleanup(this)
       write (*,*) 'ERROR in gen23_init() -- unable to allocate phase-space workspace:',&
            trim(allocation_message)
       stop 1
    endif
    this%order=0
    this%invm=0d0
    this%invm_min=0d0
    this%ETmin=0d0
    this%invm_max=0d0
    this%pp(0:3,0:maskr(this%next))=0d0
    this%p=0d0
    this%masses=m
    this%x=0d0
    this%sets=0
    ! masses of external particles
    do i=1,n
       if ((i.eq.1 .or. i.eq.2) .and. m(i).ne.0d0) then
          write (*,*) 'ERROR in gen23_init() -- ', &
               & 'incoming particles should be massless'
          write (*,*) m
          stop 1
       endif
       this%invm(ibset(0,i-1))=m(i)**2
       this%invm(ibclr(maskr(this%next),i-1))=m(i)**2
    enddo
    
    if (verbose) write (99,*) 'masses:',m(1:n)
    this%drcut=0d0
    this%ptcut=0d0
    this%sqrt_s_min=0d0
    do i=1,this%next
       if (pt_cut(i).gt.0d0) then
          this%ptcut(i)=pt_cut(i)
       endif
       do j=1,this%next
          if (dr_cut(i,j).gt.0d0) then
             this%drcut(ibset(ibset(0,i-1),j-1))=dr_cut(i,j)
          endif
          if (sqrt_s_min(i,j).gt.0d0) then
             this%sqrt_s_min(i,j)=sqrt_s_min(i,j)
          endif
       enddo
    enddo
    call setup_PS_cuts(this)

    ! Bring the colour order to a canonical order (first in the list
    ! should be particle 1, i.e., the first incoming particle).
    do i=1,this%next
       if (o(i).eq.1) then
          do j=0,this%next-1
             this%order(j+1)=o(1+mod(i+j-1,this%next))
          enddo
          exit
       endif
    enddo

    ! Put in the two sets the particles that are between the two
    ! initial state ones, assuming cyclic permutation freedom
    ! (set(:,1) contains the ones between 1 and 2; set(:,2) contains
    ! the ones between 2 and 1)
    this%sets=0
    do i=2,this%next
       if (this%order(i).eq.2) then
          do j=i+1,this%next
             this%sets(0,2)=ibset(this%sets(0,2),this%order(j)-1)
          enddo
          this%sets(1:i-2,1)=this%order(2:i-1)
          this%sets(1:this%next-i,2)=this%order(i+1:this%next)
          exit
       endif
       this%sets(0,1)=ibset(this%sets(0,1),this%order(i)-1)
    enddo


!!$    ! Put in the two sets the particles that are between the two
!!$    ! initial state ones, assuming cyclic permutation freedom
!!$    ! (set(:,1) contains the ones between 1 and 2; set(:,2) contains
!!$    ! the ones between 2 and 1)
!!$
!!$    this%sets = 0
!!$    do i=2,this%next
!!$       if (this%order(i).eq.2) then
!!$          do j=i+1, this%next
!!$             if (this%order(j).eq.this%next) cycle
!!$             this%sets(0,2)=ibset(this%sets(0,2),this%order(j)-1)
!!$          enddo
!!$          cnt1 = count(this%order(2:i-1).ne.this%next)
!!$          if (cnt1.gt.0) then
!!$             this%sets(1:cnt1,1) = pack(this%order(2:i-1),this%order(2:i-1).ne.this%next)
!!$          endif
!!$          cnt2 = count(this%order(i+1:this%next).ne.this%next)
!!$          if (cnt2.gt.0) then
!!$             this%sets(1:cnt2,2) = pack(this%order(i+1:this%next), this%order(i+1:this%next).ne.this%next)
!!$          endif
!!$          if (ubound(this%sets,1).gt.cnt2) this%sets(cnt2+1:ubound(this%sets,1),2)=0
!!$          exit
!!$       endif
!!$       if (this%order(i).ne.this%next) then
!!$          this%sets(0,1) = ibset(this%sets(0,1),this%order(i)-1)
!!$       endif
!!$    enddo
    
    ndim_extra=0
    iset(1)=popcnt(this%sets(0,1))
    iset(2)=popcnt(this%sets(0,2))
    if (verbose) then
       write (99,*) "iset:",iset(1),iset(2),"ub:",ubound(this%sets,1)
       write (99,*) "set 1:"
       do i=0,min(iset(1),ubound(this%sets,1))
          write (99,*) i,this%sets(i,1)
       enddo
       write (99,*) "set 2:"
       do i=0,min(iset(2),ubound(this%sets,1))
          write (99,*) i,this%sets(i,2)
       enddo
    endif
    if (iset(1).eq.2 .and. iset(2).gt.0) ndim_extra=ndim_extra+1
    if (iset(2).eq.2 .and. iset(1).gt.0) ndim_extra=ndim_extra+1
    if (iset(1).ge.3) ndim_extra=ndim_extra+(iset(1)-2)
    if (iset(2).ge.3) ndim_extra=ndim_extra+(iset(2)-2)
    this%ndim_extra=ndim_extra
    if (verbose) then
       write (99,*) "Additional random numbers needed:",ndim_extra
    endif
    if (present(flat)) then
       if (flat) then
          this%ip=ip_flat
          this%ip_shat=ip_flat
          this%ip_dt=ip_flat
          this%ip_mass=ip_flat
       else
          this%ip=[2d0,-1d0,-2d0]
          this%ip_shat=[2d0,-1.2d0,-2d0]
          this%ip_dt=[2d0,-1d0,-2d0]
          this%ip_mass=[2d0,-0.5d0,-2d0]
       endif
    else
       this%ip=[2d0,-1d0,-2d0]
       this%ip_shat=[2d0,-1.2d0,-2d0]
       this%ip_dt=[2d0,-1d0,-2d0]
       this%ip_mass=[2d0,-0.5d0,-2d0]
    endif
    if (verbose) then
       write (99,*) "Power in importance sampling:",this%ip
    endif
  end subroutine gen23_init


  subroutine setup_PS_cuts(this)
    ! Given the input cuts, fills the minimum (s-channel) and/or
    ! maximum (t-channel) values the invariants can be in the
    ! phase-space generation. Does not apply these cuts on
    ! invariants not used in the phase-space generation.
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    real(kind=8) :: s_cut(2)
    real(kind=8) :: mass,cut
    integer(kind=4) :: i,j,k,npart
    this%invm_min=0d0
    this%invm_max=0d0
    do k=1,maskr(this%next)
       npart=popcnt(k)
       if (btest(k,0).and.btest(k,1)) then ! both initial state particles are part of 'k'
          this%invm_min(k,1:2)=0d0 ! no cuts
       elseif (btest(k,0).or.btest(k,1)) then ! one of the initial state particles is part of 'k'
          if (npart.eq.2) then ! exaclty two particles in 'k'
             do i=1,this%next
                if (.not.btest(k,i-1)) cycle ! particle 'i' is not in combined particle 'k'
                do j=1,this%next
                   if (.not.btest(k,j-1)) cycle ! particle 'j' is not in combined particle 'k'
                   this%invm_max(k,2)=-max(this%sqrt_s_min(i,j)**2,this%ptcut(i)**2,this%ptcut(j)**2)
                   this%invm_max(maskr(this%next)-k,2)=this%invm_max(k,2) ! all but the two particles
                enddo
             enddo
          endif
       else ! only final state particles in the combined particle 'k'
          ! total mass of external particles in 'k'
          mass=0d0
          do i=0,this%next-1
             if (.not.btest(k,i)) cycle ! particle 'i' is not in combined particle 'k'
             mass=mass+sqrt(this%invm(ibset(0,i)))
          enddo
          ! total from the cuts
          cut=0d0
          do i=1,this%next-1
             if (.not.btest(k,i-1)) cycle ! particle 'i' is not in combined particle 'k'
             do j=i+1,this%next
                if (.not.btest(k,j-1)) cycle ! particle 'j' is not in combined particle 'k'
                cut=cut+max(this%sqrt_s_min(i,j)**2,2d0*this%ptcut(i)*this%ptcut(j)* &
                     (1d0-cos(this%drcut(ibset(ibset(0,i-1),j-1)))))
             enddo
          enddo
          if (npart.eq.this%next-2) then ! all final state particles are in 'k'
             this%invm_min(k,1)=mass**2
             this%invm_min(k,2)=max(cut*dble(npart)/dble(npart-1),mass**2)
          else
             this%invm_min(k,1)=mass**2
             this%invm_min(k,2)=max(cut/2d0,mass**2)
          endif
       endif
    enddo
    call setup_ETmin(this)
    call apply_bound_mode(this)
  end subroutine setup_PS_cuts

  subroutine setup_ETmin(this)
    implicit none
    ! Setup the minimum required (transverse) energy for each (combination of)
    ! final state particle(s) in the collision c.o.m. frame. Based on the
    ! masses and the ptcut (i.e., assumes that all pz=0 and pT=pTcut)
    class(phase_space_gen23),intent(inout) :: this
    integer :: i,j
    this%ETmin=0d0
    do i=1,maskr(this%next)
       if (btest(i,0).or.btest(i,1)) cycle ! skip the ones that include incoming particles
       do j=0,this%next-1
          if (btest(i,j)) this%ETmin(i,1)=this%ETmin(i,1)+sqrt(this%invm(ibset(0,j)))
          if (btest(i,j)) this%ETmin(i,2)=this%ETmin(i,2)+sqrt(this%invm(ibset(0,j))+this%ptcut(j+1)**2)
       enddo
      this%ETmin(i,1:2)=max(this%ETmin(i,1:2),sqrt(this%invm_min(i,1:2)))
    enddo
  end subroutine setup_ETmin

  subroutine apply_bound_mode(this)
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    if (.not.this%use_soft_bounds_as_actual_limits) return
    this%invm_min(:,1)=this%invm_min(:,2)
    this%invm_max(:,1)=this%invm_max(:,2)
    this%ETmin(:,1)=this%ETmin(:,2)
  end subroutine apply_bound_mode

  subroutine select_integration_bounds(var_min,var_max,var_min_eff,var_max_eff,use_soft_bounds)
    implicit none
    real(kind=8),dimension(1:2),intent(in) :: var_min,var_max
    real(kind=8),dimension(1:2),intent(out) :: var_min_eff,var_max_eff
    logical,intent(in) :: use_soft_bounds
    var_min_eff=var_min
    var_max_eff=var_max
    if (.not.use_soft_bounds) return
    var_min_eff(1)=var_min(2)
    var_max_eff(1)=var_max(2)
    if (var_min_eff(1).gt.var_max_eff(1)) var_min_eff(1)=var_max_eff(1)
    var_min_eff(2)=var_min_eff(1)
    var_max_eff(2)=var_max_eff(1)
  end subroutine select_integration_bounds

  subroutine random_to_var_inputs(power_in,var_min,var_max,power,vmin,vmax,use_soft_bounds,valid)
    implicit none
    real(kind=8),dimension(-1:1),intent(in) :: power_in
    real(kind=8),dimension(1:2),intent(in) :: var_min,var_max
    real(kind=8),dimension(3),intent(out) :: power,vmin,vmax
    logical,intent(in) :: use_soft_bounds
    logical,intent(out) :: valid
    real(kind=8),dimension(1:2) :: varmin,varmax,var_min_loc,var_max_loc
    real(kind=8) :: range_scale,range_tolerance
    valid=.false.
    power=0d0
    vmin=0d0
    vmax=0d0
    if (any(.not.ieee_is_finite(power_in)) .or. any(.not.ieee_is_finite(var_min)) .or. &
         any(.not.ieee_is_finite(var_max))) return
    call select_integration_bounds(var_min,var_max,var_min_loc,var_max_loc, &
         use_soft_bounds)
    range_scale=max(maxval(abs(var_min_loc)),maxval(abs(var_max_loc)))
    range_tolerance=128d0*epsilon(1d0)*max(spacing(0d0),range_scale)
    if (var_min_loc(1).gt.var_max_loc(1)) return
    var_min_loc(2)=max(var_min_loc(2),var_min_loc(1))
    var_max_loc(2)=min(var_max_loc(2),var_max_loc(1))
    if (var_min_loc(2).gt.var_max_loc(2)) then
       var_min_loc(2)=var_max_loc(2)
    endif
    if (var_min_loc(1).lt.0d0 .and. var_max_loc(1).le.0d0) then
       power(1:3)=power_in(1:-1:-1)
       varmin=-var_max_loc
       varmax=-var_min_loc
    elseif (var_min_loc(1).ge.0d0 .and. var_max_loc(1).gt.0d0) then
       power(1:3)=power_in(-1:1)
       varmin=var_min_loc
       varmax=var_max_loc
    else
       power(1:3)=0d0
       varmin=var_min_loc
       varmax=var_max_loc
    endif
    if (varmin(1).gt.varmax(1)) return
    if (varmin(2).lt.varmin(1)) varmin(2)=varmin(1)
    if (varmax(2).gt.varmax(1)) varmax(2)=varmax(1)
    vmin(1)=varmin(1)
    if (varmin(2)-vmin(1).lt.range_tolerance) then
       vmax(1)=vmin(1)
    else
       vmax(1)=min(varmin(2),varmax(1))
    endif
    vmin(2)=vmax(1)
    if (varmax(2)-vmin(2).lt.range_tolerance) then
       vmax(2)=vmin(2)
    else
       vmax(2)=varmax(2)
    endif
    vmin(3)=vmax(2)
    if (varmax(1)-vmin(3).lt.range_tolerance) then
       vmax(3)=vmin(3)
    else
       vmax(3)=varmax(1)
    endif
    where (vmin(1:3).le.spacing(0d0) .and. power(1:3).le.-1d0)
       power(1:3)=0d0
    end where
    valid=.true.
  end subroutine random_to_var_inputs

  subroutine random_to_var_weights(power,vmin,vmax,q)
    implicit none
    real(kind=8),dimension(3),intent(in) :: power,vmin,vmax
    real(kind=8),dimension(3),intent(out) :: q
    integer(kind=4) :: i
    real(kind=8),dimension(3) :: log_inte,scaled_inte
    real(kind=8) :: total,max_log
    logical,dimension(3) :: active
    logical :: valid_integral

    q=0d0
    log_inte=-huge(1d0)
    active=.false.
    do i=1,3
       if (vmax(i).le.vmin(i)) cycle
       call log_power_integral(vmin(i),vmax(i),power(i),log_inte(i),valid_integral)
       if (.not.valid_integral) return
       active(i)=.true.
    enddo
    if (.not.any(active)) return

    ! Match the three power laws continuously at their shared boundaries.
    ! Work in logarithms: the former direct powers overflowed for perfectly
    ! representable invariant ranges and made the probabilities depend on the
    ! arbitrary energy unit used for those invariants.
    if (active(1) .and. power(2).ne.power(1)) then
       if (vmin(2).le.0d0) return
       log_inte(1)=log_inte(1)+(power(2)-power(1))*log(vmin(2))
    endif
    if (active(3) .and. power(2).ne.power(3)) then
       if (vmin(3).le.0d0) return
       log_inte(3)=log_inte(3)+(power(2)-power(3))*log(vmin(3))
    endif
    if (any(active .and. .not.ieee_is_finite(log_inte))) return
    max_log=maxval(log_inte,mask=active)
    scaled_inte=0d0
    where (active) scaled_inte=exp(log_inte-max_log)
    total=sum(scaled_inte)
    if (.not.ieee_is_finite(total)) return
    if (total.le.0d0) return
    q=scaled_inte/total
  end subroutine random_to_var_weights

  subroutine log_power_integral(lo,hi,p,log_integral,valid)
    ! Logarithm of integral_lo^hi v**p dv without forming a dimensionful
    ! power.  Scaling both bounds by c therefore adds (p+1)*log(c), as it
    ! should, while the calculation remains finite over the double range.
    implicit none
    real(kind=8),intent(in) :: lo,hi,p
    real(kind=8),intent(out) :: log_integral
    logical,intent(out) :: valid
    real(kind=8) :: exponent,log_ratio,difference

    valid=.false.
    log_integral=0d0
    if (.not.ieee_is_finite(lo) .or. .not.ieee_is_finite(hi) .or. &
         .not.ieee_is_finite(p)) return
    if (hi.le.lo) return
    if (p.eq.0d0) then
       difference=hi-lo
       if (.not.ieee_is_finite(difference)) return
       if (difference.le.0d0) return
       log_integral=log(difference)
       valid=ieee_is_finite(log_integral)
       return
    endif
    if (hi.le.0d0) return
    exponent=p+1d0
    if (abs(exponent).lt.1d-12) then
       if (lo.le.0d0) return
       difference=log(hi)-log(lo)
       if (.not.ieee_is_finite(difference)) return
       if (difference.le.0d0) return
       log_integral=log(difference)
    elseif (exponent.gt.0d0) then
       if (lo.lt.0d0) return
       if (lo.eq.0d0) then
          difference=1d0
       else
          log_ratio=log(lo)-log(hi)
          difference=one_minus_exp(exponent*log_ratio)
       endif
       if (.not.ieee_is_finite(difference)) return
       if (difference.le.0d0) return
       log_integral=exponent*log(hi)+log(difference)-log(exponent)
    else
       if (lo.le.0d0) return
       log_ratio=log(hi)-log(lo)
       difference=one_minus_exp(exponent*log_ratio)
       if (.not.ieee_is_finite(difference)) return
       if (difference.le.0d0) return
       log_integral=exponent*log(lo)+log(difference)-log(-exponent)
    endif
    valid=ieee_is_finite(log_integral)
  end subroutine log_power_integral

  pure real(kind=8) function one_minus_exp(x)
    ! Stable 1-exp(x) for the non-positive arguments used above.  Avoid
    ! relying on the non-standard EXPM1 extension offered by some compilers.
    implicit none
    real(kind=8),intent(in) :: x
    if (abs(x).lt.1d-5) then
       one_minus_exp=-x*(1d0+x*(0.5d0+x*(1d0/6d0+x*(1d0/24d0+x/120d0))))
    else
       one_minus_exp=1d0-exp(x)
    endif
  end function one_minus_exp


  
  subroutine gen23_finalize(this)
    type(phase_space_gen23) :: this
    call gen23_cleanup(this)
  end subroutine gen23_finalize

  subroutine gen23_cleanup(this)
    implicit none
    class(phase_space_gen23),intent(inout) :: this
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
    this%use_soft_bounds_as_actual_limits=.false.
  end subroutine gen23_cleanup

  subroutine gen23_generate_momenta(this,ps)
    ! Wrapper for the routine that generates the momenta.
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    type(psv),intent(inout) :: ps
    integer(kind=4) :: i,ix,ix_e
    real(kind=8) :: ycm,sqrtshat
    real(kind=8),parameter :: random_tolerance=4096d0*epsilon(1d0)
    real(kind=8),dimension(0:3,0:maskr(this%next)) :: pp
    real(kind=8),dimension(maskr(this%next)) :: invm
    ps%jac=-47d0
    if (.not.allocated(ps%p) .or. .not.allocated(ps%x)) return
    if (size(ps%p,1).ne.4 .or. size(ps%p,2).ne.this%next .or. &
         lbound(ps%p,1).ne.0 .or. &
         size(ps%x).lt.this%ndim+this%ndim_extra) return
    ps%p=0d0
    ps%xbjrk=1d0
    if (any(.not.ieee_is_finite(ps%x(1:this%ndim+this%ndim_extra)))) return
    if (any(ps%x(1:this%ndim+this%ndim_extra).lt.-random_tolerance) .or. &
         any(ps%x(1:this%ndim+this%ndim_extra).gt.1d0+random_tolerance)) return
    ps%x(1:this%ndim+this%ndim_extra)=max(0d0, &
         min(1d0,ps%x(1:this%ndim+this%ndim_extra)))
    pp=0d0
    ps%jac=1d0
    ix=0
    ix_e=this%ndim

!!$       call generate_bw_mass(this%next-1)
!!$       this%next=this%next-1
!!$       call setup_PS_cuts(this)

    invm=this%invm
    if (this%include_pdf) then
       call generate_initial_state
       if (bad_forward_jac()) return
    else
       sqrtshat=this%sqrts
    endif
    call generate_momenta
    if (bad_forward_jac()) return
    if (ix.ne.this%ndim .or. ix_e.gt.size(ps%x) .or. .not.ieee_is_finite(ps%jac)) then
       ps%jac=-47d0
       return
    endif
!!$       this%next=this%next+1
!!$    if (ps%jac.lt.0d0) return
!!$
!!$       call decay_bw(this%next-1,this%next-1,this%next)
!!$       ! Add factors of 2*pi (since this%next was reduced by 1)
!!$       ps%jac=ps%jac/((2d0*pi)**3)
       if (debug) call test_momenta

    do i=1,this%next
       if (this%include_pdf) then
          ! Note: 'ycm' is the rapidity needed to go from lab to CM
          ! frame. Hence, here we boost from CM to lab frame with '-ycm'
          call boostz(pp(0:3,ibset(0,i-1)),-ycm,ps%p(0:3,i))
       else
          ps%p(0:3,i)=pp(0:3,ibset(0,i-1))
       endif
    enddo
    if (.not.generated_momenta_are_valid(ps%p,this%masses,ps%xbjrk,this%include_pdf)) then
       ps%p=0d0
       ps%jac=-48d0
    endif
  contains
    logical function bad_forward_jac()
      implicit none
      if (.not.ieee_is_finite(ps%jac)) then
         ps%jac=-47d0
         bad_forward_jac=.true.
         return
      endif
      bad_forward_jac=ps%jac.le.0d0
    end function bad_forward_jac

    logical function update_forward_jac(numerator,denominator)
      implicit none
      real(kind=8),intent(in) :: numerator,denominator
      real(kind=8) :: updated_jac
      update_forward_jac=.false.
      if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator) .or. &
           numerator.le.0d0 .or. denominator.le.0d0) then
         ps%jac=-47d0
         return
      endif
      if (.not.safe_phase_space_scaled_ratio(ps%jac,numerator,denominator,updated_jac)) then
         ps%jac=-47d0
         return
      endif
      ps%jac=updated_jac
      update_forward_jac=.true.
    end function update_forward_jac

    subroutine generate_bw_mass(ires)
      implicit none
      integer,intent(in) :: ires
      real(kind=8) :: smin,smax,qmass,qwidth,y
      real(kind=8),dimension(1:2) :: A,B
      smin=50d0**2
      smax=this%sqrts**2
      qmass=z_mass
      qwidth=z_width
      A=atan((qmass-smin/qmass)/qwidth)
      B=atan((qmass-smax/qmass)/qwidth)
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,B,A,y,ps%jac)
      if (bad_forward_jac()) return
      this%invm(ibset(0,ires-1))=qmass*(qmass-qwidth*tan(y))
      if (.not.update_forward_jac(qmass*qwidth,cos(y)**2)) return
    end subroutine generate_bw_mass
    subroutine decay_bw(ires,id1,id2)
      ! decay the BW
      implicit none
      integer,intent(in) :: ires,id1,id2
      integer :: i,j,k
      real(kind=8),dimension(0:3,this%next) :: pp_tmp
      if (ires.gt.id1 .or. ires.gt.id2) then
         write (*,*) 'ERROR in decay_BW',ires,id1,id2
      endif
      ! shift the already generated momenta to make place for the
      ! daughters of the resonance
      j=0
      k=0
      do i=1,this%next-2
         j=j+1
         k=k+1
         if (k.eq.id1) k=k+1
         if (k.eq.id2) k=k+1
         if (i.eq.ires) j=j+1 ! skip ires
         pp_tmp(0:3,k)=pp(0:3,ibset(0,j-1))
      enddo
      pp(0:3,ibset(0,id1-1)+ibset(0,id2-1))=pp(0:3,ibset(0,ires-1))
      invm(ibset(0,id1-1)+ibset(0,id2-1))=invm(ibset(0,ires-1))
      do i=1,this%next
         pp(0:3,ibset(0,i-1))=pp_tmp(0:3,i)
      enddo
      ! Set the decay products of the BW as massless particles
      invm(ibset(0,id1-1))=0d0
      invm(ibset(0,id2-1))=0d0
      call gens_one_step(ibset(0,id1-1),ibset(0,id2-1))
    end subroutine decay_bw
    subroutine generate_initial_state
      implicit none
      real(kind=8) :: tau
      call generate_tau(tau)
      if (bad_forward_jac()) return
      call generate_y(tau)
      if (bad_forward_jac()) return
      sqrtshat=sqrt(tau)*this%sqrts
      ps%xbjrk(1)=sqrt(tau)*exp(ycm)
      ps%xbjrk(2)=sqrt(tau)*exp(-ycm)
      if (debug) write (*,*) 'sqrtshat :',sqrtshat,ps%xbjrk(1:2),sqrtshat**2
    end subroutine generate_initial_state

    subroutine generate_tau(tau)
      implicit none
      real(kind=8),intent(out) :: tau
      real(kind=8) :: shat
      real(kind=8),dimension(1:2) :: smin,smax
      smin(1:2)=max(this%invm_min(maskr(this%next)-3,1:2),this%ETmin(maskr(this%next)-3,1:2)**2)
      smax(1:2)=this%sqrts**2
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip_shat,smin,smax,shat,ps%jac)
      if (bad_forward_jac()) return
      tau=shat/smax(1)
      if (.not.update_forward_jac(1d0,smax(1))) return
    end subroutine generate_tau

    subroutine generate_y(tau)
      implicit none
      real(kind=8),intent(in) :: tau
      real(kind=8),dimension(1:2) ::  ymin,ymax
      ymin(1:2)= log(tau)/2d0
      ymax(1:2)=-log(tau)/2d0
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,ymin,ymax,ycm,ps%jac)
      if (bad_forward_jac()) return
    end subroutine generate_y

    subroutine generate_momenta
      implicit none
      integer(kind=4) :: i,j,inext,im1
      integer(kind=4),dimension(2) :: set
      ! incoming momenta
      pp(0,ibset(0,0))=sqrtshat/2d0
      pp(0,ibset(0,1))=sqrtshat/2d0
      pp(1:2,ibset(0,0))=0d0
      pp(1:2,ibset(0,1))=0d0
      pp(3,ibset(0,0))= pp(0,ibset(0,0))
      pp(3,ibset(0,1))=-pp(0,ibset(0,1))
      ! incoming momenta when labeled from the other side
      pp(0:3,maskr(this%next)-ibset(0,0))=pp(0:3,1)
      pp(0:3,maskr(this%next)-ibset(0,1))=pp(0:3,2)
      ! invariant mass of all final state particles combined
      invm(ibset(0,0)+ibset(0,1))=sqrtshat**2
      invm(maskr(this%next)-ibset(0,0)-ibset(0,1))=sqrtshat**2
      set(1)=this%sets(0,1)
      set(2)=this%sets(0,2)
      ! Generate the central 2->2 process in case both set(1) and set(2) are not empty
      if (popcnt(set(1)).gt.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'two sets with at least two ',&
              & 'particles',popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         if (use_t_channel_at_start) then
            call gent_one_step(set(2),set(1),1)
         else
            call gens_one_step(set(2),set(1))
         endif
         if (bad_forward_jac()) return
         pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'special double t-channel (1)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t(set(1),set(2),1,2)
         if (bad_forward_jac()) return
         pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).gt.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) 'special double t-channel (2)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t(set(2),set(1),1,2)
         if (bad_forward_jac()) return
         pp(0:3,set(1)+1)=pp(0:3,2)-pp(0:3,set(2))
         invm(set(1)+1)=dot(pp(0:3,set(1)+1),pp(0:3,set(1)+1))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) '2->2 scattering with one particle in each set'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call gent_one_step(set(2),set(1),1)
         if (bad_forward_jac()) return
         pp(0:3,set(2)+2)=pp(0:3,1)-pp(0:3,set(1))
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      endif

      do i=1,2
         if (popcnt(set(i)).le.1) cycle ! at least 2 particles in a set
         inext=ibset(0,this%sets(1,i)-1)
         set(i)=set(i)-inext
         if (popcnt(set(i)).ge.2) then
            ! at least 3 particles in a set
            if (debug) write (*,*) 'At least 3 particles in a set',&
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step(set(i),inext,i)
            if (bad_forward_jac()) return
            pp(0:3,(3-i)+set(i)+inext)=pp(0:3,i)-pp(0:3,this%sets(0,(3-i)))
            invm((3-i)+set(i)+inext)=dot(pp(0:3,(3-i)+set(i)+inext),pp(0:3,(3-i)+set(i)+inext))
            pp(0:3,(3-i)+set(i))=pp(0:3,i)-pp(0:3,this%sets(0,(3-i)))-pp(0:3,inext)
            invm((3-i)+set(i))=dot(pp(0:3,(3-i)+set(i)),pp(0:3,(3-i)+set(i)))
            do j=2,popcnt(set(i))-1
               ! loop over the remaining particles in the set
               inext=ibset(0,this%sets(j,i)-1)
               im1=ibset(0,this%sets(j-1,i)-1)
               set(i)=set(i)-inext
               if (this%t_channel) then
                  call gent_one_step(inext,set(i),3-i)
               else
                  call gen23_one_step(inext,set(i),3-i,im1)
               endif
               if (bad_forward_jac()) return
            enddo
            inext=ibset(0,this%sets(j,i)-1)
            im1=ibset(0,this%sets(j-1,i)-1)
            set(i)=set(i)-inext
            if (this%t_channel) then
               call gent_one_step(inext,set(i),3-i)
            else
               call gen23_one_step(inext,set(i),3-i,im1)
            endif
            if (bad_forward_jac()) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).ne.0) then
            ! Exactly 2 particles in a set (and the other set contains at least one)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains at least one)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            im1=3-i
            pp(0:3,set(i)+inext+im1)=pp(0:3,set(i)+inext)-pp(0:3,im1)
            pp(0:3,set(i)+inext+i+im1)=pp(0:3,set(i)+inext+im1)-pp(0:3,i)
            invm(set(i)+inext+im1)=dot(pp(0:3,set(i)+inext+im1),pp(0:3,set(i)+inext+im1))
            invm(set(i)+inext+i+im1)=dot(pp(0:3,set(i)+inext+im1+i),pp(0:3,set(i)+inext+im1+i))
            if (this%t_channel) then
               call gent_one_step(set(i),inext,i)
            else
               call gen23_one_step(set(i),inext,i,im1)
            endif
            if (bad_forward_jac()) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).eq.0) then
            ! Exactly 2 particles in a set (and the other set contains none)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains none)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step(set(i),inext,i)
            if (bad_forward_jac()) return
         else
            write (*,*) 'Inconsistent sets'
            write (*,*) i,':',this%sets(:,i)
            write (*,*) 3-i,':',this%sets(:,3-i)
            ps%jac=-47d0
            return
         endif
         ! We need to get the momentum of the final particle of the set.
         pp(0:3,set(i))=pp(0:3,set(i)+inext+(3-i))+pp(0:3,(3-i))-pp(0:3,inext)
      enddo
      ! Add factors of 2*pi
      if (.not.update_forward_jac(1d0,(2d0*pi)**(3*(this%next-2)-4))) return
      ! Add flux factor
      if (.not.update_forward_jac(1d0,2d0*sqrtshat**2)) return
    end subroutine generate_momenta

    subroutine test_momenta
      ! Writes the momenta to the screen.
      implicit none
      integer(kind=4) :: i,i1,i2,i1b,i2b
      real(kind=8),dimension(0:3) :: ptot
      write (*,*) 'Momenta check:'
      if (ix.ne.this%ndim) then
         write (*,*) 'ERROR: number of random numbers used not consistent',ix,this%ndim
         stop 1
      endif
      ptot(0:3)=0d0
      do i=0,this%next-1
         ptot(0:3)=ptot(0:3)+pp(0:3,ibset(0,i))
         write (*,*) i+1,pp(0:3,ibset(0,i)),dot(pp(0,ibset(0,i)),pp(0,ibset(0,i)))
      enddo
      write (*,*) 'ptot',ptot(0:3)
      do i=1,this%next
         i1=this%order(i)
         i2=this%order(mod(i,this%next)+1)
         i1b=ibset(0,i1-1)
         i2b=ibset(0,i2-1)
         write (*,*) '(pi+pj)^2',i1,i2,invm(i1b+i2b),invm(maskr(this%next)-(i1b+i2b))&
              &,dot(pp(0:3,i1b)+pp(0:3,i2b),pp(0:3,i1b)+pp(0:3,i2b))
      enddo
      write (*,*) 'number of dimensions',this%ndim
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
      real(kind=8) :: phi,pt2,mass_tol,den_t,den_ir,den_scale,jac_den,jac_scale
      real(kind=8),dimension(1:2) :: tmin,tmax,pzmax,Eimax,yr
      logical :: soft_ok,massless_collinear_limit
      if (popcnt(i).ne.1 .or. popcnt(ir).le.1) then
         write (*,*) 'Subroutine only for i is a single particle '&
              //'and ir is more than 1',i,ir,popcnt(i),popcnt(ir)
         ps%jac=-47d0
         return
      endif
      soft_ok=.true.
      yr(1:2)=lambda(invm(ia+ib),invm(i),this%invm_min(ir,1:2))
      if (yr(1).lt.0d0) then
         ps%jac=-1d0
         return
      endif
      if (yr(2).lt.0d0) soft_ok=.false.
      yr(1)=sqrt(yr(1))
      if (soft_ok) yr(2)=sqrt(yr(2))
      tmin(1)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,1)-yr(1))/2d0
      tmax(1)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,1)+yr(1))/2d0
      if (soft_ok) then
         tmin(2)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,2)-yr(2))/2d0
         tmax(2)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,2)+yr(2))/2d0
      else
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      endif
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2) = min(tmax(1:2), this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin(1:2),this%invm_min(ir+ib,1:2))
      end where
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      pzmax(1:2)=lambda(sqrtshat**2,this%ETmin(i,1:2)**2,this%ETmin(ir,1:2)**2)
      if (pzmax(1).lt.0d0) then
         ps%jac=-1d0
         return
      endif
      if (pzmax(2).lt.0d0) soft_ok=.false.
      pzmax(1)=sqrt(pzmax(1))/(2d0*sqrtshat)
      Eimax(1)=sqrtshat-sqrt(this%ETmin(ir,1)**2+pzmax(1)**2)
      tmin(1)=max(tmin(1),invm(i)-sqrtshat*(Eimax(1)+pzmax(1)))
      tmax(1)=min(tmax(1),invm(i)-sqrtshat*(Eimax(1)-pzmax(1)))
      if (soft_ok) then
         pzmax(2)=sqrt(pzmax(2))/(2d0*sqrtshat)
         Eimax(2)=sqrtshat-sqrt(this%ETmin(ir,2)**2+pzmax(2)**2)
         tmin(2)=max(tmin(2),invm(i)-sqrtshat*(Eimax(2)+pzmax(2)))
         tmax(2)=min(tmax(2),invm(i)-sqrtshat*(Eimax(2)-pzmax(2)))
      else
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-1d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip_dt,tmin,tmax,invm(i+ia),ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) 'dt- i+ia',i+ia,invm(i+ia),tmin,tmax
      endif
      den_t=invm(i)-invm(i+ia)
      den_ir=sqrtshat**2+invm(i+ia)-invm(i)
      den_scale=max(spacing(0d0),abs(invm(i)),abs(invm(i+ia)),sqrtshat**2)
      massless_collinear_limit=abs(den_t).le.vtiny*den_scale .and. &
           abs(invm(i)).le.vtiny*den_scale .and. &
           this%ETmin(i,1).le.vtiny*sqrt(den_scale)
      if (abs(den_ir).le.vtiny*den_scale .or. &
           (abs(den_t).le.vtiny*den_scale .and. .not.massless_collinear_limit)) then
         ps%jac=-2d0
         if (debug) write (*,*) 'singular double_t constraint',den_t,den_ir
         return
      endif
      tmin(1:2)=-invm(ia+ib)-invm(i+ia)+invm(i)+this%invm_min(ir,1:2)
      if (massless_collinear_limit) then
         tmax(1:2)=0d0
      else
         tmax(1:2)=invm(i)*(invm(i)-invm(ia+ib)-invm(i+ia))/den_t
      endif
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2) = min(tmax(1:2), this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin(1:2),this%invm_min(ir+ib,1:2))
      end where
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      tmin(1:2)=max(tmin(1:2),invm(i)-sqrtshat**2*(1-this%ETmin(ir,1:2)**2/den_ir))
      if (.not.massless_collinear_limit) then
         tmax(1:2)=min(tmax(1:2),invm(i)-sqrtshat**2*(this%ETmin(i,1:2)**2/den_t))
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-2d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip_dt,tmin,tmax,invm(i+ib),ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) 'dt- i+ib',i+ib,invm(i+ib),tmin,tmax
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],phi,ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) 'dt- phi',phi
      endif
      pt2=invm(i+ia)*invm(i+ib)/invm(ia+ib)+ &
           & invm(i)**2/invm(ia+ib)-(invm(i+ia)+invm(i+ib))*invm(i)/invm(ia+ib)-invm(i)
      if (pt2.lt.-vtiny*max(spacing(0d0),abs(invm(ia+ib)))) then
         ps%jac=-3d0
         if (debug) write (*,*) "rejecting double_t point with negative pt2",pt2
         return
      elseif (pt2.lt.0d0) then
         pt2=0d0
      endif
      pp(0,i)=(-invm(i+ia)-invm(i+ib)+2d0*invm(i))/(2d0*sqrtshat)
      pp(1,i)=sqrt(pt2)*cos(phi)
      pp(2,i)=sqrt(pt2)*sin(phi)
      pp(3,i)=(invm(i+ia)-invm(i+ib))/(2d0*sqrtshat)
      pp(0,ir)=sqrtshat-pp(0,i)
      pp(1:3,ir)=-pp(1:3,i)
      invm(ir)=dot(pp(0:3,ir),pp(0:3,ir))
      mass_tol=vtiny*max(spacing(0d0),abs(invm(ia+ib)),abs(invm(i+ia)), &
           abs(invm(i+ib)),abs(invm(i)))
      if (invm(ir).lt.-mass_tol) then
         ps%jac=-4d0
         if (debug) write (*,*) "rejecting double_t point with negative remainder mass",ir,invm(ir),i
         return
      elseif (invm(ir).lt.0d0) then
         invm(ir)=0d0
      endif
      jac_den=sqrt0(lambda(invm(ir+i),0d0,0d0))
      jac_scale=max(spacing(0d0),abs(invm(ir+i)))
      if (jac_den.le.tiny_kin*jac_scale) then
         ps%jac=-4d0
         return
      endif
      if (.not.update_forward_jac(1d0,4d0*jac_den)) return
    end subroutine double_t

    subroutine gen23_one_step(i,ir,ib,im1)
      ! Generates one step using the 2->3 setup, using the invariants
      ! shat(i), s(i) and t(i) as defined in E.~Byckling and K.~Kajantie,
      ! ``Reductions of the phase-space integral in terms of simpler
      ! processes,'' Phys. Rev. 187 (1969), 2008-2016,
      ! doi:10.1103/PhysRev.187.2008.  Assumes massless incoming particles.
      implicit none
      integer(kind=4),intent(in) :: im1,i,ir,ib
      integer :: gent_status
      real(kind=8) :: tmin_S,tmax_S,smin_S,smax_S,phi1,phi2,gram4,V,sqrtGG,y,phi_rot,delta_ir,delta_scale
      real(kind=8),dimension(1:2) :: shatmin,shatmax,tmin,tmax,etminir,etmini,base,root,smin,smax,rad
      real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
      logical :: rapidity_ok
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         call generate_mass(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         call generate_mass(ir,shatmin,shatmax)
      endif
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) '23- i    ',i,invm(i)
         write (*,*) '23- ir   ',ir,invm(ir)
      endif
      call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin_S,tmax_S)
      tmin(1:2)=tmin_S
      tmax(1:2)=tmax_S
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2)=min(tmax_S,this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin_S,this%invm_min(ir+ib,1:2))
      end where
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
      call longitudinal_rapidity(pp(0:3,i+ir),y,rapidity_ok)
      if (.not.rapidity_ok) then
         ps%jac=-35d0
         return
      endif
      call boostz(pp(0,i+ir),y,piir)
      call boostz(pp(0,ib),y,pib)
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(i,1:2)**2-invm(i) .and. popcnt(i).eq.1 )
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(i,1:2)**2-invm(i)) )**2))

      elsewhere
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)))
      end where
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir,1:2)**2-invm(ir) .and. popcnt(ir).eq.1 )
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(ir,1:2)**2-invm(ir)) )**2))
      elsewhere
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)))
      end where
      base(1:2)=piir(0)**2-ETmini(1:2)**2+ETminir(1:2)**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root(1:2)=(piir(0)-ETmini(1:2)-ETminir(1:2))*(piir(0)+ETmini(1:2)-ETminir(1:2))*&
           (piir(0)-ETmini(1:2)+ETminir(1:2))*(piir(0)+ETmini(1:2)+ETminir(1:2))
      if (root(1).lt.0d0) then
         ps%jac=-33d0
         if (debug) write (*,*) 'root.lt.0d0',root
         return
      endif
      if (abs(piir(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
         ps%jac=-38d0
         if (debug) write (*,*) 'zero-energy system in gen23_one_step',piir(0)
         return
      endif
      tmin(1)=max(tmin(1),invm(ir)-pib(0)/piir(0)*(base(1)+sqrt0(root(1))))
      tmax(1)=min(tmax(1),invm(ir)-pib(0)/piir(0)*(base(1)-sqrt0(root(1))))
      if (root(2).lt.0d0) then
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      else
         tmin(2)=max(tmin(2),invm(ir)-pib(0)/piir(0)*(base(2)+sqrt0(root(2))))
         tmax(2)=min(tmax(2),invm(ir)-pib(0)/piir(0)*(base(2)-sqrt0(root(2))))
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-3d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip,tmin,tmax,invm(ir+ib),ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) '23- ir+ib',ir+ib,invm(ir+ib),tmin,tmax
      endif
      call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
           &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin_S,smax_S,V,sqrtGG)
      smin(1:2)=smin_S
      smax(1:2)=smax_S
      where (this%invm_min(i+im1,1:2).ne.0d0)
         smin(1:2)=max(smin_S,this%invm_min(i+im1,1:2))
      end where
      where (this%invm_max(i+im1,1:2).ne.0d0)
         smax(1:2)=min(smax_S,this%invm_max(i+im1,1:2))
      end where
      if (im1.gt.2) then
         ! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
         call longitudinal_rapidity(pp(0:3,im1),y,rapidity_ok)
         if (.not.rapidity_ok) then
            ps%jac=-35d0
            return
         endif
         call boostz(pp(0,i+ir),y,piirr)
         call boostz(pp(0,im1),y,pim1r)
         call boostz(pp(0,ib),y,pib)
         phi_rot=atan2(pp(2,im1),pp(1,im1))
         call rotz(piirr,-phi_rot,piir)
         call rotz(pim1r,-phi_rot,pim1)
         ! Eir > Etmin(ir) + constraint coming from t
         delta_ir=invm(ir)-invm(ir+ib)
         delta_scale=max(spacing(0d0),abs(invm(ir)),abs(invm(ir+ib)),pib(0)**2)
         if (abs(pib(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
            ps%jac=-34d0
            if (debug) write (*,*) 'unphysical Eir constraint',delta_ir,pib(0)
            return
         endif
         if (abs(delta_ir).le.vtiny*delta_scale) then
            ! For delta->0- the max has the finite limit ETmin(ir).  The
            ! same holds for a zero transverse bound from either side.
            if (delta_ir.le.0d0 .or. this%ETmin(ir,1).le.vtiny*sqrt(delta_scale)) then
               etminir(1:2)=this%ETmin(ir,1:2)
            else
               ps%jac=-34d0
               if (debug) write (*,*) 'positive singular Eir constraint',delta_ir,pib(0)
               return
            endif
         else
            etminir(1:2)=max(pib(0)*this%ETmin(ir,1:2)**2/delta_ir+delta_ir/(4d0*pib(0)),&
                 this%ETmin(ir,1:2))
         endif
         rad=(piir(0)-etminir(1:2))**2-invm(i)
         if (rad(1).lt.0d0) then
            ps%jac=-34d0
            if (debug) write (*,*) 'rad.lt.0d0 inverse',rad,piir(0),etminir,invm(i)
            return
         endif
         ! The second bound is an optional soft/pT importance-sampling
         ! sector.  It can be empty for this sequential construction even
         ! though the physical (hard) sector, stored in entry one, is valid.
         ! Do not reject the physical point in that case.
         smax(1)=min(smax(1),invm(i)+invm(im1)+2d0*(piir(0)-etminir(1))*pim1(0)+&
              2d0*sqrt(rad(1))*pim1(1))
         if (rad(2).ge.0d0) then
            smax(2)=min(smax(2),invm(i)+invm(im1)+2d0*(piir(0)-etminir(2))*pim1(0)+&
                 2d0*sqrt(rad(2))*pim1(1))
         endif

         if(invm(i).eq.0d0) then
            smin(1)=max(smin(1),2d0*this%ETmin(i,1)*(pim1(0)-pim1(1)))
            smin(2)=max(smin(2),2d0*this%ETmin(i,2)*(pim1(0)-pim1(1)*cos(this%drcut(i+im1))))
         endif
         if (rad(2).lt.0d0) smax(2)=smin(2)

      endif
      if (smin(1).ge.smax(1)) then
         ps%jac=-4d0
         if (debug) write (*,*) 'smin.ge.smax',smin,smax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip,smin,smax,invm(i+im1),ps%jac)
      if (bad_forward_jac()) return
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
           &,sqrt0(invm(i)),sqrt0(invm(ir)),pi1,ppibir1,gent_status)
      if (gent_status.ne.0) then
         ps%jac=-9d0
         return
      endif
      pr1(0:3)=pp(0:3,ir+i)-pi1(0:3)
      phi2=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
           &,invm(ir+i+ib),V,sqrtGG,0d0)
      call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi2 &
           &,sqrt0(invm(i)),sqrt0(invm(ir)),pi2,ppibir2,gent_status)
      if (gent_status.ne.0) then
         ps%jac=-9d0
         return
      endif
      pr2(0:3)=pp(0:3,ir+i)-pi2(0:3)
      if ( pi1(0)**2-pi1(3)**2.ge.this%ETmin(i,1)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir,1)**2 .and. &
           pi2(0)**2-pi2(3)**2.ge.this%ETmin(i,1)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir,1)**2 ) then
         ix_e=ix_e+1
         if (ix_e.gt.size(ps%x)) then
            ps%jac=-47d0
            return
         endif
         if(ps%x(ix_e).gt.0.5d0) then
            pp(0:3,i)=pi1(0:3)
            pp(0:3,ir)=pr1(0:3)
            pp(0:3,ib+ir)=ppibir1(0:3)
         else
            pp(0:3,i)=pi2(0:3)
            pp(0:3,ir)=pr2(0:3)
            pp(0:3,ib+ir)=ppibir2(0:3)
         endif
      elseif (pi1(0)**2-pi1(3)**2.ge.this%ETmin(i,1)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir,1)**2) then
         pp(0:3,i)=pi1(0:3)
         pp(0:3,ir)=pr1(0:3)
         pp(0:3,ib+ir)=ppibir1(0:3)
         if (.not.update_forward_jac(1d0,2d0)) return
      elseif (pi2(0)**2-pi2(3)**2.ge.this%ETmin(i,1)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir,1)**2) then
         pp(0:3,i)=pi2(0:3)
         pp(0:3,ir)=pr2(0:3)
         pp(0:3,ib+ir)=ppibir2(0:3)
         if (.not.update_forward_jac(1d0,2d0)) return
      else
         ps%jac=-19d0
         if (debug) then
            write (*,*) 'piir',pp(0:3,i+ir)
            write (*,*) 'pim1',pp(0:3,im1)
            write (*,*) '1:',phi1,(phi1+phi2)/(2d0*pi)
            write (*,*) 'i',i,this%ETmin(i,1:2),':',pi1(0:3)
            write (*,*) 'ir',ir,this%ETmin(ir,1:2),':',pr1(0:3)
            write (*,*) '2:',phi2
            write (*,*) 'i',i,this%ETmin(i,1:2),':',pi2(0:3)
            write (*,*) 'ir',ir,this%ETmin(ir,1:2),':',pr2(0:3)
            write (*,*) ''
         endif
         return
      endif
      ! Compute the Jacobian
      gram4=gram_determinant4(invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib)&
           &,invm(ir+i),invm(i+im1),invm(ir+ib+i+im1),invm(ir),invm(i)&
           &,invm(im1))
      if (gram4.ge.0d0) then 
         write (99,*) 'error, gram4 greater than or equal to zero',gram4,i,ir
         ps%jac=-5d0
         return
      endif
      if (.not.update_forward_jac(1d0,8d0*sqrt(-gram4))) return
    end subroutine gen23_one_step





    subroutine gent_one_step(i,ir,ib)
      ! One step in the usual MadGraph t-channel phase-space generation.
      implicit none
      integer(kind=4),intent(in) :: i,ir,ib
      real(kind=8) :: tmin_S,tmax_S,phi,y,jac_den,jac_scale
      real(kind=8),dimension(1:2) :: shatmin,shatmax,Eimax,tmin,tmax,etminir,base,root,etmini
      real(kind=8),dimension(0:3) :: piir,pib
      logical :: rapidity_ok
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         if (popcnt(i+ir).eq.this%next-2) then
            ! The energy of i will be
            ! Ei=(sqrtshat+(invm(i)-invm(ir))/sqrtshat)/2d0. This gives a
            ! constraint on the allowed value of invm(i), since Ei>ETmin(i)
            shatmin(1:2)=max(shatmin(1:2),max(invm(ir),this%invm_min(ir,1:2))+sqrtshat*(2d0*this%ETmin(i,1:2)-sqrtshat))
            Eimax(1:2)=sqrtshat-this%ETmin(ir,1:2) ! maximum energy for i
            if (popcnt(ir).eq.1) then
               shatmax(1:2)=min(shatmax(1:2),invm(ir)+sqrtshat*(2d0*Eimax(1:2)-sqrtshat))
            else
               shatmax(1:2)=min(shatmax,Eimax(1:2)**2)
            endif
         endif
         call generate_mass(i,shatmin,shatmax)
         if (debug) then
            write (*,*) 't - i    ',i,invm(i),shatmin,shatmax
         endif
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         if (popcnt(i+ir).eq.this%next-2) then
            ! The energy of ir will be
            ! Eir=(sqrtshat+(invm(ir)-invm(i))/sqrtshat)/2d0. This gives a
            ! constraint on the allowed value of invm(ir), since Eir>ETmin(ir)
            shatmin(1:2)=max(shatmin(1:2),invm(i)+sqrtshat*(2d0*this%ETmin(ir,1:2)-sqrtshat))
            shatmax(1:2)=min(shatmax(1:2),invm(i)+sqrtshat*(sqrtshat-2d0*max(sqrt0(invm(i)),this%ETmin(i,1:2))))
         endif
         call generate_mass(ir,shatmin,shatmax)
         if (debug) then
            write (*,*) 't - ir   ',ir,invm(ir),shatmin,shatmax
         endif
      endif
      if (bad_forward_jac()) return
      call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin_S,tmax_S)
      tmin(1:2)=tmin_S
      tmax(1:2)=tmax_S
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2)=min(tmax_S,this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin_S,this%invm_min(ir+ib,1:2))
      end where
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
      call longitudinal_rapidity(pp(0:3,i+ir),y,rapidity_ok)
      if (.not.rapidity_ok) then
         ps%jac=-35d0
         return
      endif
      call boostz(pp(0,i+ir),y,piir)
      call boostz(pp(0,ib),y,pib)
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(i,1:2)**2-invm(i) .and. popcnt(i).eq.1 )
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(i,1:2)**2-invm(i)) )**2))
      elsewhere
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)))
      end where
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir,1:2)**2-invm(ir) .and. popcnt(ir).eq.1 )
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(ir,1:2)**2-invm(ir)) )**2))
      elsewhere
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)))
      end where
      base(1:2)=piir(0)**2-ETmini(1:2)**2+ETminir(1:2)**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root(1:2)=(piir(0)-ETmini(1:2)-ETminir(1:2))*(piir(0)+ETmini(1:2)-ETminir(1:2))*&
           (piir(0)-ETmini(1:2)+ETminir(1:2))*(piir(0)+ETmini(1:2)+ETminir(1:2))
      if (root(1).lt.0d0) then
         ps%jac=-33d0
         if (debug) write (*,*) 'root.lt.0d0',root
         return
      endif
      if (abs(piir(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
         ps%jac=-38d0
         if (debug) write (*,*) 'zero-energy system in gent_one_step',piir(0)
         return
      endif
      tmin(1)=max(tmin(1),invm(ir)-pib(0)/piir(0)*(base(1)+sqrt(root(1))))
      tmax(1)=min(tmax(1),invm(ir)-pib(0)/piir(0)*(base(1)-sqrt(root(1))))
      if (root(2).lt.0d0) then
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      else
         tmin(2)=max(tmin(2),invm(ir)-pib(0)/piir(0)*(base(2)+sqrt(root(2))))
         tmax(2)=min(tmax(2),invm(ir)-pib(0)/piir(0)*(base(2)-sqrt(root(2))))
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-3d0
         if (debug) write (*,*) 'tmin.ge.tmax',tmin,tmax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip,tmin,tmax,invm(ir+ib),ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) 't- ir+ib',ir+ib,invm(ir+ib),tmin,tmax
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],phi,ps%jac)
      if (bad_forward_jac()) return
      if (debug) then
         write (*,*) 't - phi  ',i,phi,0d0,2d0*pi
      endif
      call gentcms(pp(0,ib+ir+i),pp(0,ib),invm(ib+ir),phi,sqrt0(invm(i)) &
           &,sqrt0(invm(ir)),pp(0,i),pp(0,ib+ir))
      if (bad_forward_jac()) return
      pp(0:3,ir)=pp(0:3,ib+ir+i)+pp(0:3,ib)-pp(0:3,i)
      jac_den=sqrt0(lambda(invm(ir+i),0d0,invm(ir+i+ib)))
      jac_scale=max(spacing(0d0),abs(invm(ir+i)),abs(invm(ir+i+ib)))
      if (jac_den.le.tiny_kin*jac_scale) then
         ps%jac=-8d0
         return
      endif
      if (.not.update_forward_jac(1d0,4d0*jac_den)) return
    end subroutine gent_one_step

    subroutine gens_one_step(i,ir)
      implicit none
      integer(kind=4),intent(in) :: i,ir
      real(kind=8) :: esum,costh,phi
      real(kind=8),dimension(1:2) :: shatmin,shatmax
      real(kind=8),dimension(0:3) :: p_i,p_ir
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         call generate_mass(i,shatmin,shatmax)
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         call generate_mass(ir,shatmin,shatmax)
      endif
      if (bad_forward_jac()) return
      esum=sqrt0(invm(i+ir))
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,[-1d0,-1d0],[1d0,1d0],costh,ps%jac)
      if (bad_forward_jac()) return
      ix=ix+1
      call random_to_var(ps%x(ix),ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],phi,ps%jac)
      if (bad_forward_jac()) return
      call mom2cx(esum,sqrt0(invm(i)),sqrt0(invm(ir)),costh,phi,p_i,p_ir)
      if (bad_forward_jac()) return
      call boostm(p_i,pp(0:3,i+ir),esum,pp(0:3,i))
      call boostm(p_ir,pp(0:3,i+ir),esum,pp(0:3,ir))
      if (.not.update_forward_jac(sqrt0(lambda(invm(i+ir),invm(i),invm(ir))),&
           8d0*invm(i+ir))) return
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
      real(kind=8),dimension(1:2) :: shatmin,shatmax
      where (this%invm_min(i,1:2).ne.0d0)
         shatmin(1:2)=max(shatmin(1:2),this%invm_min(i,1:2))
      end where
      where (this%invm_max(i,1:2).ne.0d0)
         shatmax(1:2)=min(shatmax(1:2),this%invm_max(i,1:2))
      end where
      if (shatmin(1).ge.shatmax(1)) then
         ps%jac=-7d0
         if (debug) write (*,*) 'shatmin.ge.shatmax',i,shatmin,shatmax
         return
      endif
      ix=ix+1
      call random_to_var(ps%x(ix),this%ip_mass,shatmin,shatmax,invm(i),ps%jac)
      if (bad_forward_jac()) return
    end subroutine generate_mass

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
      real(kind=8) :: esum,mass1,mass2,costh1,phi1,md2,ed,pp,sinth1,pp2,scale
      real(kind=8),parameter :: rZero = 0d0,rHalf = 0.5d0,rOne = 1d0,rTwo = 2d0
      p1=0d0
      p2=0d0
      scale=max(spacing(0d0),esum**2,mass1**2,mass2**2)
      if (esum.le.vtiny*sqrt(scale) .or. abs(costh1).gt.1d0+vtiny) then
         ps%jac=-8d0
         return
      endif
      md2 = (mass1-mass2)*(mass1+mass2)
      ed = md2/esum
      if ( mass1*mass2.eq.rZero ) then
         pp = (esum-abs(ed))*rHalf
         if (pp.lt.-vtiny*sqrt(scale)) then
            ps%jac=-8d0
            return
         endif
         pp=max(pp,0d0)
      else
         pp2=(md2/esum)**2-rTwo*(mass1**2+mass2**2)+esum**2
         if (pp2.lt.-vtiny*scale) then
            ps%jac=-8d0
            return
         endif
         pp=sqrt(max(pp2,0d0))*rHalf
      endif
      sinth1 = sqrt(max((rOne-costh1)*(rOne+costh1),0d0))
      p1(0) = max((esum+ed)*rHalf,rZero)
      if (esum+ed.lt.-vtiny*sqrt(scale) .or. esum-ed.lt.-vtiny*sqrt(scale)) then
         ps%jac=-8d0
         p1=0d0
         p2=0d0
         return
      endif
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
      real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,ma2,pt,pt2,pt2_tol
      real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,p1_rot
      real(kind=8),parameter :: tiny=1d-5
      p1=0d0
      pr=0d0
      ptot(0:3)=pa(0:3)+pb(0:3)
      ptotm(0)=ptot(0)
      ptotm(1:3)=-ptot(1:3)
      ma2=dot(pa,pa)
      ! determine magnitude of p1 in cms frame (from dhelas routine mom2cx)
      ESUM2 = dot(ptot,ptot)
      if (esum2 .le. 0d0) then
         ps%jac=-8d0
         return
      endif
      esum=sqrt(esum2)
      MD2=(M1-M2)*(M1+M2)
      ED=MD2/ESUM
      IF (M1*M2.EQ.0.d0) THEN
         PP2=0.25d0*(ESUM-ABS(ED))**2
      ELSE
         PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
         if(pp2.lt.0d0) then
            ps%jac=-8d0
            return
         endif
      ENDIF
      call boostm(pa,ptotm,esum,pa_cms)
      E_acms = pa_cms(0)
      p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)
      if (p_acms.le.vtiny*max(spacing(0d0),abs(E_acms),esum)) then
         ps%jac=-8d0
         return
      endif
      ! define p1 in the frame where pa_cms is aligned with the positive z axis.
      p1(0) = MAX((ESUM+ED)*0.5d0,0.d0)
      if (esum+ed.le.0d0) then
         ps%jac=-8d0
         return
      endif
      p1(3) = -(m1**2+ma2-t-2d0*p1(0)*E_acms)/(2d0*p_acms)
      pt2=pp2-p1(3)**2
      pt2_tol=tiny*max(abs(esum2),abs(pp2),abs(p1(3)**2))
      if (pt2.lt.-pt2_tol) then
         ps%jac=-8d0
         return
      elseif (pt2.lt.0d0) then
         pt2=0d0
      endif
      pt = sqrt(pt2)
      p1(1) = pt*cos(phi)
      p1(2) = pt*sin(phi)
      call rotxxx(p1,pa_cms,p1_rot)       !Rotate p1 to the pa_cms frame
      call boostm(p1_rot,ptot,esum,p1)    !boost back to lab frame
      pr(0:3)=pa(0:3)-p1(0:3)         !Return remainder of momentum
    end subroutine gentcms

    subroutine random_to_var(x,power_in,var_min,var_max,var,jac)
      ! Given a random number x, it generates var in the range var_min
      ! <= var <= var_max according to var^(power)
      implicit none
      real(kind=8),intent(in) :: x
      real(kind=8),dimension(-1:1),intent(in) :: power_in
      real(kind=8),dimension(1:2),intent(in) :: var_min,var_max
      real(kind=8),intent(out) :: var
      real(kind=8),intent(inout) :: jac
      integer(kind=4) :: k
      real(kind=8),dimension(3) :: vmin,vmax,power,q
      real(kind=8) :: xloc,updated_jac
      logical :: valid_inputs
      var=0d0
      if (.not.ieee_is_finite(jac)) then
         jac=-1d0
         return
      endif
      if (jac.le.0d0) return
      if (.not.ieee_is_finite(x)) then
         jac=-1d0
         return
      endif
      if (x.lt.-vtiny .or. x.gt.1d0+vtiny) then
         jac=-1d0
         return
      endif
      call random_to_var_inputs(power_in,var_min,var_max,power,vmin,vmax, &
           this%use_soft_bounds_as_actual_limits,valid_inputs)
      if (.not.valid_inputs) then
         jac=-1d0
         return
      endif
      call random_to_var_weights(power,vmin,vmax,q)
      if (any(.not.ieee_is_finite(q))) then
         var=vmin(2)
         jac=-1d0
         return
      endif
      if (any(q.lt.0d0) .or. sum(q(1:3)).le.0d0) then
         var=vmin(2)
         jac=-1d0
         return
      endif
      xloc=max(0d0,min(1d0,x))
      if (q(1).gt.0d0 .and. xloc.le.q(1)) then
         k = 1
         xloc = xloc / q(1)
      elseif (q(2).gt.0d0 .and. xloc.le.q(1) + q(2)) then
         k = 2
         xloc = (xloc - q(1)) / q(2)
      else
         k = 3
         if (q(3).le.0d0) then
            var=vmax(3)
            jac=-1d0
            return
         endif
         xloc = (xloc - q(1) - q(2)) / q(3)
      endif
      call random_to_var_map(xloc,power(k),vmin(k),vmax(k),var,jac)
      if (var_min(1).lt.0d0 .and. var_max(1).le.0d0) then
         var=-var
      endif
      if (.not.safe_phase_space_ratio(jac,q(k),updated_jac)) then
         var=0d0
         jac=-1d0
         return
      endif
      jac=updated_jac
      if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac)) then
         var=0d0
         jac=-1d0
         return
      endif
      if (jac.le.0d0) then
         var=0d0
         jac=-1d0
      endif
    end subroutine random_to_var

    subroutine random_to_var_map(x,power,varmin,varmax,var,jac)
      implicit none
      real(kind=8),intent(in) :: x,power,varmin,varmax
      real(kind=8),intent(out) :: var
      real(kind=8),intent(inout) :: jac
      real(kind=8) :: jac_factor,updated_jac
      real(kind=8) :: interval_scale,interval_tolerance
      logical :: valid_power_map
      var=0d0
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(power) .or. &
           .not.ieee_is_finite(varmin) .or. .not.ieee_is_finite(varmax) .or. &
           .not.ieee_is_finite(jac)) then
         jac=-1d0
         return
      endif
      interval_scale=max(abs(varmin),abs(varmax))
      interval_tolerance=128d0*epsilon(1d0)*max(spacing(0d0),interval_scale)
      if (varmax.le.varmin .or. varmax-varmin.le.interval_tolerance) then
         var=varmin
         jac=-1d0
         return
      endif
      if (x.lt.0d0 .or. x.gt.1d0) then
         jac=-1d0
         return
      endif
      if (power.eq.0d0) then
         var=varmin+x*(varmax-varmin)
         jac_factor=varmax-varmin
      elseif (power.eq.-1d0) then
         if (varmin.le.0d0) then
            jac=-1d0
            return
         endif
         jac_factor=log(varmax)-log(varmin)
         if (.not.ieee_is_finite(jac_factor)) then
            jac=-1d0
            return
         endif
         if (jac_factor.le.0d0) then
            jac=-1d0
            return
         endif
         var=exp((1d0-x)*log(varmin)+x*log(varmax))
         jac_factor=var*jac_factor
      elseif (power.eq.101d0) then
         var=varmin+(varmax-varmin)*(1d0-cos(pi*x))/2d0
         jac_factor=(varmax-varmin)*pi*sin(pi*x)/2d0
      else
         call stable_phase_space_power_map(x,power,varmin,varmax,var,&
              jac_factor,valid_power_map)
         if (.not.valid_power_map) then
            jac=-1d0
            return
         endif
      endif
      if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(jac_factor)) then
         var=0d0
         jac=-1d0
         return
      endif
      if (jac_factor.le.0d0) then
         var=0d0
         jac=-1d0
         return
      endif
      if (.not.safe_phase_space_product(jac,jac_factor,updated_jac)) then
         var=0d0
         jac=-1d0
         return
      endif
      jac=updated_jac
      if (jac.le.0d0) then
         var=0d0
         jac=-1d0
      endif
    end subroutine random_to_var_map
  end subroutine gen23_generate_momenta

  subroutine gen23_compute_x_from_momenta(this,ps)
    implicit none
    class(phase_space_gen23),intent(inout) :: this
    type(psv),intent(inout) :: ps
    real(kind=8) :: ycm,sqrtshat
    integer :: ix
    real(kind=8),dimension(0:3,0:maskr(this%next)) :: pp
    real(kind=8),dimension(maskr(this%next)) :: invm
    real(kind=8),dimension(0:3,2) :: p_tmp
    if (debug) write (*,*) 'computing x from momenta'
    ps%jac=-48d0
    if (.not.allocated(ps%p) .or. .not.allocated(ps%x)) return
    if (size(ps%p,1).ne.4 .or. size(ps%p,2).ne.this%next .or. &
         lbound(ps%p,1).ne.0 .or. size(ps%x).lt.this%ndim) return
    if (.not.generated_momenta_are_valid(ps%p,this%masses,ps%xbjrk,this%include_pdf)) return
    ps%x=0d0
    ps%jac=1d0
    ix=0
    
!!$       p_tmp(0:3,1)=ps%p(0:3,this%next-1)
!!$       p_tmp(0:3,2)=ps%p(0:3,this%next)
!!$       ps%p(0:3,this%next-1)=p_tmp(0:3,1)+p_tmp(0:3,2)
!!$       call generate_bw_mass_inverse(this%next-1)
!!$       this%next=this%next-1
       call setup_PS_cuts(this)
       
       invm=this%invm

    ! Fill the full momentum array, including all possible
    ! intermediate states:
    call fill_momentum_array
    if (bad_inverse_jac()) return
    ! get the two random number corresponding to the initial state
    if (this%include_pdf) then
       call compute_x_initial_state
       if (bad_inverse_jac()) return
    else
      sqrtshat=this%sqrts
    endif
    ! The final-state momenta configuration gives all the other random numbers
    call compute_x_final_state
    if (bad_inverse_jac()) return
    if (ix.ne.this%ndim .or. .not.ieee_is_finite(ps%jac) .or. &
         any(.not.ieee_is_finite(ps%x(1:this%ndim)))) then
       ps%jac=-47d0
       return
    endif

!!$       this%next=this%next+1
!!$       pp(0:3,ibset(0,this%next-2))=p_tmp(0:3,1)
!!$       pp(0:3,ibset(0,this%next-1))=p_tmp(0:3,2)
!!$       pp(0:3,ibset(0,this%next-1)+ibset(0,this%next-2))=p_tmp(0:3,1)+p_tmp(0:3,2)
!!$       invm(ibset(0,this%next-1)+ibset(0,this%next-2))=invm(ibset(0,this%next-2))
!!$       invm(ibset(0,this%next-2))=0d0
!!$       invm(ibset(0,this%next-1))=0d0
!!$       call decay_bw_inverse(this%next-1,this%next-1,this%next)
!!$       ! Add factors of 2*pi (since this%next was reduced by 1)
!!$       ps%jac=ps%jac/((2d0*pi)**3)
!!$       ps%p(0:3,this%next-1)=p_tmp(0:3,1)
!!$       ps%p(0:3,this%next)=p_tmp(0:3,2)
    
  contains
    logical function bad_inverse_jac()
      implicit none
      ! Inverse-status codes used below: -2/-3/-4 are singular or empty
      ! kinematic intervals, -5 is a non-negative Gram determinant,
      ! -33 is a negative reconstruction root, -35 is a light-like boost
      ! reference, -36/-37 are double-t bound failures, -38 is a zero
      ! reconstructed energy, -39 is a zero Eir denominator, -40 is an
      ! otherwise unclassified zero Jacobian, -41:-43 are inverse-variable
      ! bounds/range failures, -44 is invalid initial-state kinematics,
      ! -45/-46 are singular Eir/radial bounds, -47 is inconsistent random
      ! bookkeeping, and -48 is a non-finite generated momentum.
      if (.not.ieee_is_finite(ps%jac)) then
         ps%jac=-47d0
         bad_inverse_jac=.true.
         return
      endif
      bad_inverse_jac=ps%jac.le.0d0
      ! Preserve the diagnostic status set by the inverse step.  In
      ! particular, the multichannel caller needs to distinguish a point
      ! outside an alternative channel's bounds from a Gram-determinant or
      ! other kinematic reconstruction failure.
      if (bad_inverse_jac .and. ps%jac.eq.0d0) ps%jac=-40d0
    end function bad_inverse_jac

    logical function update_inverse_jac(numerator,denominator)
      implicit none
      real(kind=8),intent(in) :: numerator,denominator
      real(kind=8) :: updated_jac
      update_inverse_jac=.false.
      if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator) .or. &
           numerator.le.0d0 .or. denominator.le.0d0) then
         ps%jac=-47d0
         return
      endif
      if (.not.safe_phase_space_scaled_ratio(ps%jac,numerator,denominator,updated_jac)) then
         ps%jac=-47d0
         return
      endif
      ps%jac=updated_jac
      update_inverse_jac=.true.
    end function update_inverse_jac

    subroutine generate_bw_mass_inverse(ires)
      implicit none
      integer,intent(in) :: ires
      real(kind=8) :: smin,smax,qmass,qwidth,y
      real(kind=8),dimension(1:2) :: A,B
      smin=50d0**2
      smax=this%sqrts**2
      qmass=z_mass
      qwidth=z_width
      A=atan((qmass-smin/qmass)/qwidth)
      B=atan((qmass-smax/qmass)/qwidth)
      ix=ix+1
      this%invm(ibset(0,ires-1))=dot(ps%p(0:3,ires),ps%p(0:3,ires))
      y=atan((qmass-this%invm(ibset(0,ires-1))/qmass)/qwidth)
      call var_to_random(y,ip_flat,B,A,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      if (.not.update_inverse_jac(qmass*qwidth,cos(y)**2)) return
    end subroutine generate_bw_mass_inverse
    subroutine decay_bw_inverse(ires,id1,id2)
      implicit none
      integer,intent(in) :: ires,id1,id2
      integer :: i,j,k
      real(kind=8),dimension(0:3,this%next) :: pp_tmp
      if (ires.gt.id1 .or. ires.gt.id2) then
         write (*,*) 'ERROR in decay_BW',ires,id1,id2
      endif
      call gens_one_step_inverse(ibset(0,id1-1),ibset(0,id2-1))
    end subroutine decay_bw_inverse
       
    subroutine compute_x_initial_state
      implicit none
      real(kind=8) :: tau
      call compute_x_from_tau(tau)
      if (bad_inverse_jac()) return
      call compute_x_from_y(tau)
      if (bad_inverse_jac()) return
    end subroutine compute_x_initial_state
    subroutine compute_x_from_tau(tau)
      implicit none
      real(kind=8),intent(out) :: tau
      real(kind=8) :: shat
      real(kind=8),dimension(1:2) :: smin,smax
      tau=0d0
      smin(1:2)=max(this%invm_min(maskr(this%next)-3,1:2),this%ETmin(maskr(this%next)-3,1:2)**2)
      smax(1:2)=this%sqrts**2
      shat=dot(pp(0:3,3),pp(0:3,3))
      if (shat.le.0d0) then
         ps%jac=-44d0
         return
      endif
      sqrtshat=sqrt(shat)
      ix=ix+1
      call var_to_random(shat,this%ip_shat,smin,smax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      tau=shat/smax(1)
      if (.not.update_inverse_jac(1d0,smax(1))) return
    end subroutine compute_x_from_tau
    subroutine compute_x_from_y(tau)
      implicit none
      real(kind=8),intent(in) :: tau
      real(kind=8),dimension(1:2) ::  ymin,ymax
      ymin(1:2)= log(tau)/2d0
      ymax(1:2)=-log(tau)/2d0
      ix=ix+1
      call var_to_random(ycm,ip_flat,ymin,ymax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
    end subroutine compute_x_from_y
    subroutine compute_x_final_state
      implicit none
      integer(kind=4) :: i,j,inext,im1
      integer(kind=4),dimension(2) :: set
      ! invariant mass of all final state particles combined
      invm(ibset(0,0)+ibset(0,1))=sqrtshat**2
      invm(maskr(this%next)-ibset(0,0)-ibset(0,1))=sqrtshat**2
      set(1)=this%sets(0,1)
      set(2)=this%sets(0,2)
      ! Generate the central 2->2 process in case both set(1) and set(2) are not empty
      if (popcnt(set(1)).gt.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'two sets with at least two ',&
              & 'particles',popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         if (use_t_channel_at_start) then
            call gent_one_step_inverse(set(2),set(1),1)
         else
            call gens_one_step_inverse(set(2),set(1))
         endif
         if (bad_inverse_jac()) return
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).gt.1) then
         if (debug) write (*,*) 'special double t-channel (1)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t_inverse(set(1),set(2),1,2)
         if (bad_inverse_jac()) return
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      elseif (popcnt(set(1)).gt.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) 'special double t-channel (2)'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call double_t_inverse(set(2),set(1),1,2)
         if (bad_inverse_jac()) return
         invm(set(1)+1)=dot(pp(0:3,set(1)+1),pp(0:3,set(1)+1))
      elseif (popcnt(set(1)).eq.1 .and. popcnt(set(2)).eq.1) then
         if (debug) write (*,*) '2->2 scattering with one particle in each set'&
              &,popcnt(this%sets(0,1)),popcnt(this%sets(0,2))
         call gent_one_step_inverse(set(2),set(1),1)
         if (bad_inverse_jac()) return
         invm(set(2)+2)=dot(pp(0:3,set(2)+2),pp(0:3,set(2)+2))
      endif

      do i=1,2
         if (popcnt(set(i)).le.1) cycle ! at least 2 particles in a set
         inext=ibset(0,this%sets(1,i)-1)
         set(i)=set(i)-inext
         if (popcnt(set(i)).ge.2) then
            ! at least 3 particles in a set
            if (debug) write (*,*) 'At least 3 particles in a set',&
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step_inverse(set(i),inext,i)
            if (bad_inverse_jac()) return
            invm((3-i)+set(i)+inext)=dot(pp(0:3,(3-i)+set(i)+inext),pp(0:3,(3-i)+set(i)+inext))
            invm((3-i)+set(i))=dot(pp(0:3,(3-i)+set(i)),pp(0:3,(3-i)+set(i)))
            do j=2,popcnt(set(i))-1
               ! loop over the remaining particles in the set
               inext=ibset(0,this%sets(j,i)-1)
               im1=ibset(0,this%sets(j-1,i)-1)
               set(i)=set(i)-inext
               if (this%t_channel) then
                  call gent_one_step_inverse(inext,set(i),3-i)
               else
                  call gen23_one_step_inverse(inext,set(i),3-i,im1)
               endif
               if (bad_inverse_jac()) return
            enddo
            inext=ibset(0,this%sets(j,i)-1)
            im1=ibset(0,this%sets(j-1,i)-1)
            set(i)=set(i)-inext
            if (this%t_channel) then
               call gent_one_step_inverse(inext,set(i),3-i)
            else
               call gen23_one_step_inverse(inext,set(i),3-i,im1)
            endif
            if (bad_inverse_jac()) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).ne.0) then
            ! Exactly 2 particles in a set (and the other set contains at least one)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains at least one)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            im1=3-i
            invm(set(i)+inext+im1)=dot(pp(0:3,set(i)+inext+im1),pp(0:3,set(i)+inext+im1))
            invm(set(i)+inext+i+im1)=dot(pp(0:3,set(i)+inext+im1+i),pp(0:3,set(i)+inext+im1+i))
            if (this%t_channel) then
               call gent_one_step_inverse(set(i),inext,i)
            else
               call gen23_one_step_inverse(set(i),inext,i,im1)
            endif
            if (bad_inverse_jac()) return
         elseif (popcnt(set(i)).eq.1 .and. popcnt(this%sets(0,3-i)).eq.0) then
            ! Exactly 2 particles in a set (and the other set contains none)
            if (debug) write (*,*) 'Exactly 2 particles in a set (and ', &
                 & 'the other set contains none)', &
                 & popcnt(this%sets(0,i)),popcnt(this%sets(0,3-i))
            call gent_one_step_inverse(set(i),inext,i)
            if (bad_inverse_jac()) return
         else
            write (*,*) 'Inconsistent sets'
            write (*,*) i,':',this%sets(:,i)
            write (*,*) 3-i,':',this%sets(:,3-i)
            ps%jac=-48d0
            return
         endif
      enddo

      ! Add factors of 2*pi
      if (.not.update_inverse_jac(1d0,(2d0*pi)**(3*(this%next-2)-4))) return
      ! Add flux factor
      if (.not.update_inverse_jac(1d0,2d0*sqrtshat**2)) return

    end subroutine compute_x_final_state
    subroutine gent_one_step_inverse(i,ir,ib)
      implicit none
      integer(kind=4),intent(in) :: i,ir,ib
      real(kind=8) :: tmin_S,tmax_S,phi,y,esum,jac_den,jac_scale
      real(kind=8),dimension(1:2) :: shatmin,shatmax,Eimax,tmin,tmax,etminir,base,root,etmini
      real(kind=8),dimension(0:3) :: piir,pib,p_boost,pi_rot,pa_cms
      logical :: rapidity_ok
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         if (popcnt(i+ir).eq.this%next-2) then
            shatmin(1:2)=max(shatmin(1:2),max(invm(ir),this%invm_min(ir,1:2))+sqrtshat*(2d0*this%ETmin(i,1:2)-sqrtshat))
            Eimax(1:2)=sqrtshat-this%ETmin(ir,1:2) ! maximum energy for i
            if (popcnt(ir).eq.1) then
               shatmax(1:2)=min(shatmax(1:2),invm(ir)+sqrtshat*(2d0*Eimax(1:2)-sqrtshat))
            else
               shatmax(1:2)=min(shatmax(1:2),Eimax(1:2)**2)
            endif
         endif
         if (debug) write (*,*) 'generate_mass_inverse gent 1',i,ir
         call generate_mass_inverse(i,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         if (popcnt(i+ir).eq.this%next-2) then
            shatmin(1:2)=max(shatmin(1:2),invm(i)+sqrtshat*(2d0*this%ETmin(ir,1:2)-sqrtshat))
            shatmax(1:2)=min(shatmax(1:2),invm(i)+sqrtshat*(sqrtshat-2d0*max(sqrt0(invm(i)),this%ETmin(i,1:2))))
         endif
         if (debug) write (*,*) 'generate_mass_inverse gent 2',ir
         call generate_mass_inverse(ir,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin_S,tmax_S)
      tmin(1:2)=tmin_S
      tmax(1:2)=tmax_S
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2)=min(tmax_S,this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin_S,this%invm_min(ir+ib,1:2))
      end where
      ! Make sure that the t-range is compatible with the pT cut. Since t is an
      ! invariant we can compute it in any frame. Let's use the frame in which
      ! p(:,i+ir) has p_z=0, since in this frame p_z(i)=-p_z(ir). (Note that
      ! ETmin() is boost invariant in the z-direction)
      call longitudinal_rapidity(pp(0:3,i+ir),y,rapidity_ok)
      if (.not.rapidity_ok) then
         ps%jac=-35d0
         return
      endif
      call boostz(pp(0,i+ir),y,piir)
      call boostz(pp(0,ib),y,pib)
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(i,1:2)**2-invm(i) .and. popcnt(i).eq.1 )
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(i,1:2)**2-invm(i)) )**2))
      elsewhere
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)))
      end where
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir,1:2)**2-invm(ir) .and. popcnt(ir).eq.1 )
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(ir,1:2)**2-invm(ir)) )**2))
      elsewhere
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)))
      end where
      base(1:2)=piir(0)**2-ETmini(1:2)**2+ETminir(1:2)**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root(1:2)=(piir(0)-ETmini(1:2)-ETminir(1:2))*(piir(0)+ETmini(1:2)-ETminir(1:2))*&
           (piir(0)-ETmini(1:2)+ETminir(1:2))*(piir(0)+ETmini(1:2)+ETminir(1:2))
      if (root(1).lt.0d0) then
         ps%jac=-33d0
         if (debug) write (*,*) 'root.lt.0d0 in gent_one_step_inverse',root
         return
      endif
      if (abs(piir(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
         ps%jac=-38d0
         return
      endif
      tmin(1)=max(tmin(1),invm(ir)-pib(0)/piir(0)*(base(1)+sqrt(root(1))))
      tmax(1)=min(tmax(1),invm(ir)-pib(0)/piir(0)*(base(1)-sqrt(root(1))))
      if (root(2).lt.0d0) then
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      else
         tmin(2)=max(tmin(2),invm(ir)-pib(0)/piir(0)*(base(2)+sqrt(root(2))))
         tmax(2)=min(tmax(2),invm(ir)-pib(0)/piir(0)*(base(2)-sqrt(root(2))))
      endif
      invm(ir+ib)=dot(pp(0:3,ir+ib),pp(0:3,ir+ib))
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-3d0
         if (debug) write (*,*) 'tmin.ge.tmax in gent_one_step_inverse',tmin,tmax,invm(ir+ib)
         return
      endif
      if (debug) then
         write (*,*) 'ti - ir+ib',ir+ib,tmin,tmax,invm(ir+ib)
      endif
      ix=ix+1
      call var_to_random(invm(ir+ib),this%ip,tmin,tmax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      ! inverse of boosts and rotation from gentcms()
      jac_scale=max(spacing(0d0),abs(invm(i+ir)))
      if (invm(i+ir).le.tiny_kin*jac_scale) then
         ps%jac=-44d0
         return
      endif
      esum=sqrt(invm(i+ir))
      p_boost(0)=pp(0,i+ir)
      p_boost(1:3)=-pp(1:3,i+ir)
      call boostm(pp(0:3,i),p_boost,esum,pib)
      call boostm(pp(0:3,ib+ir+i),p_boost,esum,pa_cms)
      call rotxxx_inv(pib,pa_cms,pi_rot)
      phi=atan2(pi_rot(2),pi_rot(1))
      if(phi.lt.0d0) phi=phi+2d0*pi
      if (debug) then
         write (*,*) 'ti - phi',0d0,2d0*pi,phi
      endif
      ix=ix+1
      call var_to_random(phi,ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      jac_den=sqrt0(lambda(invm(ir+i),0d0,invm(ir+i+ib)))
      jac_scale=max(spacing(0d0),abs(invm(ir+i)),abs(invm(ir+i+ib)))
      if (jac_den.le.tiny_kin*jac_scale) then
         ps%jac=-44d0
         return
      endif
      if (.not.update_inverse_jac(1d0,4d0*jac_den)) return
    end subroutine gent_one_step_inverse
    subroutine gens_one_step_inverse(i,ir)
      implicit none
      integer(kind=4),intent(in) :: i,ir
      real(kind=8) :: esum,costh,phi,pnorm,jac_scale
      real(kind=8),dimension(1:2) :: shatmin,shatmax
      real(kind=8),dimension(0:3) :: p_i,p_boost
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         if (debug) write (*,*) 'generate_mass_inverse gens 1',i
         call generate_mass_inverse(i,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         if (debug) write (*,*) 'generate_mass_inverse gent 2',ir
         call generate_mass_inverse(ir,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      ! boost p(i) and p(ir) to the p(i+ir) rest frame
      jac_scale=max(spacing(0d0),abs(invm(i+ir)))
      if (invm(i+ir).le.tiny_kin*jac_scale) then
         ps%jac=-44d0
         return
      endif
      esum=sqrt(invm(i+ir))
      p_boost(0)=-pp(0,i+ir)
      p_boost(1:3)=-pp(1:3,i+ir)
      call boostm(pp(0:3,i),p_boost,esum,p_i)
      ! compute the angles from the momenta and use that to get the random numbers
      pnorm=sqrt(p_i(1)**2+p_i(2)**2+p_i(3)**2)
      if (pnorm.le.tiny_kin*max(spacing(0d0),abs(p_i(0)))) then
         ps%jac=-44d0
         return
      endif
      costh=max(-1d0,min(1d0,p_i(3)/pnorm))
      if (debug) then
         write (*,*) 'si - i',i,-1d0,1d0,costh
      endif
      ix=ix+1
      call var_to_random(costh,ip_flat,[-1d0,-1d0],[1d0,1d0],ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      phi=atan2(p_i(2),p_i(1))
      if(phi.lt.0d0) phi=phi+2d0*pi
      if (debug) then
         write (*,*) 'si - phi i',i,0d0,2d0*pi,phi
      endif
      ix=ix+1
      call var_to_random(phi,ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      ! update the Jacobian
      if (.not.update_inverse_jac(sqrt0(lambda(invm(i+ir),invm(i),invm(ir))),&
           8d0*invm(i+ir))) return
      ! compute some t-channel invariants just to make sure they are filled. 
      invm(i+1)=dot(pp(0:3,i+1),pp(0:3,i+1))
      invm(i+2)=dot(pp(0:3,i+2),pp(0:3,i+2))
      invm(ir+1)=dot(pp(0:3,ir+1),pp(0:3,ir+1))
      invm(ir+2)=dot(pp(0:3,ir+2),pp(0:3,ir+2))
    end subroutine gens_one_step_inverse
    subroutine double_t_inverse(i,ir,ia,ib)
      implicit none
      integer(kind=4),intent(in) :: i,ir,ia,ib
      real(kind=8) :: phi,mass_tol,den_t,den_ir,den_scale,jac_den,jac_scale
      real(kind=8),dimension(1:2) :: tmin,tmax,yr,Eimax,pzmax
      logical :: soft_ok,massless_collinear_limit
      if (popcnt(i).ne.1 .or. popcnt(ir).le.1) then
         write (*,*) 'Subroutine only for i is a single particle '&
              //'and ir is more than 1',i,ir,popcnt(i),popcnt(ir)
         ps%jac=-48d0
         return
      endif
      soft_ok=.true.
      yr(1:2)=lambda(invm(ia+ib),invm(i),this%invm_min(ir,1:2))
      if (yr(1).lt.0d0) then
         ps%jac=-36d0
         return
      endif
      if (yr(2).lt.0d0) soft_ok=.false.
      yr(1)=sqrt(yr(1))
      if (soft_ok) yr(2)=sqrt(yr(2))
      tmin(1)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,1)-yr(1))/2d0
      tmax(1)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,1)+yr(1))/2d0
      if (soft_ok) then
         tmin(2)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,2)-yr(2))/2d0
         tmax(2)=(-invm(ia+ib)+invm(i)+this%invm_min(ir,2)+yr(2))/2d0
      else
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      endif
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2) = min(tmax(1:2), this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin(1:2),this%invm_min(ir+ib,1:2))
      end where
      pzmax(1:2)=lambda(sqrtshat**2,this%ETmin(i,1:2)**2,this%ETmin(ir,1:2)**2)
      if (pzmax(1).lt.0d0) then
         ps%jac=-37d0
         return
      endif
      if (pzmax(2).lt.0d0) soft_ok=.false.
      pzmax(1)=sqrt(pzmax(1))/(2d0*sqrtshat)
      Eimax(1)=sqrtshat-sqrt(this%ETmin(ir,1)**2+pzmax(1)**2)
      tmin(1)=max(tmin(1),invm(i)-sqrtshat*(Eimax(1)+pzmax(1)))
      tmax(1)=min(tmax(1),invm(i)-sqrtshat*(Eimax(1)-pzmax(1)))
      if (soft_ok) then
         pzmax(2)=sqrt(pzmax(2))/(2d0*sqrtshat)
         Eimax(2)=sqrtshat-sqrt(this%ETmin(ir,2)**2+pzmax(2)**2)
         tmin(2)=max(tmin(2),invm(i)-sqrtshat*(Eimax(2)+pzmax(2)))
         tmax(2)=min(tmax(2),invm(i)-sqrtshat*(Eimax(2)-pzmax(2)))
      else
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-2d0
         if (debug) write (*,*) 'tmin.ge.tmax in double_t_inverse',tmin,tmax
         return
      endif
      invm(i+ia)=dot(pp(0:3,i+ia),pp(0:3,i+ia))
      if (debug) then
         write (*,*) 'dti- i+ia',i+ia,tmin,tmax,invm(i+ia)
      endif
      ix=ix+1
      call var_to_random(invm(i+ia),this%ip_dt,tmin,tmax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      den_t=invm(i)-invm(i+ia)
      den_ir=sqrtshat**2+invm(i+ia)-invm(i)
      den_scale=max(spacing(0d0),abs(invm(i)),abs(invm(i+ia)),sqrtshat**2)
      massless_collinear_limit=abs(den_t).le.vtiny*den_scale .and. &
           abs(invm(i)).le.vtiny*den_scale .and. &
           this%ETmin(i,1).le.vtiny*sqrt(den_scale)
      if (abs(den_ir).le.vtiny*den_scale .or. &
           (abs(den_t).le.vtiny*den_scale .and. .not.massless_collinear_limit)) then
         ps%jac=-2d0
         if (debug) write (*,*) 'singular inverse double_t constraint',den_t,den_ir
         return
      endif
      tmin(1:2)=-invm(ia+ib)-invm(i+ia)+invm(i)+this%invm_min(ir,1:2)
      if (massless_collinear_limit) then
         tmax(1:2)=0d0
      else
         tmax(1:2)=invm(i)*(invm(i)-invm(ia+ib)-invm(i+ia))/den_t
      endif
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2) = min(tmax(1:2), this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin(1:2),this%invm_min(ir+ib,1:2))
      end where
      ! Additional constraints on tmin and tmax due to pp(0,i) and pp(0,ir)
      ! being larger than ETmin(i) and ETmin(ir), respectively:
      tmin(1:2)=max(tmin(1:2),invm(i)-sqrtshat**2*(1-this%ETmin(ir,1:2)**2/den_ir))
      if (.not.massless_collinear_limit) then
         tmax(1:2)=min(tmax(1:2),invm(i)-sqrtshat**2*(this%ETmin(i,1:2)**2/den_t))
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-2d0
         if (debug) write (*,*) 'tmin.ge.tmax in double_t_inverse',tmin,tmax
         return
      endif
      invm(i+ib)=dot(pp(0:3,i+ib),pp(0:3,i+ib))
      if (debug) then
         write (*,*) 'dti- i+ib',i+ib,tmin,tmax,invm(i+ib)
      endif
      ix=ix+1
      call var_to_random(invm(i+ib),this%ip_dt,tmin,tmax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      phi=atan2(pp(2,i),pp(1,i))
      if(phi.lt.0d0) phi=phi+2d0*pi
      if (debug) then
         write (*,*) 'dti- phi',i,0d0,2d0*pi,phi
      endif
      ix=ix+1
      call var_to_random(phi,ip_flat,[0d0,0d0],[2d0*pi,2d0*pi],ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      invm(ir)=dot(pp(0:3,ir),pp(0:3,ir))
      mass_tol=vtiny*max(spacing(0d0),abs(invm(ia+ib)),abs(invm(i+ia)), &
           abs(invm(i+ib)),abs(invm(i)))
      if (invm(ir).lt.-mass_tol) then
         ps%jac=-4d0
         if (debug) write (*,*) "rejecting inverse double_t point with negative remainder mass",ir,invm(ir),i
         return
      elseif (invm(ir).lt.0d0) then
         invm(ir)=0d0
      endif
      jac_den=sqrt0(lambda(invm(ir+i),0d0,0d0))
      jac_scale=max(spacing(0d0),abs(invm(ir+i)))
      if (jac_den.le.tiny_kin*jac_scale) then
         ps%jac=-44d0
         return
      endif
      if (.not.update_inverse_jac(1d0,4d0*jac_den)) return
    end subroutine double_t_inverse
    subroutine gen23_one_step_inverse(i,ir,ib,im1)
      implicit none
      integer(kind=4),intent(in) :: im1,i,ir,ib
      integer :: gent_status
      real(kind=8) :: tmin_S,tmax_S,smin_S,smax_S,phi1,phi2,gram4,V,sqrtGG,y,phi_rot,delta_ir,delta_scale
      real(kind=8),dimension(1:2) :: shatmin,shatmax,tmin,tmax,etminir,etmini,base,root,smin,smax,rad
      real(kind=8),dimension(0:3) :: pi1,pr1,ppibir1,pi2,pr2,ppibir2,piir,pib,pim1,piirr,pim1r
      logical :: rapidity_ok
      if (popcnt(i).gt.1) then
         if (popcnt(ir).gt.1) invm(ir)=0d0 ! set this mass to zero to get the correct smax limit in shatminmax
         call shatminmax(this,i,ir,shatmin,shatmax,invm)
         if (debug) write (*,*) 'generate_mass_inverse gen23 1',i
         call generate_mass_inverse(i,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      if (popcnt(ir).gt.1) then
         call shatminmax(this,ir,i,shatmin,shatmax,invm)
         if (debug) write (*,*) 'generate_mass_inverse gen23 2',ir
         call generate_mass_inverse(ir,shatmin,shatmax)
         if (bad_inverse_jac()) return
      endif
      call tminmax(invm(ir+i),invm(ir+i+ib),invm(ir),invm(i),0d0,tmin_S,tmax_S)
      tmin(1:2)=tmin_S
      tmax(1:2)=tmax_S
      where (this%invm_max(ir+ib,1:2).ne.0d0)
         tmax(1:2)=min(tmax_S,this%invm_max(ir+ib,1:2))
      end where
      where (this%invm_min(ir+ib,1:2).ne.0d0)
         tmin(1:2)=max(tmin_S,this%invm_min(ir+ib,1:2))
      end where
      pp(0:3,i+ir)=pp(0:3,i+ir+ib)+pp(0:3,ib)
      call longitudinal_rapidity(pp(0:3,i+ir),y,rapidity_ok)
      if (.not.rapidity_ok) then
         ps%jac=-35d0
         return
      endif
      call boostz(pp(0,i+ir),y,piir)
      call boostz(pp(0,ib),y,pib)
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(i,1:2)**2-invm(i) .and. popcnt(i).eq.1 )
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(i,1:2)**2-invm(i)) )**2))
      elsewhere
         etminir(1:2)=max(this%ETmin(ir,1:2),sqrt0(invm(ir)))
      end where
      where ( piir(1)**2+piir(2)**2.lt.this%ETmin(ir,1:2)**2-invm(ir) .and. popcnt(ir).eq.1 )
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)+abs(sqrt(piir(1)**2+piir(2)**2)-&
              sqrt0(this%ETmin(ir,1:2)**2-invm(ir)) )**2))
      elsewhere
         etmini(1:2)=max(this%ETmin(i,1:2),sqrt0(invm(i)))
      end where
      base(1:2)=piir(0)**2-ETmini(1:2)**2+ETminir(1:2)**2
      ! Note, root=lambda(piir(0)**2,this%ETmin(i)**2,this%ETmin(ir)**2), but the
      ! following is more stable:
      root(1:2)=(piir(0)-ETmini(1:2)-ETminir(1:2))*(piir(0)+ETmini(1:2)-ETminir(1:2))*&
           (piir(0)-ETmini(1:2)+ETminir(1:2))*(piir(0)+ETmini(1:2)+ETminir(1:2))
      if (root(1).lt.0d0) then
         ps%jac=-33d0
         return
      endif
      if (abs(piir(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
         ps%jac=-38d0
         return
      endif
      tmin(1)=max(tmin(1),invm(ir)-pib(0)/piir(0)*(base(1)+sqrt0(root(1))))
      tmax(1)=min(tmax(1),invm(ir)-pib(0)/piir(0)*(base(1)-sqrt0(root(1))))
      if (root(2).lt.0d0) then
         tmin(2)=tmin(1)
         tmax(2)=tmin(1)
      else
         tmin(2)=max(tmin(2),invm(ir)-pib(0)/piir(0)*(base(2)+sqrt0(root(2))))
         tmax(2)=min(tmax(2),invm(ir)-pib(0)/piir(0)*(base(2)-sqrt0(root(2))))
      endif
      if (tmin(1).ge.tmax(1)) then
         ps%jac=-3d0
         return
      endif
      invm(ir+ib)=dot(pp(0:3,ir+ib),pp(0:3,ir+ib))
      if (debug) then
         write (*,*) '23i- ir+ib',ir+ib,tmin,tmax,invm(ir+ib)
      endif
      ix=ix+1
      call var_to_random(invm(ir+ib),this%ip,tmin,tmax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      call sminmax(invm(ir+i),invm(ir),invm(ir+i+im1),invm(ir+i+ib)&
           &,invm(ir+ib),invm(ir+ib+i+im1),invm(i),invm(im1),smin_S,smax_S,V,sqrtGG)
      smin(1:2)=smin_S
      smax(1:2)=smax_S
      where (this%invm_min(i+im1,1:2).ne.0d0)
         smin(1:2)=max(smin_S,this%invm_min(i+im1,1:2))
      end where
      where (this%invm_max(i+im1,1:2).ne.0d0)
         smax(1:2)=min(smax_S,this%invm_max(i+im1,1:2))
      end where
      if (im1.gt.2) then
         ! Boost and rotate in z-direction such that pp(:,im1) goes in the x-direction.
         call longitudinal_rapidity(pp(0:3,im1),y,rapidity_ok)
         if (.not.rapidity_ok) then
            ps%jac=-35d0
            return
         endif
         call boostz(pp(0,i+ir),y,piirr)
         call boostz(pp(0,im1),y,pim1r)
         call boostz(pp(0,ib),y,pib)
         phi_rot=atan2(pp(2,im1),pp(1,im1))
         call rotz(piirr,-phi_rot,piir)
         call rotz(pim1r,-phi_rot,pim1)
         ! Eir > Etmin(ir) + constraint coming from t
         delta_ir=invm(ir)-invm(ir+ib)
         delta_scale=max(spacing(0d0),abs(invm(ir)),abs(invm(ir+ib)),pib(0)**2)
         if (abs(pib(0)).le.vtiny*max(spacing(0d0),sqrtshat)) then
            ps%jac=-39d0
            if (debug) write (*,*) 'unphysical inverse Eir constraint',delta_ir,pib(0)
            return
         endif
         if (abs(delta_ir).le.vtiny*delta_scale) then
            if (delta_ir.le.0d0 .or. this%ETmin(ir,1).le.vtiny*sqrt(delta_scale)) then
               etminir(1:2)=this%ETmin(ir,1:2)
            else
               ps%jac=-45d0
               if (debug) write (*,*) 'positive singular inverse Eir constraint',delta_ir,pib(0)
               return
            endif
         else
            etminir(1:2)=max(pib(0)*this%ETmin(ir,1:2)**2/delta_ir+delta_ir/(4d0*pib(0)),&
                 this%ETmin(ir,1:2))
         endif
         rad=(piir(0)-etminir(1:2))**2-invm(i)
         if (rad(1).lt.0d0) then
            ps%jac=-46d0
            if (debug) write (*,*) 'rad.lt.0d0 inverse',rad,piir(0),etminir,invm(i)
            return
         endif
         smax(1)=min(smax(1),invm(i)+invm(im1)+2d0*(piir(0)-etminir(1))*pim1(0)+&
              2d0*sqrt(rad(1))*pim1(1))
         if (rad(2).ge.0d0) then
            smax(2)=min(smax(2),invm(i)+invm(im1)+2d0*(piir(0)-etminir(2))*pim1(0)+&
                 2d0*sqrt(rad(2))*pim1(1))
         endif
         if(invm(i).eq.0d0) then
            smin(1)=max(smin(1),2d0*this%ETmin(i,1)*(pim1(0)-pim1(1)))
            smin(2)=max(smin(2),2d0*this%ETmin(i,2)*(pim1(0)-pim1(1)*cos(this%drcut(i+im1))))
         endif
         if (rad(2).lt.0d0) smax(2)=smin(2)
      endif
      if (smin(1).ge.smax(1)) then
         ps%jac=-4d0
         return
      endif
      invm(i+im1)=dot(pp(0:3,i+im1),pp(0:3,i+im1))
      if (debug) then
         write (*,*) '23i- i+im1',i+im1,smin,smax,invm(i+im1)
      endif
      ix=ix+1
      call var_to_random(invm(i+im1),this%ip,smin,smax,ps%x(ix),ps%jac)
      if (bad_inverse_jac()) return
      ! Generate the momenta from the integration variables. Since there is an
      ! ambiguity in phi, get both of them and pick the one that passes the cuts
      ! (if it's only one). If both pass, simply pick one of the two at random
      ! with a flat prior.
      phi1=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
           &,invm(ir+i+ib),V,sqrtGG,1d0)
      call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi1 &
           &,sqrt0(invm(i)),sqrt0(invm(ir)),pi1,ppibir1,gent_status)
      if (gent_status.ne.0) then
         ps%jac=-9d0
         return
      endif
      pr1(0:3)=pp(0:3,ir+i)-pi1(0:3)
      phi2=getphifroms(invm(i+im1),invm(ir+i),invm(ir),invm(ir+i+im1)&
           &,invm(ir+i+ib),V,sqrtGG,0d0)
      call gentcms2(pp(0,ib),pp(0,ib+ir+i),pp(0,ib+ir+i+im1),invm(ir+ib),phi2 &
           &,sqrt0(invm(i)),sqrt0(invm(ir)),pi2,ppibir2,gent_status)
      if (gent_status.ne.0) then
         ps%jac=-9d0
         return
      endif
      pr2(0:3)=pp(0:3,ir+i)-pi2(0:3)
      if ( pi1(0)**2-pi1(3)**2.ge.this%ETmin(i,1)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir,1)**2 .and. &
           pi2(0)**2-pi2(3)**2.ge.this%ETmin(i,1)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir,1)**2 ) then
         continue
      elseif (pi1(0)**2-pi1(3)**2.ge.this%ETmin(i,1)**2 .and. pr1(0)**2-pr1(3)**2.ge.this%ETmin(ir,1)**2) then
         if (.not.update_inverse_jac(1d0,2d0)) return
      elseif (pi2(0)**2-pi2(3)**2.ge.this%ETmin(i,1)**2 .and. pr2(0)**2-pr2(3)**2.ge.this%ETmin(ir,1)**2) then
         if (.not.update_inverse_jac(1d0,2d0)) return
      endif
      ! Compute the Jacobian
      gram4=gram_determinant4(invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib)&
           &,invm(ir+i),invm(i+im1),invm(ir+ib+i+im1),invm(ir),invm(i)&
           &,invm(im1))
      if (gram4.ge.0d0) then 
         write (99,*) 'Warning: gram4 greater than or equal to zero in gen23_one_step_inverse',gram4,i,ir
         write (99,*) '  inverse Gram inputs i,ir,ib,im1=',i,ir,ib,im1
         write (99,*) '  shat_ip1,t_im1,t_i,shat_i,s_i,t_ip1,shat_im1,m_i_2,m_ip1_2=',&
              & invm(ir+i+im1),invm(ir+ib),invm(ir+i+ib),invm(ir+i),invm(i+im1),&
              & invm(ir+ib+i+im1),invm(ir),invm(i),invm(im1)
         ps%jac=-5d0
         return
      endif
      if (.not.update_inverse_jac(1d0,8d0*sqrt(-gram4))) return
    end subroutine gen23_one_step_inverse
    subroutine fill_momentum_array
      implicit none
      integer :: i,j
      real(kind=8),dimension(0:3) :: p
      real(kind=8) :: e_initial,pz_initial,lightcone_plus,lightcone_minus
      e_initial=ps%p(0,1)+ps%p(0,2)
      pz_initial=ps%p(3,1)+ps%p(3,2)
      lightcone_plus=e_initial+pz_initial
      lightcone_minus=e_initial-pz_initial
      if (.not.ieee_is_finite(lightcone_plus) .or. .not.ieee_is_finite(lightcone_minus)) then
         ps%jac=-44d0
         return
      endif
      if (lightcone_plus.le.0d0 .or. lightcone_minus.le.0d0) then
         ps%jac=-44d0
         return
      endif
      ycm=(log(lightcone_plus)-log(lightcone_minus))/2d0
      do i=1,maskr(this%next)
         p(0:3)=0d0
         do j=0,this%next-1
            if (btest(i,j)) then
               if (j.le.1 .and. popcnt(i).ne.1) then
                  p(0:3)=p(0:3)-ps%p(0:3,j+1)
               else
                  p(0:3)=p(0:3)+ps%p(0:3,j+1)
               endif
            endif
         enddo
         ! Keep the configured singleton masses (in particular, exact
         ! masslessness).  Compute all composite invariants before the
         ! auxiliary longitudinal boost: recomputing a nearly light-like
         ! quantity from the boosted components loses several digits through
         ! E^2-p_z^2 cancellation and can spuriously flip a Gram determinant.
         if (popcnt(i).gt.1) invm(i)=dot(p(0:3),p(0:3))
         call boostz(p(0:3),ycm,pp(0:3,i))
      enddo
    end subroutine fill_momentum_array
    subroutine var_to_random(variable,power_in,var_min,var_max,x,jac)
      ! Given a random variable var between varmin and varmax, compute
      ! the corresponding value of x between 0 and 1.
      implicit none
      real(kind=8),intent(in) :: variable
      real(kind=8),dimension(-1:1),intent(in) :: power_in
      real(kind=8),dimension(1:2),intent(in) :: var_min,var_max
      real(kind=8),intent(out) :: x
      real(kind=8),intent(inout) :: jac
      integer(kind=4) :: i,k
      real(kind=8),dimension(3) :: vmin,vmax,power,q
      real(kind=8),dimension(1:2) :: var_min_eff,var_max_eff
      real(kind=8) :: var,xloc,variable_clamped,bound_tol,range_tol,updated_jac
      logical :: found,valid_inputs
      x=0d0
      if (.not.ieee_is_finite(variable) .or. .not.ieee_is_finite(jac)) then
         jac=-41d0
         return
      endif
      if (jac.le.0d0) then
         jac=-41d0
         return
      endif
      call select_integration_bounds(var_min,var_max,var_min_eff,var_max_eff, &
           this%use_soft_bounds_as_actual_limits)
      if (var_min_eff(1).gt.var_max_eff(1)) then
         jac=-41d0
         return
      endif
      bound_tol=5d-8*max(spacing(0d0),abs(variable),abs(var_min_eff(1)),abs(var_max_eff(1)))
      if (variable.lt.var_min_eff(1)-bound_tol .or. variable.gt.var_max_eff(1)+bound_tol) then
         write (99,*) 'Warning: variable not between varmin and varmax',var_min_eff(1),variable,var_max_eff(1)
         jac=-42d0
         return
      endif
      variable_clamped=min(max(variable,var_min_eff(1)),var_max_eff(1))
      call random_to_var_inputs(power_in,var_min_eff,var_max_eff,power,vmin,vmax, &
           this%use_soft_bounds_as_actual_limits,valid_inputs)
      if (.not.valid_inputs) then
         jac=-41d0
         return
      endif
      call random_to_var_weights(power,vmin,vmax,q)
      if (any(.not.ieee_is_finite(q))) then
         jac=-41d0
         return
      endif
      if (any(q.lt.0d0) .or. sum(q).le.0d0) then
         jac=-41d0
         return
      endif
      if (var_min_eff(1).lt.0d0 .and. var_max_eff(1).le.0d0) then
         var=-variable_clamped
      else
         var=variable_clamped
      endif
      found=.false.
      k=0
      do i=1,3
         if (q(i).le.0d0) cycle
         range_tol=128d0*epsilon(1d0)*max(spacing(0d0),abs(var),abs(vmin(i)),abs(vmax(i)))
         if (var.ge.vmin(i)-range_tol .and. var.le.vmax(i)+range_tol) then
            k=i
            found=.true.
            exit
         endif
      enddo
      if (.not.found) then
         write (99,*) 'Warning: variable not in any active random-to-var range',variable,var_min_eff,var_max_eff
         jac=-43d0
         return
      endif
      call var_to_random_map(var,power(k),vmin(k),vmax(k),xloc,jac)
      x=q(k)*xloc
      if (k.gt.1) x=x+sum(q(1:k-1))
      if (.not.safe_phase_space_ratio(jac,q(k),updated_jac)) then
         x=0d0
         jac=-43d0
         return
      endif
      jac=updated_jac
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(jac)) then
         x=0d0
         jac=-43d0
         return
      endif
      if (x.lt.-vtiny .or. x.gt.1d0+vtiny .or. jac.le.0d0) then
         x=0d0
         jac=-43d0
      else
         x=max(0d0,min(1d0,x))
      endif
    end subroutine var_to_random

    subroutine var_to_random_map(var,power,varmin,varmax,x,jac)
      implicit none
      real(kind=8),intent(in) :: var,power,varmin,varmax
      real(kind=8),intent(out) :: x
      real(kind=8),intent(inout) :: jac
      real(kind=8) :: exponent,scale,ratio_power,scaled_power,log_ratio,log_scaled
      real(kind=8) :: numerator,denominator,jac_factor,interval_scale,interval_tolerance
      real(kind=8) :: updated_jac,log_jac_factor
      x=0d0
      if (.not.ieee_is_finite(var) .or. .not.ieee_is_finite(power) .or. &
           .not.ieee_is_finite(varmin) .or. .not.ieee_is_finite(varmax) .or. &
           .not.ieee_is_finite(jac)) then
         jac=-43d0
         return
      endif
      interval_scale=max(abs(varmin),abs(varmax))
      interval_tolerance=128d0*epsilon(1d0)*max(spacing(0d0),interval_scale)
      if (varmax.le.varmin) then
         jac=-43d0
         return
      endif
      if (varmax-varmin.le.interval_tolerance) then
         ! A zero-width sector is a phase-space boundary.  It has no
         ! sampling measure, so retain a finite inverse Jacobian and use its
         ! endpoint coordinate instead of classifying the point as invalid.
         x=0d0
         return
      endif
      if (power.eq.0d0) then
         x=(var-varmin)/(varmax-varmin)
         jac_factor=varmax-varmin
      elseif (power.eq.-1d0) then
         if (varmin.le.0d0 .or. var.le.0d0) then
            jac=-43d0
            return
         endif
         denominator=log(varmax)-log(varmin)
         if (.not.ieee_is_finite(denominator)) then
            jac=-43d0
            return
         endif
         if (denominator.le.0d0) then
            jac=-43d0
            return
         endif
         x=log(var/varmin)/denominator
         jac_factor=var*denominator
      elseif (power.eq.101d0) then
         x=acos(min(max(1d0-2d0*(var-varmin)/(varmax-varmin),-1d0),1d0))/pi
         jac_factor=(varmax-varmin)*pi*sin(pi*x)/2d0
      else
         if (varmin.lt.0d0 .or. varmax.le.0d0 .or. var.le.0d0) then
            jac=-43d0
            return
         endif
         exponent=1d0+power
         if (abs(exponent).le.epsilon(1d0)) then
            jac=-43d0
            return
         endif
         if (exponent.lt.0d0 .and. varmin.le.0d0) then
            jac=-43d0
            return
         endif
         if (exponent.gt.0d0) then
            scale=varmax
            log_ratio=log(varmin)-log(varmax)
            ratio_power=exp(exponent*log_ratio)
            denominator=1d0-ratio_power
            log_scaled=log(var)-log(scale)
            scaled_power=exp(exponent*log_scaled)
            numerator=scaled_power-ratio_power
         else
            scale=varmin
            log_ratio=log(varmax)-log(varmin)
            ratio_power=exp(exponent*log_ratio)
            denominator=ratio_power-1d0
            log_scaled=log(var)-log(scale)
            scaled_power=exp(exponent*log_scaled)
            numerator=scaled_power-1d0
         endif
         if (.not.ieee_is_finite(denominator)) then
            jac=-43d0
            return
         endif
         if (abs(denominator).le.128d0*epsilon(1d0)*max(spacing(0d0),abs(ratio_power),1d0)) then
            jac=-43d0
            return
         endif
         x=numerator/denominator
         log_jac_factor=log(scale)+log(abs(denominator))-&
              log(abs(exponent))-power*log_scaled
         if (.not.ieee_is_finite(log_jac_factor) .or. &
              log_jac_factor.gt.log(huge(1d0)) .or. &
              log_jac_factor.lt.log(spacing(0d0))) then
            jac=-43d0
            return
         endif
         jac_factor=exp(log_jac_factor)
      endif
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(jac_factor)) then
         x=0d0
         jac=-43d0
         return
      endif
      if (jac_factor.le.0d0) then
         x=0d0
         jac=-43d0
         return
      endif
      if (.not.safe_phase_space_product(jac,jac_factor,updated_jac)) then
         x=0d0
         jac=-43d0
         return
      endif
      jac=updated_jac
      if (jac.le.0d0) then
         x=0d0
         jac=-43d0
      endif
    end subroutine var_to_random_map
    subroutine generate_mass_inverse(i,shatmin,shatmax)
      implicit none
      integer :: i
      real(kind=8),dimension(1:2) :: shatmin,shatmax
      where (this%invm_min(i,1:2).ne.0d0)
         shatmin(1:2)=max(shatmin(1:2),this%invm_min(i,1:2))
      end where
      where (this%invm_max(i,1:2).ne.0d0)
         shatmax(1:2)=min(shatmax(1:2),this%invm_max(i,1:2))
      end where
      invm(i)=dot(pp(0:3,i),pp(0:3,i))
      if (debug) then
         write (*,*) 'mi- i',i,shatmin,shatmax,invm(i)
      endif
      ix=ix+1
      call var_to_random(invm(i),this%ip_mass,shatmin,shatmax,ps%x(ix),ps%jac)
    end subroutine generate_mass_inverse
  end subroutine gen23_compute_x_from_momenta

  real(kind=8) function dot(p1,p2)
    ! Inner product between two 4-vectors
    implicit none
    real(kind=8),intent(in),dimension(0:3) :: p1,p2
    dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
  end function dot
  subroutine shatminmax(this,j1,j2,shatmin,shatmax,invm)
    ! Determines minimum and maximum allowed s-channel invariant
    ! masses based on previously generated masses and the masses of
    ! the final-state particles.
    implicit none
    class(phase_space_gen23),intent(in) :: this
    integer(kind=4),intent(in) :: j1,j2
    real(kind=8),dimension(*),intent(in) :: invm
    real(kind=8),dimension(1:2),intent(out) :: shatmin,shatmax
    integer(kind=4) :: j
    shatmin=0d0
    do j=0,this%next-1
       if (btest(j1,j)) then
          shatmin=shatmin+sqrt0(invm(ibset(0,j)))
       endif
    enddo
    shatmin=shatmin**2
    shatmax(1:2)=(sqrt0(invm(j1+j2))-sqrt0(max(invm(j2),this%invm_min(j2,1:2))))**2
  end subroutine shatminmax
  elemental real(kind=8) function lambda(s,xa2,xb2)
    ! The usual two dimensional phase-space volume factor. See, e.g.,
    ! Eq.(A2) of E.~Byckling and K.~Kajantie, ``Reductions of the
    ! phase-space integral in terms of simpler processes,'' Phys. Rev. 187
    ! (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: xa2,xb2,S
    if (xb2.eq.0d0) then
       lambda=(s-xa2)**2
    elseif (xa2.eq.0d0) then
       lambda=(s-xb2)**2
    else
       lambda=s**2-2d0*(xa2+xb2)*s+(xa2-xb2)**2
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
    ! q^2=(q0-m)(q0+m) makes this form stable both at rest and for a
    ! nearly-rest-frame boost; it avoids division by |q_vec|^2.
    lf=(pq/(q(0)+m)+p(0))/m
    pboost(0)=(p(0)*q(0)+pq)/m
    pboost(1:3)=p(1:3)+q(1:3)*lf
    if (.not.all(ieee_is_finite(pboost))) pboost=0d0
  end subroutine boostm
  subroutine boostz(p,yb,pb)
    ! boost in the z-direction with rapidity yb
    implicit none
    real(kind=8),dimension(0:3) :: p,pb
    real(kind=8) :: yb
    pb(0)=p(0)*cosh(yb)-p(3)*sinh(yb)
    pb(1:2)=p(1:2)
    pb(3)=p(3)*cosh(yb)-p(0)*sinh(yb)
  end subroutine boostz
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
    real(kind=8),parameter :: tol=128d0*epsilon(1d0)
    real(kind=8) :: deter,matrix_scale,log_determinant
    logical :: success
    a(1:4,1)=(/ 0d0           , t_im1-shat_im1 , t_i-shat_i       , t_ip1-shat_ip1    /)
    a(1:4,2)=(/ t_im1-shat_im1, 2d0*t_im1      , t_i+t_im1-m_i_2  , t_im1+t_ip1-s_i   /)
    a(1:4,3)=(/ t_i-shat_i    , t_i+t_im1-m_i_2, 2d0*t_i          , t_i+t_ip1-m_ip1_2 /)
    a(1:4,4)=(/ t_ip1-shat_ip1, t_im1+t_ip1-s_i, t_i+t_ip1-m_ip1_2, 2d0*t_ip1         /)
    gram_determinant4=0d0
    if (.not.all(ieee_is_finite(a))) return
    matrix_scale=maxval(abs(a))
    if (matrix_scale.le.spacing(0d0)) return
    a=a/matrix_scale
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       if (.not.ieee_is_finite(deter)) return
       if (deter.eq.0d0) return
       log_determinant=log(abs(deter))+dble(n)*log(matrix_scale)-log(16d0)
       if (.not.ieee_is_finite(log_determinant)) return
       if (log_determinant.gt.log(huge(1d0))) return
       gram_determinant4=sign(exp(log_determinant),deter)
    endif
  end function gram_determinant4
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
    real(kind=8) :: t1,t2,yr,lambda_left,lambda_right
    real(kind=8) :: center,numerator,quotient,kinematic_scale
    tmin=0d0
    tmax=0d0
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(z) .or. &
         .not.ieee_is_finite(u) .or. .not.ieee_is_finite(v) .or. &
         .not.ieee_is_finite(w)) return
    center=u+w
    if (.not.ieee_is_finite(center)) return
    tmin=center
    tmax=center
    kinematic_scale=max(spacing(0d0),abs(x),abs(z),abs(u),abs(v),abs(w))
    if (abs(x).le.tiny_kin*kinematic_scale) return
    lambda_left=lambda(x,u,v)
    lambda_right=lambda(x,w,z)
    if (.not.ieee_is_finite(lambda_left) .or. &
         .not.ieee_is_finite(lambda_right)) return
    if (.not.safe_phase_space_product(lambda_left,lambda_right,yr)) return
    if (yr.le.0d0) return
    yr=sqrt(yr)
    if (.not.safe_phase_space_product(x+u-v,x+w-z,numerator)) return
    if (.not.safe_phase_space_ratio(numerator-yr,2d0*x,quotient)) return
    t1=center-quotient
    if (.not.safe_phase_space_ratio(numerator+yr,2d0*x,quotient)) return
    t2=center-quotient
    if (.not.ieee_is_finite(t1) .or. .not.ieee_is_finite(t2)) return
    tmin = min(t1,t2)
    tmax = max(t1,t2)
  end subroutine tminmax
  real(kind=8) function computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    ! Computes the determinant of V (with m_a=0) based on eq.(11) of
    ! E.~Byckling and K.~Kajantie, ``Reductions of the phase-space
    ! integral in terms of simpler processes,'' Phys. Rev. 187 (1969),
    ! 2008-2016, doi:10.1103/PhysRev.187.2008
    use LUPdecomposition
    real(kind=8),intent(in) :: shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2
    real(kind=8) :: deter,matrix_scale,log_determinant
    integer(kind=4),parameter :: n=3
    real(kind=8),dimension(n,n) :: a
    integer(kind=4),dimension(0:n) :: p
    real(kind=8),parameter :: tol=128d0*epsilon(1d0)
    logical :: success
    a(1:3,1)=(/2d0*shat_i             , shat_i-t_i    , shat_i+shat_im1-m_i_2/)
    a(1:3,2)=(/shat_i-t_i             , 0d0           , shat_im1-t_im1       /)
    a(1:3,3)=(/shat_ip1+shat_i-m_ip1_2, shat_ip1-t_ip1, 0d0                  /)
    computeV=-99d99
    if (.not.all(ieee_is_finite(a))) return
    matrix_scale=maxval(abs(a))
    if (matrix_scale.le.spacing(0d0)) return
    a=a/matrix_scale
    call LUPdecompose(a,n,tol,p,success)
    if (success) then
       call LUPdeterminant(a,p,n,deter)
       if (.not.ieee_is_finite(deter)) return
       if (deter.eq.0d0) return
       log_determinant=log(abs(deter))+dble(n)*log(matrix_scale)-log(8d0)
       if (.not.ieee_is_finite(log_determinant)) return
       if (log_determinant.gt.log(huge(1d0))) return
       computeV=-sign(exp(log_determinant),deter)
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
    real(kind=8) :: GG,s1,s2,lam,lam_scale
    V=computeV(shat_i,shat_im1,shat_ip1,t_i,t_im1,t_ip1,m_i_2,m_ip1_2)
    GG = G(t_i  , shat_ip1, shat_i  , t_ip1, m_ip1_2, 0d0) &
        *G(t_im1, shat_i  , shat_im1, t_i  , m_i_2  , 0d0)
    if (.not.ieee_is_finite(GG) .or. .not.ieee_is_finite(V)) then
       smin=shat_im1+shat_ip1
       smax=smin
       sqrtGG=0d0
       return
    endif
    if (GG.le.0d0 .or. V.eq.-99d99) then
!!$       write (*,*) 'No allowed range for s: smin=smax',GG,V
!!$       stop 1
       GG=0d0
       V=0d0
    endif
    lam=lambda(shat_i,t_i,0d0)
    lam_scale=max(spacing(0d0),abs(shat_i),abs(t_i))
    if (.not.ieee_is_finite(lam)) then
       smin=shat_im1+shat_ip1
       smax=smin
       sqrtGG=0d0
       return
    endif
    if (sqrt(abs(lam)).le.sqrt(tiny_kin)*lam_scale) then
       smin=shat_im1+shat_ip1
       smax=smin
       sqrtGG=0d0
       return
    endif
    sqrtGG=sqrt(GG)
    s1=shat_im1+shat_ip1+2d0/lam * (4d0*V + sqrtGG)
    s2=shat_im1+shat_ip1+2d0/lam * (4d0*V - sqrtGG)
    smin=min(s1,s2)
    smax=max(s1,s2)
  end subroutine sminmax
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
    real(kind=8),intent(in),dimension(0:3) :: p,q
    real(kind=8),intent(out),dimension(0:3) :: prot
    logical :: valid
    call stable_rotate_from_z_axis(p,q,prot,valid)
    if (.not.valid) prot=0d0
  end subroutine rotxxx
  subroutine rotxxx_inv(p,q,prot)
    ! Same as rotxxx, but inverse. That is, first doing
    ! rotxxx(p,q,prot) and then rotxxx_inv(prot,q,p) should give you
    ! back the original p.
    implicit none
    real(kind=8),dimension(0:3),intent(in) :: p,q
    real(kind=8),dimension(0:3),intent(out) :: prot
    logical :: valid
    call stable_rotate_to_z_axis(p,q,prot,valid)
    if (.not.valid) prot=0d0
  end subroutine rotxxx_inv
  subroutine rotz(p,phi,prot)
    implicit none
    real(kind=8),dimension(0:3) :: p,prot
    real(kind=8) :: phi
    prot(0)=p(0)
    prot(1)=p(1)*cos(phi)-p(2)*sin(phi)
    prot(2)=p(1)*sin(phi)+p(2)*cos(phi)
    prot(3)=p(3)
  end subroutine rotz
  subroutine gentcms2(pa,pb,pc,t,phi,m1,m2,p1,pr,status)
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
    integer,intent(out) :: status
    real(kind=8) :: E_acms,p_acms,esum,esum2,ed,pp2,md2,pt,pt2,pt2_tol,phi_off
    real(kind=8),dimension(0:3) :: ptot,pa_cms,ptotm,Pii,pc_cms,pc_rot,pii_rot
    real(kind=8),parameter :: tiny=1d-5
    p1=0d0
    pr=0d0
    status=0
    ptot(0:3)=pa(0:3)+pb(0:3)
    ptotm(0)=ptot(0)
    ptotm(1:3)=-ptot(1:3)
    ! determine magnitude of Pii in cms frame (from dhelas routine mom2cx)
    ESUM2 = dot(ptot,ptot)
    if (esum2 .le. 0d0) then
       status=-1
       return
    endif
    esum=sqrt(esum2)
    MD2=(M2-M1)*(M1+M2)
    ED=MD2/ESUM
    IF (M1*M2.EQ.0.d0) THEN
       PP2=0.25d0*(ESUM-ABS(ED))**2
    ELSE
       PP2=0.25d0*((MD2/ESUM)**2-2d0*(M1**2+M2**2)+ESUM**2)
       if(pp2.lt.0d0) then
          status=-2
          return
       endif
    ENDIF
    call boostm(pa,ptotm,esum,pa_cms)
    E_acms = pa_cms(0)
    p_acms = sqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)
    if (p_acms.le.vtiny*max(spacing(0d0),abs(E_acms),esum)) then
       status=-3
       return
    endif

    ! determine the offset in phi; the frame in which phi is defined
    ! is in the ptot rest-frame, with pa_cms aligned with the z-axis,
    ! and pc having zero phi angle.
    call boostm(pc,ptotm,esum,pc_cms)
    call rotxxx_inv(pc_cms,pa_cms,pc_rot)
    phi_off=atan2(pc_rot(2),pc_rot(1))

    ! define Pii in the frame where pa_cms is aligned with the positive z axis
    Pii(0) = MAX((ESUM+ED)*0.5d0,0.d0)
    if (esum+ed.le.0d0) then
       status=-4
       return
    endif
    Pii(3) = -(m2**2-t-2d0*Pii(0)*E_acms)/(2d0*p_acms)
    pt2=pp2-Pii(3)**2
    pt2_tol=tiny*max(abs(esum2),abs(pp2),abs(Pii(3)**2))
    if (pt2.lt.-pt2_tol) then
       status=-5
       return
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
  real(kind=8) function getphifroms(si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran)
    ! Given s_i (invariant mass of p_i and p_i+1, it transforms it
    ! into phi_i. Note that there are two possibilities for phi: need
    ! to pick one at random.
    ! Based on eq.(11) of E.~Byckling and K.~Kajantie, ``Reductions of
    ! the phase-space integral in terms of simpler processes,''
    ! Phys. Rev. 187 (1969), 2008-2016, doi:10.1103/PhysRev.187.2008
    implicit none
    real(kind=8),intent(in) :: si,shat_i,shat_im1,shat_ip1,t_i,V,sqrtGG,ran
    real(kind=8) :: cosphi,x,lam,lam_scale,sqrtgg_scale
    lam=lambda(shat_i,t_i,0d0)
    if (.not.ieee_is_finite(lam) .or. .not.ieee_is_finite(V) .or. &
         .not.ieee_is_finite(sqrtGG)) then
       getphifroms=0d0
       return
    endif
    lam_scale=max(spacing(0d0),abs(shat_i),abs(t_i))
    if (sqrt(abs(lam)).le.sqrt(tiny_kin)*lam_scale) then
       getphifroms=0d0
       return
    endif
    sqrtgg_scale=max(spacing(0d0),abs((si-shat_im1-shat_ip1)*0.5d0*lam),abs(4d0*V))
    if (.not.ieee_is_finite(sqrtgg_scale)) then
       getphifroms=0d0
       return
    endif
    if (abs(sqrtGG).le.tiny_kin*sqrtgg_scale) then
       getphifroms=0d0
       return
    endif
    cosphi=((si-shat_im1-shat_ip1)*0.5d0*lam-4d0*V)/sqrtGG
    if (cosphi.lt.-1d0 .or. cosphi.gt.1d0) then
       write (99,*) 'WARNING cosphi does not have a reasonable value',cosphi
       getphifroms=0d0
       return
    endif
    x=ran
    if (x.gt.0.5d0) then
       getphifroms=acos(cosphi)
    else
       getphifroms=-acos(cosphi)+2d0*pi
    endif
  end function getphifroms

end module phase_space_gen23_mod
