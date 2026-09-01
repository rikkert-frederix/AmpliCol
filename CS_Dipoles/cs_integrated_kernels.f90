module cs_integrated_kernels
  ! Universal massless pieces of the Catani--Seymour I, P and K
  ! insertion operators.  Distribution-valued kernels are returned as
  !
  !   regular(x) + plus_one/(1-x)_+
  !              + plus_log*(log((1-x)/x)/(1-x))_+
  !              + plus_log_one*(log(1-x)/(1-x))_+
  !              + delta*delta(1-x).
  !
  ! Process-dependent leading-colour correlations are deliberately not
  ! included here.  They are supplied by the integrated-history registry so
  ! that the same ordered histories and normalisations as the local dipoles
  ! are used.
  use cs_dipole_mappings, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  real(dp), parameter, public :: cs_pi = 3.1415926535897932384626433832795_dp
  real(dp), parameter, public :: cs_ca = 3.0_dp
  real(dp), parameter, public :: cs_cf_lc = cs_ca/2.0_dp
  ! For a real incoming quark reducing to a gluon Born leg, the local
  ! g->q qbar kernel contains 2*T_R and an ordered-history weight of 1/2.
  ! Changing the incoming average from the real quark to the reduced Born
  ! gluon supplies (N_c^2-1)/N_c=8/3.  The history weight is applied by
  ! integrated_beam, leaving the physical C_F=4/3 after both factors.
  real(dp), parameter, public :: cs_cf_initial_qg = &
       (cs_ca*cs_ca-1.0_dp)/cs_ca
  ! A q(real)->U(1)(Born) history has the same initial-state average as the
  ! physical-gluon parent, while its local splitting trace is smaller by
  ! one power of N_c.
  real(dp), parameter, public :: cs_u1_initial_factor = 1.0_dp/cs_ca
  real(dp), parameter, public :: cs_tr = 0.5_dp
  ! The local crossed q <- g kernel is normalized with C_F^LC, while the
  ! real and reduced-Born matrix elements carry, respectively, gluon and
  ! quark initial-state averages.  Expressed as a kernel multiplying the
  ! averaged Born contribution this gives
  !
  !   C_F^LC * N_c/(N_c^2-1) = 9/16
  !
  ! for N_c=3.  Using T_R here mixes the full-colour AP convention with the
  ! leading-colour local dipole convention and leaves a 9/8 mismatch.
  real(dp), parameter, public :: cs_tr_initial_lc = &
       cs_cf_lc*cs_ca/(cs_ca*cs_ca-1.0_dp)

  integer, parameter, public :: cs_parton_q = 1
  integer, parameter, public :: cs_parton_g = 2
  integer, parameter, public :: cs_scheme_hv = 1
  integer, parameter, public :: cs_scheme_fdh = 2
  real(dp), parameter :: kernel_coordinate_floor=sqrt(tiny(1.0_dp))

  type, public :: cs_distribution
     real(dp) :: regular = 0.0_dp
     real(dp) :: plus_one = 0.0_dp
     real(dp) :: plus_log = 0.0_dp
     real(dp) :: plus_log_one = 0.0_dp
     real(dp) :: delta = 0.0_dp
  end type cs_distribution

  public :: cs_gamma, cs_k_constant, cs_i_qg, cs_i_gg, cs_i_gg_ordered
  public :: cs_i_qqbar
  public :: cs_ap_distribution, cs_kbar_distribution
  public :: cs_ff_alpha_endpoint, cs_fi_alpha_endpoint
  public :: cs_fi_distribution, cs_fi_alpha_terms
  public :: cs_if_tilde_distribution, cs_ii_tilde_distribution
  public :: cs_if_alpha_distribution, cs_ii_alpha_distribution
  public :: cs_scale_distribution, cs_fdh_endpoint_shift

contains

  pure real(dp) function cs_gamma(parton,nf) result(value)
    integer, intent(in) :: parton,nf
    select case (parton)
    case (cs_parton_q)
       value=1.5_dp*cs_cf_lc
    case (cs_parton_g)
       value=(11.0_dp/6.0_dp)*cs_ca-(2.0_dp/3.0_dp)*cs_tr*real(nf,dp)
    case default
       value=0.0_dp
    end select
  end function cs_gamma

  pure real(dp) function cs_k_constant(parton,nf) result(value)
    integer, intent(in) :: parton,nf
    select case (parton)
    case (cs_parton_q)
       value=(3.5_dp-cs_pi**2/6.0_dp)*cs_cf_lc
    case (cs_parton_g)
       value=(67.0_dp/18.0_dp-cs_pi**2/6.0_dp)*cs_ca &
            -(10.0_dp/9.0_dp)*cs_tr*real(nf,dp)
    case default
       value=0.0_dp
    end select
  end function cs_k_constant

  pure subroutine cs_i_qg(coeff)
    ! Integral of one q -> q g splitting kernel in the HV scheme.
    real(dp), intent(out) :: coeff(-2:0)
    coeff(-2)=cs_cf_lc
    coeff(-1)=1.5_dp*cs_cf_lc
    coeff(0)=cs_cf_lc*(5.0_dp-cs_pi**2/2.0_dp)
  end subroutine cs_i_qg

  pure subroutine cs_i_gg(coeff)
    ! Integral of the symmetric g -> g g kernel.  An ordered local history
    ! must use cs_i_gg_ordered instead; the history's lc_weight is an
    ! independent colour-correlation factor.
    real(dp), intent(out) :: coeff(-2:0)
    coeff(-2)=2.0_dp*cs_ca
    coeff(-1)=(11.0_dp/3.0_dp)*cs_ca
    coeff(0)=2.0_dp*cs_ca*(50.0_dp/9.0_dp-cs_pi**2/2.0_dp)
  end subroutine cs_i_gg

  pure subroutine cs_i_gg_ordered(coeff)
    ! Integral of one ordered side of the g -> g g kernel.  Every local
    ! leading-colour history designates one colour neighbour and carries
    ! only that ordered side, for final- and initial-state emitters alike.
    ! The complementary soft pole and the other half of the collinear
    ! kernel belong to the history on the other side of the gluon.
    real(dp), intent(out) :: coeff(-2:0)
    call cs_i_gg(coeff)
    coeff=0.5_dp*coeff
  end subroutine cs_i_gg_ordered

  pure subroutine cs_i_qqbar(coeff)
    ! Integral of one g -> q qbar flavour.
    real(dp), intent(out) :: coeff(-2:0)
    coeff(-2)=0.0_dp
    coeff(-1)=-(2.0_dp/3.0_dp)*cs_tr
    coeff(0)=-(16.0_dp/9.0_dp)*cs_tr
  end subroutine cs_i_qqbar

  pure subroutine cs_ap_distribution(a,b,x,nf,kernel,info)
    ! Four-dimensional regularised Altarelli--Parisi probability P^{a b}.
    ! a is the flavour entering the incoming PDF before the splitting and b
    ! is the reduced Born incoming flavour.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info

    kernel=cs_distribution()
    info=0
    if (.not.ieee_is_finite(x)) then
       info=-20
       return
    endif
    if (nf.lt.0 .or. nf.gt.6) then
       info=-1
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif

    if (a == cs_parton_q .and. b == cs_parton_g) then
       kernel%regular=cs_cf_initial_qg*(1.0_dp+(1.0_dp-x)**2)/x
    else if (a == cs_parton_g .and. b == cs_parton_q) then
       kernel%regular=cs_tr_initial_lc*(x*x+(1.0_dp-x)**2)
    else if (a == cs_parton_q .and. b == cs_parton_q) then
       kernel%regular=-cs_cf_lc*(1.0_dp+x)
       kernel%plus_one=2.0_dp*cs_cf_lc
       kernel%delta=cs_gamma(cs_parton_q,nf)
    else if (a == cs_parton_g .and. b == cs_parton_g) then
       kernel%regular=2.0_dp*cs_ca*((1.0_dp-x)/x-1.0_dp+x*(1.0_dp-x))
       kernel%plus_one=2.0_dp*cs_ca
       kernel%delta=cs_gamma(cs_parton_g,nf)
    else
       info=-2
    endif
    if (info.eq.0 .and. .not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_ap_distribution

  pure subroutine cs_kbar_distribution(a,b,x,nf,kernel,info)
    ! MSbar finite flavour kernel Kbar^{a b}.  The factorisation-scheme
    ! kernel K_FS is zero and the colour-correlated gamma contribution is
    ! added by the caller.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info
    type(cs_distribution) :: ap
    real(dp) :: logarithm

    kernel=cs_distribution()
    info=0
    if (.not.ieee_is_finite(x)) then
       info=-20
       return
    endif
    if (nf.lt.0 .or. nf.gt.6) then
       info=-1
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    logarithm=log(1.0_dp-x)-log(x)

    if (a == cs_parton_q .and. b == cs_parton_g) then
       call cs_ap_distribution(a,b,x,nf,ap,info)
       if (info /= 0) return
       kernel%regular=ap%regular*logarithm+cs_cf_initial_qg*x
    else if (a == cs_parton_g .and. b == cs_parton_q) then
       call cs_ap_distribution(a,b,x,nf,ap,info)
       if (info /= 0) return
       kernel%regular=ap%regular*logarithm+2.0_dp*cs_tr_initial_lc*x*(1.0_dp-x)
    else if (a == cs_parton_q .and. b == cs_parton_q) then
       kernel%regular=cs_cf_lc*(-(1.0_dp+x)*logarithm+(1.0_dp-x))
       kernel%plus_log=2.0_dp*cs_cf_lc
       kernel%delta=-(5.0_dp-cs_pi**2)*cs_cf_lc
    else if (a == cs_parton_g .and. b == cs_parton_g) then
       kernel%regular=2.0_dp*cs_ca*((1.0_dp-x)/x-1.0_dp+x*(1.0_dp-x))*logarithm
       kernel%plus_log=2.0_dp*cs_ca
       kernel%delta=-(50.0_dp/9.0_dp-cs_pi**2)*cs_ca &
            +(16.0_dp/9.0_dp)*cs_tr*real(nf,dp)
    else
       info=-2
    endif
    if (info.eq.0 .and. .not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_kbar_distribution

  pure subroutine cs_ff_alpha_endpoint(alpha,primitive,correction,info)
    ! Finite change of a massless final--final endpoint relative to alpha=1.
    ! Writing it in terms of the Laurent primitive covers qg, gg and qqbar
    ! histories without independently duplicating their colour factors.
    real(dp), intent(in) :: alpha,primitive(-2:0)
    real(dp), intent(out) :: correction
    integer, intent(out) :: info
    real(dp) :: logarithm

    correction=0.0_dp
    info=0
    if (.not.ieee_is_finite(alpha) .or. .not.all(ieee_is_finite(primitive))) then
       info=-20
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (alpha >= 1.0_dp) return

    logarithm=log(alpha)
    correction=primitive(-1)*(alpha-1.0_dp-logarithm) &
         -primitive(-2)*logarithm*logarithm
    if (.not.ieee_is_finite(correction)) then
       correction=0.0_dp
       info=-20
    endif
  end subroutine cs_ff_alpha_endpoint

  pure subroutine cs_fi_alpha_endpoint(alpha,primitive,correction,info)
    ! Delta(1-x) part of the finite massless final--initial alpha change.
    ! The remaining x-dependent contribution is returned by
    ! cs_fi_alpha_terms.
    real(dp), intent(in) :: alpha,primitive(-2:0)
    real(dp), intent(out) :: correction
    integer, intent(out) :: info
    real(dp) :: logarithm

    correction=0.0_dp
    info=0
    if (.not.ieee_is_finite(alpha) .or. .not.all(ieee_is_finite(primitive))) then
       info=-20
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (alpha >= 1.0_dp) return

    logarithm=log(alpha)
    correction=-primitive(-1)*logarithm &
         -primitive(-2)*logarithm*logarithm
    if (.not.ieee_is_finite(correction)) then
       correction=0.0_dp
       info=-20
    endif
  end subroutine cs_fi_alpha_endpoint

  pure subroutine cs_fi_distribution(primitive,x,alpha,regular_gz,subtracted,info)
    ! Complete non-endpoint massless final--initial distribution.  Its action
    ! on a test function g is
    !
    !   regular_gz*g(x) + subtracted*(g(x)-g(1)).
    !
    ! At alpha=1 this is the finite FI baseline.  Restricting the local
    ! dipole to 1-x < alpha leaves only x > 1-alpha.
    real(dp), intent(in) :: primitive(-2:0),x,alpha
    real(dp), intent(out) :: regular_gz,subtracted
    integer, intent(out) :: info
    real(dp) :: one_minus_x

    regular_gz=0.0_dp
    subtracted=0.0_dp
    info=0
    if (.not.all(ieee_is_finite(primitive)) .or. .not.ieee_is_finite(x) .or. &
         .not.ieee_is_finite(alpha)) then
       info=-20
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (x <= 1.0_dp-alpha) return

    one_minus_x=1.0_dp-x
    regular_gz=2.0_dp*primitive(-2)*log(2.0_dp-x)/one_minus_x
    subtracted=-(2.0_dp*primitive(-2)*log(one_minus_x)+primitive(-1)) &
         /one_minus_x
    if (.not.ieee_is_finite(regular_gz) .or. .not.ieee_is_finite(subtracted)) then
       regular_gz=0.0_dp
       subtracted=0.0_dp
       info=-20
    endif
  end subroutine cs_fi_distribution

  pure subroutine cs_fi_alpha_terms(primitive,x,alpha,regular_gz,subtracted,info)
    ! Non-endpoint final--initial alpha change.  Its action on a test
    ! function g is
    !
    !   regular_gz*g(x) + subtracted*(g(x)-g(1)).
    !
    ! This is the distributional form obtained by changing the upper
    ! y_{ij,a} limit from one to alpha.  It is nonzero only for
    ! x < 1-alpha, exactly as the corresponding restricted local dipole.
    real(dp), intent(in) :: primitive(-2:0),x,alpha
    real(dp), intent(out) :: regular_gz,subtracted
    integer, intent(out) :: info
    real(dp) :: one_minus_x

    regular_gz=0.0_dp
    subtracted=0.0_dp
    info=0
    if (.not.all(ieee_is_finite(primitive)) .or. .not.ieee_is_finite(x) .or. &
         .not.ieee_is_finite(alpha)) then
       info=-20
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (alpha >= 1.0_dp .or. x >= 1.0_dp-alpha) return

    one_minus_x=1.0_dp-x
    regular_gz=-2.0_dp*primitive(-2)*log(2.0_dp-x)/one_minus_x
    subtracted=(2.0_dp*primitive(-2)*log(one_minus_x)+primitive(-1)) &
         /one_minus_x
    if (.not.ieee_is_finite(regular_gz) .or. .not.ieee_is_finite(subtracted)) then
       regular_gz=0.0_dp
       subtracted=0.0_dp
       info=-20
    endif
  end subroutine cs_fi_alpha_terms

  pure subroutine cs_if_tilde_distribution(a,b,x,nf,kernel,info)
    ! Spectator-dependent finite K contribution for a massless
    ! initial--final dipole at alpha=1.  It is nonzero only for diagonal
    ! splittings.  The alpha-dependent change is returned separately by
    ! cs_if_alpha_distribution.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info
    type(cs_distribution) :: ap
    real(dp) :: one_minus_x

    kernel=cs_distribution()
    info=0
    call cs_ap_distribution(a,b,x,nf,ap,info)
    if (info /= 0) return
    if (a /= b) return

    one_minus_x=1.0_dp-x
    kernel%regular=-ap%plus_one*log(2.0_dp-x)/one_minus_x
    kernel%plus_log_one=ap%plus_one
    kernel%delta=-(cs_pi**2/6.0_dp)*ap%plus_one
    if (.not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_if_tilde_distribution

  pure subroutine cs_ii_tilde_distribution(a,b,x,nf,kernel,info)
    ! Spectator-dependent finite K contribution for a massless
    ! initial--initial dipole at alpha=1.  In distributional form this is
    ! the four-dimensional AP kernel multiplied by log(1-x), together with
    ! the diagonal endpoint constant.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info
    type(cs_distribution) :: ap

    kernel=cs_distribution()
    info=0
    call cs_ap_distribution(a,b,x,nf,ap,info)
    if (info /= 0) return

    kernel%regular=ap%regular*log(1.0_dp-x)
    kernel%plus_log_one=ap%plus_one
    if (a == b) then
       kernel%delta=-(cs_pi**2/3.0_dp)*ap%plus_one
    endif
    if (.not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_ii_tilde_distribution

  pure subroutine cs_if_alpha_distribution(a,b,x,nf,alpha,kernel,info)
    ! Finite change in a massless initial--final integrated dipole relative
    ! to alpha=1.  Unlike final--initial dipoles this topology has no
    ! alpha-dependent endpoint term.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x,alpha
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info
    type(cs_distribution) :: ap
    real(dp) :: one_minus_x

    kernel=cs_distribution()
    info=0
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(alpha)) then
       info=-20
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif

    call cs_ap_distribution(a,b,x,nf,ap,info)
    if (info /= 0 .or. alpha >= 1.0_dp) return

    one_minus_x=1.0_dp-x
    kernel%regular=(ap%regular+ap%plus_one/one_minus_x)*log(alpha)
    if (ap%plus_one /= 0.0_dp) then
       kernel%regular=kernel%regular-ap%plus_one/one_minus_x*&
            log((one_minus_x+alpha)/(one_minus_x+1.0_dp))
    endif
    if (.not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_if_alpha_distribution

  pure subroutine cs_ii_alpha_distribution(a,b,x,nf,alpha,kernel,info)
    ! Finite change in an integrated initial--initial dipole when the local
    ! subtraction is restricted by
    !
    !   v_j = p_a.p_j / p_a.p_b < alpha .
    !
    ! The unresolved phase-space boundary is v_j < 1-x.  Relative to
    ! alpha=1, integration therefore adds the unregularised four-dimensional
    ! splitting kernel times log(alpha/(1-x)) for x < 1-alpha.  Expressing
    ! it through cs_ap_distribution keeps every flavour and leading-colour
    ! normalization identical to P and to the local histories.
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: x,alpha
    type(cs_distribution), intent(out) :: kernel
    integer, intent(out) :: info
    type(cs_distribution) :: ap

    kernel=cs_distribution()
    info=0
    if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(alpha)) then
       info=-20
       return
    endif
    if (x < kernel_coordinate_floor .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif

    call cs_ap_distribution(a,b,x,nf,ap,info)
    if (info /= 0) return
    if (alpha >= 1.0_dp .or. x >= 1.0_dp-alpha) return

    kernel%regular=(ap%regular+ap%plus_one/(1.0_dp-x))*&
         (log(alpha)-log(1.0_dp-x))
    if (.not.distribution_is_finite(kernel)) then
       kernel=cs_distribution()
       info=-20
    endif
  end subroutine cs_ii_alpha_distribution

  pure real(dp) function cs_fdh_endpoint_shift(parton) result(value)
    ! V_I^FDH = V_I^HV - tilde(gamma)_I.  Massive emitters have no
    ! scheme shift and are filtered by the caller.
    integer, intent(in) :: parton
    select case (parton)
    case (cs_parton_q)
       value=-0.5_dp*cs_cf_lc
    case (cs_parton_g)
       value=-cs_ca/6.0_dp
    case default
       value=0.0_dp
    end select
  end function cs_fdh_endpoint_shift

  pure subroutine cs_scale_distribution(kernel,factor)
    type(cs_distribution), intent(inout) :: kernel
    real(dp), intent(in) :: factor
    kernel%regular=factor*kernel%regular
    kernel%plus_one=factor*kernel%plus_one
    kernel%plus_log=factor*kernel%plus_log
    kernel%plus_log_one=factor*kernel%plus_log_one
    kernel%delta=factor*kernel%delta
  end subroutine cs_scale_distribution

  pure logical function distribution_is_finite(kernel)
    type(cs_distribution),intent(in) :: kernel
    distribution_is_finite=all(ieee_is_finite([kernel%regular,kernel%plus_one,&
         kernel%plus_log,kernel%plus_log_one,kernel%delta]))
  end function distribution_is_finite

end module cs_integrated_kernels
