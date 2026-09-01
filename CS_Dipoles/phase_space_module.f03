module phase_space_module
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter, public :: dp = selected_real_kind(15, 307)
  integer, parameter, public :: PS_OK = 0
  integer, parameter, public :: PS_BAD_INPUT = 1
  integer, parameter, public :: PS_NO_PHASE_SPACE = 2
  integer, parameter, public :: PS_NUMERIC_FAIL = 3

  ! Four-vector component indices: p(0:3,i) = (E, px, py, pz).
  integer, parameter, public :: P_E = 0
  integer, parameter, public :: P_X = 1
  integer, parameter, public :: P_Y = 2
  integer, parameter, public :: P_Z = 3

  ! Full scattering event column indices.
  ! Columns 1 and 2 are incoming massless particles.
  ! Columns 3..nout+2 are outgoing particles 1..nout.
  integer, parameter, public :: COL_IN1 = 1
  integer, parameter, public :: COL_IN2 = 2
  integer, parameter, public :: COL_FIRST_OUT = 3

  real(dp), parameter,private :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter, private :: ps_quartic_input_limit = &
       0.125_dp*sqrt(sqrt(huge(1.0_dp)))
  real(dp), parameter, private :: validation_rel = 1.0e-8_dp

  public :: generate_phase_space
  public :: generate_scattering_phase_space
  public :: generate_cm_scattering_phase_space
  public :: soft_deform
  public :: soft_deform_event
  public :: collinear_deform
  public :: collinear_deform_event
  public :: event_residuals
  public :: scattering_event_residuals
  public :: incoming_total
  public :: total_momentum
  public :: dot3, dot4, msq4, invariant_mass, sij

contains

  logical function finite_bounded_scalar(x)
    real(dp), intent(in) :: x
    finite_bounded_scalar = .false.
    if (.not.ieee_is_finite(x)) return
    finite_bounded_scalar = abs(x) <= ps_quartic_input_limit
  end function finite_bounded_scalar

  logical function finite_bounded_vector(x)
    real(dp), intent(in) :: x(:)
    finite_bounded_vector = all(ieee_is_finite(x))
    if (.not.finite_bounded_vector) return
    if (size(x) > 0) finite_bounded_vector = maxval(abs(x)) <= ps_quartic_input_limit
  end function finite_bounded_vector

  function norm3(v) result(r)
    real(dp), intent(in) :: v(3)
    real(dp) :: r, scale

    if (.not.all(ieee_is_finite(v))) then
       r = huge(1.0_dp)
       return
    end if
    scale = maxval(abs(v))
    if (scale == 0.0_dp) then
       r = 0.0_dp
    else if (.not.ieee_is_finite(scale)) then
       r = huge(1.0_dp)
    else
       r = scale*sqrt(sum((v/scale)**2))
    end if
  end function norm3

  function hypot2(a, b) result(r)
    real(dp), intent(in) :: a, b
    real(dp) :: r, scale

    if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b)) then
       r = huge(1.0_dp)
       return
    end if
    scale = max(abs(a), abs(b))
    if (scale == 0.0_dp) then
       r = 0.0_dp
    else if (.not.ieee_is_finite(scale)) then
       r = huge(1.0_dp)
    else
       r = scale*sqrt((a/scale)**2 + (b/scale)**2)
    end if
  end function hypot2

  function dot3(a, b) result(d)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: d
    d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function dot3

  function dot4(a, b) result(d)
    real(dp), intent(in) :: a(0:3), b(0:3)
    real(dp) :: d
    d = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function dot4

  function msq4(p) result(m2)
    real(dp), intent(in) :: p(0:3)
    real(dp) :: m2, scale, pn(0:3)

    if (.not.all(ieee_is_finite(p))) then
       m2 = huge(1.0_dp)
       return
    end if
    scale = maxval(abs(p))
    if (scale == 0.0_dp) then
       m2 = 0.0_dp
    else if (scale > sqrt(huge(1.0_dp))) then
       m2 = huge(1.0_dp)
    else
       pn = p/scale
       m2 = scale*scale*(pn(0)*pn(0) - sum(pn(1:3)*pn(1:3)))
    end if
  end function msq4

  function invariant_mass(p) result(m)
    real(dp), intent(in) :: p(0:3)
    real(dp) :: m, scale, pn(0:3), m2n

    if (.not.all(ieee_is_finite(p))) then
       m = huge(1.0_dp)
       return
    end if
    scale = maxval(abs(p))
    if (scale == 0.0_dp) then
       m = 0.0_dp
    else if (.not.ieee_is_finite(scale)) then
       m = huge(1.0_dp)
    else
       pn = p/scale
       m2n = pn(0)*pn(0) - sum(pn(1:3)*pn(1:3))
       m = scale*sqrt(max(0.0_dp, m2n))
    end if
  end function invariant_mass

  function normalized_msq4(p) result(m2n)
    real(dp), intent(in) :: p(0:3)
    real(dp) :: m2n, scale, pn(0:3)

    if (.not.all(ieee_is_finite(p))) then
       m2n = -huge(1.0_dp)
       return
    end if
    scale = maxval(abs(p))
    if (scale == 0.0_dp) then
       m2n = 0.0_dp
    else if (.not.ieee_is_finite(scale)) then
       m2n = -huge(1.0_dp)
    else
       pn = p/scale
       m2n = pn(0)*pn(0) - sum(pn(1:3)*pn(1:3))
    end if
  end function normalized_msq4

  function sij(p, i, j) result(sout)
    real(dp), intent(in) :: p(0:,:)
    integer, intent(in) :: i, j
    real(dp) :: sout
    real(dp) :: q(0:3)
    q = p(0:3, i) + p(0:3, j)
    sout = msq4(q)
  end function sij

  function safe_sqrt(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x <= 0.0_dp) then
       y = 0.0_dp
    else
       y = sqrt(x)
    end if
  end function safe_sqrt

  function kallen(a, b, c) result(lam)
    real(dp), intent(in) :: a, b, c
    real(dp) :: lam, scale, an, bn, cn

    scale = max(abs(a), abs(b), abs(c))
    if (scale == 0.0_dp) then
       lam = 0.0_dp
       return
    end if
    an = a/scale
    bn = b/scale
    cn = c/scale
    lam = scale*scale*(an*an + bn*bn + cn*cn - &
         2.0_dp*(an*bn + an*cn + bn*cn))
  end function kallen

  subroutine set_status(status, value)
    integer, intent(out), optional :: status
    integer, intent(in) :: value
    if (present(status)) status = value
  end subroutine set_status

  subroutine total_momentum(n, p, psum)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(0:3,n)
    real(dp), intent(out) :: psum(0:3)
    integer :: i
    psum = 0.0_dp
    do i = 1, n
       psum = psum + p(:, i)
    end do
  end subroutine total_momentum

  subroutine incoming_total(p, ptot)
    real(dp), intent(in) :: p(0:3,2)
    real(dp), intent(out) :: ptot(0:3)
    ptot = p(:, 1) + p(:, 2)
  end subroutine incoming_total

  function max_abs4(p) result(r)
    real(dp), intent(in) :: p(0:3)
    real(dp) :: r
    r = maxval(abs(p))
  end function max_abs4

  logical function masses_are_valid(n, mass)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n)

    masses_are_valid = n >= 1
    if (.not.masses_are_valid) return
    masses_are_valid = all(ieee_is_finite(mass))
    if (.not.masses_are_valid) return
    masses_are_valid = all(mass >= 0.0_dp) .and. maxval(mass) <= ps_quartic_input_limit
  end function masses_are_valid

  logical function final_state_is_valid(n, mass, p, total)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n), p(0:3,n), total(0:3)
    integer :: i
    real(dp) :: scale, shell_scale, shell_tol, momentum_tol, psum(0:3)

    final_state_is_valid = .false.
    if (.not.masses_are_valid(n, mass) .or. .not.all(ieee_is_finite(p)) .or. &
         .not.finite_bounded_vector(total)) return
    if (maxval(abs(p)) > ps_quartic_input_limit) return

    scale = max(maxval(abs(p)), maxval(mass), maxval(abs(total)), tiny(1.0_dp))
    momentum_tol = validation_rel*scale*real(max(1, n), dp)
    psum = 0.0_dp
    do i = 1, n
       if (p(0,i) < 0.0_dp .or. (mass(i) > 0.0_dp .and. p(0,i) <= 0.0_dp)) return
       shell_scale = max(max_abs4(p(:,i)), mass(i), tiny(1.0_dp))
       shell_tol = validation_rel*shell_scale*shell_scale
       if (abs(msq4(p(:,i)) - mass(i)*mass(i)) > shell_tol) return
       psum = psum + p(:,i)
       if (.not.finite_bounded_vector(psum)) return
    end do
    if (maxval(abs(psum - total)) > momentum_tol) return
    final_state_is_valid = .true.
  end function final_state_is_valid

  logical function scattering_event_is_valid(nout, mass, p)
    integer, intent(in) :: nout
    real(dp), intent(in) :: mass(nout), p(0:3,nout+2)
    integer :: i
    real(dp) :: incoming(0:3), beam_scale, shell_tol

    scattering_event_is_valid = .false.
    if (nout < 1) return
    if (.not.all(ieee_is_finite(p))) return
    if (maxval(abs(p)) > ps_quartic_input_limit) return
    if (p(0,1) <= 0.0_dp .or. p(0,2) <= 0.0_dp) return
    do i = 1, 2
       beam_scale = max(max_abs4(p(:,i)), tiny(1.0_dp))
       shell_tol = validation_rel*beam_scale*beam_scale
       if (abs(msq4(p(:,i))) > shell_tol) return
    end do
    incoming = p(:,1) + p(:,2)
    if (incoming(0) <= 0.0_dp .or. normalized_msq4(incoming) <= 0.0_dp) return
    scattering_event_is_valid = final_state_is_valid(nout, mass, p(:,3:nout+2), incoming)
  end function scattering_event_is_valid

  subroutine boost_one(p, beta, pb, istat)
    real(dp), intent(in) :: p(0:3)
    real(dp), intent(in) :: beta(3)
    real(dp), intent(out) :: pb(0:3)
    integer, intent(out) :: istat
    real(dp) :: b2, gamma, bp, fac

    pb = 0.0_dp
    if (.not.finite_bounded_vector(p) .or. .not.all(ieee_is_finite(beta))) then
       istat = PS_BAD_INPUT
       return
    end if
    if (maxval(abs(beta)) > 1.0_dp) then
       istat = PS_BAD_INPUT
       return
    end if

    b2 = dot3(beta, beta)
    bp = beta(1)*p(1) + beta(2)*p(2) + beta(3)*p(3)
    if (.not.ieee_is_finite(b2)) then
       istat = PS_BAD_INPUT
       return
    end if
    if (b2 < 0.0_dp) then
       istat = PS_BAD_INPUT
       return
    end if
    if (b2 <= epsilon(1.0_dp)**2) then
       pb(0) = p(0) + bp
       pb(1) = p(1) + p(0)*beta(1)
       pb(2) = p(2) + p(0)*beta(2)
       pb(3) = p(3) + p(0)*beta(3)
       if (.not.finite_bounded_vector(pb)) then
          pb = 0.0_dp
          istat = PS_NUMERIC_FAIL
          return
       end if
       istat = PS_OK
       return
    end if

    if (b2 >= 1.0_dp - 10.0_dp*epsilon(1.0_dp)) then
       pb = 0.0_dp
       istat = PS_BAD_INPUT
       return
    end if

    gamma = 1.0_dp / sqrt(1.0_dp - b2)
    fac = ((gamma - 1.0_dp)*bp/b2 + gamma*p(0))

    pb(0) = gamma*(p(0) + bp)
    pb(1) = p(1) + fac*beta(1)
    pb(2) = p(2) + fac*beta(2)
    pb(3) = p(3) + fac*beta(3)
    if (.not.finite_bounded_vector(pb)) then
       pb = 0.0_dp
       istat = PS_NUMERIC_FAIL
       return
    end if
    istat = PS_OK
  end subroutine boost_one

  subroutine boost_many(n, p, beta, pb, istat)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(0:3,n)
    real(dp), intent(in) :: beta(3)
    real(dp), intent(out) :: pb(0:3,n)
    integer, intent(out) :: istat
    integer :: i, st

    pb = 0.0_dp
    istat = PS_OK
    do i = 1, n
       call boost_one(p(:, i), beta, pb(:, i), st)
       if (st /= PS_OK) then
          pb = 0.0_dp
          istat = st
          return
       end if
    end do
  end subroutine boost_many

  function scaled_energy(n, mass, kvec, scale) result(e)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: kvec(3,n)
    real(dp), intent(in) :: scale
    real(dp) :: e, kmag, scaled_k, term
    integer :: i

    e = 0.0_dp
    do i = 1, n
       kmag = norm3(kvec(:, i))
       if (abs(scale) > 1.0_dp) then
          if (kmag > ps_quartic_input_limit/abs(scale)) then
             e = huge(1.0_dp)
             return
          end if
       end if
       scaled_k = scale*kmag
       term = hypot2(mass(i), scaled_k)
       if (.not.ieee_is_finite(term)) then
          e = huge(1.0_dp)
          return
       end if
       if (term > huge(1.0_dp) - e) then
          e = huge(1.0_dp)
          return
       end if
       e = e + term
    end do
  end function scaled_energy

  subroutine solve_spatial_scale(n, mass, kvec, target_energy, scale, istat)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: kvec(3,n)
    real(dp), intent(in) :: target_energy
    real(dp), intent(out) :: scale
    integer, intent(out) :: istat

    integer :: iter, i
    real(dp) :: sum_mass, max_knorm, lo, hi, mid, fmid, fhi, tol

    scale = 0.0_dp
    istat = PS_OK

    if (n < 1 .or. .not.finite_bounded_scalar(target_energy) .or. &
         .not.all(ieee_is_finite(mass)) .or. .not.all(ieee_is_finite(kvec))) then
       istat = PS_BAD_INPUT
       return
    end if
    if (target_energy < 0.0_dp .or. maxval(abs(mass)) > ps_quartic_input_limit .or. &
         maxval(abs(kvec)) > ps_quartic_input_limit) then
       istat = PS_BAD_INPUT
       return
    end if

    sum_mass = 0.0_dp
    max_knorm = 0.0_dp
    do i = 1, n
       if (mass(i) < 0.0_dp) then
          istat = PS_BAD_INPUT
          return
       end if
       if (mass(i) > ps_quartic_input_limit - sum_mass) then
          istat = PS_BAD_INPUT
          return
       end if
       sum_mass = sum_mass + mass(i)
       max_knorm = max(max_knorm, norm3(kvec(:, i)))
    end do

    tol = 1.0e-12_dp * max(abs(target_energy), sum_mass, tiny(1.0_dp))
    if (target_energy < sum_mass - tol) then
       istat = PS_NO_PHASE_SPACE
       return
    end if

    if (abs(target_energy - sum_mass) <= tol) then
       scale = 0.0_dp
       istat = PS_OK
       return
    end if

    if (max_knorm == 0.0_dp) then
       istat = PS_NO_PHASE_SPACE
       return
    end if

    lo = 0.0_dp
    hi = 1.0_dp
    fhi = scaled_energy(n, mass, kvec, hi) - target_energy
    do while (fhi < 0.0_dp)
       if (hi > 0.5_dp*huge(1.0_dp)) then
          istat = PS_NUMERIC_FAIL
          return
       end if
       hi = 2.0_dp * hi
       fhi = scaled_energy(n, mass, kvec, hi) - target_energy
    end do

    if (.not.ieee_is_finite(fhi)) then
       istat = PS_NUMERIC_FAIL
       return
    end if

    do iter = 1, 200
       mid = 0.5_dp*(lo + hi)
       fmid = scaled_energy(n, mass, kvec, mid) - target_energy
       if (.not.ieee_is_finite(fmid)) then
          istat = PS_NUMERIC_FAIL
          return
       end if
       if (abs(fmid) <= tol) exit
       if (fmid > 0.0_dp) then
          hi = mid
       else
          lo = mid
       end if
    end do

    scale = 0.5_dp*(lo + hi)
    istat = PS_OK
  end subroutine solve_spatial_scale

  subroutine unit_vector(vec, fallback, u)
    real(dp), intent(in) :: vec(3), fallback(3)
    real(dp), intent(out) :: u(3)
    real(dp) :: nrm

    nrm = norm3(vec)
    if (ieee_is_finite(nrm)) then
       if (nrm > 0.0_dp) then
          u = vec / nrm
          return
       end if
    end if

    nrm = norm3(fallback)
    if (ieee_is_finite(nrm)) then
       if (nrm > 0.0_dp) then
          u = fallback / nrm
          return
       end if
    end if
    u = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  end subroutine unit_vector

  subroutine perpendicular_unit(nhat, u)
    real(dp), intent(in) :: nhat(3)
    real(dp), intent(out) :: u(3)
    real(dp) :: trial(3), proj, nrm

    if (abs(nhat(1)) < 0.9_dp) then
       trial = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    else
       trial = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    end if

    proj = dot3(trial, nhat)
    u = trial - proj*nhat
    nrm = norm3(u)
    if (.not.ieee_is_finite(nrm)) then
       u = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    else if (nrm <= 32.0_dp*epsilon(1.0_dp)) then
       u = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    else
       u = u / nrm
    end if
  end subroutine perpendicular_unit

  function two_body_q(mparent, m1, m2, istat) result(q)
    real(dp), intent(in) :: mparent, m1, m2
    integer, intent(out) :: istat
    real(dp) :: q, scale, ap, a1, a2, f1, f2, f3, f4, tol

    q = 0.0_dp
    istat = PS_OK
    if (.not.finite_bounded_scalar(mparent) .or. .not.finite_bounded_scalar(m1) .or. &
         .not.finite_bounded_scalar(m2)) then
       istat = PS_BAD_INPUT
       return
    end if
    if (min(mparent, m1, m2) < 0.0_dp) then
       istat = PS_BAD_INPUT
       return
    end if

    scale = max(mparent, m1, m2)
    tol = 1.0e-12_dp*max(scale, tiny(1.0_dp))

    if (mparent < m1 + m2 - tol) then
       istat = PS_NO_PHASE_SPACE
       return
    end if

    if (mparent == 0.0_dp) then
       if (m1 == 0.0_dp .and. m2 == 0.0_dp) then
          q = 0.0_dp
          istat = PS_OK
       else
          istat = PS_NO_PHASE_SPACE
       end if
       return
    end if

    scale = max(mparent, m1, m2, tiny(1.0_dp))
    ap = mparent/scale
    a1 = m1/scale
    a2 = m2/scale
    f1 = ap - a1 - a2
    f2 = ap + a1 + a2
    f3 = ap - a1 + a2
    f4 = ap + a1 - a2
    if (min(f1, f2, f3, f4) < -1.0e-12_dp) then
       istat = PS_NO_PHASE_SPACE
       return
    end if
    q = scale*sqrt(max(0.0_dp, f1*f2*f3*f4))/(2.0_dp*ap)
    if (.not.finite_bounded_scalar(q)) then
       q = 0.0_dp
       istat = PS_NUMERIC_FAIL
    end if
  end function two_body_q

  subroutine split_pair_from_cluster(cvec, m1, m2, qrel, ehat, nhat, pi_out, pj_out, istat)
    real(dp), intent(in) :: cvec(0:3), m1, m2, qrel, ehat(3), nhat(3)
    real(dp), intent(out) :: pi_out(0:3), pj_out(0:3)
    integer, intent(out) :: istat

    real(dp) :: mcluster, kabs, cth, ei_star, ej_star, e_over_m, k_over_m
    real(dp) :: ei_lab, ppar_i, trans(3)

    istat = PS_OK
    pi_out = 0.0_dp
    pj_out = 0.0_dp

    if (.not.finite_bounded_vector(cvec) .or. .not.finite_bounded_scalar(m1) .or. &
         .not.finite_bounded_scalar(m2) .or. .not.finite_bounded_scalar(qrel) .or. &
         .not.all(ieee_is_finite(ehat)) .or. .not.all(ieee_is_finite(nhat))) then
       istat = PS_BAD_INPUT
       return
    end if
    if (min(m1, m2, qrel) < 0.0_dp .or. cvec(0) <= 0.0_dp) then
       istat = PS_BAD_INPUT
       return
    end if
    if (abs(norm3(ehat) - 1.0_dp) > 1.0e-6_dp .or. &
         abs(norm3(nhat) - 1.0_dp) > 1.0e-6_dp) then
       istat = PS_BAD_INPUT
       return
    end if

    kabs = norm3(cvec(1:3))
    cth = max(-1.0_dp, min(1.0_dp, dot3(ehat, nhat)))
    trans = ehat - cth*nhat

    if (m1 == 0.0_dp .and. m2 == 0.0_dp) then
       ! Stable massless formula. It remains well conditioned when cvec^2 -> 0.
       ei_lab = 0.5_dp*(cvec(0) + kabs*cth)
       ppar_i = 0.5_dp*(kabs + cvec(0)*cth)
       pi_out(0) = ei_lab
       pi_out(1:3) = ppar_i*nhat + qrel*trans
       pj_out = cvec - pi_out
       if (.not.finite_bounded_vector(pi_out) .or. .not.finite_bounded_vector(pj_out)) then
          pi_out = 0.0_dp
          pj_out = 0.0_dp
          istat = PS_NUMERIC_FAIL
       end if
       return
    end if

    mcluster = invariant_mass(cvec)
    if (.not.finite_bounded_scalar(mcluster)) then
       istat = PS_NO_PHASE_SPACE
       return
    end if
    if (mcluster <= 0.0_dp) then
       istat = PS_NO_PHASE_SPACE
       return
    end if

    ei_star = hypot2(m1, qrel)
    ej_star = hypot2(m2, qrel)
    if (abs((ei_star + ej_star) - mcluster) > &
         1.0e-8_dp*max(mcluster, ei_star + ej_star, tiny(1.0_dp))) then
       ! The construction assumes that cvec has the invariant mass implied by qrel.
       istat = PS_NUMERIC_FAIL
       return
    end if

    if (mcluster < max(cvec(0), kabs)/ps_quartic_input_limit) then
       istat = PS_NUMERIC_FAIL
       return
    end if
    e_over_m = cvec(0) / mcluster
    k_over_m = kabs / mcluster
    if (.not.ieee_is_finite(e_over_m) .or. .not.ieee_is_finite(k_over_m)) then
       istat = PS_NUMERIC_FAIL
       return
    end if
    if (max(e_over_m, k_over_m) > ps_quartic_input_limit) then
       istat = PS_NUMERIC_FAIL
       return
    end if

    ei_lab = e_over_m*ei_star + k_over_m*qrel*cth
    ppar_i = k_over_m*ei_star + e_over_m*qrel*cth

    pi_out(0) = ei_lab
    pi_out(1:3) = ppar_i*nhat + qrel*trans
    pj_out = cvec - pi_out
    if (.not.finite_bounded_vector(pi_out) .or. .not.finite_bounded_vector(pj_out)) then
       pi_out = 0.0_dp
       pj_out = 0.0_dp
       istat = PS_NUMERIC_FAIL
    end if
  end subroutine split_pair_from_cluster

  subroutine generate_phase_space(n, mass, ptot, p, status)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: ptot(0:3)
    real(dp), intent(out) :: p(0:3,n)
    integer, intent(out), optional :: status

    real(dp), allocatable :: q(:,:), qr(:,:), pcm(:,:), kvec(:,:)
    real(dp) :: ecm, sum_mass, beta(3), qsum(0:3), qmass, scale
    real(dp) :: r(4), cth, sth, phi, e, prod, tol
    integer :: i, st

    call set_status(status, PS_OK)
    p = 0.0_dp

    if (n < 1 .or. .not.masses_are_valid(n, mass) .or. &
         .not.finite_bounded_vector(ptot)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    if (ptot(0) <= 0.0_dp .or. normalized_msq4(ptot) <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    ecm = invariant_mass(ptot)
    if (.not.finite_bounded_scalar(ecm)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (ecm <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    sum_mass = 0.0_dp
    do i = 1, n
       if (mass(i) > ps_quartic_input_limit - sum_mass) then
          call set_status(status, PS_BAD_INPUT)
          return
       end if
       sum_mass = sum_mass + mass(i)
    end do
    tol = 1.0e-12_dp * max(ecm, sum_mass, tiny(1.0_dp))

    if (ecm < sum_mass - tol) then
       call set_status(status, PS_NO_PHASE_SPACE)
       return
    end if

    if (n == 1) then
       if (abs(ecm - mass(1)) > tol) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       p(:, 1) = ptot
       if (.not.final_state_is_valid(n, mass, p, ptot)) then
          p = 0.0_dp
          call set_status(status, PS_NUMERIC_FAIL)
          return
       end if
       call set_status(status, PS_OK)
       return
    end if

    allocate(q(0:3,n), qr(0:3,n), pcm(0:3,n), kvec(3,n))

    if (abs(ecm - sum_mass) <= tol) then
       pcm = 0.0_dp
       do i = 1, n
          pcm(0, i) = max(0.0_dp, mass(i))
       end do
       beta = ptot(1:3) / ptot(0)
       call boost_many(n, pcm, beta, p, st)
       if (st /= PS_OK) then
          p = 0.0_dp
       else if (.not.final_state_is_valid(n, mass, p, ptot)) then
          p = 0.0_dp
          st = PS_NUMERIC_FAIL
       end if
       call set_status(status, st)
       return
    end if

    do i = 1, n
       call random_number(r)
       cth = 2.0_dp*r(1) - 1.0_dp
       phi = 2.0_dp*pi*r(2)
       prod = max(r(3)*r(4), tiny(1.0_dp))
       e = -log(prod)
       sth = sqrt(max(0.0_dp, 1.0_dp - cth*cth))
       q(0, i) = e
       q(1, i) = e*sth*cos(phi)
       q(2, i) = e*sth*sin(phi)
       q(3, i) = e*cth
    end do

    qsum = 0.0_dp
    do i = 1, n
       qsum = qsum + q(:, i)
    end do

    qmass = invariant_mass(qsum)
    if (.not.ieee_is_finite(qmass)) then
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if
    if (qmass <= 0.0_dp) then
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if

    beta = -qsum(1:3) / qsum(0)
    call boost_many(n, q, beta, qr, st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    q = (ecm/qmass) * qr
    if (.not.all(ieee_is_finite(q))) then
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if
    if (maxval(abs(q)) > ps_quartic_input_limit) then
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if
    do i = 1, n
       kvec(:, i) = q(1:3, i)
    end do

    call solve_spatial_scale(n, mass, kvec, ecm, scale, st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    pcm = 0.0_dp
    do i = 1, n
       pcm(1:3, i) = scale*q(1:3, i)
       pcm(0, i) = hypot2(mass(i), norm3(pcm(1:3, i)))
    end do

    beta = ptot(1:3) / ptot(0)
    call boost_many(n, pcm, beta, p, st)
    if (st /= PS_OK) then
       p = 0.0_dp
    else if (.not.final_state_is_valid(n, mass, p, ptot)) then
       p = 0.0_dp
       st = PS_NUMERIC_FAIL
    end if
    call set_status(status, st)
  end subroutine generate_phase_space


  subroutine generate_scattering_phase_space(nout, mass, pinitial, p, status)
    integer, intent(in) :: nout
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: pinitial(0:3,2)
    real(dp), intent(out) :: p(0:3,nout+2)
    integer, intent(out), optional :: status

    real(dp) :: ptot(0:3), tol
    integer :: st

    call set_status(status, PS_OK)
    p = 0.0_dp

    if (nout < 1 .or. .not.masses_are_valid(nout, mass) .or. &
         .not.all(ieee_is_finite(pinitial))) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (maxval(abs(pinitial)) > ps_quartic_input_limit) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    if (pinitial(0,1) <= 0.0_dp .or. pinitial(0,2) <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    tol = validation_rel*max(max_abs4(pinitial(:,1)), max_abs4(pinitial(:,2)), &
         tiny(1.0_dp))**2
    if (abs(msq4(pinitial(:,1))) > tol .or. abs(msq4(pinitial(:,2))) > tol) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    ptot = pinitial(:, 1) + pinitial(:, 2)
    if (ptot(0) <= 0.0_dp .or. normalized_msq4(ptot) <= 0.0_dp .or. &
         .not.finite_bounded_vector(ptot)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    p(:, 1:2) = pinitial
    call generate_phase_space(nout, mass, ptot, p(:, 3:nout+2), st)
    if (st == PS_OK .and. .not.scattering_event_is_valid(nout, mass, p)) st = PS_NUMERIC_FAIL
    call set_status(status, st)
    if (st /= PS_OK) p = 0.0_dp
  end subroutine generate_scattering_phase_space

  subroutine generate_cm_scattering_phase_space(nout, mass, sqrts, p, status)
    integer, intent(in) :: nout
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: sqrts
    real(dp), intent(out) :: p(0:3,nout+2)
    integer, intent(out), optional :: status

    real(dp) :: pinitial(0:3,2), ebeam

    p = 0.0_dp
    if (nout < 1 .or. .not.masses_are_valid(nout, mass) .or. &
         .not.finite_bounded_scalar(sqrts)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (sqrts <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    ebeam = 0.5_dp*sqrts
    pinitial = 0.0_dp
    pinitial(0,1) = ebeam
    pinitial(3,1) = ebeam
    pinitial(0,2) = ebeam
    pinitial(3,2) = -ebeam

    call generate_scattering_phase_space(nout, mass, pinitial, p, status)
  end subroutine generate_cm_scattering_phase_space

  subroutine soft_deform(n, mass, pin, isoft, lambda, pout, status)
    integer, intent(in) :: n, isoft
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: pin(0:3,n)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,n)
    integer, intent(out), optional :: status

    integer :: i, a, c, st, nspec
    integer, allocatable :: idx(:)
    real(dp), allocatable :: pcm(:,:), pout_cm(:,:), krest(:,:), knew(:,:)
    real(dp), allocatable :: mm(:), kvec(:,:)
    real(dp) :: ptot(0:3), p_cm_tot(0:3), beta(3), beta_back(3), beta_r(3)
    real(dp) :: r0(0:3), rnew(0:3), mnew, scale, tol

    call set_status(status, PS_OK)
    pout = 0.0_dp

    if (n < 2 .or. isoft < 1 .or. isoft > n .or. &
         .not.finite_bounded_scalar(lambda) .or. &
         .not.masses_are_valid(n, mass) .or. .not.all(ieee_is_finite(pin))) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (lambda < 0.0_dp .or. maxval(abs(pin)) > ps_quartic_input_limit) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    call total_momentum(n, pin, ptot)
    if (.not.final_state_is_valid(n, mass, pin, ptot) .or. ptot(0) <= 0.0_dp .or. &
         normalized_msq4(ptot) <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    nspec = n - 1
    allocate(idx(nspec), pcm(0:3,n), pout_cm(0:3,n))
    allocate(krest(0:3,nspec), knew(0:3,nspec), mm(nspec), kvec(3,nspec))

    c = 0
    do i = 1, n
       if (i /= isoft) then
          c = c + 1
          idx(c) = i
          mm(c) = mass(i)
       end if
    end do

    beta = -ptot(1:3) / ptot(0)
    beta_back = ptot(1:3) / ptot(0)
    call boost_many(n, pin, beta, pcm, st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    p_cm_tot = 0.0_dp
    p_cm_tot(0) = invariant_mass(ptot)

    pout_cm = 0.0_dp
    pout_cm(1:3, isoft) = lambda * pcm(1:3, isoft)
    pout_cm(0, isoft) = hypot2(mass(isoft), norm3(pout_cm(1:3, isoft)))

    r0 = p_cm_tot - pcm(:, isoft)
    rnew = p_cm_tot - pout_cm(:, isoft)
    mnew = invariant_mass(rnew)

    tol = 1.0e-10_dp * max(mnew, mass(idx(1)), tiny(1.0_dp))

    if (nspec == 1) then
       if (abs(mnew - mass(idx(1))) > tol) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       pout_cm(:, idx(1)) = rnew
    else
       if (r0(0) <= 0.0_dp .or. normalized_msq4(r0) <= 0.0_dp) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       beta_r = -r0(1:3) / r0(0)
       do a = 1, nspec
          call boost_one(pcm(:, idx(a)), beta_r, krest(:, a), st)
          if (st /= PS_OK) then
             call set_status(status, st)
             return
          end if
          kvec(:, a) = krest(1:3, a)
       end do

       call solve_spatial_scale(nspec, mm, kvec, mnew, scale, st)
       if (st /= PS_OK) then
          call set_status(status, st)
          return
       end if

       do a = 1, nspec
          knew(1:3, a) = scale*krest(1:3, a)
          knew(0, a) = hypot2(mm(a), norm3(knew(1:3, a)))
       end do

       if (rnew(0) <= 0.0_dp) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       beta_r = rnew(1:3) / rnew(0)
       do a = 1, nspec
          call boost_one(knew(:, a), beta_r, pout_cm(:, idx(a)), st)
          if (st /= PS_OK) then
             call set_status(status, st)
             return
          end if
       end do
    end if

    call boost_many(n, pout_cm, beta_back, pout, st)
    if (st == PS_OK .and. .not.final_state_is_valid(n, mass, pout, ptot)) then
       pout = 0.0_dp
       st = PS_NUMERIC_FAIL
    end if
    call set_status(status, st)
  end subroutine soft_deform

  subroutine collinear_deform(n, mass, pin, i1, i2, lambda, pout, status)
    integer, intent(in) :: n, i1, i2
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: pin(0:3,n)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,n)
    integer, intent(out), optional :: status

    integer :: i, a, c, st, nspec
    integer, allocatable :: idx(:)
    real(dp), allocatable :: pcm(:,:), pout_cm(:,:), krest(:,:)
    real(dp), allocatable :: mm(:)
    real(dp) :: ptot(0:3), p_cm_tot(0:3), beta(3), beta_back(3), beta_r(3)
    real(dp) :: c0(0:3), r0(0:3), cnew(0:3), rnew(0:3), pstar_i(0:3)
    real(dp) :: m0, mr, q0, qnew, mnew, ecm, outer_p, outer_e_c, outer_e_r
    real(dp) :: lam_outer, lam_scale, lam_tol, nhat(3), ehat(3), fallback(3), tol

    call set_status(status, PS_OK)
    pout = 0.0_dp

    if (n < 2 .or. i1 < 1 .or. i1 > n .or. i2 < 1 .or. i2 > n .or. &
         i1 == i2 .or. .not.finite_bounded_scalar(lambda) .or. &
         .not.masses_are_valid(n, mass) .or. .not.all(ieee_is_finite(pin))) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (lambda < 0.0_dp .or. maxval(abs(pin)) > ps_quartic_input_limit) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    call total_momentum(n, pin, ptot)
    if (.not.final_state_is_valid(n, mass, pin, ptot) .or. ptot(0) <= 0.0_dp .or. &
         normalized_msq4(ptot) <= 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    nspec = n - 2
    if (nspec == 0 .and. abs(lambda - 1.0_dp) > 1.0e-14_dp) then
       call set_status(status, PS_NO_PHASE_SPACE)
       return
    end if

    allocate(idx(max(1, nspec)), pcm(0:3,n), pout_cm(0:3,n))
    allocate(krest(0:3,max(1, nspec)), mm(max(1, nspec)))

    c = 0
    do i = 1, n
       if (i /= i1 .and. i /= i2) then
          c = c + 1
          idx(c) = i
          mm(c) = mass(i)
       end if
    end do

    beta = -ptot(1:3) / ptot(0)
    beta_back = ptot(1:3) / ptot(0)
    call boost_many(n, pin, beta, pcm, st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    ecm = invariant_mass(ptot)
    p_cm_tot = 0.0_dp
    p_cm_tot(0) = ecm

    c0 = pcm(:, i1) + pcm(:, i2)
    r0 = p_cm_tot - c0
    m0 = invariant_mass(c0)
    mr = invariant_mass(r0)
    tol = 1.0e-10_dp * max(ecm, m0, mr, tiny(1.0_dp))

    q0 = two_body_q(m0, mass(i1), mass(i2), st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    fallback = pcm(1:3, i1)
    if (norm3(fallback) == 0.0_dp) fallback = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    call unit_vector(c0(1:3), fallback, nhat)

    if (m0 > 0.0_dp .and. q0 > 0.0_dp) then
       beta_r = -c0(1:3) / c0(0)
       call boost_one(pcm(:, i1), beta_r, pstar_i, st)
       if (st /= PS_OK) then
          call set_status(status, st)
          return
       end if
       call unit_vector(pstar_i(1:3), nhat, ehat)
    else
       call unit_vector(pcm(1:3, i1), nhat, ehat)
    end if

    if (nspec > 1) then
       if (r0(0) <= 0.0_dp .or. normalized_msq4(r0) <= 0.0_dp) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       beta_r = -r0(1:3) / r0(0)
       do a = 1, nspec
          call boost_one(pcm(:, idx(a)), beta_r, krest(:, a), st)
          if (st /= PS_OK) then
             call set_status(status, st)
             return
          end if
       end do
    end if

    if (q0 > 0.0_dp) then
       if (lambda > ps_quartic_input_limit/q0) then
          call set_status(status, PS_NUMERIC_FAIL)
          return
       end if
    end if
    qnew = lambda*q0
    mnew = hypot2(mass(i1), qnew) + hypot2(mass(i2), qnew)
    if (.not.finite_bounded_scalar(mnew)) then
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if

    if (ecm < mnew + mr - tol) then
       call set_status(status, PS_NO_PHASE_SPACE)
       return
    end if

    lam_outer = kallen(ecm*ecm, mnew*mnew, mr*mr)
    lam_scale = max(ecm, mnew, mr, tiny(1.0_dp))
    lam_tol = 1.0e-10_dp*lam_scale**4
    if (.not.ieee_is_finite(lam_outer)) then
       call set_status(status, PS_NO_PHASE_SPACE)
       return
    end if
    if (lam_outer < -lam_tol) then
       call set_status(status, PS_NO_PHASE_SPACE)
       return
    end if
    outer_p = sqrt(max(0.0_dp, lam_outer)) / (2.0_dp*ecm)
    outer_e_c = (ecm*ecm + mnew*mnew - mr*mr) / (2.0_dp*ecm)
    outer_e_r = ecm - outer_e_c

    cnew(0) = outer_e_c
    cnew(1:3) = outer_p*nhat
    rnew(0) = outer_e_r
    rnew(1:3) = -outer_p*nhat

    pout_cm = 0.0_dp

    if (nspec == 0) then
       ! Only possible for lambda = 1, already checked above.
    else if (nspec == 1) then
       pout_cm(:, idx(1)) = rnew
    else
       if (rnew(0) <= 0.0_dp) then
          call set_status(status, PS_NO_PHASE_SPACE)
          return
       end if
       beta_r = rnew(1:3) / rnew(0)
       do a = 1, nspec
          call boost_one(krest(:, a), beta_r, pout_cm(:, idx(a)), st)
          if (st /= PS_OK) then
             call set_status(status, st)
             return
          end if
       end do
    end if

    call split_pair_from_cluster(cnew, mass(i1), mass(i2), qnew, ehat, nhat, &
         pout_cm(:, i1), pout_cm(:, i2), st)
    if (st /= PS_OK) then
       call set_status(status, st)
       return
    end if

    call boost_many(n, pout_cm, beta_back, pout, st)
    if (st == PS_OK .and. .not.final_state_is_valid(n, mass, pout, ptot)) then
       pout = 0.0_dp
       st = PS_NUMERIC_FAIL
    end if
    call set_status(status, st)
  end subroutine collinear_deform

  subroutine soft_deform_event(nout, mass, pin, isoft_col, lambda, pout, status)
    integer, intent(in) :: nout, isoft_col
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: pin(0:3,nout+2)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,nout+2)
    integer, intent(out), optional :: status

    integer :: st, ifinal

    pout = 0.0_dp
    if (nout < 1 .or. isoft_col < 3 .or. isoft_col > nout + 2 .or. &
         .not.finite_bounded_scalar(lambda) .or. &
         .not.scattering_event_is_valid(nout, mass, pin)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (lambda < 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    ifinal = isoft_col - 2
    pout(:, 1:2) = pin(:, 1:2)
    call soft_deform(nout, mass, pin(:, 3:nout+2), ifinal, lambda, &
         pout(:, 3:nout+2), st)
    if (st == PS_NO_PHASE_SPACE) then
       ! At fixed incoming momentum a two-body final state with a massive
       ! spectator cannot generally reach a soft limit.  Vary s-hat instead.
       call soft_deform_variable_shat_event(nout, mass, pin, ifinal, lambda, pout, st)
    end if
    if (st /= PS_OK) then
       pout = 0.0_dp
       call set_status(status, st)
       return
    end if
    if (.not.scattering_event_is_valid(nout, mass, pout)) then
       pout = 0.0_dp
       call set_status(status, PS_NUMERIC_FAIL)
       return
    end if
    call set_status(status, PS_OK)
  end subroutine soft_deform_event

  subroutine soft_deform_variable_shat_event(nout, mass, pin, ifinal, lambda, pout, status)
    integer, intent(in) :: nout, ifinal
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: pin(0:3,nout+2)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,nout+2)
    integer, intent(out) :: status

    integer :: i, a, c, st, nspec
    integer, allocatable :: idx(:)
    real(dp), allocatable :: pcm(:,:), pout_cm(:,:), krest(:,:)
    real(dp) :: ptot(0:3), pnew(0:3), p_cm_tot(0:3), beta(3), beta_back(3), beta_r(3)
    real(dp) :: r0(0:3), rnew(0:3), m_spec, tol, e1, e2

    status = PS_OK
    pout = 0.0_dp
    if (nout < 2 .or. ifinal < 1 .or. ifinal > nout .or. &
         .not.finite_bounded_scalar(lambda) .or. &
         .not.scattering_event_is_valid(nout, mass, pin)) then
       status = PS_BAD_INPUT
       return
    end if
    if (lambda < 0.0_dp) then
       status = PS_BAD_INPUT
       return
    end if

    nspec = nout - 1
    allocate(idx(nspec), pcm(0:3,nout), pout_cm(0:3,nout), krest(0:3,nspec))
    c = 0
    do i = 1, nout
       if (i /= ifinal) then
          c = c + 1
          idx(c) = i
       end if
    end do

    ptot = pin(:, 1) + pin(:, 2)
    if (ptot(0) <= 0.0_dp .or. normalized_msq4(ptot) <= 0.0_dp) then
       status = PS_BAD_INPUT
       return
    end if
    beta = -ptot(1:3) / ptot(0)
    beta_back = -beta
    call boost_many(nout, pin(:, 3:nout+2), beta, pcm, st)
    if (st /= PS_OK) then
       status = st
       return
    end if

    p_cm_tot = 0.0_dp
    p_cm_tot(0) = invariant_mass(ptot)
    pout_cm = 0.0_dp
    pout_cm(1:3, ifinal) = lambda*pcm(1:3, ifinal)
    pout_cm(0, ifinal) = hypot2(mass(ifinal), norm3(pout_cm(1:3, ifinal)))

    r0 = p_cm_tot - pcm(:, ifinal)
    tol = 1.0e-10_dp*max(p_cm_tot(0), max_abs4(r0), tiny(1.0_dp))
    if (r0(0) <= 0.0_dp .or. normalized_msq4(r0) < -1.0e-10_dp) then
       status = PS_NO_PHASE_SPACE
       return
    end if
    m_spec = invariant_mass(r0)
    if (nspec > 1) then
       if (normalized_msq4(r0) <= 0.0_dp) then
          status = PS_NO_PHASE_SPACE
          return
       end if
       beta_r = -r0(1:3) / r0(0)
       do a = 1, nspec
          call boost_one(pcm(:, idx(a)), beta_r, krest(:, a), st)
          if (st /= PS_OK) then
             status = st
             return
          end if
       end do
    end if

    rnew = 0.0_dp
    rnew(1:3) = -pout_cm(1:3, ifinal)
    rnew(0) = hypot2(m_spec, norm3(rnew(1:3)))
    if (nspec == 1) then
       pout_cm(:,idx(1)) = rnew
    else
       beta_r = rnew(1:3) / rnew(0)
       do a = 1, nspec
          call boost_one(krest(:, a), beta_r, pout_cm(:, idx(a)), st)
          if (st /= PS_OK) then
             status = st
             return
          end if
       end do
    end if

    call boost_many(nout, pout_cm, beta_back, pout(:, 3:nout+2), st)
    if (st /= PS_OK) then
       status = st
       return
    end if
    call total_momentum(nout, pout(:, 3:nout+2), pnew)
    tol = 1.0e-10_dp*max(max_abs4(pnew), tiny(1.0_dp))
    if (abs(pnew(1)) > tol .or. abs(pnew(2)) > tol) then
       status = PS_BAD_INPUT
       return
    end if
    e1 = 0.5_dp*(pnew(0) + pnew(3))
    e2 = 0.5_dp*(pnew(0) - pnew(3))
    if (e1 <= 0.0_dp .or. e2 <= 0.0_dp) then
       status = PS_NO_PHASE_SPACE
       return
    end if
    pout(:, 1) = (/ e1, 0.0_dp, 0.0_dp, e1 /)
    pout(:, 2) = (/ e2, 0.0_dp, 0.0_dp, -e2 /)
    if (.not.scattering_event_is_valid(nout, mass, pout)) then
       pout = 0.0_dp
       status = PS_NUMERIC_FAIL
    end if
  end subroutine soft_deform_variable_shat_event

  subroutine collinear_deform_event(nout, mass, pin, icol1, icol2, lambda, pout, status)
    integer, intent(in) :: nout, icol1, icol2
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: pin(0:3,nout+2)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,nout+2)
    integer, intent(out), optional :: status

    integer :: st, if1, if2, iinit, ifinal
    logical :: c1_initial, c2_initial, c1_final, c2_final

    pout = 0.0_dp
    if (nout < 1 .or. icol1 < 1 .or. icol1 > nout + 2 .or. &
         icol2 < 1 .or. icol2 > nout + 2 .or. icol1 == icol2 .or. &
         .not.finite_bounded_scalar(lambda) .or. &
         .not.scattering_event_is_valid(nout, mass, pin)) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if
    if (lambda < 0.0_dp) then
       call set_status(status, PS_BAD_INPUT)
       return
    end if

    c1_initial = (icol1 == 1 .or. icol1 == 2)
    c2_initial = (icol2 == 1 .or. icol2 == 2)
    c1_final = (icol1 >= 3)
    c2_final = (icol2 >= 3)

    if (c1_final .and. c2_final) then
       if1 = icol1 - 2
       if2 = icol2 - 2
       pout(:, 1:2) = pin(:, 1:2)
       call collinear_deform(nout, mass, pin(:, 3:nout+2), if1, if2, lambda, &
            pout(:, 3:nout+2), st)
       if (st /= PS_OK) then
          pout = 0.0_dp
          call set_status(status, st)
          return
       end if
       pout(:, 1:2) = pin(:, 1:2)
       if (.not.scattering_event_is_valid(nout, mass, pout)) then
          pout = 0.0_dp
          call set_status(status, PS_NUMERIC_FAIL)
          return
       end if
       call set_status(status, PS_OK)
       return
    end if

    if (c1_initial .and. c2_final) then
       iinit = icol1
       ifinal = icol2
       call collinear_initial_final_deform(nout, mass, pin, iinit, ifinal, lambda, pout, st)
       if (st /= PS_OK) pout = 0.0_dp
       call set_status(status, st)
       return
    end if

    if (c2_initial .and. c1_final) then
       iinit = icol2
       ifinal = icol1
       call collinear_initial_final_deform(nout, mass, pin, iinit, ifinal, lambda, pout, st)
       if (st /= PS_OK) pout = 0.0_dp
       call set_status(status, st)
       return
    end if

    ! Two initial-state particles cannot be brought to a final-state collinear limit.
    call set_status(status, PS_BAD_INPUT)
  end subroutine collinear_deform_event

  subroutine collinear_initial_final_deform(nout, mass, pin, iinit_col, ifinal_col, lambda, pout, status)
    integer, intent(in) :: nout, iinit_col, ifinal_col
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: pin(0:3,nout+2)
    real(dp), intent(in) :: lambda
    real(dp), intent(out) :: pout(0:3,nout+2)
    integer, intent(out) :: status

    integer :: i, a, c, st, nspec, ifinal
    integer, allocatable :: idx(:)
    real(dp), allocatable :: pcm(:,:), pout_cm(:,:), krest(:,:), knew(:,:)
    real(dp), allocatable :: mm(:), kvec(:,:)
    real(dp) :: ptot(0:3), p_cm_tot(0:3)
    real(dp) :: beta(3), beta_back(3), beta_r(3)
    real(dp) :: f0(0:3), fnew(0:3), r0(0:3), rnew(0:3)
    real(dp) :: nhat(3), uperp(3), kvec_f(3), kperp(3)
    real(dp) :: pabs, kpar, kperp_abs, theta0, theta_new, mnew, scale

    status = PS_OK
    pout = 0.0_dp

    if (nout < 1 .or. iinit_col < 1 .or. iinit_col > 2 .or. &
         ifinal_col < 3 .or. ifinal_col > nout + 2 .or. &
         .not.finite_bounded_scalar(lambda) .or. &
         .not.scattering_event_is_valid(nout, mass, pin)) then
       status = PS_BAD_INPUT
       return
    end if
    if (lambda < 0.0_dp) then
       status = PS_BAD_INPUT
       return
    end if

    ifinal = ifinal_col - 2
    nspec = nout - 1

    if (nspec == 0) then
       if (abs(lambda - 1.0_dp) <= 1.0e-14_dp) then
          pout = pin
          status = PS_OK
       else
          status = PS_NO_PHASE_SPACE
       end if
       return
    end if

    ptot = pin(:, 1) + pin(:, 2)

    allocate(idx(max(1, nspec)), pcm(0:3,nout+2), pout_cm(0:3,nout+2))
    allocate(krest(0:3,max(1, nspec)), knew(0:3,max(1, nspec)))
    allocate(mm(max(1, nspec)), kvec(3,max(1, nspec)))

    c = 0
    do i = 1, nout
       if (i /= ifinal) then
          c = c + 1
          idx(c) = i + 2
          mm(c) = mass(i)
       end if
    end do

    beta = -ptot(1:3) / ptot(0)
    beta_back = ptot(1:3) / ptot(0)
    call boost_many(nout+2, pin, beta, pcm, st)
    if (st /= PS_OK) then
       status = st
       return
    end if

    p_cm_tot = 0.0_dp
    p_cm_tot(0) = invariant_mass(ptot)

    call unit_vector(pcm(1:3, iinit_col), (/ 0.0_dp, 0.0_dp, 1.0_dp /), nhat)

    f0 = pcm(:, ifinal_col)
    kvec_f = f0(1:3)
    pabs = norm3(kvec_f)

    if (pabs == 0.0_dp) then
       if (abs(lambda - 1.0_dp) <= 1.0e-14_dp) then
          pout = pin
          status = PS_OK
       else
          status = PS_NO_PHASE_SPACE
       end if
       return
    end if

    kpar = dot3(kvec_f, nhat)
    kperp = kvec_f - kpar*nhat
    kperp_abs = norm3(kperp)
    if (kperp_abs > 0.0_dp) then
       uperp = kperp / kperp_abs
    else
       call perpendicular_unit(nhat, uperp)
    end if

    theta0 = atan2(kperp_abs, kpar)
    theta_new = min(pi, lambda*theta0)

    fnew(0) = f0(0)
    fnew(1:3) = pabs*cos(theta_new)*nhat + pabs*sin(theta_new)*uperp

    r0 = p_cm_tot - f0
    rnew = p_cm_tot - fnew
    mnew = invariant_mass(rnew)

    pout_cm = 0.0_dp
    pout_cm(:, 1:2) = pcm(:, 1:2)
    pout_cm(:, ifinal_col) = fnew

    if (nspec == 1) then
       if (rnew(0) <= 0.0_dp .or. normalized_msq4(rnew) < -1.0e-10_dp) then
          status = PS_NO_PHASE_SPACE
          return
       end if
       if (abs(mnew - mm(1)) > &
            1.0e-8_dp*max(mnew, mm(1), tiny(1.0_dp))) then
          status = PS_NO_PHASE_SPACE
          return
       end if
       pout_cm(:, idx(1)) = rnew
    else
       if (r0(0) <= 0.0_dp .or. normalized_msq4(r0) <= 0.0_dp .or. &
            rnew(0) <= 0.0_dp .or. normalized_msq4(rnew) < -1.0e-10_dp) then
          status = PS_NO_PHASE_SPACE
          return
       end if

       beta_r = -r0(1:3) / r0(0)
       do a = 1, nspec
          call boost_one(pcm(:, idx(a)), beta_r, krest(:, a), st)
          if (st /= PS_OK) then
             status = st
             return
          end if
          kvec(:, a) = krest(1:3, a)
       end do

       call solve_spatial_scale(nspec, mm, kvec, mnew, scale, st)
       if (st /= PS_OK) then
          status = st
          return
       end if

       do a = 1, nspec
          knew(1:3, a) = scale*krest(1:3, a)
          knew(0, a) = hypot2(mm(a), norm3(knew(1:3, a)))
       end do

       beta_r = rnew(1:3) / rnew(0)
       do a = 1, nspec
          call boost_one(knew(:, a), beta_r, pout_cm(:, idx(a)), st)
          if (st /= PS_OK) then
             status = st
             return
          end if
       end do
    end if

    call boost_many(nout+2, pout_cm, beta_back, pout, st)
    if (st /= PS_OK) then
       status = st
       return
    end if

    ! Keep the beam momenta bitwise equal to the input convention.
    pout(:, 1:2) = pin(:, 1:2)
    if (.not.scattering_event_is_valid(nout, mass, pout)) then
       pout = 0.0_dp
       status = PS_NUMERIC_FAIL
       return
    end if
    status = PS_OK
  end subroutine collinear_initial_final_deform

  subroutine event_residuals(n, mass, p, total, momentum_residual, mass_residual)
    integer, intent(in) :: n
    real(dp), intent(in) :: mass(n)
    real(dp), intent(in) :: p(0:3,n), total(0:3)
    real(dp), intent(out) :: momentum_residual, mass_residual
    integer :: i
    real(dp) :: psum(0:3), dm

    call total_momentum(n, p, psum)
    momentum_residual = maxval(abs(psum - total))
    mass_residual = 0.0_dp
    do i = 1, n
       dm = abs(msq4(p(:, i)) - mass(i)*mass(i))
       if (dm > mass_residual) mass_residual = dm
    end do
  end subroutine event_residuals

  subroutine scattering_event_residuals(nout, mass, p, momentum_residual, &
       final_mass_residual, incoming_mass_residual)
    integer, intent(in) :: nout
    real(dp), intent(in) :: mass(nout)
    real(dp), intent(in) :: p(0:3,nout+2)
    real(dp), intent(out) :: momentum_residual, final_mass_residual
    real(dp), intent(out), optional :: incoming_mass_residual

    integer :: i
    real(dp) :: psum_in(0:3), psum_out(0:3), dm, im

    psum_in = p(:, 1) + p(:, 2)
    call total_momentum(nout, p(:, 3:nout+2), psum_out)
    momentum_residual = maxval(abs(psum_out - psum_in))

    final_mass_residual = 0.0_dp
    do i = 1, nout
       dm = abs(msq4(p(:, i+2)) - mass(i)*mass(i))
       if (dm > final_mass_residual) final_mass_residual = dm
    end do

    if (present(incoming_mass_residual)) then
       im = max(abs(msq4(p(:, 1))), abs(msq4(p(:, 2))))
       incoming_mass_residual = im
    end if
  end subroutine scattering_event_residuals

end module phase_space_module
