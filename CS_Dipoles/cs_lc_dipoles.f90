module cs_lc_spin_dipoles
  use cs_dipole_mappings, only: dp, dot4, new_index,cs_dot4_scale,&
       cs_value_is_resolved,cs_roundoff_tolerance,cs_normalize_unit_interval,&
       cs_quartic_input_limit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(dp), parameter :: pi_dp = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: ca = 3.0_dp
  real(dp), parameter :: cf_lc = ca / 2.0_dp
  real(dp), parameter :: tr = 0.5_dp
  real(dp), parameter :: tr_u1 = tr / ca
  real(dp), parameter :: tiny_dip = tiny(1.0_dp)

  integer, parameter :: ch_none = 0
  integer, parameter :: ch_i_qg = 11
  integer, parameter :: ch_i_gq = 12
  integer, parameter :: ch_i_qq = 13
  integer, parameter :: ch_i_gg = 14

contains

  logical function is_initial(idx)
    integer, intent(in) :: idx
    is_initial = idx <= 2
  end function is_initial


  pure logical function dipole_status_is_numerical(status) result(is_numerical)
    ! Status values caused by a phase-space boundary or a denominator that
    ! cannot be resolved reliably.  Invalid flavours, arguments, and
    ! unsupported kernels are deliberately excluded: those are programming
    ! errors and must remain fatal to the caller.
    integer, intent(in) :: status
    select case (status)
    case (-20,-17:-10,-101,-103,-104,-106,-108,-201,-202,-203,-301,-302,-303)
       is_numerical=.true.
    case default
       is_numerical=.false.
    end select
  end function dipole_status_is_numerical


  logical function is_gluon(f)
    integer, intent(in) :: f
    is_gluon = (f == 21 .or. f == 99)
  end function is_gluon


  logical function is_u1_gluon(f)
    integer, intent(in) :: f
    is_u1_gluon = (f == 99)
  end function is_u1_gluon


  logical function is_quark(f)
    integer, intent(in) :: f
    is_quark = (f >= 1 .and. f <= 6) .or. (f <= -1 .and. f >= -6)
  end function is_quark


  logical function same_quark_flavour(f1,f2)
    integer, intent(in) :: f1,f2
    same_quark_flavour=.false.
    if (.not.is_quark(f1) .or. .not.is_quark(f2)) return
    same_quark_flavour=f1 == f2 .or. f1 == -f2
  end function same_quark_flavour


  logical function is_q_qbar_pair(f1, f2)
    integer, intent(in) :: f1, f2
    is_q_qbar_pair = .false.
    if (.not.is_quark(f1) .or. .not.is_quark(f2)) return
    is_q_qbar_pair = f1 == -f2
  end function is_q_qbar_pair


  real(dp) function vector_trace(f)
    integer, intent(in) :: f
    if (is_u1_gluon(f)) then
       vector_trace = tr_u1
    else
       vector_trace = tr
    end if
  end function vector_trace


  real(dp) function gcontra(mu, nu) result(g)
    integer, intent(in) :: mu, nu

    if (mu /= nu) then
       g = 0.0_dp
    else if (mu == 0) then
       g = 1.0_dp
    else
       g = -1.0_dp
    end if
  end function gcontra


  real(dp) function minus_gcontra(mu, nu) result(g)
    integer, intent(in) :: mu, nu
    g = -gcontra(mu, nu)
  end function minus_gcontra


  subroutine lower_complex_vector(v, vlower)
    complex(dp), intent(in) :: v(0:3)
    complex(dp), intent(out) :: vlower(0:3)

    vlower(0) = v(0)
    vlower(1) = -v(1)
    vlower(2) = -v(2)
    vlower(3) = -v(3)
  end subroutine lower_complex_vector


  subroutine tensor_to_helicity(vten, eps, vhel)
    ! vhel(a,b) = eps_a^*_mu V^{mu nu} eps_b_nu.
    ! eps is assumed to contain contravariant polarization vectors eps^mu.

    real(dp), intent(in) :: vten(0:3,0:3)
    complex(dp), intent(in) :: eps(0:3,2)
    complex(dp), intent(out) :: vhel(2,2)

    complex(dp) :: elow(0:3,2)
    integer :: a, b, mu, nu

    do a = 1, 2
       call lower_complex_vector(eps(:,a), elow(:,a))
    end do

    vhel = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do a = 1, 2
       do b = 1, 2
          do mu = 0, 3
             do nu = 0, 3
                vhel(a,b) = vhel(a,b) + conjg(elow(mu,a)) * vten(mu,nu) * elow(nu,b)
             end do
          end do
       end do
    end do
  end subroutine tensor_to_helicity


  real(dp) function contract_rho_v(rho, vhel) result(val)
    ! rho(lambda,lambda') = sum_other_helicities M(lambda) conjg(M(lambda')).
    ! Then M^* V M = sum rho(lambda,lambda') V(lambda',lambda).

    complex(dp), intent(in) :: rho(2,2), vhel(2,2)
    complex(dp) :: cval
    integer :: a, b

    cval = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do a = 1, 2
       do b = 1, 2
          cval = cval + rho(a,b) * vhel(b,a)
       end do
    end do

    val = real(cval, dp)
  end function contract_rho_v


  subroutine zero_tensor(vten)
    real(dp), intent(out) :: vten(0:3,0:3)
    vten = 0.0_dp
  end subroutine zero_tensor


  subroutine add_minus_g_term(vten, coeff)
    real(dp), intent(inout) :: vten(0:3,0:3)
    real(dp), intent(in) :: coeff
    integer :: mu, nu

    do mu = 0, 3
       do nu = 0, 3
          vten(mu,nu) = vten(mu,nu) + coeff * minus_gcontra(mu,nu)
       end do
    end do
  end subroutine add_minus_g_term


  subroutine add_outer_term(vten, coeff, r)
    real(dp), intent(inout) :: vten(0:3,0:3)
    real(dp), intent(in) :: coeff, r(0:3)
    integer :: mu, nu

    do mu = 0, 3
       do nu = 0, 3
          vten(mu,nu) = vten(mu,nu) + coeff * r(mu) * r(nu)
       end do
    end do
  end subroutine add_outer_term


  subroutine scalar_to_vhel(scalar, vhel)
    real(dp), intent(in) :: scalar
    complex(dp), intent(out) :: vhel(2,2)

    vhel = cmplx(0.0_dp, 0.0_dp, kind=dp)
    vhel(1,1) = cmplx(scalar, 0.0_dp, kind=dp)
    vhel(2,2) = cmplx(scalar, 0.0_dp, kind=dp)
  end subroutine scalar_to_vhel


  subroutine final_velocity_factors(p_i,p_j,p_k,mass_parent,mass_spectator,vreal,vborn,info)
    ! Massive final-final relative velocities.  Evaluate the Kallen
    ! functions with scale-aware boundary handling and validate the two
    ! scalar-product denominators before forming either ratio.
    real(dp), intent(in) :: p_i(0:3),p_j(0:3),p_k(0:3)
    real(dp), intent(in) :: mass_parent,mass_spectator
    real(dp), intent(out) :: vreal,vborn
    integer, intent(out) :: info
    real(dp) :: qvec(0:3),pairvec(0:3)
    real(dp) :: q2,pij2,mk2,mp2,q2_scale,pij2_scale
    real(dp) :: lambda_real,lambda_born,lambda_real_scale,lambda_born_scale
    real(dp) :: denom_real,denom_born,denom_real_scale,denom_born_scale
    logical :: valid

    vreal=0.0_dp
    vborn=0.0_dp
    info=0
    if (.not.all(ieee_is_finite(p_i)) .or. .not.all(ieee_is_finite(p_j)) .or. &
         .not.all(ieee_is_finite(p_k)) .or. .not.ieee_is_finite(mass_parent) .or. &
         .not.ieee_is_finite(mass_spectator)) then
       info=-20
       return
    endif
    if (mass_parent.lt.0.0_dp .or. mass_spectator.lt.0.0_dp) then
       info=-20
       return
    endif

    qvec=p_i+p_j+p_k
    pairvec=p_i+p_j
    q2=dot4(qvec,qvec)
    pij2=dot4(pairvec,pairvec)
    mk2=mass_spectator*mass_spectator
    mp2=mass_parent*mass_parent
    q2_scale=cs_dot4_scale(qvec,qvec)
    pij2_scale=cs_dot4_scale(pairvec,pairvec)
    if (.not.ieee_is_finite(q2) .or. .not.ieee_is_finite(pij2) .or. &
         .not.ieee_is_finite(q2_scale) .or. .not.ieee_is_finite(pij2_scale)) then
       info=-20
       return
    endif
    if (q2.le.cs_roundoff_tolerance(q2_scale)) then
       info=-20
       return
    endif
    if (pij2.lt.-cs_roundoff_tolerance(pij2_scale)) then
       info=-20
       return
    endif

    lambda_real=q2*q2+pij2*pij2+mk2*mk2-&
         2.0_dp*(q2*pij2+q2*mk2+pij2*mk2)
    lambda_born=q2*q2+mp2*mp2+mk2*mk2-&
         2.0_dp*(q2*mp2+q2*mk2+mp2*mk2)
    lambda_real_scale=abs(q2*q2)+abs(pij2*pij2)+abs(mk2*mk2)+&
         2.0_dp*(abs(q2*pij2)+abs(q2*mk2)+abs(pij2*mk2))
    lambda_born_scale=abs(q2*q2)+abs(mp2*mp2)+abs(mk2*mk2)+&
         2.0_dp*(abs(q2*mp2)+abs(q2*mk2)+abs(mp2*mk2))
    if (.not.ieee_is_finite(lambda_real) .or. .not.ieee_is_finite(lambda_born) .or. &
         .not.ieee_is_finite(lambda_real_scale) .or. &
         .not.ieee_is_finite(lambda_born_scale)) then
       info=-20
       return
    endif
    if (lambda_real.le.cs_roundoff_tolerance(lambda_real_scale) .or. &
         lambda_born.lt.-cs_roundoff_tolerance(lambda_born_scale)) then
       info=-20
       return
    endif

    denom_real=q2-pij2-mk2
    denom_born=q2-mp2-mk2
    denom_real_scale=abs(q2)+abs(pij2)+mk2
    denom_born_scale=abs(q2)+mp2+mk2
    if (denom_real.le.cs_roundoff_tolerance(denom_real_scale) .or. &
         denom_born.le.cs_roundoff_tolerance(denom_born_scale)) then
       info=-20
       return
    endif
    vreal=sqrt(lambda_real)/denom_real
    vborn=sqrt(max(0.0_dp,lambda_born))/denom_born
    if (.not.ieee_is_finite(vreal) .or. .not.ieee_is_finite(vborn)) then
       vreal=0.0_dp
       vborn=0.0_dp
       info=-20
       return
    endif
    call cs_normalize_unit_interval(vreal,valid)
    if (.not.valid) then
       vreal=0.0_dp
       vborn=0.0_dp
       info=-20
       return
    endif
    call cs_normalize_unit_interval(vborn,valid)
    if (.not.valid) then
       vreal=0.0_dp
       vborn=0.0_dp
       info=-20
       return
    endif
    if (vreal.le.cs_roundoff_tolerance(1.0_dp)) then
       vreal=0.0_dp
       vborn=0.0_dp
       info=-20
    endif
  end subroutine final_velocity_factors


  subroutine cs_lc_dipole_spinrho(p, flav_real, flav_born, ijk, alpha_s, rho, eps_parent, dip, info, lc_weight, &
       mass_real, mass_parent)
    ! Leading-colour, spin-correlated Catani-Seymour dipole.  The optional
    ! masses activate the massive final-state kernels.
    !
    ! p(:,1), p(:,2) are incoming physical momenta with positive energy.
    ! j is the unresolved final-state parton and must satisfy j > 2.
    !
    ! rho(lambda,lambda') is the Born spin-density matrix of the mapped
    ! emitter/parent leg:
    !
    !   rho(lambda,lambda') = sum_h M(lambda,h) conjg(M(lambda',h)).
    !
    ! eps_parent(:,lambda) are contravariant polarization vectors for the
    ! mapped parent leg. They are used only if the mapped parent is a gluon.
    !
    ! lc_weight replaces -T_spectator . T_emitter / T_emitter^2.
    ! Use zero for non-leading-colour-connected spectator choices.
    !
    ! The result includes the local 8*pi*alpha_s or 16*pi*alpha_s factors
    ! in the Catani-Seymour splitting kernels, but no flux, PDF, symmetry,
    ! phase-space, or measurement-function factors.

    real(dp), intent(in) :: p(0:,:)
    integer, intent(in) :: flav_real(:), flav_born(:)
    integer, dimension(3),intent(in) :: ijk
    real(dp), intent(in) :: alpha_s
    complex(dp), intent(in) :: rho(2,2)
    complex(dp), intent(in) :: eps_parent(0:3,2)
    real(dp), intent(out) :: dip
    integer, intent(out) :: info
    real(dp), intent(in), optional :: lc_weight
    real(dp), intent(in), optional :: mass_real(:), mass_parent

    integer :: n, ni, nk, istat,i,j,k
    logical :: iini, kini
    real(dp) :: wt, pref, vcontract_alt
    real(dp) :: sij, sik, sjk, dotij, sij_parent
    real(dp) :: dotij_scale,parent_scale,sik_scale,sjk_scale
    real(dp) :: topology_den,topology_scale
    real(dp) :: mi, mj, mk, mparent
    real(dp) :: x, y, z, u
    complex(dp) :: vhel(2,2)
    real(dp) :: vcontract
    logical :: use_masses,valid

    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
    
    dip = 0.0_dp
    istat = 0
    info = 0

    if (present(lc_weight)) then
       wt = lc_weight
    else
       wt = 1.0_dp
    end if

    n = size(p, 2)

    use_masses = present(mass_real) .or. present(mass_parent)
    if (use_masses) then
       if (.not.(present(mass_real) .and. present(mass_parent))) then
          istat = -2
          goto 900
       endif
       if (size(mass_real) /= n) then
          istat = -2
          goto 900
       endif
       if (.not.all(ieee_is_finite(mass_real)) .or. .not.ieee_is_finite(mass_parent)) then
          istat=-20
          goto 900
       endif
       if (any(mass_real < 0.0_dp) .or. mass_parent < 0.0_dp) then
          istat = -2
          goto 900
       endif
    else
       mi = 0.0_dp
       mj = 0.0_dp
       mk = 0.0_dp
       mparent = 0.0_dp
    endif

    if (size(p, 1) /= 4) then
       istat = -1
       goto 900
    end if

    if (size(flav_real) /= n .or. size(flav_born) /= n - 1) then
       istat = -2
       goto 900
    end if
    if (.not.all(ieee_is_finite(p)) .or. .not.ieee_is_finite(alpha_s) .or. &
         .not.ieee_is_finite(wt)) then
       istat=-20
       goto 900
    endif
    if (.not.all(ieee_is_finite(real(rho,dp))) .or. &
         .not.all(ieee_is_finite(aimag(rho))) .or. &
         .not.all(ieee_is_finite(real(eps_parent,dp))) .or. &
         .not.all(ieee_is_finite(aimag(eps_parent)))) then
       istat=-20
       goto 900
    endif
    if (maxval(abs(p)).gt.cs_quartic_input_limit) then
       istat=-20
       goto 900
    endif
    if (use_masses) then
       if (max(maxval(abs(mass_real)),abs(mass_parent)).gt.cs_quartic_input_limit) then
          istat=-20
          goto 900
       endif
    endif
    if (alpha_s.lt.0.0_dp) then
       istat=-2
       goto 900
    endif

    if (i < 1 .or. i > n .or. j < 1 .or. j > n .or. k < 1 .or. k > n) then
       istat = -3
       goto 900
    end if

    if (i == j .or. i == k .or. j == k) then
       istat = -4
       goto 900
    end if

    if (j <= 2) then
       istat = -5
       goto 900
    end if

    if (use_masses) then
       mi = mass_real(i)
       mj = mass_real(j)
       mk = mass_real(k)
       mparent = mass_parent
       if (mj > 100.0_dp*epsilon(1.0_dp)*max(tiny(1.0_dp),abs(mparent))) then
          istat = -7
          goto 900
       endif
    endif

    ni = new_index(i, j)
    nk = new_index(k, j)

    if (ni <= 0 .or. nk <= 0) then
       istat = -6
       goto 900
    end if

    if (abs(wt) <= tiny_dip) then
       dip = 0.0_dp
       goto 900
    end if

    iini = is_initial(i)
    kini = is_initial(k)

    dotij = dot4(p(:,i), p(:,j))
    dotij_scale=cs_dot4_scale(p(:,i),p(:,j))
    sij = 2.0_dp*dotij
    ! Evaluate the parent virtuality from masses and p_i.p_j.  Forming
    ! (p_i+p_j)^2 directly loses the same small invariant through a second
    ! large cancellation in a collinear configuration.
    sij_parent=mi*mi+mj*mj+2.0_dp*dotij-mparent*mparent
    parent_scale=mi*mi+mj*mj+mparent*mparent+2.0_dp*dotij_scale

    if (.not.cs_value_is_resolved(dotij,dotij_scale) .or. &
         .not.cs_value_is_resolved(sij_parent,parent_scale)) then
       istat = -10
       goto 900
    end if

    if (.not. iini .and. .not. kini) then
       ! FF: final emitter, final unresolved, final spectator.

       sik = 2.0_dp * dot4(p(:,i), p(:,k))
       sjk = 2.0_dp * dot4(p(:,j), p(:,k))
       sik_scale=2.0_dp*cs_dot4_scale(p(:,i),p(:,k))
       sjk_scale=2.0_dp*cs_dot4_scale(p(:,j),p(:,k))

       topology_den=sik+sjk
       topology_scale=sik_scale+sjk_scale
       if (.not.cs_value_is_resolved(topology_den,topology_scale)) then
          istat = -11
          goto 900
       endif

       z=sik/topology_den
       call cs_normalize_unit_interval(z,valid)
       if (.not.valid) then
          istat=-11
          goto 900
       endif
       topology_den=sij+sik+sjk
       topology_scale=2.0_dp*dotij_scale+sik_scale+sjk_scale
       if (.not.cs_value_is_resolved(topology_den,topology_scale)) then
          istat=-11
          goto 900
       endif
       y=sij/topology_den
       call cs_normalize_unit_interval(y,valid)
       if (.not.valid) then
          istat=-11
          goto 900
       endif

       call final_splitting_matrix(.false., alpha_s, flav_real(i), flav_real(j), flav_born(ni), &
            z, y_dummy=y, x_dummy=0.0_dp, p_i=p(:,i), p_j=p(:,j), &
            p_k=p(:,k), mass_i=mi, mass_j=mj, mass_k=mk, mass_parent=mparent, &
            eps_parent=eps_parent, vhel=vhel, info=istat)
       
       if (istat /= 0) goto 900

       pref = wt / sij_parent
       vcontract = contract_rho_v(rho, vhel)
       dip = pref * vcontract

    else if (.not. iini .and. kini) then

       ! FI: final emitter, final unresolved, initial spectator.

       sik = 2.0_dp * dot4(p(:,i), p(:,k))
       sjk = 2.0_dp * dot4(p(:,j), p(:,k))
       sik_scale=2.0_dp*cs_dot4_scale(p(:,i),p(:,k))
       sjk_scale=2.0_dp*cs_dot4_scale(p(:,j),p(:,k))

       topology_den=sik+sjk
       topology_scale=sik_scale+sjk_scale
       if (.not.cs_value_is_resolved(topology_den,topology_scale)) then
          istat = -12
          goto 900
       end if

       z=sik/topology_den
       x=(topology_den-sij)/topology_den
       call cs_normalize_unit_interval(z,valid)
       if (.not.valid) then
          istat=-12
          goto 900
       endif
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-13
          goto 900
       endif

       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat = -13
          goto 900
       end if

       call final_splitting_matrix(.true., alpha_s, flav_real(i), flav_real(j), flav_born(ni), &
            z, y_dummy=0.0_dp, x_dummy=x, p_i=p(:,i), p_j=p(:,j), &
            p_k=p(:,k), mass_i=mi, mass_j=mj, mass_k=mk, mass_parent=mparent, &
            eps_parent=eps_parent, vhel=vhel, info=istat)
       if (istat /= 0) goto 900

       pref = wt / (sij_parent * x)
       vcontract = contract_rho_v(rho, vhel)
       dip = pref * vcontract

    else if (iini .and. .not. kini) then

       ! IF: initial emitter, final unresolved, final spectator.

       sik = 2.0_dp * dot4(p(:,i), p(:,k))
       sjk = 2.0_dp * dot4(p(:,j), p(:,k))
       sik_scale=2.0_dp*cs_dot4_scale(p(:,i),p(:,k))

       topology_den=sij+sik
       topology_scale=2.0_dp*dotij_scale+sik_scale
       if (.not.cs_value_is_resolved(topology_den,topology_scale)) then
          istat = -14
          goto 900
       end if

       x=(topology_den-sjk)/topology_den
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-15
          goto 900
       endif

       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat = -15
          goto 900
       end if

       u=sij/topology_den
       call cs_normalize_unit_interval(u,valid)
       if (.not.valid) then
          istat=-14
          goto 900
       endif

       call initial_final_splitting_matrix(alpha_s, flav_real(i), flav_real(j), flav_born(ni), &
            x, u, p_emit=p(:,i), p_unres=p(:,j), p_spec=p(:,k), mass_spec=mk, eps_parent=eps_parent, &
            vhel=vhel, info=istat)
       if (istat /= 0) goto 900

       pref = wt / (sij_parent * x)
       vcontract = contract_rho_v(rho, vhel)
       dip = pref * vcontract

    else

       ! II: initial emitter, final unresolved, initial spectator.

       sik = 2.0_dp * dot4(p(:,i), p(:,k))
       sjk = 2.0_dp * dot4(p(:,j), p(:,k))
       sik_scale=2.0_dp*cs_dot4_scale(p(:,i),p(:,k))

       if (.not.cs_value_is_resolved(sik,sik_scale)) then
          istat = -16
          goto 900
       end if

       x = (sik - sij - sjk) / sik
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-17
          goto 900
       endif

       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat = -17
          goto 900
       end if

       call initial_initial_splitting_matrix(alpha_s, flav_real(i), flav_real(j), flav_born(ni), &
            x, p_emit=p(:,i), p_unres=p(:,j), p_spec=p(:,k), eps_parent=eps_parent, &
            vhel=vhel, info=istat)
       if (istat /= 0) goto 900

       pref = wt / (sij_parent * x)
       vcontract = contract_rho_v(rho, vhel)
       dip = pref * vcontract

    end if

900 continue
    if (istat.eq.0 .and. .not.ieee_is_finite(dip)) then
       dip=0.0_dp
       istat=-20
    endif
    info = istat
  end subroutine cs_lc_dipole_spinrho


  subroutine final_splitting_matrix(is_fi, alpha_s, fi, fj, fp, z, y_dummy, x_dummy, p_i, p_j, p_k, &
       mass_i, mass_j, mass_k, mass_parent, eps_parent, vhel, info)
    ! Final-state emitter kernels for FF and FI.
    ! If is_fi = .false., use y_dummy as y_ij,k.
    ! If is_fi = .true.,  use x_dummy as x_ij,a.

    logical, intent(in) :: is_fi
    real(dp), intent(in) :: alpha_s
    integer, intent(in) :: fi, fj, fp
    real(dp), intent(in) :: z, y_dummy, x_dummy
    real(dp), intent(in) :: p_i(0:3), p_j(0:3), p_k(0:3)
    real(dp), intent(in) :: mass_i, mass_j, mass_k, mass_parent
    complex(dp), intent(in) :: eps_parent(0:3,2)
    complex(dp), intent(out) :: vhel(2,2)
    integer, intent(out) :: info

    real(dp) :: zi, zj, zq, denom_i, denom_j
    real(dp) :: scalar, aterm, coeff, dotij, y, vreal, vtilde
    real(dp) :: dotij_scale
    real(dp) :: zim, zjm
    real(dp) :: r(0:3), vten(0:3,0:3)
    logical :: massive
    integer :: velocity_info

    info = 0
    vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
    zi = z
    zj = 1.0_dp - z
    dotij = dot4(p_i, p_j)
    dotij_scale=cs_dot4_scale(p_i,p_j)
    massive = mass_i > 0.0_dp .or. mass_j > 0.0_dp .or. &
         mass_k > 0.0_dp .or. mass_parent > 0.0_dp

    if (.not.cs_value_is_resolved(dotij,dotij_scale)) then
       info = -101
       return
    end if

    if (is_quark(fp)) then

       ! Parent quark: q -> q g.

       if (fi == fp .and. is_gluon(fj)) then
          zq = zi
       else if (is_gluon(fi) .and. fj == fp) then
          zq = zj
       else
          info = -102
          return
       end if

       if (is_fi) then
          denom_i = 1.0_dp - zq + (1.0_dp - x_dummy)
       else
          denom_i = 1.0_dp - zq * (1.0_dp - y_dummy)
       end if

       if (.not.cs_value_is_resolved(denom_i,1.0_dp+abs(zq)+abs(x_dummy)+abs(y_dummy))) then
          info = -103
          return
       end if

       if (massive .and. .not.is_fi) then
          call final_velocity_factors(p_i,p_j,p_k,mass_parent,mass_k,vreal,vtilde,velocity_info)
          if (velocity_info.ne.0) then
             info = -106
             return
          endif
          scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * &
               (2.0_dp/denom_i - (vtilde/vreal)*(1.0_dp + zq + mass_parent*mass_parent/dotij))
       else if (massive .and. is_fi) then
          scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * &
               (2.0_dp/denom_i - (1.0_dp + zq) - mass_parent*mass_parent/dotij)
       else
          scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * (2.0_dp/denom_i - (1.0_dp + zq))
       endif
       call scalar_to_vhel(scalar, vhel)

    else

       ! Parent gluon: g -> g g or g -> q qbar.

       if (is_gluon(fi) .and. is_gluon(fj)) then

          if (is_u1_gluon(fp)) then
             ! The U(1) gluon has no triple-vector coupling.
             info = -105
             return
          end if

          if (is_fi) then
             denom_i = 1.0_dp - zi + (1.0_dp - x_dummy)
             denom_j = 1.0_dp - zj + (1.0_dp - x_dummy)
          else
             denom_i = 1.0_dp - zi * (1.0_dp - y_dummy)
             denom_j = 1.0_dp - zj * (1.0_dp - y_dummy)
          end if

          if (.not.cs_value_is_resolved(denom_i,1.0_dp+abs(zi)+abs(x_dummy)+abs(y_dummy)) .or. &
               .not.cs_value_is_resolved(denom_j,1.0_dp+abs(zj)+abs(x_dummy)+abs(y_dummy))) then
             info = -104
             return
          end if

          if (massive .and. (mass_i > 0.0_dp .or. mass_j > 0.0_dp .or. &
               mass_parent > 0.0_dp)) then
             info = -107
             return
          endif
          ! Ordered LC dipoles designate leg j as unresolved.  The
          ! complementary i-soft pole is supplied by the dipole where i is
          ! the unresolved leg; keeping it here would over-subtract soft
          ! gluons in a colour-ordered sum.
          aterm = 1.0_dp/denom_i - 1.0_dp
          if (massive .and. .not.is_fi) then
             call final_velocity_factors(p_i,p_j,p_k,mass_parent,mass_k,vreal,vtilde,velocity_info)
             if (velocity_info.ne.0) then
                info = -107
                return
             endif
          else
             vreal = 1.0_dp
          endif
          zim = zi - 0.5_dp*(1.0_dp-vreal)
          zjm = zj - 0.5_dp*(1.0_dp-vreal)
          r = zim*p_i - zjm*p_j
          coeff = 1.0_dp / (vreal*dotij)

          call zero_tensor(vten)
          call add_minus_g_term(vten, 16.0_dp*pi_dp*alpha_s*ca*aterm)
          call add_outer_term(vten, 8.0_dp*pi_dp*alpha_s*ca*coeff, r)
          call tensor_to_helicity(vten, eps_parent, vhel)

       else if (is_q_qbar_pair(fi, fj)) then

          if (massive .and. (mass_i > 0.0_dp .or. mass_j > 0.0_dp .or. &
               mass_parent > 0.0_dp)) then
             info = -108
             return
          endif
          if (massive .and. .not.is_fi) then
             call final_velocity_factors(p_i,p_j,p_k,mass_parent,mass_k,vreal,vtilde,velocity_info)
             if (velocity_info.ne.0) then
                info = -108
                return
             endif
          else
             vreal = 1.0_dp
          endif
          zim = zi - 0.5_dp*(1.0_dp-vreal)
          zjm = zj - 0.5_dp*(1.0_dp-vreal)
          r = zim*p_i - zjm*p_j
          coeff = -2.0_dp / (vreal*dotij)

          call zero_tensor(vten)
          if (is_u1_gluon(fp)) then
             call add_minus_g_term(vten, 8.0_dp*pi_dp*alpha_s*tr_u1/vreal)
             call add_outer_term(vten, 8.0_dp*pi_dp*alpha_s*tr_u1*coeff, r)
          else
             call add_minus_g_term(vten, 8.0_dp*pi_dp*alpha_s*tr/vreal)
             call add_outer_term(vten, 8.0_dp*pi_dp*alpha_s*tr*coeff, r)
          end if
          call tensor_to_helicity(vten, eps_parent, vhel)

       else
          info = -105
          return
       end if
    end if
    if (info.eq.0) then
       if (.not.all(ieee_is_finite(real(vhel,dp))) .or. &
            .not.all(ieee_is_finite(aimag(vhel)))) then
          vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
          info=-20
       endif
    endif
  end subroutine final_splitting_matrix


  subroutine initial_final_splitting_matrix(alpha_s, fa, fj, fp, x, u, p_emit, p_unres, p_spec, mass_spec, &
       eps_parent, vhel, info)
    ! Initial-state emitter, final-state spectator kernels.
    ! p_emit is the real incoming emitter a, p_unres is emitted final j,
    ! p_spec is final spectator k.

    real(dp), intent(in) :: alpha_s, x, u
    integer, intent(in) :: fa, fj, fp
    real(dp), intent(in) :: p_emit(0:3), p_unres(0:3), p_spec(0:3)
    real(dp), intent(in) :: mass_spec
    complex(dp), intent(in) :: eps_parent(0:3,2)
    complex(dp), intent(out) :: vhel(2,2)
    integer, intent(out) :: info

    integer :: ch
    real(dp) :: scalar, denom, dotjk, dotjk_scale, coeff, aterm, parent_tr
    real(dp) :: r(0:3), vten(0:3,0:3)

    info = 0
    vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
    if (mass_spec < 0.0_dp) then
       info = -205
       return
    endif
    call initial_channel(fa, fj, fp, ch)

    select case (ch)

    case (ch_i_qg)

       denom = 1.0_dp - x + u

       if (.not.cs_value_is_resolved(denom,1.0_dp+abs(x)+abs(u))) then
          info = -201
          return
       end if

       scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * (2.0_dp/denom - (1.0_dp + x))
       call scalar_to_vhel(scalar, vhel)

    case (ch_i_gq)

       ! In the crossed initial-state q <- g channel, the reduced Born
       ! amplitude carries one fewer power of N_c than the real ordered
       ! amplitude.  The leading-colour normalization is therefore C_F^LC,
       ! not T_R.
       scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * (1.0_dp - 2.0_dp*x*(1.0_dp - x))
       call scalar_to_vhel(scalar, vhel)

    case (ch_i_qq)

       dotjk = dot4(p_unres, p_spec)
       dotjk_scale=cs_dot4_scale(p_unres,p_spec)

       if (x.le.cs_roundoff_tolerance(1.0_dp) .or. &
            u.le.cs_roundoff_tolerance(1.0_dp) .or. &
            1.0_dp-u.le.cs_roundoff_tolerance(1.0_dp) .or. &
            .not.cs_value_is_resolved(dotjk,dotjk_scale)) then
          info = -202
          return
       end if

       r = p_unres/u - p_spec/(1.0_dp - u)
       coeff = ((1.0_dp - x)/x) * (2.0_dp*u*(1.0_dp - u)/dotjk)
       parent_tr = vector_trace(fp)

       call zero_tensor(vten)
       call add_minus_g_term(vten, 8.0_dp*pi_dp*alpha_s*(2.0_dp*parent_tr)*x)
       call add_outer_term(vten, 8.0_dp*pi_dp*alpha_s*(2.0_dp*parent_tr)*coeff, r)
       call tensor_to_helicity(vten, eps_parent, vhel)

    case (ch_i_gg)

       if (is_u1_gluon(fp) .or. is_u1_gluon(fa) .or. is_u1_gluon(fj)) then
          info = -204
          return
       end if

       dotjk = dot4(p_unres, p_spec)
       dotjk_scale=cs_dot4_scale(p_unres,p_spec)
       denom = 1.0_dp - x + u

       if (x.le.cs_roundoff_tolerance(1.0_dp) .or. &
            u.le.cs_roundoff_tolerance(1.0_dp) .or. &
            1.0_dp-u.le.cs_roundoff_tolerance(1.0_dp) .or. &
            .not.cs_value_is_resolved(dotjk,dotjk_scale) .or. &
            .not.cs_value_is_resolved(denom,1.0_dp+abs(x)+abs(u))) then
          info = -203
          return
       end if

       r = p_unres/u - p_spec/(1.0_dp - u)
       aterm = 1.0_dp/denom - 1.0_dp + x*(1.0_dp - x)
       coeff = ((1.0_dp - x)/x) * (u*(1.0_dp - u)/dotjk)

       call zero_tensor(vten)
       call add_minus_g_term(vten, 16.0_dp*pi_dp*alpha_s*ca*aterm)
       call add_outer_term(vten, 16.0_dp*pi_dp*alpha_s*ca*coeff, r)
       call tensor_to_helicity(vten, eps_parent, vhel)

    case default

       info = -204
       return

    end select
    if (info.eq.0) then
       if (.not.all(ieee_is_finite(real(vhel,dp))) .or. &
            .not.all(ieee_is_finite(aimag(vhel)))) then
          vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
          info=-20
       endif
    endif
  end subroutine initial_final_splitting_matrix


  subroutine initial_initial_splitting_matrix(alpha_s, fa, fj, fp, x, p_emit, p_unres, p_spec, eps_parent, vhel, info)
    ! Initial-state emitter, initial-state spectator kernels.
    ! p_emit is real incoming emitter a, p_unres is emitted final j,
    ! p_spec is incoming spectator b.

    real(dp), intent(in) :: alpha_s, x
    integer, intent(in) :: fa, fj, fp
    real(dp), intent(in) :: p_emit(0:3), p_unres(0:3), p_spec(0:3)
    complex(dp), intent(in) :: eps_parent(0:3,2)
    complex(dp), intent(out) :: vhel(2,2)
    integer, intent(out) :: info

    integer :: ch
    real(dp) :: scalar, dotab, dotja, dotjb, coeff, aterm, parent_tr
    real(dp) :: dotab_scale,dotja_scale,dotjb_scale
    real(dp) :: r(0:3), vten(0:3,0:3)

    info = 0
    vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
    call initial_channel(fa, fj, fp, ch)

    select case (ch)

    case (ch_i_qg)

       if (1.0_dp-x.le.cs_roundoff_tolerance(1.0_dp)) then
          info = -301
          return
       end if

       scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * (2.0_dp/(1.0_dp - x) - (1.0_dp + x))
       call scalar_to_vhel(scalar, vhel)

    case (ch_i_gq)

       scalar = 8.0_dp*pi_dp*alpha_s*cf_lc * (1.0_dp - 2.0_dp*x*(1.0_dp - x))
       call scalar_to_vhel(scalar, vhel)

    case (ch_i_qq)

       dotab = dot4(p_emit, p_spec)
       dotja = dot4(p_unres, p_emit)
       dotjb = dot4(p_unres, p_spec)
       dotab_scale=cs_dot4_scale(p_emit,p_spec)
       dotja_scale=cs_dot4_scale(p_unres,p_emit)
       dotjb_scale=cs_dot4_scale(p_unres,p_spec)

       if (x.le.cs_roundoff_tolerance(1.0_dp) .or. &
            .not.cs_value_is_resolved(dotab,dotab_scale) .or. &
            .not.cs_value_is_resolved(dotja,dotja_scale) .or. &
            .not.cs_value_is_resolved(dotjb,dotjb_scale)) then
          info = -302
          return
       end if

       r = p_unres - (dotja/dotab)*p_spec
       coeff = ((1.0_dp - x)/x) * (2.0_dp*dotab/(dotja*dotjb))
       parent_tr = vector_trace(fp)

       call zero_tensor(vten)
       call add_minus_g_term(vten, 8.0_dp*pi_dp*alpha_s*(2.0_dp*parent_tr)*x)
       call add_outer_term(vten, 8.0_dp*pi_dp*alpha_s*(2.0_dp*parent_tr)*coeff, r)
       call tensor_to_helicity(vten, eps_parent, vhel)

    case (ch_i_gg)

       if (is_u1_gluon(fp) .or. is_u1_gluon(fa) .or. is_u1_gluon(fj)) then
          info = -304
          return
       end if

       dotab = dot4(p_emit, p_spec)
       dotja = dot4(p_unres, p_emit)
       dotjb = dot4(p_unres, p_spec)
       dotab_scale=cs_dot4_scale(p_emit,p_spec)
       dotja_scale=cs_dot4_scale(p_unres,p_emit)
       dotjb_scale=cs_dot4_scale(p_unres,p_spec)

       if (x.le.cs_roundoff_tolerance(1.0_dp) .or. &
            1.0_dp-x.le.cs_roundoff_tolerance(1.0_dp) .or. &
            .not.cs_value_is_resolved(dotab,dotab_scale) .or. &
            .not.cs_value_is_resolved(dotja,dotja_scale) .or. &
            .not.cs_value_is_resolved(dotjb,dotjb_scale)) then
          info = -303
          return
       end if

       r = p_unres - (dotja/dotab)*p_spec
       aterm = x/(1.0_dp - x) + x*(1.0_dp - x)
       coeff = ((1.0_dp - x)/x) * (dotab/(dotja*dotjb))

       call zero_tensor(vten)
       call add_minus_g_term(vten, 16.0_dp*pi_dp*alpha_s*ca*aterm)
       call add_outer_term(vten, 16.0_dp*pi_dp*alpha_s*ca*coeff, r)
       call tensor_to_helicity(vten, eps_parent, vhel)

    case default

       info = -304
       return

    end select
    if (info.eq.0) then
       if (.not.all(ieee_is_finite(real(vhel,dp))) .or. &
            .not.all(ieee_is_finite(aimag(vhel)))) then
          vhel=cmplx(0.0_dp,0.0_dp,kind=dp)
          info=-20
       endif
    endif
  end subroutine initial_initial_splitting_matrix


  subroutine initial_channel(fa, fj, fp, ch)
    ! Infer initial-state splitting channel.
    !
    ! fa: real initial-state flavour
    ! fj: unresolved final-state flavour
    ! fp: mapped Born initial-state parent flavour

    integer, intent(in) :: fa, fj, fp
    integer, intent(out) :: ch

    ch = ch_none

    if (is_quark(fp) .and. fa == fp .and. is_gluon(fj)) then

       ! parent q, real incoming q, unresolved g
       ch = ch_i_qg

    else if (is_gluon(fa) .and. same_quark_flavour(fp,fj)) then

       ! parent q, real incoming g, unresolved q/qbar
       ch = ch_i_gq

    else if (is_gluon(fp) .and. same_quark_flavour(fa,fj)) then

       ! parent g, real incoming q/qbar, unresolved q/qbar
       ch = ch_i_qq

    else if (is_gluon(fp) .and. is_gluon(fa) .and. is_gluon(fj)) then

       ! parent g, real incoming g, unresolved g
       ch = ch_i_gg

    end if
  end subroutine initial_channel

end module cs_lc_spin_dipoles
