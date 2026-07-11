module cs_lc_spin_dipoles
  use cs_dipole_mappings, only: dp, dot4, new_index
  implicit none

  real(dp), parameter :: pi_dp = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: ca = 3.0_dp
  real(dp), parameter :: cf_lc = ca / 2.0_dp
  real(dp), parameter :: tr = 0.5_dp
  real(dp), parameter :: tr_u1 = tr / ca
  real(dp), parameter :: tiny_dip = 1.0d-30

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
    is_quark = .not.is_gluon(f)
  end function is_quark


  logical function is_q_qbar_pair(f1, f2)
    integer, intent(in) :: f1, f2
    is_q_qbar_pair = (is_quark(f1) .and. f1 == -f2)
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
    integer, intent(out), optional :: info
    real(dp), intent(in), optional :: lc_weight
    real(dp), intent(in), optional :: mass_real(:), mass_parent

    integer :: n, ni, nk, istat,i,j,k
    logical :: iini, kini
    real(dp) :: wt, pref, vcontract_alt
    real(dp) :: sij, sik, sjk, dotij, sij_parent
    real(dp) :: mi, mj, mk, mparent
    real(dp) :: x, y, z, u
    complex(dp) :: vhel(2,2)
    real(dp) :: vcontract
    logical :: use_masses

    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
    
    dip = 0.0_dp
    istat = 0

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
       if (size(mass_real) /= n .or. any(mass_real < 0.0_dp) .or. mass_parent < 0.0_dp) then
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
       if (mj > 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,mparent)) then
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

    sij = 2.0_dp * dot4(p(:,i), p(:,j))
    dotij = dot4(p(:,i), p(:,j))
    sij_parent = dot4(p(:,i)+p(:,j),p(:,i)+p(:,j))-mparent*mparent

    if (abs(sij_parent) <= tiny_dip .or. abs(dotij) <= tiny_dip) then
       istat = -10
       goto 900
    end if

    if (.not. iini .and. .not. kini) then
       ! FF: final emitter, final unresolved, final spectator.

       sik = 2.0_dp * dot4(p(:,i), p(:,k))
       sjk = 2.0_dp * dot4(p(:,j), p(:,k))

       if (abs(sik + sjk) <= tiny_dip .or. abs(sij + sik + sjk) <= tiny_dip) then
          istat = -11
          goto 900
       end if

       z = sik / (sik + sjk)
       y = sij / (sij + sik + sjk)

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

       if (abs(sik + sjk) <= tiny_dip) then
          istat = -12
          goto 900
       end if

       z = sik / (sik + sjk)
       x = (sik + sjk - sij) / (sik + sjk)

       if (abs(x) <= tiny_dip) then
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

       if (abs(sij + sik) <= tiny_dip) then
          istat = -14
          goto 900
       end if

       x = (sij + sik - sjk) / (sij + sik)

       if (abs(x) <= tiny_dip) then
          istat = -15
          goto 900
       end if

       u = sij / (sij + sik)

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

       if (abs(sik) <= tiny_dip) then
          istat = -16
          goto 900
       end if

       x = (sik - sij - sjk) / sik

       if (abs(x) <= tiny_dip) then
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
    if (present(info)) info = istat
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
    real(dp) :: q2, pij2, mk2, mp2, lambda_real, lambda_born
    real(dp) :: zim, zjm
    real(dp) :: r(0:3), vten(0:3,0:3)
    logical :: massive

    info = 0
    zi = z
    zj = 1.0_dp - z
    dotij = dot4(p_i, p_j)
    massive = mass_i > 100.0_dp*epsilon(1.0_dp) .or. mass_j > 100.0_dp*epsilon(1.0_dp) .or. &
         mass_k > 100.0_dp*epsilon(1.0_dp) .or. mass_parent > 100.0_dp*epsilon(1.0_dp)

    if (abs(dotij) <= tiny_dip) then
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

       if (abs(denom_i) <= tiny_dip) then
          info = -103
          return
       end if

       if (massive .and. .not.is_fi) then
          q2 = dot4(p_i+p_j+p_k,p_i+p_j+p_k)
          pij2 = dot4(p_i+p_j,p_i+p_j)
          mk2 = mass_k*mass_k
          mp2 = mass_parent*mass_parent
          lambda_real = q2*q2 + pij2*pij2 + mk2*mk2 - 2.0_dp*(q2*pij2 + q2*mk2 + pij2*mk2)
          lambda_born = q2*q2 + mp2*mp2 + mk2*mk2 - 2.0_dp*(q2*mp2 + q2*mk2 + mp2*mk2)
          if (q2 <= tiny_dip .or. lambda_real <= tiny_dip .or. lambda_born < -tiny_dip) then
             info = -106
             return
          endif
          vreal = sqrt(max(0.0_dp,lambda_real))/(q2-pij2-mk2)
          vtilde = sqrt(max(0.0_dp,lambda_born))/(q2-mp2-mk2)
          if (abs(vreal) <= tiny_dip) then
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

          if (abs(denom_i) <= tiny_dip .or. abs(denom_j) <= tiny_dip) then
             info = -104
             return
          end if

          if (massive .and. (mass_i > 100.0_dp*epsilon(1.0_dp) .or. &
               mass_j > 100.0_dp*epsilon(1.0_dp) .or. mass_parent > 100.0_dp*epsilon(1.0_dp))) then
             info = -107
             return
          endif
          ! Ordered LC dipoles designate leg j as unresolved.  The
          ! complementary i-soft pole is supplied by the dipole where i is
          ! the unresolved leg; keeping it here would over-subtract soft
          ! gluons in a colour-ordered sum.
          aterm = 1.0_dp/denom_i - 1.0_dp
          if (massive .and. .not.is_fi) then
             q2 = dot4(p_i+p_j+p_k,p_i+p_j+p_k)
             pij2 = dot4(p_i+p_j,p_i+p_j)
             mk2 = mass_k*mass_k
             lambda_real = q2*q2 + pij2*pij2 + mk2*mk2 - 2.0_dp*(q2*pij2 + q2*mk2 + pij2*mk2)
             if (q2 <= tiny_dip .or. lambda_real <= tiny_dip) then
                info = -107
                return
             endif
             vreal = sqrt(lambda_real)/(q2-pij2-mk2)
             if (abs(vreal) <= tiny_dip) then
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

          if (massive .and. (mass_i > 100.0_dp*epsilon(1.0_dp) .or. &
               mass_j > 100.0_dp*epsilon(1.0_dp) .or. mass_parent > 100.0_dp*epsilon(1.0_dp))) then
             info = -108
             return
          endif
          if (massive .and. .not.is_fi) then
             q2 = dot4(p_i+p_j+p_k,p_i+p_j+p_k)
             pij2 = dot4(p_i+p_j,p_i+p_j)
             mk2 = mass_k*mass_k
             lambda_real = q2*q2 + pij2*pij2 + mk2*mk2 - 2.0_dp*(q2*pij2 + q2*mk2 + pij2*mk2)
             if (q2 <= tiny_dip .or. lambda_real <= tiny_dip) then
                info = -108
                return
             endif
             vreal = sqrt(lambda_real)/(q2-pij2-mk2)
             if (abs(vreal) <= tiny_dip) then
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
    real(dp) :: scalar, denom, dotjk, coeff, aterm, parent_tr
    real(dp) :: r(0:3), vten(0:3,0:3)

    info = 0
    if (mass_spec < 0.0_dp) then
       info = -205
       return
    endif
    call initial_channel(fa, fj, fp, ch)

    select case (ch)

    case (ch_i_qg)

       denom = 1.0_dp - x + u

       if (abs(denom) <= tiny_dip) then
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

       if (abs(x) <= tiny_dip .or. abs(u) <= tiny_dip .or. abs(1.0_dp - u) <= tiny_dip .or. &
           abs(dotjk) <= tiny_dip) then
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
       denom = 1.0_dp - x + u

       if (abs(x) <= tiny_dip .or. abs(u) <= tiny_dip .or. abs(1.0_dp - u) <= tiny_dip .or. &
           abs(dotjk) <= tiny_dip .or. abs(denom) <= tiny_dip) then
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
    real(dp) :: r(0:3), vten(0:3,0:3)

    info = 0
    call initial_channel(fa, fj, fp, ch)

    select case (ch)

    case (ch_i_qg)

       if (abs(1.0_dp - x) <= tiny_dip) then
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

       if (abs(x) <= tiny_dip .or. abs(dotab) <= tiny_dip .or. &
           abs(dotja) <= tiny_dip .or. abs(dotjb) <= tiny_dip) then
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

       if (abs(x) <= tiny_dip .or. abs(1.0_dp - x) <= tiny_dip .or. &
           abs(dotab) <= tiny_dip .or. abs(dotja) <= tiny_dip .or. abs(dotjb) <= tiny_dip) then
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

    else if (is_quark(fp) .and. is_gluon(fa) .and. is_quark(fj) .and. abs(fp) == abs(fj)) then

       ! parent q, real incoming g, unresolved q/qbar
       ch = ch_i_gq

    else if (is_gluon(fp) .and. is_quark(fa) .and. is_quark(fj) .and. abs(fa) == abs(fj)) then

       ! parent g, real incoming q/qbar, unresolved q/qbar
       ch = ch_i_qq

    else if (is_gluon(fp) .and. is_gluon(fa) .and. is_gluon(fj)) then

       ! parent g, real incoming g, unresolved g
       ch = ch_i_gg

    end if
  end subroutine initial_channel

end module cs_lc_spin_dipoles
