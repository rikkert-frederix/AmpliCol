module cs_dipole_mappings
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: tiny_kin = 1.0d-30

contains

  function dot4(p, q) result(d)
    real(dp), intent(in) :: p(0:), q(0:)
    real(dp) :: d

    d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot4


  logical function is_initial(idx)
    integer, intent(in) :: idx

    is_initial = idx <= 2
  end function is_initial


  function new_index(old_index, removed_index) result(idx)
    integer, intent(in) :: old_index, removed_index
    integer :: idx

    if (old_index < removed_index) then
       idx = old_index
    else if (old_index > removed_index) then
       idx = old_index - 1
    else
       idx = 0
    end if
  end function new_index


  subroutine boost_K_to_Ktilde(qvec, Kvec, Ktvec, qtilde, info)
    real(dp), intent(in)  :: qvec(0:), Kvec(0:), Ktvec(0:)
    real(dp), intent(out) :: qtilde(0:)
    integer, intent(out)  :: info

    real(dp) :: Ksq, Qsum(0:3), Qsq
    real(dp) :: acoef, bcoef

    info = 0

    Qsum = Kvec + Ktvec
    Ksq  = dot4(Kvec, Kvec)
    Qsq  = dot4(Qsum, Qsum)

    if (abs(Ksq) <= tiny_kin .or. abs(Qsq) <= tiny_kin) then
       qtilde = 0.0_dp
       info = -20
       return
    end if

    acoef = 2.0_dp * dot4(qvec, Qsum) / Qsq
    bcoef = 2.0_dp * dot4(qvec, Kvec) / Ksq

    qtilde = qvec - acoef*Qsum + bcoef*Ktvec
  end subroutine boost_K_to_Ktilde


  subroutine cs_map(p, ijk, ptilde, info, xout, yout, mass_real, mass_parent)
    ! Catani-Seymour momentum mappings.  The original massless interface is
    ! retained; supplying masses enables the massive final-state mappings.
    !
    ! Convention:
    !   p(0:3,1:nexternal)
    !
    !   p(:,1) and p(:,2) are the two incoming momenta.
    !   All other particles are final-state momenta.
    !
    !   Metric: p.q = p0 q0 - p1 q1 - p2 q2 - p3 q3.
    !
    ! Arguments:
    !   i = emitter
    !   j = unresolved emitted final-state parton
    !   k = spectator
    !
    ! The dipole type is inferred automatically:
    !
    !   i >  2, k >  2  -> FF
    !   i >  2, k <= 2  -> FI
    !   i <= 2, k >  2  -> IF
    !   i <= 2, k <= 2  -> II
    !
    ! The unresolved parton j is removed from ptilde.
    !
    ! Optional outputs:
    !   xout = CS x variable for FI, IF, II mappings
    !   yout = CS y variable for FF mappings
    !
    ! info:
    !    0   success
    !   -1   bad array sizes
    !   -2   bad indices
    !   -3   i,j,k not distinct
    !   -4   unresolved j is not final-state
    !   -10  singular FF denominator
    !   -11  singular FI denominator
    !   -12  singular IF denominator
    !   -13  singular II denominator
    !   -20  singular II Lorentz transform

    real(dp), intent(in)  :: p(0:,:)
    integer, dimension(3),intent(in)  :: ijk
    real(dp), intent(out) :: ptilde(0:,:)
    integer, intent(out), optional :: info
    real(dp), intent(out), optional :: xout, yout
    real(dp), intent(in), optional :: mass_real(:), mass_parent

    integer :: n, nm, i, j, k
    integer :: l, ni, nj, nk, nl
    integer :: istat, locinfo

    logical :: i_is_initial, k_is_initial

    real(dp) :: x, y, den, omx, omy
    real(dp) :: Kvec(0:3), Ktvec(0:3), tmp(0:3)
    real(dp) :: mi2, mj2, mk2, mij2, q2, pij2, lambda_real, lambda_born
    real(dp) :: map_scale, qdotk, qvec(0:3), pairvec(0:3)
    logical :: use_masses

    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
    
    istat = 0
    ptilde = 0.0_dp

    if (present(xout)) xout = 1.0_dp
    if (present(yout)) yout = 0.0_dp

    n  = size(p, 2)
    nm = size(ptilde, 2)

    use_masses = present(mass_real) .or. present(mass_parent)
    if (use_masses) then
       if (.not.(present(mass_real) .and. present(mass_parent))) then
          istat = -5
          goto 900
       endif
       if (size(mass_real) /= n .or. any(mass_real < 0.0_dp) .or. mass_parent < 0.0_dp) then
          istat = -5
          goto 900
       endif
    endif

    if (size(p,1) /= 4 .or. size(ptilde,1) /= 4 .or. nm /= n-1) then
       istat = -1
       goto 900
    end if

    if (i < 1 .or. i > n .or. j < 1 .or. j > n .or. k < 1 .or. k > n) then
       istat = -2
       goto 900
    end if

    if (i == j .or. i == k .or. j == k) then
       istat = -3
       goto 900
    end if

    ! In Catani-Seymour real-emission mappings, j is the unresolved
    ! emitted parton and must be final-state.
    if (j <= 2) then
       istat = -4
       goto 900
    end if

    ni = new_index(i, j)
    nj = new_index(j, j)
    nk = new_index(k, j)

    if (ni == 0 .or. nk == 0 .or. nj /= 0) then
       istat = -3
       goto 900
    end if

    if (use_masses) then
       mi2 = mass_real(i)*mass_real(i)
       mj2 = mass_real(j)*mass_real(j)
       mk2 = mass_real(k)*mass_real(k)
       mij2 = mass_parent*mass_parent
       if (mass_real(j) > 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,mass_parent)) then
          istat = -5
          goto 900
       endif
    else
       mi2 = 0.0_dp
       mj2 = 0.0_dp
       mk2 = 0.0_dp
       mij2 = 0.0_dp
    endif

    i_is_initial = is_initial(i)
    k_is_initial = is_initial(k)

    ! Default compact copy: copy all momenta except unresolved j.
    do l = 1, n
       if (l /= j) then
          nl = new_index(l, j)
          ptilde(:,nl) = p(:,l)
       end if
    end do


    if (.not. i_is_initial .and. .not. k_is_initial) then

       ! FF: final emitter, final unresolved, final spectator.
       !
       ! y_ij,k = p_i.p_j / (p_i.p_j + p_i.p_k + p_j.p_k)
       !
       ! ptilde_k  = p_k / (1 - y)
       ! ptilde_ij = p_i + p_j - y/(1-y) p_k

       if (use_masses .and. (mi2 > 0.0_dp .or. mk2 > 0.0_dp .or. mij2 > 0.0_dp)) then
          qvec = p(:,i) + p(:,j) + p(:,k)
          pairvec = p(:,i) + p(:,j)
          q2 = dot4(qvec,qvec)
          pij2 = dot4(pairvec,pairvec)
          lambda_real = q2*q2 + pij2*pij2 + mk2*mk2 - 2.0_dp*(q2*pij2 + q2*mk2 + pij2*mk2)
          lambda_born = q2*q2 + mij2*mij2 + mk2*mk2 - 2.0_dp*(q2*mij2 + q2*mk2 + mij2*mk2)
          if (q2 <= tiny_kin .or. lambda_real <= tiny_kin .or. lambda_born < -tiny_kin) then
             istat = -10
             goto 900
          endif
          map_scale = sqrt(max(0.0_dp,lambda_born)/lambda_real)
          qdotk = dot4(qvec,p(:,k))
          ptilde(:,nk) = map_scale*(p(:,k) - qdotk/q2*qvec) + &
               (q2 + mk2 - mij2)/(2.0_dp*q2)*qvec
          ptilde(:,ni) = qvec - ptilde(:,nk)
          y = dot4(p(:,i),p(:,j))/(dot4(p(:,i),p(:,j)) + dot4(p(:,i),p(:,k)) + dot4(p(:,j),p(:,k)))
          if (present(yout)) yout = y
          goto 900
       endif

       den = dot4(p(:,i), p(:,j)) + dot4(p(:,i), p(:,k)) + dot4(p(:,j), p(:,k))

       if (abs(den) <= tiny_kin) then
          istat = -10
          goto 900
       end if

       y = dot4(p(:,i), p(:,j)) / den
       if (present(yout)) yout = y

       omy = 1.0_dp - y

       if (abs(omy) <= tiny_kin) then
          istat = -10
          goto 900
       end if

       ptilde(:,ni) = p(:,i) + p(:,j) - (y/omy)*p(:,k)
       ptilde(:,nk) = p(:,k) / omy


    else if (.not. i_is_initial .and. k_is_initial) then

       ! FI: final emitter, final unresolved, initial spectator.
       !
       ! x_ij,k = 1 - p_i.p_j / [p_k.(p_i+p_j)]
       !
       ! ptilde_k  = x p_k
       ! ptilde_ij = p_i + p_j - (1-x) p_k

       den = dot4(p(:,k), p(:,i)) + dot4(p(:,k), p(:,j))

       if (abs(den) <= tiny_kin) then
          istat = -11
          goto 900
       end if

       if (use_masses .and. (mi2 > 0.0_dp .or. mij2 > 0.0_dp)) then
          x = (den - dot4(p(:,i),p(:,j)) + 0.5_dp*(mij2-mi2-mj2))/den
       else
          x = (den - dot4(p(:,i), p(:,j))) / den
       endif
       if (present(xout)) xout = x

       omx = 1.0_dp - x

       ptilde(:,nk) = x*p(:,k)
       ptilde(:,ni) = p(:,i) + p(:,j) - omx*p(:,k)


    else if (i_is_initial .and. .not. k_is_initial) then

       ! IF: initial emitter, final unresolved, final spectator.
       !
       ! x_ij,k = 1 - p_j.p_k / [p_i.(p_j+p_k)]
       !
       ! ptilde_i = x p_i
       ! ptilde_k = p_j + p_k - (1-x) p_i

       den = dot4(p(:,i), p(:,j)) + dot4(p(:,i), p(:,k))

       if (abs(den) <= tiny_kin) then
          istat = -12
          goto 900
       end if

       x = (den - dot4(p(:,j), p(:,k))) / den
       if (present(xout)) xout = x

       omx = 1.0_dp - x

       ptilde(:,ni) = x*p(:,i)
       ptilde(:,nk) = p(:,j) + p(:,k) - omx*p(:,i)


    else

       ! II: initial emitter, final unresolved, initial spectator.
       !
       ! x_ij,k = 1 - [p_j.p_i + p_j.p_k] / [p_i.p_k]
       !
       ! ptilde_i = x p_i
       ! ptilde_k = p_k
       !
       ! All remaining final-state momenta are transformed by the Lorentz
       ! transformation mapping
       !
       !   K      = p_i + p_k - p_j
       !
       ! to
       !
       !   Ktilde = x p_i + p_k.

       den = dot4(p(:,i), p(:,k))

       if (abs(den) <= tiny_kin) then
          istat = -13
          goto 900
       end if

       x = (den - dot4(p(:,j), p(:,i)) - dot4(p(:,j), p(:,k))) / den
       if (present(xout)) xout = x

       ptilde(:,ni) = x*p(:,i)
       ptilde(:,nk) = p(:,k)

       Kvec  = p(:,i) + p(:,k) - p(:,j)
       Ktvec = x*p(:,i) + p(:,k)

       do l = 1, n
          if (l /= i .and. l /= j .and. l /= k) then
             nl = new_index(l, j)

             call boost_K_to_Ktilde(p(:,l), Kvec, Ktvec, tmp, locinfo)

             if (locinfo /= 0) then
                istat = locinfo
                goto 900
             end if

             ptilde(:,nl) = tmp
          end if
       end do

    end if

900 continue
    if (present(info)) info = istat

  end subroutine cs_map

end module cs_dipole_mappings
