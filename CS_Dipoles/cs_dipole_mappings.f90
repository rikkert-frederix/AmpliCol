module cs_dipole_mappings
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: tiny_kin = tiny(1.0_dp)
  real(dp), parameter :: cs_roundoff_safety = 128.0_dp
  real(dp), parameter :: cs_quartic_input_limit = &
       0.125_dp*sqrt(sqrt(huge(1.0_dp)))

contains

  function dot4(p, q) result(d)
    real(dp), intent(in) :: p(0:), q(0:)
    real(dp) :: d

    d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot4


  pure real(dp) function cs_dot4_scale(p,q) result(scale)
    ! Absolute scale of the terms entering a Minkowski product.  This
    ! measures how much cancellation is required to obtain dot4(p,q).
    real(dp), intent(in) :: p(0:),q(0:)
    scale=abs(p(0)*q(0))+abs(p(1)*q(1))+abs(p(2)*q(2))+abs(p(3)*q(3))
  end function cs_dot4_scale


  pure real(dp) function cs_roundoff_tolerance(scale) result(tolerance)
    ! A scale-aware lower bound below which an invariant cannot be trusted
    ! in double precision.  The safety factor covers the several arithmetic
    ! operations and subsequent ratios used by the mappings and kernels.
    real(dp), intent(in) :: scale
    tolerance=max(tiny_kin,cs_roundoff_safety*epsilon(1.0_dp)*abs(scale))
  end function cs_roundoff_tolerance


  pure logical function cs_value_is_resolved(value,scale) result(resolved)
    real(dp), intent(in) :: value,scale
    resolved=.false.
    if (.not.ieee_is_finite(value) .or. .not.ieee_is_finite(scale)) return
    resolved=abs(value).gt.cs_roundoff_tolerance(scale)
  end function cs_value_is_resolved


  pure subroutine cs_normalize_unit_interval(value,valid)
    ! Accept only physical unit-interval variables, while snapping the
    ! scale-sized excursions produced by floating-point roundoff.
    real(dp), intent(inout) :: value
    logical, intent(out) :: valid
    real(dp) :: tolerance

    valid=.false.
    if (.not.ieee_is_finite(value)) return
    tolerance=cs_roundoff_safety*epsilon(1.0_dp)*max(1.0_dp,abs(value))
    if (value.lt.-tolerance .or. value.gt.1.0_dp+tolerance) return
    value=max(0.0_dp,min(1.0_dp,value))
    valid=.true.
  end subroutine cs_normalize_unit_interval


  logical function is_initial(idx)
    integer, intent(in) :: idx

    is_initial = idx <= 2
  end function is_initial


  subroutine cs_dipole_cut_variable(p, ijk, mass_real, mass_parent, cut_variable, info)
    ! Return the variable restricted by an alpha cut.  The variable must
    ! vanish in the unresolved limit selected by the dipole.
    real(dp), intent(in) :: p(0:,:), mass_real(:), mass_parent
    integer, intent(in) :: ijk(3)
    real(dp), intent(out) :: cut_variable
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(dp) :: den, den_scale, x, mi2, mj2, mk2, parent2
    real(dp) :: q2, q2_scale, mui2, muk, muk2, yplus, endpoint_den
    real(dp) :: dij, dik, djk, dij_scale, dik_scale, djk_scale
    logical :: valid

    cut_variable=0.0_dp
    info=0
    n=size(p,2)
    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
    if (size(p,1) /= 4 .or. size(mass_real) /= n) then
       info=-1
       return
    endif
    if (.not.all(ieee_is_finite(p)) .or. .not.all(ieee_is_finite(mass_real)) .or. &
         .not.ieee_is_finite(mass_parent)) then
       info=-20
       return
    endif
    if (maxval(abs(p)).gt.cs_quartic_input_limit .or. &
         max(maxval(abs(mass_real)),abs(mass_parent)).gt.cs_quartic_input_limit) then
       info=-20
       return
    endif
    if (i < 1 .or. i > n .or. j < 1 .or. j > n .or. k < 1 .or. k > n) then
       info=-2
       return
    endif
    if (i == j .or. i == k .or. j == k) then
       info=-3
       return
    endif
    if (j <= 2) then
       info=-4
       return
    endif
    if (any(mass_real < 0.0_dp) .or. mass_parent < 0.0_dp) then
       info=-4
       return
    endif
    if (mass_real(j) > 100.0_dp*epsilon(1.0_dp)* &
         max(tiny(1.0_dp),abs(mass_parent))) then
       info=-5
       return
    endif

    mi2=mass_real(i)*mass_real(i)
    mj2=mass_real(j)*mass_real(j)
    mk2=mass_real(k)*mass_real(k)
    parent2=mass_parent*mass_parent
    dij=dot4(p(:,i),p(:,j))
    dij_scale=cs_dot4_scale(p(:,i),p(:,j))
    if (i > 2 .and. k > 2) then
       dik=dot4(p(:,i),p(:,k))
       djk=dot4(p(:,j),p(:,k))
       dik_scale=cs_dot4_scale(p(:,i),p(:,k))
       djk_scale=cs_dot4_scale(p(:,j),p(:,k))
       den=dij+dik+djk
       den_scale=dij_scale+dik_scale+djk_scale
       if (.not.cs_value_is_resolved(den,den_scale)) then
          info=-10
          return
       endif
       cut_variable=dij/den
       if (parent2 > 0.0_dp .or. mk2 > 0.0_dp) then
          ! For a massive FF dipole alpha is defined relative to the
          ! kinematic endpoint: y_ij,k < alpha*y_+.  This keeps alpha=1 as
          ! the unrestricted case and matches the massive integrated
          ! finite terms.
          q2=parent2+mk2+2.0_dp*den
          q2_scale=parent2+mk2+2.0_dp*den_scale
          if (.not.ieee_is_finite(q2) .or. .not.ieee_is_finite(q2_scale)) then
             info=-20
             return
          endif
          if (q2.le.cs_roundoff_tolerance(q2_scale)) then
             info=-14
             return
          endif
          mui2=parent2/q2
          muk2=mk2/q2
          muk=sqrt(muk2)
          endpoint_den=1.0_dp-mui2-muk2
          if (endpoint_den.le.cs_roundoff_tolerance(1.0_dp+abs(mui2)+abs(muk2))) then
             info=-14
             return
          endif
          yplus=1.0_dp-2.0_dp*muk*(1.0_dp-muk)/&
               endpoint_den
          if (.not.ieee_is_finite(yplus)) then
             info=-20
             return
          endif
          if (yplus.le.cs_roundoff_tolerance(1.0_dp+abs(yplus))) then
             info=-14
             return
          endif
          cut_variable=cut_variable/yplus
       endif
    elseif (i > 2) then
       dik=dot4(p(:,k),p(:,i))
       djk=dot4(p(:,k),p(:,j))
       dik_scale=cs_dot4_scale(p(:,k),p(:,i))
       djk_scale=cs_dot4_scale(p(:,k),p(:,j))
       den=dik+djk
       den_scale=dik_scale+djk_scale
       if (.not.cs_value_is_resolved(den,den_scale)) then
          info=-11
          return
       endif
       x=(den-dij+0.5_dp*(parent2-mi2-mj2))/den
       cut_variable=1.0_dp-x
    elseif (k > 2) then
       dik=dot4(p(:,i),p(:,k))
       dik_scale=cs_dot4_scale(p(:,i),p(:,k))
       den=dij+dik
       den_scale=dij_scale+dik_scale
       if (.not.cs_value_is_resolved(den,den_scale)) then
          info=-12
          return
       endif
       cut_variable=dij/den
    else
       den=dot4(p(:,i),p(:,k))
       den_scale=cs_dot4_scale(p(:,i),p(:,k))
       if (.not.cs_value_is_resolved(den,den_scale)) then
          info=-13
          return
       endif
       ! Initial--initial restriction variable
       !
       !   v_j = p_i.p_j / p_i.p_k .
       !
       ! Do not use v_j/(1+v_j) here.  That is a different variable and
       ! makes the local restriction inconsistent with the analytically
       ! integrated initial--initial dipole.
       cut_variable=dij/den
    endif
    call cs_normalize_unit_interval(cut_variable,valid)
    if (.not.valid) then
       cut_variable=0.0_dp
       info=-20
    endif
  end subroutine cs_dipole_cut_variable


  integer function cs_dipole_topology(ijk)
    integer, intent(in) :: ijk(3)
    if (ijk(1) > 2) then
       if (ijk(3) > 2) then
          cs_dipole_topology=1 ! FF
       else
          cs_dipole_topology=2 ! FI
       endif
    else
       if (ijk(3) > 2) then
          cs_dipole_topology=3 ! IF
       else
          cs_dipole_topology=4 ! II
       endif
    endif
  end function cs_dipole_topology


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

    real(dp) :: Ksq, Ksq_scale, Qsum(0:3), Qsq, Qsq_scale
    real(dp) :: acoef, bcoef

    info = 0
    qtilde = 0.0_dp

    if (size(qvec).ne.4 .or. size(Kvec).ne.4 .or. size(Ktvec).ne.4 .or. &
         size(qtilde).ne.4) then
       info=-1
       return
    endif
    if (.not.all(ieee_is_finite(qvec)) .or. .not.all(ieee_is_finite(Kvec)) .or. &
         .not.all(ieee_is_finite(Ktvec))) then
       info=-20
       return
    endif

    Qsum = Kvec + Ktvec
    Ksq  = dot4(Kvec, Kvec)
    Qsq  = dot4(Qsum, Qsum)
    Ksq_scale=cs_dot4_scale(Kvec,Kvec)
    Qsq_scale=cs_dot4_scale(Qsum,Qsum)

    if (.not.cs_value_is_resolved(Ksq,Ksq_scale) .or. &
         .not.cs_value_is_resolved(Qsq,Qsq_scale)) then
       info = -20
       return
    end if
    if (Ksq.le.0.0_dp .or. Qsq.le.0.0_dp) then
       info=-20
       return
    endif

    acoef = 2.0_dp * dot4(qvec, Qsum) / Qsq
    bcoef = 2.0_dp * dot4(qvec, Kvec) / Ksq

    qtilde = qvec - acoef*Qsum + bcoef*Ktvec
    if (.not.all(ieee_is_finite(qtilde))) then
       qtilde=0.0_dp
       info=-20
    endif
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
    !   -20  non-finite or numerically unresolved kinematics

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

    real(dp) :: x, y, den, den_scale, omx, omy
    real(dp) :: Kvec(0:3), Ktvec(0:3), tmp(0:3)
    real(dp) :: mi2, mj2, mk2, mij2, q2, pij2, lambda_real, lambda_born
    real(dp) :: q2_scale, lambda_real_scale, lambda_born_scale
    real(dp) :: map_scale, qdotk, qvec(0:3), pairvec(0:3)
    real(dp) :: dij, dik, djk, dij_scale, dik_scale, djk_scale
    logical :: use_masses,valid

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
       if (size(mass_real) /= n) then
          istat = -5
          goto 900
       endif
       if (.not.all(ieee_is_finite(mass_real)) .or. .not.ieee_is_finite(mass_parent)) then
          istat=-20
          goto 900
       endif
       if (any(mass_real < 0.0_dp) .or. mass_parent < 0.0_dp) then
          istat = -5
          goto 900
       endif
    endif

    if (size(p,1) /= 4 .or. size(ptilde,1) /= 4 .or. nm /= n-1) then
       istat = -1
       goto 900
    end if
    if (.not.all(ieee_is_finite(p))) then
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
       if (mass_real(j) > 100.0_dp*epsilon(1.0_dp)* &
            max(tiny(1.0_dp),abs(mass_parent))) then
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
    dij=dot4(p(:,i),p(:,j))
    dij_scale=cs_dot4_scale(p(:,i),p(:,j))

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

       dik=dot4(p(:,i),p(:,k))
       djk=dot4(p(:,j),p(:,k))
       dik_scale=cs_dot4_scale(p(:,i),p(:,k))
       djk_scale=cs_dot4_scale(p(:,j),p(:,k))
       den=dij+dik+djk
       den_scale=dij_scale+dik_scale+djk_scale
       if (.not.cs_value_is_resolved(den,den_scale)) then
          istat = -10
          goto 900
       endif
       y=dij/den
       call cs_normalize_unit_interval(y,valid)
       if (.not.valid) then
          istat=-10
          goto 900
       endif

       if (use_masses .and. (mi2 > 0.0_dp .or. mk2 > 0.0_dp .or. mij2 > 0.0_dp)) then
          qvec = p(:,i) + p(:,j) + p(:,k)
          pairvec = p(:,i) + p(:,j)
          q2 = dot4(qvec,qvec)
          q2_scale=cs_dot4_scale(qvec,qvec)
          pij2 = dot4(pairvec,pairvec)
          lambda_real = q2*q2 + pij2*pij2 + mk2*mk2 - 2.0_dp*(q2*pij2 + q2*mk2 + pij2*mk2)
          lambda_born = q2*q2 + mij2*mij2 + mk2*mk2 - 2.0_dp*(q2*mij2 + q2*mk2 + mij2*mk2)
          lambda_real_scale=abs(q2*q2)+abs(pij2*pij2)+abs(mk2*mk2)+&
               2.0_dp*(abs(q2*pij2)+abs(q2*mk2)+abs(pij2*mk2))
          lambda_born_scale=abs(q2*q2)+abs(mij2*mij2)+abs(mk2*mk2)+&
               2.0_dp*(abs(q2*mij2)+abs(q2*mk2)+abs(mij2*mk2))
          if (.not.ieee_is_finite(q2) .or. .not.ieee_is_finite(q2_scale)) then
             istat=-20
             goto 900
          endif
          if (q2.le.cs_roundoff_tolerance(q2_scale)) then
             istat = -10
             goto 900
          endif
          if (.not.ieee_is_finite(lambda_real) .or. .not.ieee_is_finite(lambda_born) .or. &
               .not.ieee_is_finite(lambda_real_scale) .or. &
               .not.ieee_is_finite(lambda_born_scale)) then
             istat=-20
             goto 900
          endif
          if (lambda_real.le.cs_roundoff_tolerance(lambda_real_scale) .or. &
               lambda_born.lt.-cs_roundoff_tolerance(lambda_born_scale)) then
             istat=-10
             goto 900
          endif
          map_scale = sqrt(max(0.0_dp,lambda_born)/lambda_real)
          qdotk = dot4(qvec,p(:,k))
          ptilde(:,nk) = map_scale*(p(:,k) - qdotk/q2*qvec) + &
               (q2 + mk2 - mij2)/(2.0_dp*q2)*qvec
          ptilde(:,ni) = qvec - ptilde(:,nk)
          if (present(yout)) yout = y
          goto 900
       endif

       if (present(yout)) yout = y

       omy = 1.0_dp - y

       if (.not.cs_value_is_resolved(omy,1.0_dp+abs(y))) then
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

       dik=dot4(p(:,k),p(:,i))
       djk=dot4(p(:,k),p(:,j))
       dik_scale=cs_dot4_scale(p(:,k),p(:,i))
       djk_scale=cs_dot4_scale(p(:,k),p(:,j))
       den=dik+djk
       den_scale=dik_scale+djk_scale

       if (.not.cs_value_is_resolved(den,den_scale)) then
          istat = -11
          goto 900
       end if

       if (use_masses .and. (mi2 > 0.0_dp .or. mij2 > 0.0_dp)) then
          x = (den - dot4(p(:,i),p(:,j)) + 0.5_dp*(mij2-mi2-mj2))/den
       else
          x = (den - dij) / den
       endif
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-11
          goto 900
       endif
       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat=-11
          goto 900
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

       dik=dot4(p(:,i),p(:,k))
       dik_scale=cs_dot4_scale(p(:,i),p(:,k))
       den=dij+dik
       den_scale=dij_scale+dik_scale

       if (.not.cs_value_is_resolved(den,den_scale)) then
          istat = -12
          goto 900
       end if

       x = (den - dot4(p(:,j), p(:,k))) / den
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-12
          goto 900
       endif
       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat=-12
          goto 900
       endif
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

       den=dot4(p(:,i),p(:,k))
       den_scale=cs_dot4_scale(p(:,i),p(:,k))

       if (.not.cs_value_is_resolved(den,den_scale)) then
          istat = -13
          goto 900
       end if

       x = (den - dij - dot4(p(:,j), p(:,k))) / den
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then
          istat=-13
          goto 900
       endif
       if (x.le.cs_roundoff_tolerance(1.0_dp)) then
          istat=-13
          goto 900
       endif
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
    if (istat.eq.0) then
       if (.not.all(ieee_is_finite(ptilde))) istat=-20
       if (present(xout)) then
          if (.not.ieee_is_finite(xout)) istat=-20
       endif
       if (present(yout)) then
          if (.not.ieee_is_finite(yout)) istat=-20
       endif
    endif
    if (istat.ne.0) then
       ptilde=0.0_dp
       if (present(xout)) xout=1.0_dp
       if (present(yout)) yout=0.0_dp
    endif
    if (present(info)) info = istat

  end subroutine cs_map


  subroutine cs_born_pushback_weight(p, ptilde, ijk, mass_real, mass_parent, weight, info)
    ! Reciprocal four-dimensional CS radiation measure used to embed a
    ! mapped Born contribution in an (m+1)-body integration.  The II
    ! expression uses vbar=vtilde/(1-x), so its (1-x) conversion is part
    ! of the Jacobian.  No local-dipole 1/x factor belongs here.  For a
    ! massive emitter or spectator this is only the local measure: the
    ! mass-dependent radiation-domain volume is not included, so those
    ! mappings are not used as unit-normalized Born-recycling histories.
    real(dp), intent(in) :: p(0:,:), ptilde(0:,:), mass_real(:), mass_parent
    integer, intent(in) :: ijk(3)
    real(dp), intent(out) :: weight
    integer, intent(out) :: info
    integer :: i,j,k,ni,nk,n
    real(dp) :: den,den_scale,x,y,jac,pi
    real(dp) :: mi,mj,mk,mp,q(0:3),q2,q2_scale,a,lam
    real(dp) :: dij,dik,djk,dij_scale,dik_scale,djk_scale
    logical :: valid

    weight=0.0_dp
    info=0
    pi=acos(-1.0_dp)
    i=ijk(1); j=ijk(2); k=ijk(3)
    n=size(p,2)
    if (size(p,1).ne.4 .or. size(ptilde,1).ne.4 .or. n.lt.3 .or. &
         size(ptilde,2).ne.n-1 .or. size(mass_real).ne.n) then
       info=-1; return
    endif
    if (.not.all(ieee_is_finite(p)) .or. .not.all(ieee_is_finite(ptilde)) .or. &
         .not.all(ieee_is_finite(mass_real)) .or. .not.ieee_is_finite(mass_parent)) then
       info=-20; return
    endif
    if (max(maxval(abs(p)),maxval(abs(ptilde))).gt.cs_quartic_input_limit .or. &
         max(maxval(abs(mass_real)),abs(mass_parent)).gt.cs_quartic_input_limit) then
       info=-20; return
    endif
    if (i.lt.1 .or. i.gt.n .or. j.lt.1 .or. j.gt.n .or. k.lt.1 .or. k.gt.n) then
       info=-2; return
    endif
    if (i.eq.j .or. i.eq.k .or. j.eq.k) then
       info=-3; return
    endif
    if (j.le.2) then
       info=-4; return
    endif
    if (any(mass_real.lt.0.0_dp) .or. mass_parent.lt.0.0_dp) then
       info=-4; return
    endif
    ni=new_index(i,j); nk=new_index(k,j)
    if (ni.lt.1 .or. ni.gt.n-1 .or. nk.lt.1 .or. nk.gt.n-1) then
       info=-5; return
    endif

    if (i.gt.2 .and. k.gt.2) then
       dij=dot4(p(:,i),p(:,j)); dik=dot4(p(:,i),p(:,k)); djk=dot4(p(:,j),p(:,k))
       dij_scale=cs_dot4_scale(p(:,i),p(:,j))
       dik_scale=cs_dot4_scale(p(:,i),p(:,k))
       djk_scale=cs_dot4_scale(p(:,j),p(:,k))
       den=dij+dik+djk
       den_scale=dij_scale+dik_scale+djk_scale
       if (.not.cs_value_is_resolved(den,den_scale) .or. den.le.0.0_dp) then
          info=-10; return
       endif
       y=dij/den
       call cs_normalize_unit_interval(y,valid)
       if (.not.valid) then; info=-10; return; endif
       mi=mass_real(i); mj=mass_real(j); mk=mass_real(k); mp=mass_parent
       if (mi.gt.0.0_dp .or. mj.gt.0.0_dp .or. mk.gt.0.0_dp .or. mp.gt.0.0_dp) then
          q=p(:,i)+p(:,j)+p(:,k); q2=dot4(q,q)
          q2_scale=cs_dot4_scale(q,q)
          if (.not.ieee_is_finite(q2) .or. .not.ieee_is_finite(q2_scale)) then
             info=-20; return
          endif
          if (q2.le.cs_roundoff_tolerance(q2_scale)) then; info=-11; return; endif
          mi=mi/sqrt(q2); mj=mj/sqrt(q2); mk=mk/sqrt(q2); mp=mp/sqrt(q2)
          a=1.0_dp-mi*mi-mj*mj-mk*mk
          lam=1.0_dp+mp**4+mk**4-2.0_dp*(mp*mp+mk*mk+mp*mp*mk*mk)
          if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(lam)) then
             info=-20; return
          endif
          if (a.le.cs_roundoff_tolerance(1.0_dp+mi*mi+mj*mj+mk*mk) .or. &
               lam.le.cs_roundoff_tolerance(1.0_dp+mp**4+mk**4) .or. &
               1.0_dp-y.le.cs_roundoff_tolerance(1.0_dp+abs(y))) then
             info=-12; return
          endif
          jac=q2/(16.0_dp*pi*pi)*a*a/sqrt(lam)*(1.0_dp-y)
       else
          jac=2.0_dp*dot4(ptilde(:,ni),ptilde(:,nk))/(16.0_dp*pi*pi)*(1.0_dp-y)
       endif
    elseif (i.gt.2) then
       jac=2.0_dp*dot4(ptilde(:,ni),p(:,k))/(16.0_dp*pi*pi)
    elseif (k.gt.2) then
       jac=2.0_dp*dot4(ptilde(:,nk),p(:,i))/(16.0_dp*pi*pi)
    else
       den=dot4(p(:,i),p(:,k))
       den_scale=cs_dot4_scale(p(:,i),p(:,k))
       if (.not.cs_value_is_resolved(den,den_scale)) then; info=-13; return; endif
       x=(den-dot4(p(:,j),p(:,i))-dot4(p(:,j),p(:,k)))/den
       call cs_normalize_unit_interval(x,valid)
       if (.not.valid) then; info=-13; return; endif
       jac=2.0_dp*den/(16.0_dp*pi*pi)*(1.0_dp-x)
    endif
    if (.not.ieee_is_finite(jac)) then; info=-20; return; endif
    if (jac.le.tiny_kin) then; info=-20; return; endif
    weight=1.0_dp/jac
    if (.not.ieee_is_finite(weight)) then
       weight=0.0_dp
       info=-20
    endif
  end subroutine cs_born_pushback_weight

end module cs_dipole_mappings
