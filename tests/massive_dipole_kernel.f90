program massive_dipole_kernel
  use cs_lc_spin_dipoles
  use cs_dipole_mappings, only: dot4,cs_dot4_scale,cs_value_is_resolved
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none

  real(dp) :: p(0:3,5), masses(5), dip,dotij,dot_scale,nan_value
  real(dp) :: pscaled(0:3,5),masses_scaled(5),dip_reference
  complex(dp) :: rho(2,2), eps(0:3,2)
  integer :: info,minimum_integer
  integer :: flav_q(5), born_q(4), flav_g(5), born_g(4)

  rho = cmplx(0.0_dp,0.0_dp,kind=dp)
  rho(1,1) = cmplx(1.0_dp,0.0_dp,kind=dp)
  rho(2,2) = cmplx(1.0_dp,0.0_dp,kind=dp)
  eps = cmplx(0.0_dp,0.0_dp,kind=dp)
  eps(1,1) = cmplx(1.0_dp,0.0_dp,kind=dp)
  eps(2,2) = cmplx(1.0_dp,0.0_dp,kind=dp)

  flav_q = [1,-1,6,21,1]
  born_q = [1,-1,6,1]
  p = 0.0_dp
  p(:,1) = [500.0_dp,0.0_dp,0.0_dp,500.0_dp]
  p(:,2) = [500.0_dp,0.0_dp,0.0_dp,-500.0_dp]
  p(:,3) = [350.0_dp,100.0_dp,0.0_dp,sqrt(350.0_dp**2-100.0_dp**2-173.0_dp**2)]
  p(:,4) = [100.0_dp,0.0_dp,80.0_dp,60.0_dp]
  p(:,5) = [250.0_dp,-40.0_dp,30.0_dp,sqrt(250.0_dp**2-40.0_dp**2-30.0_dp**2-5.0_dp**2)]
  masses = [0.0_dp,0.0_dp,173.0_dp,0.0_dp,5.0_dp]

  call cs_lc_dipole_spinrho(p,flav_q,born_q,[3,4,5],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=173.0_dp)
  call check_status(info,'massive FF quark emitter')
  dip_reference=dip
  call check_scaled_kernel(1.0e-20_dp,dip_reference,'small-scale massive FF kernel')
  call check_scaled_kernel(1.0e20_dp,dip_reference,'large-scale massive FF kernel')

  call cs_lc_dipole_spinrho(p,flav_q,born_q,[3,4,1],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=173.0_dp)
  call check_status(info,'massive FI quark emitter')

  flav_q = [1,-1,1,21,6]
  born_q = [1,-1,1,6]
  p(:,3) = [300.0_dp,100.0_dp,0.0_dp,sqrt(300.0_dp**2-100.0_dp**2)]
  masses = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,5.0_dp]
  call cs_lc_dipole_spinrho(p,flav_q,born_q,[1,4,5],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'massive IF spectator')

  flav_g = [1,-1,21,21,21]
  born_g = [1,-1,21,21]
  call cs_lc_dipole_spinrho(p,flav_g,born_g,[3,4,5],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'massive FF gluon spectator')

  ! Regression for a gen23 point that rounded an initial--final collinear
  ! invariant to zero.  It must be identified as numerically unresolved,
  ! never interpreted as a valid zero-valued dipole.
  flav_q=[-1,2,21,21,24]
  born_q=[-1,2,21,24]
  masses=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,80.419002445756163_dp]
  p(:,1)=[5236.9248837217456_dp,0.0_dp,0.0_dp,5236.9248837217456_dp]
  p(:,2)=[5496.9281402961697_dp,0.0_dp,0.0_dp,-5496.9281402961697_dp]
  p(:,3)=[260.13903026870383_dp,-1.7763568394002505e-15_dp,&
       -7.1054273576010019e-15_dp,260.13903026870383_dp]
  p(:,4)=[5496.5924210813537_dp,11.056786979018097_dp,&
       46.478195495017438_dp,-5496.3847907778609_dp]
  p(:,5)=[4977.1215726678574_dp,-11.056786979018096_dp,&
       -46.478195495017431_dp,4976.2425039347327_dp]
  dotij=dot4(p(:,1),p(:,3))
  dot_scale=cs_dot4_scale(p(:,1),p(:,3))
  call check_true(.not.cs_value_is_resolved(dotij,dot_scale),&
       'reproduced collinear invariant is unresolved')
  call cs_lc_dipole_spinrho(p,flav_q,born_q,[1,3,4],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=0.0_dp)
  call check_expected_status(info,-10,'unresolved initial-final invariant')
  call check_true(dip.eq.0.0_dp,'unresolved dipole is flagged before use')
  call check_true(dipole_status_is_numerical(info),&
       'unresolved invariant status is classified as numerical')
  call check_true(.not.dipole_status_is_numerical(-102),&
       'unsupported flavour status remains fatal')
  call check_true(dipole_status_is_numerical(-20),&
       'ill-resolved boost status is classified as numerical')
  call check_true(is_quark(1) .and. is_quark(-6),&
       'physical quark flavours are recognized')
  call check_true(.not.is_quark(0) .and. .not.is_quark(11) .and. .not.is_quark(21),&
       'non-parton flavours are not classified as quarks')
  call check_true(.not.is_q_qbar_pair(11,-11),&
       'a lepton pair is not classified as a q-qbar splitting')
  minimum_integer=-huge(0)-1
  call check_true(.not.is_quark(minimum_integer),&
       'minimum integer is not classified as a quark')
  call check_true(.not.is_q_qbar_pair(1,minimum_integer),&
       'minimum integer is rejected before q-qbar sign comparison')

  ! A nearby, representable collinear point remains evaluable.
  p(:,3)=[260.13903026870383_dp,1.0e-2_dp,0.0_dp,&
       sqrt(260.13903026870383_dp**2-1.0e-4_dp)]
  call cs_lc_dipole_spinrho(p,flav_q,born_q,[1,3,4],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'resolved nearby initial-final invariant')

  nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
  p(0,3)=nan_value
  call cs_lc_dipole_spinrho(p,flav_q,born_q,[1,3,4],1.0_dp/(4.0_dp*pi_dp),rho,eps,dip,info, &
       mass_real=masses,mass_parent=0.0_dp)
  call check_expected_status(info,-20,'non-finite dipole input')
  call check_true(dip.eq.0.0_dp,'non-finite dipole input clears output')

  write (*,'(a)') 'dipole kernel regression passed'

contains

  subroutine check_scaled_kernel(scale,reference,label)
    real(dp),intent(in) :: scale,reference
    character(len=*),intent(in) :: label
    real(dp) :: scaled_value,tolerance
    integer :: local_info

    pscaled=scale*p
    masses_scaled=scale*masses
    call cs_lc_dipole_spinrho(pscaled,flav_q,born_q,[3,4,5],1.0_dp/(4.0_dp*pi_dp),&
         rho,eps,scaled_value,local_info,mass_real=masses_scaled,mass_parent=scale*173.0_dp)
    call check_status(local_info,label)
    tolerance=2.0e-10_dp*max(1.0_dp,abs(reference))
    call check_true(abs(scaled_value*scale*scale-reference).le.tolerance,&
         trim(label)//' violates scale covariance')
  end subroutine check_scaled_kernel

  subroutine check_status(status,label)
    integer,intent(in) :: status
    character(len=*),intent(in) :: label
    if (status /= 0) then
       write (*,*) 'kernel failed:',trim(label),status
       stop 1
    endif
  end subroutine check_status

  subroutine check_expected_status(status,expected,label)
    integer,intent(in) :: status,expected
    character(len=*),intent(in) :: label
    if (status.ne.expected) then
       write (*,*) 'unexpected kernel status:',trim(label),status,expected
       stop 1
    endif
  end subroutine check_expected_status

  subroutine check_true(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if (.not.condition) then
       write (*,*) 'failed:',trim(label)
       stop 1
    endif
  end subroutine check_true

end program massive_dipole_kernel
