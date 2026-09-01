program massive_dipole_mapping
  use cs_dipole_mappings, only: dp, cs_map, cs_dipole_cut_variable, &
       cs_born_pushback_weight, dot4
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  implicit none

  real(dp) :: p(0:3,5), pt(0:3,4), pdegenerate(0:3,5), masses(5), mp
  integer :: info
  real(dp) :: cut_variable,expected_cut,raw_y,q2,mui2,muk,muk2,yplus,pushback_weight,nan_value

  p = 0.0_dp
  p(:,1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
  p(:,2) = [500.0_dp, 0.0_dp, 0.0_dp,-500.0_dp]
  p(:,3) = [350.0_dp, 100.0_dp, 0.0_dp, sqrt(350.0_dp**2-100.0_dp**2-173.0_dp**2)]
  p(:,4) = [100.0_dp, 0.0_dp, 80.0_dp, 60.0_dp]
  p(:,5) = [250.0_dp, -40.0_dp, 30.0_dp, sqrt(250.0_dp**2-40.0_dp**2-30.0_dp**2-5.0_dp**2)]
  masses = [0.0_dp, 0.0_dp, 173.0_dp, 0.0_dp, 5.0_dp]

  ! FF: massive emitter and massive spectator.
  mp = masses(3)
  call cs_map(p,[3,4,5],pt,info,mass_real=masses,mass_parent=mp)
  call check_status(info,'FF')
  call cs_born_pushback_weight(p,pt,[3,4,5],masses,mp,pushback_weight,info)
  call check_positive_weight(info,pushback_weight,'FF pushback')
  call check_mass(pt(:,3),mp,'FF parent')
  call check_mass(pt(:,4),masses(5),'FF spectator')
  call check_vector(pt(:,3)+pt(:,4)-p(:,3)-p(:,4)-p(:,5),'FF momentum')

  ! FI: massive final-state emitter and massless initial spectator.
  call cs_map(p,[3,4,1],pt,info,mass_real=masses,mass_parent=mp)
  call check_status(info,'FI')
  call cs_born_pushback_weight(p,pt,[3,4,1],masses,mp,pushback_weight,info)
  call check_positive_weight(info,pushback_weight,'FI pushback')
  call check_mass(pt(:,3),mp,'FI parent')
  call check_mass(pt(:,1),0.0_dp,'FI initial spectator')

  ! IF: massless initial-state emitter and massive final-state spectator.
  call cs_map(p,[1,4,5],pt,info,mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'IF')
  call cs_born_pushback_weight(p,pt,[1,4,5],masses,0.0_dp,pushback_weight,info)
  call check_positive_weight(info,pushback_weight,'IF pushback')
  call check_mass(pt(:,4),masses(5),'IF spectator')
  call check_mass(pt(:,1),0.0_dp,'IF initial parent')

  call cs_dipole_cut_variable(p,[3,4,5],masses,mp,cut_variable,info)
  call check_cut(info,cut_variable,'FF')
  raw_y=dot4(p(:,3),p(:,4))/(dot4(p(:,3),p(:,4))+&
       dot4(p(:,3),p(:,5))+dot4(p(:,4),p(:,5)))
  q2=mp*mp+masses(5)*masses(5)+2.0_dp*&
       (dot4(p(:,3),p(:,4))+dot4(p(:,3),p(:,5))+dot4(p(:,4),p(:,5)))
  mui2=mp*mp/q2
  muk2=masses(5)*masses(5)/q2
  muk=sqrt(muk2)
  yplus=1.0_dp-2.0_dp*muk*(1.0_dp-muk)/(1.0_dp-mui2-muk2)
  call check_close(cut_variable,raw_y/yplus,'FF uses y/yplus')
  call cs_dipole_cut_variable(p,[3,4,1],masses,mp,cut_variable,info)
  call check_cut(info,cut_variable,'FI')
  call cs_dipole_cut_variable(p,[1,4,5],masses,0.0_dp,cut_variable,info)
  call check_cut(info,cut_variable,'IF')
  call cs_dipole_cut_variable(p,[1,4,2],masses,0.0_dp,cut_variable,info)
  call check_cut(info,cut_variable,'II')
  expected_cut=dot4(p(:,1),p(:,4))/dot4(p(:,1),p(:,2))
  call check_close(cut_variable,expected_cut,'II uses v_j')

  call cs_map(p,[1,4,2],pt,info,mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'II')
  call cs_born_pushback_weight(p,pt,[1,4,2],masses,0.0_dp,pushback_weight,info)
  call check_positive_weight(info,pushback_weight,'II pushback')

  call cs_born_pushback_weight(p,pt,[0,4,2],masses,0.0_dp,pushback_weight,info)
  call check_expected_status(info,-2,'pushback index bounds')
  call cs_born_pushback_weight(p,pt,[1,4,1],masses,0.0_dp,pushback_weight,info)
  call check_expected_status(info,-3,'pushback distinct indices')
  pdegenerate=p
  pdegenerate(:,2)=pdegenerate(:,1)
  call cs_born_pushback_weight(pdegenerate,pt,[1,4,2],masses,0.0_dp,pushback_weight,info)
  call check_expected_status(info,-13,'pushback zero II invariant')

  ! The maps are homogeneous in the momentum unit.  Exercise scales far
  ! below the former absolute 1e-30 invariant floor and far above the
  ! ordinary collider scale.
  call check_scaled_mapping(1.0e-20_dp,[3,4,5],mp,'small-scale FF')
  call check_scaled_mapping(1.0e-20_dp,[3,4,1],mp,'small-scale FI')
  call check_scaled_mapping(1.0e-20_dp,[1,4,5],0.0_dp,'small-scale IF')
  call check_scaled_mapping(1.0e-20_dp,[1,4,2],0.0_dp,'small-scale II')
  call check_scaled_mapping(1.0e20_dp,[3,4,5],mp,'large-scale FF')
  call check_scaled_mapping(1.0e20_dp,[3,4,1],mp,'large-scale FI')
  call check_scaled_mapping(1.0e20_dp,[1,4,5],0.0_dp,'large-scale IF')
  call check_scaled_mapping(1.0e20_dp,[1,4,2],0.0_dp,'large-scale II')

  ! Non-finite and non-physical variables must be rejected without exposing
  ! NaNs or a partially filled mapped point.
  nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
  pdegenerate=p
  pdegenerate(0,3)=nan_value
  call cs_map(pdegenerate,[3,4,5],pt,info,mass_real=masses,mass_parent=mp)
  call check_expected_status(info,-20,'mapping rejects NaN')
  call check_zero_matrix(pt,'failed mapping clears output')
  call cs_dipole_cut_variable(pdegenerate,[3,4,5],masses,mp,cut_variable,info)
  call check_expected_status(info,-20,'alpha variable rejects NaN')
  call cs_born_pushback_weight(pdegenerate,pt,[3,4,5],masses,mp,pushback_weight,info)
  call check_expected_status(info,-20,'pushback rejects NaN')

  pdegenerate=p
  pdegenerate(:,4)=[2000.0_dp,0.0_dp,0.0_dp,-2000.0_dp]
  call cs_dipole_cut_variable(pdegenerate,[1,4,2],masses,0.0_dp,cut_variable,info)
  call check_expected_status(info,-20,'alpha variable rejects value above one')
  call cs_map(pdegenerate,[1,4,2],pt,info,mass_real=masses,mass_parent=0.0_dp)
  call check_expected_status(info,-13,'mapping rejects negative x')

  write (*,'(a)') 'massive dipole mapping regression passed'

contains

  subroutine check_scaled_mapping(scale,indices,parent_mass,label)
    real(dp),intent(in) :: scale,parent_mass
    integer,intent(in) :: indices(3)
    character(len=*),intent(in) :: label
    real(dp) :: pscaled(0:3,5),masses_scaled(5)
    real(dp) :: mapped_reference(0:3,4),mapped_scaled(0:3,4)
    real(dp) :: cut_reference,cut_scaled,weight_reference,weight_scaled
    real(dp) :: difference,normalisation
    integer :: local_info

    pscaled=scale*p
    masses_scaled=scale*masses
    call cs_map(p,indices,mapped_reference,local_info,mass_real=masses,mass_parent=parent_mass)
    call check_status(local_info,trim(label)//' reference map')
    call cs_map(pscaled,indices,mapped_scaled,local_info,mass_real=masses_scaled,&
         mass_parent=scale*parent_mass)
    call check_status(local_info,trim(label)//' scaled map')
    difference=maxval(abs(mapped_scaled/scale-mapped_reference))
    normalisation=max(1.0_dp,maxval(abs(mapped_reference)))
    if (difference.gt.2.0e-10_dp*normalisation) then
       write (*,*) 'scaled mapping mismatch:',trim(label),difference,normalisation
       stop 1
    endif

    call cs_dipole_cut_variable(p,indices,masses,parent_mass,cut_reference,local_info)
    call check_status(local_info,trim(label)//' reference alpha variable')
    call cs_dipole_cut_variable(pscaled,indices,masses_scaled,scale*parent_mass,cut_scaled,local_info)
    call check_status(local_info,trim(label)//' scaled alpha variable')
    call check_close(cut_scaled,cut_reference,trim(label)//' alpha scale invariance')

    call cs_born_pushback_weight(p,mapped_reference,indices,masses,parent_mass,weight_reference,local_info)
    call check_positive_weight(local_info,weight_reference,trim(label)//' reference pushback')
    call cs_born_pushback_weight(pscaled,mapped_scaled,indices,masses_scaled,scale*parent_mass,&
         weight_scaled,local_info)
    call check_positive_weight(local_info,weight_scaled,trim(label)//' scaled pushback')
    call check_close(weight_scaled*scale*scale,weight_reference,trim(label)//' pushback scaling')
  end subroutine check_scaled_mapping

  subroutine check_status(status,label)
    integer,intent(in) :: status
    character(len=*),intent(in) :: label
    if (status /= 0) then
       write (*,*) 'mapping failed:',trim(label),status
       stop 1
    endif
  end subroutine check_status

  subroutine check_mass(q,mass,label)
    real(dp),intent(in) :: q(0:3),mass
    character(len=*),intent(in) :: label
    if (abs(dot4(q,q)-mass*mass) > 1.0e-8_dp*max(1.0_dp,mass*mass)) then
       write (*,*) 'mass-shell check failed:',trim(label),dot4(q,q),mass*mass
       stop 1
    endif
  end subroutine check_mass

  subroutine check_vector(q,label)
    real(dp),intent(in) :: q(0:3)
    character(len=*),intent(in) :: label
    if (maxval(abs(q)) > 1.0e-8_dp) then
       write (*,*) 'momentum check failed:',trim(label),q
       stop 1
    endif
  end subroutine check_vector

  subroutine check_zero_matrix(q,label)
    real(dp),intent(in) :: q(0:,:)
    character(len=*),intent(in) :: label
    if (.not.all(ieee_is_finite(q))) then
       write (*,*) 'zero-matrix check is non-finite:',trim(label)
       stop 1
    endif
    if (maxval(abs(q)).gt.0.0_dp) then
       write (*,*) 'zero-matrix check failed:',trim(label),maxval(abs(q))
       stop 1
    endif
  end subroutine check_zero_matrix

  subroutine check_cut(status,cut,label)
    integer,intent(in) :: status
    real(dp),intent(in) :: cut
    character(len=*),intent(in) :: label
    if (status /= 0 .or. cut < 0.0_dp .or. cut > 1.0_dp) then
       write (*,*) 'alpha-variable check failed:',trim(label),status,cut
       stop 1
    endif
  end subroutine check_cut

  subroutine check_close(value,expected,label)
    real(dp),intent(in) :: value,expected
    character(len=*),intent(in) :: label
    if (abs(value-expected) > 1.0e-13_dp*max(1.0_dp,abs(expected))) then
       write (*,*) 'value check failed:',trim(label),value,expected
       stop 1
    endif
  end subroutine check_close

  subroutine check_positive_weight(status,value,label)
    integer,intent(in) :: status
    real(dp),intent(in) :: value
    character(len=*),intent(in) :: label
    if (status.ne.0 .or. value.le.0.0_dp) then
       write (*,*) 'pushback-weight check failed:',trim(label),status,value
       stop 1
    endif
  end subroutine check_positive_weight

  subroutine check_expected_status(status,expected,label)
    integer,intent(in) :: status,expected
    character(len=*),intent(in) :: label
    if (status.ne.expected) then
       write (*,*) 'status check failed:',trim(label),status,expected
       stop 1
    endif
  end subroutine check_expected_status

end program massive_dipole_mapping
