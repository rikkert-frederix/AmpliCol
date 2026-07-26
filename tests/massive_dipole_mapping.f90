program massive_dipole_mapping
  use cs_dipole_mappings, only: dp, cs_map, cs_dipole_cut_variable, dot4
  implicit none

  real(dp) :: p(0:3,5), pt(0:3,4), masses(5), mp
  integer :: info
  real(dp) :: cut_variable,expected_cut,raw_y,q2,mui2,muk,muk2,yplus

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
  call check_mass(pt(:,3),mp,'FF parent')
  call check_mass(pt(:,4),masses(5),'FF spectator')
  call check_vector(pt(:,3)+pt(:,4)-p(:,3)-p(:,4)-p(:,5),'FF momentum')

  ! FI: massive final-state emitter and massless initial spectator.
  call cs_map(p,[3,4,1],pt,info,mass_real=masses,mass_parent=mp)
  call check_status(info,'FI')
  call check_mass(pt(:,3),mp,'FI parent')
  call check_mass(pt(:,1),0.0_dp,'FI initial spectator')

  ! IF: massless initial-state emitter and massive final-state spectator.
  call cs_map(p,[1,4,5],pt,info,mass_real=masses,mass_parent=0.0_dp)
  call check_status(info,'IF')
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

  write (*,'(a)') 'massive dipole mapping regression passed'

contains

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

end program massive_dipole_mapping
