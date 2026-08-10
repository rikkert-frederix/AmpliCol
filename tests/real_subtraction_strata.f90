program test_real_subtraction_strata
  use real_subtraction_strata
  implicit none
  logical :: alpha_active(3),mapped_pass(3)
  real(kind=8) :: regular,migration,distance
  real(kind=8) :: mapped_margins(3)
  integer :: stratum

  alpha_active=[.true.,.true.,.false.]
  mapped_pass=[.true.,.true.,.false.]
  stratum=classify_real_subtraction_stratum(.true.,alpha_active,mapped_pass)
  call require(stratum.eq.real_stratum_regular,'matching measurements were not regular')

  mapped_pass=[.true.,.false.,.false.]
  stratum=classify_real_subtraction_stratum(.true.,alpha_active,mapped_pass)
  call require(stratum.eq.real_stratum_migration,'active measurement mismatch was not isolated')
  call split_real_subtraction_weight(-17d0,stratum,regular,migration)
  call require(regular.eq.0d0 .and. migration.eq.-17d0,'migration split changed the point weight')
  call require(regular+migration.eq.-17d0,'migration split does not close pointwise')

  mapped_pass=[.true.,.true.,.true.]
  stratum=classify_real_subtraction_stratum(.true.,alpha_active,mapped_pass)
  call require(stratum.eq.real_stratum_regular,'inactive mismatch affected the classification')
  call split_real_subtraction_weight(23d0,stratum,regular,migration)
  call require(regular.eq.23d0 .and. migration.eq.0d0,'regular split changed the point weight')
  call require(regular+migration.eq.23d0,'regular split does not close pointwise')

  mapped_pass=[.true.,.false.,.false.]
  mapped_margins=[0.20d0,2.0d0,0.01d0]
  distance=migration_pt_distance(-0.40d0,mapped_margins,.true.,alpha_active,mapped_pass)
  call require(abs(distance-0.40d0).lt.1d-14,'wrong migration threshold distance')

  mapped_margins=[0.20d0,-2.0d0,0.01d0]
  distance=migration_pt_distance(-0.40d0,mapped_margins,.true.,alpha_active,mapped_pass)
  call require(distance.lt.0d0,'non-pT migration received a pT-threshold distance')

  mapped_pass=[.true.,.true.,.false.]
  distance=migration_pt_distance(-0.40d0,mapped_margins,.true.,alpha_active,mapped_pass)
  call require(distance.lt.0d0,'regular point received a migration threshold distance')

  write(*,'(a)') 'real subtraction strata: PASS'

contains

  subroutine require(condition,message)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: message
    if (.not.condition) then
       write(*,'(a)') 'FAIL: '//trim(message)
       stop 1
    endif
  end subroutine require

end program test_real_subtraction_strata
