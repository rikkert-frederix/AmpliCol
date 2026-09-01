program phase_space_module_safety
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use phase_space_module
  implicit none

  integer, parameter :: nout = 4
  real(dp), parameter :: scale_low = 1.0e-20_dp
  real(dp), parameter :: scale_high = 1.0e20_dp
  real(dp) :: mass(nout), mass_low(nout), mass_high(nout)
  real(dp) :: pbase(0:3,nout+2), plow(0:3,nout+2), phigh(0:3,nout+2)
  real(dp) :: pnext(0:3,nout+2), qbase(0:3,nout+2)
  real(dp) :: qlow(0:3,nout+2), qhigh(0:3,nout+2), bad(0:3,nout+2)
  real(dp) :: mass2(2), p2(0:3,4), q2(0:3,4), nan_value
  real(dp) :: res_p, res_m, res_in, original_s, deformed_s
  integer :: status

  mass = (/ 0.0_dp, 0.0_dp, 5.0_dp, 10.0_dp /)
  mass_low = scale_low*mass
  mass_high = scale_high*mass

  call set_test_seed()
  call generate_cm_scattering_phase_space(nout, mass, 1000.0_dp, pbase, status)
  call require(status == PS_OK, 'base phase-space generation failed')
  call require_valid(nout, mass, pbase, 'base phase-space point')

  call set_test_seed()
  call generate_cm_scattering_phase_space(nout, mass_low, scale_low*1000.0_dp, &
       plow, status)
  call require(status == PS_OK, 'low-scale phase-space generation failed')
  call require_valid(nout, mass_low, plow, 'low-scale phase-space point')
  call require(relative_difference(plow/scale_low, pbase) < 2.0e-10_dp, &
       'phase-space generation is not covariant at low scale')

  call set_test_seed()
  call generate_cm_scattering_phase_space(nout, mass_high, scale_high*1000.0_dp, &
       phigh, status)
  call require(status == PS_OK, 'high-scale phase-space generation failed')
  call require_valid(nout, mass_high, phigh, 'high-scale phase-space point')
  call require(relative_difference(phigh/scale_high, pbase) < 2.0e-10_dp, &
       'phase-space generation is not covariant at high scale')

  call set_test_seed()
  call generate_cm_scattering_phase_space(nout, mass, 1000.0_dp, qbase, status)
  call require(status == PS_OK, 'first RNG progression point failed')
  call generate_cm_scattering_phase_space(nout, mass, 1000.0_dp, pnext, status)
  call require(status == PS_OK, 'second RNG progression point failed')
  call require(relative_difference(pnext, qbase) > 1.0e-8_dp, &
       'phase-space generator resets the global RNG on each call')

  call soft_deform_event(nout, mass, pbase, COL_FIRST_OUT, 1.0e-4_dp, qbase, status)
  call require(status == PS_OK, 'base soft deformation failed')
  call require_valid(nout, mass, qbase, 'base soft deformation')
  call soft_deform_event(nout, mass_low, plow, COL_FIRST_OUT, 1.0e-4_dp, qlow, status)
  call require(status == PS_OK, 'low-scale soft deformation failed')
  call soft_deform_event(nout, mass_high, phigh, COL_FIRST_OUT, 1.0e-4_dp, qhigh, status)
  call require(status == PS_OK, 'high-scale soft deformation failed')
  call require(relative_difference(qlow/scale_low, qbase) < 2.0e-8_dp, &
       'soft deformation is not covariant at low scale')
  call require(relative_difference(qhigh/scale_high, qbase) < 2.0e-8_dp, &
       'soft deformation is not covariant at high scale')

  call collinear_deform_event(nout, mass, pbase, COL_FIRST_OUT, &
       COL_FIRST_OUT+1, 1.0e-4_dp, qbase, status)
  call require(status == PS_OK, 'base final-final collinear deformation failed')
  call require_valid(nout, mass, qbase, 'base final-final collinear deformation')
  call collinear_deform_event(nout, mass_low, plow, COL_FIRST_OUT, &
       COL_FIRST_OUT+1, 1.0e-4_dp, qlow, status)
  call require(status == PS_OK, 'low-scale final-final collinear deformation failed')
  call collinear_deform_event(nout, mass_high, phigh, COL_FIRST_OUT, &
       COL_FIRST_OUT+1, 1.0e-4_dp, qhigh, status)
  call require(status == PS_OK, 'high-scale final-final collinear deformation failed')
  call require(relative_difference(qlow/scale_low, qbase) < 2.0e-8_dp, &
       'final-final deformation is not covariant at low scale')
  call require(relative_difference(qhigh/scale_high, qbase) < 2.0e-8_dp, &
       'final-final deformation is not covariant at high scale')

  call collinear_deform_event(nout, mass, pbase, COL_IN1, COL_FIRST_OUT, &
       1.0e-4_dp, qbase, status)
  call require(status == PS_OK, 'base initial-final collinear deformation failed')
  call require_valid(nout, mass, qbase, 'base initial-final collinear deformation')
  call collinear_deform_event(nout, mass_low, plow, COL_IN1, COL_FIRST_OUT, &
       1.0e-4_dp, qlow, status)
  call require(status == PS_OK, 'low-scale initial-final collinear deformation failed')
  call collinear_deform_event(nout, mass_high, phigh, COL_IN1, COL_FIRST_OUT, &
       1.0e-4_dp, qhigh, status)
  call require(status == PS_OK, 'high-scale initial-final collinear deformation failed')
  call require(relative_difference(qlow/scale_low, qbase) < 2.0e-8_dp, &
       'initial-final deformation is not covariant at low scale')
  call require(relative_difference(qhigh/scale_high, qbase) < 2.0e-8_dp, &
       'initial-final deformation is not covariant at high scale')

  mass2 = 0.0_dp
  call set_test_seed()
  call generate_cm_scattering_phase_space(2, mass2, 1000.0_dp, p2, status)
  call require(status == PS_OK, 'massless two-body generation failed')
  original_s = msq4(p2(:,1) + p2(:,2))
  call soft_deform_event(2, mass2, p2, COL_FIRST_OUT, 0.1_dp, q2, status)
  call require(status == PS_OK, 'massless two-body variable-shat soft deformation failed')
  call require_valid(2, mass2, q2, 'massless two-body soft deformation')
  deformed_s = msq4(q2(:,1) + q2(:,2))
  call require(deformed_s > 0.0_dp .and. deformed_s < original_s, &
       'two-body soft fallback did not lower the incoming invariant mass')

  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  bad = pbase
  bad(0,3) = nan_value
  call soft_deform_event(nout, mass, bad, COL_FIRST_OUT, 0.5_dp, qbase, status)
  call require(status == PS_BAD_INPUT .and. all(qbase == 0.0_dp), &
       'soft deformation did not reject NaN input cleanly')
  call collinear_deform_event(nout, mass, pbase, COL_IN1, COL_FIRST_OUT, &
       nan_value, qbase, status)
  call require(status == PS_BAD_INPUT .and. all(qbase == 0.0_dp), &
       'collinear deformation did not reject a NaN parameter cleanly')
  call generate_cm_scattering_phase_space(nout, mass, nan_value, qbase, status)
  call require(status == PS_BAD_INPUT .and. all(qbase == 0.0_dp), &
       'phase-space generation did not reject NaN energy cleanly')

  bad = plow
  bad(0,1) = 1.1_dp*bad(0,1)
  call soft_deform_event(nout, mass_low, bad, COL_FIRST_OUT, 0.5_dp, qlow, status)
  call require(status == PS_BAD_INPUT .and. all(qlow == 0.0_dp), &
       'low-scale off-shell incoming momentum was accepted')

  print *, 'phase_space_module safety regression: PASS'

contains

  subroutine set_test_seed()
    integer :: i, nseed
    integer, allocatable :: seed(:)

    call random_seed(size=nseed)
    allocate(seed(nseed))
    do i = 1, nseed
       seed(i) = 104729 + 7919*i
    end do
    call random_seed(put=seed)
  end subroutine set_test_seed

  subroutine require_valid(n, masses, p, label)
    integer, intent(in) :: n
    real(dp), intent(in) :: masses(n), p(0:3,n+2)
    character(len=*), intent(in) :: label
    real(dp) :: local_scale

    call scattering_event_residuals(n, masses, p, res_p, res_m, res_in)
    local_scale = max(maxval(abs(p)), maxval(masses), tiny(1.0_dp))
    call require(all(p(0,:) >= 0.0_dp), trim(label)//' has negative energy')
    call require(res_p <= 2.0e-8_dp*local_scale*real(max(1,n),dp), &
         trim(label)//' violates momentum conservation')
    call require(res_m <= 2.0e-8_dp*local_scale*local_scale, &
         trim(label)//' has an off-shell final momentum')
    call require(res_in <= 2.0e-8_dp*local_scale*local_scale, &
         trim(label)//' has an off-shell incoming momentum')
  end subroutine require_valid

  real(dp) function relative_difference(a, b)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: scale

    scale = max(maxval(abs(a)), maxval(abs(b)), tiny(1.0_dp))
    relative_difference = maxval(abs(a-b))/scale
  end function relative_difference

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not.condition) then
       write (*,'(a)') 'FAIL: '//trim(message)
       error stop 1
    end if
  end subroutine require

end program phase_space_module_safety
