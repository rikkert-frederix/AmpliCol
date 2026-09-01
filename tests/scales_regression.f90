program scales_regression
  use common, only: phys_model
  use run_parameters, only: reset_run_parameters,z_mass,alphaS_MZ
  use scales, only: set_scale,alphas_Q,newton1
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  implicit none

  real(kind=8) :: value,output,nan,p(0:3,4),cmass,bmass
  integer :: info
  integer :: flavours(4)
  common /qmass/ cmass,bmass

  call reset_run_parameters()
  call phys_model%init_part()

  value=alphas_Q(z_mass,2,alphaS_MZ,info)
  call assert_true(info.eq.0 .and. ieee_is_finite(value) .and. &
       abs(value-alphaS_MZ).le.1d-12,'alpha_s(MZ) normalization')
  value=alphas_Q(2d0,2,alphaS_MZ,info)
  call assert_true(info.eq.0 .and. ieee_is_finite(value) .and. value.gt.alphaS_MZ,&
       'finite low-scale alpha_s evolution')
  value=alphas_Q(1d-100,2,alphaS_MZ,info)
  call assert_true(info.eq.-20 .and. value.eq.0d0,'Landau-region scale rejected')
  value=alphas_Q(z_mass,4,alphaS_MZ,info)
  call assert_true(info.eq.-4 .and. value.eq.0d0,'unsupported alpha_s loop order rejected')
  nan=ieee_value(0d0,ieee_quiet_nan)
  value=alphas_Q(nan,2,alphaS_MZ,info)
  call assert_true(info.eq.-20 .and. value.eq.0d0,'non-finite alpha_s scale rejected')
  value=alphas_Q(z_mass,2,huge(1d0),info)
  call assert_true(info.eq.-20 .and. value.eq.0d0,'unsafe alpha_s input rejected')

  call newton1(0d0,alphaS_MZ,output,3,5,info)
  call assert_true(info.eq.0 .and. abs(output-alphaS_MZ).le.1d-12,&
       'zero-separation Newton evolution')
  call newton1(-100d0,alphaS_MZ,output,2,5,info)
  call assert_true(info.eq.-20 .and. output.eq.0d0,'invalid Newton domain rejected')

  bmass=cmass
  value=alphas_Q(z_mass,2,alphaS_MZ,info)
  call assert_true(info.eq.-4 .and. value.eq.0d0,'invalid threshold ordering rejected')
  cmass=1.42d0
  bmass=4.7d0
  value=alphas_Q(z_mass,2,alphaS_MZ,info)
  call assert_true(info.eq.0 .and. abs(value-alphaS_MZ).le.1d-12,&
       'threshold cache recovers after invalid input')

  p=0d0
  p(:,1)=[100d0,0d0,0d0,100d0]
  p(:,2)=[100d0,0d0,0d0,-100d0]
  p(:,3)=[100d0,30d0,0d0,0d0]
  p(:,4)=[100d0,-40d0,0d0,0d0]
  flavours=[1,-1,23,23]
  call set_scale(3,4,p,flavours,value,info)
  call assert_true(info.eq.0 .and. abs(value-200d0).le.1d-12,&
       'partonic invariant-mass scale')
  call set_scale(4,4,p,flavours,value,info)
  call assert_true(info.eq.-4 .and. value.eq.0d0,'minimum-jet scale without a jet rejected')
  p(0,3)=nan
  call set_scale(1,4,p,flavours,value,info)
  call assert_true(info.eq.-20 .and. value.eq.0d0,'non-finite scale momentum rejected')

  write(*,'(a)') 'scale regression tests: PASS'

contains

  subroutine assert_true(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if (.not.condition) then
       write(*,*) 'FAIL: ',trim(label)
       stop 1
    endif
  end subroutine assert_true

end program scales_regression
