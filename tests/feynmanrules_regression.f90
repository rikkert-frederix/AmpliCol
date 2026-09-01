program feynmanrules_regression
  use FeynmanRules
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none

  real(kind=8),parameter :: tolerance=2d-14
  real(kind=8) :: p(0:3),sqrt_mass,nan_value
  real(kind=8) :: momenta(0:3,2)
  real(kind=8) :: real_vector(4)
  complex(kind=8) :: spinor(4),expected(4),vector(4),weyl(2),scalar(1),reference(4)
  complex(kind=8) :: contraction

  p=[4d0,0d0,0d0,0d0]
  sqrt_mass=2d0
  call ext_fermion_outflow(p,1,1,spinor,4d0)
  expected=[(0d0,0d0),(-2d0,0d0),(0d0,0d0),(-2d0,0d0)]
  call assert_complex_close(spinor,expected,'positive-energy outflow rest spinor, helicity +')
  call ext_fermion_outflow(p,-1,1,spinor,4d0)
  expected=[(2d0,0d0),(0d0,0d0),(2d0,0d0),(0d0,0d0)]
  call assert_complex_close(spinor,expected,'positive-energy outflow rest spinor, helicity -')
  call ext_fermion_inflow(p,1,1,spinor,4d0)
  expected=[(0d0,0d0),(-2d0,0d0),(0d0,0d0),(2d0,0d0)]
  call assert_complex_close(spinor,expected,'positive-energy inflow rest spinor, helicity +')
  call ext_fermion_inflow(p,-1,1,spinor,4d0)
  expected=[(2d0,0d0),(0d0,0d0),(-2d0,0d0),(0d0,0d0)]
  call assert_complex_close(spinor,expected,'positive-energy inflow rest spinor, helicity -')

  p(0)=-4d0
  call ext_fermion_outflow(p,1,1,spinor,4d0)
  expected=[(2d0,0d0),(0d0,0d0),(-2d0,0d0),(0d0,0d0)]
  call assert_complex_close(spinor,expected,'negative-energy outflow rest spinor')
  call ext_fermion_inflow(p,-1,1,spinor,4d0)
  expected=[(0d0,0d0),(2d0,0d0),(0d0,0d0),(2d0,0d0)]
  call assert_complex_close(spinor,expected,'negative-energy inflow rest spinor')

  ! A nonzero mass must never be reclassified as massless based on its units.
  p=[1d-12,0d0,0d0,0d0]
  sqrt_mass=1d-6
  call ext_fermion_outflow(p,-1,1,spinor,1d-12)
  expected=[cmplx(sqrt_mass,0d0,kind=8),(0d0,0d0),&
       cmplx(sqrt_mass,0d0,kind=8),(0d0,0d0)]
  call assert_complex_close(spinor,expected,'tiny nonzero mass uses massive spinor')

  call reset_feynman_numerical_status()
  p=[10d0,3d0,4d0,sqrt(75d0)]
  call ext_massless_vector_cmplx(p,1,1,vector)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(vector),&
       'regular external massless vector')

  call reset_feynman_numerical_status()
  p=[2d0,0d0,0d0,0d0]
  call ext_massive_vector(p,0,1,vector,2d0)
  expected=[(0d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)]
  call assert_true(feynman_numerical_status_ok(),'massive-vector rest status')
  call assert_complex_close(vector,expected,'massive-vector rest wavefunction')

  call reset_feynman_numerical_status()
  p=0d0
  call ext_massless_vector_cmplx(p,1,1,vector)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'zero-energy external massless vector rejected')
  call assert_true(all(vector.eq.(0d0,0d0)),&
       'zero-energy external massless-vector output cleared')

  nan_value=ieee_value(0d0,ieee_quiet_nan)
  call reset_feynman_numerical_status()
  p=[nan_value,0d0,0d0,0d0]
  call ext_fermion_outflow(p,1,1,spinor,0d0)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite external fermion momentum rejected')
  call assert_true(all(spinor.eq.(0d0,0d0)),&
       'non-finite external fermion output cleared')

  call reset_feynman_numerical_status()
  call ext_fermion_inflow_weyl(p,-1,1,weyl,-1)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite external Weyl momentum rejected')
  call assert_true(all(weyl.eq.(0d0,0d0)),&
       'non-finite external Weyl output cleared')

  call reset_feynman_numerical_status()
  call ext_scalar(p,1,scalar)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite external scalar momentum rejected')
  call assert_true(all(scalar.eq.(0d0,0d0)),&
       'non-finite external scalar output cleared')

  call reset_feynman_numerical_status()
  p=[2d0,0d0,0d0,0d0]
  vector=(1d0,0d0)
  call MasslessVectorPropagator(vector,p)
  expected=(0d0,-0.25d0)
  call assert_true(feynman_numerical_status_ok(),'regular massless-vector propagator status')
  call assert_complex_close(vector,expected,'regular massless-vector propagator value')

  call reset_feynman_numerical_status()
  real_vector=1d0
  call MasslessVectorPropagator_real(real_vector,p)
  call assert_true(feynman_numerical_status_ok(),'regular real massless-vector propagator status')
  call assert_true(maxval(abs(real_vector-0.25d0)).le.tolerance,&
       'regular real massless-vector propagator value')

  call reset_feynman_numerical_status()
  p=[1d0,1d0,0d0,0d0]
  vector=(1d0,0d0)
  call MasslessVectorPropagator(vector,p)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'singular massless-vector propagator rejected')
  call assert_true(all(vector.eq.(0d0,0d0)),'singular massless-vector output cleared')

  call reset_feynman_numerical_status()
  p=[3d0,0d0,0d0,0d0]
  vector=(1d0,0d0)
  call MassiveVectorPropagator(vector,p,2d0,0.5d0)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(vector),&
       'regular massive-vector propagator')

  call reset_feynman_numerical_status()
  vector=(1d0,0d0)
  call FermionPropagator(vector,p,2d0,0.5d0)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(vector),&
       'regular fermion propagator')
  call reset_feynman_numerical_status()
  vector=(1d0,0d0)
  call AntifermionPropagator(vector,p,2d0,0.5d0)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(vector),&
       'regular antifermion propagator')

  call reset_feynman_numerical_status()
  weyl=(1d0,0d0)
  call FermionPropagator_weyl(weyl,p,1)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(weyl),&
       'regular Weyl-fermion propagator')
  call reset_feynman_numerical_status()
  weyl=(1d0,0d0)
  call AntifermionPropagator_weyl(weyl,p,-1)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(weyl),&
       'regular Weyl-antifermion propagator')

  call reset_feynman_numerical_status()
  scalar=(1d0,0d0)
  call ScalarPropagator(scalar,p,2d0,0.5d0)
  call assert_true(feynman_numerical_status_ok() .and. complex_vector_is_finite(scalar),&
       'regular scalar propagator')

  call reset_feynman_numerical_status()
  p=[2d0,0d0,0d0,0d0]
  scalar=(1d0,0d0)
  call ScalarPropagator(scalar,p,2d0,0d0)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'zero-width massive pole rejected')
  call assert_true(all(scalar.eq.(0d0,0d0)),'massive-pole scalar output cleared')

  call reset_feynman_numerical_status()
  p=[nan_value,0d0,0d0,0d0]
  vector=(1d0,0d0)
  call MasslessVectorPropagator(vector,p)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite propagator momentum rejected')

  call reset_feynman_numerical_status()
  p=[2d0,0d0,0d0,0d0]
  vector=cmplx(huge(1d0),0d0,kind=8)
  call MasslessVectorPropagator(vector,p)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'unsafe propagator current rejected before multiplication')

  call reset_feynman_numerical_status()
  momenta=0d0
  momenta(0,1)=nan_value
  call validate_feynman_momenta(momenta)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite evaluator momentum array rejected')

  call reset_feynman_numerical_status()
  spinor=(1d0,0d0)
  vector=(1d0,0d0)
  vector(1)=cmplx(nan_value,0d0,kind=8)
  contraction=ContractFermionCurrents(spinor,0,vector,0)
  call assert_true(.not.feynman_numerical_status_ok(),&
       'non-finite terminal contraction current rejected')
  call assert_true(contraction.eq.(0d0,0d0),&
       'invalid current contraction output cleared')

  ! Real-gluon kernels must be exact restrictions of their complex counterparts.
  spinor=[cmplx(1d0,-2d0,kind=8),cmplx(-3d0,4d0,kind=8),&
       cmplx(5d0,6d0,kind=8),cmplx(-7d0,-8d0,kind=8)]
  real_vector=[0.5d0,-1.25d0,2.5d0,-0.75d0]
  vector=cmplx(real_vector,0d0,kind=8)
  call AntiquarkColourFlowVectorToAntiquark(spinor,vector,reference)
  call AntiquarkColourFlowVectorToAntiquark_real(spinor,real_vector,expected)
  call assert_complex_close(expected,reference,&
       'real antiquark-vector kernel agrees with complex kernel')
  call ColourFlowVectorAntiquarkToAntiquark(vector,spinor,reference)
  call ColourFlowVectorAntiquarkToAntiquark_real(real_vector,spinor,expected)
  call assert_complex_close(expected,reference,&
       'real vector-antiquark kernel agrees with complex kernel')

  write(*,'(a)') 'Feynman-rule numerical regression tests: PASS'

contains

  subroutine assert_true(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if (.not.condition) then
       write(*,*) 'FAIL: ',trim(label)
       stop 1
    endif
  end subroutine assert_true

  subroutine assert_complex_close(actual,reference,label)
    complex(kind=8),intent(in) :: actual(:),reference(:)
    character(len=*),intent(in) :: label
    real(kind=8) :: scale
    call assert_true(complex_vector_is_finite(actual),trim(label)//' is finite')
    scale=max(1d0,maxval(abs(reference)))
    call assert_true(maxval(abs(actual-reference)).le.tolerance*scale,label)
  end subroutine assert_complex_close

  logical function complex_vector_is_finite(values)
    complex(kind=8),intent(in) :: values(:)
    complex_vector_is_finite=all(ieee_is_finite(real(values,kind=8))) .and. &
         all(ieee_is_finite(aimag(values)))
  end function complex_vector_is_finite

end program feynmanrules_regression
