program fd_gauge_components
  use FeynmanRules
  implicit none
  integer,parameter :: dp=kind(1d0)
  real(kind=dp),parameter :: tol=2d-12
  real(kind=dp),dimension(0:3) :: p,p_on
  complex(kind=dp),dimension(4) :: source,legacy,mapped,unitary_external
  complex(kind=dp),dimension(5) :: fd_current,fd_external
  complex(kind=dp),dimension(4) :: v1,v2,v3,direct,factorised,legacy_factorised,tmp_vector
  complex(kind=dp),dimension(6) :: tensor12,tensor23
  complex(kind=dp),dimension(16) :: auxiliary
  real(kind=dp) :: mass,width

  mass=80.419002445756163d0
  width=2.0476d0
  p=[250d0,30d0,-45d0,70d0]
  source=[cmplx(1.2d0,-0.4d0,dp),cmplx(-0.7d0,0.2d0,dp),&
       cmplx(0.3d0,1.1d0,dp),cmplx(-1.4d0,-0.6d0,dp)]

  legacy=source
  call GluonPropagator_mass(legacy,p,mass,width)
  call fd_lift_massive_current(source,p,mass,fd_current)
  call GluonPropagator_mass_fd(fd_current,p,mass,width)
  call fd_massive_to_unitary(fd_current,p,mass,mapped)
  call assert_close('massive propagator factorisation',mapped,legacy,tol)

  p_on(1:3)=[30d0,40d0,25d0]
  p_on(0)=sqrt(mass**2+sum(p_on(1:3)**2))
  call ext_gluon_mass(p_on,0,1,unitary_external,mass)
  call ext_gluon_mass_fd(p_on,0,1,fd_external,mass)
  call fd_massive_to_unitary(fd_external,p_on,mass,mapped)
  call assert_close('longitudinal external wavefunction',mapped,unitary_external,tol)
  if (abs(fd_external(5)+cmplx(0d0,1d0,dp)).gt.tol) then
     write (*,*) 'FD component test failed: longitudinal Goldstone component',fd_external(5)
     stop 1
  endif

  call ext_gluon_mass(p_on,1,1,unitary_external,mass)
  call ext_gluon_mass_fd(p_on,1,1,fd_external,mass)
  call fd_massive_to_unitary(fd_external,p_on,mass,mapped)
  call assert_close('transverse external wavefunction',mapped,unitary_external,tol)
  if (abs(fd_external(5)).gt.tol) then
     write (*,*) 'FD component test failed: transverse Goldstone component',fd_external(5)
     stop 1
  endif

  v1=[cmplx(1d0,2d0,dp),cmplx(2d0,-1d0,dp),&
       cmplx(3d0,0.5d0,dp),cmplx(-1d0,1d0,dp)]
  v2=[cmplx(2d0,-0.5d0,dp),cmplx(-1d0,2d0,dp),&
       cmplx(0.5d0,-1d0,dp),cmplx(3d0,1d0,dp)]
  v3=[cmplx(-1d0,1d0,dp),cmplx(1d0,0.2d0,dp),&
       cmplx(2d0,-2d0,dp),cmplx(0.3d0,0.7d0,dp)]
  call FourGluon(v1,v2,v3,direct)
  call TwoGluonToTensor(v1,v2,tensor12)
  call TensorGluontoGluon(tensor12,v3,legacy_factorised)
  call TwoGluonToTensor(v2,v3,tensor23)
  call GluonTensortoGluon(v1,tensor23,tmp_vector)
  legacy_factorised=legacy_factorised+tmp_vector
  call assert_close('legacy six-component four-vector decomposition',legacy_factorised,direct,tol)
  call FDTwoVectorToAux(v1,v2,auxiliary)
  call FDAuxVectorToVector(auxiliary,v3,factorised)
  call assert_close('16-component four-vector decomposition',factorised,direct,tol)

  write (*,*) 'FD-gauge component tests passed'

contains

  subroutine assert_close(label,actual,expected,threshold)
    implicit none
    character(len=*),intent(in) :: label
    complex(kind=dp),dimension(:),intent(in) :: actual,expected
    real(kind=dp),intent(in) :: threshold
    real(kind=dp) :: scale,error
    scale=max(1d0,maxval(abs(expected)))
    error=maxval(abs(actual-expected))/scale
    if (error.gt.threshold) then
       write (*,*) 'FD component test failed: ',trim(label)
       write (*,*) 'relative error:',error
       write (*,*) 'actual:',actual
       write (*,*) 'expected:',expected
       stop 1
    endif
  end subroutine assert_close

end program fd_gauge_components
