program phase_space_onebody_regression
  use phase_space_base
  use phase_space_onebody_mod
  implicit none
  integer,parameter :: dp=kind(1d0),n=3
  real(kind=dp),parameter :: pi=3.1415926535897932384626433832795d0
  real(kind=dp),parameter :: sqrts=14000d0,mass=125d0
  type(phase_space_onebody) :: phase_space
  type(psv) :: point
  integer,dimension(n) :: order
  real(kind=dp),dimension(n) :: masses,pt_cut,rap_cut
  real(kind=dp),dimension(n,n) :: dr_cut,sqrt_s_min
  real(kind=dp),dimension(0:3,n) :: generated
  real(kind=dp) :: tau,ymax,expected_jac,expected_x

  order=[1,2,3]
  masses=[0d0,0d0,mass]
  pt_cut=-1d0
  rap_cut=-1d0
  dr_cut=-1d0
  sqrt_s_min=-1d0
  call phase_space%init(sqrts,n,masses,order,pt_cut,rap_cut,dr_cut,&
       sqrt_s_min,.false.,.true.)
  if (phase_space%ndim.ne.1 .or. phase_space%ndim_extra.ne.0) then
     write (*,*) 'Unexpected one-body phase-space dimension',&
          phase_space%ndim,phase_space%ndim_extra
     stop 1
  endif
  if (.not.phase_space%can_invert_momenta) then
     write (*,*) 'One-body phase space must support inverse maps'
     stop 1
  endif

  allocate(point%x(1),point%p(0:3,n))
  expected_x=0.37123456789d0
  point%x(1)=expected_x
  call phase_space%generate_momenta(point)
  tau=(mass/sqrts)**2
  ymax=-0.5d0*log(tau)
  expected_jac=2d0*ymax*pi/(mass**2*sqrts**2)
  call assert_close('Jacobian',point%jac,expected_jac,2d-14)
  call assert_close('Bjorken product',point%xbjrk(1)*point%xbjrk(2),tau,2d-14)
  call assert_close('final mass',minkowski_square(point%p(:,3)),mass**2,2d-14)
  if (maxval(abs(point%p(:,3)-point%p(:,1)-point%p(:,2))).gt.1d-12) then
     write (*,*) 'One-body momenta do not conserve four-momentum'
     stop 1
  endif

  generated=point%p
  point%x=-1d0
  point%xbjrk=-1d0
  point%jac=-1d0
  call phase_space%compute_x_from_momenta(point)
  call assert_close('inverse random coordinate',point%x(1),expected_x,2d-13)
  call assert_close('inverse Jacobian',point%jac,expected_jac,2d-14)
  if (maxval(abs(point%p-generated)).gt.0d0) then
     write (*,*) 'Inverse map changed the input momenta'
     stop 1
  endif

  call phase_space%cleanup()
  write (*,'(a)') 'One-body phase-space regression passed'

contains

  real(kind=dp) function minkowski_square(p)
    real(kind=dp),dimension(0:3),intent(in) :: p
    minkowski_square=p(0)**2-p(1)**2-p(2)**2-p(3)**2
  end function minkowski_square

  subroutine assert_close(label,actual,expected,tolerance)
    character(len=*),intent(in) :: label
    real(kind=dp),intent(in) :: actual,expected,tolerance
    if (abs(actual-expected).gt.tolerance*max(1d0,abs(expected))) then
       write (*,*) trim(label),' mismatch:',actual,expected
       stop 1
    endif
  end subroutine assert_close

end program phase_space_onebody_regression
