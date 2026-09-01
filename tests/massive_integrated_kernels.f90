program massive_integrated_kernels_test
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use cs_integrated_kernels
  use cs_massive_integrated_kernels
  implicit none

  real(kind=8), parameter :: tol=2d-10
  real(kind=8) :: c(-2:0),cref(-2:0),ell,q2,mi,mk,szone,sx,x,alpha,alpha_near_one
  real(kind=8) :: mu2,mu_ren,mu_fac,expected,gz,g1,value,nan
  real(kind=8),parameter :: scale_factors(2)=[1d-20,1d20]
  type(cs_convolution_kernel) :: kernel,kernel_base
  type(cs_distribution) :: alpha_kernel
  integer :: info,a,b,split,iscale

  ell=0.37d0
  q2=1.0d6
  mi=173d0
  mk=173d0
  alpha=1d0

  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,mi,mk,q2,&
       ell,alpha,cs_scheme_hv,c,info)
  call assert_status(info,'massive-emitter massive-spectator FF')
  call assert_finite(c,'massive-emitter massive-spectator FF')

  ! Roundoff at the alpha=1 boundary must not turn an analytically
  ! non-negative square-root argument into a rejected phase-space point.
  alpha_near_one=nearest(1d0,-1d0)
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,mi,mk,q2,&
       ell,alpha_near_one,cs_scheme_hv,c,info)
  call assert_status(info,'massive FF alpha boundary')
  call assert_finite(c,'massive FF alpha boundary')
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,0d0,mk,q2,&
       ell,alpha_near_one,cs_scheme_hv,c,info)
  call assert_status(info,'massive-spectator quark FF alpha boundary')
  call assert_finite(c,'massive-spectator quark FF alpha boundary')
  call cs_massive_ff_endpoint(cs_parton_g,cs_massive_split_gg,0d0,mk,q2,&
       ell,alpha_near_one,cs_scheme_hv,c,info)
  call assert_status(info,'massive-spectator gluon FF alpha boundary')
  call assert_finite(c,'massive-spectator gluon FF alpha boundary')

  ! Fixed values independently cross-checked against the original MadDipole
  ! finiteterms.f implementation.
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,mi,mk,q2,&
       ell,0.2d0,cs_scheme_hv,c,info)
  call assert_vector_close(c,[0d0,-3.67979151661138282d0,&
       -20.4243030206532836d0],tol,'massive FF reference value')
  cref=c
  do iscale=1,size(scale_factors)
     call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,&
          scale_factors(iscale)*mi,scale_factors(iscale)*mk,&
          scale_factors(iscale)**2*q2,ell,0.2d0,cs_scheme_hv,c,info)
     call assert_status(info,'massive FF scale covariance')
     call assert_vector_close(c,cref,tol,'massive FF scale covariance')
  enddo

  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,mi,0d0,q2,&
       ell,0.2d0,cs_scheme_hv,c,info)
  call assert_status(info,'massive-emitter massless-spectator FF alpha')
  call assert_finite(c,'massive-emitter massless-spectator FF alpha')
  call assert_vector_close(c,[0.75d0,-1.08611650529588033d0,&
       -13.8129815733969448d0],tol,'massive-emitter FF reference value')

  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,0d0,mk,q2,&
       ell,0.2d0,cs_scheme_fdh,c,info)
  call assert_status(info,'massless-quark massive-spectator FF')
  call assert_finite(c,'massless-quark massive-spectator FF')
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,0d0,mk,q2,&
       ell,0.2d0,cs_scheme_hv,cref,info)
  call assert_vector_close(cref,[0.75d0,-0.336116505295877666d0,&
       -8.43442575974700759d0],tol,'massive-spectator FF reference value')
  call assert_close(c(0)-cref(0),-0.75d0,tol,&
       'massive-spectator quark FDH shift')

  do split=cs_massive_split_gg,cs_massive_split_qqbar
     call cs_massive_ff_endpoint(cs_parton_g,split,0d0,mk,q2,ell,0.2d0,&
          cs_scheme_hv,c,info)
     call assert_status(info,'massless-gluon massive-spectator FF')
     call assert_finite(c,'massless-gluon massive-spectator FF')
     if (split.eq.cs_massive_split_gg) then
        call assert_vector_close(c,[1.5d0,0.282187967977586496d0,&
             -13.7349355000855802d0],tol,'massive-spectator gg reference value')
     else
        call assert_vector_close(c,[0d0,-0.333333333333333481d0,&
             -1.38587766799676149d0],tol,&
             'massive-spectator qqbar reference value')
     endif
  enddo

  szone=4.0d5
  mi=173d0
  alpha=0.2d0
  call cs_massive_fi_endpoint(mi,szone,ell,alpha,c,info)
  call assert_status(info,'massive FI endpoint')
  call assert_finite(c,'massive FI endpoint')
  call assert_vector_close(c,[0d0,-2.49718825314288573d0,&
       -8.95900540362021758d0],tol,'massive FI endpoint reference value')
  cref=c
  do iscale=1,size(scale_factors)
     call cs_massive_fi_endpoint(scale_factors(iscale)*mi,&
          scale_factors(iscale)**2*szone,ell,alpha,c,info)
     call assert_status(info,'massive FI endpoint scale covariance')
     call assert_vector_close(c,cref,tol,'massive FI endpoint scale covariance')
  enddo
  mu2=mi*mi/szone
  expected=cs_cf_lc*(1d0+log(mu2/(1d0+mu2)))
  call assert_close(c(-2),0d0,tol,'massive FI double pole')
  call assert_close(c(-1),expected,tol,'massive FI single pole')

  x=0.91d0
  ! The fixed Born invariant is x times the real-emission invariant.
  sx=szone/x
  call cs_massive_fi_convolution(mi,sx,szone,x,alpha,kernel,info)
  call assert_status(info,'massive FI convolution')
  call assert_kernel_finite(kernel,'massive FI convolution')
  call assert_kernel_close(kernel,[-19.746497205397716d0,&
       55.4930722920641486d0,55.4930722920641486d0,0d0],tol,&
       'massive FI convolution reference value')
  kernel_base=kernel
  do iscale=1,size(scale_factors)
     call cs_massive_fi_convolution(scale_factors(iscale)*mi,&
          scale_factors(iscale)**2*sx,scale_factors(iscale)**2*szone,&
          x,alpha,kernel,info)
     call assert_status(info,'massive FI convolution scale covariance')
     call assert_kernel_close(kernel,[kernel_base%regular,kernel_base%plus_z,&
          kernel_base%plus_one,kernel_base%delta],tol,&
          'massive FI convolution scale covariance')
  enddo
  call assert_close(kernel%plus_z,kernel%plus_one,tol,&
       'massive FI plus endpoint')

  mk=173d0
  mu_ren=sqrt(2.5d5)
  mu_fac=mu_ren
  do a=cs_parton_q,cs_parton_g
     do b=cs_parton_q,cs_parton_g
        call cs_massive_if_endpoint(a,b,mk,szone,ell,cs_scheme_hv,5,c,info)
        call assert_status(info,'massive IF endpoint')
        call assert_finite(c,'massive IF endpoint')
        call cs_massive_if_convolution(a,b,mk,sx,szone,x,&
             mu_ren,mu_fac,alpha,5,kernel,info)
        call assert_status(info,'massive IF convolution')
        call assert_kernel_finite(kernel,'massive IF convolution')
     enddo
  enddo

  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,ell,&
       cs_scheme_hv,5,c,info)
  call assert_vector_close(c,[1.5d0,0.108233297501044490d0,&
       0.590895114799644627d0],tol,'massive IF endpoint reference value')
  mu_ren=sqrt(exp(ell)*szone)
  mu_fac=mu_ren
  call cs_massive_if_convolution(cs_parton_q,cs_parton_q,mk,sx,szone,x,&
       mu_ren,mu_fac,alpha,5,kernel,info)
  call assert_kernel_close(kernel,[2.0689712525664170d0,&
       -171.73912070350647d0,-175.268224965703638d0,0d0],tol,&
       'massive IF convolution reference value')
  kernel_base=kernel
  do iscale=1,size(scale_factors)
     call cs_massive_if_convolution(cs_parton_q,cs_parton_q,&
          scale_factors(iscale)*mk,scale_factors(iscale)**2*sx,&
          scale_factors(iscale)**2*szone,x,scale_factors(iscale)*mu_ren,&
          scale_factors(iscale)*mu_fac,alpha,5,kernel,info)
     call assert_status(info,'massive IF convolution scale covariance')
     call assert_kernel_close(kernel,[kernel_base%regular,kernel_base%plus_z,&
          kernel_base%plus_one,kernel_base%delta],tol,&
          'massive IF convolution scale covariance')
  enddo
  call cs_massive_if_convolution(cs_parton_q,cs_parton_g,mk,sx,szone,x,&
       mu_ren,mu_fac,alpha,5,kernel,info)
  ! This AmpliCol value includes the quark/gluon initial-state averaging
  ! ratio; integrated_beam subsequently applies the ordered 1/2 history
  ! weight.
  call assert_kernel_close(kernel,[-10.008572924824385d0,0d0,0d0,0d0],&
       tol,'massive IF qg reference value')
  call cs_massive_if_convolution(cs_parton_g,cs_parton_q,mk,sx,szone,x,&
       mu_ren,mu_fac,alpha,5,kernel,info)
  call assert_kernel_close(kernel,[-1.9271629744640093d0,0d0,0d0,0d0],&
       tol,'massive IF gq reference value')
  call cs_massive_if_endpoint(cs_parton_g,cs_parton_g,mk,szone,ell,&
       cs_scheme_hv,5,c,info)
  call assert_vector_close(c,[3d0,0.216466595002089979d0,&
       1.18179022959928925d0],tol,'massive IF gg endpoint reference value')
  call cs_massive_if_convolution(cs_parton_g,cs_parton_g,mk,sx,szone,x,&
       mu_ren,mu_fac,alpha,5,kernel,info)
  call assert_kernel_close(kernel,[0.26974149760112098d0,&
       -343.47824140701294d0,-350.536449931407276d0,0d0],tol,&
       'massive IF gg reference value')

  ! The massive IF reference formula contains a non-cancelling diagonal
  ! delta term gamma_i log(mu_R^2/mu_F^2).  All regular and plus terms are
  ! factorised at mu_F.
  call cs_massive_if_convolution(cs_parton_q,cs_parton_q,mk,sx,szone,x,&
       mu_fac,mu_fac,alpha,5,kernel_base,info)
  call cs_massive_if_convolution(cs_parton_q,cs_parton_q,mk,sx,szone,x,&
       2d0*mu_fac,mu_fac,alpha,5,kernel,info)
  call assert_close(kernel%regular,kernel_base%regular,tol,'massive IF qq muF regular')
  call assert_close(kernel%plus_z,kernel_base%plus_z,tol,'massive IF qq muF plus-z')
  call assert_close(kernel%plus_one,kernel_base%plus_one,tol,'massive IF qq muF plus-one')
  call assert_close(kernel%delta,cs_gamma(cs_parton_q,5)*log(4d0),tol,&
       'massive IF qq muR/muF delta')
  call cs_massive_if_convolution(cs_parton_g,cs_parton_g,mk,sx,szone,x,&
       2d0*mu_fac,mu_fac,alpha,3,kernel,info)
  call assert_close(kernel%delta,cs_gamma(cs_parton_g,3)*log(4d0),tol,&
       'massive IF gg muR/muF delta')
  call cs_massive_if_convolution(cs_parton_q,cs_parton_g,mk,sx,szone,x,&
       2d0*mu_fac,mu_fac,alpha,5,kernel,info)
  call assert_close(kernel%delta,0d0,tol,'massive IF off-diagonal muR/muF delta')

  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,ell,&
       cs_scheme_hv,5,cref,info)
  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,ell,&
       cs_scheme_fdh,5,c,info)
  call assert_close(c(0)-cref(0),-0.75d0,tol,'massive IF quark FDH shift')
  call cs_massive_if_endpoint(cs_parton_g,cs_parton_g,mk,szone,ell,&
       cs_scheme_hv,5,cref,info)
  call cs_massive_if_endpoint(cs_parton_g,cs_parton_g,mk,szone,ell,&
       cs_scheme_fdh,5,c,info)
  call assert_close(c(0)-cref(0),-0.5d0,tol,'massive IF gluon FDH shift')

  ! The alpha-dependent regular term has a smooth massless-spectator
  ! limit.  This checks its LC normalization against the established
  ! massless IF kernels for every crossed flavour channel.
  do a=cs_parton_q,cs_parton_g
     do b=cs_parton_q,cs_parton_g
        call cs_massive_if_convolution(a,b,1d-3,sx,szone,x,&
             mu_ren,mu_fac,alpha,5,kernel,info)
        call assert_status(info,'massive IF alpha massless limit')
        call cs_massive_if_convolution(a,b,1d-3,sx,szone,x,&
             mu_ren,mu_fac,1d0,5,kernel_base,info)
        call assert_status(info,'massive IF alpha-one massless limit')
        call cs_if_alpha_distribution(a,b,x,5,alpha,alpha_kernel,info)
        call assert_status(info,'massless IF alpha comparison')
        call assert_close(kernel%regular-kernel_base%regular,&
             alpha_kernel%regular,2d-8,'massive IF alpha massless limit')
     enddo
  enddo

  ! The generalized plus distribution must retain distinct z and endpoint
  ! values for a massive spectator.
  call cs_massive_if_convolution(cs_parton_q,cs_parton_q,mk,sx,szone,x,&
       sqrt(2.5d5),sqrt(2.5d5),alpha,5,kernel,info)
  call assert_true(abs(kernel%plus_z-kernel%plus_one).gt.1d-10,&
       'massive IF has distinct plus endpoint')
  gz=1.7d0
  g1=0.9d0
  value=cs_apply_convolution(kernel,gz,g1)
  expected=(kernel%regular+kernel%plus_z)*gz-kernel%plus_one*g1+kernel%delta*g1
  call assert_close(value,expected,tol,'massive convolution action')
  kernel=cs_convolution_kernel()
  kernel%regular=huge(1d0)
  value=cs_apply_convolution(kernel,2d0,1d0,info)
  call assert_true(info.eq.-20 .and. value.eq.0d0,&
       'unsafe massive convolution action rejected before overflow')

  call cs_massive_fi_endpoint(0d0,szone,ell,alpha,c,info)
  call assert_true(info.ne.0,'massless input rejected by massive FI kernel')
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,mi,mk,&
       (mi+mk)**2,ell,alpha,cs_scheme_hv,c,info)
  call assert_true(info.ne.0,'FF threshold rejected')

  nan=ieee_value(0d0,ieee_quiet_nan)
  call cs_massive_ff_endpoint(cs_parton_q,cs_massive_split_qg,nan,mk,q2,&
       ell,alpha,cs_scheme_hv,c,info)
  call assert_true(info.eq.-20 .and. all(c.eq.0d0),'non-finite massive FF input rejected')
  call cs_massive_fi_endpoint(mi,nan,ell,alpha,c,info)
  call assert_true(info.eq.-20 .and. all(c.eq.0d0),'non-finite massive FI endpoint input rejected')
  call cs_massive_fi_convolution(mi,sx,szone,nan,alpha,kernel,info)
  call assert_true(info.eq.-20,'non-finite massive FI convolution input rejected')
  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,nan,&
       cs_scheme_hv,5,c,info)
  call assert_true(info.eq.-20 .and. all(c.eq.0d0),'non-finite massive IF endpoint input rejected')
  call cs_massive_if_convolution(cs_parton_q,cs_parton_q,mk,sx,szone,x,&
       mu_ren,nan,alpha,5,kernel,info)
  call assert_true(info.eq.-20,'non-finite massive IF convolution input rejected')
  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,ell,999,5,c,info)
  call assert_true(info.ne.0,'unknown massive endpoint scheme rejected')
  call cs_massive_if_endpoint(cs_parton_q,cs_parton_q,mk,szone,ell,&
       cs_scheme_hv,7,c,info)
  call assert_true(info.ne.0,'invalid active-flavour count rejected')

  write(*,'(a)') 'massive integrated kernel tests: PASS'

contains

  subroutine assert_status(status,label)
    integer, intent(in) :: status
    character(len=*), intent(in) :: label
    if (status.ne.0) then
       write(*,*) 'FAIL status: ',trim(label),status
       stop 1
    endif
  end subroutine assert_status

  subroutine assert_finite(values,label)
    real(kind=8), intent(in) :: values(-2:0)
    character(len=*), intent(in) :: label
    if (.not.all(ieee_is_finite(values))) then
       write(*,*) 'FAIL non-finite: ',trim(label),values
       stop 1
    endif
  end subroutine assert_finite

  subroutine assert_kernel_finite(k,label)
    type(cs_convolution_kernel), intent(in) :: k
    character(len=*), intent(in) :: label
    if (.not.all(ieee_is_finite([k%regular,k%plus_z,k%plus_one,k%delta]))) then
       write(*,*) 'FAIL non-finite kernel: ',trim(label)
       stop 1
    endif
  end subroutine assert_kernel_finite

  subroutine assert_kernel_close(k,want,relative_tolerance,label)
    type(cs_convolution_kernel), intent(in) :: k
    real(kind=8), intent(in) :: want(4),relative_tolerance
    character(len=*), intent(in) :: label
    real(kind=8) :: actual(4)
    actual=[k%regular,k%plus_z,k%plus_one,k%delta]
    if (any(abs(actual-want).gt.relative_tolerance*max(1d0,abs(want)))) then
       write(*,*) 'FAIL: ',trim(label)
       write(*,*) actual
       write(*,*) want
       stop 1
    endif
  end subroutine assert_kernel_close

  subroutine assert_close(actual,want,relative_tolerance,label)
    real(kind=8), intent(in) :: actual,want,relative_tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-want).gt.relative_tolerance*max(1d0,abs(want))) then
       write(*,*) 'FAIL: ',trim(label),actual,want
       stop 1
    endif
  end subroutine assert_close

  subroutine assert_vector_close(actual,want,relative_tolerance,label)
    real(kind=8), intent(in) :: actual(-2:0),want(-2:0),relative_tolerance
    character(len=*), intent(in) :: label
    if (any(abs(actual-want).gt.relative_tolerance*max(1d0,abs(want)))) then
       write(*,*) 'FAIL: ',trim(label)
       write(*,*) actual
       write(*,*) want
       stop 1
    endif
  end subroutine assert_vector_close

  subroutine assert_true(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
       write(*,*) 'FAIL: ',trim(label)
       stop 1
    endif
  end subroutine assert_true

end program massive_integrated_kernels_test
