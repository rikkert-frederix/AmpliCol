program integrated_kernels_test
  use cs_integrated_kernels
  implicit none
  real(kind=8), parameter :: tol=1d-13
  real(kind=8) :: alpha,c(-2:0),x,correction,regular_gz,subtracted
  real(kind=8) :: lf,current_regular,current_plus,expected
  type(cs_distribution) :: p,k,kt,ka
  integer :: info

  call cs_i_qg(c)
  call assert_close(c(-2),cs_cf_lc,'qg double pole')
  call assert_close(c(-1),1.5d0*cs_cf_lc,'qg single pole')
  call assert_close(c(0),cs_cf_lc*(5d0-cs_pi**2/2d0),'qg finite')

  call cs_i_gg(c)
  call assert_close(c(-2),2d0*cs_ca,'gg double pole')
  call assert_close(c(-1),(11d0/3d0)*cs_ca,'gg single pole')

  call cs_i_gg_ordered(c)
  call assert_close(c(-2),cs_ca,'ordered final-state gg double pole')
  call assert_close(c(-1),(11d0/6d0)*cs_ca,'ordered final-state gg single pole')
  call assert_close(c(0),cs_ca*(50d0/9d0-cs_pi**2/2d0),&
       'ordered final-state gg finite')

  call cs_i_qqbar(c)
  call assert_close(c(-2),0d0,'qqbar double pole')
  call assert_close(c(-1),-(2d0/3d0)*cs_tr,'qqbar single pole')

  alpha=0.1d0
  call cs_i_qg(c)
  call cs_ff_alpha_endpoint(alpha,c,correction,info)
  call assert_true(info.eq.0,'FF alpha endpoint status')
  call assert_close(correction,c(-1)*(alpha-1d0-log(alpha))-&
       c(-2)*log(alpha)**2,'FF alpha endpoint')
  call cs_fi_alpha_endpoint(alpha,c,correction,info)
  call assert_true(info.eq.0,'FI alpha endpoint status')
  call assert_close(correction,-c(-1)*log(alpha)-c(-2)*log(alpha)**2,&
       'FI alpha endpoint')
  call cs_ff_alpha_endpoint(1d0,c,correction,info)
  call assert_close(correction,0d0,'FF alpha endpoint at one')
  call cs_fi_alpha_endpoint(1d0,c,correction,info)
  call assert_close(correction,0d0,'FI alpha endpoint at one')

  call cs_i_qqbar(c)
  call cs_ff_alpha_endpoint(alpha,c,correction,info)
  call assert_close(correction,c(-1)*(alpha-1d0-log(alpha)),&
       'FF alpha qqbar linear term')
  call cs_fi_alpha_endpoint(alpha,c,correction,info)
  call assert_close(correction,-c(-1)*log(alpha),&
       'FI alpha qqbar endpoint')

  x=0.37d0
  call cs_ap_distribution(cs_parton_q,cs_parton_q,x,5,p,info)
  call assert_true(info.eq.0,'Pqq status')
  call assert_close(p%regular,-cs_cf_lc*(1d0+x),'Pqq regular')
  call assert_close(p%plus_one,2d0*cs_cf_lc,'Pqq plus')
  call assert_close(p%delta,1.5d0*cs_cf_lc,'Pqq delta')

  call cs_kbar_distribution(cs_parton_g,cs_parton_g,x,5,k,info)
  call assert_true(info.eq.0,'Kgg status')
  call assert_close(k%plus_log,2d0*cs_ca,'Kgg plus-log')
  call assert_close(k%delta,-(50d0/9d0-cs_pi**2)*cs_ca+(16d0/9d0)*cs_tr*5d0,'Kgg delta')

  call cs_ap_distribution(cs_parton_g,cs_parton_q,x,5,p,info)
  call assert_true(info.eq.0,'P g<-q status')
  call assert_close(p%regular,cs_tr_initial_lc*(x*x+(1d0-x)**2),&
       'P g<-q LC initial-average normalization')

  call cs_kbar_distribution(cs_parton_g,cs_parton_q,x,5,k,info)
  call assert_true(info.eq.0,'K g<-q status')
  call assert_close(k%regular,p%regular*log((1d0-x)/x)+&
       2d0*cs_tr_initial_lc*x*(1d0-x),'K g<-q LC initial-average normalization')

  call cs_if_alpha_distribution(cs_parton_q,cs_parton_q,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_q,cs_parton_q,x,5,p,info)
  call assert_true(info.eq.0,'IF alpha qq status')
  call assert_close(ka%regular,(p%regular+p%plus_one/(1d0-x))*log(alpha)-&
       p%plus_one/(1d0-x)*log((1d0-x+alpha)/(2d0-x)),&
       'IF alpha qq kernel')
  call cs_if_alpha_distribution(cs_parton_g,cs_parton_q,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_g,cs_parton_q,x,5,p,info)
  call assert_true(info.eq.0,'IF alpha gq status')
  call assert_close(ka%regular,p%regular*log(alpha),&
       'IF alpha off-diagonal kernel')
  call cs_if_alpha_distribution(cs_parton_q,cs_parton_q,x,1d0,ka,info)
  call assert_close(ka%regular,0d0,'IF alpha at one')

  call cs_i_qg(c)
  call cs_fi_alpha_terms(c,x,alpha,regular_gz,subtracted,info)
  call assert_true(info.eq.0,'FI alpha terms status')
  call assert_close(regular_gz,-2d0*c(-2)*log(2d0-x)/(1d0-x),&
       'FI alpha regular test-function term')
  call assert_close(subtracted,(2d0*c(-2)*log(1d0-x)+c(-1))/(1d0-x),&
       'FI alpha subtracted term')
  call cs_fi_alpha_terms(c,0.95d0,alpha,regular_gz,subtracted,info)
  call assert_close(regular_gz,0d0,'FI alpha support regular')
  call assert_close(subtracted,0d0,'FI alpha support subtracted')
  call cs_fi_alpha_terms(c,x,1d0,regular_gz,subtracted,info)
  call assert_close(regular_gz,0d0,'FI alpha regular at one')
  call assert_close(subtracted,0d0,'FI alpha subtracted at one')

  call cs_fi_distribution(c,x,1d0,regular_gz,subtracted,info)
  call assert_true(info.eq.0,'FI complete distribution status')
  call assert_close(regular_gz,2d0*c(-2)*log(2d0-x)/(1d0-x),&
       'FI alpha-one regular baseline')
  call assert_close(subtracted,-(2d0*c(-2)*log(1d0-x)+c(-1))/(1d0-x),&
       'FI alpha-one plus baseline')
  call cs_fi_distribution(c,x,alpha,regular_gz,subtracted,info)
  call assert_close(regular_gz,0d0,'FI restricted support regular')
  call assert_close(subtracted,0d0,'FI restricted support plus')
  call cs_fi_distribution(c,0.95d0,alpha,regular_gz,subtracted,info)
  call assert_true(abs(regular_gz).gt.0d0,'FI retained support regular')
  call assert_true(abs(subtracted).gt.0d0,'FI retained support plus')

  call cs_if_tilde_distribution(cs_parton_q,cs_parton_q,x,5,kt,info)
  call assert_true(info.eq.0,'IF tilde qq status')
  call assert_close(kt%regular,-2d0*cs_cf_lc*log(2d0-x)/(1d0-x),&
       'IF tilde qq log(2-x)')
  call assert_close(kt%plus_log_one,2d0*cs_cf_lc,'IF tilde qq plus-log')
  call assert_close(kt%delta,-cs_pi**2*cs_cf_lc/3d0,'IF tilde qq delta')
  call cs_if_tilde_distribution(cs_parton_g,cs_parton_q,x,5,kt,info)
  call assert_close(kt%regular,0d0,'IF tilde off-diagonal regular')
  call assert_close(kt%delta,0d0,'IF tilde off-diagonal delta')

  call cs_ii_tilde_distribution(cs_parton_q,cs_parton_q,x,5,kt,info)
  call assert_true(info.eq.0,'II tilde qq status')
  call assert_close(kt%regular,-cs_cf_lc*(1d0+x)*log(1d0-x),&
       'II tilde qq regular')
  call assert_close(kt%plus_log_one,2d0*cs_cf_lc,'II tilde qq plus-log')
  call assert_close(kt%delta,-2d0*cs_pi**2*cs_cf_lc/3d0,'II tilde qq delta')

  ! Distribution-level comparison with the massless finiteii q->q formula.
  lf=0.41d0
  call cs_ap_distribution(cs_parton_q,cs_parton_q,x,5,p,info)
  call cs_kbar_distribution(cs_parton_q,cs_parton_q,x,5,k,info)
  call cs_ii_tilde_distribution(cs_parton_q,cs_parton_q,x,5,kt,info)
  current_regular=-(lf-log(x))*p%regular+k%regular+kt%regular
  expected=cs_cf_lc*((1d0+x)*lf-2d0*(1d0+x)*log(1d0-x)+(1d0-x))
  call assert_close(current_regular,expected,'assembled II qq regular')
  current_plus=-(lf-log(x))*p%plus_one/(1d0-x)+&
       k%plus_log*log((1d0-x)/x)/(1d0-x)+&
       kt%plus_log_one*log(1d0-x)/(1d0-x)
  expected=cs_cf_lc*(-2d0*lf+4d0*log(1d0-x))/(1d0-x)
  call assert_close(current_plus,expected,'assembled II qq plus')

  call cs_ii_alpha_distribution(cs_parton_g,cs_parton_q,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_g,cs_parton_q,x,5,p,info)
  call assert_true(info.eq.0,'II alpha g<-q status')
  call assert_close(ka%regular,p%regular*log(alpha/(1d0-x)),&
       'II alpha g<-q matches P normalization')
  call cs_ii_alpha_distribution(cs_parton_q,cs_parton_q,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_q,cs_parton_q,x,5,p,info)
  call assert_close(ka%regular,(p%regular+p%plus_one/(1d0-x))*&
       log(alpha/(1d0-x)),'II alpha qq unregularised kernel')
  call cs_ii_alpha_distribution(cs_parton_q,cs_parton_g,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_q,cs_parton_g,x,5,p,info)
  call assert_close(ka%regular,p%regular*log(alpha/(1d0-x)),&
       'II alpha qg unregularised kernel')
  call cs_ii_alpha_distribution(cs_parton_g,cs_parton_g,x,alpha,ka,info)
  call cs_ap_distribution(cs_parton_g,cs_parton_g,x,5,p,info)
  call assert_close(ka%regular,(p%regular+p%plus_one/(1d0-x))*&
       log(alpha/(1d0-x)),'II alpha gg unregularised kernel')
  call cs_ii_alpha_distribution(cs_parton_q,cs_parton_q,0.95d0,alpha,ka,info)
  call assert_close(ka%regular,0d0,'II alpha support')
  call cs_ii_alpha_distribution(cs_parton_q,cs_parton_q,x,1d0,ka,info)
  call assert_close(ka%regular,0d0,'II alpha one')

  call assert_close(cs_fdh_endpoint_shift(cs_parton_q),-0.5d0*cs_cf_lc,'FDH quark shift')
  call assert_close(cs_fdh_endpoint_shift(cs_parton_g),-cs_ca/6d0,'FDH gluon shift')

  write(*,'(a)') 'integrated kernel tests: PASS'

contains

  subroutine assert_close(value,expected,label)
    real(kind=8), intent(in) :: value,expected
    character(len=*), intent(in) :: label
    if (abs(value-expected).gt.tol*max(1d0,abs(expected))) then
       write(*,*) 'FAIL: ',trim(label),value,expected
       stop 1
    endif
  end subroutine assert_close

  subroutine assert_true(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
       write(*,*) 'FAIL: ',trim(label)
       stop 1
    endif
  end subroutine assert_true

end program integrated_kernels_test
