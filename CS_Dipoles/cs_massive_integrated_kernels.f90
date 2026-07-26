module cs_massive_integrated_kernels
  ! Massive Catani--Seymour endpoint and convolution kernels for a
  ! massless unresolved parton.  The formulae follow the conventions of
  ! Utilities/Docs/finiteterms.f, with the overall alpha_s/(2*pi), the
  ! Born contribution, and the colour-correlation history weight removed.
  !
  ! Unlike cs_distribution, cs_convolution_kernel permits the value
  ! multiplying g(z) in a plus distribution to differ from the value
  ! multiplying g(1).  This is required by the massive-spectator kernels.
  use cs_dipole_mappings, only: dp
  use cs_integrated_kernels, only: cs_pi,cs_ca,cs_cf_lc,cs_cf_initial_qg,&
       cs_tr_initial_lc,&
       cs_parton_q,cs_parton_g,cs_scheme_hv,cs_scheme_fdh
  implicit none
  private

  integer, parameter, public :: cs_massive_split_qg=1
  integer, parameter, public :: cs_massive_split_gg=2
  integer, parameter, public :: cs_massive_split_qqbar=3

  type, public :: cs_convolution_kernel
     real(dp) :: regular=0.0_dp
     real(dp) :: plus_z=0.0_dp
     real(dp) :: plus_one=0.0_dp
     real(dp) :: delta=0.0_dp
  end type cs_convolution_kernel

  public :: cs_massive_ff_endpoint
  public :: cs_massive_fi_endpoint,cs_massive_fi_convolution
  public :: cs_massive_if_endpoint,cs_massive_if_convolution
  public :: cs_apply_convolution,cs_dilog

contains

  pure real(dp) function cs_apply_convolution(kernel,gz,g1) result(value)
    type(cs_convolution_kernel), intent(in) :: kernel
    real(dp), intent(in) :: gz,g1
    value=(kernel%regular+kernel%plus_z)*gz-kernel%plus_one*g1+kernel%delta*g1
  end function cs_apply_convolution

  pure real(dp) function cs_kallen(x,y,z) result(value)
    real(dp), intent(in) :: x,y,z
    value=x*x+y*y+z*z-2.0_dp*(x*y+x*z+y*z)
  end function cs_kallen

  pure real(dp) function cs_dilog(x) result(value)
    ! CERNLIB DDILOG, expressed in modern Fortran.  All arguments used by
    ! the massive kernels are real and lie on the real branch.
    real(dp), intent(in) :: x
    real(dp), parameter :: c(0:19)=[&
         0.42996693560813697_dp, 0.40975987533077105_dp,&
        -0.01858843665014592_dp, 0.00145751084062268_dp,&
        -0.00014304184442340_dp, 0.00001588415541880_dp,&
        -0.00000190784959387_dp, 0.00000024195180854_dp,&
        -0.00000003193341274_dp, 0.00000000434545080_dp,&
        -0.00000000060578480_dp, 0.00000000008612098_dp,&
        -0.00000000001244332_dp, 0.00000000000182256_dp,&
        -0.00000000000027007_dp, 0.00000000000004042_dp,&
        -0.00000000000000610_dp, 0.00000000000000093_dp,&
        -0.00000000000000014_dp, 0.00000000000000002_dp]
    real(dp) :: t,y,s,a,h,alfa,b0,b1,b2
    integer :: i

    if (x == 1.0_dp) then
       value=cs_pi*cs_pi/6.0_dp
       return
    elseif (x == -1.0_dp) then
       value=-cs_pi*cs_pi/12.0_dp
       return
    endif
    t=-x
    if (t <= -2.0_dp) then
       y=-1.0_dp/(1.0_dp+t)
       s=1.0_dp
       a=-cs_pi*cs_pi/3.0_dp+0.5_dp*(log(-t)**2-log(1.0_dp+1.0_dp/t)**2)
    elseif (t < -1.0_dp) then
       y=-1.0_dp-t
       s=-1.0_dp
       a=log(-t)
       a=-cs_pi*cs_pi/6.0_dp+a*(a+log(1.0_dp+1.0_dp/t))
    elseif (t <= -0.5_dp) then
       y=-(1.0_dp+t)/t
       s=1.0_dp
       a=log(-t)
       a=-cs_pi*cs_pi/6.0_dp+a*(-0.5_dp*a+log(1.0_dp+t))
    elseif (t < 0.0_dp) then
       y=-t/(1.0_dp+t)
       s=-1.0_dp
       a=0.5_dp*log(1.0_dp+t)**2
    elseif (t <= 1.0_dp) then
       y=t
       s=1.0_dp
       a=0.0_dp
    else
       y=1.0_dp/t
       s=-1.0_dp
       a=cs_pi*cs_pi/6.0_dp+0.5_dp*log(t)**2
    endif
    h=2.0_dp*y-1.0_dp
    alfa=2.0_dp*h
    b1=0.0_dp
    b2=0.0_dp
    do i=19,0,-1
       b0=c(i)+alfa*b1-b2
       b2=b1
       b1=b0
    enddo
    value=-(s*(b0-h*b2)+a)
  end function cs_dilog

  pure subroutine laurent_from_finite(f0,f1,f2,factual,normalisation,coeff)
    ! If F(L)=A L^2+B L+C is the finite coefficient after expanding the
    ! universal dimensional factor, the corresponding poles are 2A and B.
    real(dp), intent(in) :: f0,f1,f2,factual,normalisation
    real(dp), intent(out) :: coeff(-2:0)
    real(dp) :: a,b
    a=0.5_dp*(f2-2.0_dp*f1+f0)
    b=f1-f0-a
    coeff(-2)=normalisation*2.0_dp*a
    coeff(-1)=normalisation*b
    coeff(0)=normalisation*factual
  end subroutine laurent_from_finite

  pure subroutine cs_massive_ff_endpoint(parton,split,mi,mk,q2,ell,alpha,scheme,coeff,info)
    integer, intent(in) :: parton,split,scheme
    real(dp), intent(in) :: mi,mk,q2,ell,alpha
    real(dp), intent(out) :: coeff(-2:0)
    integer, intent(out) :: info
    real(dp) :: f0,f1,f2,fa,norm

    coeff=0.0_dp
    info=0
    if (mi < 0.0_dp .or. mk < 0.0_dp .or. q2 <= 0.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (q2 <= (mi+mk)**2) then
       info=-4
       return
    endif

    if (parton == cs_parton_q .and. split == cs_massive_split_qg) then
       norm=cs_cf_lc
    elseif (parton == cs_parton_g .and. split == cs_massive_split_gg .and. mi == 0.0_dp) then
       ! One ordered side.  The independent history colour weight is
       ! applied by integrated_dipoles.
       norm=0.5_dp*cs_ca
    elseif (parton == cs_parton_g .and. split == cs_massive_split_qqbar .and. mi == 0.0_dp) then
       ! finiteff carries the conventional T_R/C_A normalization.  Restore
       ! the one-flavour g -> q qbar primitive before the independent
       ! history colour weight is applied.
       norm=cs_ca
    else
       info=-2
       return
    endif

    call ff_finite_base(parton,split,mi,mk,q2,0.0_dp,alpha,scheme,f0,info)
    if (info /= 0) return
    call ff_finite_base(parton,split,mi,mk,q2,1.0_dp,alpha,scheme,f1,info)
    if (info /= 0) return
    call ff_finite_base(parton,split,mi,mk,q2,2.0_dp,alpha,scheme,f2,info)
    if (info /= 0) return
    call ff_finite_base(parton,split,mi,mk,q2,ell,alpha,scheme,fa,info)
    if (info /= 0) return
    call laurent_from_finite(f0,f1,f2,fa,norm,coeff)
  end subroutine cs_massive_ff_endpoint

  pure subroutine ff_finite_base(parton,split,mi,mk,q2,ell,alpha,scheme,value,info)
    integer, intent(in) :: parton,split,scheme
    real(dp), intent(in) :: mi,mk,q2,ell,alpha
    real(dp), intent(out) :: value
    integer, intent(out) :: info
    real(dp) :: mui,muk,mui2,muk2,sq,qik2,lam,v,rhoi,rhok,rho,rs
    real(dp) :: yp,xp,xm,x,a,b,c,d,yl,e,term

    value=0.0_dp
    info=0
    sq=sqrt(q2)
    mui2=mi*mi/q2
    muk2=mk*mk/q2
    mui=mi/sq
    muk=mk/sq
    qik2=q2-mi*mi-mk*mk

    if (parton == cs_parton_q) then
       if (mi > 0.0_dp .and. mk > 0.0_dp) then
          lam=cs_kallen(1.0_dp,mui2,muk2)
          if (lam <= 0.0_dp .or. qik2 <= 0.0_dp) then
             info=-4
             return
          endif
          v=sqrt(lam)/(1.0_dp-mui2-muk2)
          rhoi=sqrt((1.0_dp-v+2.0_dp*mui2/(1.0_dp-mui2-muk2))/&
               (1.0_dp+v+2.0_dp*mui2/(1.0_dp-mui2-muk2)))
          rhok=sqrt((1.0_dp-v+2.0_dp*muk2/(1.0_dp-mui2-muk2))/&
               (1.0_dp+v+2.0_dp*muk2/(1.0_dp-mui2-muk2)))
          rho=sqrt((1.0_dp-v)/(1.0_dp+v))
          value=0.5_dp*(6.0_dp+2.0_dp*ell+&
               2.0_dp*mk*((4.0_dp*mk-2.0_dp*sq)/qik2+1.0_dp/(mk-sq))-&
               2.0_dp*cs_pi**2/(3.0_dp*v)-&
               4.0_dp*log(((sq-mk)**2-mi*mi)/q2)+&
               2.0_dp*log(mi*(1.0_dp-mk/sq)/sq)-&
               4.0_dp*mi*mi*log(mi/(sq-mk))/qik2-&
               (-4.0_dp*log(rho*rho)*log(1.0_dp+rho*rho)+&
               log(rhoi*rhoi)**2+log(rhok*rhok)**2-&
               4.0_dp*log(rho)*(ell+2.0_dp*log(q2/qik2)))/(2.0_dp*v)-&
               2.0_dp*(-2.0_dp*cs_dilog(rho*rho)+cs_dilog(1.0_dp-rhoi*rhoi)+&
               cs_dilog(1.0_dp-rhok*rhok))/v)

          if (alpha < 1.0_dp) then
             a=2.0_dp*muk/(1.0_dp-mui2-muk2)
             b=2.0_dp*(1.0_dp-muk)/(1.0_dp-mui2-muk2)
             c=2.0_dp*(1.0_dp-muk)*muk/(1.0_dp-mui2-muk2)
             d=(1.0_dp-mui2-muk2)/2.0_dp
             xp=(-mui2+(1.0_dp-muk)**2+sqrt(lam))/(1.0_dp-mui2-muk2)
             xm=(-mui2+(1.0_dp-muk)**2-sqrt(lam))/(1.0_dp-mui2-muk2)
             yp=1.0_dp-2.0_dp*(1.0_dp-muk)*muk/(1.0_dp-mui2-muk2)
             term=((4.0_dp*mui2*muk2)/&
                  ((mui2-(1.0_dp-muk)**2)*(1.0_dp-mui2-muk2))+1.0_dp/yp-alpha*yp)*&
                  (yp-alpha*yp)
             if (term < 0.0_dp) then
                info=-5
                return
             endif
             x=yp-alpha*yp+sqrt(term)
             term=1.0_dp/(1.0_dp-muk)-2.0_dp*(2.0_dp-2.0_dp*mui2-muk)/&
                  (1.0_dp-mui2-muk2)+1.5_dp*(1.0_dp+alpha*yp)+&
                  mui2*(1.0_dp-alpha*yp)/(2.0_dp*(mui2+alpha*(1.0_dp-mui2-muk2)*yp))-&
                  2.0_dp*log(alpha*(1.0_dp-mui2-muk2)*yp/&
                  (-mui2+(1.0_dp-muk)**2))+&
                  (1.0_dp+mui2-muk2)*log((mui2+alpha*(1.0_dp-mui2-muk2)*yp)/&
                  (1.0_dp-muk)**2)/(2.0_dp*(1.0_dp-mui2-muk2))
             term=term+2.0_dp*(-log(b)*log(a*(-b+xm)/((a+b)*xm))+&
                  log(b-x)*log((a+x)*(-b+xm)/((a+b)*(-x+xm)))+&
                  log((c+xm)/(a+xm))*log((-x+xm)/xm)+&
                  log(a*(b-xp))*log(xp)+0.5_dp*log((a+x)/a)*&
                  log(a*(a+x)*(a+xp)**2)-&
                  log(c)*log((a-c)*xp/(a*(c+xp)))+&
                  log(d)*log((a+x)*xm*xp/(a*(-x+xm)*(-x+xp)))-&
                  log((a+x)*(b-xp))*log(-x+xp)+&
                  log(c+x)*log((a-c)*(-x+xp)/((a+x)*(c+xp)))-&
                  cs_dilog(b/(a+b))+cs_dilog(c/(-a+c))+&
                  cs_dilog((b-x)/(a+b))-cs_dilog((c+x)/(-a+c))+&
                  cs_dilog(b/(b-xm))-cs_dilog((b-x)/(b-xm))-&
                  cs_dilog(xm/(a+xm))+cs_dilog(xm/(c+xm))+&
                  cs_dilog((-x+xm)/(a+xm))-cs_dilog((-x+xm)/(c+xm))+&
                  cs_dilog(a/(a+xp))-cs_dilog((a+x)/(a+xp))-&
                  cs_dilog(xp/(-b+xp))-cs_dilog(c/(c+xp))+&
                  cs_dilog((c+x)/(c+xp))+cs_dilog((-x+xp)/(-b+xp)))/v
             value=value+term
          endif
       elseif (mi > 0.0_dp .and. mk == 0.0_dp) then
          value=(72.0_dp+6.0_dp*ell*(4.0_dp+ell)-11.0_dp*cs_pi**2+&
               24.0_dp*mui2*log(mui2)/(-1.0_dp+mui2)+&
               6.0_dp*(4.0_dp*log(1.0_dp-mui2)**2+&
               (2.0_dp+2.0_dp*ell-log(mui2))*log(mui2)-&
               4.0_dp*log(1.0_dp-mui2)*(2.0_dp+ell+log(mui2)))-&
               24.0_dp*cs_dilog(1.0_dp-mui2))/24.0_dp
          if (alpha < 1.0_dp) then
             value=value-2.0_dp*log(alpha)+&
                  2.0_dp*log(alpha+(1.0_dp-alpha)*mui2)/(1.0_dp-mui2)+&
                  0.5_dp*(-2.0_dp+3.0_dp*alpha-&
                  alpha/(alpha+(1.0_dp-alpha)*mui2)-&
                  (3.0_dp-mui2)*log(alpha+(1.0_dp-alpha)*mui2)/(1.0_dp-mui2))+&
                  2.0_dp*(-log(alpha)*log(mui2)-cs_dilog((-1.0_dp+mui2)/mui2)+&
                  cs_dilog(alpha*(-1.0_dp+mui2)/mui2))
          endif
       elseif (mi == 0.0_dp .and. mk > 0.0_dp) then
          rs=0.0_dp
          if (scheme == cs_scheme_fdh) rs=-0.5_dp
          value=(36.0_dp*ell*(1.0_dp+muk)+6.0_dp*ell**2*(1.0_dp+muk)-&
               11.0_dp*cs_pi**2+24.0_dp*(5.0_dp+rs)+&
               muk*(-11.0_dp*cs_pi**2+24.0_dp*(2.0_dp+rs))-&
               36.0_dp*(1.0_dp+muk)*log((1.0_dp-muk)**2)-&
               6.0_dp*(1.0_dp+muk)*(-4.0_dp*log(1.0_dp-muk2)**2+&
               log(muk2)*(-2.0_dp*ell+log(muk2))+&
               4.0_dp*log(1.0_dp-muk2)*(ell+log(muk2)))-&
               24.0_dp*(1.0_dp+muk)*cs_dilog(1.0_dp-muk2))/&
               (24.0_dp*(1.0_dp+muk))
          if (alpha < 1.0_dp) then
             yp=(1.0_dp-muk2)/(1.0_dp+muk2)
             term=(1.0_dp-alpha)*(1.0_dp-alpha*yp*yp)
             if (term < 0.0_dp) then
                info=-5
                return
             endif
             xp=yp*(1.0_dp-alpha)+sqrt(term)
             value=value-1.5_dp*((1.0_dp-alpha)*yp+log(alpha))-&
                  2.0_dp*log((1.0_dp-xp+yp)/(1.0_dp+yp))**2+&
                  log((1.0_dp+2.0_dp*xp*yp-yp*yp)/&
                  ((1.0_dp+xp-yp)*(1.0_dp-xp+yp)))**2+&
                  4.0_dp*(log((1.0_dp+xp-yp)/(1.0_dp-yp))*log((1.0_dp+yp)/2.0_dp)+&
                  log((1.0_dp+yp)/(2.0_dp*yp))*&
                  log((1.0_dp+2.0_dp*xp*yp-yp*yp)/(1.0_dp-yp*yp))-&
                  cs_dilog((1.0_dp-yp)/2.0_dp)+cs_dilog((1.0_dp+xp-yp)/2.0_dp)+&
                  cs_dilog((1.0_dp-yp)/(1.0_dp+yp))-&
                  cs_dilog((1.0_dp+2.0_dp*xp*yp-yp*yp)/(1.0_dp+yp)**2))
          endif
       else
          info=-2
       endif

    elseif (parton == cs_parton_g .and. mi == 0.0_dp .and. mk > 0.0_dp) then
       if (split == cs_massive_split_qqbar) then
          e=(-1.0_dp+muk)*(-6.0_dp*ell*(1.0_dp+muk)-4.0_dp*(4.0_dp+muk))+&
               12.0_dp*(-1.0_dp+muk2)*log(1.0_dp-muk)+&
               12.0_dp*muk2*log(2.0_dp*muk/(1.0_dp+muk))
          value=e/(18.0_dp*cs_ca*(-1.0_dp+muk2))
          if (alpha < 1.0_dp) then
             value=value-(-(-1.0_dp+alpha)*(-1.0_dp+muk)**2+log(alpha)+&
                  muk2*(log(4.0_dp)-log(alpha)+2.0_dp*log(muk)-&
                  2.0_dp*log(1.0_dp+muk)-&
                  2.0_dp*log(1.0_dp+alpha-2.0_dp*alpha/(1.0_dp+muk))))/&
                  (3.0_dp*(-1.0_dp+muk2)*cs_ca)
          endif
       elseif (split == cs_massive_split_gg) then
          rs=0.0_dp
          if (scheme == cs_scheme_fdh) rs=-1.0_dp/6.0_dp
          e=200.0_dp+66.0_dp*ell+9.0_dp*ell**2-132.0_dp*muk/(1.0_dp+muk)-&
               15.0_dp*cs_pi**2-132.0_dp*log(1.0_dp-muk)-&
               24.0_dp*muk2*log(2.0_dp*muk/(1.0_dp+muk))/(-1.0_dp+muk2)+&
               66.0_dp*log(1.0_dp-muk2)+36.0_dp*rs+&
               9.0_dp*(4.0_dp*(ell-log(muk))*log(muk)-&
               2.0_dp*(ell+2.0_dp*log(muk))*log(1.0_dp-muk2)+&
               log(1.0_dp-muk2)**2)-36.0_dp*cs_dilog(1.0_dp-muk2)
          value=e/18.0_dp
          if (alpha < 1.0_dp) then
             term=(-1.0_dp+muk)**2*(alpha**2*(-1.0_dp+muk)**2+&
                  (1.0_dp+muk)**2-2.0_dp*alpha*(1.0_dp+muk2))
             if (term < 0.0_dp .and. abs(term) > 100.0_dp*epsilon(1.0_dp)) then
                info=-5
                return
             endif
             yl=1.0_dp+alpha*(-1.0_dp+muk)**2-muk2-sqrt(max(0.0_dp,term))
             e=11.0_dp*(-2.0_dp+2.0_dp*muk+yl)**2/&
                  ((-1.0_dp+muk2)*(-2.0_dp+yl))-44.0_dp*log(2.0_dp-2.0_dp*muk)-&
                  22.0_dp*log(muk)+24.0_dp*log(2.0_dp/(1.0_dp+muk))*&
                  (log(2.0_dp/(1.0_dp+muk))+2.0_dp*log(1.0_dp+muk))+&
                  2.0_dp*((-11.0_dp+15.0_dp*muk2)*log(2.0_dp*muk)+&
                  4.0_dp*muk2*(-log(-8.0_dp*(-1.0_dp+muk)*muk2)+&
                  log((-2.0_dp+yl)**2+4.0_dp*muk2*(-1.0_dp+yl)))+&
                  (11.0_dp-15.0_dp*muk2)*log(2.0_dp-yl))/(-1.0_dp+muk2)+&
                  22.0_dp*log(2.0_dp-2.0_dp*muk2-yl)+22.0_dp*log(yl)-&
                  12.0_dp*(4.0_dp*log(1.0_dp-yl/2.0_dp)*&
                  log(-yl/(-1.0_dp+muk2))-log(-yl/(-1.0_dp+muk2))**2+&
                  log(-2.0_dp*(-2.0_dp+2.0_dp*muk2+yl)/&
                  ((-1.0_dp+muk2)*(-2.0_dp+yl)))**2+&
                  2.0_dp*log(-yl/(-1.0_dp+muk2))*&
                  (log(-2.0_dp*(-2.0_dp+2.0_dp*muk2+yl)/&
                  ((-1.0_dp+muk2)*(-2.0_dp+yl)))-&
                  2.0_dp*log(1.0_dp+yl/(-2.0_dp+2.0_dp*muk2))))+&
                  48.0_dp*cs_dilog(1.0_dp-muk)-48.0_dp*cs_dilog(1.0_dp/(1.0_dp+muk))-&
                  48.0_dp*cs_dilog(yl/2.0_dp)+48.0_dp*cs_dilog(yl/(2.0_dp-2.0_dp*muk2))
             value=value-e/6.0_dp
          endif
       else
          info=-2
       endif
    else
       info=-2
    endif
  end subroutine ff_finite_base

  pure subroutine cs_massive_fi_endpoint(mi,szone,ell,alpha,coeff,info)
    real(dp), intent(in) :: mi,szone,ell,alpha
    real(dp), intent(out) :: coeff(-2:0)
    integer, intent(out) :: info
    real(dp) :: f0,f1,f2,fa
    coeff=0.0_dp
    info=0
    if (mi <= 0.0_dp .or. szone <= 0.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    call fi_endpoint_base(mi,szone,0.0_dp,alpha,f0)
    call fi_endpoint_base(mi,szone,1.0_dp,alpha,f1)
    call fi_endpoint_base(mi,szone,2.0_dp,alpha,f2)
    call fi_endpoint_base(mi,szone,ell,alpha,fa)
    call laurent_from_finite(f0,f1,f2,fa,cs_cf_lc,coeff)
  end subroutine cs_massive_fi_endpoint

  pure subroutine fi_endpoint_base(mi,szone,ell,alpha,value)
    real(dp), intent(in) :: mi,szone,ell,alpha
    real(dp), intent(out) :: value
    real(dp) :: mu2
    mu2=mi*mi/szone
    value=(1.0_dp+log(mu2/(1.0_dp+mu2)))*ell-&
         2.0_dp*cs_dilog(-mu2)-cs_pi**2/3.0_dp+2.0_dp+&
         0.5_dp*log(mu2)**2+0.5_dp*log(1.0_dp+mu2)**2-&
         2.0_dp*log(mu2)*log(1.0_dp+mu2)+log(mu2)
    if (alpha < 1.0_dp) value=value+2.0_dp*log(alpha)*&
         (log((1.0_dp+mu2)/mu2)-1.0_dp)
  end subroutine fi_endpoint_base

  pure subroutine cs_massive_fi_convolution(mi,s,szone,x,alpha,kernel,info)
    real(dp), intent(in) :: mi,s,szone,x,alpha
    type(cs_convolution_kernel), intent(out) :: kernel
    integer, intent(out) :: info
    real(dp) :: mu2,mu2one
    kernel=cs_convolution_kernel()
    info=0
    if (mi <= 0.0_dp .or. s <= 0.0_dp .or. szone <= 0.0_dp) then
       info=-1
       return
    endif
    if (x <= 0.0_dp .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    if (x <= 1.0_dp-alpha) return
    mu2=mi*mi/s
    mu2one=mi*mi/szone
    kernel%regular=cs_cf_lc*((1.0_dp-x)/(2.0_dp*(1.0_dp-x+mu2)**2)+&
         2.0_dp*log((2.0_dp-x+mu2)*mu2one/&
         ((1.0_dp+mu2one)*(1.0_dp-x+mu2)))/(1.0_dp-x))
    kernel%plus_z=cs_cf_lc*2.0_dp*(log((1.0_dp+mu2one)/mu2one)-1.0_dp)/(1.0_dp-x)
    kernel%plus_one=kernel%plus_z
  end subroutine cs_massive_fi_convolution

  pure subroutine cs_massive_if_endpoint(a,b,mk,szone,ell,scheme,nf,coeff,info)
    integer, intent(in) :: a,b,scheme,nf
    real(dp), intent(in) :: mk,szone,ell
    real(dp), intent(out) :: coeff(-2:0)
    integer, intent(out) :: info
    real(dp) :: f0,f1,f2,fa,norm
    coeff=0.0_dp
    info=0
    if (mk <= 0.0_dp .or. szone <= 0.0_dp) then
       info=-1
       return
    endif
    if (nf < 0) then
       info=-1
       return
    endif
    if (a == cs_parton_q .and. b == cs_parton_q) then
       norm=cs_cf_lc
    elseif (a == cs_parton_g .and. b == cs_parton_g) then
       norm=cs_ca
    elseif ((a == cs_parton_q .and. b == cs_parton_g) .or. &
         (a == cs_parton_g .and. b == cs_parton_q)) then
       return
    else
       info=-2
       return
    endif
    call if_endpoint_base(a,b,mk,szone,0.0_dp,scheme,f0)
    call if_endpoint_base(a,b,mk,szone,1.0_dp,scheme,f1)
    call if_endpoint_base(a,b,mk,szone,2.0_dp,scheme,f2)
    call if_endpoint_base(a,b,mk,szone,ell,scheme,fa)
    call laurent_from_finite(f0,f1,f2,fa,norm,coeff)
  end subroutine cs_massive_if_endpoint

  pure subroutine if_endpoint_base(a,b,mk,szone,ell,scheme,value)
    integer, intent(in) :: a,b,scheme
    real(dp), intent(in) :: mk,szone,ell
    real(dp), intent(out) :: value
    real(dp) :: mu2,rs
    mu2=mk*mk/szone
    rs=0.0_dp
    if (scheme == cs_scheme_fdh) then
       if (a == cs_parton_q) rs=-0.5_dp
       if (a == cs_parton_g) rs=-1.0_dp/6.0_dp
    endif
    if (a == cs_parton_q .and. b == cs_parton_q) then
       value=(2.0_dp*ell**2-cs_pi**2+4.0_dp*rs+&
            2.0_dp*log(1.0_dp+mu2)*(2.0_dp*ell+log(1.0_dp+mu2))+&
            8.0_dp*cs_dilog(1.0_dp/(1.0_dp+mu2)))/4.0_dp
    elseif (a == cs_parton_g .and. b == cs_parton_g) then
       value=(6.0_dp*ell**2-3.0_dp*cs_pi**2+12.0_dp*rs+&
            6.0_dp*log(1.0_dp+mu2)*(2.0_dp*ell+log(1.0_dp+mu2))+&
            24.0_dp*cs_dilog(1.0_dp/(1.0_dp+mu2)))/12.0_dp
    else
       value=0.0_dp
    endif
  end subroutine if_endpoint_base

  pure subroutine cs_massive_if_convolution(a,b,mk,s,szone,x,mu2_over_s,&
       mu2_over_szone,alpha,nf,kernel,info)
    integer, intent(in) :: a,b,nf
    real(dp), intent(in) :: mk,s,szone,x,mu2_over_s,mu2_over_szone,alpha
    type(cs_convolution_kernel), intent(out) :: kernel
    integer, intent(out) :: info
    real(dp) :: m2,m2one,zp,l,lone,omx,core

    kernel=cs_convolution_kernel()
    info=0
    if (mk <= 0.0_dp .or. s <= 0.0_dp .or. szone <= 0.0_dp) then
       info=-1
       return
    endif
    if (nf < 0) then
       info=-1
       return
    endif
    if (x <= 0.0_dp .or. x >= 1.0_dp) then
       info=-1
       return
    endif
    if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
       info=-3
       return
    endif
    m2=mk*mk/s
    m2one=mk*mk/szone
    omx=1.0_dp-x
    zp=omx/(omx+m2)
    l=log(mu2_over_s)
    lone=log(mu2_over_szone)

    if (a == cs_parton_q .and. b == cs_parton_q) then
       kernel%regular=-cs_cf_lc*(-((-1.0_dp+x*x)*l)+&
            (-1.0_dp+x*x)*log(omx)-2.0_dp*log(2.0_dp-x)+&
            (-omx)*(-omx+(1.0_dp+x)*log(omx/(omx+m2))))/(-omx)
       kernel%plus_z=cs_cf_lc*2.0_dp*(l-2.0_dp*log(omx))/(-omx)
       kernel%plus_one=cs_cf_lc*2.0_dp*(lone-2.0_dp*log(omx))/(-omx)
       kernel%plus_z=kernel%plus_z+cs_cf_lc*2.0_dp*&
            log((2.0_dp-x)/(2.0_dp-x+m2))/omx
       kernel%plus_one=kernel%plus_one+cs_cf_lc*2.0_dp*&
            log(1.0_dp/(1.0_dp+m2one))/omx
       if (zp > alpha) then
          kernel%regular=kernel%regular-cs_cf_lc*(-(1.0_dp+x)*log(zp/alpha)+&
               2.0_dp*log((1.0_dp+alpha-x)*zp/&
               (alpha*(1.0_dp-x+zp)))/omx)
       endif
    elseif (a == cs_parton_q .and. b == cs_parton_g) then
       core=x*x-l*(2.0_dp+(-2.0_dp+x)*x)+&
            (2.0_dp+(-2.0_dp+x)*x)*log(omx)+&
            (2.0_dp+(-2.0_dp+x)*x)*log(omx/(omx+m2))
       kernel%regular=cs_cf_initial_qg*core/x
       if (zp > alpha) then
          kernel%regular=kernel%regular-cs_cf_initial_qg*&
               (2.0_dp*m2*log((1.0_dp-zp)/(1.0_dp-alpha))/x+&
               (1.0_dp+omx*omx)*log(zp/alpha)/x)
       endif
    elseif (a == cs_parton_g .and. b == cs_parton_q) then
       core=-l+2.0_dp*(1.0_dp+l)*x-2.0_dp*(1.0_dp+l)*x*x+&
            (1.0_dp+2.0_dp*(-1.0_dp+x)*x)*log(omx)+&
            (1.0_dp+2.0_dp*(-1.0_dp+x)*x)*log(omx/(omx+m2))
       kernel%regular=cs_tr_initial_lc*core
       if (zp > alpha) kernel%regular=kernel%regular-cs_tr_initial_lc*&
            (omx*omx+x*x)*log(zp/alpha)
    elseif (a == cs_parton_g .and. b == cs_parton_g) then
       core=l-3.0_dp*l*x+3.0_dp*l*x*x-2.0_dp*l*x**3+l*x**4-&
            (-1.0_dp+x)*(-1.0_dp+x*(2.0_dp+(-1.0_dp+x)*x))*log(omx)+&
            x*log(2.0_dp-x)-m2*log(m2/(1.0_dp+m2-x))+&
            m2*x*log(m2/(1.0_dp+m2-x))-&
            (-1.0_dp+x)*(-1.0_dp+x*(2.0_dp+(-1.0_dp+x)*x))*&
            log(omx/(omx+m2))
       kernel%regular=2.0_dp*cs_ca*core/((-1.0_dp+x)*x)
       kernel%plus_z=2.0_dp*cs_ca*(l-2.0_dp*log(omx))/(-omx)
       kernel%plus_one=2.0_dp*cs_ca*(lone-2.0_dp*log(omx))/(-omx)
       kernel%plus_z=kernel%plus_z+2.0_dp*cs_ca*&
            log((2.0_dp-x)/(2.0_dp-x+m2))/omx
       kernel%plus_one=kernel%plus_one+2.0_dp*cs_ca*&
            log(1.0_dp/(1.0_dp+m2one))/omx
       if (zp > alpha) then
          kernel%regular=kernel%regular+cs_ca*(-2.0_dp*m2*&
               log((1.0_dp-zp)/(1.0_dp-alpha))/x-&
               2.0_dp*(-1.0_dp+omx/x+omx*x)*log(zp/alpha)+&
               2.0_dp*log(alpha*(1.0_dp-x+zp)/&
               ((1.0_dp+alpha-x)*zp))/omx)
       endif
    else
       info=-2
    endif
  end subroutine cs_massive_if_convolution

end module cs_massive_integrated_kernels
