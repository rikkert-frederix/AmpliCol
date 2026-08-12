module FeynmanRules
! Generic rules are named for their Lorentz fields.  Rules with explicit
! quark/gluon names below are the fixed-normalization QCD colour-flow kernels.
contains
  subroutine ext_massless_vector_cmplx(p,ihel,idum,wf)
    ! External massless-vector wavefunction. From HELAS.
    implicit none
    integer :: ihel,idum
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(4) :: wf
    real(kind=8),parameter :: rZero=0d0,sqh=sqrt(0.5d0)
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8) :: hel,pp,pt,pzpt
    if (p(0).eq.0d0) then
       write (*,*) 'Cannot generate an external massless vector with zero energy'
       write (*,*) p
       stop 1
    elseif (p(0).gt.0d0) then
       hel = dble(ihel)
       pp = p(0)
       pt = sqrt(p(1)**2+p(2)**2)
       wf(1) = cZero
       wf(4) = dcmplx( hel*pt/pp*sqh )
       if ( pt.ne.rZero ) then
          pzpt = p(3)/(pp*pt)*sqh*hel
          wf(2) = dcmplx( -p(1)*pzpt , -p(2)/pt*sqh )
          wf(3) = dcmplx( -p(2)*pzpt ,  p(1)/pt*sqh )
       else
          wf(2) = dcmplx( -hel*sqh )
          wf(3) = dcmplx( rZero , sign(sqh,p(3)) )
       endif
    else
       hel = dble(-ihel)
       pp = -p(0)
       pt = sqrt(p(1)**2+p(2)**2)
       wf(1) = cZero
       wf(4) = dcmplx( hel*pt/pp*sqh )
       if ( pt.ne.rZero ) then
          pzpt = -p(3)/(pp*pt)*sqh*hel
          wf(2) = dcmplx( p(1)*pzpt ,  p(2)/pt*sqh )
          wf(3) = dcmplx( p(2)*pzpt , -p(1)/pt*sqh )
       else
          wf(2) = dcmplx( -hel*sqh )
          wf(3) = dcmplx( rZero , -sign(sqh,p(3)) )
       endif
    endif
  end subroutine ext_massless_vector_cmplx


  subroutine ext_massive_vector(p,nhel,nsv,wf,vmass)
    implicit none
    real(kind=8),dimension(0:3) :: p(0:3)
    complex(kind=8),dimension(4) :: wf
    real(kind=8) :: vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
    integer nhel,nsv,nsvahl
    real(kind=8),parameter :: rZero=0d0, rHalf=0.5d0, rOne=1d0, rTwo=2d0
    sqh = dsqrt(rHalf)
    hel = dble(nhel)
    nsvahl = nsv*abs(nhel)
    pt2 = p(1)**2+p(2)**2
    pp = min(p(0),dsqrt(pt2+p(3)**2))
    pt = min(pp,dsqrt(pt2))
    if ( vmass.ne.rZero ) then
       hel0 = rOne-dabs(hel)
       if ( pp.eq.rZero ) then
          wf(1) = dcmplx( rZero )
          wf(2) = dcmplx(-hel*sqh )
          wf(3) = dcmplx( rZero , nsvahl*sqh )
          wf(4) = dcmplx( hel0 )
       else
          emp = p(0)/(vmass*pp)
          wf(1) = dcmplx( hel0*pp/vmass )
          wf(4) = dcmplx( hel0*p(3)*emp+hel*pt/pp*sqh )
          if ( pt.ne.rZero ) then
             pzpt = p(3)/(pp*pt)*sqh*hel
             wf(2) = dcmplx( hel0*p(1)*emp-p(1)*pzpt , -nsvahl*p(2)/pt*sqh       )
             wf(3) = dcmplx( hel0*p(2)*emp-p(2)*pzpt ,  nsvahl*p(1)/pt*sqh       )
          else
             wf(2) = dcmplx( -hel*sqh )
             wf(3) = dcmplx( rZero , nsvahl*sign(sqh,p(3)) )
          endif
       endif
    else
       pp = p(0)
       pt = sqrt(p(1)**2+p(2)**2)
       wf(1) = dcmplx( rZero )
       wf(4) = dcmplx( hel*pt/pp*sqh )
       if ( pt.ne.rZero ) then
          pzpt = p(3)/(pp*pt)*sqh*hel
          wf(2) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
          wf(3) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
       else
          wf(2) = dcmplx( -hel*sqh )
          wf(3) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
       endif
    endif
  end subroutine ext_massive_vector

  subroutine ext_fermion_outflow(p,nhel,idum,wf,fmass)
    ! Flowing-out fermion number: final-state fermion (p(0)>0) or
    ! initial-state antifermion (p(0)<0).
    implicit none
    integer :: nhel,idum
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(4) :: wf
    complex(kind=8), dimension(2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8) :: sqp0p3
    real(kind=8), parameter :: tiny=1d-8
    real(kind=8) :: fmass, pp,pp3,lim
    real(kind=8) :: omega(2),sfomeg(2),sf(2)
    integer :: im,ip,nsf,nh

    lim=tiny

    if (p(0).gt.0d0) then
       ! outgoing final state momenta
       if (abs(fmass).lt.lim) then
          if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
             sqp0p3 = 0d0
          else
             sqp0p3 = dsqrt(max(p(0)+p(3),rZero))
          end if
          chi(1) = dcmplx( sqp0p3 )
          if ( sqp0p3.eq.rZero ) then
             chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
          else
             chi(2) = dcmplx( nhel*p(1), -p(2) )/sqp0p3
          endif
          if ( nhel.eq.1 ) then
             wf(1) = chi(1)
             wf(2) = chi(2)
             wf(3) = cZero
             wf(4) = cZero
          else
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(2)
             wf(4) = chi(1)
          endif
       else
          nsf=+1
          nh=nsf*nhel
          pp = abs(dsqrt(p(1)**2+p(2)**2+p(3)**2))
          sf(1) = dble(1+nsf+(1-nsf)*nh)*0.5d0
          sf(2) = dble(1+nsf-(1-nsf)*nh)*0.5d0
          omega(1) = dsqrt(p(0)+pp)
          omega(2) = fmass/omega(1)
          ip = (3+nh)/2
          im = (3-nh)/2
          sfomeg(1) = sf(1)*omega(ip)
          sfomeg(2) = sf(2)*omega(im)
          pp3 = max(pp+p(3),rZero)
          chi(1) = dcmplx( dsqrt(pp3*0.5d0/pp) )
          if ( pp3.eq.rZero ) then
             chi(2) = dcmplx(-nh )
          else
             chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
          endif
          wf(1) = sfomeg(2)*chi(im)
          wf(2) = sfomeg(2)*chi(ip)
          wf(3) = sfomeg(1)*chi(im)
          wf(4) = sfomeg(1)*chi(ip)
       endif
    else
       ! "outgoing" initial state momenta
       if (abs(fmass).lt.lim) then
          if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).gt.0d0) then
             sqp0p3 = 0d0
          else
             sqp0p3 = -dsqrt(max(-(p(0)+p(3)),rZero))
          end if
          chi(1) = dcmplx( sqp0p3 )
          if ( sqp0p3.eq.rZero ) then
             chi(2) = dcmplx(-nhel )*dsqrt(rTwo*abs(p(0)))
          else
             chi(2) = dcmplx( -nhel*(-p(1)), -(-p(2)) )/sqp0p3
          endif
          if ( -nhel.eq.1 ) then
             wf(1) = chi(1)
             wf(2) = chi(2)
             wf(3) = cZero
             wf(4) = cZero
          else
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(2)
             wf(4) = chi(1)
          endif
       else
          nsf=-1
          nh=nsf*nhel
          pp = abs(dsqrt(p(1)**2+p(2)**2+p(3)**2))
          sf(1) = dble(1+nsf+(1-nsf)*nh)*0.5d0
          sf(2) = dble(1+nsf-(1-nsf)*nh)*0.5d0
          omega(1) = dsqrt(abs(p(0))+pp)
          omega(2) = fmass/omega(1)
          ip = (3+nh)/2
          im = (3-nh)/2
          sfomeg(1) = sf(1)*omega(ip)
          sfomeg(2) = sf(2)*omega(im)
          pp3 = max(pp+(-p(3)),rZero)
          chi(1) = dcmplx( dsqrt(pp3*0.5d0/pp) )
          if ( pp3.eq.rZero ) then
             chi(2) = dcmplx(-nh )
          else
             chi(2) = dcmplx( nh*(-p(1)) , -(-p(2)) )/dsqrt(rTwo*pp*pp3)
          endif
          wf(1) = sfomeg(2)*chi(im)
          wf(2) = sfomeg(2)*chi(ip)
          wf(3) = sfomeg(1)*chi(im)
          wf(4) = sfomeg(1)*chi(ip)
       endif
    endif
  end subroutine ext_fermion_outflow


  subroutine ext_fermion_inflow(p,nhel,idum,wf,fmass)
    ! Flowing-in fermion number: final-state antifermion (p(0)>0) or
    ! initial-state fermion (p(0)<0).
    implicit none
    integer :: nhel,idum
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(4) :: wf
    complex(kind=8), dimension(2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8) :: sqp0p3
    real(kind=8), parameter :: tiny=1d-8
    real(kind=8) :: fmass, pp,pp3,lim
    real(kind=8) :: omega(2),sfomeg(2),sf(2)
    integer :: im,ip,nsf,nh

    lim=tiny

    if(p(0).gt.0d0) then
       ! outgoing final state momenta
       if (abs(fmass).lt.lim) then
          if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
             sqp0p3 = 0d0
          else
             sqp0p3 = -dsqrt(max(p(0)+p(3),rZero))
          end if
          chi(1) = dcmplx( sqp0p3 )
          if ( sqp0p3.eq.rZero ) then
             chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
          else
             chi(2) = dcmplx(-nhel*p(1), p(2) )/sqp0p3
          endif
          if ( -nhel.eq.1 ) then
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(1)
             wf(4) = chi(2)
          else
             wf(1) = chi(2)
             wf(2) = chi(1)
             wf(3) = cZero
             wf(4) = cZero
          endif
       else
          nsf=-1
          nh=nsf*nhel
          pp = abs(dsqrt(p(1)**2+p(2)**2+p(3)**2))
          sf(1) = dble(1+nsf+(1-nsf)*nh)*0.5d0
          sf(2) = dble(1+nsf-(1-nsf)*nh)*0.5d0
          omega(1) = dsqrt(p(0)+pp)
          omega(2) = fmass/omega(1)
          ip = (3+nh)/2
          im = (3-nh)/2
          sfomeg(1) = sf(1)*omega(ip)
          sfomeg(2) = sf(2)*omega(im)
          pp3 = max(pp+p(3),rZero)
          chi(1) = dcmplx( dsqrt(pp3*0.5d0/pp) )
          if ( pp3.eq.rZero ) then
             chi(2) = dcmplx(-nh )
          else
             chi(2) = dcmplx( nh*p(1) , p(2) )/dsqrt(rTwo*pp*pp3)
          endif
          wf(1) = sfomeg(1)*chi(im)
          wf(2) = sfomeg(1)*chi(ip)
          wf(3) = sfomeg(2)*chi(im)
          wf(4) = sfomeg(2)*chi(ip)
       endif
    else
       ! "outgoing" initial state momenta
       if (abs(fmass).lt.lim) then
          if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).gt.0d0) then
             sqp0p3 = 0d0
          else
             sqp0p3 = dsqrt(max(-(p(0)+p(3)),rZero))
          end if
          chi(1) = dcmplx( sqp0p3 )
          if ( sqp0p3.eq.rZero ) then
             chi(2) = dcmplx( -nhel )*dsqrt(rTwo*abs(p(0)))
          else
             chi(2) = dcmplx( nhel*(-p(1)), (-p(2)) )/sqp0p3
          endif
          if ( nhel.eq.1 ) then
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(1)
             wf(4) = chi(2)
          else
             wf(1) = chi(2)
             wf(2) = chi(1)
             wf(3) = cZero
             wf(4) = cZero
          endif
       else
          nsf=+1
          nh=nsf*nhel
          pp = abs(dsqrt(p(1)**2+p(2)**2+p(3)**2))
          sf(1) = dble(1+nsf+(1-nsf)*nh)*0.5d0
          sf(2) = dble(1+nsf-(1-nsf)*nh)*0.5d0
          omega(1) = dsqrt(abs(p(0))+pp)
          omega(2) = fmass/omega(1)
          ip = (3+nh)/2
          im = (3-nh)/2
          sfomeg(1) = sf(1)*omega(ip)
          sfomeg(2) = sf(2)*omega(im)
          pp3 = max(pp+(-p(3)),rZero)
          chi(1) = dcmplx( dsqrt(pp3*0.5d0/pp) )
          if ( pp3.eq.rZero ) then
             chi(2) = dcmplx(-nh )
          else
             chi(2) = dcmplx( nh*(-p(1)) , (-p(2)) )/dsqrt(rTwo*pp*pp3)
          endif
          wf(1) = sfomeg(1)*chi(im)
          wf(2) = sfomeg(1)*chi(ip)
          wf(3) = sfomeg(2)*chi(im)
          wf(4) = sfomeg(2)*chi(ip)
       endif
    endif
  end subroutine ext_fermion_inflow

  subroutine ext_fermion_outflow_weyl(p,nhel,idum,wf,chirality)
    implicit none
    integer :: nhel,idum,chirality
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(2) :: wf
    complex(kind=8), dimension(2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    real(kind=8) :: sqp0p3

    wf(1:2)=(0d0,0d0)
    if (p(0).gt.0d0) then
       if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
          sqp0p3 = 0d0
       else
          sqp0p3 = dsqrt(max(p(0)+p(3),rZero))
       end if
       chi(1) = dcmplx( sqp0p3 )
       if ( sqp0p3.eq.rZero ) then
          chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
       else
          chi(2) = dcmplx( nhel*p(1), -p(2) )/sqp0p3
       endif
       if (nhel.eq.1 .and. chirality.eq.1) then
          wf(1:2)=chi(1:2)
       elseif (nhel.eq.-1 .and. chirality.eq.-1) then
          wf(1)=chi(2)
          wf(2)=chi(1)
       endif
    else
       if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).gt.0d0) then
          sqp0p3 = 0d0
       else
          sqp0p3 = -dsqrt(max(-(p(0)+p(3)),rZero))
       end if
       chi(1) = dcmplx( sqp0p3 )
       if ( sqp0p3.eq.rZero ) then
          chi(2) = dcmplx(-nhel )*dsqrt(rTwo*abs(p(0)))
       else
          chi(2) = dcmplx( -nhel*(-p(1)), -(-p(2)) )/sqp0p3
       endif
       if (nhel.eq.-1 .and. chirality.eq.1) then
          wf(1:2)=chi(1:2)
       elseif (nhel.eq.1 .and. chirality.eq.-1) then
          wf(1)=chi(2)
          wf(2)=chi(1)
       endif
    endif
  end subroutine ext_fermion_outflow_weyl

  subroutine ext_fermion_inflow_weyl(p,nhel,idum,wf,chirality)
    implicit none
    integer :: nhel,idum,chirality
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(2) :: wf
    complex(kind=8), dimension(2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    real(kind=8) :: sqp0p3

    wf(1:2)=(0d0,0d0)
    if(p(0).gt.0d0) then
       if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
          sqp0p3 = 0d0
       else
          sqp0p3 = -dsqrt(max(p(0)+p(3),rZero))
       end if
       chi(1) = dcmplx( sqp0p3 )
       if ( sqp0p3.eq.rZero ) then
          chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
       else
          chi(2) = dcmplx(-nhel*p(1), p(2) )/sqp0p3
       endif
       if (nhel.eq.1 .and. chirality.eq.1) then
          wf(1)=chi(2)
          wf(2)=chi(1)
       elseif (nhel.eq.-1 .and. chirality.eq.-1) then
          wf(1:2)=chi(1:2)
       endif
    else
       if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).gt.0d0) then
          sqp0p3 = 0d0
       else
          sqp0p3 = dsqrt(max(-(p(0)+p(3)),rZero))
       end if
       chi(1) = dcmplx( sqp0p3 )
       if ( sqp0p3.eq.rZero ) then
          chi(2) = dcmplx( -nhel )*dsqrt(rTwo*abs(p(0)))
       else
          chi(2) = dcmplx( nhel*(-p(1)), (-p(2)) )/sqp0p3
       endif
       if (nhel.eq.-1 .and. chirality.eq.1) then
          wf(1)=chi(2)
          wf(2)=chi(1)
       elseif (nhel.eq.1 .and. chirality.eq.-1) then
          wf(1:2)=chi(1:2)
       endif
    endif
  end subroutine ext_fermion_inflow_weyl

  subroutine ext_scalar(p,idum,wf)
    implicit none
    integer :: idum
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(1) :: wf
    real(kind=8),parameter :: rOne=1d0
    wf(1)=(rOne,0d0)
  end subroutine ext_scalar

  subroutine ThreeGluon(wf1,pwf1,wf2,pwf2,wf)
    ! Colour-ordered three-gluon interaction
    implicit none
    complex(kind=8),dimension(4) :: wf1,wf2,wf
    real(kind=8),dimension(0:3) :: pwf1,pwf2
    complex(kind=8),parameter :: prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3
    TMP1 = (wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4))
    TMP2 = (wf1(1)*pwf2(0)-wf1(2)*pwf2(1)-wf1(3)*pwf2(2)-wf1(4)*pwf2(3))
    TMP3 = (wf2(1)*pwf1(0)-wf2(2)*pwf1(1)-wf2(3)*pwf1(2)-wf2(4)*pwf1(3))
    wf(1:4) = prefact*(TMP1*(pwf1(0:3)-pwf2(0:3))+2d0*(TMP2*wf2(1:4)-TMP3*wf1(1:4)))
  end subroutine ThreeGluon
  subroutine FourGluon(wf1,wf2,wf3,wf)
    ! Colour-ordered four-gluon interaction
    implicit none
    complex(kind=8),dimension(4) :: wf1,wf2,wf3,wf
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    complex(kind=8) :: TMP1,TMP2,TMP3
    TMP1 = (wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4))
    TMP2 = (wf1(1)*wf3(1)-wf1(2)*wf3(2)-wf1(3)*wf3(3)-wf1(4)*wf3(4))
    TMP3 = (wf2(1)*wf3(1)-wf2(2)*wf3(2)-wf2(3)*wf3(3)-wf2(4)*wf3(4))
    wf(1:4) = prefact*(2d0*wf2(1:4)*TMP2-wf1(1:4)*TMP3-wf3(1:4)*TMP1)
  end subroutine FourGluon
  subroutine TwoGluonToAuxTensor(wfg1,wfg2,wfT)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfg2
    complex(kind=8),dimension(6) :: wfT
    complex(kind=8),parameter :: prefact=(0d0,1d0)
    ! Since it is an anti-symmetric 4x4 tensor, take only the upper-right triangle.
    wfT(1)=(wfg1(1)*wfg2(2)-wfg1(2)*wfg2(1))! * prefact
    wfT(2)=(wfg1(1)*wfg2(3)-wfg1(3)*wfg2(1))! * prefact
    wfT(3)=(wfg1(1)*wfg2(4)-wfg1(4)*wfg2(1))! * prefact
    wfT(4)=(wfg1(2)*wfg2(3)-wfg1(3)*wfg2(2))! * prefact
    wfT(5)=(wfg1(2)*wfg2(4)-wfg1(4)*wfg2(2))! * prefact
    wfT(6)=(wfg1(3)*wfg2(4)-wfg1(4)*wfg2(3))! * prefact
  end subroutine TwoGluonToAuxTensor
  subroutine AuxTensorGluonToGluon(wfT1,wfg2,wfg)
    implicit none
    complex(kind=8),dimension(4) :: wfg2,wfg
    complex(kind=8),dimension(6) :: wfT1
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact
  end subroutine AuxTensorGluonToGluon
  subroutine GluonAuxTensorToGluon(wfg1,wfT2,wfg)
    implicit none 
    complex(kind=8),dimension(4) :: wfg1,wfg
    complex(kind=8),dimension(6) :: wfT2
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact
  end subroutine GluonAuxTensorToGluon
  subroutine QuarkColourFlowVectorToQuark(wfq1,wfg2,wfq)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfg2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    wfq(1)=prefact*(TMP1*wfq1(3)+TMP3*wfq1(4)) 
    wfq(2)=prefact*(TMP2*wfq1(4)+TMP4*wfq1(3)) 
    wfq(3)=prefact*(TMP2*wfq1(1)-TMP3*wfq1(2))  
    wfq(4)=prefact*(TMP1*wfq1(2)-TMP4*wfq1(1))  
  end subroutine QuarkColourFlowVectorToQuark

  subroutine ColourFlowVectorQuarkToQuark(wfg1,wfq2,wfq)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    wfq(1)=prefact*(TMP1*wfq2(3)+TMP3*wfq2(4)) ! sl1 ! minus sign
    wfq(2)=prefact*(TMP2*wfq2(4)+TMP4*wfq2(3)) ! sl2
    wfq(3)=prefact*(TMP2*wfq2(1)-TMP3*wfq2(2)) ! sr1
    wfq(4)=prefact*(TMP1*wfq2(2)-TMP4*wfq2(1)) ! sr2
  end subroutine ColourFlowVectorQuarkToQuark

  subroutine VectorFermionToFermion(wfg1,wfq2,wfq,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    ! L
    wfq(1)=prefact*(TMP1*wfq2(3)+TMP3*wfq2(4))*coupl(1) ! sl1 ! minus sign
    wfq(2)=prefact*(TMP2*wfq2(4)+TMP4*wfq2(3))*coupl(1) ! sl2
    ! R
    wfq(3)=prefact*(TMP2*wfq2(1)-TMP3*wfq2(2))*coupl(2) ! sr1
    wfq(4)=prefact*(TMP1*wfq2(2)-TMP4*wfq2(1))*coupl(2) ! sr2
  end subroutine VectorFermionToFermion

  subroutine VectorAntifermionToAntifermion(wfg1,wfq2,wfq,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    ! L
    wfq(1)=prefact*(TMP2*wfq2(3)-TMP4*wfq2(4))*coupl(2) !sr1
    wfq(2)=prefact*(TMP1*wfq2(4)-TMP3*wfq2(3))*coupl(2) !sr2
    ! R
    wfq(3)=prefact*(TMP1*wfq2(1)+TMP4*wfq2(2))*coupl(1) !sl1
    wfq(4)=prefact*(TMP2*wfq2(2)+TMP3*wfq2(1))*coupl(1) !sl2
  end subroutine VectorAntifermionToAntifermion

  subroutine AntiquarkColourFlowVectorToAntiquark(wfq1,wfg2,wfq)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfg2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    wfq(1)=prefact*(TMP2*wfq1(3)-TMP4*wfq1(4)) !sr1
    wfq(2)=prefact*(TMP1*wfq1(4)-TMP3*wfq1(3)) !sr2
    wfq(3)=prefact*(TMP1*wfq1(1)+TMP4*wfq1(2)) !sl1
    wfq(4)=prefact*(TMP2*wfq1(2)+TMP3*wfq1(1)) !sl2
  end subroutine AntiquarkColourFlowVectorToAntiquark

  subroutine ColourFlowVectorAntiquarkToAntiquark(wfg1,wfq2,wfq)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    wfq(1)=prefact*(TMP2*wfq2(3)-TMP4*wfq2(4)) !sr1 
    wfq(2)=prefact*(TMP1*wfq2(4)-TMP3*wfq2(3)) !sr2
    wfq(3)=prefact*(TMP1*wfq2(1)+TMP4*wfq2(2)) !sl1
    wfq(4)=prefact*(TMP2*wfq2(2)+TMP3*wfq2(1)) !sl2
  end subroutine ColourFlowVectorAntiquarkToAntiquark

  subroutine QuarkAntiquarkToColourFlowU1Vector(wfq1,wfq2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq2,wfg
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    TMP1=wfq1(3)*wfq2(1)+wfq1(2)*wfq2(4)
    TMP2=wfq1(4)*wfq2(2)+wfq1(1)*wfq2(3)
    TMP3=wfq1(2)*wfq2(3)-wfq1(4)*wfq2(1)
    TMP4=wfq1(1)*wfq2(4)-wfq1(3)*wfq2(2)
    wfg(1)=( TMP1 + TMP2 )*prefact*coupl(1)
    wfg(2)=( TMP4 + TMP3 )*prefact*coupl(1)
    wfg(3)=(-TMP4 + TMP3 )*cImag*prefact*coupl(1)
    wfg(4)=(-TMP1 + TMP2 )*prefact*coupl(1)
  end subroutine QuarkAntiquarkToColourFlowU1Vector
  subroutine AntiquarkQuarkToGluon(wfq1,wfq2,wfg)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq2,wfg
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfq2(3)*wfq1(1)+wfq2(2)*wfq1(4)
    TMP2=wfq2(4)*wfq1(2)+wfq2(1)*wfq1(3)
    TMP3=wfq2(2)*wfq1(3)-wfq2(4)*wfq1(1)
    TMP4=wfq2(1)*wfq1(4)-wfq2(3)*wfq1(2)
    wfg(1)=( TMP1 + TMP2 )*prefact
    wfg(2)=( TMP4 + TMP3 )*prefact
    wfg(3)=(-TMP4 + TMP3 )*cImag*prefact
    wfg(4)=(-TMP1 + TMP2 )*prefact
  end subroutine AntiquarkQuarkToGluon

  subroutine VectorVectorToVector(wf1,pwf1,wf2,pwf2,wf,coupl)
    ! Yang-Mills three-vector interaction with an explicit model coupling.
    implicit none
    complex(kind=8),dimension(4) :: wf1,wf2,wf
    real(kind=8),dimension(0:3) :: pwf1,pwf2
    complex(kind=8),parameter :: prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    TMP1 = wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4)
    TMP2 = wf1(1)*(2d0*pwf2(0)+pwf1(0))-wf1(2)*(2d0*pwf2(1)+pwf1(1))-wf1(3)*(2d0*pwf2(2)+pwf1(2))-wf1(4)*(2d0*pwf2(3)+pwf1(3))
    TMP3 = wf2(1)*(-2d0*pwf1(0)-pwf2(0))-wf2(2)*(-2d0*pwf1(1)-pwf2(1))-wf2(3)*(-2d0*pwf1(2)-pwf2(2))-wf2(4)*(-2d0*pwf1(3)-pwf2(3))
    TMP4 = prefact*coupl(1)
    wf(1:4) = TMP4*(TMP1*(pwf1(0:3)-pwf2(0:3))+TMP2*wf2(1:4)+TMP3*wf1(1:4))
  end subroutine VectorVectorToVector

  subroutine FermionVectorToFermion(wfq1,wfg2,wfq,coupl) ! from fvoxxx.f
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfg2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    TMP5=prefact*coupl(1) ! L
    wfq(1)=TMP5*(TMP1*wfq1(3)+TMP3*wfq1(4))
    wfq(2)=TMP5*(TMP2*wfq1(4)+TMP4*wfq1(3))
    TMP5=prefact*coupl(2) ! R
    wfq(3)=TMP5*(TMP2*wfq1(1)-TMP3*wfq1(2))
    wfq(4)=TMP5*(TMP1*wfq1(2)-TMP4*wfq1(1))
  end subroutine FermionVectorToFermion

  subroutine QuarkColourFlowVectorToQuark_weyl(wfq1,wfg2,wfq,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq1,wfq
    complex(kind=8),dimension(4) :: wfg2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    if (chirality.eq.1) then
       wfq(1)=prefact*(TMP2*wfq1(1)-TMP3*wfq1(2))
       wfq(2)=prefact*(TMP1*wfq1(2)-TMP4*wfq1(1))
    elseif (chirality.eq.-1) then
       wfq(1)=prefact*(TMP1*wfq1(1)+TMP3*wfq1(2))
       wfq(2)=prefact*(TMP2*wfq1(2)+TMP4*wfq1(1))
    else
       call QuarkColourFlowVectorToQuark(wfq1,wfg2,wfq)
    endif
  end subroutine QuarkColourFlowVectorToQuark_weyl

  subroutine ColourFlowVectorQuarkToQuark_weyl(wfg1,wfq2,wfq,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4) :: wfg1
    complex(kind=8),dimension(*) :: wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    if (chirality.eq.1) then
       wfq(1)=prefact*(TMP2*wfq2(1)-TMP3*wfq2(2))
       wfq(2)=prefact*(TMP1*wfq2(2)-TMP4*wfq2(1))
    elseif (chirality.eq.-1) then
       wfq(1)=prefact*(TMP1*wfq2(1)+TMP3*wfq2(2))
       wfq(2)=prefact*(TMP2*wfq2(2)+TMP4*wfq2(1))
    else
       call ColourFlowVectorQuarkToQuark(wfg1,wfq2,wfq)
    endif
  end subroutine ColourFlowVectorQuarkToQuark_weyl

  subroutine FermionVectorToFermion_weyl(wfq1,wfg2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq1,wfq
    complex(kind=8),dimension(4) :: wfg2
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    if (chirality.eq.1) then
       TMP5=prefact*coupl(2)
       wfq(1)=TMP5*(TMP2*wfq1(1)-TMP3*wfq1(2))
       wfq(2)=TMP5*(TMP1*wfq1(2)-TMP4*wfq1(1))
    elseif (chirality.eq.-1) then
       TMP5=prefact*coupl(1)
       wfq(1)=TMP5*(TMP1*wfq1(1)+TMP3*wfq1(2))
       wfq(2)=TMP5*(TMP2*wfq1(2)+TMP4*wfq1(1))
    else
       call FermionVectorToFermion(wfq1,wfg2,wfq,coupl)
    endif
  end subroutine FermionVectorToFermion_weyl

  subroutine VectorFermionToFermion_weyl(wfg1,wfq2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4) :: wfg1
    complex(kind=8),dimension(*) :: wfq2,wfq
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    if (chirality.eq.1) then
       TMP5=prefact*coupl(2)
       wfq(1)=TMP5*(TMP2*wfq2(1)-TMP3*wfq2(2))
       wfq(2)=TMP5*(TMP1*wfq2(2)-TMP4*wfq2(1))
    elseif (chirality.eq.-1) then
       TMP5=prefact*coupl(1)
       wfq(1)=TMP5*(TMP1*wfq2(1)+TMP3*wfq2(2))
       wfq(2)=TMP5*(TMP2*wfq2(2)+TMP4*wfq2(1))
    else
       call VectorFermionToFermion(wfg1,wfq2,wfq,coupl)
    endif
  end subroutine VectorFermionToFermion_weyl

  subroutine AntifermionVectorToAntifermion(wfq1,wfg2,wfq,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfg2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    real(kind=8),dimension(2) :: coupl
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    TMP5=prefact*coupl(2)
    wfq(1)=TMP5*(TMP2*wfq1(3)-TMP4*wfq1(4)) !sr1
    wfq(2)=TMP5*(TMP1*wfq1(4)-TMP3*wfq1(3)) !sr2
    TMP5=prefact*coupl(1)
    wfq(3)=TMP5*(TMP1*wfq1(1)+TMP4*wfq1(2)) !sl1
    wfq(4)=TMP5*(TMP2*wfq1(2)+TMP3*wfq1(1)) !sl2
  end subroutine AntifermionVectorToAntifermion

  subroutine AntiquarkColourFlowVectorToAntiquark_weyl(wfq1,wfg2,wfq,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq1,wfq
    complex(kind=8),dimension(4) :: wfg2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    if (chirality.eq.1) then
       wfq(1)=prefact*(TMP1*wfq1(1)+TMP4*wfq1(2))
       wfq(2)=prefact*(TMP2*wfq1(2)+TMP3*wfq1(1))
    elseif (chirality.eq.-1) then
       wfq(1)=prefact*(TMP2*wfq1(1)-TMP4*wfq1(2))
       wfq(2)=prefact*(TMP1*wfq1(2)-TMP3*wfq1(1))
    else
       call AntiquarkColourFlowVectorToAntiquark(wfq1,wfg2,wfq)
    endif
  end subroutine AntiquarkColourFlowVectorToAntiquark_weyl

  subroutine ColourFlowVectorAntiquarkToAntiquark_weyl(wfg1,wfq2,wfq,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4) :: wfg1
    complex(kind=8),dimension(*) :: wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    if (chirality.eq.1) then
       wfq(1)=prefact*(TMP1*wfq2(1)+TMP4*wfq2(2))
       wfq(2)=prefact*(TMP2*wfq2(2)+TMP3*wfq2(1))
    elseif (chirality.eq.-1) then
       wfq(1)=prefact*(TMP2*wfq2(1)-TMP4*wfq2(2))
       wfq(2)=prefact*(TMP1*wfq2(2)-TMP3*wfq2(1))
    else
       call ColourFlowVectorAntiquarkToAntiquark(wfg1,wfq2,wfq)
    endif
  end subroutine ColourFlowVectorAntiquarkToAntiquark_weyl

  subroutine AntifermionVectorToAntifermion_weyl(wfq1,wfg2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq1,wfq
    complex(kind=8),dimension(4) :: wfg2
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=wfg2(2)+cImag*wfg2(3)
    TMP4=wfg2(2)-cImag*wfg2(3)
    if (chirality.eq.1) then
       TMP5=prefact*coupl(1)
       wfq(1)=TMP5*(TMP1*wfq1(1)+TMP4*wfq1(2))
       wfq(2)=TMP5*(TMP2*wfq1(2)+TMP3*wfq1(1))
    elseif (chirality.eq.-1) then
       TMP5=prefact*coupl(2)
       wfq(1)=TMP5*(TMP2*wfq1(1)-TMP4*wfq1(2))
       wfq(2)=TMP5*(TMP1*wfq1(2)-TMP3*wfq1(1))
    else
       call AntifermionVectorToAntifermion(wfq1,wfg2,wfq,coupl)
    endif
  end subroutine AntifermionVectorToAntifermion_weyl

  subroutine VectorAntifermionToAntifermion_weyl(wfg1,wfq2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4) :: wfg1
    complex(kind=8),dimension(*) :: wfq2,wfq
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4,TMP5
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    if (chirality.eq.1) then
       TMP5=prefact*coupl(1)
       wfq(1)=TMP5*(TMP1*wfq2(1)+TMP4*wfq2(2))
       wfq(2)=TMP5*(TMP2*wfq2(2)+TMP3*wfq2(1))
    elseif (chirality.eq.-1) then
       TMP5=prefact*coupl(2)
       wfq(1)=TMP5*(TMP2*wfq2(1)-TMP4*wfq2(2))
       wfq(2)=TMP5*(TMP1*wfq2(2)-TMP3*wfq2(1))
    else
       call VectorAntifermionToAntifermion(wfg1,wfq2,wfq,coupl)
    endif
  end subroutine VectorAntifermionToAntifermion_weyl

  subroutine WeylToDirac(wf_weyl,chirality,wf_dirac)
    ! Embed a compact massless-fermion current in the four-component
    ! representation used by the massive-fermion Feynman rules.  This is
    ! needed when a flavour-changing vector vertex connects a massless and
    ! a massive fermion, for example b-W-t.
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*),intent(in) :: wf_weyl
    complex(kind=8),dimension(4),intent(out) :: wf_dirac

    wf_dirac=(0d0,0d0)
    if (chirality.eq.1) then
       wf_dirac(1:2)=wf_weyl(1:2)
    elseif (chirality.eq.-1) then
       wf_dirac(3:4)=wf_weyl(1:2)
    elseif (chirality.eq.0) then
       wf_dirac(1:4)=wf_weyl(1:4)
    else
       write (*,*) 'Invalid Weyl chirality in WeylToDirac',chirality
       stop 1
    endif
  end subroutine WeylToDirac

  subroutine FermionVectorToFermion_mixed(wfq1,wfg2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*),intent(in) :: wfq1
    complex(kind=8),dimension(4),intent(in) :: wfg2
    complex(kind=8),dimension(4),intent(out) :: wfq
    real(kind=8),dimension(2),intent(in) :: coupl
    complex(kind=8),dimension(4) :: wfq1_dirac

    call WeylToDirac(wfq1,chirality,wfq1_dirac)
    call FermionVectorToFermion(wfq1_dirac,wfg2,wfq,coupl)
  end subroutine FermionVectorToFermion_mixed

  subroutine AntifermionVectorToAntifermion_mixed(wfq1,wfg2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*),intent(in) :: wfq1
    complex(kind=8),dimension(4),intent(in) :: wfg2
    complex(kind=8),dimension(4),intent(out) :: wfq
    real(kind=8),dimension(2),intent(in) :: coupl
    complex(kind=8),dimension(4) :: wfq1_dirac

    call WeylToDirac(wfq1,chirality,wfq1_dirac)
    call AntifermionVectorToAntifermion(wfq1_dirac,wfg2,wfq,coupl)
  end subroutine AntifermionVectorToAntifermion_mixed

  subroutine VectorFermionToFermion_mixed(wfg1,wfq2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4),intent(in) :: wfg1
    complex(kind=8),dimension(*),intent(in) :: wfq2
    complex(kind=8),dimension(4),intent(out) :: wfq
    real(kind=8),dimension(2),intent(in) :: coupl
    complex(kind=8),dimension(4) :: wfq2_dirac

    call WeylToDirac(wfq2,chirality,wfq2_dirac)
    call VectorFermionToFermion(wfg1,wfq2_dirac,wfq,coupl)
  end subroutine VectorFermionToFermion_mixed

  subroutine VectorAntifermionToAntifermion_mixed(wfg1,wfq2,wfq,coupl,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(4),intent(in) :: wfg1
    complex(kind=8),dimension(*),intent(in) :: wfq2
    complex(kind=8),dimension(4),intent(out) :: wfq
    real(kind=8),dimension(2),intent(in) :: coupl
    complex(kind=8),dimension(4) :: wfq2_dirac

    call WeylToDirac(wfq2,chirality,wfq2_dirac)
    call VectorAntifermionToAntifermion(wfg1,wfq2_dirac,wfq,coupl)
  end subroutine VectorAntifermionToAntifermion_mixed

  complex(kind=8) function ContractFermionCurrents(wf_terminal,chirality_terminal, &
       wf_external,chirality_external)
    ! The terminal current has not received a fermion propagator.  Its compact
    ! Weyl representation therefore occupies the opposite Dirac block from a
    ! propagated/external compact current of the same chirality.
    implicit none
    integer,intent(in) :: chirality_terminal,chirality_external
    complex(kind=8),dimension(*),intent(in) :: wf_terminal,wf_external

    if (chirality_terminal.ne.0 .and. chirality_external.ne.0) then
       ContractFermionCurrents=sum(wf_terminal(1:2)*wf_external(1:2))
    elseif (chirality_terminal.eq.0 .and. chirality_external.eq.0) then
       ContractFermionCurrents=sum(wf_terminal(1:4)*wf_external(1:4))
    elseif (chirality_terminal.eq.0) then
       if (chirality_external.eq.1) then
          ContractFermionCurrents=sum(wf_terminal(1:2)*wf_external(1:2))
       elseif (chirality_external.eq.-1) then
          ContractFermionCurrents=sum(wf_terminal(3:4)*wf_external(1:2))
       else
          write (*,*) 'Invalid external chirality in ContractFermionCurrents',chirality_external
          stop 1
       endif
    else
       if (chirality_terminal.eq.1) then
          ContractFermionCurrents=sum(wf_terminal(1:2)*wf_external(3:4))
       elseif (chirality_terminal.eq.-1) then
          ContractFermionCurrents=sum(wf_terminal(1:2)*wf_external(1:2))
       else
          write (*,*) 'Invalid terminal chirality in ContractFermionCurrents',chirality_terminal
          stop 1
       endif
    endif
  end function ContractFermionCurrents

  subroutine VectorVectorToAuxTensor(wfg1,wfg2,wfT,coupl)
    ! This vertex includes all factors needed to make the non-propagating
    ! auxiliary-tensor connection the identity.
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfg2
    complex(kind=8),dimension(6) :: wfT
    real(kind=8),dimension(2) :: coupl
    ! Since it is an anti-symmetric 4x4 tensor, take only the upper-right triangle.
    wfT(1)=(wfg1(1)*wfg2(2)-wfg1(2)*wfg2(1))*coupl(1)
    wfT(2)=(wfg1(1)*wfg2(3)-wfg1(3)*wfg2(1))*coupl(1)
    wfT(3)=(wfg1(1)*wfg2(4)-wfg1(4)*wfg2(1))*coupl(1)
    wfT(4)=(wfg1(2)*wfg2(3)-wfg1(3)*wfg2(2))*coupl(1)
    wfT(5)=(wfg1(2)*wfg2(4)-wfg1(4)*wfg2(2))*coupl(1)
    wfT(6)=(wfg1(3)*wfg2(4)-wfg1(4)*wfg2(3))*coupl(1)
  end subroutine VectorVectorToAuxTensor

  subroutine AuxTensorVectorToVector(wfT1,wfg2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg2,wfg
    complex(kind=8),dimension(6) :: wfT1
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact*coupl(1)
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact*coupl(1)
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact*coupl(1)
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact*coupl(1)
  end subroutine AuxTensorVectorToVector

  subroutine VectorAuxTensorToVector(wfg1,wfT2,wfg,coupl)
    implicit none 
    complex(kind=8),dimension(4) :: wfg1,wfg
    complex(kind=8),dimension(6) :: wfT2
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact*coupl(1)
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact*coupl(1)
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact*coupl(1)
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact*coupl(1)
  end subroutine VectorAuxTensorToVector


  subroutine FermionScalarToFermion(wfq1,wfs2,wfq,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq
    complex(kind=8),dimension(1) :: wfs2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    wfq(1:4)=-prefact*coupl(1)*wfs2(1)*wfq1(1:4)
  end subroutine FermionScalarToFermion

  subroutine VectorVectorToScalar(wfg1,wfg2,wfs,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfg2
    complex(kind=8),dimension(1) :: wfs
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    complex(kind=8) :: TMP
    TMP = wfg1(1)*wfg2(1)-wfg1(2)*wfg2(2)-wfg1(3)*wfg2(3)-wfg1(4)*wfg2(4)
    wfs(1)= prefact*coupl(1)*TMP
  end subroutine VectorVectorToScalar

  subroutine ScalarVectorToVector(wfs1,wfg2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg2,wfg
    complex(kind=8),dimension(1) :: wfs1
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1:4)= prefact*coupl(1)*wfs1(1)*(wfg2(1:4))
  end subroutine ScalarVectorToVector

  subroutine VectorScalarToVector(wfg1,wfs2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfg
    complex(kind=8),dimension(1) :: wfs2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1:4)= prefact*coupl(1)*wfs2(1)*(wfg1(1:4))
  end subroutine VectorScalarToVector

  subroutine ScalarScalarToScalar(wfs1,wfs2,wfs,coupl)
    implicit none
    complex(kind=8),dimension(1) :: wfs1,wfs2,wfs
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(2) :: coupl
    complex(kind=8) :: TMP
    TMP=(1d0,0d0)
    if (coupl(2).eq.-10d0) TMP=(0d0,1d0)
    wfs(1)= prefact*TMP*coupl(1)*wfs1(1)*wfs2(1)
  end subroutine ScalarScalarToScalar

  subroutine LeptonAntileptonToVector(wfq1,wfq2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq2,wfg,wfg_temp
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    ! L
    wfg_temp(1)=( wfq1(3)*wfq2(1)+wfq1(4)*wfq2(2) )*prefact*coupl(1)
    wfg_temp(2)=(-wfq1(4)*wfq2(1)-wfq1(3)*wfq2(2) )*prefact*coupl(1)
    wfg_temp(3)=(-wfq1(4)*wfq2(1)+wfq1(3)*wfq2(2) )*cImag*prefact*coupl(1)
    wfg_temp(4)=(-wfq1(3)*wfq2(1)+wfq1(4)*wfq2(2) )*prefact*coupl(1)
    ! R
    wfg(1)=( wfq1(1)*wfq2(3)+wfq1(2)*wfq2(4) )*prefact*coupl(2)
    wfg(2)=( wfq1(1)*wfq2(4)+wfq1(2)*wfq2(3) )*prefact*coupl(2)
    wfg(3)=(-wfq1(1)*wfq2(4)+wfq1(2)*wfq2(3) )*cImag*prefact*coupl(2)
    wfg(4)=( wfq1(1)*wfq2(3)-wfq1(2)*wfq2(4) )*prefact*coupl(2)
    ! add
    wfg(1:4)=wfg(1:4)+wfg_temp(1:4)
  end subroutine LeptonAntileptonToVector

  subroutine AntileptonLeptonToVector(wfq1,wfq2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq2,wfg,wfg_temp
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    real(kind=8),dimension(2) :: coupl
    ! L
    wfg_temp(1)=( wfq2(3)*wfq1(1)+wfq2(4)*wfq1(2) )*prefact*coupl(1)
    wfg_temp(2)=(-wfq2(4)*wfq1(1)-wfq2(3)*wfq1(2) )*prefact*coupl(1)
    wfg_temp(3)=(-wfq2(4)*wfq1(1)+wfq2(3)*wfq1(2) )*cImag*prefact*coupl(1)
    wfg_temp(4)=(-wfq2(3)*wfq1(1)+wfq2(4)*wfq1(2) )*prefact*coupl(1)
    ! R
    wfg(1)=( wfq2(1)*wfq1(3)+wfq2(2)*wfq1(4) )*prefact*coupl(2)
    wfg(2)=( wfq2(1)*wfq1(4)+wfq2(2)*wfq1(3) )*prefact*coupl(2)
    wfg(3)=(-wfq2(1)*wfq1(4)+wfq2(2)*wfq1(3) )*cImag*prefact*coupl(2)
    wfg(4)=( wfq2(1)*wfq1(3)-wfq2(2)*wfq1(4) )*prefact*coupl(2)
    ! add
    wfg(1:4)=wfg(1:4)+wfg_temp(1:4)
  end subroutine AntileptonLeptonToVector

  subroutine QuarkAntiquarkToColourFlowU1Vector_weyl(wfq1,wfq2,wfg,coupl,chirality1,chirality2)
    implicit none
    integer,intent(in) :: chirality1,chirality2
    complex(kind=8),dimension(*) :: wfq1,wfq2
    complex(kind=8),dimension(4) :: wfg
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: q1,q2,q3,q4,a1,a2,a3,a4,TMP1,TMP2,TMP3,TMP4,TMP5
    TMP5=prefact*coupl(1)
    wfg(1:4)=(0d0,0d0)
    if (chirality1.eq.-1 .and. chirality2.eq.1) then
       wfg(1)=TMP5*(wfq1(1)*wfq2(1)+wfq1(2)*wfq2(2))
       wfg(2)=-TMP5*(wfq1(1)*wfq2(2)+wfq1(2)*wfq2(1))
       wfg(3)=cImag*TMP5*(wfq1(1)*wfq2(2)-wfq1(2)*wfq2(1))
       wfg(4)=TMP5*(-wfq1(1)*wfq2(1)+wfq1(2)*wfq2(2))
       return
    elseif (chirality1.eq.1 .and. chirality2.eq.-1) then
       wfg(1)=TMP5*(wfq1(2)*wfq2(2)+wfq1(1)*wfq2(1))
       wfg(2)=TMP5*(wfq1(1)*wfq2(2)+wfq1(2)*wfq2(1))
       wfg(3)=cImag*TMP5*(-wfq1(1)*wfq2(2)+wfq1(2)*wfq2(1))
       wfg(4)=TMP5*(-wfq1(2)*wfq2(2)+wfq1(1)*wfq2(1))
       return
    elseif (chirality1.ne.0 .and. chirality2.ne.0) then
       return
    elseif (chirality1.eq.0 .and. chirality2.eq.0) then
       call QuarkAntiquarkToColourFlowU1Vector(wfq1,wfq2,wfg,coupl)
       return
    endif

    if (chirality1.eq.0) then
       q1=wfq1(1); q2=wfq1(2); q3=wfq1(3); q4=wfq1(4)
    elseif (chirality1.eq.1) then
       q1=wfq1(1); q2=wfq1(2); q3=(0d0,0d0); q4=(0d0,0d0)
    else
       q1=(0d0,0d0); q2=(0d0,0d0); q3=wfq1(1); q4=wfq1(2)
    endif
    if (chirality2.eq.0) then
       a1=wfq2(1); a2=wfq2(2); a3=wfq2(3); a4=wfq2(4)
    elseif (chirality2.eq.1) then
       a1=wfq2(1); a2=wfq2(2); a3=(0d0,0d0); a4=(0d0,0d0)
    else
       a1=(0d0,0d0); a2=(0d0,0d0); a3=wfq2(1); a4=wfq2(2)
    endif
    TMP1=q3*a1+q2*a4
    TMP2=q4*a2+q1*a3
    TMP3=q2*a3-q4*a1
    TMP4=q1*a4-q3*a2
    wfg(1)=( TMP1 + TMP2 )*TMP5
    wfg(2)=( TMP4 + TMP3 )*TMP5
    wfg(3)=(-TMP4 + TMP3 )*cImag*TMP5
    wfg(4)=(-TMP1 + TMP2 )*TMP5
  end subroutine QuarkAntiquarkToColourFlowU1Vector_weyl

  subroutine AntiquarkQuarkToGluon_weyl(wfq1,wfq2,wfg,chirality1,chirality2)
    implicit none
    integer,intent(in) :: chirality1,chirality2
    complex(kind=8),dimension(*) :: wfq1,wfq2
    complex(kind=8),dimension(4) :: wfg
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: q1,q2,q3,q4,a1,a2,a3,a4,TMP1,TMP2,TMP3,TMP4
    wfg(1:4)=(0d0,0d0)
    if (chirality1.eq.1 .and. chirality2.eq.-1) then
       wfg(1)=prefact*(wfq2(1)*wfq1(1)+wfq2(2)*wfq1(2))
       wfg(2)=-prefact*(wfq2(1)*wfq1(2)+wfq2(2)*wfq1(1))
       wfg(3)=cImag*prefact*(wfq2(1)*wfq1(2)-wfq2(2)*wfq1(1))
       wfg(4)=prefact*(-wfq2(1)*wfq1(1)+wfq2(2)*wfq1(2))
       return
    elseif (chirality1.eq.-1 .and. chirality2.eq.1) then
       wfg(1)=prefact*(wfq2(2)*wfq1(2)+wfq2(1)*wfq1(1))
       wfg(2)=prefact*(wfq2(1)*wfq1(2)+wfq2(2)*wfq1(1))
       wfg(3)=cImag*prefact*(-wfq2(1)*wfq1(2)+wfq2(2)*wfq1(1))
       wfg(4)=prefact*(-wfq2(2)*wfq1(2)+wfq2(1)*wfq1(1))
       return
    elseif (chirality1.ne.0 .and. chirality2.ne.0) then
       return
    elseif (chirality1.eq.0 .and. chirality2.eq.0) then
       call AntiquarkQuarkToGluon(wfq1,wfq2,wfg)
       return
    endif

    if (chirality2.eq.0) then
       q1=wfq2(1); q2=wfq2(2); q3=wfq2(3); q4=wfq2(4)
    elseif (chirality2.eq.1) then
       q1=wfq2(1); q2=wfq2(2); q3=(0d0,0d0); q4=(0d0,0d0)
    else
       q1=(0d0,0d0); q2=(0d0,0d0); q3=wfq2(1); q4=wfq2(2)
    endif
    if (chirality1.eq.0) then
       a1=wfq1(1); a2=wfq1(2); a3=wfq1(3); a4=wfq1(4)
    elseif (chirality1.eq.1) then
       a1=wfq1(1); a2=wfq1(2); a3=(0d0,0d0); a4=(0d0,0d0)
    else
       a1=(0d0,0d0); a2=(0d0,0d0); a3=wfq1(1); a4=wfq1(2)
    endif
    TMP1=q3*a1+q2*a4
    TMP2=q4*a2+q1*a3
    TMP3=q2*a3-q4*a1
    TMP4=q1*a4-q3*a2
    wfg(1)=( TMP1 + TMP2 )*prefact
    wfg(2)=( TMP4 + TMP3 )*prefact
    wfg(3)=(-TMP4 + TMP3 )*cImag*prefact
    wfg(4)=(-TMP1 + TMP2 )*prefact
  end subroutine AntiquarkQuarkToGluon_weyl

  subroutine LeptonAntileptonToVector_weyl(wfq1,wfq2,wfg,coupl,chirality1,chirality2)
    implicit none
    integer,intent(in) :: chirality1,chirality2
    complex(kind=8),dimension(*) :: wfq1,wfq2
    complex(kind=8),dimension(4) :: wfg
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: l1,l2,l3,l4,a1,a2,a3,a4,TMP5
    wfg(1:4)=(0d0,0d0)
    if (chirality1.eq.-1 .and. chirality2.eq.1) then
       TMP5=prefact*coupl(1)
       wfg(1)=TMP5*(wfq1(1)*wfq2(1)+wfq1(2)*wfq2(2))
       wfg(2)=-TMP5*(wfq1(2)*wfq2(1)+wfq1(1)*wfq2(2))
       wfg(3)=cImag*TMP5*(-wfq1(2)*wfq2(1)+wfq1(1)*wfq2(2))
       wfg(4)=TMP5*(-wfq1(1)*wfq2(1)+wfq1(2)*wfq2(2))
       return
    elseif (chirality1.eq.1 .and. chirality2.eq.-1) then
       TMP5=prefact*coupl(2)
       wfg(1)=TMP5*(wfq1(1)*wfq2(1)+wfq1(2)*wfq2(2))
       wfg(2)=TMP5*(wfq1(1)*wfq2(2)+wfq1(2)*wfq2(1))
       wfg(3)=cImag*TMP5*(-wfq1(1)*wfq2(2)+wfq1(2)*wfq2(1))
       wfg(4)=TMP5*(wfq1(1)*wfq2(1)-wfq1(2)*wfq2(2))
       return
    elseif (chirality1.ne.0 .and. chirality2.ne.0) then
       return
    elseif (chirality1.eq.0 .and. chirality2.eq.0) then
       call LeptonAntileptonToVector(wfq1,wfq2,wfg,coupl)
       return
    endif

    if (chirality1.eq.0) then
       l1=wfq1(1); l2=wfq1(2); l3=wfq1(3); l4=wfq1(4)
    elseif (chirality1.eq.1) then
       l1=wfq1(1); l2=wfq1(2); l3=(0d0,0d0); l4=(0d0,0d0)
    else
       l1=(0d0,0d0); l2=(0d0,0d0); l3=wfq1(1); l4=wfq1(2)
    endif
    if (chirality2.eq.0) then
       a1=wfq2(1); a2=wfq2(2); a3=wfq2(3); a4=wfq2(4)
    elseif (chirality2.eq.1) then
       a1=wfq2(1); a2=wfq2(2); a3=(0d0,0d0); a4=(0d0,0d0)
    else
       a1=(0d0,0d0); a2=(0d0,0d0); a3=wfq2(1); a4=wfq2(2)
    endif
    wfg(1)=prefact*(coupl(1)*(l3*a1+l4*a2)+coupl(2)*(l1*a3+l2*a4))
    wfg(2)=prefact*(coupl(1)*(-l4*a1-l3*a2)+coupl(2)*(l1*a4+l2*a3))
    wfg(3)=cImag*prefact*(coupl(1)*(-l4*a1+l3*a2)+coupl(2)*(-l1*a4+l2*a3))
    wfg(4)=prefact*(coupl(1)*(-l3*a1+l4*a2)+coupl(2)*(l1*a3-l2*a4))
  end subroutine LeptonAntileptonToVector_weyl

  subroutine AntileptonLeptonToVector_weyl(wfq1,wfq2,wfg,coupl,chirality1,chirality2)
    implicit none
    integer,intent(in) :: chirality1,chirality2
    complex(kind=8),dimension(*) :: wfq1,wfq2
    complex(kind=8),dimension(4) :: wfg
    real(kind=8),dimension(2) :: coupl
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: l1,l2,l3,l4,a1,a2,a3,a4,TMP5
    wfg(1:4)=(0d0,0d0)
    if (chirality1.eq.1 .and. chirality2.eq.-1) then
       TMP5=prefact*coupl(1)
       wfg(1)=TMP5*(wfq2(1)*wfq1(1)+wfq2(2)*wfq1(2))
       wfg(2)=-TMP5*(wfq2(2)*wfq1(1)+wfq2(1)*wfq1(2))
       wfg(3)=cImag*TMP5*(-wfq2(2)*wfq1(1)+wfq2(1)*wfq1(2))
       wfg(4)=TMP5*(-wfq2(1)*wfq1(1)+wfq2(2)*wfq1(2))
       return
    elseif (chirality1.eq.-1 .and. chirality2.eq.1) then
       TMP5=prefact*coupl(2)
       wfg(1)=TMP5*(wfq2(1)*wfq1(1)+wfq2(2)*wfq1(2))
       wfg(2)=TMP5*(wfq2(1)*wfq1(2)+wfq2(2)*wfq1(1))
       wfg(3)=cImag*TMP5*(-wfq2(1)*wfq1(2)+wfq2(2)*wfq1(1))
       wfg(4)=TMP5*(wfq2(1)*wfq1(1)-wfq2(2)*wfq1(2))
       return
    elseif (chirality1.ne.0 .and. chirality2.ne.0) then
       return
    elseif (chirality1.eq.0 .and. chirality2.eq.0) then
       call AntileptonLeptonToVector(wfq1,wfq2,wfg,coupl)
       return
    endif

    if (chirality2.eq.0) then
       l1=wfq2(1); l2=wfq2(2); l3=wfq2(3); l4=wfq2(4)
    elseif (chirality2.eq.1) then
       l1=wfq2(1); l2=wfq2(2); l3=(0d0,0d0); l4=(0d0,0d0)
    else
       l1=(0d0,0d0); l2=(0d0,0d0); l3=wfq2(1); l4=wfq2(2)
    endif
    if (chirality1.eq.0) then
       a1=wfq1(1); a2=wfq1(2); a3=wfq1(3); a4=wfq1(4)
    elseif (chirality1.eq.1) then
       a1=wfq1(1); a2=wfq1(2); a3=(0d0,0d0); a4=(0d0,0d0)
    else
       a1=(0d0,0d0); a2=(0d0,0d0); a3=wfq1(1); a4=wfq1(2)
    endif
    wfg(1)=prefact*(coupl(1)*(l3*a1+l4*a2)+coupl(2)*(l1*a3+l2*a4))
    wfg(2)=prefact*(coupl(1)*(-l4*a1-l3*a2)+coupl(2)*(l1*a4+l2*a3))
    wfg(3)=cImag*prefact*(coupl(1)*(-l4*a1+l3*a2)+coupl(2)*(-l1*a4+l2*a3))
    wfg(4)=prefact*(coupl(1)*(-l3*a1+l4*a2)+coupl(2)*(l1*a3-l2*a4))
  end subroutine AntileptonLeptonToVector_weyl
  
  subroutine MasslessVectorPropagator(wfg,p)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfg
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    propagator=-cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
    wfg(1:4)=wfg(1:4)*propagator
  end subroutine MasslessVectorPropagator

  subroutine MassiveVectorPropagator(wfg,p,vm,vw)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfg
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8) :: vm,vw
    complex(kind=8) :: TMP
    propagator=-cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2-vm**2+cImag*vm*vw)
    TMP=(p(0)*wfg(1)-p(1)*wfg(2)-p(2)*wfg(3)-p(3)*wfg(4))/vm**2
    wfg(1)=(wfg(1)-p(0)*TMP)*propagator
    wfg(2)=(wfg(2)-p(1)*TMP)*propagator
    wfg(3)=(wfg(3)-p(2)*TMP)*propagator
    wfg(4)=(wfg(4)-p(3)*TMP)*propagator
  end subroutine MassiveVectorPropagator

  subroutine FermionPropagator(wfq,p,fm,fw)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfq
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: prefact
    complex(kind=8),dimension(1:4) :: tmp_p,tmp_val
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8) :: fm,fw
    prefact=cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2-fm**2+cImag*fm*fw)
    tmp_val(1:4)=wfq(1:4)
    tmp_p(1)=(p(0)+p(3))
    tmp_p(2)=(p(0)-p(3))
    tmp_p(3)=(p(1)+cImag*p(2))
    tmp_p(4)=(p(1)-cImag*p(2))
    wfq(1)=(tmp_p(1)*tmp_val(3)+tmp_p(3)*tmp_val(4)+fm*tmp_val(1))*prefact
    wfq(2)=(tmp_p(2)*tmp_val(4)+tmp_p(4)*tmp_val(3)+fm*tmp_val(2))*prefact
    wfq(3)=(tmp_p(2)*tmp_val(1)-tmp_p(3)*tmp_val(2)+fm*tmp_val(3))*prefact
    wfq(4)=(tmp_p(1)*tmp_val(2)-tmp_p(4)*tmp_val(1)+fm*tmp_val(4))*prefact
  end subroutine FermionPropagator

  subroutine FermionPropagator_weyl(wfq,p,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: prefact,tmp1,tmp2,tmp3,tmp4,val1,val2
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    prefact=cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
    tmp1=(p(0)+p(3))
    tmp2=(p(0)-p(3))
    tmp3=(p(1)+cImag*p(2))
    tmp4=(p(1)-cImag*p(2))
    val1=wfq(1)
    val2=wfq(2)
    if (chirality.eq.1) then
       wfq(1)=(tmp1*val1+tmp3*val2)*prefact
       wfq(2)=(tmp2*val2+tmp4*val1)*prefact
    elseif (chirality.eq.-1) then
       wfq(1)=(tmp2*val1-tmp3*val2)*prefact
       wfq(2)=(tmp1*val2-tmp4*val1)*prefact
    else
       write (*,*) 'FermionPropagator_weyl called with zero chirality'
       stop 1
    endif
  end subroutine FermionPropagator_weyl

  subroutine AntifermionPropagator(wfq,p,fm,fw)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfq
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: prefact
    complex(kind=8),dimension(1:4) :: tmp_p,tmp_val
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8) :: fm,fw
    prefact=cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2-fm**2+cImag*fm*fw)
    tmp_val(1:4)=wfq(1:4)
    tmp_p(1)=-(p(0)+p(3))
    tmp_p(2)=-(p(0)-p(3))
    tmp_p(3)=-(p(1)+cImag*p(2))
    tmp_p(4)=-(p(1)-cImag*p(2))
    wfq(1)=(tmp_p(2)*tmp_val(3)-tmp_p(4)*tmp_val(4)+fm*tmp_val(1))*prefact
    wfq(2)=(tmp_p(1)*tmp_val(4)-tmp_p(3)*tmp_val(3)+fm*tmp_val(2))*prefact
    wfq(3)=(tmp_p(1)*tmp_val(1)+tmp_p(4)*tmp_val(2)+fm*tmp_val(3))*prefact
    wfq(4)=(tmp_p(2)*tmp_val(2)+tmp_p(3)*tmp_val(1)+fm*tmp_val(4))*prefact
  end subroutine AntifermionPropagator

  subroutine AntifermionPropagator_weyl(wfq,p,chirality)
    implicit none
    integer,intent(in) :: chirality
    complex(kind=8),dimension(*) :: wfq
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: prefact,tmp1,tmp2,tmp3,tmp4,val1,val2
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    prefact=cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
    tmp1=-(p(0)+p(3))
    tmp2=-(p(0)-p(3))
    tmp3=-(p(1)+cImag*p(2))
    tmp4=-(p(1)-cImag*p(2))
    val1=wfq(1)
    val2=wfq(2)
    if (chirality.eq.1) then
       wfq(1)=(tmp2*val1-tmp4*val2)*prefact
       wfq(2)=(tmp1*val2-tmp3*val1)*prefact
    elseif (chirality.eq.-1) then
       wfq(1)=(tmp1*val1+tmp4*val2)*prefact
       wfq(2)=(tmp2*val2+tmp3*val1)*prefact
    else
       write (*,*) 'AntifermionPropagator_weyl called with zero chirality'
       stop 1
    endif
  end subroutine AntifermionPropagator_weyl

  subroutine ScalarPropagator(wfs,p,sm,sw)
    implicit none
    complex(kind=8),dimension(1),intent(inout) :: wfs
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8) :: sm,sw
    propagator=cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2-sm**2+cImag*sm*sw)
    wfs(1)=wfs(1)*propagator
  end subroutine ScalarPropagator
end module FeynmanRules
