module FeynmanRules
contains
  subroutine ext_gluon_real(p,ihel,ifinal,wf)
    implicit none
    integer ihel,ifinal
    real(kind=8), dimension(0:3) :: p
    real(kind=8), dimension(4) :: wf
    complex(kind=8),dimension(4) :: wf0,wf1
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8),parameter :: sqh=sqrt(0.5d0)
    call ext_gluon_cmplx(p, 1,ifinal,wf1)
    call ext_gluon_cmplx(p,-1,ifinal,wf0)
    if (ihel.eq.1) then
       wf(1:4)=dble(cImag*(wf1(1:4)+wf0(1:4)))*sqh
    elseif (ihel.eq.-1) then
       wf(1:4)=-dble(wf1(1:4)-wf0(1:4))*sqh
    endif
  end subroutine ext_gluon_real
  subroutine ext_gluon_cmplx(p,ihel,idum,wf)
    ! External gluon wavefunction. From HELAS.
    implicit none
    integer :: ihel,idum
    real(kind=8), dimension(0:3) :: p
    complex(kind=8), dimension(4) :: wf
    real(kind=8),parameter :: rZero=0d0,sqh=sqrt(0.5d0)
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8) :: hel,pp,pt,pzpt
    if (p(0).eq.0d0) then
       write (*,*) 'Cannot generate external gluon with zero energy'
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
       hel = dble(ihel)
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
  
  end subroutine ext_gluon_cmplx


  subroutine ext_gluon_mass(p,nhel,nsv,wf,vmass)
      implicit none
      double complex wf(4)
      double precision p(0:3),vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
      integer nhel,nsv,nsvahl

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )

      sqh = dsqrt(rHalf)
      hel = dble(nhel)
      nsvahl = nsv*dabs(hel)
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


  end subroutine ext_gluon_mass

  subroutine ext_quark(p,nhel,idum,wf,fmass)
  ! flowing-out fermion number, i.e., final state quark (p(0)>0) or initial
  ! state anti-quark (p(0)<0)
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
             ! oxxx, nsf=+1, nhel=+1
             wf(1) = chi(1)
             wf(2) = chi(2)
             wf(3) = cZero
             wf(4) = cZero
          else
             ! oxxx, nsf=+1, nhel=-1
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(2)
             wf(4) = chi(1)
          endif
       else
          write (*,*) 'CHECK THIS'
          stop 1
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
             ! oxxx, nsf=-1, nhel=+1
             wf(1) = chi(1)
             wf(2) = chi(2)
             wf(3) = cZero
             wf(4) = cZero
          else
             ! oxxx, nsf=-1, nhel=-1
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(2)
             wf(4) = chi(1)
          endif
       else
          write (*,*) 'CHECK THIS'
          stop 1
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
  end subroutine ext_quark


  subroutine ext_antiquark(p,nhel,idum,wf,fmass)
    ! flowing-in fermion number, i.e., final state anti-quark (p(0)>0), or
    ! initial state quark (p(0)<0)
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
             ! ixxx, nsf=-1, nhel=-1
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(1)
             wf(4) = chi(2)
          else
             ! ixxx, nsf=-1, nhel=+1
             wf(1) = chi(2)
             wf(2) = chi(1)
             wf(3) = cZero
             wf(4) = cZero
          endif
       else
          write (*,*) 'CHECK THIS'
          stop 1
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
             ! ixxx, nsf=+1, nhel=+1
             wf(1) = cZero
             wf(2) = cZero
             wf(3) = chi(1)
             wf(4) = chi(2)
          else
             ! ixxx, nsf=+1, nhel=-1
             wf(1) = chi(2)
             wf(2) = chi(1)
             wf(3) = cZero
             wf(4) = cZero
          endif
       else
          write (*,*) 'CHECK THIS'
          stop 1
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
  end subroutine ext_antiquark

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
  subroutine ThreeGluon_real(wf1,pwf1,wf2,pwf2,wf)
    ! Colour-ordered three-gluon interaction
    implicit none
    real(kind=8),dimension(4) :: wf1,wf2,wf
    real(kind=8),dimension(0:3) :: pwf1,pwf2
    real(kind=8),parameter :: prefact=1d0/sqrt(2d0)
    real(kind=8) :: TMP1,TMP2,TMP3
    TMP1 = (wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4))
    TMP2 = (wf1(1)*pwf2(0)-wf1(2)*pwf2(1)-wf1(3)*pwf2(2)-wf1(4)*pwf2(3))
    TMP3 = (wf2(1)*pwf1(0)-wf2(2)*pwf1(1)-wf2(3)*pwf1(2)-wf2(4)*pwf1(3))
    wf(1:4) = prefact*(TMP1*(pwf1(0:3)-pwf2(0:3))+2d0*(TMP2*wf2(1:4)-TMP3*wf1(1:4)))
  end subroutine ThreeGluon_Real
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
  subroutine TwoGluontoTensor(wfg1,wfg2,wfT)
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
!!$    do i=1,4
!!$       wfT(1:4,i)=(wfg1(1:4)*wfg2(i)-wfg2(1:4)*wfg1(i))
!!$    enddo
  end subroutine TwoGluontoTensor
  subroutine TwoGluontoTensor_real(wfg1,wfg2,wfT)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
    implicit none
    real(kind=8),dimension(4) :: wfg1,wfg2
    real(kind=8),dimension(6) :: wfT
    wfT(1)=(wfg1(1)*wfg2(2)-wfg1(2)*wfg2(1))
    wfT(2)=(wfg1(1)*wfg2(3)-wfg1(3)*wfg2(1))
    wfT(3)=(wfg1(1)*wfg2(4)-wfg1(4)*wfg2(1))
    wfT(4)=(wfg1(2)*wfg2(3)-wfg1(3)*wfg2(2))
    wfT(5)=(wfg1(2)*wfg2(4)-wfg1(4)*wfg2(2))
    wfT(6)=(wfg1(3)*wfg2(4)-wfg1(4)*wfg2(3))
  end subroutine TwoGluontoTensor_real
  subroutine TensorGluontoGluon(wfT1,wfg2,wfg)
    implicit none
    complex(kind=8),dimension(4) :: wfg2,wfg
    complex(kind=8),dimension(6) :: wfT1
!!$    complex(kind=8),parameter :: prefact=(-0.5d0,0d0)
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact
!!$    do i=1,4
!!$       wfg(i)=((wfT1(1,i)*wfg2(1)-wfT1(2,i)*wfg2(2)-wfT1(3,i)*wfg2(3)-wfT1(4,i)*wfg2(4))- &
!!$               (wfT1(i,1)*wfg2(1)-wfT1(i,2)*wfg2(2)-wfT1(i,3)*wfg2(3)-wfT1(i,4)*wfg2(4)))*0.25d0
!!$    enddo
  end subroutine TensorGluontoGluon
  subroutine TensorGluontoGluon_real(wfT1,wfg2,wfg)
    implicit none
    real(kind=8),dimension(4) :: wfg2,wfg
    real(kind=8),dimension(6) :: wfT1
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact
  end subroutine TensorGluontoGluon_Real
  subroutine GluonTensortoGluon(wfg1,wfT2,wfg)
    implicit none 
    complex(kind=8),dimension(4) :: wfg1,wfg
    complex(kind=8),dimension(6) :: wfT2
!!$    complex(kind=8),parameter :: prefact=(-0.5d0,0d0)
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact
!!$    do i=1,4
!!$       wfg(i)=-((wfg1(1)*wfT2(1,i)-wfg1(2)*wfT2(2,i)-wfg1(3)*wfT2(3,i)-wfg1(4)*wfT2(4,i))- &
!!$               (wfg1(1)*wfT2(i,1)-wfg1(2)*wfT2(i,2)-wfg1(3)*wfT2(i,3)-wfg1(4)*wfT2(i,4)))*0.25d0
!!$    enddo
  end subroutine GluonTensortoGluon
  subroutine GluonTensortoGluon_real(wfg1,wfT2,wfg)
    implicit none 
    real(kind=8),dimension(4) :: wfg1,wfg
    real(kind=8),dimension(6) :: wfT2
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact
  end subroutine GluonTensortoGluon_Real







  subroutine QuarkGluontoQuark(wfq1,wfg2,wfq) ! from fvoxxx.f
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
  end subroutine QuarkGluontoQuark

  subroutine QuarkGluontoQuark_real(wfq1,wfg2,wfq)
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq
    real(kind=8),dimension(4) :: wfg2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8) :: TMP1,TMP2
    complex(kind=8) :: TMP3,TMP4
    TMP1=wfg2(1)+wfg2(4)
    TMP2=wfg2(1)-wfg2(4)
    TMP3=dcmplx(wfg2(2),wfg2(3))
    TMP4=dcmplx(wfg2(2),-wfg2(3))
    wfq(1)=prefact*(TMP1*wfq1(3)+TMP3*wfq1(4)) !sl1
    wfq(2)=prefact*(TMP2*wfq1(4)+TMP4*wfq1(3)) !sl2
    wfq(3)=prefact*(TMP2*wfq1(1)-TMP3*wfq1(2)) !sr1
    wfq(4)=prefact*(TMP1*wfq1(2)-TMP4*wfq1(1)) !sr2
  end subroutine QuarkGluontoQuark_real

  subroutine GluonQuarktoQuark(wfg1,wfq2,wfq) ! from fvoxxx.f
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
  end subroutine GluonQuarktoQuark

  subroutine GluonQuarktoQuark_real(wfg1,wfq2,wfq) ! from fvoxxx.f
    implicit none
    complex(kind=8),dimension(4) :: wfq2,wfq
    real(kind=8),dimension(4) :: wfg1
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8) :: TMP1,TMP2
    complex(kind=8) :: TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=dcmplx(wfg1(2),wfg1(3))
    TMP4=dcmplx(wfg1(2),-wfg1(3))
    wfq(1)=prefact*(TMP1*wfq2(3)+TMP3*wfq2(4)) ! sl1 ! minus sign
    wfq(2)=prefact*(TMP2*wfq2(4)+TMP4*wfq2(3)) ! sl2
    wfq(3)=prefact*(TMP2*wfq2(1)-TMP3*wfq2(2)) ! sr1
    wfq(4)=prefact*(TMP1*wfq2(2)-TMP4*wfq2(1)) ! sr2
  end subroutine GluonQuarktoQuark_real

  subroutine AquarkGluontoAquark(wfq1,wfg2,wfq) ! TV from fvixxx.f
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
  end subroutine AquarkGluontoAquark

  subroutine GluonAquarktoAquark(wfg1,wfq2,wfq) ! TV from fvixxx.f
    implicit none
    complex(kind=8),dimension(4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfg1(1)+wfg1(4)
    TMP2=wfg1(1)-wfg1(4)
    TMP3=wfg1(2)+cImag*wfg1(3)
    TMP4=wfg1(2)-cImag*wfg1(3)
    wfq(1)=prefact*(TMP2*wfq2(3)-TMP4*wfq2(4)) !sr1 ! minus sign!
    wfq(2)=prefact*(TMP1*wfq2(4)-TMP3*wfq2(3)) !sr2
    wfq(3)=prefact*(TMP1*wfq2(1)+TMP4*wfq2(2)) !sl1
    wfq(4)=prefact*(TMP2*wfq2(2)+TMP3*wfq2(1)) !sl2
  end subroutine GluonAquarktoAquark

  subroutine QuarkAquarktoGluon(wfq1,wfq2,wfg) ! TV from jioxxx.f
    implicit none
    complex(kind=8),dimension(4) :: wfq1,wfq2,wfg
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8) :: TMP1,TMP2,TMP3,TMP4
    TMP1=wfq1(3)*wfq2(1)+wfq1(2)*wfq2(4)
    TMP2=wfq1(4)*wfq2(2)+wfq1(1)*wfq2(3)
    TMP3=wfq1(2)*wfq2(3)-wfq1(4)*wfq2(1)
    TMP4=wfq1(1)*wfq2(4)-wfq1(3)*wfq2(2)
    wfg(1)=( TMP1 + TMP2 )*prefact
    wfg(2)=( TMP4 + TMP3 )*prefact
    wfg(3)=(-TMP4 + TMP3 )*cImag*prefact
    wfg(4)=(-TMP1 + TMP2 )*prefact
  end subroutine QuarkAquarktoGluon
  subroutine AquarkQuarktoGluon(wfq1,wfq2,wfg) ! TV from jioxxx.f
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
  end subroutine AquarkQuarktoGluon


  subroutine ThreeGluon_coupl(wf1,pwf1,wf2,pwf2,wf,coupl)
    ! Colour-ordered three-gluon interaction
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
  end subroutine ThreeGluon_Coupl


  subroutine QuarkGluontoQuark_coupl(wfq1,wfg2,wfqcoupl) ! from fvoxxx.f
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
  end subroutine QuarkGluontoQuark_coupl

  subroutine AquarkGluontoAquark_coupl(wfq1,wfg2,wfq,coupl) ! TV from fvixxx.f
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
  end subroutine AquarkGluontoAquark_coupl

  subroutine TwoGluontoTensor_coupl(wfg1,wfg2,wfT,coupl)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
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
  end subroutine TwoGluontoTensor_coupl

  subroutine TensorGluontoGluon_coupl(wfT1,wfg2,wfg,coupl)
    implicit none
    complex(kind=8),dimension(4) :: wfg2,wfg
    complex(kind=8),dimension(6) :: wfT1
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact*coupl(1)
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact*coupl(1)
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact*coupl(1)
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact*coupl(1)
  end subroutine TensorGluontoGluon_coupl

  subroutine GluonTensortoGluon_coupl(wfg1,wfT2,wfg,coupl)
    implicit none 
    complex(kind=8),dimension(4) :: wfg1,wfg
    complex(kind=8),dimension(6) :: wfT2
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    real(kind=8),dimension(2) :: coupl
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact*coupl(1)
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact*coupl(1)
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact*coupl(1)
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact*coupl(1)
  end subroutine GluonTensortoGluon_coupl









  subroutine GluonPropagator(wfg,p)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfg
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    propagator=-cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
    wfg(1:4)=wfg(1:4)*propagator
  end subroutine GluonPropagator

  subroutine GluonPropagator_real(wfg,p)
    implicit none
    real(kind=8),dimension(1:4),intent(inout) :: wfg
    real(kind=8),dimension(0:3),intent(in) :: p
    real(kind=8) :: propagator
    propagator=1d0/(p(0)**2-p(1)**2-p(2)**2-p(3)**2)
    wfg(1:4)=wfg(1:4)*propagator
  end subroutine GluonPropagator_Real

  subroutine GluonPropagator_mass(wfg,p,vm,vw)
    implicit none
    complex(kind=8),dimension(1:4),intent(inout) :: wfg
    real(kind=8),dimension(0:3),intent(in) :: p
    complex(kind=8) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8) :: vm,vw
    complex(kind=8) :: TMP
    propagator=-cImag/(p(0)**2-p(1)**2-p(2)**2-p(3)**2-vm**2+cImag*vm*vw)
!!$    TMP=(p(0)*wfg(1)-p(1)*wfg(2)-p(2)*wfg(3)-p(3)*wfg(4))/(vm**2-cImag*vm*vw)
    TMP=(p(0)*wfg(1)-p(1)*wfg(2)-p(2)*wfg(3)-p(3)*wfg(4))/vm**2
    wfg(1)=(wfg(1)-p(0)*TMP)*propagator
    wfg(2)=(wfg(2)-p(1)*TMP)*propagator
    wfg(3)=(wfg(3)-p(2)*TMP)*propagator
    wfg(4)=(wfg(4)-p(3)*TMP)*propagator
  end subroutine GluonPropagator_mass

  subroutine QuarkPropagator(wfq,p,fm,fw)
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
  end subroutine QuarkPropagator

  subroutine AquarkPropagator(wfq,p,fm,fw)
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
  end subroutine AquarkPropagator

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
