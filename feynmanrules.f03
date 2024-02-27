module FeynmanRules
contains
  subroutine ext_gluon_real(np,p,ihel,ifinal,wf)
    implicit none
    integer,intent(in) :: np
    integer ihel,ifinal
    real(kind=8), dimension(np,0:3) :: p
    real(kind=8), dimension(np,4) :: wf
    complex(kind=8),dimension(np,4) :: wf0,wf1
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    real(kind=8),parameter :: sqh=sqrt(0.5d0)
    call ext_gluon_cmplx(np,p,1,ifinal,wf1)
    call ext_gluon_cmplx(np,p,0,ifinal,wf0)
    if (ihel.eq.1) then
       wf(1:np,1:4)=dble(cImag*(wf1(1:np,1:4)+wf0(1:np,1:4)))*sqh
    elseif (ihel.eq.0) then
       wf(1:np,1:4)=-dble(wf1(1:np,1:4)-wf0(1:np,1:4))*sqh
    endif
  end subroutine ext_gluon_real
  subroutine ext_gluon_cmplx(np,p,ihel,ifinal,wf)
    ! External gluon wavefunction. From HELAS.
    implicit none
    integer,intent(in) :: np
    integer :: ihel,ifinal
    real(kind=8), dimension(np,0:3) :: p
    complex(kind=8), dimension(np,4) :: wf
    real(kind=8),parameter :: rZero=0d0,sqh=sqrt(0.5d0)
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8) :: hel
    real(kind=8),dimension(np) :: pp,pt,pzpt
!!$    hel = dble(2*ihel-1)
!!$    pp = p(0)
!!$    pt = sqrt(p(1)**2+p(2)**2)
!!$    wf(1) = cZero
!!$    wf(4) = dcmplx( hel*pt/pp*sqh )
!!$    if ( pt.ne.rZero ) then
!!$       pzpt = p(3)/(pp*pt)*sqh*hel
!!$       wf(2) = dcmplx( -p(1)*pzpt , -ifinal*p(2)/pt*sqh )
!!$       wf(3) = dcmplx( -p(2)*pzpt ,  ifinal*p(1)/pt*sqh )
!!$    else
!!$       wf(2) = dcmplx( -hel*sqh )
!!$       wf(3) = dcmplx( rZero , ifinal*sign(sqh,p(3)) )
!!$    endif

    if (p(1,0).gt.0d0) then
       hel = dble(2*ihel-1)
       pp(1:np) = p(1:np,0)
       pt(1:np) = sqrt(p(1:np,1)**2+p(1:np,2)**2)
       wf(1:np,1) = dcmplx( rZero )
       wf(1:np,4) = dcmplx( hel*pt(1:np)/pp(1:np)*sqh )
       if ( pt(1).ne.rZero ) then
          pzpt(1:np) = p(1:np,3)/(pp(1:np)*pt(1:np))*sqh*hel
          wf(1:np,2) = dcmplx( -p(1:np,1)*pzpt(1:np) , -ifinal*p(1:np,2)/pt(1:np)*sqh )
          wf(1:np,3) = dcmplx( -p(1:np,2)*pzpt(1:np) ,  ifinal*p(1:np,1)/pt(1:np)*sqh )
       else
          wf(1:np,2) = dcmplx( -hel*sqh )
          wf(1:np,3) = dcmplx( rZero , ifinal*sign(sqh,p(1:np,3)) )
       endif
    else
!!$       hel = -dble(2*ihel-1)
       hel = dble(2*ihel-1)
       pp(1:np) = -p(1:np,0)
       pt(1:np) = sqrt(p(1:np,1)**2+p(1:np,2)**2)
       wf(1:np,1) = dcmplx( rZero )
       wf(1:np,4) = dcmplx( hel*pt(1:np)/pp(1:np)*sqh )
       if ( pt(1).ne.rZero ) then
          pzpt(1:np) = -p(1:np,3)/(pp(1:np)*pt(1:np))*sqh*hel
          wf(1:np,2) = dcmplx( p(1:np,1)*pzpt(1:np) , -ifinal*p(1:np,2)/pt(1:np)*sqh )
          wf(1:np,3) = dcmplx( p(1:np,2)*pzpt(1:np) ,  ifinal*p(1:np,1)/pt(1:np)*sqh )
       else
          wf(1:np,2) = dcmplx( -hel*sqh )
          wf(1:np,3) = dcmplx( rZero , ifinal*sign(sqh,p(1:np,3)) )
       endif
    endif
  
  end subroutine ext_gluon_cmplx
  subroutine ext_quark(np,p,ihel,ifinal,wf)
  ! flowing-out fermion number, i.e., final state quark (p(0)>0) or initial
  ! state anti-quark (p(0)<0)
    implicit none
    integer,intent(in) :: np
    integer :: ihel,ifinal
    real(kind=8), dimension(np,0:3) :: p
    complex(kind=8), dimension(np,4) :: wf
    complex(kind=8), dimension(np,2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8),dimension(np) :: sqp0p3
    integer :: nhel
    if (p(1,0).gt.0d0) then
       ! outgoing final state momenta
       nhel = 2*ihel-1
       if(p(1,1).eq.0d0.and.p(1,2).eq.0d0.and.p(1,3).lt.0d0) then
          sqp0p3(1:np) = 0d0
       else
          sqp0p3(1:np) = dsqrt(max(p(1:np,0)+p(1:np,3),rZero))
       end if
       chi(1:np,1) = dcmplx( sqp0p3(1:np) )
       if ( sqp0p3(1).eq.rZero ) then
          chi(1:np,2) = dcmplx(-nhel )*dsqrt(rTwo*p(1:np,0))
       else
          chi(1:np,2) = dcmplx( nhel*p(1:np,1), -p(1:np,2) )/sqp0p3(1:np)
       endif
       if ( nhel.eq.1 ) then
          wf(1:np,1) = chi(1:np,1)
          wf(1:np,2) = chi(1:np,2)
          wf(1:np,3) = cZero
          wf(1:np,4) = cZero
       else
          wf(1:np,1) = cZero
          wf(1:np,2) = cZero
          wf(1:np,3) = chi(1:np,2)
          wf(1:np,4) = chi(1:np,1)
       endif
    else
       ! "outgoing" initial state momenta
       nhel = (2*ihel-1)
       if(p(1,1).eq.0d0.and.p(1,2).eq.0d0.and.p(1,3).gt.0d0) then
          sqp0p3(1:np) = 0d0
       else
          sqp0p3(1:np) = -dsqrt(max(-(p(1:np,0)+p(1:np,3)),rZero))
       end if
       chi(1:np,1) = dcmplx( sqp0p3(1:np) )
       if ( sqp0p3(1).eq.rZero ) then
          chi(1:np,2) = dcmplx(-nhel )*dsqrt(rTwo*abs(p(1:np,0)))
       else
          chi(1:np,2) = dcmplx( -nhel*p(1:np,1), -p(1:np,2) )/sqp0p3(1:np)
       endif
       if ( -nhel.eq.1 ) then
          wf(1:np,1) = cZero
          wf(1:np,2) = cZero
          wf(1:np,3) = chi(1:np,2)
          wf(1:np,4) = chi(1:np,1)
       else
          wf(1:np,1) = chi(1:np,1)
          wf(1:np,2) = chi(1:np,2)
          wf(1:np,3) = cZero
          wf(1:np,4) = cZero
       endif
    endif
  end subroutine ext_quark
  subroutine ext_antiquark(np,p,ihel,ifinal,wf)
  ! flowing-in fermion number, i.e., final state anti-quark (p(0)>0), or
  ! initial state quark (p(0)<0)
    implicit none
    integer,intent(in) :: np
    integer :: ihel,ifinal
    real(kind=8), dimension(np,0:3) :: p
    complex(kind=8), dimension(np,4) :: wf
    complex(kind=8), dimension(np,2) :: chi
    real(kind=8),parameter :: rzero=0d0,rTwo=2d0
    complex(kind=8),parameter :: cZero=(0d0,0d0)
    real(kind=8),dimension(np) :: sqp0p3
    integer :: nhel
    if(p(1,0).gt.0d0) then
! outgoing final state momenta
       nhel = (2*ihel-1)
       if(p(1,1).eq.0d0.and.p(1,2).eq.0d0.and.p(1,3).lt.0d0) then
          sqp0p3(1:np) = 0d0
       else
          sqp0p3(1:np) = -dsqrt(max(p(1:np,0)+p(1:np,3),rZero))
       end if
       chi(1:np,1) = dcmplx( sqp0p3(1:np) )
       if ( sqp0p3(1).eq.rZero ) then
          chi(1:np,2) = dcmplx(-nhel )*dsqrt(rTwo*p(1:np,0))
       else
          chi(1:np,2) = dcmplx(-nhel*p(1:np,1), p(1:np,2) )/sqp0p3(1:np)
       endif
       if ( -nhel.eq.1 ) then
          wf(1:np,1) = cZero
          wf(1:np,2) = cZero
          wf(1:np,3) = chi(1:np,1)
          wf(1:np,4) = chi(1:np,2)
       else
          wf(1:np,1) = chi(1:np,2)
          wf(1:np,2) = chi(1:np,1)
          wf(1:np,3) = cZero
          wf(1:np,4) = cZero
       endif
    else
! "outgoing" initial state momenta
       nhel = 2*ihel-1
       if(p(1,1).eq.0d0.and.p(1,2).eq.0d0.and.p(1,3).gt.0d0) then
          sqp0p3(1:np) = 0d0
       else
          sqp0p3(1:np) = dsqrt(max(-(p(1:np,0)+p(1:np,3)),rZero))
       end if
       chi(1:np,1) = dcmplx( sqp0p3(1:np) )
       if ( sqp0p3(1).eq.rZero ) then
          chi(1:np,2) = dcmplx( -nhel )*dsqrt(rTwo*abs(p(1:np,0)))
       else
          chi(1:np,2) = dcmplx( nhel*p(1:np,1), p(1:np,2) )/sqp0p3(1:np)
       endif
       if ( nhel.eq.1 ) then
          wf(1:np,1) = chi(1:np,2)
          wf(1:np,2) = chi(1:np,1)
          wf(1:np,3) = cZero
          wf(1:np,4) = cZero
       else
          wf(1:np,1) = cZero
          wf(1:np,2) = cZero
          wf(1:np,3) = chi(1:np,1)
          wf(1:np,4) = chi(1:np,2)
       endif
    endif
  end subroutine ext_antiquark
  subroutine ThreeGluon(np,wf1,pwf1,wf2,pwf2,wf)
    ! Colour-ordered three-gluon interaction
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wf1,wf2,wf
    real(kind=8),dimension(np,0:3) :: pwf1,pwf2
    complex(kind=8),parameter :: prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8),dimension(np) :: TMP1,TMP2,TMP3
    integer :: i
    TMP1(1:np) = (wf1(1:np,1)*wf2(1:np,1)-wf1(1:np,2)*wf2(1:np,2)-wf1(1:np,3)*wf2(1:np,3)-wf1(1:np,4)*wf2(1:np,4))
    TMP2(1:np) = (wf1(1:np,1)*pwf2(1:np,0)-wf1(1:np,2)*pwf2(1:np,1)-wf1(1:np,3)*pwf2(1:np,2)-wf1(1:np,4)*pwf2(1:np,3))
    TMP3(1:np) = (wf2(1:np,1)*pwf1(1:np,0)-wf2(1:np,2)*pwf1(1:np,1)-wf2(1:np,3)*pwf1(1:np,2)-wf2(1:np,4)*pwf1(1:np,3))
    do i=1,4
       wf(1:np,i) = prefact*(TMP1(1:np)*(pwf1(1:np,i-1)-pwf2(1:np,i-1))+2d0*(TMP2(1:np)*wf2(1:np,i)-TMP3(1:np)*wf1(1:np,i)))
    enddo
  end subroutine ThreeGluon
  subroutine ThreeGluon_real(np,wf1,pwf1,wf2,pwf2,wf)
    ! Colour-ordered three-gluon interaction
    implicit none
    integer,intent(in) :: np
    real(kind=8),dimension(np,4) :: wf1,wf2,wf
    real(kind=8),dimension(np,0:3) :: pwf1,pwf2
    real(kind=8),parameter :: prefact=1d0/sqrt(2d0)
    real(kind=8),dimension(np) :: TMP1,TMP2,TMP3
    integer :: i
    TMP1(1:np) = (wf1(1:np,1)*wf2(1:np,1)-wf1(1:np,2)*wf2(1:np,2)-wf1(1:np,3)*wf2(1:np,3)-wf1(1:np,4)*wf2(1:np,4))
    TMP2(1:np) = (wf1(1:np,1)*pwf2(1:np,0)-wf1(1:np,2)*pwf2(1:np,1)-wf1(1:np,3)*pwf2(1:np,2)-wf1(1:np,4)*pwf2(1:np,3))
    TMP3(1:np) = (wf2(1:np,1)*pwf1(1:np,0)-wf2(1:np,2)*pwf1(1:np,1)-wf2(1:np,3)*pwf1(1:np,2)-wf2(1:np,4)*pwf1(1:np,3))
    do i=1,4
       wf(1:np,i) = prefact*(TMP1(1:np)*(pwf1(1:np,i-1)-pwf2(1:np,i-1))+2d0*(TMP2(1:np)*wf2(1:np,i)-TMP3(1:np)*wf1(1:np,i)))
    enddo
  end subroutine ThreeGluon_Real
  subroutine FourGluon(np,wf1,wf2,wf3,wf)
    ! Colour-ordered four-gluon interaction
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wf1,wf2,wf3,wf
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    complex(kind=8),dimension(np) :: TMP1,TMP2,TMP3
    integer :: i
    TMP1(1:np) = (wf1(1:np,1)*wf2(1:np,1)-wf1(1:np,2)*wf2(1:np,2)-wf1(1:np,3)*wf2(1:np,3)-wf1(1:np,4)*wf2(1:np,4))
    TMP2(1:np) = (wf1(1:np,1)*wf3(1:np,1)-wf1(1:np,2)*wf3(1:np,2)-wf1(1:np,3)*wf3(1:np,3)-wf1(1:np,4)*wf3(1:np,4))
    TMP3(1:np) = (wf2(1:np,1)*wf3(1:np,1)-wf2(1:np,2)*wf3(1:np,2)-wf2(1:np,3)*wf3(1:np,3)-wf2(1:np,4)*wf3(1:np,4))
    do i=1,4
       wf(1:np,i) = prefact*(2d0*wf2(1:np,i)*TMP2(1:np)-wf1(1:np,i)*TMP3(1:np)-wf3(1:np,i)*TMP1(1:np))
    enddo
  end subroutine FourGluon
  subroutine TwoGluontoTensor(np,wfg1,wfg2,wfT)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfg1,wfg2
    complex(kind=8),dimension(np,6) :: wfT
    complex(kind=8),parameter :: prefact=(0d0,1d0)
    ! Since it is an anti-symmetric 4x4 tensor, take only the upper-right triangle.
    wfT(1:np,1)=(wfg1(1:np,1)*wfg2(1:np,2)-wfg1(1:np,2)*wfg2(1:np,1))! * prefact
    wfT(1:np,2)=(wfg1(1:np,1)*wfg2(1:np,3)-wfg1(1:np,3)*wfg2(1:np,1))! * prefact
    wfT(1:np,3)=(wfg1(1:np,1)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,1))! * prefact
    wfT(1:np,4)=(wfg1(1:np,2)*wfg2(1:np,3)-wfg1(1:np,3)*wfg2(1:np,2))! * prefact
    wfT(1:np,5)=(wfg1(1:np,2)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,2))! * prefact
    wfT(1:np,6)=(wfg1(1:np,3)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,3))! * prefact
!!$    do i=1,4
!!$       wfT(1:4,i)=(wfg1(1:4)*wfg2(i)-wfg2(1:4)*wfg1(i))
!!$    enddo
  end subroutine TwoGluontoTensor
  subroutine TwoGluontoTensor_real(np,wfg1,wfg2,wfT)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
    implicit none
    integer,intent(in) :: np
    real(kind=8),dimension(np,4) :: wfg1,wfg2
    real(kind=8),dimension(np,6) :: wfT
    wfT(1:np,1)=(wfg1(1:np,1)*wfg2(1:np,2)-wfg1(1:np,2)*wfg2(1:np,1))
    wfT(1:np,2)=(wfg1(1:np,1)*wfg2(1:np,3)-wfg1(1:np,3)*wfg2(1:np,1))
    wfT(1:np,3)=(wfg1(1:np,1)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,1))
    wfT(1:np,4)=(wfg1(1:np,2)*wfg2(1:np,3)-wfg1(1:np,3)*wfg2(1:np,2))
    wfT(1:np,5)=(wfg1(1:np,2)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,2))
    wfT(1:np,6)=(wfg1(1:np,3)*wfg2(1:np,4)-wfg1(1:np,4)*wfg2(1:np,3))
  end subroutine TwoGluontoTensor_real
  subroutine TensorGluontoGluon(np,wfT1,wfg2,wfg)
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfg2,wfg
    complex(kind=8),dimension(np,6) :: wfT1
!!$    complex(kind=8),parameter :: prefact=(-0.5d0,0d0)
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1:np,1)=(wfT1(1:np,1)*wfg2(1:np,2)+wfT1(1:np,2)*wfg2(1:np,3)+wfT1(1:np,3)*wfg2(1:np,4))*prefact
    wfg(1:np,2)=(wfT1(1:np,1)*wfg2(1:np,1)+wfT1(1:np,4)*wfg2(1:np,3)+wfT1(1:np,5)*wfg2(1:np,4))*prefact
    wfg(1:np,3)=(wfT1(1:np,2)*wfg2(1:np,1)-wfT1(1:np,4)*wfg2(1:np,2)+wfT1(1:np,6)*wfg2(1:np,4))*prefact
    wfg(1:np,4)=(wfT1(1:np,3)*wfg2(1:np,1)-wfT1(1:np,5)*wfg2(1:np,2)-wfT1(1:np,6)*wfg2(1:np,3))*prefact
!!$    do i=1,4
!!$       wfg(i)=((wfT1(1,i)*wfg2(1)-wfT1(2,i)*wfg2(2)-wfT1(3,i)*wfg2(3)-wfT1(4,i)*wfg2(4))- &
!!$               (wfT1(i,1)*wfg2(1)-wfT1(i,2)*wfg2(2)-wfT1(i,3)*wfg2(3)-wfT1(i,4)*wfg2(4)))*0.25d0
!!$    enddo
  end subroutine TensorGluontoGluon
  subroutine TensorGluontoGluon_real(np,wfT1,wfg2,wfg)
    implicit none
    integer,intent(in) :: np
    real(kind=8),dimension(np,4) :: wfg2,wfg
    real(kind=8),dimension(np,6) :: wfT1
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1:np,1)=(wfT1(1:np,1)*wfg2(1:np,2)+wfT1(1:np,2)*wfg2(1:np,3)+wfT1(1:np,3)*wfg2(1:np,4))*prefact
    wfg(1:np,2)=(wfT1(1:np,1)*wfg2(1:np,1)+wfT1(1:np,4)*wfg2(1:np,3)+wfT1(1:np,5)*wfg2(1:np,4))*prefact
    wfg(1:np,3)=(wfT1(1:np,2)*wfg2(1:np,1)-wfT1(1:np,4)*wfg2(1:np,2)+wfT1(1:np,6)*wfg2(1:np,4))*prefact
    wfg(1:np,4)=(wfT1(1:np,3)*wfg2(1:np,1)-wfT1(1:np,5)*wfg2(1:np,2)-wfT1(1:np,6)*wfg2(1:np,3))*prefact
  end subroutine TensorGluontoGluon_Real
  subroutine GluonTensortoGluon(np,wfg1,wfT2,wfg)
    implicit none 
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfg1,wfg
    complex(kind=8),dimension(np,6) :: wfT2
!!$    complex(kind=8),parameter :: prefact=(-0.5d0,0d0)
    complex(kind=8),parameter :: prefact=(0d0,0.5d0)
    wfg(1:np,1)=(-wfg1(1:np,2)*wfT2(1:np,1)-wfg1(1:np,3)*wfT2(1:np,2)-wfg1(1:np,4)*wfT2(1:np,3))*prefact
    wfg(1:np,2)=(-wfg1(1:np,1)*wfT2(1:np,1)-wfg1(1:np,3)*wfT2(1:np,4)-wfg1(1:np,4)*wfT2(1:np,5))*prefact
    wfg(1:np,3)=(-wfg1(1:np,1)*wfT2(1:np,2)+wfg1(1:np,2)*wfT2(1:np,4)-wfg1(1:np,4)*wfT2(1:np,6))*prefact
    wfg(1:np,4)=(-wfg1(1:np,1)*wfT2(1:np,3)+wfg1(1:np,2)*wfT2(1:np,5)+wfg1(1:np,3)*wfT2(1:np,6))*prefact
!!$    do i=1,4
!!$       wfg(i)=-((wfg1(1)*wfT2(1,i)-wfg1(2)*wfT2(2,i)-wfg1(3)*wfT2(3,i)-wfg1(4)*wfT2(4,i))- &
!!$               (wfg1(1)*wfT2(i,1)-wfg1(2)*wfT2(i,2)-wfg1(3)*wfT2(i,3)-wfg1(4)*wfT2(i,4)))*0.25d0
!!$    enddo
  end subroutine GluonTensortoGluon
  subroutine GluonTensortoGluon_real(np,wfg1,wfT2,wfg)
    implicit none 
    integer,intent(in) :: np
    real(kind=8),dimension(np,4) :: wfg1,wfg
    real(kind=8),dimension(np,6) :: wfT2
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1:np,1)=(-wfg1(1:np,2)*wfT2(1:np,1)-wfg1(1:np,3)*wfT2(1:np,2)-wfg1(1:np,4)*wfT2(1:np,3))*prefact
    wfg(1:np,2)=(-wfg1(1:np,1)*wfT2(1:np,1)-wfg1(1:np,3)*wfT2(1:np,4)-wfg1(1:np,4)*wfT2(1:np,5))*prefact
    wfg(1:np,3)=(-wfg1(1:np,1)*wfT2(1:np,2)+wfg1(1:np,2)*wfT2(1:np,4)-wfg1(1:np,4)*wfT2(1:np,6))*prefact
    wfg(1:np,4)=(-wfg1(1:np,1)*wfT2(1:np,3)+wfg1(1:np,2)*wfT2(1:np,5)+wfg1(1:np,3)*wfT2(1:np,6))*prefact
  end subroutine GluonTensortoGluon_Real
  subroutine GluonQuarktoQuark(np,wfg1,wfq2,wfq)
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfg1,wfq2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8),dimension(np) :: TMP1,TMP2,TMP3,TMP4
    TMP1(1:np)=wfg1(1:np,1)+wfg1(1:np,4)
    TMP2(1:np)=wfg1(1:np,1)-wfg1(1:np,4)
    TMP3(1:np)=wfg1(1:np,2)+cImag*wfg1(1:np,3)
    TMP4(1:np)=wfg1(1:np,2)-cImag*wfg1(1:np,3)
    wfq(1:np,1)=prefact*(TMP1(1:np)*wfq2(1:np,3)+TMP3(1:np)*wfq2(1:np,4))
    wfq(1:np,2)=prefact*(TMP2(1:np)*wfq2(1:np,4)+TMP4(1:np)*wfq2(1:np,3))
    wfq(1:np,3)=prefact*(TMP2(1:np)*wfq2(1:np,1)-TMP3(1:np)*wfq2(1:np,2))
    wfq(1:np,4)=prefact*(TMP1(1:np)*wfq2(1:np,2)-TMP4(1:np)*wfq2(1:np,1))
  end subroutine GluonQuarktoQuark
  subroutine QuarkGluontoQuark(np,wfq1,wfg2,wfq)
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfq1,wfg2,wfq
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    complex(kind=8),dimension(np) :: TMP1,TMP2,TMP3,TMP4
    TMP1(1:np)=wfg2(1:np,1)+wfg2(1:np,4)
    TMP2(1:np)=wfg2(1:np,1)-wfg2(1:np,4)
    TMP3(1:np)=wfg2(1:np,2)+cImag*wfg2(1:np,3)
    TMP4(1:np)=wfg2(1:np,2)-cImag*wfg2(1:np,3)
    wfq(1:np,1)=prefact*(TMP1(1:np)*wfq1(1:np,3)+TMP3(1:np)*wfq1(1:np,4))
    wfq(1:np,2)=prefact*(TMP2(1:np)*wfq1(1:np,4)+TMP4(1:np)*wfq1(1:np,3))
    wfq(1:np,3)=prefact*(TMP2(1:np)*wfq1(1:np,1)-TMP3(1:np)*wfq1(1:np,2))
    wfq(1:np,4)=prefact*(TMP1(1:np)*wfq1(1:np,2)-TMP4(1:np)*wfq1(1:np,1))
  end subroutine QuarkGluontoQuark
  subroutine QuarkGluontoQuark_real(np,wfq1,wfg2,wfq)
    implicit none
    integer,intent(in) :: np
    complex(kind=8),dimension(np,4) :: wfq1,wfq
    real(kind=8),dimension(np,4) :: wfg2
    complex(kind=8), parameter :: cImag=(0d0,1d0),prefact=(0d0,1d0)/sqrt(2d0)
    real(kind=8),dimension(np) :: TMP1,TMP2
    complex(kind=8),dimension(np) :: TMP3,TMP4
    TMP1(1:np)=wfg2(1:np,1)+wfg2(1:np,4)
    TMP2(1:np)=wfg2(1:np,1)-wfg2(1:np,4)
    TMP3(1:np)=dcmplx(wfg2(1:np,2),wfg2(1:np,3))
    TMP4(1:np)=dcmplx(wfg2(1:np,2),-wfg2(1:np,3))
    wfq(1:np,1)=prefact*(TMP1(1:np)*wfq1(1:np,3)+TMP3(1:np)*wfq1(1:np,4))
    wfq(1:np,2)=prefact*(TMP2(1:np)*wfq1(1:np,4)+TMP4(1:np)*wfq1(1:np,3))
    wfq(1:np,3)=prefact*(TMP2(1:np)*wfq1(1:np,1)-TMP3(1:np)*wfq1(1:np,2))
    wfq(1:np,4)=prefact*(TMP1(1:np)*wfq1(1:np,2)-TMP4(1:np)*wfq1(1:np,1))
  end subroutine QuarkGluontoQuark_Real
  subroutine GluonPropagator(np,wfg,nhel,p)
    implicit none
    integer,intent(in) :: np
    integer,intent(in) :: nhel
    complex(kind=8),dimension(1:np,1:4,nhel),intent(inout) :: wfg
    real(kind=8),dimension(1:np,0:3),intent(in) :: p
    complex(kind=8),dimension(1:np) :: propagator
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    integer :: ihel,i
    propagator(1:np)=-cImag/(p(1:np,0)**2-p(1:np,1)**2-p(1:np,2)**2-p(1:np,3)**2)
    do ihel=1,nhel
       do i=1,4
          wfg(1:np,i,ihel)=wfg(1:np,i,ihel)*propagator(1:np)
       enddo
    enddo
  end subroutine GluonPropagator
  subroutine GluonPropagator_real(np,wfg,nhel,p)
    implicit none
    integer,intent(in) :: np
    integer,intent(in) :: nhel
    real(kind=8),dimension(1:np,1:4,nhel),intent(inout) :: wfg
    real(kind=8),dimension(1:np,0:3),intent(in) :: p
    real(kind=8),dimension(1:np) :: propagator
    integer :: ihel,i
    propagator(1:np)=1d0/(p(1:np,0)**2-p(1:np,1)**2-p(1:np,2)**2-p(1:np,3)**2)
    do ihel=1,nhel
       do i=1,4
          wfg(1:np,i,ihel)=wfg(1:np,i,ihel)*propagator(1:np)
       enddo
    enddo
  end subroutine GluonPropagator_Real
  subroutine QuarkPropagator(np,wfq,nhel,p)
    implicit none
    integer,intent(in) :: np
    integer,intent(in) :: nhel
    complex(kind=8),dimension(1:np,1:4,nhel),intent(inout) :: wfq
    real(kind=8),dimension(1:np,0:3),intent(in) :: p
    complex(kind=8),dimension(1:np) :: prefact
    complex(kind=8),dimension(1:np,1:4) :: tmp_p,tmp_val
    complex(kind=8),parameter :: cImag=(0d0,1d0)
    integer :: ih
    prefact(1:np)=cImag/(p(1:np,0)**2-p(1:np,1)**2-p(1:np,2)**2-p(1:np,3)**2)
    do ih=1,nhel
       tmp_val(1:np,1:4)=wfq(1:np,1:4,ih)
       tmp_p(1:np,1)=p(1:np,0)+p(1:np,3)
       tmp_p(1:np,2)=p(1:np,0)-p(1:np,3)
       tmp_p(1:np,3)=p(1:np,1)+cImag*p(1:np,2)
       tmp_p(1:np,4)=p(1:np,1)-cImag*p(1:np,2)
       wfq(1:np,1,ih)=(tmp_p(1:np,1)*tmp_val(1:np,3)+tmp_p(1:np,3)*tmp_val(1:np,4))*prefact(1:np)
       wfq(1:np,2,ih)=(tmp_p(1:np,2)*tmp_val(1:np,4)+tmp_p(1:np,4)*tmp_val(1:np,3))*prefact(1:np)
       wfq(1:np,3,ih)=(tmp_p(1:np,2)*tmp_val(1:np,1)-tmp_p(1:np,3)*tmp_val(1:np,2))*prefact(1:np)
       wfq(1:np,4,ih)=(tmp_p(1:np,1)*tmp_val(1:np,2)-tmp_p(1:np,4)*tmp_val(1:np,1))*prefact(1:np)
    enddo
  end subroutine QuarkPropagator
end module FeynmanRules
