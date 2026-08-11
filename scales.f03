module scales
  use common
  use particles
  use cuts
contains
  subroutine set_scale(choice,n,p,ipdg,scale)
    implicit none
    integer,intent(in) :: choice,n
    integer,dimension(n),intent(in) :: ipdg
    real(kind=8),dimension(0:3,n),intent(in) :: p
    real(kind=8),intent(out) :: scale
    real(kind=8),dimension(0:3) :: ptot
    integer :: i
    if (choice.eq.0) then
       ! z-mass
       scale=phys_model%get_mass(23)
    elseif (choice.eq.1) then
       ! H_T
       scale=0d0
       do i=3,n
          scale=scale+sqrt(max(0d0,(p(0,i)+p(3,i))*(p(0,i)-p(3,i))))
       enddo
    elseif (choice.eq.2) then
       ! H_T/2d0
       scale=0d0
       do i=3,n
          scale=scale+sqrt(max(0d0,(p(0,i)+p(3,i))*(p(0,i)-p(3,i))))
       enddo
       scale=scale/2d0
    elseif (choice.eq.3) then
       ! sqrt(s-hat)
       scale=sqrt(sumdot(p(0,1),p(0,2)))
    elseif (choice.eq.4) then
       ! min jet-pT
       scale=sqrts
       do i=3,n
          if (phys_model%is_jet(ipdg(i))) scale=min(scale,pT(p(0,i)))
       enddo
       if (scale.eq.sqrts) then
          write (*,*) 'Found no jet. Scale choice not valid',choice
       endif
    elseif (choice.eq.5) then
       ! invariant mass non-jet system
       ptot=0d0
       do i=3,n
          if (phys_model%is_jet(ipdg(i))) cycle
          ptot(0:3)=ptot(0:3)+p(0:3,i)
       enddo
       scale=sqrt(dot(ptot,ptot))
    endif
    if (scale.le.1d0) then
       write (*,*) 'Found scale smaller than 1 GeV. Scale choice not valid',choice
       stop 1
    endif
  end subroutine set_scale



!C-----------------------------------------------------------------------------
!C
  real(kind=8) FUNCTION ALPHAS_Q(Q,nloop,asmz)
!c
!c     Evaluation of strong coupling constant alpha_S
!c     Author: R.K. Ellis
!c
!c     q -- scale at which alpha_s is to be evaluated
!c
!c-- common block alfas.inc
!c     asmz -- value of alpha_s at the mass of the Z-boson
!c     nloop -- the number of loops (1,2, or 3) at which beta 
!c
!c     function is evaluated to determine running.
!c     the values of the cmass and the bmass should be set
!c     in common block qmass.
!C-----------------------------------------------------------------------------
    IMPLICIT NONE
!c
    DOUBLE PRECISION Q,T,AMZ0,AMB,AMC
    DOUBLE PRECISION AS_OUT
    INTEGER NLOOP0,NF3,NF4,NF5
    PARAMETER(NF5=5,NF4=4,NF3=3)
!C
    REAL*8       CMASS,BMASS
    COMMON/QMASS/CMASS,BMASS
    DATA CMASS,BMASS/1.42D0,4.7D0/  ! HEAVY QUARK MASSES FOR THRESHOLDS
!C
!C
    real*8 ZMASS0
    SAVE AMZ0,NLOOP0,AMB,AMC,ZMASS0
    DATA AMZ0,NLOOP0,ZMASS0/0D0,0,0D0/
    integer nloop
    real*8 asmz

    IF (Q .LE. 0D0) THEN 
       WRITE(6,*) 'q .le. 0 in alphas'
       WRITE(6,*) 'q= ',Q
       STOP
    ENDIF
    IF (asmz .LE. 0D0) THEN 
       WRITE(6,*) 'asmz .le. 0 in alphas',asmz
       STOP
       asmz=0.1185D0
    ENDIF
    IF (CMASS .LE. 0.3D0) THEN 
       WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
       STOP
       CMASS=1.42D0
    ENDIF
    IF (BMASS .LE. 0D0) THEN 
       WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
       WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
       STOP
       BMASS=4.7D0
    ENDIF
!c--- establish value of coupling at b- and c-mass and save
    IF ((asmz .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0) .OR. (z_mass .NE. ZMASS0)) THEN
       AMZ0=asmz
       NLOOP0=NLOOP
       ZMASS0=z_mass
       T=2D0*DLOG(BMASS/z_mass)
       CALL NEWTON1(T,asmz,AMB,NLOOP,NF5)
       T=2D0*DLOG(CMASS/BMASS)
       CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
    ENDIF

!c--- evaluate strong coupling at scale q
    IF (Q  .LT. BMASS) THEN
       IF (Q  .LT. CMASS) THEN
          T=2D0*DLOG(Q/CMASS)
          CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
       ELSE
          T=2D0*DLOG(Q/BMASS)
          CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
       ENDIF
    ELSE
       T=2D0*DLOG(Q/z_mass)
       CALL NEWTON1(T,asmz,AS_OUT,NLOOP,NF5)
    ENDIF
    ALPHAS_Q=AS_OUT
    RETURN
  END FUNCTION ALPHAS_Q


  SUBROUTINE NEWTON1(T,A_IN,A_OUT,NLOOP,NF)
!C     Author: R.K. Ellis
!c---  calculate a_out using nloop beta-function evolution 
!c---  with nf flavours, given starting value as-in
!c---  given as_in and logarithmic separation between 
!c---  input scale and output scale t.
!c---  Evolution is performed using Newton's method,
!c---  with a precision given by tol.
    IMPLICIT NONE
    INTEGER :: NLOOP,NF
    REAL*8 :: T,A_IN,A_OUT,AS,F2,F3,F,FP,DELTA
!C---     B0=(11.-2.*NF/3.)/4./PI
    REAL*8,dimension(3:5),parameter :: B0=[0.716197243913527D0,0.66314559621623D0,0.61009394851893D0] 
!C---     C1=(102.D0-38.D0/3.D0*NF)/4.D0/PI/(11.D0-2.D0/3.D0*NF)
    REAL*8,dimension(3:5),parameter :: C1=[.565884242104515D0,0.49019722472304D0,0.40134724779695D0]
!C---     C2=(2857.D0/2.D0-5033*NF/18.D0+325*NF**2/54)
!C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*NF)
    REAL*8,dimension(3:5),parameter :: C2=[0.453013579178645D0,0.30879037953664D0,0.14942733137107D0]
!C---     DEL=SQRT(4*C2-C1**2)
    REAL*8,dimension(3:5),parameter :: DEL=[1.22140465909230D0,0.99743079911360D0,0.66077962451190D0]
    real*8,parameter :: TOL=5d-4
    F2(AS)=1D0/AS+C1(NF)*LOG((C1(NF)*AS)/(1D0+C1(NF)*AS))
    F3(AS)=1D0/AS+0.5D0*C1(NF) &
         & *LOG((C2(NF)*AS**2)/(1D0+C1(NF)*AS+C2(NF)*AS**2)) &
         & -(C1(NF)**2-2D0*C2(NF))/DEL(NF) &
         & *ATAN((2D0*C2(NF)*AS+C1(NF))/DEL(NF))
    A_OUT=A_IN/(1D0+A_IN*B0(NF)*T)
    IF (NLOOP .EQ. 1) RETURN
    A_OUT=A_IN/(1D0+B0(NF)*A_IN*T+C1(NF)*A_IN*LOG(1D0+A_IN*B0(NF)*T))
    IF (A_OUT .LT. 0D0) AS=0.3D0
30  AS=A_OUT
    IF (NLOOP .EQ. 2) THEN
       F=B0(NF)*T+F2(A_IN)-F2(AS)
       FP=1D0/(AS**2*(1D0+C1(NF)*AS))
    ENDIF
    IF (NLOOP .EQ. 3) THEN
       F=B0(NF)*T+F3(A_IN)-F3(AS)
       FP=1D0/(AS**2*(1D0+C1(NF)*AS+C2(NF)*AS**2))
    ENDIF
    A_OUT=AS-F/FP
    DELTA=ABS(F/FP/AS)
    IF (DELTA .GT. TOL) GO TO 30
    RETURN
  END SUBROUTINE NEWTON1
end module scales
