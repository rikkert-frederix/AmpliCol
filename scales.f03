module scales
  use common
  use particles
  use cuts
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
contains
  subroutine set_scale(choice,n,p,ipdg,scale,status)
    implicit none
    integer,intent(in) :: choice,n
    integer,dimension(n),intent(in) :: ipdg
    real(kind=8),dimension(0:3,n),intent(in) :: p
    real(kind=8),intent(out) :: scale
    integer,intent(out),optional :: status
    real(kind=8),parameter :: momentum_limit=0.125d0*sqrt(huge(1d0))
    real(kind=8),dimension(0:3) :: ptot
    real(kind=8) :: mass_squared,mass_scale,tolerance
    integer :: i
    logical :: found_jet
    scale=0d0
    if (present(status)) status=0
    if (n.lt.2) then
       call scale_failure(-20,'invalid momentum input')
       return
    endif
    if (.not.all(ieee_is_finite(p))) then
       call scale_failure(-20,'invalid momentum input')
       return
    endif
    if (maxval(abs(p)).gt.momentum_limit) then
       call scale_failure(-20,'invalid momentum input')
       return
    endif
    if (choice.eq.0) then
       ! z-mass
       scale=phys_model%get_mass(23)
    elseif (choice.eq.1) then
       ! H_T
       scale=0d0
       do i=3,n
          scale=scale+pT(p(:,i))
       enddo
    elseif (choice.eq.2) then
       ! H_T/2d0
       scale=0d0
       do i=3,n
          scale=scale+pT(p(:,i))
       enddo
       scale=scale/2d0
    elseif (choice.eq.3) then
       ! sqrt(s-hat)
       ptot=p(:,1)+p(:,2)
       call invariant_mass_from_vector(ptot,scale,mass_squared,mass_scale,tolerance)
       if (mass_squared.lt.-tolerance) then
          call scale_failure(-20,'spacelike incoming total')
          return
       endif
    elseif (choice.eq.4) then
       ! min jet-pT
       scale=sqrts
       found_jet=.false.
       do i=3,n
          if (phys_model%is_jet(ipdg(i))) then
             found_jet=.true.
             scale=min(scale,pT(p(:,i)))
          endif
       enddo
       if (.not.found_jet) then
          call scale_failure(-4,'minimum-jet-pT scale has no jet')
          return
       endif
    elseif (choice.eq.5) then
       ! invariant mass non-jet system
       ptot=0d0
       do i=3,n
          if (phys_model%is_jet(ipdg(i))) cycle
          ptot(0:3)=ptot(0:3)+p(0:3,i)
       enddo
       call invariant_mass_from_vector(ptot,scale,mass_squared,mass_scale,tolerance)
       if (mass_squared.lt.-tolerance) then
          call scale_failure(-20,'spacelike non-jet total')
          return
       endif
    else
       call scale_failure(-4,'unknown scale choice')
       return
    endif
    if (.not.ieee_is_finite(scale)) then
       call scale_failure(-20,'non-finite scale')
       return
    endif
    if (scale.le.1d0) then
       call scale_failure(-20,'non-finite or sub-GeV scale')
       return
    endif
  contains
    subroutine invariant_mass_from_vector(momentum,mass,mass2,arithmetic_scale,roundoff)
      real(kind=8),intent(in) :: momentum(0:3)
      real(kind=8),intent(out) :: mass,mass2,arithmetic_scale,roundoff
      mass2=dot(momentum,momentum)
      arithmetic_scale=sum(abs(momentum*momentum))
      roundoff=4096d0*epsilon(1d0)*max(tiny(1d0),arithmetic_scale)
      if (.not.ieee_is_finite(mass2)) then
         mass=0d0
         return
      endif
      if (mass2.ge.-roundoff) then
         mass=sqrt(max(0d0,mass2))
      else
         mass=0d0
      endif
    end subroutine invariant_mass_from_vector

    subroutine scale_failure(code,message)
      integer,intent(in) :: code
      character(len=*),intent(in) :: message
      scale=0d0
      if (present(status)) then
         status=code
      else
         write(*,*) 'ERROR: '//trim(message)//' for scale choice',choice
         stop 1
      endif
    end subroutine scale_failure
  end subroutine set_scale



!C-----------------------------------------------------------------------------
!C
  real(kind=8) FUNCTION ALPHAS_Q(Q,nloop,asmz,status)
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
    DOUBLE PRECISION,intent(in) :: Q,asmz
    INTEGER,intent(in) :: nloop
    INTEGER,intent(out),optional :: status
    DOUBLE PRECISION :: T,AMZ0,AMB,AMC,AMB_NEW,AMC_NEW
    DOUBLE PRECISION :: AS_OUT
    INTEGER NLOOP0,NF3,NF4,NF5
    PARAMETER(NF5=5,NF4=4,NF3=3)
!C
    REAL*8       CMASS,BMASS
    COMMON/QMASS/CMASS,BMASS
    DATA CMASS,BMASS/1.42D0,4.7D0/  ! HEAVY QUARK MASSES FOR THRESHOLDS
!C
!C
    real*8 ZMASS0
    real(kind=8),parameter :: scale_limit=0.125d0*sqrt(huge(1d0))
    real(kind=8),parameter :: coupling_limit=0.125d0*huge(1d0)**0.125d0
    real(kind=8) :: CMASS0,BMASS0
    integer :: newton_status
    SAVE AMZ0,NLOOP0,AMB,AMC,ZMASS0,CMASS0,BMASS0
    DATA AMZ0,NLOOP0,ZMASS0,CMASS0,BMASS0/0D0,0,0D0,0D0,0D0/

    ALPHAS_Q=0d0
    if (present(status)) status=0
    IF (.not.ieee_is_finite(Q)) THEN
       call alphas_failure(-20,'invalid evolution scale')
       return
    ENDIF
    IF (Q.LE.0D0 .or. Q.GT.scale_limit) THEN
       call alphas_failure(-20,'invalid evolution scale')
       return
    ENDIF
    IF (.not.ieee_is_finite(asmz)) THEN
       call alphas_failure(-20,'invalid alpha_s(MZ)')
       return
    ENDIF
    IF (asmz.LE.0D0 .or. asmz.GT.coupling_limit) THEN
       call alphas_failure(-20,'invalid alpha_s(MZ)')
       return
    ENDIF
    IF (NLOOP.LT.1 .or. NLOOP.GT.3) THEN
       call alphas_failure(-4,'unsupported loop order')
       return
    ENDIF
    IF (.not.ieee_is_finite(CMASS) .or. .not.ieee_is_finite(BMASS) .or. &
         .not.ieee_is_finite(z_mass)) THEN
       call alphas_failure(-4,'invalid heavy-quark or Z-mass thresholds')
       return
    ENDIF
    IF (CMASS.LE.0.3D0 .or. BMASS.LE.CMASS .or. z_mass.LE.BMASS) THEN
       call alphas_failure(-4,'invalid heavy-quark or Z-mass thresholds')
       return
    ENDIF
!c--- establish value of coupling at b- and c-mass and save
    IF ((asmz.NE.AMZ0) .OR. (NLOOP.NE.NLOOP0) .OR. (z_mass.NE.ZMASS0) .OR. &
         (CMASS.NE.CMASS0) .OR. (BMASS.NE.BMASS0)) THEN
       T=2D0*DLOG(BMASS/z_mass)
       CALL NEWTON1(T,asmz,AMB_NEW,NLOOP,NF5,newton_status)
       if (newton_status.ne.0) then
          call alphas_failure(newton_status,'failed evolution to the bottom threshold')
          return
       endif
       T=2D0*DLOG(CMASS/BMASS)
       CALL NEWTON1(T,AMB_NEW,AMC_NEW,NLOOP,NF4,newton_status)
       if (newton_status.ne.0) then
          call alphas_failure(newton_status,'failed evolution to the charm threshold')
          return
       endif
       AMB=AMB_NEW
       AMC=AMC_NEW
       AMZ0=asmz
       NLOOP0=NLOOP
       ZMASS0=z_mass
       CMASS0=CMASS
       BMASS0=BMASS
    ENDIF

!c--- evaluate strong coupling at scale q
    IF (Q  .LT. BMASS) THEN
       IF (Q  .LT. CMASS) THEN
          T=2D0*DLOG(Q/CMASS)
          CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3,newton_status)
       ELSE
          T=2D0*DLOG(Q/BMASS)
          CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4,newton_status)
       ENDIF
    ELSE
       T=2D0*DLOG(Q/z_mass)
       CALL NEWTON1(T,asmz,AS_OUT,NLOOP,NF5,newton_status)
    ENDIF
    if (newton_status.ne.0) then
       call alphas_failure(-20,'alpha_s evolution did not produce a finite positive result')
       return
    endif
    if (.not.ieee_is_finite(AS_OUT)) then
       call alphas_failure(-20,'alpha_s evolution did not produce a finite positive result')
       return
    endif
    if (AS_OUT.le.0d0 .or. AS_OUT.gt.coupling_limit) then
       call alphas_failure(-20,'alpha_s evolution did not produce a finite positive result')
       return
    endif
    ALPHAS_Q=AS_OUT
    RETURN
  contains
    subroutine alphas_failure(code,message)
      integer,intent(in) :: code
      character(len=*),intent(in) :: message
      ALPHAS_Q=0d0
      if (present(status)) then
         status=code
      else
         write(6,*) 'ERROR: '//trim(message)//' in alphas_Q'
         stop 1
      endif
    end subroutine alphas_failure
  END FUNCTION ALPHAS_Q


  SUBROUTINE NEWTON1(T,A_IN,A_OUT,NLOOP,NF,STATUS)
!C     Author: R.K. Ellis
!c---  calculate a_out using nloop beta-function evolution 
!c---  with nf flavours, given starting value as-in
!c---  given as_in and logarithmic separation between 
!c---  input scale and output scale t.
!c---  Evolution is performed using Newton's method,
!c---  with a precision given by tol.
    IMPLICIT NONE
    INTEGER,intent(in) :: NLOOP,NF
    REAL*8,intent(in) :: T,A_IN
    REAL*8,intent(out) :: A_OUT
    INTEGER,intent(out),optional :: STATUS
    REAL*8 :: AS,F,FP,DELTA,DENOM,ONE_LOOP,STEP,CANDIDATE
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
    real*8,parameter :: coupling_floor=sqrt(tiny(1d0))
    real*8,parameter :: coupling_limit=0.125d0*huge(1d0)**0.125d0
    integer,parameter :: max_iterations=100
    integer :: iteration
    logical :: candidate_is_valid

    A_OUT=0d0
    if (present(STATUS)) STATUS=0
    if (.not.ieee_is_finite(T) .or. .not.ieee_is_finite(A_IN)) then
       call newton_failure(-20,'invalid input')
       return
    endif
    if (A_IN.lt.coupling_floor .or. A_IN.gt.coupling_limit .or. &
         NLOOP.lt.1 .or. NLOOP.gt.3 .or. NF.lt.3 .or. NF.gt.5) then
       call newton_failure(-20,'invalid input')
       return
    endif
    DENOM=1D0+A_IN*B0(NF)*T
    if (.not.ieee_is_finite(DENOM)) then
       call newton_failure(-20,'scale is below the perturbative domain')
       return
    endif
    if (DENOM.le.0d0) then
       call newton_failure(-20,'scale is below the perturbative domain')
       return
    endif
    A_OUT=A_IN/DENOM
    if (NLOOP.eq.1) then
       if (.not.valid_coupling(A_OUT)) call newton_failure(-20,'invalid one-loop result')
       return
    endif
    ONE_LOOP=DENOM
    DENOM=1D0+B0(NF)*A_IN*T+C1(NF)*A_IN*LOG(ONE_LOOP)
    if (.not.ieee_is_finite(DENOM)) then
       call newton_failure(-20,'invalid higher-loop starting value')
       return
    endif
    if (DENOM.le.0d0) then
       call newton_failure(-20,'invalid higher-loop starting value')
       return
    endif
    A_OUT=A_IN/DENOM
    if (.not.valid_coupling(A_OUT)) then
       call newton_failure(-20,'invalid higher-loop starting value')
       return
    endif
    do iteration=1,max_iterations
       AS=A_OUT
       IF (NLOOP.eq.2) THEN
          F=B0(NF)*T+beta_integral2(A_IN,NF)-beta_integral2(AS,NF)
          FP=1D0/(AS**2*(1D0+C1(NF)*AS))
       ELSE
          F=B0(NF)*T+beta_integral3(A_IN,NF)-beta_integral3(AS,NF)
          FP=1D0/(AS**2*(1D0+C1(NF)*AS+C2(NF)*AS**2))
       ENDIF
       if (.not.ieee_is_finite(F) .or. .not.ieee_is_finite(FP)) then
          call newton_failure(-20,'non-finite Newton residual')
          return
       endif
       if (FP.eq.0d0) then
          call newton_failure(-20,'zero Newton derivative')
          return
       endif
       STEP=F/FP
       if (.not.ieee_is_finite(STEP)) then
          call newton_failure(-20,'non-finite Newton step')
          return
       endif
       DELTA=ABS(STEP/AS)
       CANDIDATE=AS-STEP
       candidate_is_valid=valid_coupling(CANDIDATE)
       if (.not.candidate_is_valid) CANDIDATE=0.5d0*AS
       A_OUT=CANDIDATE
       if (candidate_is_valid) then
          if (ieee_is_finite(DELTA)) then
             if (DELTA.le.TOL) return
          endif
       endif
    enddo
    call newton_failure(-20,'Newton iteration did not converge')
    RETURN
  contains
    pure real(kind=8) function beta_integral2(value,nf_local)
      real(kind=8),intent(in) :: value
      integer,intent(in) :: nf_local
      beta_integral2=1D0/value+C1(nf_local)*&
           LOG((C1(nf_local)*value)/(1D0+C1(nf_local)*value))
    end function beta_integral2

    pure real(kind=8) function beta_integral3(value,nf_local)
      real(kind=8),intent(in) :: value
      integer,intent(in) :: nf_local
      beta_integral3=1D0/value+0.5D0*C1(nf_local)*&
           LOG((C2(nf_local)*value**2)/(1D0+C1(nf_local)*value+C2(nf_local)*value**2))-&
           (C1(nf_local)**2-2D0*C2(nf_local))/DEL(nf_local)*&
           ATAN((2D0*C2(nf_local)*value+C1(nf_local))/DEL(nf_local))
    end function beta_integral3

    pure logical function valid_coupling(value)
      real(kind=8),intent(in) :: value
      valid_coupling=.false.
      if (.not.ieee_is_finite(value)) return
      valid_coupling=value.ge.coupling_floor .and. value.le.coupling_limit
    end function valid_coupling

    subroutine newton_failure(code,message)
      integer,intent(in) :: code
      character(len=*),intent(in) :: message
      A_OUT=0d0
      if (present(STATUS)) then
         STATUS=code
      else
         write(6,*) 'ERROR: '//trim(message)//' in NEWTON1'
         stop 1
      endif
    end subroutine newton_failure
  END SUBROUTINE NEWTON1
end module scales
