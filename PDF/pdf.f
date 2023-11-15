      subroutine PDF_eval(ih,ipdgs,x,xmu,pdfs)
!      
! ih = 1 for proton, -1 for anti-proton and 0 for lepton collisions
! ipdgs is a logical array (-6:7) that tells if we should compute the
!     PDF for that parton (with 0=gluon and 7=photon, and the rest
!     anti-quarks).
! x is Bjorken x
! xmu is factorisation scale in GeV (NOT the scale squared!)
! pdfs are the values of the PDFs. Only computed are the ones for which
!     ipdgs(..) is true.
!
      implicit none
      DOUBLE  PRECISION x,xmu,pdfs(-6:7)
      INTEGER IH,ipdg,ipart
      logical ipdgs(-6:7)
      if (ih.eq.0) then
c     Lepton collisions (no PDF). 
         pdfs=1d0
         return
      endif
      ! make sure we have a reasonable Bjorken x.
      if (x.lt.0d0 .or. x.gt.1d0) then
         pdfs=0d0
         return
      endif

      if (ih.ne.1) then
         write (*,*) 'PDFs only for proton collisions (for '/
     $        /'anti-proton collisions needs flipping of ipdgs)'
         stop 1
      endif

c The actual call to the PDFs
      call fdist(ipdgs,x,xmu,pdfs(-6))
      return
      end

      
      subroutine fdist(ipdgs,x,xmu,fx)
      implicit none
      integer ih
      double precision fx(-6:7),x,xmu,nnfx(-6:7)
      logical ipdgs(-6:7)
      fx(-6:7)=0d0
      call NNevolvePDF(ipdgs,x,xmu,nnfx)
      fx(-6:7)=nnfx(-6:7)/x
      return
      end
      
  

      subroutine PDF_initialise
      implicit none
      real*8 asmz
      call NNPDFDriver('NNPDF23nlo_as_0119_qed_mem0.grid')      
      call NNinitPDF(0)
      asmz=0.119d0
      write (*,*) 'NNPDF 23 NLO (alpha_s 0.119) QED set initialised'
      return
      end


      INTEGER FUNCTION NEXTUNOPEN()
C     *****************************************************************
C     ***
C     Returns an unallocated FORTRAN i/o unit.
C     *****************************************************************
C     ***

      LOGICAL EX
C     
      DO 10 N = 10, 300
      INQUIRE (UNIT=N, OPENED=EX)
      IF (.NOT. EX) THEN
        NEXTUNOPEN = N
        RETURN
      ENDIF
 10   CONTINUE
      STOP ' There is no available I/O unit. '
C     *************************
      END

      SUBROUTINE OPENDATA(TABLEFILE)
C     *****************************************************************
C     ***
C     generic subroutine to open the table files in the right
C      directories
C     *****************************************************************
C     ***
      IMPLICIT NONE
C     
      CHARACTER TABLEFILE*(*),UP*3,LIB*4,DIR*4,TEMPNAME*100
      DATA DIR/'PDF/'/
      INTEGER IU,NEXTUNOPEN,I
      EXTERNAL NEXTUNOPEN
      COMMON/IU/IU
      CHARACTER*300 TEMPNAME2, PATH
      CHARACTER*25 UPBUFF
      INTEGER POS, FINE2

      IU=NEXTUNOPEN()
C     Try in the current directory (for cluster use)
 5    TEMPNAME=TABLEFILE
      OPEN(IU,FILE=TEMPNAME,STATUS='old',ERR=10)
      RETURN
C     then try PDF directory
 10   TEMPNAME=DIR//TABLEFILE
      OPEN(IU,FILE=TEMPNAME,STATUS='old',ERR=30)
      RETURN

 30   CONTINUE
      PRINT*,'table for the pdf NOT found  !!!'
      PRINT*,TEMPNAME
      stop 1
      
      RETURN
      END


      
