      subroutine PDF_eval(ih,ipdgs,x,xmu,pdfs)
      use pdf_internal_interface, only: fdist
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
      DOUBLE PRECISION,INTENT(IN) :: x,xmu
      DOUBLE PRECISION,INTENT(OUT) :: pdfs(-6:7)
      INTEGER,INTENT(IN) :: IH
      LOGICAL,INTENT(IN) :: ipdgs(-6:7)
      pdfs=0d0
      if (ih.eq.0) then
c     Lepton collisions (no PDF). 
         pdfs=1d0
         return
      endif
      ! make sure we have a reasonable Bjorken x.
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(xmu)) then
         return
      endif
      if (x.le.0d0 .or. x.gt.1d0 .or. xmu.le.0d0
     $     .or. xmu.gt.dsqrt(huge(1d0))) then
         return
      endif

      if (ih.ne.1) then
         write (*,*) 'PDFs only for proton collisions (for '/
     $        /'anti-proton collisions needs flipping of ipdgs)'
         stop 1
      endif

c The actual call to the PDFs
      call fdist(ipdgs,x,xmu,pdfs(-6:7))
      return
      end

      
      subroutine fdist(ipdgs,x,xmu,fx)
      use pdf_internal_interface, only: NNevolvePDF
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      double precision,intent(out) :: fx(-6:7)
      double precision,intent(in) :: x,xmu
      double precision nnfx(-6:7)
      integer i
      logical,intent(in) :: ipdgs(-6:7)
      fx(-6:7)=0d0
      if (.not.ieee_is_finite(x) .or. .not.ieee_is_finite(xmu))
     $     return
      if (x.le.0d0 .or. x.gt.1d0 .or. xmu.le.0d0
     $     .or. xmu.gt.dsqrt(huge(1d0))) return
      call NNevolvePDF(ipdgs,x,xmu,nnfx)
      if (.not.all(ieee_is_finite(nnfx))) then
         write (*,*) 'Internal PDF returned a non-finite density'
         stop 1
      endif
      do i=-6,7
         if (dabs(nnfx(i)).gt.huge(1d0)*x) then
            fx(i)=dsign(huge(1d0),nnfx(i))
         else
            fx(i)=nnfx(i)/x
         endif
      enddo
      return
      end
      
  

      subroutine PDF_initialise(gridfilename,irep)
      use pdf_internal_interface, only: NNPDFDriver,NNinitPDF
      implicit none
      character(len=*),intent(in) :: gridfilename
      integer,intent(in) :: irep
      call NNPDFDriver(gridfilename)
      call NNinitPDF(irep)
      write (*,*) 'Internal PDF grid initialised: ',gridfilename
      write (99,*) 'Internal PDF grid initialised: ',gridfilename
      return
      end


      INTEGER FUNCTION NEXTUNOPEN()
C     *****************************************************************
C     ***
C     Returns an unallocated FORTRAN i/o unit.
C     *****************************************************************
C     ***

      LOGICAL EX
      INTEGER N,IOS
C     
      NEXTUNOPEN = -1
      DO 10 N = 10, 300
      INQUIRE (UNIT=N, OPENED=EX, IOSTAT=IOS)
      IF (IOS .NE. 0) GOTO 10
      IF (.NOT. EX) THEN
        NEXTUNOPEN = N
        RETURN
      ENDIF
 10   CONTINUE
      WRITE (*,*) 'There is no available legacy PDF I/O unit (10:300)'
      STOP 1
C     *************************
      END

      SUBROUTINE OPENDATA(TABLEFILE)
      USE PDF_INTERNAL_INTERFACE, ONLY: NEXTUNOPEN
C     *****************************************************************
C     ***
C     generic subroutine to open the table files in the right
C      directories
C     *****************************************************************
C     ***
      IMPLICIT NONE
C     
      CHARACTER TABLEFILE*(*),DIR*4
      CHARACTER(LEN=LEN(TABLEFILE)+4) TEMPNAME
      CHARACTER*256 IOMESSAGE
      DATA DIR/'PDF/'/
      INTEGER IU,IOS
      COMMON/IU/IU

      IU=NEXTUNOPEN()
C     Try in the current directory (for cluster use)
      TEMPNAME=TRIM(TABLEFILE)
      OPEN(UNIT=IU,FILE=TRIM(TEMPNAME),STATUS='old',ACTION='read',
     $     IOSTAT=IOS,IOMSG=IOMESSAGE)
      IF (IOS .NE. 0) GOTO 10
      RETURN
C     then try PDF directory
 10   TEMPNAME=DIR//TABLEFILE
      OPEN(UNIT=IU,FILE=TRIM(TEMPNAME),STATUS='old',ACTION='read',
     $     IOSTAT=IOS,IOMSG=IOMESSAGE)
      IF (IOS .NE. 0) GOTO 30
      RETURN

 30   CONTINUE
      PRINT*,'PDF table not found: ',TRIM(TABLEFILE)
      PRINT*,'Last path tried: ',TRIM(TEMPNAME)
      PRINT*,'I/O error: ',TRIM(IOMESSAGE)
      stop 1
      
      RETURN
      END


      
