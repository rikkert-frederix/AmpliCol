c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nwgt
      character*(*) weights_info(*)
      integer j,kk,l
      character*3 bb(5)
      data bb/'inc','020','050','100','250'/
      double precision pi
      parameter (pi=3.14159265358979312d0)
c      
      call HwU_inithist(nwgt,weights_info)
c
      call HwU_book(1,'rates',3,0d0,3d0)
      call HwU_book(2,'rates ratios',3,0d0,3d0)
      
      call HwU_book(3,'ratio NLC/LC narrow',100,0.7d0,1.4d0)
      call HwU_book(4,'ratio full/LC narrow',100,0.7d0,1.4d0)
      call HwU_book(5,'ratio full/NLC narrow',100,0.95d0,1.05d0)

      call HwU_book(6,'ratio log10(NLC/LC) narrow',100,-0.5d0,0.5d0)
      call HwU_book(7,'ratio log10(full/LC) narrow',100,-0.5d0,0.5d0)
      call HwU_book(8,'ratio log10(full/NLC) narrow',100,-0.05d0,0.05d0)

      call HwU_book(9,'ratio NLC/LC',100,-25d0,25d0)
      call HwU_book(10,'ratio full/LC',100,-25d0,25d0)
      call HwU_book(11,'ratio full/NLC',100,-25d0,25d0)

      call HwU_book(12,'ratio log10(NLC/LC)',100,-2d0,2d0)
      call HwU_book(13,'ratio log10(full/LC)',100,-2d0,2d0)
      call HwU_book(14,'ratio log10(full/NLC)',100,-2d0,2d0)

      call HwU_book(15,'ratio log10(NLC/LC) wide',100,-6d0,6d0)
      call HwU_book(16,'ratio log10(full/LC) wide',100,-6d0,6d0)
      call HwU_book(17,'ratio log10(full/NLC) wide',100,-6d0,6d0)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision dummy
      call HwU_write_file
      return                
      end

      subroutine HwU_write_file
      implicit none
      double precision :: xnorm
      open (unit=99,file='events.HwU',status='unknown')
      xnorm=1d0
      call HwU_output(99,xnorm)
      close (99)
      return
      end subroutine HwU_write_file
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(nexternal,p,evt_wgt_LC,evt_wgt_NLC
     $     ,evt_wgt_full)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nexternal
      double precision p(0:3,nexternal),evt_wgt_LC(1),evt_wgt_NLC(1)
     $     ,evt_wgt_full(1),www(1)

      ! total rates
      call HwU_fill(1,0.5d0,evt_wgt_LC)
      call HwU_fill(1,1.5d0,evt_wgt_NLC)
      call HwU_fill(1,2.5d0,evt_wgt_full)

!      call HwU_fill(2,0.5d0,evt_wgt_NLC/evt_wgt_full)
!      call HwU_fill(2,1.5d0,evt_wgt_NLC/evt_wgt_LC)
!      call HwU_fill(2,2.5d0,evt_wgt_LC/evt_wgt_full)

      call HwU_fill(2,0.5d0,evt_wgt_LC)
      call HwU_fill(2,1.5d0,evt_wgt_NLC)
      call HwU_fill(2,2.5d0,evt_wgt_full)

      ! weight distribution
      www(1)=evt_wgt_LC(1)
      call HwU_fill(3,evt_wgt_NLC(1)/evt_wgt_LC(1),www)
      call HwU_fill(4,evt_wgt_full(1)/evt_wgt_LC(1),www)
      call HwU_fill(5,evt_wgt_full(1)/evt_wgt_NLC(1),www)
      if (evt_wgt_NLC(1).gt.0d0) call HwU_fill(6,log10(evt_wgt_NLC(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(7,log10(evt_wgt_full(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(8,log10(evt_wgt_full(1)
     $     /evt_wgt_NLC(1)),www)
      call HwU_fill(9,evt_wgt_NLC(1)/evt_wgt_LC(1),www)
      call HwU_fill(10,evt_wgt_full(1)/evt_wgt_LC(1),www)
      call HwU_fill(11,evt_wgt_full(1)/evt_wgt_NLC(1),www)
      if (evt_wgt_NLC(1).gt.0d0) call HwU_fill(12,log10(evt_wgt_NLC(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(13,log10(evt_wgt_full(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(14,log10(evt_wgt_full(1)
     $     /evt_wgt_NLC(1)),www)
      
      if (evt_wgt_NLC(1).gt.0d0) call HwU_fill(15,log10(evt_wgt_NLC(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(16,log10(evt_wgt_full(1)
     $     /evt_wgt_LC(1)),www)
      if (evt_wgt_full(1).gt.0d0) call HwU_fill(17,log10(evt_wgt_full(1)
     $     /evt_wgt_NLC(1)),www)
      
      call HwU_add_points
      return
      end


      function getpseudorap(p)
      implicit none
      real*8 getpseudorap,en,ptx,pty,pl,tiny,pt,eta,th,p(0:3)
      parameter (tiny=1.d-5)
c
      en=p(0)
      ptx=p(1)
      pty=p(2)
      pl=p(3)
      pt=sqrt(ptx**2+pty**2)
      if(pt.lt.tiny.and.abs(pl).lt.tiny)then
        eta=sign(1.d0,pl)*1.d8
      else
        th=atan2(pt,pl)
        eta=-log(tan(th/2.d0))
      endif
      getpseudorap=eta
      return
      end

      function getpt(p)
      implicit none
      real*8 getpt,p(0:3)
      getpt=dsqrt(p(1)**2+p(2)**2)
      return
      end

      double precision function getdelphi(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),denom,temp
      DENOM = SQRT(P1(1)**2 + P1(2)**2) * SQRT(P2(1)**2 + P2(2)**2)
      TEMP = MAX(-0.99999999D0, (P1(1)*P2(1) + P1(2)*P2(2)) / DENOM)
      TEMP = MIN( 0.99999999D0, TEMP)
      getdelphi = ACOS(TEMP)
      return
      end

      double precision function getdelR(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision getdelphi,getpseudorap
      external getdelphi,getpseudorap
      getdelR = sqrt(max((getdelphi(P1,P2))**2+(getpseudorap(p1)
     $     -getpseudorap(p2))**2,0d0))
      return
      end

      double precision function getinvm(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),p(0:3)
      integer i
      do i=0,3
         p(i)=p1(i)+p2(i)
      enddo
      getinvm=sqrt(max(p(0)**2-p(1)**2-p(2)**2-p(3)**2,0d0))
      return
      end
      
      double precision function getinvm3(p1,p2,p3)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p(0:3)
      integer i
      do i=0,3
         p(i)=p1(i)+p2(i)+p3(i)
      enddo
      getinvm3=sqrt(max(p(0)**2-p(1)**2-p(2)**2-p(3)**2,0d0))
      return
      end
      
