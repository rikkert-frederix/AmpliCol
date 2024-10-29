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
      character*8 HwUtype(3)
      data HwUtype/'|T@LC   ','|T@NLC  ','|T@full '/
      double precision pi,tiny
      parameter (pi=3.14159265358979312d0,tiny=1d-6)
c      
      call HwU_inithist(nwgt,weights_info)
c
      call HwU_book(1,'weights LC to full ',500,0d0-tiny,5d0-tiny)
      call HwU_book(2,'weights LC to NLC  ',500,0d0-tiny,5d0-tiny)
      call HwU_book(3,'weights NLC to full',500,0d0-tiny,5d0-tiny)
      
      call HwU_book(4,'rates',3,0d0,3d0)

      call HwU_book(5,'subprocesses',100,0.5d0,100.5d0)

      call HwU_book(6,'pt j1',100,25d0,125d0)
      call HwU_book(7,'pt j2',100,25d0,125d0)
      call HwU_book(8,'pt j3',100,25d0,125d0)
      call HwU_book(9,'pt j4',100,25d0,125d0)

      call HwU_book(10,'eta j1',70,-7d0,+7d0)
      call HwU_book(11,'eta j2',70,-7d0,+7d0)
      call HwU_book(12,'eta j3',70,-7d0,+7d0)
      call HwU_book(13,'eta j4',70,-7d0,+7d0)
      
      call HwU_book(14,'dr j1j2',100,0d0,15d0)
      call HwU_book(15,'dr j1j3',100,0d0,15d0)
      call HwU_book(16,'dr j1j4',100,0d0,15d0)
      call HwU_book(17,'dr j2j3',100,0d0,15d0)
      call HwU_book(18,'dr j2j4',100,0d0,15d0)
      call HwU_book(19,'dr j3j4',100,0d0,15d0)

      call HwU_book(20,'m j1j2',100,0d0,100d0)
      call HwU_book(21,'m j1j3',100,0d0,100d0)
      call HwU_book(22,'m j1j4',100,0d0,100d0)
      call HwU_book(23,'m j2j3',100,0d0,100d0)
      call HwU_book(24,'m j2j4',100,0d0,100d0)
      call HwU_book(25,'m j3j4',100,0d0,100d0)

      call HwU_book(26,'m j1j2 l',100,0d0,7500d0)
      call HwU_book(27,'m j1j3 l',100,0d0,7500d0)
      call HwU_book(28,'m j1j4 l',100,0d0,7500d0)
      call HwU_book(29,'m j2j3 l',100,0d0,7500d0)
      call HwU_book(30,'m j2j4 l',100,0d0,7500d0)
      call HwU_book(31,'m j3j4 l',100,0d0,7500d0)

      call HwU_book(32,'pt j1 l',100,0d0,1500d0)
      call HwU_book(33,'pt j2 l',100,0d0,1500d0)
      call HwU_book(34,'pt j3 l',100,0d0,1500d0)
      call HwU_book(35,'pt j4 l',100,0d0,1500d0)

      call HwU_book(36,'pt j1 nqq0',100,25d0,125d0)
      call HwU_book(37,'pt j2 nqq0',100,25d0,125d0)
      call HwU_book(38,'pt j3 nqq0',100,25d0,125d0)
      call HwU_book(39,'pt j4 nqq0',100,25d0,125d0)

      call HwU_book(40,'eta j1 nqq0',70,-7d0,+7d0)
      call HwU_book(41,'eta j2 nqq0',70,-7d0,+7d0)
      call HwU_book(42,'eta j3 nqq0',70,-7d0,+7d0)
      call HwU_book(43,'eta j4 nqq0',70,-7d0,+7d0)
      
      call HwU_book(44,'dr j1j2 nqq0',100,0d0,15d0)
      call HwU_book(45,'dr j1j3 nqq0',100,0d0,15d0)
      call HwU_book(46,'dr j1j4 nqq0',100,0d0,15d0)
      call HwU_book(47,'dr j2j3 nqq0',100,0d0,15d0)
      call HwU_book(48,'dr j2j4 nqq0',100,0d0,15d0)
      call HwU_book(49,'dr j3j4 nqq0',100,0d0,15d0)

      call HwU_book(50,'m j1j2 nqq0',100,0d0,100d0)
      call HwU_book(51,'m j1j3 nqq0',100,0d0,100d0)
      call HwU_book(52,'m j1j4 nqq0',100,0d0,100d0)
      call HwU_book(53,'m j2j3 nqq0',100,0d0,100d0)
      call HwU_book(54,'m j2j4 nqq0',100,0d0,100d0)
      call HwU_book(55,'m j3j4 nqq0',100,0d0,100d0)

      call HwU_book(56,'m j1j2 l nqq0',100,0d0,7500d0)
      call HwU_book(57,'m j1j3 l nqq0',100,0d0,7500d0)
      call HwU_book(58,'m j1j4 l nqq0',100,0d0,7500d0)
      call HwU_book(59,'m j2j3 l nqq0',100,0d0,7500d0)
      call HwU_book(60,'m j2j4 l nqq0',100,0d0,7500d0)
      call HwU_book(61,'m j3j4 l nqq0',100,0d0,7500d0)

      call HwU_book(62,'pt j1 l nqq0',100,0d0,1500d0)
      call HwU_book(63,'pt j2 l nqq0',100,0d0,1500d0)
      call HwU_book(64,'pt j3 l nqq0',100,0d0,1500d0)
      call HwU_book(65,'pt j4 l nqq0',100,0d0,1500d0)

      
      call HwU_book(66,'pt j1 nqq1',100,25d0,125d0)
      call HwU_book(67,'pt j2 nqq1',100,25d0,125d0)
      call HwU_book(68,'pt j3 nqq1',100,25d0,125d0)
      call HwU_book(69,'pt j4 nqq1',100,25d0,125d0)

      call HwU_book(70,'eta j1 nqq1',70,-7d0,+7d0)
      call HwU_book(71,'eta j2 nqq1',70,-7d0,+7d0)
      call HwU_book(72,'eta j3 nqq1',70,-7d0,+7d0)
      call HwU_book(73,'eta j4 nqq1',70,-7d0,+7d0)
      
      call HwU_book(74,'dr j1j2 nqq1',100,0d0,15d0)
      call HwU_book(75,'dr j1j3 nqq1',100,0d0,15d0)
      call HwU_book(76,'dr j1j4 nqq1',100,0d0,15d0)
      call HwU_book(77,'dr j2j3 nqq1',100,0d0,15d0)
      call HwU_book(78,'dr j2j4 nqq1',100,0d0,15d0)
      call HwU_book(79,'dr j3j4 nqq1',100,0d0,15d0)

      call HwU_book(80,'m j1j2 nqq1',100,0d0,100d0)
      call HwU_book(81,'m j1j3 nqq1',100,0d0,100d0)
      call HwU_book(82,'m j1j4 nqq1',100,0d0,100d0)
      call HwU_book(83,'m j2j3 nqq1',100,0d0,100d0)
      call HwU_book(84,'m j2j4 nqq1',100,0d0,100d0)
      call HwU_book(85,'m j3j4 nqq1',100,0d0,100d0)

      call HwU_book(86,'m j1j2 l nqq1',100,0d0,7500d0)
      call HwU_book(87,'m j1j3 l nqq1',100,0d0,7500d0)
      call HwU_book(88,'m j1j4 l nqq1',100,0d0,7500d0)
      call HwU_book(89,'m j2j3 l nqq1',100,0d0,7500d0)
      call HwU_book(90,'m j2j4 l nqq1',100,0d0,7500d0)
      call HwU_book(91,'m j3j4 l nqq1',100,0d0,7500d0)

      call HwU_book(92,'pt j1 l nqq1',100,0d0,1500d0)
      call HwU_book(93,'pt j2 l nqq1',100,0d0,1500d0)
      call HwU_book(94,'pt j3 l nqq1',100,0d0,1500d0)
      call HwU_book(95,'pt j4 l nqq1',100,0d0,1500d0)
      
      call HwU_book(96,'pt j1 nqq2',100,25d0,125d0)
      call HwU_book(97,'pt j2 nqq2',100,25d0,125d0)
      call HwU_book(98,'pt j3 nqq2',100,25d0,125d0)
      call HwU_book(99,'pt j4 nqq2',100,25d0,125d0)

      call HwU_book(100,'eta j1 nqq2',70,-7d0,+7d0)
      call HwU_book(101,'eta j2 nqq2',70,-7d0,+7d0)
      call HwU_book(102,'eta j3 nqq2',70,-7d0,+7d0)
      call HwU_book(103,'eta j4 nqq2',70,-7d0,+7d0)
     
      call HwU_book(104,'dr j1j2 nqq2',100,0d0,15d0)
      call HwU_book(105,'dr j1j3 nqq2',100,0d0,15d0)
      call HwU_book(106,'dr j1j4 nqq2',100,0d0,15d0)
      call HwU_book(107,'dr j2j3 nqq2',100,0d0,15d0)
      call HwU_book(108,'dr j2j4 nqq2',100,0d0,15d0)
      call HwU_book(109,'dr j3j4 nqq2',100,0d0,15d0)

      call HwU_book(110,'m j1j2 nqq2',100,0d0,100d0)
      call HwU_book(111,'m j1j3 nqq2',100,0d0,100d0)
      call HwU_book(112,'m j1j4 nqq2',100,0d0,100d0)
      call HwU_book(113,'m j2j3 nqq2',100,0d0,100d0)
      call HwU_book(114,'m j2j4 nqq2',100,0d0,100d0)
      call HwU_book(115,'m j3j4 nqq2',100,0d0,100d0)
      
      call HwU_book(116,'m j1j2 l nqq2',100,0d0,7500d0)
      call HwU_book(117,'m j1j3 l nqq2',100,0d0,7500d0)
      call HwU_book(118,'m j1j4 l nqq2',100,0d0,7500d0)
      call HwU_book(119,'m j2j3 l nqq2',100,0d0,7500d0)
      call HwU_book(120,'m j2j4 l nqq2',100,0d0,7500d0)
      call HwU_book(121,'m j3j4 l nqq2',100,0d0,7500d0)

      call HwU_book(122,'pt j1 l nqq2',100,0d0,1500d0)
      call HwU_book(123,'pt j2 l nqq2',100,0d0,1500d0)
      call HwU_book(124,'pt j3 l nqq2',100,0d0,1500d0)
      call HwU_book(125,'pt j4 l nqq2',100,0d0,1500d0)
      
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
      character*140 outfile
      common /c_to_analysis/outfile
      open (unit=99,file=outfile,status='unknown')
      xnorm=1d0
      call HwU_output(99,xnorm)
      close (99)
      return
      end subroutine HwU_write_file
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(nexternal,p,ipdg,evt_wgt_LC,evt_wgt_NLC
     $     ,evt_wgt_full)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nexternal,ipdg(1:nexternal)
      double precision p(0:3,nexternal),evt_wgt_LC(1),evt_wgt_NLC(1)
     $     ,evt_wgt_full(1),www(1),ptj1,ptj2,ptj3,ptj4,eta1,eta2,eta3
     $     ,eta4,dr12,dr13,dr14,dr23,dr24,dr34,m12,m13,m14,m23,m24,m34
      double precision pjet(0:3,nexternal-2),max
      integer nqq,iflav(2),i,imax,j,l
      logical filled(3:nexternal)
      double precision pt,eta,dr,m
      external pt,eta,dr,m
      ! weight distribution
      www(1)=1d0
      call HwU_fill(1,evt_wgt_full(1)/evt_wgt_LC(1),www)
      call HwU_fill(2,evt_wgt_NLC(1)/evt_wgt_LC(1),www)
      call HwU_fill(3,evt_wgt_full(1)/evt_wgt_NLC(1),www)
      
      ! total rates
      call HwU_fill(4,0.5d0,evt_wgt_LC)
      call HwU_fill(4,1.5d0,evt_wgt_NLC)
      call HwU_fill(4,2.5d0,evt_wgt_full)

      iflav(1:2)=0
      nqq=0
      do i=1,nexternal
         if (abs(ipdg(i)).ge.1 .and. abs(ipdg(i)).lt.6) then
            nqq=nqq+1
            if (iflav(1).eq.0) then
               iflav(1)=abs(ipdg(i))
            elseif (iflav(1).ne.abs(ipdg(i))) then
               iflav(2)=abs(ipdg(i))
            endif
         endif
      enddo
      nqq=nqq/2

      ! order jets in pT
      filled(3:nexternal)=.false.
      do i=1,nexternal-2
         max=0d0
         do j=3,nexternal
            if (filled(j)) cycle
            if (pt(p(0,j)).gt.max) then
               max=pt(p(0,j))
               imax=j
            endif
         enddo
         pjet(0:3,i)=p(0:3,imax)
         filled(imax)=.true.
      enddo
         
      ptj1=pt(pjet(0,1))
      ptj2=pt(pjet(0,2))
      eta1=eta(pjet(0,1))
      eta2=eta(pjet(0,2))
      dr12=dr(pjet(0,1),pjet(0,2))
      m12=m(pjet(0,1),pjet(0,2))
      if (nexternal.ge.5) then
         ptj3=pt(pjet(0,3))
         eta3=eta(pjet(0,3))
         dr13=dr(pjet(0,1),pjet(0,3))
         dr23=dr(pjet(0,2),pjet(0,3))
         m13=m(pjet(0,1),pjet(0,3))
         m23=m(pjet(0,2),pjet(0,3))
      else
         ptj3=-1d0
         eta3=-100d0
         dr13=-1d0
         dr23=-1d0
         m13=-1d0
         m23=-1d0
      endif
      if (nexternal.ge.6) then
         ptj4=pt(pjet(0,4))
         eta4=eta(pjet(0,4))
         dr14=dr(pjet(0,1),pjet(0,4))
         dr24=dr(pjet(0,2),pjet(0,4))
         dr34=dr(pjet(0,3),pjet(0,4))
         m14=m(pjet(0,1),pjet(0,4))
         m24=m(pjet(0,2),pjet(0,4))
         m34=m(pjet(0,3),pjet(0,4))
      else
         ptj4=-1d0
         eta4=-100d0
         dr14=-1d0
         dr24=-1d0
         dr34=-1d0
         m14=-1d0
         m24=-1d0
         m34=-1d0
      endif
      
      
! subprocesses
      if (nqq.eq.0) then
! gg -> n g
         call HwU_fill(5,1d0,evt_wgt_full)
      elseif (nqq.eq.1) then
         if (ipdg(1).eq.21 .and. ipdg(2).eq.21) then
! gg -> qqbar +g
            call HwU_fill(5,3d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! qqbar -> +g
            call HwU_fill(5,4d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! qbarq -> +g
            call HwU_fill(5,5d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).eq.21)) then
! qg -> q +g
            call HwU_fill(5,6d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! gq -> q +g
            call HwU_fill(5,7d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).eq.21)) then
! qbarg -> qbar +g
            call HwU_fill(5,8d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! gbarq -> qbar +g
            call HwU_fill(5,9d0,evt_wgt_full)
         else
            write (*,*) 'unknown process #1',ipdg
            write (*,*) nqq,iflav
            stop 1
         endif
      elseif(nqq.eq.2 .and. iflav(1).gt.0 .and. iflav(2).gt.0
     $                              .and. iflav(1).ne.iflav(2)) then
         if (ipdg(1).eq.21 .and. ipdg(2).eq.21) then
! gg -> qqbar QQbar +g
            call HwU_fill(5,11d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6) .and.
     $           (ipdg(1).eq.-ipdg(2))) then
! qqbar -> QQbar +g
            call HwU_fill(5,12d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6) .and.
     $           (ipdg(1).eq.-ipdg(2))) then
! qbarq -> QQbar +g
            call HwU_fill(5,13d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).eq.21)) then
! qg -> q QQbar +g
            call HwU_fill(5,14d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! gq -> q QQbar +g
            call HwU_fill(5,15d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).eq.21)) then
! qbarg -> qbar QQbar +g
            call HwU_fill(5,16d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! gqbar -> qbar QQbar +g
            call HwU_fill(5,17d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6) .and.
     $           (ipdg(1).ne.ipdg(2))) then
! qQ -> qQ +g
            call HwU_fill(5,18d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6) .and.
     $           (ipdg(1).ne.-ipdg(2))) then
! qbarQ -> q Qbar +g
            call HwU_fill(5,19d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6) .and.
     $           (ipdg(1).ne.-ipdg(2))) then
! qQbar -> q Qbar +g
            call HwU_fill(5,20d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6) .and.
     $           (ipdg(1).ne.ipdg(2))) then
! qbarQbar -> qbar Qbar +g
            call HwU_fill(5,21d0,evt_wgt_full)
         else
            write (*,*) 'unknown process #2',ipdg
            write (*,*) nqq,iflav
            stop 1
         endif
      elseif(nqq.eq.2) then
         if (ipdg(1).eq.21 .and. ipdg(2).eq.21) then
! gg -> qqbar qqbar +g
            call HwU_fill(5,23d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! qqbar -> qqbar +g
            call HwU_fill(5,24d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! qbarq -> qqbar +g
            call HwU_fill(5,25d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).eq.21)) then
! qg -> q qqbar +g
            call HwU_fill(5,26d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! gq -> q qqbar +g
            call HwU_fill(5,27d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).eq.21)) then
! qbarg -> qbar qqbar +g
            call HwU_fill(5,28d0,evt_wgt_full)
         elseif ((ipdg(1).eq.21) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! gqbar -> qbar qqbar +g
            call HwU_fill(5,29d0,evt_wgt_full)
         elseif ((ipdg(1).ge.1 .and. ipdg(1).le.6) .and.
     $           (ipdg(2).ge.1 .and. ipdg(2).le.6)) then
! qq -> q q +g
            call HwU_fill(5,30d0,evt_wgt_full)
         elseif ((ipdg(1).le.-1 .and. ipdg(1).ge.-6) .and.
     $           (ipdg(2).le.-1 .and. ipdg(2).ge.-6)) then
! qbarqbar -> qbar qbar +g
            call HwU_fill(5,31d0,evt_wgt_full)
         else
            write (*,*) 'unknown process #3',ipdg
            write (*,*) nqq,iflav
            stop 1
         endif
      endif
      if (nqq.eq.1) then
         call HwU_fill(5,dble(32+iflav(1)),evt_wgt_full)
      elseif (nqq.eq.2 .and. iflav(2).eq.0) then
         call HwU_fill(5,dble(38+iflav(1)),evt_wgt_full)
      endif


      do i=1,4
         if (i.eq.2 .and. nqq.ne.0) cycle
         if (i.eq.3 .and. nqq.ne.1) cycle
         if (i.eq.4 .and. nqq.ne.2) cycle
         l=(i-1)*30
         call HwU_fill(l+6,ptj1,evt_wgt_full)
         call HwU_fill(l+7,ptj2,evt_wgt_full)
         call HwU_fill(l+8,ptj3,evt_wgt_full)
         call HwU_fill(l+9,ptj4,evt_wgt_full)

         call HwU_fill(l+10,eta1,evt_wgt_full)
         call HwU_fill(l+11,eta2,evt_wgt_full)
         call HwU_fill(l+12,eta3,evt_wgt_full)
         call HwU_fill(l+13,eta4,evt_wgt_full)
         
         call HwU_fill(l+14,dr12,evt_wgt_full)
         call HwU_fill(l+15,dr13,evt_wgt_full)
         call HwU_fill(l+16,dr14,evt_wgt_full)
         call HwU_fill(l+17,dr23,evt_wgt_full)
         call HwU_fill(l+18,dr24,evt_wgt_full)
         call HwU_fill(l+19,dr34,evt_wgt_full)

         call HwU_fill(l+20,m12,evt_wgt_full)
         call HwU_fill(l+21,m13,evt_wgt_full)
         call HwU_fill(l+22,m14,evt_wgt_full)
         call HwU_fill(l+23,m23,evt_wgt_full)
         call HwU_fill(l+24,m24,evt_wgt_full)
         call HwU_fill(l+25,m34,evt_wgt_full)
         
         call HwU_fill(l+26,m12,evt_wgt_full)
         call HwU_fill(l+27,m13,evt_wgt_full)
         call HwU_fill(l+28,m14,evt_wgt_full)
         call HwU_fill(l+29,m23,evt_wgt_full)
         call HwU_fill(l+30,m24,evt_wgt_full)
         call HwU_fill(l+31,m34,evt_wgt_full)

         call HwU_fill(l+32,ptj1,evt_wgt_full)
         call HwU_fill(l+33,ptj2,evt_wgt_full)
         call HwU_fill(l+34,ptj3,evt_wgt_full)
         call HwU_fill(l+35,ptj4,evt_wgt_full)
      enddo
      
      call HwU_add_points
      return
      end

      double precision function pt(p1)
      implicit none
      double precision p1(0:3)
      pt=sqrt(p1(1)**2+p1(2)**2)
      end
      double precision function eta(p1)
      implicit none
      double precision p1(0:3),theta
      theta=acos(p1(3)/sqrt(p1(1)**2+p1(2)**2+p1(3)**2))
      eta=-log(dtan(theta/2d0))
      end
      double precision function dr(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),delta_phi,eta
      external delta_phi,eta
      dr=sqrt(delta_phi(p1,p2)**2+(eta(p1)-eta(p2))**2)
      end
      double precision function delta_phi(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3),denom,pt
      external pt
      denom=pt(p1)*pt(p2)
      delta_phi=acos((p1(1)*p2(1)+p1(2)*p2(2))/denom)
      end
      double precision function m(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      m=sqrt((p1(0)+p2(0))**2-(p1(1)+p2(1))**2-
     $       (p1(2)+p2(2))**2-(p1(3)+p2(3))**2)
      end
