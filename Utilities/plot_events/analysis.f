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
      integer nexternal,ipdg(nexternal)
      double precision p(0:3,nexternal),evt_wgt_LC(1),evt_wgt_NLC(1)
     $     ,evt_wgt_full(1),www(1)

      ! weight distribution
      www(1)=1d0
      call HwU_fill(1,evt_wgt_full(1)/evt_wgt_LC(1),www)
      call HwU_fill(2,evt_wgt_NLC(1)/evt_wgt_LC(1),www)
      call HwU_fill(3,evt_wgt_full(1)/evt_wgt_NLC(1),www)
      
      ! total rates
      call HwU_fill(4,0.5d0,evt_wgt_LC)
      call HwU_fill(4,1.5d0,evt_wgt_NLC)
      call HwU_fill(4,2.5d0,evt_wgt_full)

      call HwU_add_points
      return
      end

