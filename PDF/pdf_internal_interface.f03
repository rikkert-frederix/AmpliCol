module pdf_internal_state
  implicit none
  private
  integer,parameter,public :: max_pdf_members=100,max_pdf_flavours=14
  integer,parameter,public :: max_pdf_x_points=100,max_pdf_q2_points=60
  integer,public :: nfl=0,nx=0,nq2=0,mem=-1,rep=-1
  real(kind=8),public :: alphas=0d0
  real(kind=8),public :: xgrid(max_pdf_x_points),logxgrid(max_pdf_x_points)
  real(kind=8),public :: q2grid(max_pdf_q2_points),logq2grid(max_pdf_q2_points)
  real(kind=8),public :: pdfgrid(0:max_pdf_members,max_pdf_flavours,&
       max_pdf_x_points,max_pdf_q2_points)
  logical,public :: hasphoton=.false.,pdf_grid_loaded=.false.
end module pdf_internal_state

module pdf_internal_interface
  implicit none
  private
  public :: PDF_eval,PDF_initialise,fdist,NNPDFDriver,NNinitPDF,NNevolvePDF
  public :: NEXTUNOPEN,OPENDATA,readPDFSet,lh_polin2,lh_polint

  interface
     subroutine PDF_eval(ih,requested,x,scale,pdfs)
       implicit none
       integer,intent(in) :: ih
       logical,intent(in) :: requested(-6:7)
       real(kind=8),intent(in) :: x,scale
       real(kind=8),intent(out) :: pdfs(-6:7)
     end subroutine PDF_eval

     subroutine PDF_initialise(grid_filename,member)
       implicit none
       character(len=*),intent(in) :: grid_filename
       integer,intent(in) :: member
     end subroutine PDF_initialise

     subroutine fdist(requested,x,scale,pdfs)
       implicit none
       logical,intent(in) :: requested(-6:7)
       real(kind=8),intent(in) :: x,scale
       real(kind=8),intent(out) :: pdfs(-6:7)
     end subroutine fdist

     subroutine NNPDFDriver(grid_filename)
       implicit none
       character(len=*),intent(in) :: grid_filename
     end subroutine NNPDFDriver

     subroutine NNinitPDF(member)
       implicit none
       integer,intent(in) :: member
     end subroutine NNinitPDF

     subroutine NNevolvePDF(requested,x,scale,pdfs)
       implicit none
       logical,intent(in) :: requested(-6:7)
       real(kind=8),intent(in) :: x,scale
       real(kind=8),intent(out) :: pdfs(-6:7)
     end subroutine NNevolvePDF

     integer function NEXTUNOPEN()
       implicit none
     end function NEXTUNOPEN

     subroutine OPENDATA(table_file)
       implicit none
       character(len=*),intent(in) :: table_file
     end subroutine OPENDATA

     subroutine readPDFSet(grid_filename)
       implicit none
       character(len=*),intent(in) :: grid_filename
     end subroutine readPDFSet

     subroutine lh_polin2(x1a,x2a,values,m,n,x1,x2,value,error_estimate)
       implicit none
       integer,intent(in) :: m,n
       real(kind=8),intent(in) :: x1a(m),x2a(n),values(m,n),x1,x2
       real(kind=8),intent(out) :: value,error_estimate
     end subroutine lh_polin2

     subroutine lh_polint(xa,values,n,x,value,error_estimate)
       implicit none
       integer,intent(in) :: n
       real(kind=8),intent(in) :: xa(n),values(n),x
       real(kind=8),intent(out) :: value,error_estimate
     end subroutine lh_polint
  end interface
end module pdf_internal_interface
