! Minimal LHAGlue-compatible entry points for builds using the bundled PDF grid.
! The ordinary build continues to use LHAPDF. Select this adapter explicitly
! with `make PDF_BACKEND=internal ...`; the standalone color probe also links
! it so that an otherwise-unused PDF symbol does not impose LHAPDF on the probe.

subroutine InitPDFsetbyname(setname)
  implicit none
  character(len=*),intent(in) :: setname
  call PDF_initialise
end subroutine InitPDFsetbyname

subroutine initPDF(member)
  implicit none
  integer,intent(in) :: member
end subroutine initPDF

subroutine setlhaparm(parameter)
  implicit none
  character(len=*),intent(in) :: parameter
end subroutine setlhaparm

subroutine evolvePDF(x,q,fxq)
  implicit none
  real(kind=8),intent(in) :: x,q
  real(kind=8),intent(out) :: fxq(-6:7)
  logical :: ipdgs(-6:7)
  ipdgs=.true.
  call PDF_eval(1,ipdgs,x,q,fxq)
  fxq=fxq*x
end subroutine evolvePDF

real(kind=8) function alphaspdf(q)
  use scales, only: alphas_Q
  implicit none
  real(kind=8),intent(in) :: q
  alphaspdf=alphas_Q(q,2,0.119d0)
end function alphaspdf
