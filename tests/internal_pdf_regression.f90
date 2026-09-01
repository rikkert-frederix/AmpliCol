program internal_pdf_regression
  use pdf_internal_interface, only: PDF_initialise,PDF_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  character(len=1024) :: grid_filename,member_text
  logical :: requested(-6:7)
  real(kind=8) :: pdfs(-6:7)
  integer :: member,ios

  grid_filename='PDF/NNPDF23nlo_as_0119_qed_mem0.grid'
  member=0
  if (command_argument_count().ge.1) call get_command_argument(1,grid_filename)
  if (command_argument_count().ge.2) then
     call get_command_argument(2,member_text)
     read(member_text,*,iostat=ios) member
     if (ios.ne.0) error stop 'invalid member argument in internal-PDF regression'
  endif

  ! The legacy driver writes initialization details to the run-log unit.
  open(unit=99,status='scratch',action='write')
  call PDF_initialise(trim(grid_filename),member)

  requested=.true.
  call PDF_eval(1,requested,0.1d0,91.2d0,pdfs)
  if (.not.all(ieee_is_finite(pdfs))) error stop 'internal PDF returned non-finite values'
  if (all(pdfs.eq.0d0)) error stop 'internal PDF returned an empty table'

  requested=.false.
  requested(0)=.true.
  call PDF_eval(1,requested,0.2d0,20d0,pdfs)
  if (.not.ieee_is_finite(pdfs(0))) error stop 'internal gluon PDF is non-finite'
  if (any(pdfs(-6:-1).ne.0d0) .or. any(pdfs(1:7).ne.0d0)) &
       error stop 'internal PDF evaluated unrequested flavours'

  call PDF_eval(1,requested,0d0,20d0,pdfs)
  if (any(pdfs.ne.0d0)) error stop 'out-of-support internal PDF was not zeroed'
  call PDF_eval(1,requested,ieee_value(0d0,ieee_quiet_nan),20d0,pdfs)
  if (any(pdfs.ne.0d0)) error stop 'NaN-x internal PDF request was not zeroed'
  call PDF_eval(1,requested,0.2d0,ieee_value(0d0,ieee_quiet_nan),pdfs)
  if (any(pdfs.ne.0d0)) error stop 'NaN-scale internal PDF request was not zeroed'
  close(99)
  write(*,'(a)') 'Internal-PDF regression: PASS'
end program internal_pdf_regression
