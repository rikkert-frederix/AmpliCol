program test_integration_histograms
  use integration_histograms
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  implicit none
  integer :: nintegrals(2)
  real(kind=8) :: value,error

  nintegrals=[1,1]
  call histogram_initialize(nintegrals,.true.,'tests/integration_histograms.HwU')
  call histogram_book(1,'correlated observable [pb/bin]',2,0d0,2d0)
  call histogram_book(2,'mapped observable [pb/bin]',2,0d0,2d0)

  ! Leaf 1, point 1: a real event and its counterevent must be squared only
  ! after cancellation, giving a point value of one rather than two
  ! independent contributions of ten and minus nine.
  call histogram_begin_point(1,1)
  call histogram_fill(1,0.5d0,10d0,4d0)
  call histogram_fill(1,0.5d0,-9d0,0d0)
  call histogram_fill(2,0.5d0,10d0,0d0)
  call histogram_fill(2,1.5d0,-9d0,0d0)
  call histogram_commit_point()

  call histogram_begin_point(1,1)
  call histogram_fill(1,0.5d0,3d0,2d0)
  call histogram_commit_point()

  ! Leaf 2 has a different number of points.  Its sample mean must be added
  ! to leaf 1's sample mean, not folded into one global point average.
  call histogram_begin_point(2,1)
  call histogram_fill(1,0.5d0,5d0,0d0)
  call histogram_commit_point()
  call histogram_finalize_iteration()

  call histogram_get_bin(1,1,.true.,value,error)
  call assert_close(value,7d0,1d-12,'stratified NLO mean')
  call assert_close(error,sqrt(0.5d0),1d-12,'correlated NLO uncertainty')

  call histogram_get_bin(1,1,.false.,value,error)
  call assert_close(value,3d0,1d-12,'Born mean')
  call assert_close(error,sqrt(0.5d0),1d-12,'Born uncertainty')

  call histogram_get_bin(1,2,.true.,value,error)
  call assert_close(value,0d0,1d-12,'empty bin value')
  call assert_close(error,0d0,1d-12,'empty bin uncertainty')

  call histogram_get_bin(2,1,.true.,value,error)
  call assert_close(value,5d0,1d-12,'real observable bin')
  call assert_close(error,sqrt(12.5d0),1d-12,'real observable bin uncertainty')
  call histogram_get_bin(2,2,.true.,value,error)
  call assert_close(value,-4.5d0,1d-12,'mapped counterevent bin')
  call assert_close(error,sqrt(10.125d0),1d-12,'mapped counterevent bin uncertainty')

  ! A second iteration has different statistics in each leaf.  Combining it
  ! must reproduce the moments of the complete per-leaf samples.
  call histogram_begin_point(1,1)
  call histogram_fill(1,0.5d0,4d0,5d0)
  call histogram_commit_point()

  call histogram_begin_point(2,1)
  call histogram_fill(1,0.5d0,1d0,0d0)
  call histogram_commit_point()
  call histogram_begin_point(2,1)
  call histogram_fill(1,0.5d0,3d0,0d0)
  call histogram_commit_point()
  call histogram_finalize_iteration()

  call histogram_get_bin(1,1,.true.,value,error)
  call assert_close(value,17d0/3d0,1d-12,'combined-iteration NLO mean')
  call assert_close(error,sqrt(38d0/27d0),1d-12,'combined-iteration NLO uncertainty')

  call histogram_get_bin(1,1,.false.,value,error)
  call assert_close(value,11d0/3d0,1d-12,'combined-iteration Born mean')
  call assert_close(error,sqrt(14d0/27d0),1d-12,'combined-iteration Born uncertainty')

  ! One invalid fill discards the complete correlated point, including valid
  ! counterevents that were staged before the failure.
  call histogram_initialize([1],.true.,'tests/integration_histograms.HwU')
  call histogram_book(1,'invalid-point guard [pb/bin]',1,0d0,1d0)
  call histogram_begin_point(1,1)
  call histogram_fill(1,0.5d0,7d0,3d0)
  call histogram_fill(1,0.5d0,ieee_value(0d0,ieee_quiet_nan),0d0)
  call histogram_commit_point()
  call histogram_finalize_iteration()
  call histogram_get_bin(1,1,.true.,value,error)
  call assert_close(value,0d0,1d-12,'invalid point is discarded atomically')
  if (.not.ieee_is_finite(error)) then
     write(*,*) 'FAIL: invalid histogram point contaminated uncertainty'
     stop 1
  endif

  write(*,*) 'integration histogram tests passed'

contains

  subroutine assert_close(actual,expected,tolerance,label)
    real(kind=8),intent(in) :: actual,expected,tolerance
    character(len=*),intent(in) :: label
    if (abs(actual-expected).gt.tolerance) then
       write(*,*) 'FAIL: ',trim(label),actual,expected
       stop 1
    endif
  end subroutine assert_close

end program test_integration_histograms
