program coupling_order_pruning_library_regression
  use amp1_1_lib
  implicit none

  integer,parameter :: dp=kind(1d0),n=6,n_amps=1,n_sectors=3
  real(kind=dp),dimension(0:3,n) :: p
  complex(kind=dp),dimension(n_amps) :: expected,actual
  complex(kind=dp),dimension(n_amps,n_sectors) :: expected_by_order,actual_by_order

  open(unit=14,file='Library/amp1_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p
  read(14) expected
  read(14) expected_by_order
  close(14)
  if (any(expected.ne.cmplx(0d0,0d0,kind=dp)) .or.&
       any(expected_by_order.ne.cmplx(0d0,0d0,kind=dp))) then
     write (*,*) 'Empty pruned library stored a nonzero reference coefficient'
     stop 1
  endif

  call check_point(p)
  call check_point(1.219d0*p)
  write (*,'(a)') 'Empty coupling-order pruning library regression passed'

contains

  subroutine check_point(momentum)
    implicit none
    real(kind=dp),dimension(0:3,n),intent(in) :: momentum

    actual=cmplx(7d0,-3d0,kind=dp)
    actual_by_order=cmplx(-5d0,2d0,kind=dp)
    call evaluate_amp1_1(momentum,actual)
    call evaluate_amp1_1_by_order(momentum,actual_by_order)
    if (any(actual.ne.cmplx(0d0,0d0,kind=dp)) .or.&
         any(actual_by_order.ne.cmplx(0d0,0d0,kind=dp))) then
       write (*,*) 'Generated locally empty amplitude did not return exact zero'
       stop 1
    endif
  end subroutine check_point

end program coupling_order_pruning_library_regression
