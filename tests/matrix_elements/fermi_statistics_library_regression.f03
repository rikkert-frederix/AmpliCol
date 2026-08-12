program fermi_statistics_library_regression
  use amp1_1_lib
  implicit none
  integer,parameter :: dp=kind(1d0)
  real(kind=dp),dimension(0:3,6) :: p
  complex(kind=dp),dimension(1) :: expected,actual
  complex(kind=dp),dimension(1,1) :: expected_by_order,actual_by_order
  real(kind=dp) :: relative_difference

  open(unit=14,file='Library/amp1_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p
  read(14) expected
  read(14) expected_by_order
  close(14)
  call evaluate_amp1_1(p,actual)
  call evaluate_amp1_1_by_order(p,actual_by_order)
  relative_difference=abs(actual(1)-expected(1))/&
       max(1d-30,abs(actual(1))+abs(expected(1)))
  if (relative_difference.gt.1d-10) then
     write (*,*) 'Generated library dropped a Fermi-statistics sign'
     write (*,*) 'expected/actual:',expected(1),actual(1)
     stop 1
  endif
  relative_difference=abs(actual_by_order(1,1)-expected_by_order(1,1))/&
       max(1d-30,abs(actual_by_order(1,1))+abs(expected_by_order(1,1)))
  if (relative_difference.gt.1d-10 .or.&
       abs(actual(1)-actual_by_order(1,1)).gt.1d-10*max(1d-30,abs(actual(1)))) then
     write (*,*) 'Generated library dropped a coupling-order sector'
     write (*,*) 'expected/actual:',expected_by_order(1,1),actual_by_order(1,1)
     stop 1
  endif
  write (*,'(a,es12.4)') 'Generated-library regression passed; relative difference=',&
       relative_difference
end program fermi_statistics_library_regression
