program coupling_order_library_regression
  use amp1_1_lib
  implicit none

  integer,parameter :: dp=kind(1d0),n=6,n_amps=6,n_sectors=3
  real(kind=dp),dimension(0:3,n) :: p
  complex(kind=dp),dimension(n_amps) :: expected,actual
  complex(kind=dp),dimension(n_amps,n_sectors) :: expected_by_order,actual_by_order
  real(kind=dp) :: relative_difference

  open(unit=14,file='Library/amp1_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p
  read(14) expected
  read(14) expected_by_order
  close(14)

  call evaluate_amp1_1(p,actual)
  call evaluate_amp1_1_by_order(p,actual_by_order)
  relative_difference=maxval(abs(actual_by_order-expected_by_order)/&
       max(1d-30,abs(actual_by_order),abs(expected_by_order)))
  if (relative_difference.gt.1d-10) then
     write (*,*) 'Generated multi-root library changed a coupling-order sector'
     write (*,*) 'maximum relative difference:',relative_difference
     stop 1
  endif
  relative_difference=maxval(abs(actual-sum(actual_by_order,dim=2))/&
       max(1d-30,abs(actual),abs(sum(actual_by_order,dim=2))))
  if (relative_difference.gt.1d-10) then
     write (*,*) 'Generated multi-root library sectors do not reconstruct the amplitude'
     write (*,*) 'maximum relative difference:',relative_difference
     stop 1
  endif
  write (*,'(a,es12.4)') 'Multi-root coupling-order library regression passed; difference=',&
       relative_difference
end program coupling_order_library_regression
