program mixed_spinor_library_regression
  use amp1_1_lib, only: evaluate_amp1_1
  use amp2_1_lib, only: evaluate_amp2_1
  implicit none
  integer,parameter :: dp=kind(1d0),n=4
  real(kind=dp),dimension(0:3,n) :: p
  complex(kind=dp),dimension(1) :: expected,actual

  open(unit=14,file='Library/amp1_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p
  read(14) expected
  close(14)
  call evaluate_amp1_1(p,actual)
  call compare_amplitudes('b g > t W-',expected(1),actual(1))

  open(unit=14,file='Library/amp2_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p
  read(14) expected
  close(14)
  call evaluate_amp2_1(p,actual)
  call compare_amplitudes('b~ g > t~ W+',expected(1),actual(1))
  write (*,'(a)') 'Mixed-spinor generated-library regression passed'

contains

  subroutine compare_amplitudes(label,dynamic,generated)
    implicit none
    character(len=*),intent(in) :: label
    complex(kind=dp),intent(in) :: dynamic,generated
    real(kind=dp) :: relative_difference

    relative_difference=abs(generated-dynamic)/&
         max(1d-30,abs(generated)+abs(dynamic))
    if (abs(dynamic).lt.1d-12 .or. relative_difference.gt.1d-11) then
       write (*,*) trim(label),': generated amplitude differs from dynamic'
       write (*,*) 'dynamic/generated:',dynamic,generated
       write (*,*) 'relative difference:',relative_difference
       stop 1
    endif
    write (*,'(a,a,a,es12.4)') 'LIBRARY_MIXED_SPINOR[',trim(label),&
         '] relative difference=',relative_difference
  end subroutine compare_amplitudes

end program mixed_spinor_library_regression
