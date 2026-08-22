program heft_library_regression
  use amp1_1_lib, only: evaluate_amp1_1
  use amp2_1_lib, only: evaluate_amp2_1
  use amp3_1_lib, only: evaluate_amp3_1
  implicit none
  integer,parameter :: dp=kind(1d0)
  real(kind=dp),parameter :: gs=1.2177157847767197d0
  real(kind=dp),parameter :: gew=sqrt(8d0*acos(-1d0)/132.507d0)
  real(kind=dp),parameter :: gheft=5.248659993424408d-5
  real(kind=dp),dimension(0:3,5) :: p5
  real(kind=dp),dimension(0:3,4) :: p4
  real(kind=dp),dimension(0:3,3) :: p3
  complex(kind=dp),dimension(1) :: expected,generated

  open(unit=14,file='Library/amp1_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p5
  read(14) expected
  close(14)
  call evaluate_amp1_1(p5,generated,gs,gew,gheft)
  call assert_close('generated gg > ggh',generated(1),expected(1))
  call assert_close('fixed gg > ggh coefficient',generated(1),&
       cmplx(-1.12707793396816619d-4,-2.27460135830594135d-4,kind=dp))

  open(unit=14,file='Library/amp2_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p4
  read(14) expected
  close(14)
  call evaluate_amp2_1(p4,generated,gs,gew,gheft)
  call assert_close('generated ttbar > gh',generated(1),expected(1))
  call assert_close('fixed ttbar > gh interference coefficient',generated(1),&
       cmplx(1.55455111506466086d-3,-1.0426671423761962d0,kind=dp))

  open(unit=14,file='Library/amp3_1_lib.data',form='unformatted',&
       access='stream',status='old',action='read')
  read(14) p3
  read(14) expected
  close(14)
  call evaluate_amp3_1(p3,generated,gs,gew,gheft)
  call assert_close('generated scalar-terminal gg > h',generated(1),expected(1))
  call assert_close('fixed scalar-terminal gg > h',generated(1),&
       cmplx(0d0,-7812.5d0*gheft,kind=dp))

  write (*,'(a)') 'HEFT generated-library regression passed'

contains

  subroutine assert_close(label,value,reference)
    implicit none
    character(len=*),intent(in) :: label
    complex(kind=dp),intent(in) :: value,reference
    real(kind=dp) :: difference,scale

    difference=abs(value-reference)
    scale=max(1d-30,abs(value)+abs(reference))
    if (difference.gt.1d-11*scale .and. difference.gt.1d-13) then
       write (*,*) trim(label),' mismatch:',value,reference
       write (*,*) 'absolute/relative difference:',difference,difference/scale
       stop 1
    endif
  end subroutine assert_close

end program heft_library_regression
