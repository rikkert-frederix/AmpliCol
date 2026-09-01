program random_number_regression
  use random_number_interface, only: ntuple,ranmar,rmarin
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer(kind=8),parameter :: maximum_seed=904866561_8
  integer(kind=8) :: iseed
  common /to_seed/ iseed
  character(len=32) :: mode
  real(kind=8) :: value
  real(kind=8),dimension(8) :: first,second
  integer :: i

  mode='max'
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  open(unit=99,status='scratch',action='write')
  select case (trim(mode))
  case ('max')
     iseed=maximum_seed
     call ntuple(value,0d0,1d0,1)
     if (iseed.ne.maximum_seed) error stop 'RNG initialization changed the shared seed'
     if (.not.ieee_is_finite(value) .or. value.le.0d0 .or. value.ge.1d0) &
          error stop 'maximum-seed RNG value is outside (0,1)'
     write(*,'(a,1x,es24.16,1x,i0)') 'RNG_MAX',value,iseed
  case ('repeat')
     call rmarin(1802,9373)
     do i=1,size(first)
        call ranmar(first(i))
     enddo
     call rmarin(1802,9373)
     do i=1,size(second)
        call ranmar(second(i))
     enddo
     if (any(first.ne.second)) error stop 'RNG reinitialization is not reproducible'
     if (any(first.lt.0d0) .or. any(first.ge.1d0)) &
          error stop 'RNG value is outside [0,1)'
     write(*,'(a)') 'RNG_REPEAT PASS'
  case ('uninitialized')
     call ranmar(value)
     error stop 'uninitialized RANMAR unexpectedly succeeded'
  case ('invalid-seed')
     call rmarin(31329,0)
     error stop 'invalid RANMAR seed unexpectedly succeeded'
  case ('wide-interval')
     iseed=1_8
     call ntuple(value,-huge(1d0),huge(1d0),1)
     error stop 'overflowing RNG interval unexpectedly succeeded'
  case default
     error stop 'unknown random-number regression mode'
  end select
  close(99)
end program random_number_regression
