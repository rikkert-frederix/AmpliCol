program math_functions_regression
  use math_functions
  implicit none
  character(len=32) :: mode
  integer :: value,index
  integer(kind=8) :: value8
  integer,dimension(3) :: permutation,next_permutation

  mode='success'
  value=0
  value8=0_8
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  select case (trim(mode))
  case ('success')
     if (checked_multiply(0,huge(1),'zero-left').ne.0 .or. &
          checked_multiply(huge(1),0,'zero-right').ne.0) &
          error stop 'zero-factor multiplication failed'
     if (checked_multiply8(0_8,huge(1_8),'zero-left-8').ne.0_8 .or. &
          checked_multiply8(huge(1_8),0_8,'zero-right-8').ne.0_8) &
          error stop 'zero-factor 64-bit multiplication failed'
     if (checked_add(huge(1)-1,1,'boundary-add').ne.huge(1)) &
          error stop 'default-integer boundary addition failed'
     if (checked_add8(huge(1_8)-1_8,1_8,'boundary-add-8').ne.huge(1_8)) &
          error stop '64-bit boundary addition failed'
     if (factorial(12).ne.479001600) error stop 'default factorial failed'
     if (factorial8(20).ne.2432902008176640000_8) error stop '64-bit factorial failed'
     if (checked_integer_power(3,10,'small power').ne.59049) &
          error stop 'integer power failed'
     if (checked_integer_power(0,huge(1),'zero power base').ne.0 .or. &
          checked_integer_power(1,huge(1),'unit power base').ne.1) &
          error stop 'constant-base power did not handle a boundary exponent'
     if (permutation_count(huge(1),0).ne.1) &
          error stop 'zero-length boundary permutation count failed'
     if (permutation_count(4,3).ne.24) error stop 'partial permutation count failed'
     permutation=[1,2,3]
     do index=1,24
        if (any(permutation.lt.1) .or. any(permutation.gt.4) .or. &
             count(permutation.eq.permutation(1)).ne.1 .or. &
             count(permutation.eq.permutation(2)).ne.1 .or. &
             count(permutation.eq.permutation(3)).ne.1) &
             error stop 'partial permutation successor is invalid'
        if (index.lt.24) then
           call get_next_iperm(3,permutation,next_permutation,4)
           permutation=next_permutation
        endif
     enddo
     if (any(permutation.ne.[4,3,2])) error stop 'partial permutation sequence ended incorrectly'
     if (ifindloc([2,4],2,4).ne.2 .or. ifindloc([2,4],2,3).ne.3) &
          error stop 'integer lookup failed'
     write(*,'(a)') 'Math-functions regression: PASS'
  case ('multiply-overflow')
     value=checked_multiply(huge(1),2,'expected overflow')
  case ('multiply8-overflow')
     value8=checked_multiply8(huge(1_8),2_8,'expected overflow')
  case ('add-overflow')
     value=checked_add(huge(1),1,'expected overflow')
  case ('add8-overflow')
     value8=checked_add8(huge(1_8),1_8,'expected overflow')
  case ('power-overflow')
     value=checked_integer_power(2,bit_size(1),'expected overflow')
  case ('factorial-overflow')
     value=factorial(13)
  case ('negative')
     value=checked_multiply(-1,0,'expected negative rejection')
  case ('final-permutation')
     permutation=[4,3,2]
     call get_next_iperm(3,permutation,next_permutation,4)
  case default
     error stop 'unknown math-functions regression mode'
  end select
  if (value.eq.-huge(1) .or. value8.eq.-huge(1_8)) write(*,*) 'unreachable'
end program math_functions_regression
