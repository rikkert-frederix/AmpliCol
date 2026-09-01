program bitset_regression
  use bitset_mod, only: bitset
  implicit none
  type(bitset) :: first,second,combined,restored
  character(len=32) :: mode
  character(len=256) :: message
  integer :: iunit,ios,i
  integer(kind=8) :: raw_word,value

  mode='roundtrip'
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  select case (trim(mode))
  case ('roundtrip')
     call first%init(130)
     do i=1,130
        if (any(i.eq.[1,64,65,128,130])) call first%set_bit(i)
     enddo
     if (first%count_bits().ne.5) error stop 'wrong bit count across word boundaries'
     call second%init(130)
     call second%set_bit(64)
     call second%set_bit(130)
     combined=first .and. second
     if (combined%count_bits().ne.2) error stop 'bitset AND returned the wrong mask'
     call first%clear_bit(65)
     if (first%test_bit(65) .or. first%count_bits().ne.4) &
          error stop 'bit clearing failed'

     open(newunit=iunit,status='scratch',access='stream',form='unformatted',&
          action='readwrite')
     call first%bitset_write_unformatted(iunit,ios,message)
     if (ios.ne.0) error stop 'could not serialize bitset'
     rewind(iunit)
     call restored%bitset_read_unformatted(iunit,ios,message)
     if (ios.ne.0) error stop 'could not deserialize bitset'
     do i=1,130
        if (restored%test_bit(i).neqv.first%test_bit(i)) &
             error stop 'bitset round trip changed a bit'
     enddo
     close(iunit)

     call first%init(0)
     if (first%count_bits().ne.0 .or. first%bitset_to_integer().ne.0_8) &
          error stop 'initialized empty bitset is inconsistent'

     ! A failed read must leave the old value untouched.
     call restored%init(5)
     call restored%set_bit(3)
     open(newunit=iunit,status='scratch',access='stream',form='unformatted',&
          action='readwrite')
     write(iunit) 130
     rewind(iunit)
     call restored%bitset_read_unformatted(iunit,ios,message)
     if (ios.eq.0 .or. restored%n_bits.ne.5 .or. .not.restored%test_bit(3)) &
          error stop 'truncated bitset read was not atomic'
     close(iunit)

     ! Non-canonical padding bits are masked on input.
     raw_word=-1_8
     open(newunit=iunit,status='scratch',access='stream',form='unformatted',&
          action='readwrite')
     write(iunit) 1
     write(iunit) raw_word
     rewind(iunit)
     call restored%bitset_read_unformatted(iunit,ios,message)
     if (ios.ne.0 .or. restored%count_bits().ne.1) then
        write(*,*) 'Padding diagnostic:',ios,restored%n_bits,&
             restored%count_bits(),restored%bits
        error stop 'bitset padding was not canonicalized'
     endif
     close(iunit)
     write(*,'(a)') 'Bitset regression: PASS'
  case ('mismatch')
     call first%init(1)
     call second%init(2)
     combined=first .and. second
     error stop 'mismatched bitset AND unexpectedly succeeded'
  case ('out-of-range')
     call first%init(5)
     call first%set_bit(6)
     error stop 'out-of-range bit update unexpectedly succeeded'
  case ('malformed')
     first%n_bits=5
     i=first%count_bits()
     error stop 'malformed bitset unexpectedly succeeded'
  case ('too-large-integer')
     call first%init(64)
     value=first%bitset_to_integer()
     error stop 'oversized integer conversion unexpectedly succeeded'
  case default
     error stop 'unknown bitset regression mode'
  end select
end program bitset_regression
