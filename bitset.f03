module bitset_mod
  implicit none
  private
  integer, parameter :: int_kind = selected_int_kind(18)  ! large enough integer kind
  integer, parameter :: max_bitset_bits = 10000000
  integer, parameter :: bits_per_word = int(bit_size(0_int_kind),kind=kind(0))

  public :: bitset,max_bitset_bits

  type :: bitset
    integer(kind=int_kind), allocatable :: bits(:)
    integer :: n_bits = 0
  contains
    procedure :: init
    procedure :: set_bit
    procedure :: clear_bit
    procedure :: test_bit
    procedure :: count_bits
    procedure :: bitset_to_integer

    procedure, pass :: bitset_write_unformatted
    procedure, pass :: bitset_read_unformatted

    ! Overloaded bitwise AND operator
    procedure, pass(lhs) :: bitset_and
    generic :: operator(.and.) => bitset_and
  end type bitset

contains
  integer function words_for_bits(n)
    integer, intent(in) :: n
    integer :: word_bits

    word_bits = bits_per_word
    words_for_bits = 1
    if (n > 0) words_for_bits = 1 + (n-1)/word_bits
  end function words_for_bits

  subroutine init(this, n)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: n
    integer :: n_ints, stat
    character(len=256) :: message

    if (n < 0 .or. n > max_bitset_bits) then
      write (*,*) 'Error: invalid or unsupported bitset length:',n
      stop 1
    end if

    this%n_bits = n
    if (allocated(this%bits)) deallocate(this%bits)

    ! Calculate number of integers needed to hold n bits
    n_ints = words_for_bits(n)

    allocate(this%bits(n_ints),stat=stat,errmsg=message)
    if (stat /= 0) then
      this%n_bits = 0
      write (*,*) 'Error allocating bitset storage:',trim(message)
      stop 1
    end if
    this%bits = 0_int_kind
  end subroutine init


  subroutine set_bit(this, pos)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    if (this%n_bits < 0 .or. this%n_bits > max_bitset_bits .or. &
         .not.allocated(this%bits)) then
      write (*,*) 'Error: malformed bitset in set_bit'
      stop 1
    end if
    if (size(this%bits) /= words_for_bits(this%n_bits)) then
      write (*,*) 'Error: malformed bitset in set_bit'
      stop 1
    end if
    if (pos < 1 .or. pos > this%n_bits) then
      write (*,*) 'Error: bit position is outside the bitset in set_bit:',&
           pos,this%n_bits
      stop 1
    end if

    idx = (pos - 1) / bits_per_word + 1
    bitpos = mod(pos - 1, bits_per_word)

    this%bits(idx) = ibset(this%bits(idx), bitpos)
  end subroutine set_bit


  subroutine clear_bit(this, pos)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    if (this%n_bits < 0 .or. this%n_bits > max_bitset_bits .or. &
         .not.allocated(this%bits)) then
      write (*,*) 'Error: malformed bitset in clear_bit'
      stop 1
    end if
    if (size(this%bits) /= words_for_bits(this%n_bits)) then
      write (*,*) 'Error: malformed bitset in clear_bit'
      stop 1
    end if
    if (pos < 1 .or. pos > this%n_bits) then
      write (*,*) 'Error: bit position is outside the bitset in clear_bit:',&
           pos,this%n_bits
      stop 1
    end if

    idx = (pos - 1) / bits_per_word + 1
    bitpos = mod(pos - 1, bits_per_word)

    this%bits(idx) = ibclr(this%bits(idx), bitpos)
  end subroutine clear_bit


  logical function test_bit(this, pos)
    class(bitset), intent(in) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    test_bit = .false.
    if (this%n_bits == 0 .and. .not.allocated(this%bits)) return
    if (this%n_bits < 0 .or. this%n_bits > max_bitset_bits .or. &
         .not.allocated(this%bits)) then
      write (*,*) 'Error: malformed bitset in test_bit'
      stop 1
    end if
    if (size(this%bits) /= words_for_bits(this%n_bits)) then
      write (*,*) 'Error: malformed bitset in test_bit'
      stop 1
    end if
    if (pos < 1 .or. pos > this%n_bits) then
      write (*,*) 'Error: bit position is outside the bitset in test_bit:',&
           pos,this%n_bits
      stop 1
    end if

    idx = (pos - 1) / bits_per_word + 1
    bitpos = mod(pos - 1, bits_per_word)

    test_bit = btest(this%bits(idx), bitpos)
  end function test_bit


  integer function count_bits(this)
    class(bitset), intent(in) :: this
    integer :: i, total

    total = 0
    if (this%n_bits == 0 .and. .not.allocated(this%bits)) then
      count_bits = 0
      return
    end if
    if (this%n_bits < 0 .or. this%n_bits > max_bitset_bits .or. &
         .not.allocated(this%bits)) then
      write (*,*) 'Error: malformed bitset in count_bits'
      stop 1
    end if
    if (size(this%bits) /= words_for_bits(this%n_bits)) then
      write (*,*) 'Error: malformed bitset in count_bits'
      stop 1
    end if
    do i = 1, size(this%bits)
      total = total + popcnt(this%bits(i))
    end do
    count_bits = total
  end function count_bits


  ! Overloaded bitwise AND operator: c = a .and. b
  function bitset_and(lhs, rhs) result(c)
    class(bitset), intent(in) :: lhs
    class(bitset), intent(in) :: rhs
    type(bitset) :: c
    integer :: i

    if (lhs%n_bits < 0 .or. lhs%n_bits > max_bitset_bits .or. &
         rhs%n_bits < 0 .or. rhs%n_bits > max_bitset_bits) then
      write (*,*) 'Error: malformed operand in bitset AND'
      stop 1
    end if
    if ((lhs%n_bits.ne.0 .and. .not.allocated(lhs%bits)) .or. &
         (rhs%n_bits.ne.0 .and. .not.allocated(rhs%bits))) then
      write (*,*) 'Error: malformed operand in bitset AND'
      stop 1
    endif
    if (lhs%n_bits /= rhs%n_bits) then
      write (*,*) 'Error: bitset AND operands have different lengths:',&
           lhs%n_bits,rhs%n_bits
      stop 1
    end if
    if (lhs%n_bits.eq.0 .and. .not.allocated(lhs%bits) .and. &
         .not.allocated(rhs%bits)) then
      call c%init(0)
      return
    endif
    if (.not.allocated(lhs%bits) .or. .not.allocated(rhs%bits)) then
      write (*,*) 'Error: malformed operand in bitset AND'
      stop 1
    end if
    if (size(lhs%bits) /= words_for_bits(lhs%n_bits) .or. &
         size(rhs%bits) /= words_for_bits(rhs%n_bits)) then
      write (*,*) 'Error: malformed operand in bitset AND'
      stop 1
    end if

    call c%init(lhs%n_bits)
    do i = 1, size(lhs%bits)
      c%bits(i) = iand(lhs%bits(i), rhs%bits(i))
    end do
  end function bitset_and



  subroutine bitset_write_unformatted(this, unit, iostat, iomsg)
    class(bitset), intent(in) :: this
    integer, intent(in) :: unit
    integer, optional, intent(out) :: iostat
    character(len=*), optional, intent(out) :: iomsg
    integer :: stat
    character(len=256) :: message
    if (this%n_bits < 0 .or. this%n_bits > max_bitset_bits .or. &
         .not.allocated(this%bits)) then
       if (present(iostat)) iostat = 1
       if (present(iomsg)) iomsg = 'Cannot write malformed bitset'
       return
    end if
    if (size(this%bits) /= words_for_bits(this%n_bits)) then
       if (present(iostat)) iostat = 1
       if (present(iomsg)) iomsg = 'Cannot write malformed bitset'
       return
    end if
    ! Write n_bits first so we know how many bits on read
    write(unit, iostat=stat) this%n_bits
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error writing bitset length'
       return
    end if
    
    ! Write the bits array
    write(unit, iostat=stat, iomsg=message) this%bits
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error writing bitset contents: '//trim(message)
       return
    end if

    if (present(iostat)) iostat = 0
    if (present(iomsg)) iomsg = ''
  end subroutine bitset_write_unformatted

  
  subroutine bitset_read_unformatted(this, unit, iostat, iomsg)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: unit
    integer, optional, intent(out) :: iostat
    character(len=*), optional, intent(out) :: iomsg

    integer :: stat, n, remainder, allocation_status
    integer(kind=int_kind),allocatable :: new_bits(:)
    character(len=256) :: message

    ! Read number of bits
    read(unit, iostat=stat) n
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error reading bitset length'
       return
    end if

    if (n < 0 .or. n > max_bitset_bits) then
       if (present(iostat)) iostat = 1
       if (present(iomsg)) iomsg = 'Invalid or unsupported bitset length'
       return
    end if

    allocate(new_bits(words_for_bits(n)),stat=allocation_status,errmsg=message)
    if (allocation_status /= 0) then
       if (present(iostat)) iostat = allocation_status
       if (present(iomsg)) iomsg = 'Error allocating bitset input: '//trim(message)
       return
    endif

    ! Read into temporary storage so a truncated or malformed record leaves
    ! the previous bitset unchanged.
    read(unit, iostat=stat, iomsg=message) new_bits
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error reading bitset bits: '//trim(message)
       return
    end if

    remainder=mod(n,bits_per_word)
    if (remainder /= 0) new_bits(size(new_bits))=&
         ibits(new_bits(size(new_bits)),0,remainder)
    if (allocated(this%bits)) deallocate(this%bits)
    call move_alloc(new_bits,this%bits)
    this%n_bits=n

    if (present(iostat)) iostat = 0
    if (present(iomsg)) iomsg = ''
  end subroutine bitset_read_unformatted

  function bitset_to_integer(this) result(val)
    class(bitset), intent(in) :: this
    integer(kind=8) :: val
    integer :: total_bits, i

    ! Check if the bitset fits in a single integer
    total_bits = this%n_bits
    if (total_bits < 0 .or. total_bits > max_bitset_bits) then
       print *, 'Error: malformed bitset length:',total_bits
       stop 1
    endif
    if (total_bits.eq.0 .and. .not.allocated(this%bits)) then
       val=0_8
       return
    endif
    if (.not.allocated(this%bits)) then
       print *, 'Error: malformed bitset storage'
       stop 1
    endif
    if (size(this%bits).ne.words_for_bits(total_bits)) then
       print *, 'Error: malformed bitset storage'
       stop 1
    endif
    if (total_bits > int(bit_size(val),kind=kind(total_bits))-1) then
       print *, "Error: bitset too large to fit in a positive signed integer of kind=", int_kind
       stop 1
    end if

    val = 0_8
    do i = 1, total_bits
       if (test_bit(this, i)) then
          val = ibset(val, i - 1)
       end if
    end do
  end function bitset_to_integer
  
end module bitset_mod
