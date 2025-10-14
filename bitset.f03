module bitset_mod
  implicit none
  private
  integer, parameter :: int_kind = selected_int_kind(18)  ! large enough integer kind

  public :: bitset

  type :: bitset
    integer(kind=int_kind), allocatable :: bits(:)
    integer :: n_bits = 0
  contains
    procedure :: init
    procedure :: set_bit
    procedure :: clear_bit
    procedure :: test_bit
    procedure :: count_bits

    procedure, pass :: bitset_write_unformatted
    procedure, pass :: bitset_read_unformatted

    ! Overloaded bitwise AND operator
    procedure, pass(lhs) :: bitset_and
    generic :: operator(.and.) => bitset_and
  end type bitset

contains
  subroutine init(this, n)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: n
    integer :: n_ints

    this%n_bits = n
    if (allocated(this%bits)) deallocate(this%bits)

    ! Calculate number of integers needed to hold n bits
    n_ints = (n + bit_size(this%bits(1))-1) / bit_size(this%bits(1))
    if (n_ints < 1) n_ints = 1

    allocate(this%bits(n_ints))
    this%bits = 0_int_kind
  end subroutine init


  subroutine set_bit(this, pos)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    if (pos < 1 .or. pos > this%n_bits) return

    idx = (pos - 1) / bit_size(this%bits(1)) + 1
    bitpos = mod(pos - 1, bit_size(this%bits(1)))

    this%bits(idx) = ibset(this%bits(idx), bitpos)
  end subroutine set_bit


  subroutine clear_bit(this, pos)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    if (pos < 1 .or. pos > this%n_bits) return

    idx = (pos - 1) / bit_size(this%bits(1)) + 1
    bitpos = mod(pos - 1, bit_size(this%bits(1)))

    this%bits(idx) = ibclr(this%bits(idx), bitpos)
  end subroutine clear_bit


  logical function test_bit(this, pos)
    class(bitset), intent(in) :: this
    integer, intent(in) :: pos
    integer :: idx, bitpos

    test_bit = .false.
    if (pos < 1 .or. pos > this%n_bits) return

    idx = (pos - 1) / bit_size(this%bits(1)) + 1
    bitpos = mod(pos - 1, bit_size(this%bits(1)))

    test_bit = btest(this%bits(idx), bitpos)
  end function test_bit


  integer function count_bits(this)
    class(bitset), intent(in) :: this
    integer :: i, total

    total = 0
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

    if (lhs%n_bits /= rhs%n_bits) then
      ! bitsizes must match; here we just return empty
      call c%init(0)
      return
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
    ! Write n_bits first so we know how many bits on read
    write(unit, iostat=stat) this%n_bits
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error writing bitset length'
       return
    end if
    
    ! Write the bits array
    write(unit) this%bits

    if (present(iostat)) iostat = 0
  end subroutine bitset_write_unformatted

  
  subroutine bitset_read_unformatted(this, unit, iostat, iomsg)
    class(bitset), intent(inout) :: this
    integer, intent(in) :: unit
    integer, optional, intent(out) :: iostat
    character(len=*), optional, intent(out) :: iomsg

    integer :: stat, n

    ! Read number of bits
    read(unit, iostat=stat) n
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error reading bitset length'
       return
    end if

    call this%init(n)

    ! Read the bits array
    read(unit, iostat=stat) this%bits
    if (stat /= 0) then
       if (present(iostat)) iostat = stat
       if (present(iomsg)) iomsg = 'Error reading bitset bits'
       return
    end if

    if (present(iostat)) iostat = 0
  end subroutine bitset_read_unformatted
end module bitset_mod

