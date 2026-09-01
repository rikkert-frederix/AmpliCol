module random_number_state
  implicit none
  logical :: ranmar_initialized=.false.
end module random_number_state

module random_number_interface
  implicit none
  private
  public :: ran2,ntuple,get_offset,get_base,get_moffset,ranmar,rmarin

  ! The production generator and a few standalone utilities provide RAN2 as a
  ! legacy external function.  Keep that implementation replaceable while
  ! giving every modern caller an explicit, shared interface.
  interface
     function ran2() result(value)
       implicit none
       real(kind=8) :: value
     end function ran2
     subroutine ntuple(value,lower,upper,configuration)
       implicit none
       real(kind=8),intent(out) :: value
       real(kind=8),intent(in) :: lower,upper
       integer,intent(in) :: configuration
     end subroutine ntuple
     subroutine get_offset(offset)
       implicit none
       integer,intent(out) :: offset
     end subroutine get_offset
     subroutine get_base(seed)
       implicit none
       integer(kind=8),intent(out) :: seed
     end subroutine get_base
     subroutine get_moffset(offset)
       implicit none
       integer,intent(out) :: offset
     end subroutine get_moffset
     subroutine ranmar(value)
       implicit none
       real(kind=8),intent(out) :: value
     end subroutine ranmar
     subroutine rmarin(first_seed,second_seed)
       implicit none
       integer,intent(in) :: first_seed,second_seed
     end subroutine rmarin
  end interface
end module random_number_interface
