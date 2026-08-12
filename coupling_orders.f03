module coupling_orders
  implicit none
  private

  integer,parameter,public :: coupling_order_mode_automatic=0
  integer,parameter,public :: coupling_order_mode_explicit=1
  integer,parameter,public :: coupling_order_unbounded=-1

  type,public :: coupling_order_selection_type
     integer :: mode=coupling_order_mode_automatic
     integer :: as_min2=coupling_order_unbounded
     integer :: as_max2=coupling_order_unbounded
     integer :: aew_min2=coupling_order_unbounded
     integer :: aew_max2=coupling_order_unbounded
     integer :: resolved_as2=coupling_order_unbounded
   contains
     procedure :: is_valid => coupling_order_selection_is_valid
     procedure :: allows => coupling_order_selection_allows
  end type coupling_order_selection_type

  type(coupling_order_selection_type),save,public :: coupling_order_selection

  public :: reset_coupling_order_selection,set_coupling_order_selection
  public :: resolve_automatic_coupling_order

contains

  subroutine reset_coupling_order_selection()
    coupling_order_selection%mode=coupling_order_mode_automatic
    coupling_order_selection%as_min2=coupling_order_unbounded
    coupling_order_selection%as_max2=coupling_order_unbounded
    coupling_order_selection%aew_min2=coupling_order_unbounded
    coupling_order_selection%aew_max2=coupling_order_unbounded
    coupling_order_selection%resolved_as2=coupling_order_unbounded
  end subroutine reset_coupling_order_selection

  subroutine set_coupling_order_selection(mode,as_min2,as_max2,aew_min2,aew_max2,ok)
    integer,intent(in) :: mode,as_min2,as_max2,aew_min2,aew_max2
    logical,intent(out) :: ok
    type(coupling_order_selection_type) :: candidate

    candidate%mode=mode
    candidate%as_min2=as_min2
    candidate%as_max2=as_max2
    candidate%aew_min2=aew_min2
    candidate%aew_max2=aew_max2
    candidate%resolved_as2=coupling_order_unbounded
    ok=candidate%is_valid()
    if (ok) coupling_order_selection=candidate
  end subroutine set_coupling_order_selection

  subroutine resolve_automatic_coupling_order(as2,ok)
    integer,intent(in) :: as2
    logical,intent(out),optional :: ok
    logical :: valid

    valid=coupling_order_selection%mode.eq.coupling_order_mode_automatic .and. as2.ge.0
    if (valid) coupling_order_selection%resolved_as2=as2
    if (present(ok)) ok=valid
  end subroutine resolve_automatic_coupling_order

  logical function coupling_order_selection_is_valid(this)
    class(coupling_order_selection_type),intent(in) :: this

    coupling_order_selection_is_valid=.false.
    if (this%mode.ne.coupling_order_mode_automatic .and. &
         this%mode.ne.coupling_order_mode_explicit) return
    if (this%as_min2.lt.coupling_order_unbounded .or. &
         this%as_max2.lt.coupling_order_unbounded .or. &
         this%aew_min2.lt.coupling_order_unbounded .or. &
         this%aew_max2.lt.coupling_order_unbounded .or. &
         this%resolved_as2.lt.coupling_order_unbounded) return
    if (this%as_min2.ne.coupling_order_unbounded .and. &
         this%as_max2.ne.coupling_order_unbounded) then
       if (this%as_min2.gt.this%as_max2) return
    endif
    if (this%aew_min2.ne.coupling_order_unbounded .and. &
         this%aew_max2.ne.coupling_order_unbounded) then
       if (this%aew_min2.gt.this%aew_max2) return
    endif
    if (this%mode.eq.coupling_order_mode_automatic) then
       if (this%as_min2.ne.coupling_order_unbounded .or. &
            this%as_max2.ne.coupling_order_unbounded .or. &
            this%aew_min2.ne.coupling_order_unbounded .or. &
            this%aew_max2.ne.coupling_order_unbounded) return
    else
       if (this%resolved_as2.ne.coupling_order_unbounded) return
    endif
    coupling_order_selection_is_valid=.true.
  end function coupling_order_selection_is_valid

  logical function coupling_order_selection_allows(this,as2,aew2)
    class(coupling_order_selection_type),intent(in) :: this
    integer,intent(in) :: as2,aew2

    coupling_order_selection_allows=.false.
    if (as2.lt.0 .or. aew2.lt.0) return
    if (this%mode.eq.coupling_order_mode_automatic) then
       ! Before resolution all pairs are candidates for discovering the global
       ! maximum. Afterwards automatic mode is one exact squared-order slice.
       coupling_order_selection_allows=this%resolved_as2.eq.coupling_order_unbounded .or. &
            as2.eq.this%resolved_as2
       return
    endif
    if (this%as_min2.ne.coupling_order_unbounded .and. as2.lt.this%as_min2) return
    if (this%as_max2.ne.coupling_order_unbounded .and. as2.gt.this%as_max2) return
    if (this%aew_min2.ne.coupling_order_unbounded .and. aew2.lt.this%aew_min2) return
    if (this%aew_max2.ne.coupling_order_unbounded .and. aew2.gt.this%aew_max2) return
    coupling_order_selection_allows=.true.
  end function coupling_order_selection_allows

end module coupling_orders
