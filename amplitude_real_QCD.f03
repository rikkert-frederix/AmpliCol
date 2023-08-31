module amplitude_mod
  implicit none
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(1000) :: vertices,order
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
  end type interaction
  type amplitude
     type(current),dimension(1000) :: current_list
     type(interaction),dimension(1000) :: interaction_list
     complex(kind=16),dimension(:),allocatable :: amps
     integer :: n_cur,n_vert
   contains
     procedure :: init_OneOrder
  end type amplitude
contains
  subroutine init_OneOrder(this,n,part,order)
    implicit none
    class(amplitude) :: this
    integer::n
    integer,dimension(n)::part,order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,i
    ! consistency checks:
    ! count number of qqbar pairs (should be event
    ! in colour order, make sure that first and last are quark and anti-quark
    ! also, in the middle, anti-quark should become just before quark.

    this%n_cur=0
    this%n_vert=0
    do isize=1,n-1
       if(isize.eq.1) then
          ! external currents
          do nc=1,n
             this%n_cur=this%n_cur+1
             this%current_list(this%n_cur)%order(1)=order(nc)
             if (nc.le.2 .and. abs(part(nc)).le.6) then
                this%current_list(this%n_cur)%type=anti_current(part(order(nc))) ! switch quark <--> anti-quark for initial states
             else
                this%current_list(this%n_cur)%type=part(order(nc))
             endif
             this%current_list(this%n_cur)%bin=ibset(0,order(nc)-1)
             this%current_list(this%n_cur)%nhel=2
          enddo
       else
          do nc=1,n-isize
             do isplit=1,isize-1
                n1=isplit
                n2=isize-isplit
                bin1=0 ; do i=nc,nc+isplit-1       ; bin1=ibset(bin1,order(i)-1) ; enddo
                bin2=0 ; do i=nc+isplit,nc+isize-1 ; bin2=ibset(bin2,order(i)-1) ; enddo
                do ic1=1,this%n_cur
                   if (this%current_list(ic1)%bin.ne.bin1) cycle
                   do ic2=1,this%n_cur
                      if (this%current_list(ic2)%bin.ne.bin2) cycle
                      call add_if_allowed_threevertex()
                   enddo
                enddo
             enddo
          enddo
       endif
    enddo
  contains
    integer function anti_current(itype)
      implicit none
      integer :: itype
      if (abs(itype).le.6) then
         anti_current=-itype
      else
         anti_current=itype
      endif
    end function anti_current
    subroutine add_if_allowed_threevertex()
      implicit none
      if (this%current_list(ic1)%type.eq.21 .and. this%current_list(ic2)%type.eq.21) then
         ! add the gluon-gluon to gluon vertex
         call add_vertex(0,21)
         ! add the gluon-gluon to tensor vertex
         call add_vertex(1,-21)
      elseif (this%current_list(ic1)%type.eq.-21 .and. this%current_list(ic2)%type.eq.21) then
         ! add a tensor-gluon to gluon vertex
         call add_vertex(2,21)
      elseif (this%current_list(ic1)%type.eq.21 .and. this%current_list(ic2)%type.eq.-21) then
         ! add a gluon-tensor to gluon vertex
         call add_vertex(3,21)
      elseif (this%current_list(ic1)%type.eq.21 .and. &
           (this%current_list(ic2)%type.ge.1 .and. this%current_list(ic2)%type.le.6)) then
         ! add a gluon-quark to quark vertex
         call add_vertex(4,this%current_list(ic2)%type)
      elseif (this%current_list(ic1)%type.eq.21 .and. &
           (this%current_list(ic2)%type.le.-1 .and. this%current_list(ic2)%type.ge.-6)) then
         ! add a gluon-antiquark to antiquark vertex
         call add_vertex(5,this%current_list(ic2)%type)
      elseif ((this%current_list(ic1)%type.ge.1 .and. this%current_list(ic1)%type.le.6) .and. &
           this%current_list(ic2)%type.eq.21) then
         ! add a quark-gluon to quark vertex
         call add_vertex(6,this%current_list(ic1)%type)
      elseif ((this%current_list(ic1)%type.le.-1 .and. this%current_list(ic1)%type.ge.-6) .and. &
           this%current_list(ic2)%type.eq.21) then
         ! add a antiquark-gluon to antiquark vertex
         call add_vertex(7,this%current_list(ic1)%type)
      elseif ((this%current_list(ic1)%type.ge.1 .and. this%current_list(ic1)%type.le.6) .and. &
           (this%current_list(ic2)%type.eq.-this%current_list(ic1)%type)) then
         ! add a quark-antiquark to gluon vertex
         call add_vertex(8,21)
      elseif ((this%current_list(ic1)%type.le.-1 .and. this%current_list(ic1)%type.ge.-6) .and. &
           (this%current_list(ic2)%type.eq.-this%current_list(ic1)%type)) then
         ! add a antiquark-quark to gluon vertex
         call add_vertex(9,21)
      endif
    end subroutine add_if_allowed_threevertex
    subroutine add_vertex(itype,ctype)
      implicit none
      integer :: itype,ctype
      if (isize.eq.n-1 .and. ctype.ne.anti_current(this%current_list(n)%type)) return
      this%n_vert=this%n_vert+1
      this%interaction_list(this%n_vert)%type=itype
      this%interaction_list(this%n_vert)%currents(1)=ic1
      this%interaction_list(this%n_vert)%currents(2)=ic2
      call add_current(ctype)
    end subroutine add_vertex
    subroutine add_current(itype)
      implicit none
      integer :: itype
      integer :: i
      do i=1,this%n_cur
         if (bin1+bin2.eq.this%current_list(i)%bin .and. itype.eq.this%current_list(i)%type) then
            this%current_list(this%n_cur)%n_vert=this%current_list(this%n_cur)%n_vert+1
            this%current_list(this%n_cur)%vertices(this%current_list(this%n_cur)%n_vert)=this%n_vert
            return
         endif
      enddo
      this%n_cur=this%n_cur+1
      this%current_list(this%n_cur)%order(1:isize)=[this%current_list(ic1)%order(1:n1),this%current_list(ic2)%order(1:n2)]
      this%current_list(this%n_cur)%type=itype
      this%current_list(this%n_cur)%bin=bin1+bin2
      this%current_list(this%n_cur)%nhel=this%current_list(ic1)%nhel*this%current_list(ic2)%nhel
      this%current_list(this%n_cur)%vertices(1)=this%n_vert
      this%current_list(this%n_cur)%n_vert=1
    end subroutine add_current
  end subroutine init_OneOrder
end module amplitude_mod
