module amplitude_mod
  implicit none
  type current
     integer :: type,bin,nhel
     integer,dimension(:),allocatable :: vertices,vertex_sign,order
  end type current
  type interaction
     integer :: type
     integer,dimension(:),allocatable :: currents
  end type interaction
  type amplitude
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=16),dimension(:),allocatable :: amps
  end type amplitude
contains
  subroutine init_OneOrder(n,part,order)
    implicit none
    integer::n
    integer,dimension(n)::part,order
    ! consistency checks:
    ! count number of qqbar pairs (should be event
    ! in colour order, make sure that first and last are quark and anti-quark
    ! also, in the middle, anti-quark should become just before quark.
    
    
    n_cur=-1
    do isize=1,n-1
       if(isize.eq.1) then
          ! external currents
          do nc=1,next
             n_cur=n_cur+1
             current_list(n_cur)%order(1)=order(nc)
             if (nc.le.2 .and. abs(part(nc)).le.6) then
                current_list(n_cur)%type=-part(nc) ! switch quark <--> anti-quark for initial states
             else
                current_list(n_cur)%type=part(nc)
             endif
             current_list(n_cur)%bin=ibset(0,order(nc)-1)
             current_list(n_cur)%nhel=2
          enddo
       else
          do nc=1,next-isize
             do isplit=1,isize-1
                n1=isplit
                n2=isize-isplit
                bin1=sum(current_list(order(1:nc))%bin)
                bin2=sum(current_list(order(nc+1:nc+isize-1))%bin)
                do ic1=1,n_cur
                   if (current_list(ic1)%bin.ne.bin1) cycle
                   do ic2=1,n_cur
                      if (current_list(ic2)%bin.ne.bin2) cycle
                      call add_if_allowed_threevertex()
                   enddo
                enddo
             enddo
          enddo
       endif
    enddo
  contains
    subroutine add_if_allowed_threevertex()
      implicit none
      if (current_list(ic1)%type.eq.21 .and. current_list(ic2)%type.eq.21) then
         ! add the gluon-gluon to gluon vertex
         call add_vertex(0)
         call add_current(21)
         if (isize.ne.n-1) then
            ! add a gluon-gluon to tensor vertex
            call add_vertex(1)
            call add_current(-21)
         endif
      elseif (current_list(ic1)%type.eq.-21 .and. current_list(ic2)%type.eq.21) then
         ! add a tensor-gluon to gluon vertex
         call add_vertex(2)
         call add_current(21)
      elseif (current_list(ic1)%type.eq.21 .and. current_list(ic2)%type.eq.-21) then
         ! add a gluon-tensor to gluon vertex
         call add_vertex(3)
         call add_current(21)
      elseif (current_list(ic1)%type.eq.21 .and. (current_list(ic2)%type.ge.1 .and. current_list(ic2)%type.le.6)) then
         ! add a gluon-quark to quark vertex
         call add_vertex(4)
         call add_current(current_list(ic2)%type)
      elseif (current_list(ic1)%type.eq.21 .and. (current_list(ic2)%type.le.-1 .and. current_list(ic2)%type.ge.-6)) then
         ! add a gluon-antiquark to antiquark vertex
         call add_vertex(5)
         call add_current(current_list(ic2)%type)
      elseif ((current_list(ic1)%type.ge.1 .and. current_list(ic1)%type.le.6) .and. current_list(ic2)%type.eq.21) then
         ! add a quark-gluon to quark vertex
         call add_vertex(6)
         call add_current(current_list(ic1)%type)
      elseif ((current_list(ic1)%type.le.-1 .and. current_list(ic1)%type.ge.-6) .and. current_list(ic2)%type.eq.21) then
         ! add a antiquark-gluon to antiquark vertex
         call add_vertex(7)
         call add_current(current_list(ic1)%type)
            
      elseif ((current_list(ic1)%type.ge.1 .and. current_list(ic1)%type.le.6) .and. (current_list(ic2)%type.eq.-current_list(ic1)%type)) then
         ! add a quark-antiquark to gluon vertex
         call add_vertex(8)
         call add_current(21)
      elseif ((current_list(ic1)%type.le.-1 .and. current_list(ic1)%type.ge.-6) .and. (current_list(ic2)%type.eq.-current_list(ic1)%type)) then
         ! add a antiquark-quark to gluon vertex
         call add_vertex(9)
         call add_current(21)
            
    end subroutine add_if_allowed_threevertex
  end subroutine init_OneOrder
end module amplitude_mod
