module amplitude_mod
  implicit none
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     integer,dimension(:,:),allocatable :: val
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
  end type interaction
  type amplitude
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=16),dimension(:),allocatable :: amps
     integer :: n_cur,n_vert
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end
   contains
     procedure :: init_OneOrder
  end type amplitude
contains
  subroutine init_OneOrder(this,n,part,order)
    implicit none
    class(amplitude) :: this
    integer::n
    integer,dimension(n)::part,order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,ivert,i
    ! consistency checks:
    ! - count number of qqbar pairs (should be even)
    ! - in colour order, make sure that first and last are quark and anti-quark
    ! - also, in the middle, anti-quark should become just before quark.

    allocate(this%current_list(1000))
    allocate(this%interaction_list(1000))
    do nc=1,1000
       allocate(this%current_list(nc)%vertices(100))
       allocate(this%current_list(nc)%order(n-1))
    enddo
    allocate(this%n_cur_start(n-1))
    allocate(this%n_cur_end(n-1))
    allocate(this%n_vert_start(2:n-1))
    allocate(this%n_vert_end(2:n-1))
 
    this%n_cur=0
    this%n_vert=0
    do isize=1,n-1
       this%n_cur_start(isize)=this%n_cur+1
       if (isize.ge.2) this%n_vert_start(isize)=this%n_vert+1
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
             this%current_list(this%n_cur)%n_vert=0
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
       this%n_cur_end(isize)=this%n_cur
       if (isize.ge.2) this%n_vert_end(isize)=this%n_vert
    enddo
    ! All done. But there could be currents that are not needed. Filter them out
    call filter_needed()
  contains
    subroutine filter_needed()
      implicit none
      logical,dimension(:),allocatable :: is_needed_cur,is_needed_ver
      integer,dimension(:),allocatable :: where_to_cur,where_to_ver
      integer :: to_skip
      allocate(is_needed_cur(this%n_cur))
      allocate(is_needed_ver(this%n_vert))
      allocate(where_to_cur(this%n_cur))
      allocate(where_to_ver(this%n_vert))
      ! assume nothing is needed
      is_needed_cur(:)=.false.
      is_needed_ver(:)=.false.
      ! loop through the list backward: the closing current is needed, hence
      ! all currents that lead into that closing current are also needed. 
      do nc=this%n_cur,1,-1
         if (is_needed_cur(nc) .or. nc.ge.this%n_cur_start(n-1) .or. nc.le.this%n_cur_end(1)) then
            is_needed_cur(nc)=.true.
            do ivert=1,this%current_list(nc)%n_vert
               is_needed_ver(this%current_list(nc)%vertices(ivert))=.true.
               is_needed_cur(this%interaction_list(this%current_list(nc)%vertices(ivert))%currents(1))=.true.
               is_needed_cur(this%interaction_list(this%current_list(nc)%vertices(ivert))%currents(2))=.true.
            enddo
         endif
      enddo
      ! now we know which ones we can skip. Determine where to move the remaining currents to
      to_skip=0
      do nc=1,this%n_cur
         if (.not. is_needed_cur(nc)) then
            to_skip=to_skip+1
            cycle
         endif
         where_to_cur(nc)=nc-to_skip
      enddo
      to_skip=0
      do ivert=1,this%n_vert
         if (.not. is_needed_ver(ivert)) then
            to_skip=to_skip+1
            cycle
         endif
         where_to_ver(ivert)=ivert-to_skip
      enddo
      ! do the actual shifting of the currents in the list
      do nc=1,this%n_cur
         if (.not.is_needed_cur(nc)) cycle
         this%current_list(where_to_cur(nc))=this%current_list(nc)
         do ivert=1,this%current_list(where_to_cur(nc))%n_vert
            this%current_list(where_to_cur(nc))%vertices(ivert)= &
                 where_to_ver(this%current_list(where_to_cur(nc))%vertices(ivert))
         enddo
      enddo
      ! do the actual shifting of the interactions in the list
      do ivert=1,this%n_vert
         if (.not.is_needed_ver(ivert)) cycle
         this%interaction_list(where_to_ver(ivert))=this%interaction_list(ivert)
         this%interaction_list(where_to_ver(ivert))%currents(1:2)= &
              where_to_cur(this%interaction_list(where_to_ver(ivert))%currents(1:2))
      enddo
      ! and also the shifting of the auxiliary arrays and variables
      do isize=1,n-1
         this%n_cur_start(isize)=where_to_cur(this%n_cur_start(isize))
         this%n_cur_end(isize)=where_to_cur(this%n_cur_end(isize))
         if (isize.ge.2) then
            this%n_vert_start(isize)=where_to_ver(this%n_vert_start(isize))
            this%n_vert_end(isize)=where_to_ver(this%n_vert_end(isize))
         endif
      enddo
      this%n_cur=where_to_cur(this%n_cur)
      this%n_vert=where_to_ver(this%n_vert)
      deallocate(is_needed_ver)
      deallocate(is_needed_cur)
      deallocate(where_to_ver)
      deallocate(where_to_cur)
    end subroutine filter_needed
    
    integer function anti_current(ctype)
      implicit none
      integer :: ctype
      if (abs(ctype).le.6) then
         anti_current=-ctype
      else
         anti_current=ctype
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
           (this%current_list(ic2)%type.eq.anti_current(this%current_list(ic1)%type))) then
         ! add a quark-antiquark to gluon vertex
         call add_vertex(8,21)
      elseif ((this%current_list(ic1)%type.le.-1 .and. this%current_list(ic1)%type.ge.-6) .and. &
           (this%current_list(ic2)%type.eq.anti_current(this%current_list(ic1)%type))) then
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
    subroutine add_current(ctype)
      implicit none
      integer :: ctype
      integer :: i
      ! Check if this interaction can be added to an existing current
      do i=1,this%n_cur
         if (bin1+bin2.eq.this%current_list(i)%bin .and. ctype.eq.this%current_list(i)%type) then
            this%current_list(i)%n_vert=this%current_list(i)%n_vert+1
            this%current_list(i)%vertices(this%current_list(i)%n_vert)=this%n_vert
            return
         endif
      enddo
      ! Need a new current
      this%n_cur=this%n_cur+1
      this%current_list(this%n_cur)%order(1:isize)=[this%current_list(ic1)%order(1:n1),this%current_list(ic2)%order(1:n2)]
      this%current_list(this%n_cur)%type=ctype
      this%current_list(this%n_cur)%bin=bin1+bin2
      this%current_list(this%n_cur)%nhel=this%current_list(ic1)%nhel*this%current_list(ic2)%nhel
      this%current_list(this%n_cur)%vertices(1)=this%n_vert
      this%current_list(this%n_cur)%n_vert=1
    end subroutine add_current
  end subroutine init_OneOrder
!!$  subroutine evaluate_OneOrder(this,n)
!!$    implicit none
!!$    class(amplitude) :: this
!!$    integer :: n
!!$    do isize=1,n-1
!!$       if (isize.eq.1) then
!!$          ! fill the external wave_functions
!!$          do ic=this%n_cur_start(isize),this%n_cur_end(isize)
!!$             do ihel=1,this%current_list(ic)%nhel
!!$                if (this%current_list(ic)%type.eq.21) then
!!$                   call ext_gluon(pp(0,this%current_list(ic)%bin),ihel,1,this%current_list(ic)%val(1:4,ihel))
!!$                elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
!!$                   call ext_quark(pp(0,this%current_list(ic)%bin),ihel,1,this%current_list(ic)%val(1:4,ihel))
!!$                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
!!$                   call ext_antiquark(pp(0,this%current_list(ic)%bin),ihel,1,this%current_list(ic)%val(1:4,ihel))
!!$                endif
!!$             enddo
!!$          enddo
!!$          cycle
!!$       endif
!!$
!!$
!!$       
!!$       do ivert=this%n_vert_start(isize),this%n_vert_end(isize)
!!$          if (this%vert_type(ivert).eq.0) then
!!$             call threeGluon(current(1:4,this%vert_cur(1,ivert)),pp(0:3,this%cur_bin(this%vert_cur(1,ivert))),&
!!$                             current(1:4,this%vert_cur(2,ivert)),pp(0:3,this%cur_bin(this%vert_cur(2,ivert))),&
!!$                             current_out(1:4,ivert))
!!$          elseif(this%vert_type(ivert).eq.1) then
!!$             call TwoGluonToTensor(current(1:4,this%vert_cur(1,ivert)),&
!!$                                   current(1:4,this%vert_cur(2,ivert)),&
!!$                                   current_out(1:6,ivert))
!!$          elseif(this%vert_type(ivert).eq.2) then
!!$             call TensorGluontoGluon(current(1:6,this%vert_cur(1,ivert)),&
!!$                                     current(1:4,this%vert_cur(2,ivert)),&
!!$                                     current_out(1:4,ivert))
!!$          elseif(this%vert_type(ivert).eq.3) then
!!$             call GluonTensortoGluon(current(1:4,this%vert_cur(1,ivert)),&
!!$                                     current(1:6,this%vert_cur(2,ivert)),&
!!$                                     current_out(1:4,ivert))
!!$          elseif(this%vert_type(ivert).eq.99) then
!!$             call FourGluon(current(1:4,this%vert_cur(1,ivert)),&
!!$                            current(1:4,this%vert_cur(2,ivert)),&
!!$                            current(1:4,this%vert_cur(3,ivert)),&
!!$                            current_out(1:4,ivert))
!!$          elseif()
!!$             ...
!!$
!!$
!!$             
!!$          else
!!$             write (*,*) 'Unknown vertex type',ivert,this%vert_type(ivert)
!!$             stop 1
!!$          endif
!!$       enddo
!!$       ! compute the currents by combining the interactions
!!$       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
!!$          if (this%cur_type(ic).eq.0) then ! gluon current
!!$             current(1:4,ic)=0d0
!!$             do ivert=1,this%cur_n_vert(ic)
!!$                if (.not.this%cur_vert_sign(ivert,ic)) then
!!$                   current(1:4,ic)=current(1:4,ic)+current_out(1:4,this%cur_vertices(ivert,ic))
!!$                else
!!$                   current(1:4,ic)=current(1:4,ic)-current_out(1:4,this%cur_vertices(ivert,ic))
!!$                endif
!!$             enddo
!!$             ! include the gluon propagator
!!$             if (isize.ne.this%next-1)  then
!!$                propagator=1d0/(pp(0,this%cur_bin(ic))**2-pp(1,this%cur_bin(ic))**2- &
!!$                                pp(2,this%cur_bin(ic))**2-pp(3,this%cur_bin(ic))**2)
!!$                current(1:4,ic)=current(1:4,ic)*propagator
!!$             endif
!!$          elseif (this%cur_type(ic).eq.-1) then ! tensor current
!!$             current(1:6,ic)=0d0
!!$             do ivert=1,this%cur_n_vert(ic)
!!$                if (.not.this%cur_vert_sign(ivert,ic)) then
!!$                   current(1:6,ic)=current(1:6,ic)+current_out(1:6,this%cur_vertices(ivert,ic))
!!$                else
!!$                   current(1:6,ic)=current(1:6,ic)-current_out(1:6,this%cur_vertices(ivert,ic))
!!$                endif
!!$             enddo
!!$          else
!!$             write (*,*) 'Unknown current type',ic,this%cur_type(ic)
!!$             stop 1
!!$          endif
!!$       enddo
!!$    enddo
!!$  end subroutine evaluate_OneOrder
end module amplitude_mod
