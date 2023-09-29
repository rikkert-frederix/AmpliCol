module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:,:),allocatable :: val
     real(kind=8),dimension(0:3) :: pp
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
     complex(kind=8),dimension(:,:),allocatable :: val
  end type interaction
  type amplitude_QCD
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     integer :: n_cur,n_vert,imode,nColOrd,n_qqbar
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end,helmap

     integer,dimension(:),allocatable :: col_value_LC,col_value_NLC,col_value_full
     integer,dimension(:,:),allocatable :: perm,col_index_LC,row_index_LC,col_index_NLC,row_index_NLC, &
          col_index_full,row_index_full
   contains
     procedure :: init,evaluate,init_col
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,part,order)
    implicit none
    class(amplitude_QCD) :: this
    integer::n,imode
    integer,dimension(n)::part,order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,iv,i,max_cur,max_vert,nperm
    if (imode.eq.1) then
       write (*,*) 'Initialising amplitude for:'
       write (*,*) '   - all polarisation/helicity configurations'
       write (*,*) '   - a single colour order'
    elseif (imode.eq.2) then
       write (*,*) 'Initialising amplitude for:'
       write (*,*) '   - a single polarisation/helicity configuration'
       write (*,*) '   - all colour orders'
    else
       write (*,*) 'ERROR unknown operation mode',imode
       stop 1
    endif
    this%imode=imode

    if (this%imode.eq.2) call define_canonical_color_order()

    call check_input_consistency()

    if (this%imode.eq.1) then
       call set_max_cur()
       call set_max_vert()
       this%nColOrd=1
    elseif (this%imode.eq.2) then
       write (*,*) 'WARNING: need to set max_cur and max_vert'
       max_cur=100000
       max_vert=100000
    endif
    allocate(this%current_list(max_cur))
    allocate(this%interaction_list(max_vert))
    allocate(this%n_cur_start(n-1))
    allocate(this%n_cur_end(n-1))
    allocate(this%n_vert_start(2:n-1))
    allocate(this%n_vert_end(2:n-1))

    this%n_cur=0
    this%n_vert=0
    do isize=1,n-1
       this%n_cur_start(isize)=this%n_cur+1
       if (isize.ge.2) this%n_vert_start(isize)=this%n_vert+1
       if (isize.eq.1) then
          ! external currents
          do nc=1,n
             this%n_cur=this%n_cur+1
             allocate(this%current_list(this%n_cur)%order(isize))
             this%current_list(this%n_cur)%order(1)=order(nc)
             if (order(nc).le.2 .and. abs(part(order(nc))).le.6) then
                this%current_list(this%n_cur)%type=anti_current(part(order(nc))) ! switch quark <--> anti-quark for initial states
             else
                this%current_list(this%n_cur)%type=part(order(nc))
             endif
             this%current_list(this%n_cur)%bin=ibset(0,order(nc)-1)
             if (this%imode.eq.1) then
                this%current_list(this%n_cur)%nhel=2
             elseif (this%imode.eq.2) then
                this%current_list(this%n_cur)%nhel=1
             endif
             this%current_list(this%n_cur)%n_vert=0
          enddo
       else
          if (this%imode.eq.1) then
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
          elseif (this%imode.eq.2) then
             do isplit=1,isize-1
                n1=isplit
                n2=isize-isplit
                do ic1=this%n_cur_start(n1),this%n_cur_end(n1)
                   do ic2=this%n_cur_start(n2),this%n_cur_end(n2)
                      call add_if_allowed_threevertex()
                   enddo
                enddo
             enddo
          endif
       endif
       this%n_cur_end(isize)=this%n_cur
       if (isize.ge.2) this%n_vert_end(isize)=this%n_vert
    enddo
    if (this%n_vert.gt.max_vert) then
       write (*,*) 'ERROR: too many interactions: max_vert not set correctly',max_vert,this%n_vert
       stop 1
    endif
    if (this%n_cur.gt.max_cur) then
       write (*,*) 'ERROR: too many currents: max_cur not set correctly',max_cur,this%n_cur
       stop 1
    endif
    
    ! All done. But there could be currents that are not needed. Filter them out
    call filter_dead_trees()
    write (*,*) 'Total number of currents and vertices',this%n_cur,this%n_vert
    ! create the helicity map
    if (this%imode.eq.1) call create_helicity_map()
    ! allocate and fill the colour orders
    if (imode.eq.2) then
       nperm=this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
       allocate(this%perm(1:n-1,1:nperm*2))
       do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
          this%perm(1:n-1,nc-this%n_cur_start(n-1)+1)=this%current_list(nc)%order(1:n-1)
       enddo
       do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
          this%perm(1:n-1,nperm+(nc-this%n_cur_start(n-1)+1))=this%current_list(nc)%order(n-1:1:-1)
       enddo
    endif
  contains
    subroutine define_canonical_color_order()
      use math_functions
      implicit none
      integer :: iord,i
      ! define a canonical colour order
       order(1:n)=0
       if (all(part(1:n).eq.21)) then
          do i=1,n
             order(i)=i
          enddo
          this%nColOrd=factorial(n-1)
       else
          iord=1
          do i=1,n
             if (part(i).eq.21) then
                iord=iord+1
                order(iord)=i
             elseif (i.gt.2 .and. part(i).ge.1 .and. part(i).le.6) then
                order(1)=i
             elseif (i.lt.2 .and. part(i).le.-1 .and. part(i).ge.-6) then
                order(1)=i
             elseif (i.gt.2 .and. part(i).le.-1 .and. part(i).ge.-6) then
                order(n)=i
             elseif (i.gt.2 .and. part(i).ge.1 .and. part(i).le.6) then
                order(n)=i
             endif
          enddo
          this%nColOrd=factorial(n-2)
       endif
       if (any(order(1:n).eq.0)) then
          write (*,*) 'ERROR: canonical order not found'
          stop 1
       endif
     end subroutine define_canonical_color_order
    subroutine check_input_consistency()
      implicit none
      integer,dimension(6) :: quark_flav
      integer :: i,j
      this%n_qqbar=0
      quark_flav=0
      do i=1,n
         if (i.le.2) then
            if (part(i).ne.21) then
               quark_flav(abs(part(i)))=quark_flav(abs(part(i)))-sign(1,part(i))
               if (part(i).lt.0) this%n_qqbar=this%n_qqbar+1
            endif
         else
            if (part(i).ne.21) then
               quark_flav(abs(part(i)))=quark_flav(abs(part(i)))+sign(1,part(i))
               if (part(i).gt.0) this%n_qqbar=this%n_qqbar+1
            endif
         endif
      enddo
      if (any(quark_flav(:).ne.0)) then
         write (*,*) 'ERROR: inconsistent quark flavours',part(1:n)
         stop 1
      endif
      if (this%n_qqbar.gt.1) then
         write (*,*) 'ERROR: code only working for 0, or 1 qqbar pairs',this%n_qqbar
         write (*,*) part
         stop 1
      endif
      if (any(order(:).gt.n) .or. any(order(:).lt.1)) then
         write (*,*) 'ERROR: inconsistent colour order. An element is too large or too small',order
         stop 1
      endif
      do i=1,n-1
         do j=i+1,n
            if (order(i).eq.order(j)) then
               write (*,*) 'ERROR: inconsistent colour order. An element appears twice',order
               stop 1
            endif
         enddo
      enddo
      if (this%n_qqbar.gt.0) then
         if (order(1).le.2) then
            if (.not.(part(order(1)).le.-1 .and. part(order(1)).ge.-6)) then
               write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
               write (*,*) order
               write (*,*) part
               stop 1
            endif
         else
            if (.not.(part(order(1)).ge.1 .and. part(order(1)).le.6)) then
               write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
               write (*,*) order
               write (*,*) part
               stop 1
            endif
         endif
         if (order(n).le.2) then
            if (.not.(part(order(n)).ge.1 .and. part(order(n)).le.6)) then
               write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
               write (*,*) order
               write (*,*) part
               stop 1
            endif
         else
            if (.not.(part(order(n)).le.-1 .and. part(order(n)).ge.-6)) then
               write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
               write (*,*) order
               write (*,*) part
               stop 1
            endif
         endif
      endif
      if (this%n_qqbar.ge.2) then
         do i=1,n
            if (order(i).eq.1 .or. order(i).eq.n) cycle
            if (order(i).eq.2) then
               if (part(order(i)).gt.1 .and. part(order(i)).lt.6) then
                  ! next should be a quark
                  if (.not.(part(order(i+1)).gt.1 .and. part(order(i+1)).lt.6)) then
                     write (*,*) 'ERROR: in the colour order, after an initial state quark should come a final state quak'
                     write (*,*) order
                     write (*,*) part
                     stop 1
                  endif
               endif
            else
               if (part(order(i)).lt.-1 .and. part(order(i)).gt.-6) then
                  ! next should be a quark
                  if (.not.(part(order(i+1)).gt.1 .and. part(order(i+1)).lt.6)) then
                     write (*,*) 'ERROR: in the colour order, after a final state anti-quark should come a quark'
                     write (*,*) order
                     write (*,*) part
                     stop 1
                  endif
               endif
            endif
         enddo
      endif
    end subroutine check_input_consistency
    subroutine set_max_cur()
      ! rough upper bound for the maximum number of currents
      implicit none
      max_cur=0
      do isize=1,n-1
         if (isize.eq.1 .or. isize.eq.n-1) then
            max_cur=max_cur+(n-isize)
         else
            max_cur=max_cur+(n-isize)*2
         endif
      enddo
      max_cur=max_cur+1
    end subroutine set_max_cur
    subroutine set_max_vert()
      ! rough upper bound on the maximum number of interactions
      implicit none
      max_vert=0
      do isize=2,n-1
         if (isize.eq.2) then
            max_vert=max_vert+(n-isize)*3
         else
            max_vert=max_vert+isize*(n-isize)*3
         endif
      enddo
    end subroutine set_max_vert
    subroutine create_helicity_map()
      implicit none
      integer :: nhel,ih
      nhel=product(this%current_list(1:n)%nhel)
      allocate(this%helmap(nhel))
      do ih=1,nhel
         this%helmap(ih)=0
         do i=1,n
            if (btest(ih-1,i-1)) this%helmap(ih)=ibset(this%helmap(ih),order(i)-1)
         enddo
         this%helmap(ih)=this%helmap(ih)+1
      enddo
    end subroutine create_helicity_map
    subroutine filter_dead_trees()
      ! some currents can be removed, since the "tree" starting from some of
      ! the initial state particles might lead to a dead end with no possible
      ! interactions for that current and the remaining external
      ! particles. Hence, they do not need to be computed since they can not
      ! lead to a valid Feynman diagram. To filter them out, one starts at the
      ! end, and goes backwards through the list and see if there are any
      ! currents that were not needed (i.e., they are not the input to a
      ! vertex that is used anywhere).
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
      ! loop through the list backward: if we got to the end, i.e.,
      ! nc.ge.this%n_cur_start(n-1), we now that it is a valid tree. This
      ! means that all inputs to that final current are also needed. By moving
      ! backwards through the list, we can filter out all the branches of the
      ! tree that are needed. (The nc.le.this%n_cur_end(1) is needed only to
      ! make sure that the current that corresponds to the n'th final state
      ! particle is marked as needed, since that current does not enter any of
      ! the trees: it is only needed to close the tree and get the amplitudes.)
      do nc=this%n_cur,1,-1
         if (is_needed_cur(nc) .or. nc.ge.this%n_cur_start(n-1) .or. nc.le.this%n_cur_end(1)) then
            is_needed_cur(nc)=.true.
            do iv=1,this%current_list(nc)%n_vert
               is_needed_ver(this%current_list(nc)%vertices(iv))=.true.
               is_needed_cur(this%interaction_list(this%current_list(nc)%vertices(iv))%currents(1))=.true.
               is_needed_cur(this%interaction_list(this%current_list(nc)%vertices(iv))%currents(2))=.true.
            enddo
         endif
      enddo
      ! now we know which ones we can skip. Determine where to move the remaining currents
      to_skip=0
      do nc=1,this%n_cur
         if (.not. is_needed_cur(nc)) then
            to_skip=to_skip+1
            cycle
         endif
         where_to_cur(nc)=nc-to_skip
      enddo
      to_skip=0
      do iv=1,this%n_vert
         if (.not. is_needed_ver(iv)) then
            to_skip=to_skip+1
            cycle
         endif
         where_to_ver(iv)=iv-to_skip
      enddo
      ! do the actual shifting of the currents in the list
      do nc=1,this%n_cur
         if (.not.is_needed_cur(nc)) cycle
         this%current_list(where_to_cur(nc))=this%current_list(nc)
         do iv=1,this%current_list(where_to_cur(nc))%n_vert
            this%current_list(where_to_cur(nc))%vertices(iv)= &
                 where_to_ver(this%current_list(where_to_cur(nc))%vertices(iv))
         enddo
      enddo
      ! do the actual shifting of the interactions in the list
      do iv=1,this%n_vert
         if (.not.is_needed_ver(iv)) cycle
         this%interaction_list(where_to_ver(iv))=this%interaction_list(iv)
         this%interaction_list(where_to_ver(iv))%currents(1:2)= &
              where_to_cur(this%interaction_list(where_to_ver(iv))%currents(1:2))
      enddo
      ! and also the shifting of the auxiliary arrays and variables
      do isize=1,n-1
         do nc=this%n_cur_start(isize),this%n_cur
            if (where_to_cur(nc).ne.0) then
               this%n_cur_start(isize)=where_to_cur(nc)
               exit
            endif
         enddo
         do nc=this%n_cur_end(isize),1,-1
            if (where_to_cur(nc).ne.0) then
               this%n_cur_end(isize)=where_to_cur(nc)
               exit
            endif
         enddo
         if (isize.ge.2) then
            do iv=this%n_vert_start(isize),this%n_vert
               if (where_to_ver(iv).ne.0) then
                  this%n_vert_start(isize)=where_to_ver(iv)
                  exit
               endif
            enddo
            do iv=this%n_vert_end(isize),1,-1
               if (where_to_ver(iv).ne.0) then
                  this%n_vert_end(isize)=where_to_ver(iv)
                  exit
               endif
            enddo
         endif
      enddo
      do nc=this%n_cur,1,-1
         if (where_to_cur(nc).ne.0) then
            this%n_cur=where_to_cur(nc)
            exit
         endif
      enddo
      do iv=this%n_vert,1,-1
         if (where_to_ver(iv).ne.0) then
            this%n_vert=where_to_ver(iv)
            exit
         endif
      enddo
      deallocate(is_needed_ver)
      deallocate(is_needed_cur)
      deallocate(where_to_ver)
      deallocate(where_to_cur)
    end subroutine filter_dead_trees
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
      ! check if we should consider the current combination, and if
      ! so, and the corresponding vertices to the list.
      implicit none
      if (this%imode.eq.2) then
         if (this%n_qqbar.eq.1) then
            ! if quark is in there, it should be the very first particle
            if (any(this%current_list(ic1)%order(1:n1).eq.order(1)) .and. &
                       this%current_list(ic1)%order(1).ne.order(1)) return
            if (any(this%current_list(ic2)%order(1:n2).eq.order(1))) return
            ! anti-quark should not be part of it, since it will close the current
            if (any(this%current_list(ic1)%order(1:n1).eq.order(n))) return
            if (any(this%current_list(ic2)%order(1:n2).eq.order(n))) return
         elseif (this%n_qqbar.eq.0) then
            ! final gluon should not be part of it, since it will close the current
            if (any(this%current_list(ic1)%order(1:n1).eq.order(n))) return
            if (any(this%current_list(ic2)%order(1:n2).eq.order(n))) return
         else
            write (*,*) 'Only implemented for 0 or 1 qqbar pair',this%n_qqbar
            stop
         endif
         ! check that all particles are different in the two currents:
         if (popcnt(ieor(this%current_list(ic1)%bin,this%current_list(ic2)%bin)).ne.isize) return
         if (use_symmetry) then
            ! For the gluon (and tensor) currents only consider one ordering;
            ! the other will be obtained from symmetry.
            if (this%n_qqbar.eq.0 .or. (this%n_qqbar.eq.1 .and. this%current_list(ic1)%order(1).ne.order(1))) then
               if (maxval(this%current_list(ic1)%order(1:n1)).ge.maxval(this%current_list(ic2)%order(1:n2))) return
            endif
         endif
      endif
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
      if (isize.eq.n-1 .and. ctype.ne.anti_current(this%current_list(n)%type)) return ! dead tree. Filter already here
      this%n_vert=this%n_vert+1
      this%interaction_list(this%n_vert)%type=itype
      this%interaction_list(this%n_vert)%currents(1)=ic1
      this%interaction_list(this%n_vert)%currents(2)=ic2
      call add_all_currents(ctype)
    end subroutine add_vertex
    subroutine add_all_currents(ctype)
      implicit none
      logical :: vertex_sign
      integer,dimension(isize,8) :: ip
      integer :: i,cur_bin,ctype
      if (.not.use_symmetry) then
         cur_bin=this%current_list(ic1)%bin+this%current_list(ic2)%bin
         call add_current(.false.,cur_bin,[this%current_list(ic1)%order(1:n1),this%current_list(ic2)%order(1:n2)],ctype)
         return
      endif
      ! Need to consider the 8 possible permutations (1, 2 or 4 permutations will actually be a valid order)
      ip(1:isize,1)=[this%current_list(ic1)%order(1:n1   ),this%current_list(ic2)%order(1:n2   )]
      ip(1:isize,2)=[this%current_list(ic2)%order(1:n2   ),this%current_list(ic1)%order(1:n1   )]
      ip(1:isize,3)=[this%current_list(ic1)%order(n1:1:-1),this%current_list(ic2)%order(1:n2   )]
      ip(1:isize,4)=[this%current_list(ic2)%order(1:n2   ),this%current_list(ic1)%order(n1:1:-1)]
      ip(1:isize,5)=[this%current_list(ic1)%order(1:n1   ),this%current_list(ic2)%order(n2:1:-1)]
      ip(1:isize,6)=[this%current_list(ic2)%order(n2:1:-1),this%current_list(ic1)%order(1:n1   )]
      ip(1:isize,7)=[this%current_list(ic1)%order(n1:1:-1),this%current_list(ic2)%order(n2:1:-1)]
      ip(1:isize,8)=[this%current_list(ic2)%order(n2:1:-1),this%current_list(ic1)%order(n1:1:-1)]
      do i=1,8
         if (n1.eq.1 .and. (i.eq.3 .or. i.eq.4 .or. i.eq.7 .or. i.eq.8)) cycle
         if (n2.eq.1 .and. (i.eq.5 .or. i.eq.6 .or. i.eq.7 .or. i.eq.8)) cycle
         if (valid_current_order(ip(1:isize,i))) then
            if (i.eq.1 .or. &
                 (i.eq.3 .and. mod(n1,2).eq.1)    .or. (i.eq.4 .and. mod(n1,2).eq.0) .or. &
                 (i.eq.5 .and. mod(n2,2).eq.1)    .or. (i.eq.6 .and. mod(n2,2).eq.0) .or. &
                 (i.eq.7 .and. mod(isize,2).eq.0) .or. (i.eq.8 .and. mod(isize,2).eq.1)) then
               vertex_sign=.false. ! no extra sign needed
            else
               vertex_sign=.true.  ! permutation requires a minus sign
            endif
            cur_bin=this%current_list(ic1)%bin+this%current_list(ic2)%bin
            call add_current(vertex_sign,cur_bin,ip(1:isize,i),ctype)
         endif
      enddo
    end subroutine add_all_currents
    logical function valid_current_order(ip)
      ! Checks that ip(1:isize) is an order for a current to be considered:
      ! the smallest number needs to come before the largest number in this
      ! list.
      implicit none
      integer :: i,maxi,mini,min_loc,max_loc
      integer,dimension(isize) :: ip
      if (this%n_qqbar.eq.1 .and. (any(ip(2:isize).eq.order(1)))) then
         ! if there is a quark, it can only be at the first position
         valid_current_order=.false.
         return
      endif
      if (this%n_qqbar.eq.1 .and. ip(1).eq.order(1)) then
         ! if there is a quark, and it is part of the current (it must be at
         ! position 1), then it is a valid order, since no symmetry can be
         ! used.
         valid_current_order=.true.
         return
      endif
      ! This must be an all-gluon (or tensor) current. Here we take only one
      ! single order. Define it such that smallest label comes before the
      ! biggest. This must be compatible with what orders are skipped in
      ! 'add_if_allowed_threevertex()'.
      maxi=0
      mini=100
      do i=1,isize
         if (ip(i).gt.maxi) then
            maxi=ip(i)
            max_loc=i
         endif
         if (ip(i).lt.mini) then
            mini=ip(i)
            min_loc=i
         endif
      enddo
      if (min_loc.gt.max_loc) then
         valid_current_order=.false.
         return
      endif
      valid_current_order=.true.
    end function valid_current_order
    subroutine add_current(vertex_sign,cur_bin,ip,ctype)
      implicit none
      logical :: vertex_sign
      integer,dimension(isize) :: ip
      integer :: ctype,cur_bin
      integer :: i
      ! Check if this interaction can be added to an existing current
      do i=1,this%n_cur
         if (ctype.ne.this%current_list(i)%type) cycle
         if (cur_bin.ne.this%current_list(i)%bin) cycle
         if (any(this%current_list(i)%order(1:isize).ne.ip(1:isize))) cycle
         this%current_list(i)%n_vert=this%current_list(i)%n_vert+1
         this%current_list(i)%vertices(this%current_list(i)%n_vert)=this%n_vert
         this%current_list(i)%vertex_sign(this%current_list(i)%n_vert)=vertex_sign
         return
      enddo
      ! Need a new current
      this%n_cur=this%n_cur+1
      allocate(this%current_list(this%n_cur)%order(isize))
      this%current_list(this%n_cur)%order(1:isize)=ip(1:isize)
      this%current_list(this%n_cur)%type=ctype
      this%current_list(this%n_cur)%bin=cur_bin
      this%current_list(this%n_cur)%nhel=this%current_list(ic1)%nhel*this%current_list(ic2)%nhel
      if (ctype.eq.21) then
         allocate(this%current_list(this%n_cur)%vertices(5*(isize-1)))
         allocate(this%current_list(this%n_cur)%vertex_sign(5*(isize-1)))
      elseif (ctype.eq.-21) then
         allocate(this%current_list(this%n_cur)%vertices(isize-1))
         allocate(this%current_list(this%n_cur)%vertex_sign(isize-1))
      else
         allocate(this%current_list(this%n_cur)%vertices(2*(isize-1)))
         allocate(this%current_list(this%n_cur)%vertex_sign(2*(isize-1)))
      endif
      this%current_list(this%n_cur)%vertices(1)=this%n_vert
      this%current_list(this%n_cur)%vertex_sign(1)=vertex_sign
      this%current_list(this%n_cur)%n_vert=1
    end subroutine add_current
  end subroutine init
  subroutine evaluate(this,n,p,hel)
    use FeynmanRules
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,hel
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih1,ih2,ih,ih_in,ip
    if (.not. allocated(this%current_list(1)%val)) then
       do ic=1,this%n_cur
          if (this%current_list(ic)%type.eq.-21) then
             allocate(this%current_list(ic)%val(1:6,1:this%current_list(ic)%nhel))
          else
             allocate(this%current_list(ic)%val(1:4,1:this%current_list(ic)%nhel))
          endif
       enddo
       do iv=1,this%n_vert
          if (this%interaction_list(iv)%type.eq.1) then
             allocate(this%interaction_list(iv)%val(1:6,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                          this%current_list(this%interaction_list(iv)%currents(2))%nhel))
          else
             allocate(this%interaction_list(iv)%val(1:4,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                          this%current_list(this%interaction_list(iv)%currents(2))%nhel))
          endif
       enddo
       if (this%imode.eq.1) then
          allocate(this%amps(1:this%current_list(this%n_cur)%nhel*this%current_list(n)%nhel))
       elseif (this%imode.eq.2) then
          allocate(this%amps(1:this%nColOrd))
       endif
    endif


    do isize=1,n-1
       if (isize.eq.1) then
          ! fill the external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)
             if (this%current_list(ic)%order(1).le.2) then
                this%current_list(ic)%pp(0:3)=-p(0:3,this%current_list(ic)%order(1))
             else
                this%current_list(ic)%pp(0:3)=p(0:3,this%current_list(ic)%order(1))
             endif
             do ih=1,this%current_list(ic)%nhel
                if (this%current_list(ic)%nhel.eq.1) then
                   if (btest(hel-1,this%current_list(ic)%order(1)-1)) then
                      ih_in=2
                   else
                      ih_in=1
                   endif
                else
                   ih_in=ih
                endif
                if (this%current_list(ic)%type.eq.21) then
                   call ext_gluon_cmplx(this%current_list(ic)%pp(0:3),ih_in-1,1,this%current_list(ic)%val(1:4,ih))
!!$                   call ext_gluon_real(this%current_list(ic)%pp(0:3),ih_in-1,1,this%current_list(ic)%val(1:4,ih))
                elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
                   call ext_quark(this%current_list(ic)%pp(0:3),ih_in-1,1,this%current_list(ic)%val(1:4,ih))
                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                   call ext_antiquark(this%current_list(ic)%pp(0:3),ih_in-1,1,this%current_list(ic)%val(1:4,ih))
                endif
             enddo
          enddo
          cycle
       endif
       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          do ih1=1,this%current_list(this%interaction_list(iv)%currents(1))%nhel
             do ih2=1,this%current_list(this%interaction_list(iv)%currents(2))%nhel
                ih=(ih2-1)*this%current_list(this%interaction_list(iv)%currents(1))%nhel+ih1
                if (this%interaction_list(iv)%type.eq.0) then
                   call threeGluon(this%current_list(this%interaction_list(iv)%currents(1))%val(1:4,ih1),&
                                        this%current_list(this%interaction_list(iv)%currents(1))%pp(0:3),&
                                   this%current_list(this%interaction_list(iv)%currents(2))%val(1:4,ih2),&
                                        this%current_list(this%interaction_list(iv)%currents(2))%pp(0:3),&
                                        this%interaction_list(iv)%val(1:4,ih))
                elseif(this%interaction_list(iv)%type.eq.1) then
                   call TwoGluonToTensor(this%current_list(this%interaction_list(iv)%currents(1))%val(1:4,ih1),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val(1:4,ih2),&
                                         this%interaction_list(iv)%val(1:6,ih))
                elseif(this%interaction_list(iv)%type.eq.2) then
                   call TensorGluontoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val(1:6,ih1),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val(1:4,ih2),&
                                           this%interaction_list(iv)%val(1:4,ih))
                elseif(this%interaction_list(iv)%type.eq.3) then
                   call GluonTensortoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val(1:4,ih1),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val(1:6,ih2),&
                                           this%interaction_list(iv)%val(1:4,ih))
                elseif(this%interaction_list(iv)%type.eq.6) then
                   call QuarkGluontoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val(1:4,ih1),&
                                          this%current_list(this%interaction_list(iv)%currents(2))%val(1:4,ih2),&
                                          this%interaction_list(iv)%val(1:4,ih))
                else
                   write (*,*) 'Unknown vertex type: not yet implemented',iv,this%interaction_list(iv)%type
                   stop 1
                endif
             enddo
          enddo
       enddo
       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          call compute_momentum_current()
          if (this%current_list(ic)%type.eq.21) then
             call combine_interactions(4)
             ! a gluon current
             if (isize.ne.n-1)  then
                call include_gluon_propagator()
             endif
          elseif ((this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6)) then
             ! a quark current
             call combine_interactions(4)
             if (isize.ne.n-1)  then
                call include_quark_propagator()
             endif
          elseif (this%current_list(ic)%type.eq.-21) then
             ! the non-propagating tensor current
             call combine_interactions(6)
          else
             write (*,*) 'Unknown current type',ic,this%current_list(ic)%type
             stop 1
          endif
       enddo
    enddo

    call compute_amps_from_currents
  contains
    subroutine compute_amps_from_currents
      implicit none
      if (this%imode.eq.1) then
         do ih1=1,this%current_list(this%n_cur)%nhel ! helicities for combined current of particles 1 to n-1 (in the colour order)
            do ih2=1,this%current_list(n)%nhel       ! and for the current for particle n
               ih=(ih2-1)*this%current_list(this%n_cur)%nhel+ih1
               this%amps(this%helmap(ih))=sum(this%current_list(this%n_cur)%val(1:4,ih1)*this%current_list(n)%val(1:4,ih2))
            enddo
         enddo
      elseif (this%imode.eq.2) then
         do ic=this%n_cur_start(n-1),this%n_cur_end(n-1)
            this%amps(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val(1:4,1)*this%current_list(n)%val(1:4,1))
         enddo
         if (use_symmetry .and. this%n_qqbar.eq.0) then
            do ic=this%n_cur_end(n-1)-this%n_cur_start(n-1)+2, (this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*2
               ip=ic-(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)
               if (mod(n,2).eq.1) then
                  this%amps(ic)=-this%amps(ip)
               else
                  this%amps(ic)=this%amps(ip)
               endif
            enddo
         endif
      endif
    end subroutine compute_amps_from_currents
    subroutine compute_momentum_current()
      implicit none
      this%current_list(ic)%pp(0:3)=&
           this%current_list(this%interaction_list(this%current_list(ic)%vertices(1))%currents(1))%pp(0:3)+&
           this%current_list(this%interaction_list(this%current_list(ic)%vertices(1))%currents(2))%pp(0:3)
    end subroutine compute_momentum_current
    subroutine combine_interactions(dim)
      implicit none
      integer :: dim,iv
      this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)=(0d0,0d0)
      do iv=1,this%current_list(ic)%n_vert
         if (this%current_list(ic)%vertex_sign(iv))then
            this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)=&
                 this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)-&
                 this%interaction_list(this%current_list(ic)%vertices(iv))%val(1:dim,1:this%current_list(ic)%nhel)
         else
            this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)=&
                 this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)+&
                 this%interaction_list(this%current_list(ic)%vertices(iv))%val(1:dim,1:this%current_list(ic)%nhel)
         endif
      enddo
    end subroutine combine_interactions
    subroutine include_gluon_propagator()
      implicit none
      call GluonPropagator(this%current_list(ic)%val,this%current_list(ic)%nhel,this%current_list(ic)%pp)
    end subroutine include_gluon_propagator
    subroutine include_quark_propagator()
      implicit none
      call QuarkPropagator(this%current_list(ic)%val,this%current_list(ic)%nhel,this%current_list(ic)%pp)
    end subroutine include_quark_propagator
  end subroutine evaluate


  subroutine init_col(this,n,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer :: col_acc,n,i,jperm,iperm,col_fac,imax,max_val,nperm,nw
    integer,dimension(:),allocatable :: ic,ir,iper,jper
    real(kind=4) :: tBefore,tAfter

    if (any(this%current_list(1:n)%type.ne.21)) then
       write (*,*) 'ERROR: colour factor computation assumes all-gluon amplitudes',this%current_list(1:n)%type
       stop 1
    endif
    
    call cpu_time(tBefore)
    write (*,'(a,i3,a)') ' Setting up colour factor (col_acc =',col_acc,')...'
    allocate(iper(1:n))
    allocate(jper(1:n))
    nperm=factorial(n-1)
    if (col_acc.ge.2) allocate(this%col_value_full((n+1)/2))
    if (col_acc.ge.1) allocate(this%col_value_NLC(2))
    allocate(this%col_value_LC(1))
    allocate(ic(max((n+1)/2,2)))
    allocate(ir(max((n+1)/2,2)))
    if (col_acc.ge.2) allocate(this%col_index_full(nperm**2,(n+1)/2))
    if (col_acc.ge.2) allocate(this%row_index_full(0:nperm,(n+1)/2))
    if (col_acc.ge.1) allocate(this%col_index_NLC(nperm**2,2))
    if (col_acc.ge.1) allocate(this%row_index_NLC(0:nperm,2))
    allocate(this%col_index_LC(nperm,1))
    allocate(this%row_index_LC(0:nperm,1))
    if (col_acc.ge.2) this%row_index_full(0,:)=0
    if (col_acc.ge.1) this%row_index_NLC(0,:)=0
    this%row_index_LC(0,:)=0
    do i=1,(n+1)/2
       if (col_acc.gt.2) this%col_value_full(i)=3**(n-2*(i-1))
       if (col_acc.ge.1 .and. i.le.2) this%col_value_NLC(i)=3**(n-2*(i-1))
       if (i.le.1) this%col_value_LC(i)=3**(n-2*(i-1))
    enddo
    ic=0
    ir=0
    do iperm=1,nperm
       nw=iperm
       iper(1:n)=[this%perm(1:n-1,nw),n]
       do jperm=iperm,nperm
          nw=jperm
          jper(1:n)=[this%perm(1:n-1,nw),n]
          call compute_color_factor(col_acc,n,iper,jper,col_fac,.true.)
          if (col_fac.eq.0) cycle
          do i=1,(n+1)/2
             if (col_acc.ge.2) then
                if (this%col_value_full(i).eq.col_fac) exit
             elseif (col_acc.ge.1 .and. i.le.2) then
                if (this%col_value_NLC(i).eq.col_fac) exit
             elseif (col_acc.ge.0 .and. i.le.1) then
                if (this%col_value_LC(i).eq.col_fac) exit
             endif
          enddo
          if (col_acc.eq.1 .and. i.gt.2) cycle
          if (col_acc.eq.0 .and. i.gt.1) cycle
          ic(i)=ic(i)+1
          ir(i)=ir(i)+1
          if (col_acc.ge.2) this%col_index_full(ic(i),i)=jperm
          if (col_acc.ge.1.and.i.le.2) this%col_index_NLC(ic(i),i)=jperm
          if (i.le.1) this%col_index_LC(ic(i),i)=jperm
       enddo
       if (col_acc.ge.2) this%row_index_full(iperm,:)=ir(:)
       if (col_acc.ge.1) this%row_index_NLC(iperm,1:2)=ir(1:2)
       this%row_index_LC(iperm,1)=ir(1)
    enddo
    
    ! remove the one with the most entries.
    if (col_acc.ge.2) then
       imax=0
       max_val=0
       do i=1,(n+1)/2
          if (this%row_index_full(nperm,i).gt.max_val) then
             max_val=max(this%row_index_full(nperm,i),max_val)
             imax=i
          endif
       enddo
       if (all(this%row_index_full(nperm,:).ne.0)) then
          this%row_index_full(:,imax)=0
          this%col_value_full(:)=this%col_value_full(:)-this%col_value_full(imax)
       endif
    endif
    call cpu_time(tAfter)
    write (*,*) '... colour setup in',tAfter-tBefore,'seconds'
  contains
    subroutine compute_color_factor(col_acc,n,iper,jper,col_fac,color_flow)
      use color_algebra
      implicit none
      integer :: i,col_fac,n,acc,col_acc
      integer,dimension(n) :: iper,jper
      logical :: color_flow
      real(kind=16) :: col_factor
      if (color_flow) then
         call color_flow_factor(n,jper,iper,col_fac)
         if (col_fac.ge.n-2*col_acc) then
            col_fac=3**col_fac
         else
            col_fac=0
         endif
      else
         if (col_acc.eq.0) then ! LC
            if (all(iper.eq.jper)) then
               col_fac=3**n
            else
               col_fac=0
            endif
         elseif (col_acc.eq.1) then ! NLC
            if (all(iper.eq.jper)) then
               col_fac = 3**n - n * 3**(n-2)
            else
               call check_NLC(n,jper,iper,acc)
               col_fac=acc*3**(n-2)
            endif
         else ! NNLC and beyond
            Tr(0,0,0)=1 ! one term
            Tr(0,0,1)=2 ! that term is a product of two terms
            Tr(0,1,1)=n ! both terms in the product are a trace with next matrices
            Tr(0,2,1)=n
            Tr(1:n,1,1)=iper(1:n) ! the order of the matrices in each term
            Tr(1:n,2,1)=jper(1:n)
            call Tr_complex_conjugate(2,1) ! take the complex conjugate of the jperm term
            coef(1)=(1d0,0d0)
            coef_Nc(:,:)=0
            coef_Nc(0,1)=1
            ! compute the colour factor by simplifying the colour string
            call Tr_full_simplify(col_factor) 
            col_fac=0
            do i=n,max(n-2*col_acc,0),-1 ! do not include any Nc
               ! contributions with negative
               ! powers, since they must cancel.
               col_fac=col_fac+coef_nc(i,0)*3**i
            enddo
         endif
      endif
    end subroutine compute_color_factor
  end subroutine init_col

  
end module amplitude_QCD_mod
