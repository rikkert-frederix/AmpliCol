module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: use_mom_dict=.false.
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
     real(kind=8),dimension(0:3) :: pp
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
  end type interaction
  type amplitude_QCD
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: mom_dict
     integer :: n_cur,n_vert,imode,nColOrd,n_qqbar
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end,helmap

     integer,dimension(:),allocatable :: col_value_LC,col_value_NLC,col_value_full
     integer,dimension(:,:),allocatable :: perm,col_index_LC,row_index_LC,col_index_NLC,row_index_NLC, &
          col_index_full,row_index_full
     integer,dimension(:,:,:),allocatable :: row_index,col_index
     integer,dimension(:),allocatable :: n_col_vals
     real(kind=8),dimension(:,:),allocatable :: diff_col_vals
   contains
     procedure :: init,evaluate,init_col,init_col2
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,part,order)
    use math_functions
    implicit none
    class(amplitude_QCD) :: this
    integer::n,imode,j,ic
    integer,dimension(n)::part,order,stand_order,max_order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,iv,i,max_cur,max_vert,nperm
    integer,dimension(:),allocatable :: permutations_dict2
    integer(kind=8),dimension(:),allocatable :: permutations_dict1
    integer(kind=8) :: val,max_val,max_cur_lab
    integer :: jperm,iperm
    logical decompose_4vert
    real(kind=8) :: tAfter,tBefore
    integer(kind=8),dimension(:),allocatable :: current_dict

    decompose_4vert =.true.

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
       call cpu_time(tBefore)
       call set_max_cur()
       call set_max_vert()
       this%nColOrd=1
       allocate(current_dict(max_cur)) 
       call create_current_dict()
       call cpu_time(tAfter)
       write (*,*) '   dictionary created ',tAfter-tBefore
    elseif (this%imode.eq.2) then
       write (*,*) 'WARNING: need to set max_cur and max_vert'
       if (n.lt.7) then
       max_cur=100000
       max_vert=100000
       else
       max_cur=1000000
       max_vert=1000000
       endif
       call cpu_time(tBefore)
       allocate(current_dict(max_cur)) 
       call create_current_dict()
       call cpu_time(tAfter)
       write (*,*) '   dictionary created ',tAfter-tBefore
    endif

    if (use_mom_dict) then
          max_cur_lab=0
          do i=1,n-1
            max_order(i)=n-1-i+1
         enddo
         do i=1,n-1
            max_cur_lab=max_cur_lab+int(max_order(n-1+1-i),kind=8)*int(n-1,kind=8)**int(i-1,kind=8)
         enddo
       allocate(this%mom_dict(max_cur_lab,0:3))
    endif

    do i=1,n
    stand_order(i) = i
    enddo
    stand_order=order

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
             if (order(nc).le.2 .and. abs(part(order(nc))).le.6) then ! initial quark states
                this%current_list(this%n_cur)%type=anti_current(part(order(nc))) ! switch quark <--> anti-quark for initial states
             else
                this%current_list(this%n_cur)%type=part(order(nc))
             endif
             this%current_list(this%n_cur)%bin=ibset(0,order(nc)-1) ! give binary label
             if (this%imode.eq.1) then
                this%current_list(this%n_cur)%nhel=2 ! all possible helicities !!! MASSLESS ONLY
             elseif (this%imode.eq.2) then
                this%current_list(this%n_cur)%nhel=1 ! only one helicity
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
       if (this%n_qqbar.eq.0) then
          if (use_symmetry) then
             allocate(this%perm(1:n-1,1:nperm*2))
          else
             allocate(this%perm(1:n-1,1:nperm*2))
          endif
          do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
             this%perm(1:n-1,nc-this%n_cur_start(n-1)+1)=this%current_list(nc)%order(1:n-1)
          enddo
          if (use_symmetry) then
             do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
                this%perm(1:n-1,nperm+(nc-this%n_cur_start(n-1)+1))=this%current_list(nc)%order(n-1:1:-1)
             enddo
          endif
       elseif (this%n_qqbar.eq.1) then
          allocate(this%perm(1:n-2,1:nperm))
          do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
             this%perm(1:n-2,nc-this%n_cur_start(n-1)+1)=this%current_list(nc)%order(2:n-1)
          enddo
       endif
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
             elseif (i.le.2 .and. part(i).le.-1 .and. part(i).ge.-6) then
                order(1)=i
             elseif (i.gt.2 .and. part(i).le.-1 .and. part(i).ge.-6) then
                order(n)=i
             elseif (i.le.2 .and. part(i).ge.1 .and. part(i).le.6) then
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
               !stop 1
            endif
         else
            if (.not.(part(order(1)).ge.1 .and. part(order(1)).le.6)) then
               write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
               write (*,*) order
               write (*,*) part
               !stop 1
            endif
         endif
         if (order(n).le.2) then
            if (.not.(part(order(n)).ge.1 .and. part(order(n)).le.6)) then
               write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
               write (*,*) order
               write (*,*) part
               !stop 1
            endif
         else
            if (.not.(part(order(n)).le.-1 .and. part(order(n)).ge.-6)) then
               write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
               write (*,*) order
               write (*,*) part
               !stop 1
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
            if (this%n_qqbar.eq.0) max_cur=max_cur+(n-isize)*2
            if (this%n_qqbar.eq.1) max_cur=max_cur+((n-isize-1)*2+1)
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
      ! For imode.eq.1 (summing over helicities), the
      ! this%amps(1:nhel) array contains the helicities (in binary
      ! format) of the external particles. That means for element
      ! 'ihel', the helicities of the external particles are:
      !
      ! do i=1,next
      !    if (btest(ihel-1,i-1)) then
      !       hel(i)=1   <--- positive helicity
      !    else
      !       hel(i)=0   <--- negative helicity
      !    endif
      ! enddo
      !
      ! However, in the way the this%amps() are constructed, the
      ! labeling is according to the colour order, i.e., in the above
      ! loop the i's are over the colour order positions. The helmap
      ! compensates for this, such that this%amps(this%helmap(1:nhel))
      ! contains the order according to the external particle labels.
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
            if (any(this%current_list(ic2)%order(1:n2).eq.order(1))) then
                    return
            endif        
            ! anti-quark should not be part of it, since it will close the current
            if (any(this%current_list(ic1)%order(1:n1).eq.order(n))) return
            if (any(this%current_list(ic2)%order(1:n2).eq.order(n))) return
            if (abs(this%current_list(ic1)%type).ge.1 .and. abs(this%current_list(ic1)%type).le.6) then
              if (this%current_list(ic2)%type.eq.-21) return
            endif
            if (abs(this%current_list(ic2)%type).ge.1 .and. abs(this%current_list(ic2)%type).le.6) then
               if (this%current_list(ic1)%type.eq.-21) return
            endif

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
         if (use_symmetry .and. imode.ne.1) then
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
      if (.not.use_symmetry .or. imode.eq.1) then
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
      integer,dimension(isize) :: ip ! permutation of the current
      integer :: ctype,cur_bin, ic,ik
      integer(kind=8) :: val
      integer :: i,check

!!$      if (ctype.eq.21) then
!!$        ! gluon current
!!$         call get_value(ip,0,val)
!!$      elseif (ctype.eq.-21) then
!!$        ! tensor current
!!$         call get_value(ip,-1,val)
!!$      elseif (ctype.ge.1 .and. ctype.le.6) then
!!$        ! quark current
!!$         call get_value(ip,1,val)
!!$      endif
!!$      call solve_dict(val,ic)

!!$      if ((this%current_list(ic)%n_vert.ne.0).and.(this%current_list(ic)%type.eq.ctype)) then
!!$          if (all(this%current_list(ic)%order(1:isize).eq.ip)) then
!!$          this%current_list(ic)%n_vert=this%current_list(ic)%n_vert+1
!!$          this%current_list(ic)%vertices(this%current_list(ic)%n_vert)=this%n_vert
!!$          this%current_list(ic)%vertex_sign(this%current_list(ic)%n_vert)=vertex_sign
!!$          return
!!$       endif
!!$      endif

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
  
!!$      if (this%imode.eq.1) then  ! temporary hack
              ik = this%n_cur
!!$      elseif (this%imode.eq.2) then
!!$              ik=ic
!!$      endif

      allocate(this%current_list(ik)%order(isize))
      this%current_list(ik)%order(1:isize)=ip(1:isize)
      this%current_list(ik)%type=ctype
      this%current_list(ik)%bin=cur_bin

      this%current_list(ik)%nhel=this%current_list(ic1)%nhel*this%current_list(ic2)%nhel
      if (ctype.eq.21) then
         allocate(this%current_list(ik)%vertices(5*(isize-1)))
         allocate(this%current_list(ik)%vertex_sign(5*(isize-1)))
      elseif (ctype.eq.-21) then
         allocate(this%current_list(ik)%vertices(isize-1))
         allocate(this%current_list(ik)%vertex_sign(isize-1))
      else
         allocate(this%current_list(ik)%vertices(2*(isize-1)))
         allocate(this%current_list(ik)%vertex_sign(2*(isize-1)))
      endif
      this%current_list(ik)%vertices(1)=this%n_vert
      this%current_list(ik)%vertex_sign(1)=vertex_sign
      this%current_list(ik)%n_vert=1
    end subroutine add_current

    subroutine create_current_dict()
      ! create an ordered dictionary that uniquely gives every current
      ! a label. This can be used to quickly find, (O(logN)), a
      ! current in the list of currents
      implicit none
      integer :: size,i,key,j,factor
      integer(kind=8) :: val
      integer,dimension(:),allocatable :: ips_in,ips,ips_re
      integer*2 compar
      integer(kind=8),dimension(max_cur) :: temp

      key=0
      if (imode.eq.2) then
      factor=1
      do isize=1,n-1
         if (isize.eq.1) then
           size=n-isize+1
           factor=n-1
         else
           factor=factor*(n-isize)
           size=factor
         endif
         allocate(ips_in(1:isize))
         allocate(ips_re(1:isize))
         do i=1,isize
            ips_in(i)=i
         enddo
         allocate(ips(1:isize))
         do i=1,size

          if (this%n_qqbar.eq.1) then
           do j=1,isize
            if (ips_in(j).eq.order(n)) then
                    ips_re(j)=n
            elseif (ips_in(j).eq.n) then
                    ips_re(j)=order(n)
            else 
                    ips_re(j)=ips_in(j)
            endif
           enddo
          else
            ips_re=ips_in
          endif
            if (valid_current_order(ips_re))then 
               if (any(ips_re==order(1)) .and. this%n_qqbar.ge.1) then
                 key=key+1
                 call get_value(ips_re,1,val) ! add the quark
                 current_dict(key)=val
               else
                 key=key+1
                 call get_value(ips_re,0,val) ! add the gluon
                 current_dict(key)=val
                 if (isize.ne.1 .and. isize.ne.n-1 .and. decompose_4vert) then
                  key=key+1
                  call get_value(ips_re,-1,val) ! add the tensor
                  current_dict(key)=val
                 endif
               endif
               endif
               if (isize.eq.1) call get_next_iperm(isize,ips_in,ips,n)
               if (isize.gt.1) then
                 call get_next_iperm(isize,ips_in,ips,n-1)
               endif
               ips_in=ips
         enddo
         deallocate(ips_in)
         deallocate(ips_re)
         deallocate(ips)
      enddo
      do i=key+1,max_cur
         current_dict(i)=10000000
      enddo

      !call bubble_test(current_dict,max_cur, temp)
      !current_dict=temp
      !write(*,*) current_dict(1:40)
      !stop 1

      elseif (imode.eq.1) then
         do isize=1,n-1
            if (isize.eq.1) then
               size=n-isize+1 ! also include the external closing current in the dictionary
            else
               size=n-isize
            endif
            allocate(ips_in(1:isize))
            do i=1,size
               do j=1,isize
                  ips_in(j)=order(i+j-1)
                  !ips_in(j)=i+j-1 ! for standard order
               enddo
               !if (.not. valid_current_order(ips_in)) cycle
               if (any(ips_in==order(1)) .and. this%n_qqbar.ge.1) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add the quark
                  current_dict(key)=val
               else
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.1 .and. isize.ne.n-1) then
                     key=key+1
                     call get_value(ips_in,-1,val) ! add the tensor
                     current_dict(key)=val
                  endif
               endif
            enddo
            deallocate(ips_in)
         enddo
      endif
    end subroutine create_current_dict

    subroutine bubble_test(vec,len,ret_vec)
      implicit none
      integer :: len
      integer(kind=8), dimension(len) :: vec,ret_vec
      integer(kind=8) :: temp, bubble, lsup, j

      lsup=len
      do while (lsup .gt. 1)
         bubble = 0 !bubble in the greatest element out of order
          do j = 1, (len-1)
           if (vec(j).gt.vec(j+1)) then
            temp = vec(j)
            vec(j) = vec(j+1)
            vec(j+1) = temp
            bubble = j
           endif
          enddo
         lsup = bubble
       enddo

       ret_vec=vec
      
    end subroutine

    subroutine get_next_iperm(ip,ips_in,ips,n)
     implicit none
     integer :: ip,n,i_up,i,j
     integer,dimension(ip) :: ips,ips_in
     logical :: found
     found=.false.
     ips(1:ip)=ips_in(1:ip)
     do i_up=ip,1,-1
       do while (ips(i_up).lt.n)
          ips(i_up)=ips(i_up)+1
          if (any(ips(1:i_up-1).eq.ips(i_up))) cycle
          found=.true.
          exit
       enddo
       if (found) exit
     enddo
     do i=i_up+1,ip
       do j=1,n
          if (any(ips(1:i).eq.j)) then
             continue
          else
             ips(i)=j
             exit
          endif
       enddo
     enddo
    end subroutine get_next_iperm

    subroutine get_value(ips,itype,val)
      ! Give every current a unique value. This is based on the
      ! (external) particles that are part of the current as well as
      ! the current type.
      implicit none
      integer,dimension(isize) :: ips
      integer :: j,itype
      integer(kind=8) :: val
      val=0
      ! Give a unique identifier based on the external
      ! particles. Simply convert the list to an integer with base
      ! equal to the number of external particles.
      do j=1,isize
         val=val+int(ips(isize+1-j),kind=8)*int(n,kind=8)**int(j-1,kind=8)
      enddo
      ! Take the types into account (we have only 3 types (gluon,
      ! tensor and quark), so multiply by three (and add one for the
      ! tensor and two for quark))
      val=val*int(3,kind=8) ! gluon
      if (itype.eq.-1) then
         val=val+int(1,kind=8) ! tensor
      endif
      if (itype.eq.1) then
         val=val+int(2,kind=8) ! quark
      endif
    end subroutine get_value

    subroutine solve_dict(val,key)
      ! Given the value 'val', find the corresponding key in the
      ! 'current_dict' dictionary. Use a binary search
      ! algorithm. (This only works if the dictionary values are
      ! ordered, and all values only appear once).
      implicit none
      integer :: key,left,middle,right,new_middle
      integer(kind=8) :: val
      left=1
      right=max_cur
     
      if (imode.gt.2) then
      do while (left.le.right)
         middle=(right+left)/2
         if (current_dict(middle).eq.val) then
            key=middle
            return
         elseif(current_dict(middle).gt.val) then
            right=middle-1
         else
            left=middle+1
         endif
      enddo

      elseif (imode.le.2) then
      do i=1,max_cur
         if (current_dict(i) .eq. val) then
            key=i
            return
         endif
      enddo
      endif
      write (*,*) 'value not found in current dictionary',val

      stop 1
    end subroutine solve_dict

  end subroutine init


  subroutine evaluate(this,n,p,hel)
    use FeynmanRules
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,hel
    real(kind=8),dimension(0:3,n) :: p
    real(kind=8),dimension(:,:),allocatable :: mom_dict
    integer :: ic,iv,isize,ih1,ih2,ih,ih_in,ip,i
    integer :: max_cur_lab
    integer :: count_ext,count_vert,count_comb
    real :: tBefore,tAfter
    integer :: ifinal

    if (.not. allocated(this%current_list(1)%val_c) .and. .not.allocated(this%current_list(1)%val_r)) then
       do ic=1,this%n_cur
          if (this%current_list(ic)%type.eq.-21) then
             if (use_real_gluons) then
                allocate(this%current_list(ic)%val_r(1:6,1:this%current_list(ic)%nhel))
             else
                allocate(this%current_list(ic)%val_c(1:6,1:this%current_list(ic)%nhel))
             endif
          elseif (this%current_list(ic)%type.eq.21 .and. use_real_gluons) then
             allocate(this%current_list(ic)%val_r(1:4,1:this%current_list(ic)%nhel))
          else
             allocate(this%current_list(ic)%val_c(1:4,1:this%current_list(ic)%nhel))
          endif
       enddo
       do iv=1,this%n_vert
          if (this%interaction_list(iv)%type.eq.1) then
             if (use_real_gluons) then
                allocate(this%interaction_list(iv)%val_r(1:6,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                             this%current_list(this%interaction_list(iv)%currents(2))%nhel))
             else
                allocate(this%interaction_list(iv)%val_c(1:6,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                             this%current_list(this%interaction_list(iv)%currents(2))%nhel))
             endif
          elseif ((this%interaction_list(iv)%type.eq.0 .or. &
                   this%interaction_list(iv)%type.eq.2 .or. &
                   this%interaction_list(iv)%type.eq.3) .and. use_real_gluons  ) then
             allocate(this%interaction_list(iv)%val_r(1:4,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                          this%current_list(this%interaction_list(iv)%currents(2))%nhel))
          else
             allocate(this%interaction_list(iv)%val_c(1:4,1:this%current_list(this%interaction_list(iv)%currents(1))%nhel* &
                                                          this%current_list(this%interaction_list(iv)%currents(2))%nhel))
          endif
       enddo
       if (this%imode.eq.1) then
          if (use_real_gluons .and. this%n_qqbar.eq.0) then
             allocate(this%amps_r(1:this%current_list(this%n_cur)%nhel*this%current_list(n)%nhel))
          else
             allocate(this%amps(1:this%current_list(this%n_cur)%nhel*this%current_list(n)%nhel))
          endif
       elseif (this%imode.eq.2) then
          if (use_real_gluons .and. this%n_qqbar.eq.0) then
             allocate(this%amps_r(1:this%nColOrd))
          else
             allocate(this%amps(1:this%nColOrd))
          endif
       endif
    endif

    count_ext = 0
    count_vert=0
    count_comb=0

    do isize=1,n-1
       if (isize.eq.1) then
          ! fill the external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)   
             if (this%current_list(ic)%order(1).le.2) then
                this%current_list(ic)%pp(0:3)=-p(0:3,this%current_list(ic)%order(1))
                ifinal=-1
             else
                this%current_list(ic)%pp(0:3)=p(0:3,this%current_list(ic)%order(1))
                ifinal=1
             endif

             do ih=1,this%current_list(ic)%nhel
                count_ext = count_ext + 1
                if (this%current_list(ic)%nhel.eq.1) then
                   if (btest(hel-1,this%current_list(ic)%order(1)-1)) then
                      ih_in=1  ! + helicity 
                   else
                      ih_in=0  ! - helicity 
                   endif
                else
                   ih_in=ih-1
                endif
                if (this%current_list(ic)%type.eq.21) then
                   if (use_real_gluons) then
                           call ext_gluon_real(this%current_list(ic)%pp(0:3),ih_in,ifinal,this%current_list(ic)%val_r(1:4,ih))
                   else
                           call ext_gluon_cmplx(this%current_list(ic)%pp(0:3),ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                      
                   endif
                elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
                        call ext_quark(this%current_list(ic)%pp(0:3),ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                        call ext_antiquark(this%current_list(ic)%pp(0:3),ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                endif
             enddo
          enddo

          if (use_mom_dict) then
             call fill_mom_dict()
          endif
          cycle
       endif

       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          count_vert = count_vert + 1
          do ih1=1,this%current_list(this%interaction_list(iv)%currents(1))%nhel
             do ih2=1,this%current_list(this%interaction_list(iv)%currents(2))%nhel
                ih=(ih2-1)*this%current_list(this%interaction_list(iv)%currents(1))%nhel+ih1

                if (this%interaction_list(iv)%type.eq.0) then
                   if (use_real_gluons) then
                      call threeGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4,ih1),&
                                                this%current_list(this%interaction_list(iv)%currents(1))%pp(0:3),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                this%current_list(this%interaction_list(iv)%currents(2))%pp(0:3),&
                                                this%interaction_list(iv)%val_r(1:4,ih))
                   else
                      call threeGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                           this%current_list(this%interaction_list(iv)%currents(1))%pp(0:3),&
                                      this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%pp(0:3),&
                                           this%interaction_list(iv)%val_c(1:4,ih))
                   endif

                elseif(this%interaction_list(iv)%type.eq.1) then
                   if (use_real_gluons) then
                      call TwoGluonToTensor_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4,ih1),&
                                                 this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                 this%interaction_list(iv)%val_r(1:6,ih))
                   else
                      call TwoGluonToTensor(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                            this%interaction_list(iv)%val_c(1:6,ih))
                   endif
                elseif(this%interaction_list(iv)%type.eq.2) then
                   if (use_real_gluons) then
                      call TensorGluontoGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:6,ih1),&
                                                   this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                   this%interaction_list(iv)%val_r(1:4,ih))
                   else
                      call TensorGluontoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:6,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))
                   endif
                elseif(this%interaction_list(iv)%type.eq.3) then
                   if (use_real_gluons) then
                       call GluonTensortoGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4,ih1),&
                                                    this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:6,ih2),&
                                                    this%interaction_list(iv)%val_r(1:4,ih))
                    else
                       call GluonTensortoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:6,ih2),&
                                               this%interaction_list(iv)%val_c(1:4,ih))
                    endif
                 elseif(this%interaction_list(iv)%type.eq.6) then
                    if (use_real_gluons) then
                       call QuarkGluontoQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                                   this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                   this%interaction_list(iv)%val_c(1:4,ih))
                    else
                       call QuarkGluontoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))
                    endif
                else
                   write (*,*) 'Unknown vertex type: not yet implemented',iv,this%interaction_list(iv)%type
                   stop 1
                endif
             enddo
          enddo
       enddo

       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          count_comb = count_comb + 1
          if (use_mom_dict) then
            this%current_list(ic)%pp(0:3) = this%mom_dict(this%current_list(ic)%bin,0:3)
          else
            call compute_momentum_current()
          endif

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
               if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
                  this%amps_r(this%helmap(ih))=sum(this%current_list(this%n_cur)%val_r(1:4,ih1)*this%current_list(n)%val_r(1:4,ih2))
               else
                  this%amps(this%helmap(ih))=sum(this%current_list(this%n_cur)%val_c(1:4,ih1)*this%current_list(n)%val_c(1:4,ih2))
               endif
            enddo
         enddo

      elseif (this%imode.eq.2) then
         do ic=this%n_cur_start(n-1),this%n_cur_end(n-1)
            if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
               this%amps_r(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_r(1:4,1)*this%current_list(n)%val_r(1:4,1))
            else
               this%amps(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_c(1:4,1)*this%current_list(n)%val_c(1:4,1))
            endif
         enddo

         if (use_symmetry .and. this%n_qqbar.eq.0) then
            do ic=this%n_cur_end(n-1)-this%n_cur_start(n-1)+2, (this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*2
               ip=ic-(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)
               if (use_real_gluons .and. this%n_qqbar.eq.0) then
                  if (mod(n,2).eq.1) then
                     this%amps_r(ic)=-this%amps_r(ip)
                  else
                     this%amps_r(ic)=this%amps_r(ip)
                  endif
               else
                  if (mod(n,2).eq.1) then
                     this%amps(ic)=-this%amps(ip)
                  else
                     this%amps(ic)=this%amps(ip)
                  endif
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
      if (use_real_gluons) then
         this%current_list(ic)%val_r(1:dim,1:this%current_list(ic)%nhel)=0d0
         do iv=1,this%current_list(ic)%n_vert
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_r(1:dim,1:this%current_list(ic)%nhel)=&
                    this%current_list(ic)%val_r(1:dim,1:this%current_list(ic)%nhel)-&
                    this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim,1:this%current_list(ic)%nhel)
            else
               this%current_list(ic)%val_r(1:dim,1:this%current_list(ic)%nhel)=&
                    this%current_list(ic)%val_r(1:dim,1:this%current_list(ic)%nhel)+&
                    this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim,1:this%current_list(ic)%nhel)
            endif
         enddo

      else
         this%current_list(ic)%val_c(1:dim,1:this%current_list(ic)%nhel)=(0d0,0d0)
         do iv=1,this%current_list(ic)%n_vert
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_c(1:dim,1:this%current_list(ic)%nhel)=&
                    this%current_list(ic)%val_c(1:dim,1:this%current_list(ic)%nhel)-&
                    this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim,1:this%current_list(ic)%nhel)
            else
               this%current_list(ic)%val_c(1:dim,1:this%current_list(ic)%nhel)=&
                    this%current_list(ic)%val_c(1:dim,1:this%current_list(ic)%nhel)+&
                    this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim,1:this%current_list(ic)%nhel)
            endif
         enddo
      endif
    end subroutine combine_interactions
    subroutine include_gluon_propagator()
      implicit none
      if (use_real_gluons) then
         call GluonPropagator_real(this%current_list(ic)%val_r,this%current_list(ic)%nhel,this%current_list(ic)%pp)
      else
         call GluonPropagator(this%current_list(ic)%val_c,this%current_list(ic)%nhel,this%current_list(ic)%pp)
      endif
    end subroutine include_gluon_propagator

    subroutine include_quark_propagator()
      implicit none
      call QuarkPropagator(this%current_list(ic)%val_c,this%current_list(ic)%nhel,this%current_list(ic)%pp)
    end subroutine include_quark_propagator

    subroutine fill_mom_dict()
      implicit none
      integer ic
      this%mom_dict=0d0
      do ic=this%n_cur_start(2),this%n_cur_end(n-1)
          if (.not.(all(this%mom_dict(this%current_list(ic)%bin,:).eq.0d0))) cycle
          this%mom_dict(this%current_list(ic)%bin,0:3)=&
          this%current_list(this%interaction_list(this%current_list(ic)%vertices(1))%currents(1))%pp(0:3)+&
          this%current_list(this%interaction_list(this%current_list(ic)%vertices(1))%currents(2))%pp(0:3)
          this%current_list(ic)%pp(0:3)=this%mom_dict(this%current_list(ic)%bin,0:3)
      enddo

   end subroutine fill_mom_dict

     
  end subroutine evaluate

  subroutine init_col2(this,n,part,order,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer :: col_acc,n
    integer,dimension(n) :: part,order,iper,jper
    integer :: iperm,jperm,nperm,ival,iacc
    integer,dimension(1:3) :: n_vals
    integer,parameter :: max_vals=100
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(max_vals,1:3) :: diff_vals
    integer,dimension(:,:),allocatable :: ic,ir
    logical colour_flow
    if (this%n_qqbar.eq.0) then
       nperm=factorial(n-1)
    else
       nperm=factorial(n-2)
    endif
    if (this%n_qqbar.eq.0) then
       colour_flow=.true.
    elseif (this%n_qqbar.eq.1) then
       colour_flow=.false.
    endif
! first check a single row in the colour matrix to determine how many
! different colour factors there are
    n_vals(1:3)=0
    iperm=1
    if (this%n_qqbar.eq.0) then
       iper(1:n)=[this%perm(1:n-1,iperm),n]
    elseif (this%n_qqbar.eq.1) then
       iper(1:n)=[order(1),this%perm(1:n-2,iperm),order(n)]
    endif
    do jperm=iperm,nperm
       if (this%n_qqbar.eq.0) then
          jper(1:n)=[this%perm(1:n-1,jperm),n]
       elseif (this%n_qqbar.eq.1) then
          jper(1:n)=[order(1),this%perm(1:n-2,jperm),order(n)]
       endif
       call compute_color_factor(col_acc,n,iper,jper,col_fac,colour_flow)
       if (iperm.ne.jperm) col_fac(1:3)=col_fac(1:3)*2d0 ! include a factor 2 for the off-diagonal terms
       do iacc=1,3
          if (col_fac(iacc).eq.0d0) cycle
          do ival=1,n_vals(iacc)
             if (col_fac(iacc).eq.diff_vals(ival,iacc)) exit
          enddo
          if (ival.ge.max_vals) then
             write (*,*) 'Too many different colour factors. Increase max_vals',&
                  ival,n_vals(1:3),max_vals
             stop 1
          elseif (ival.eq.n_vals(iacc)+1) then
             ! new colour factor
             n_vals(iacc)=n_vals(iacc)+1
             diff_vals(n_vals(iacc),iacc)=col_fac(iacc)
          endif
       enddo
    enddo
    write (*,*) 'A single row in the colour matrix has',n_vals(1:3),&
         ' different colour factors at LC, NLC and full colour, respectively'
    
! Allocate the arrays now that we now their sizes
    allocate(ic(1:maxval(n_vals(1:3)),1:3))
    allocate(ir(1:maxval(n_vals(1:3)),1:3))
    allocate(this%col_index(1:nperm**2,1:maxval(n_vals(1:3)),1:3))
    allocate(this%row_index(0:nperm,1:maxval(n_vals(1:3)),1:3))
    this%row_index(0,1:maxval(n_vals(1:3)),1:3)=0
    allocate(this%n_col_vals(1:3))
    this%n_col_vals(1:3)=n_vals(1:3)
    allocate(this%diff_col_vals(1:maxval(n_vals(1:3)),1:3))
    do iacc=1,3
       this%diff_col_vals(1:n_vals(iacc),iacc)=diff_vals(1:n_vals(iacc),iacc)
    enddo
    
! Compute all the colour factors and fill the col_index and row_index arrays
    ic=0
    ir=0
    do iperm=1,nperm
       if (this%n_qqbar.eq.0) then
          iper(1:n)=[this%perm(1:n-1,iperm),n]
       elseif (this%n_qqbar.eq.1) then
          iper(1:n)=[order(1),this%perm(1:n-2,iperm),order(n)]
       endif
       do jperm=iperm,nperm ! only include upper triangle (i.e., loop starts at iperm instead of 1)
          if (this%n_qqbar.eq.0) then
             jper(1:n)=[this%perm(1:n-1,jperm),n]
          elseif (this%n_qqbar.eq.1) then
             jper(1:n)=[order(1),this%perm(1:n-2,jperm),order(n)]
          endif
          call compute_color_factor(col_acc,n,iper,jper,col_fac,colour_flow)
          if (iperm.ne.jperm) col_fac(1:3)=col_fac(1:3)*2d0 ! include a factor 2 for the off-diagonal terms
          do iacc=1,3
             if (col_fac(iacc).eq.0d0) cycle
             do ival=1,n_vals(iacc)
                if (col_fac(iacc).eq.diff_vals(ival,iacc)) exit
             enddo
             ic(ival,iacc)=ic(ival,iacc)+1
             ir(ival,iacc)=ir(ival,iacc)+1
             this%col_index(ic(ival,iacc),ival,iacc)=jperm
          enddo
       enddo
       do iacc=1,3
          this%row_index(iperm,1:n_vals(iacc),iacc)=ir(1:n_vals(iacc),iacc)
       enddo
    enddo
  contains
    subroutine compute_color_factor(col_acc,n,iper,jper,col_fac,color_flow)
      use color_algebra
      implicit none
      integer :: i,n,acc,col_acc,color_fac
      real(kind=8),dimension(1:3) :: col_fac
      integer,dimension(n) :: iper,jper
      logical :: color_flow
      real(kind=16) :: col_factor
      if (color_flow) then
         if (this%n_qqbar.ne.0) then
            write (*,*) 'Can only compute color-flow colour factor for all-gluon processes'
            stop 1
         endif
         col_fac=0d0
         call color_flow_factor(n,jper,iper,color_fac)
         if (color_fac.ge.n-2*min(col_acc,0)) then
            col_fac(1)=dble(3**color_fac)
         endif
         if(color_fac.ge.n-2*min(col_acc,1)) then
            col_fac(2)=dble(3**color_fac)
         endif
         if(color_fac.ge.n-2*col_acc) then
            col_fac(3)=dble(3**color_fac)
         endif
      else
         col_fac(1:3)=0d0
         if (col_acc.ge.0) then ! LC
            if (this%n_qqbar.eq.0) then
               if (all(iper.eq.jper)) then
                  col_fac(1)=dble(3**n)
               endif
            elseif (this%n_qqbar.eq.1) then
               if (all(iper.eq.jper)) then
                  col_fac(1)=dble(3**(n-1))
               endif
            endif
         endif
         if (col_acc.ge.1) then ! NLC
            if (this%n_qqbar.eq.0) then
               if (all(iper.eq.jper)) then
                  col_fac(2) = dble(3**n - n * 3**(n-2))
               else
                  call check_NLC(n,jper,iper,acc)
                  col_fac(2)=dble(acc*3**(n-2))
               endif
            elseif (this%n_qqbar.eq.1) then
               if (all(iper.eq.jper)) then
                  col_fac(2) = dble(3**(n-1) - (n-2) * 3**(n-3))
               else
                  call check_NLC_1qqbar(n,jper(2:n-1),iper(2:n-1),acc)
                  col_fac(2)=dble(acc*3**(n-3))
               endif
            endif
         endif
         if (col_acc.ge.2) then
            call Tr_allocate(n)
            if (this%n_qqbar.eq.0) then
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
               col_fac(3)=0d0
!!$               do i=n,n-2*min(col_acc,n),-1 
!!$               do i=-1,n-2*min(col_acc,n),-1 
               do i=n,max(n-2*col_acc,0),-1  ! do not include any Nc
                                             ! contributions with negative
                                             ! powers, since they must cancel.
                  if (i.ge.0) then
                     col_fac(3)=col_fac(3)+dble(coef_nc(i,0)*3**i)
                  else
                     col_fac(3)=col_fac(3)+coef_nc(i,0)*3d0**i
                  endif
               enddo
            elseif (this%n_qqbar.eq.1) then
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=1 ! that term is single string of matrices
               Tr(0,1,1)=2*(n-2)
               Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
               Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
               coef(1)=(1d0,0d0)
               coef_Nc(:,:)=0
               coef_Nc(0,1)=1
               call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
               col_fac(3)=0d0
               do i=n,n-2*min(col_acc,n),-1
                  if (i.ge.0) then
                     col_fac(3)=col_fac(3)+dble(coef_nc(i,0)*3**i)
                  else
                     col_fac(3)=col_fac(3)+coef_nc(i,0)*3d0**i
                  endif
               enddo
            endif
            call Tr_deallocate
         endif
      endif
    end subroutine compute_color_factor
  end subroutine init_col2

  subroutine init_col(this,n,part,order,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer :: col_acc,n,i,jperm,iperm,col_fac,imax,max_val,nperm,nw
    integer,dimension(:),allocatable :: ic,ir,iper,jper
    integer,dimension(n) :: part,order
    real(kind=4) :: tBefore,tAfter
    integer :: n_nlc,add

    !if (any(this%current_list(1:n)%type.ne.21)) then
    !   write (*,*) 'ERROR: colour factor computation assumes all-gluon amplitudes',this%current_list(1:n)%type
    !   stop 1
    !endif

    call cpu_time(tBefore)
    write (*,'(a,i3,a)') ' Setting up colour factor (col_acc =',col_acc,')...'
    allocate(iper(1:n))
    allocate(jper(1:n))
    
    if (this%n_qqbar.eq.0) then
       nperm=factorial(n-1)
       n_nlc = 2 ! 1 for LC, 2 for NLC
       add = 0
    else
       nperm=factorial(n-2)
       n_nlc = 3 ! 1 for LC, 2 for +NLC, 3 for -NLC
       add = 1
    endif

    if (col_acc.ge.2) allocate(this%col_value_full((n+1)/2))
    allocate(this%col_value_LC(1))
    if (col_acc.ge.1) allocate(this%col_value_NLC(n_nlc))
    allocate(ic(max((n+1)/2,n_nlc)))
    allocate(ir(max((n+1)/2,n_nlc)))
    if (col_acc.ge.2) allocate(this%col_index_full(nperm**2,(n+1)/2))
    if (col_acc.ge.2) allocate(this%row_index_full(0:nperm,(n+1)/2))
    if (col_acc.ge.1) allocate(this%col_index_NLC(nperm**2,n_nlc))
    if (col_acc.ge.1) allocate(this%row_index_NLC(0:nperm,n_nlc))
    allocate(this%col_index_LC(nperm,1))
    allocate(this%row_index_LC(0:nperm,1))
    if (col_acc.ge.2) this%row_index_full(0,:)=0
    if (col_acc.ge.1) this%row_index_NLC(0,:)=0
    this%row_index_LC(0,:)=0
    
    if (this%n_qqbar.eq.0) then
       do i=1,(n+1)/2+add
          if (col_acc.gt.2)              this%col_value_full(i)=3**(n-2*(i-1))
          if (col_acc.ge.1 .and. i.le.2) this%col_value_NLC(i)=3**(n-2*(i-1))
          if (i.le.1)                    this%col_value_LC(i)=3**(n-2*(i-1))
       enddo
    elseif (this%n_qqbar.eq.1) then
       do i=1,(n+1)/2+add
          if (col_acc.gt.2)              this%col_value_full(i)=3**(n-1-2*(i-1))
          if (col_acc.ge.1 .and. i.le.3) then
             if (i.eq.1) this%col_value_NLC(i)=3**(n-1-2*(i-1))-(n-2)*3**(n-1-2*((i+1)-1))
             if (i.eq.2) this%col_value_NLC(i)=3**(n-1-2*(i-1))   ! +NLC contributions
             if (i.eq.3) this%col_value_NLC(i)=-3**(n-1-2*((i-1)-1))  ! -NLC contributions
          endif
          if (i.le.1) this%col_value_LC(i)=3**(n-1-2*(i-1))
       enddo
    endif

    !write(*,*) 'LC',this%col_value_NLC(1)
    !write(*,*) 'NLC',this%col_value_NLC(2)
    !write(*,*) 'NLC',this%col_value_NLC(3)

    ic=0
    ir=0
    do iperm=1,nperm
       !write(*,*) 'iperm',iperm
       nw=iperm

       if (this%n_qqbar.eq.0) then
          iper(1:n)=[this%perm(1:n-1,nw),n]
       elseif (this%n_qqbar.eq.1) then
          iper(1:n)=[order(1),this%perm(1:n-2,nw),order(n)]
       endif

       !write(*,*) 'iper',iper
       do jperm=iperm,nperm
          nw=jperm
          if (this%n_qqbar.eq.0) then
             jper(1:n)=[this%perm(1:n-1,nw),n]
            !write(*,*) 'jper',jper
             call compute_color_factor(col_acc,n,iper,jper,col_fac,.true.)
          elseif (this%n_qqbar.eq.1) then
             jper(1:n)=[order(1),this%perm(1:n-2,nw),order(n)]
             !write(*,*) 'jper',jper
             call compute_color_factor(col_acc,n,iper,jper,col_fac,.false.)
          endif

          !write(*,*) 'colfac',col_fac

          if (col_fac.eq.0) cycle

          do i=1,(n+1)/2+add
             if (col_acc.ge.2) then
                if (this%col_value_full(i).eq.col_fac) exit
             elseif (col_acc.ge.1 .and. i.le.n_nlc) then
                if (this%col_value_NLC(i).eq.col_fac) exit
             elseif (col_acc.ge.0 .and. i.le.1) then
                if (this%col_value_LC(i).eq.col_fac) exit
             endif
          enddo

          if (col_acc.eq.1 .and. i.gt.n_nlc) cycle
          if (col_acc.eq.0 .and. i.gt.1) cycle

          ic(i)=ic(i)+1
          ir(i)=ir(i)+1

          if (col_acc.ge.2)                this%col_index_full(ic(i),i)=jperm
          if (col_acc.ge.1.and.i.le.n_nlc) this%col_index_NLC(ic(i),i)=jperm
          if (i.le.1)                      this%col_index_LC(ic(i),i)=jperm
       enddo



       if (col_acc.ge.2) this%row_index_full(iperm,:)=ir(:)
       if (this%n_qqbar.eq.0) then
            if (col_acc.ge.1) this%row_index_NLC(iperm,1:2)=ir(1:2)
       elseif (this%n_qqbar.eq.1) then
            if (col_acc.ge.1) this%row_index_NLC(iperm,1:3)=ir(1:3)
       endif
       this%row_index_LC(iperm,1)=ir(1)
    enddo

    ! remove the one with the most entries.
    if (col_acc.ge.2 .and. this%n_qqbar.eq.0) then
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
          if (this%n_qqbar.eq.0) then
            if (all(iper.eq.jper)) then
               col_fac=3**n
            else
               col_fac=0
            endif
          elseif (this%n_qqbar.eq.1) then
            if (all(iper.eq.jper)) then
               col_fac=3**(n-1)
            else
               col_fac=0
            endif
          endif
         elseif (col_acc.eq.1) then ! NLC
          if (this%n_qqbar.eq.0) then
            if (all(iper.eq.jper)) then
               col_fac = 3**n - n * 3**(n-2)
            else
               call check_NLC(n,jper,iper,acc)
               col_fac=acc*3**(n-2)
            endif

          elseif (this%n_qqbar.eq.1) then
            if (all(iper.eq.jper)) then
               col_fac = 3**(n-1) - (n-2) * 3**(n-3)
            else
               call check_NLC_1qqbar(n,jper(2:n-1),iper(2:n-1),acc)
               col_fac=acc*3**(n-3)
            endif
          endif
         else ! NNLC and beyond
          if (this%n_qqbar.eq.0) then
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
           elseif (this%n_qqbar.eq.1) then
              call Tr_allocate(n)
              Tr(0,0,0)=1 ! one term
              Tr(0,0,1)=1 ! that term is single string of matrices
              Tr(0,1,1)=2*(n-2)
              Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
              Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
              coef(1)=(1d0,0d0)
              coef_Nc(:,:)=0
              coef_Nc(0,1)=1
              call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
              !col_fac=0
              !do i=n,max(n-2*col_acc,0),-1 ! do not include any Nc
               ! contributions with negative
               ! powers, since they must cancel.
              ! col_fac=col_fac+coef_nc(i,0)*3**i
              !enddo
              col_fac=dble(col_factor)
           endif
         endif
      endif
    end subroutine compute_color_factor
  end subroutine init_col

  
end module amplitude_QCD_mod
