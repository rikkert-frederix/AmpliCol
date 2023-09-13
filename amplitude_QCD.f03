module amplitude_mod
  implicit none
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     complex(kind=8),dimension(:,:),allocatable :: val
     real(kind=8),dimension(0:3) :: pp
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
     complex(kind=8),dimension(:,:),allocatable :: val
  end type interaction
  type amplitude
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     integer :: n_cur,n_vert
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end
   contains
     procedure :: init_OneOrder,evaluate_OneOrder
  end type amplitude
contains
  subroutine init_OneOrder(this,n,part,order)
    implicit none
    class(amplitude) :: this
    integer::n
    integer,dimension(n)::part,order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,iv,i
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
       if (isize.eq.1) then
          ! external currents
          do nc=1,n
             this%n_cur=this%n_cur+1
             this%current_list(this%n_cur)%order(1)=order(nc)
             if (order(nc).le.2 .and. abs(part(order(nc))).le.6) then
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
    call filter_dead_trees()
  contains
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
  subroutine evaluate_OneOrder(this,n,p)
    use FeynmanRules
    implicit none
    class(amplitude) :: this
    integer :: n
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih1,ih2,ih
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
                if (this%current_list(ic)%type.eq.21) then
                   call ext_gluon_cmplx(this%current_list(ic)%pp(0:3),ih-1,1,this%current_list(ic)%val(1:4,ih))
                elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
                   call ext_quark(this%current_list(ic)%pp(0:3),2*ih-3,1,this%current_list(ic)%val(1:4,ih))
                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                   call ext_antiquark(this%current_list(ic)%pp(0:3),2*ih-3,1,this%current_list(ic)%val(1:4,ih))
                endif
                write (*,*) ic,ih,this%current_list(ic)%val(1:4,ih),this%current_list(ic)%order(1)
             enddo
          enddo
          cycle
       endif
       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          do ih1=1,this%current_list(this%interaction_list(iv)%currents(1))%nhel
             do ih2=1,this%current_list(this%interaction_list(iv)%currents(2))%nhel
                ih=(ih1-1)*this%current_list(this%interaction_list(iv)%currents(2))%nhel+ih2
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

    if (.not.allocated(this%amps))allocate(this%amps(1:this%current_list(this%n_cur)%nhel*this%current_list(n)%nhel))
    call compute_amps_from_currents
    
  contains
    subroutine compute_amps_from_currents
      implicit none
      do ih1=1,this%current_list(this%n_cur)%nhel
         do ih2=1,this%current_list(n)%nhel
            ih=(ih1-1)*this%current_list(n)%nhel+ih2
            this%amps(ih)=sum(this%current_list(this%n_cur)%val(1:4,ih1)*this%current_list(n)%val(1:4,ih2))
         enddo
      enddo
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
         this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)=&
              this%current_list(ic)%val(1:dim,1:this%current_list(ic)%nhel)+&
              this%interaction_list(this%current_list(ic)%vertices(iv))%val(1:dim,1:this%current_list(ic)%nhel)
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
  end subroutine evaluate_OneOrder
  
end module amplitude_mod
