module ext_wfs
  implicit none
  integer :: next
  real(kind=8),dimension(:,:),allocatable :: p_ext,wf_ext
end module ext_wfs

module amplitude_mod
  implicit none
  type amplitude
     integer :: next,nperm,isize,col_acc
     integer,dimension(:),allocatable :: ncomb,nsplit2,nsplit3,istart,inverted,nhel
     integer,dimension(:,:),allocatable :: ips,colmap,helmap
     integer,dimension(:,:,:),allocatable :: imap2,imap3
     real(kind=8),dimension(:,:),allocatable :: amps
     logical :: sum_hel
     integer,dimension(:,:),allocatable :: col_index,row_index
     integer,dimension(:),allocatable :: col_value
   contains
     procedure :: init,evaluate,init_onlycol,evaluate_order,init_CSR,init_onlycol_CSR
  end type amplitude
  type amplitude_cache
     integer :: next,nperm,col_acc
     integer,dimension(:),allocatable :: col_value_LC,col_value_NLC,col_value_full
     integer,dimension(:,:),allocatable :: perm,col_index_LC,row_index_LC,col_index_NLC,row_index_NLC, &
          col_index_full,row_index_full
     real(kind=8),dimension(:),allocatable :: amps
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,cur_type,cur_n_vert,vert_type,n_vert_start,n_vert_end,cur_bin
     integer,dimension(:,:),allocatable :: vert_cur,cur_vertices
     logical,dimension(:,:),allocatable :: cur_vert_sign
   contains
     procedure :: setup_imap_cache,setup_colmap_cache,setup_colmap_cache_NLC,evaluate_cache
  end type amplitude_cache
  type col_amp
     real(kind=8),dimension(:,:,:),allocatable :: wf,pp
     integer,dimension(:),allocatable :: order
     integer :: n
   contains
     procedure :: evaluate_order_v3
  end type col_amp
  type(col_amp),dimension(:),allocatable,public :: col_amp_list
  private
  public :: amplitude,amplitude_cache,col_amp
  public :: iperm_encode,iperm_decode,factorial,factorial_dble,factorial8,evaluate_order_v2,evaluate_order_v3
contains
  subroutine init(this,next,col_acc,sum_hel,order)
    implicit none
    class(amplitude) :: this
    integer :: next,col_acc,i,j
    integer,dimension(next),optional :: order
    integer,dimension(next) :: o
    logical :: sum_hel
    this%next=next
    this%sum_hel=sum_hel
    if (present(order)) then
       this%nperm=1
       ! Bring the colour order to a canonical order (final in the
       ! list should be particle 'next' such that the momenta get
       ! assigned correctly).
       o=order
       do i=1,next
          if (o(i).eq.next) then
             do j=0,next-1
                order(j+1)=o(1+mod(i+j,next))
             enddo
             exit
          endif
       enddo
    else
       this%nperm=factorial(next-1)
    endif
    this%col_acc=col_acc
    write (*,*) 'Setup colmap'
    call setup_colmap(this,next,order)
    write (*,*) 'Setup imap'
    call setup_imap(this,order)
    write (*,*) 'Setup helmap' 
    call setup_helmap(this,next)
    if (.not. allocated(this%amps)) allocate(this%amps(0:this%nhel(this%isize+1)-1,this%nperm))
  end subroutine init

  subroutine init_CSR(this,next,col_acc,sum_hel,order)
    implicit none
    class(amplitude) :: this
    integer :: next,col_acc,i,j
    integer,dimension(next),optional :: order
    integer,dimension(next) :: o
    logical :: sum_hel
    this%next=next
    this%sum_hel=sum_hel
    if (present(order)) then
       this%nperm=1
       ! Bring the colour order to a canonical order (final in the
       ! list should be particle 'next' such that the momenta get
       ! assigned correctly).
       o=order
       do i=1,next
          if (o(i).eq.next) then
             do j=0,next-1
                order(j+1)=o(1+mod(i+j,next))
             enddo
             exit
          endif
       enddo
    else
       this%nperm=factorial(next-1)
    endif
    this%col_acc=col_acc
    write (*,*) 'Setup colmap'
    call setup_colmap_CSR(this,next,order)
    write (*,*) 'Setup imap'
    call setup_imap(this,order)
    write (*,*) 'Setup helmap' 
    call setup_helmap(this,next)
    if (.not. allocated(this%amps)) allocate(this%amps(0:this%nhel(this%isize+1)-1,this%nperm))
  end subroutine init_CSR

  subroutine init_onlycol_CSR(this,next,col_acc,sum_hel,order)
    implicit none
    class(amplitude) :: this
    integer :: next,col_acc,i,j
    integer,dimension(next),optional :: order
    integer,dimension(next) :: o
    logical :: sum_hel
    this%next=next
    this%sum_hel=sum_hel
    if (present(order)) then
       this%nperm=1
       ! Bring the colour order to a canonical order (final in the
       ! list should be particle 'next' such that the momenta get
       ! assigned correctly).
       o=order
       do i=1,next
          if (o(i).eq.next) then
             do j=0,next-1
                order(j+1)=o(1+mod(i+j,next))
             enddo
             exit
          endif
       enddo
    else
       this%nperm=factorial(next-1)
    endif
    this%col_acc=col_acc
    write (*,*) 'Setup colmap' 
    call setup_colmap_CSR(this,next,order)
  end subroutine init_onlycol_CSR

  subroutine init_onlycol(this,next,col_acc,sum_hel,order)
    implicit none
    class(amplitude) :: this
    integer :: next,col_acc,i,j
    integer,dimension(next),optional :: order
    integer,dimension(next) :: o
    logical :: sum_hel
    this%next=next
    this%sum_hel=sum_hel
    if (present(order)) then
       this%nperm=1
       ! Bring the colour order to a canonical order (final in the
       ! list should be particle 'next' such that the momenta get
       ! assigned correctly).
       o=order
       do i=1,next
          if (o(i).eq.next) then
             do j=0,next-1
                order(j+1)=o(1+mod(i+j,next))
             enddo
             exit
          endif
       enddo
    else
       this%nperm=factorial(next-1)
    endif
    this%col_acc=col_acc
    write (*,*) 'Setup colmap' 
    call setup_colmap(this,next,order)
  end subroutine init_onlycol

  subroutine setup_colmap_cache(this,col_acc)
    use color_algebra
    implicit none
    class(amplitude_cache) :: this
    integer :: col_acc,n,i,jperm,iperm,col_fac,imax,max_val,nperm,nw
    integer,dimension(:),allocatable :: ic,ir,iper,jper
    real(kind=4) :: tBefore,tAfter

    call cpu_time(tBefore)
    write (*,'(a,i3,a)') ' Setting up colour factor (col_acc =',col_acc,')...'
    n=this%next
    allocate(iper(1:n))
    allocate(jper(1:n))
    nperm=factorial(n-1)
    this%col_acc=col_acc
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
       if (iperm.le.nperm/2) then
          nw=iperm
          iper(1:n)=[this%perm(1:n-1,nw),n]
       else
          nw=iperm-nperm/2
          iper(1:n)=[this%perm(n-1:1:-1,nw),n]
       endif
       do jperm=iperm,nperm
          if (jperm.le.nperm/2) then
             nw=jperm
             jper(1:n)=[this%perm(1:n-1,nw),n]
          else
             nw=jperm-nperm/2
             jper(1:n)=[this%perm(n-1:1:-1,nw),n]
          endif
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
    
    ! remove the largest one.
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
  end subroutine setup_colmap_cache
  
  subroutine setup_colmap_cache_NLC(this,col_acc)
    use color_algebra
    implicit none
    class(amplitude_cache) :: this
    integer :: col_acc,n,i,jperm,iperm,col_fac,nperm,nw,Q1_start,Q1_end,Q2_start,Q2_end,max_col_index
    integer,dimension(:),allocatable :: ic,ir,iper,jper,permutations_dict2
    integer(kind=8),dimension(:),allocatable :: permutations_dict1
    integer(kind=8) :: val
    real(kind=4) :: tBefore,tAfter

    call cpu_time(tBefore)
    write (*,'(a,i3,a)') ' Setting up colour factor (col_acc =',col_acc,')...'
    n=this%next
    allocate(iper(1:n))
    allocate(jper(1:n))
    nperm=factorial(n-1)
    allocate(permutations_dict1(nperm)) ! ordered dictionary of all the unique numbers for the permutations
    allocate(permutations_dict2(nperm)) ! dictionary that relates the keys of permutations_dict1 to the iperm
    call create_permutations_dict()
    this%col_acc=col_acc
    if (col_acc.ge.1) allocate(this%col_value_NLC(2))
    allocate(this%col_value_LC(1))
    allocate(ic(2))
    allocate(ir(2))
    ! max_col_index is equal to the number of NLC terms. This is
    ! therefore equal to the number of rows in the colour matrix times
    ! the number of non-zero NLC terms in each row. For the
    ! colour-flow basis (which is slightly less efficient than the
    ! fundamental basis) this is equal to the polynomial below. In
    ! fact, since we are only considering the upper-right triangle of
    ! the colour matrix, we can devide the polynomial by two.
    max_col_index=nperm*(n*(2-n-2*n**2+n**3)/24)/2
    if (col_acc.ge.1) allocate(this%col_index_NLC(max_col_index,2))
    if (col_acc.ge.1) allocate(this%row_index_NLC(0:nperm,2))
    allocate(this%col_index_LC(nperm,1))
    allocate(this%row_index_LC(0:nperm,1))
    if (col_acc.ge.1) this%row_index_NLC(0,:)=0
    this%row_index_LC(0,:)=0
    this%col_value_LC(1)= 3**n ! colour factor for the LC contributions
    if (col_acc.ge.1) then ! NLC
       this%col_value_NLC(1)=3**n ! colour factor for the LC contributions
       this%col_value_NLC(2)=3**(n-2) ! colour factor for the NLC contributions
    endif
    ic=0
    ir=0
    do iperm=1,nperm
       if (iperm.le.nperm/2) then
          nw=iperm
          iper(1:n)=[this%perm(1:n-1,nw),n]
       else
          nw=iperm-nperm/2
          iper(1:n)=[this%perm(n-1:1:-1,nw),n]
       endif
       ! LC contribution:
       i=1
       ic(i)=ic(i)+1
       ir(i)=ir(i)+1
       if (col_acc.ge.1) this%col_index_NLC(ic(i),i)=iperm
       this%col_index_LC(ic(i),i)=iperm
       ! NLC contributions (make the block interchange and check if it is a
       ! valid one): Tr[R.Q1.S.Q2.P]*Tr[R.Q2.S.Q1.P]
       jperm=0
       do Q1_start=1,n-2
          do Q1_end=Q1_start,n-2
             do Q2_start=Q1_end+1,n-1
                do Q2_end=Q2_start,n-1
                   jper(1:n)=[iper(1:Q1_start-1),iper(Q2_start:Q2_end),iper(Q1_end+1:Q2_start-1),&
                              iper(Q1_start:Q1_end),iper(Q2_end+1:n)]
                   call compute_color_factor(col_acc,n,iper,jper,col_fac,.true.)
                   if (col_fac.eq.0) cycle
                   call get_value(jper(1:n-1),val)
                   call solve_dict(val,jperm)
                   if (jperm.gt.iperm) cycle
                   do i=1,2
                      if (col_acc.ge.1) then
                         if (this%col_value_NLC(i).eq.col_fac) exit
                      elseif (col_acc.ge.0 .and. i.le.1) then
                         if (this%col_value_LC(i).eq.col_fac) exit
                      endif
                   enddo
                   if (col_acc.eq.0 .and. i.gt.1) cycle
                   ic(i)=ic(i)+1
                   ir(i)=ir(i)+1
                   if (col_acc.ge.1.and.i.le.2) this%col_index_NLC(ic(i),i)=jperm
                   if (i.le.1) this%col_index_LC(ic(i),i)=jperm
                enddo
             enddo
          enddo
       enddo
       if (col_acc.ge.1) this%row_index_NLC(iperm,1:2)=ir(1:2)
       this%row_index_LC(iperm,1)=ir(1)
    enddo
    if (any(ic(:).gt.max_col_index)) then
       write (*,*) 'ERROR: more NLC terms than expected',ic(:),max_col_index
       stop 1
    endif
    call cpu_time(tAfter)
    write (*,*) '... colour setup in',tAfter-tBefore,'seconds'
  contains
    subroutine create_permutations_dict()
      implicit none
      integer :: nw
      do iperm=1,nperm
         if (iperm.le.nperm/2) then
            nw=iperm
            call get_value(this%perm(1:n-1,nw),val)
         else
            nw=iperm-nperm/2
            call get_value(this%perm(n-1:1:-1,nw),val)
         endif
         call find_dict_pos_and_insert()
      enddo
    end subroutine create_permutations_dict
    subroutine find_dict_pos_and_insert()
      ! find the location (i.e. 'key') in the permutations array (which should
      ! be ordered in the val's) to where to insert the current val, and
      ! insert it there. Also update the permutations_dict2, which has as
      ! the same key, but its value is the current 'iperm' being considered.
      implicit none
      integer :: left,right,middle,key
      if (iperm.eq.1) then
         permutations_dict1(iperm)=val
         permutations_dict2(iperm)=iperm
      elseif (iperm.eq.2) then
         if (val.gt.permutations_dict1(iperm-1)) then
            permutations_dict1(iperm)=val
            permutations_dict2(iperm)=iperm
         else
            permutations_dict1(iperm)=permutations_dict1(iperm-1)
            permutations_dict2(iperm)=permutations_dict2(iperm-1)
            permutations_dict1(iperm-1)=val
            permutations_dict2(iperm-1)=iperm
         endif
      else
         left=1
         right=iperm-1
         do
            middle=(right+left)/2
            if (permutations_dict1(middle).lt.val) then
               if (middle.eq.iperm-1) then
                  ! add it at the final position
                  key=middle+1
                  permutations_dict1(key)=val
                  permutations_dict2(key)=iperm
                  return
               elseif (permutations_dict1(middle+1).gt.val) then
                  key=middle+1
                  ! shift all above the key by one position
                  permutations_dict1(key+1:iperm)=permutations_dict1(key:iperm-1)
                  permutations_dict2(key+1:iperm)=permutations_dict2(key:iperm-1)
                  ! add the value
                  permutations_dict1(key)=val
                  permutations_dict2(key)=iperm
                  return
               endif
               left=middle+1
            else
               right=middle
            endif
         enddo
      endif
    end subroutine find_dict_pos_and_insert
    subroutine get_value(ips,val)
      ! Give every permutation 'ips' a unique value.
      implicit none
      integer,dimension(n-1) :: ips
      integer :: j
      integer(kind=8) :: val
      val=0
      ! Give a unique identifier based on colour order. Simply convert the
      ! list to an integer with base equal to the number of external
      ! particles.
      do j=1,n-1
         val=val+int(ips(n-j),kind=8)*int(n,kind=8)**int(j-1,kind=8)
      enddo
    end subroutine get_value
    subroutine solve_dict(val,key)
      ! Given the value 'val', find the corresponding key in the
      ! 'permutations_dict1'. Use that key to return the val of the
      ! 'permutations_dict2' dictionary, which corresponds to the
      ! 'iperm' of the orderx. Use a binary search algorithm. (This
      ! only works if the dictionary values are ordered, and all
      ! values only appear once).
      implicit none
      integer :: key,left,middle,right
      integer(kind=8) :: val
      left=1
      right=nperm
      do while (left.le.right)
         middle=(right+left)/2
         if (permutations_dict1(middle).eq.val) then
            key=permutations_dict2(middle)
            return
         elseif(permutations_dict1(middle).gt.val) then
            right=middle-1
         else
            left=middle+1
         endif
      enddo
      write (*,*) 'value not found in permutations dictionary',val
      stop 1
    end subroutine solve_dict
  end subroutine setup_colmap_cache_NLC
  
  subroutine evaluate_cache(this,p,hel)
    implicit none
    class(amplitude_cache) :: this
    integer,dimension(this%next) :: hel
    real(kind=8),dimension(0:3,this%next) :: p
    integer :: isize,ic,ihel,ivert,iperm,ip,i
    real(kind=8),dimension(:,:),allocatable :: current,pp,current_out
!!$    real(kind=8),dimension(:,:,:),allocatable :: 
    real(kind=8) :: propagator

    if (.not.allocated(current)) allocate(current(1:6,0:this%n_cur_end(this%next-1)))
    if (.not.allocated(pp)) allocate(pp(0:3,maskr(this%next-1)))
    if (.not.allocated(current_out)) allocate(current_out(1:6,1:this%n_vert_end(this%next-1)))

    ! Setup the momenta for all intermediate (and external) particles (use binary labeling)
    do ip=1,maskr(this%next-1)
       pp(0:3,ip)=0d0
       do i=1,2 ! treat incoming momenta as out-going
          if (btest(ip,i-1)) pp(0:3,ip)=pp(0:3,ip)-p(0:3,i)
       enddo
       do i=3,this%next-1
          if (btest(ip,i-1)) pp(0:3,ip)=pp(0:3,ip)+p(0:3,i)
       enddo
    enddo
    
    do isize=1,this%next-1
       if (isize.eq.1) then
          ! this final wave function
          call v_ext(p(0:3,this%next),hel(this%next),1,current(1,0))
          ! fill the other external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)
             if (this%cur_type(ic).ne.0) then
                write (*,*) 'external particle is not a gluon'
                stop
             endif
             ihel=hel(ic)
             call v_ext(pp(0,this%cur_bin(ic)),ihel,1,current(1,ic))
          enddo
          cycle
       endif
       ! compute the interactions
       do ivert=this%n_vert_start(isize),this%n_vert_end(isize)
          if (this%vert_type(ivert).eq.0) then
             call threeGluon(current(1:4,this%vert_cur(1,ivert)),pp(0:3,this%cur_bin(this%vert_cur(1,ivert))),&
                             current(1:4,this%vert_cur(2,ivert)),pp(0:3,this%cur_bin(this%vert_cur(2,ivert))),&
                             current_out(1:4,ivert))
          elseif(this%vert_type(ivert).eq.1) then
             call TwoGluonToTensor(current(1:4,this%vert_cur(1,ivert)),&
                                   current(1:4,this%vert_cur(2,ivert)),&
                                   current_out(1:6,ivert))
          elseif(this%vert_type(ivert).eq.2) then
             call TensorGluontoGluon(current(1:6,this%vert_cur(1,ivert)),&
                                     current(1:4,this%vert_cur(2,ivert)),&
                                     current_out(1:4,ivert))
          elseif(this%vert_type(ivert).eq.3) then
             call GluonTensortoGluon(current(1:4,this%vert_cur(1,ivert)),&
                                     current(1:6,this%vert_cur(2,ivert)),&
                                     current_out(1:4,ivert))
          elseif(this%vert_type(ivert).eq.99) then
             call FourGluon(current(1:4,this%vert_cur(1,ivert)),&
                            current(1:4,this%vert_cur(2,ivert)),&
                            current(1:4,this%vert_cur(3,ivert)),&
                            current_out(1:4,ivert))
          else
             write (*,*) 'Unknown vertex type',ivert,this%vert_type(ivert)
             stop 1
          endif
       enddo
       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          if (this%cur_type(ic).eq.0) then ! gluon current
             current(1:4,ic)=0d0
             do ivert=1,this%cur_n_vert(ic)
                if (.not.this%cur_vert_sign(ivert,ic)) then
                   current(1:4,ic)=current(1:4,ic)+current_out(1:4,this%cur_vertices(ivert,ic))
                else
                   current(1:4,ic)=current(1:4,ic)-current_out(1:4,this%cur_vertices(ivert,ic))
                endif
             enddo
             ! include the gluon propagator
             if (isize.ne.this%next-1)  then
                propagator=1d0/(pp(0,this%cur_bin(ic))**2-pp(1,this%cur_bin(ic))**2- &
                                pp(2,this%cur_bin(ic))**2-pp(3,this%cur_bin(ic))**2)
                current(1:4,ic)=current(1:4,ic)*propagator
             endif
          elseif (this%cur_type(ic).eq.-1) then ! tensor current
             current(1:6,ic)=0d0
             do ivert=1,this%cur_n_vert(ic)
                if (.not.this%cur_vert_sign(ivert,ic)) then
                   current(1:6,ic)=current(1:6,ic)+current_out(1:6,this%cur_vertices(ivert,ic))
                else
                   current(1:6,ic)=current(1:6,ic)-current_out(1:6,this%cur_vertices(ivert,ic))
                endif
             enddo
          else
             write (*,*) 'Unknown current type',ic,this%cur_type(ic)
             stop 1
          endif
       enddo
    enddo

       
    this%nperm=this%n_cur_end(this%next-1)-this%n_cur_start(this%next-1)+1
    
    if (.not.allocated(this%amps)) allocate(this%amps(1:this%nperm*2))
    this%amps(1:this%nperm)=0d0
    do iperm=1,this%nperm    ! permutations
       ip=iperm+this%n_cur_start(this%next-1)-1
       this%amps(iperm)=this%amps(iperm)+ &
            (current(1,ip)*current(1,0)+ &
             current(2,ip)*current(2,0)+ &
             current(3,ip)*current(3,0)+ &
             current(4,ip)*current(4,0))
    enddo
    do iperm=this%nperm+1,this%nperm*2    ! permutations
       ip=iperm-this%nperm
       if (mod(this%next,2).eq.1) then
          this%amps(iperm)=-this%amps(ip)
       else
          this%amps(iperm)=this%amps(ip)
       endif
    enddo

!!$    write (*,*) '3vert',this%amps(1:this%nperm)
!!$    stop
  end subroutine evaluate_cache

  subroutine evaluate(this,p,hel)
    implicit none
    class(amplitude) :: this
    integer,dimension(this%next),optional :: hel
    integer,parameter :: zero=0
    integer :: i,ih,j,isplit,iperm,icount,ih1,ih2,ih3,ip,ji,k,ihel
    real(kind=8),dimension(:,:),allocatable :: pp
    real(kind=8),dimension(0:3,this%next) :: p
    real(kind=8) :: propagator
    real(kind=8),dimension(4) :: wfout
    ! The wavefunctions. The first index is the Lorentz index (goes from 1 to
    ! 4), the second index is for the helicities (goes from 1 to 2^x, where x is
    ! the total number of external gluons combined into this wavefunction) and
    ! the 3rd argument is the wavefunction number (goes from 1 to isize).
    real(kind=8),dimension(:,:,:),allocatable :: wf

    if (.not. allocated(wf)) allocate(wf(4,maskr(this%next-1)+1,0:this%isize))
    if (.not. allocated(pp)) allocate(pp(0:3,0:this%isize))

    ! Main loop over the wavefunction number. This computes all the external and
    ! all off-shell intermediate wavefunctions from already computed ones.
    do i=0,this%isize
       if (i.eq.0) then
          ! The 'final' wavefunction that will be used in the end to close the
          ! amplitude.
          pp(0:3,0)=p(0:3,this%next)
          do ih=1,this%nhel(i)
             if (present(hel)) then
                ihel=hel(this%next)
             else
                ihel=ih-1
             endif
             call v_ext(pp(0,0),ihel,1,wf(1,ih,0))
          enddo
       elseif (i.lt.this%istart(2)) then
          ! The wavefunctions of the external particles but the last.
          if (this%ips(1,i).le.2) then
             pp(0:3,i)=-p(0:3,this%ips(1,i)) ! treat all momenta as outgoing
          else
             pp(0:3,i)=p(0:3,this%ips(1,i))
          endif
          do ih=1,this%nhel(i)
             if (present(hel)) then
                ihel=hel(this%ips(1,i))
             else
                ihel=ih-1
             endif
             call v_ext(pp(0,i),ihel,1,wf(1,ih,i))
          enddo
       else
          if (this%inverted(i).ne.0) then
             ! Already computed the wavefunction with opposite ordering.
             ! Do some bit swapping to map the helicity states correctly.
             do k=0,this%nhel(i)-1
                ji=0
                do j=1,this%ncomb(i)
                   call mvbits(k,j-1,1,ji,this%ncomb(i)-j)
                enddo
                wf(1:4,k+1,i)=sign(1,this%inverted(i))*wf(1:4,ji+1,abs(this%inverted(i)))
             enddo
             pp(0:3,i)=pp(0:3,abs(this%inverted(i)))
             cycle
          endif
          ! from now on in the loop, create new wavefunctions from existing
          ! ones. To understand to which combination the current loop iterator
          ! 'i' belongs, one has to look at the 'setup_imap' subroutine. Start
          ! by initialising the wavefunction to zero.
          wf(1:4,1:this%nhel(i),i)=0d0
          ! combine two wavefunctions to form a third. Loop over the possible splits
          do isplit=1,this%nsplit2(i)
             if (isplit.eq.1) then
                ! define momenta of the newly formed wavefunction. Do this only
                ! for the very first split since this is independent from the
                ! actual split used.
                pp(0:3,i)=pp(0:3,this%imap2(1,isplit,i))+pp(0:3,this%imap2(2,isplit,i))
             endif
             ! double-loop over the helicities of the daughter wavefunctions.
             do ih2=1,this%nhel(this%imap2(2,isplit,i))
                do ih1=1,this%nhel(this%imap2(1,isplit,i))
                   call ThreeGluon(wf(1,ih1,this%imap2(1,isplit,i)),pp(0,this%imap2(1,isplit,i)), &
                                   wf(1,ih2,this%imap2(2,isplit,i)),pp(0,this%imap2(2,isplit,i)), &
                                   wfout)
                   ! helicity label:
                   icount=(ih2-1)*this%nhel(this%imap2(1,isplit,i))+ih1
                   ! add the computed wavefunction to the current wave function
                   wf(1:4,icount,i)=wf(1:4,icount,i)+wfout(1:4)
                enddo
             enddo
          enddo
          ! combine three wavefunctions to form a third
          do isplit=1,this%nsplit3(i)
             ! triple-loop over the helicities of the daughter wavefunctions
             do ih3=1,this%nhel(this%imap3(3,isplit,i))
                do ih2=1,this%nhel(this%imap3(2,isplit,i))
                   icount=(ih3-1)*this%nhel(this%imap3(2,isplit,i))* &
                                  this%nhel(this%imap3(1,isplit,i))+ &
                          (ih2-1)*this%nhel(this%imap3(1,isplit,i))
                   do ih1=1,this%nhel(this%imap3(1,isplit,i))
                      call FourGluon(wf(1,ih1,this%imap3(1,isplit,i)), &
                                     wf(1,ih2,this%imap3(2,isplit,i)), &
                                     wf(1,ih3,this%imap3(3,isplit,i)), &
                                     wfout)
                      ! helicity label:
                      ! add the computed wavefunction to the current wave function
                      wf(1:4,icount+ih1,i)=wf(1:4,icount+ih1,i)+wfout(1:4)
                   enddo
                enddo
             enddo
          enddo
          
          ! compute and include the propagator (except when we are ready to
          ! close with the final wavefunction to form the amplitudes)
          if (i.lt.this%istart(this%next-1))  then
             propagator=1d0/(pp(0,i)**2-pp(1,i)**2-pp(2,i)**2-pp(3,i)**2)
             wf(1:4,1:this%nhel(i),i)= &
                  wf(1:4,1:this%nhel(i),i)*propagator
          endif
       endif
    enddo
    
    ! All wavefunctions have now been computed.
    ! Multiply by the final wavefunction to form the amplitudes:
    this%amps(:,1:this%nperm)=0d0
    do iperm=1,this%nperm    ! permutations
       ip=iperm+this%istart(this%next-1)-1
       do j=1,this%nhel(ip)     ! helicities of combined wavefunctions
          do i=1,this%nhel(0)   ! helicities of final wavefunction
             ih=j-1
             if (i.eq.2) ih=ibset(ih,this%next-1) ! helicity label for amplitude
             this%amps(ih,iperm)=this%amps(ih,iperm)+ &
                           (wf(1,j,ip)*wf(1,i,0)+ &
                            wf(2,j,ip)*wf(2,i,0)+ &
                            wf(3,j,ip)*wf(3,i,0)+ &
                            wf(4,j,ip)*wf(4,i,0))
          enddo
       enddo
    enddo
    deallocate(wf)
    deallocate(pp)
  end subroutine evaluate
   
  subroutine v_ext(p,ihel,ifinal,wf)
    implicit none
    integer ihel,ifinal
    real(kind=8), dimension(0:3) :: p
    real(kind=8), dimension(4) :: wf
    complex*16,dimension(4) :: wf0,wf1
    complex*16,parameter :: cImag=(0d0,1d0)
    call v_ext1(p,1,ifinal,wf1)
    call v_ext1(p,0,ifinal,wf0)
    if (ihel.eq.1) then
       wf(1:4)=dble(cImag*(wf1(1:4)+wf0(1:4)))*sqrt(0.5d0)
    elseif (ihel.eq.0) then
       wf(1:4)=-dble(wf1(1:4)-wf0(1:4))*sqrt(0.5d0)
    endif
  end subroutine v_ext

  subroutine v_ext1(p,ihel,ifinal,wf)
    ! External gluon wavefunction. From HELAS.
    implicit none
    integer :: ihel,ifinal
    real(kind=8), dimension(0:3) :: p
    complex*16, dimension(4) :: wf
    real(kind=8),parameter :: rzero=0d0,sqh=sqrt(0.5d0)
    real(kind=8) :: hel,pp,pt,pzpt
    hel = dble(2*ihel-1)
    pp = p(0)
    pt = sqrt(p(1)**2+p(2)**2)
    wf(1) = dcmplx( rZero )
    wf(4) = dcmplx( hel*pt/pp*sqh )
    if ( pt.ne.rZero ) then
       pzpt = p(3)/(pp*pt)*sqh*hel
       wf(2) = dcmplx( -p(1)*pzpt , -ifinal*p(2)/pt*sqh )
       wf(3) = dcmplx( -p(2)*pzpt ,  ifinal*p(1)/pt*sqh )
    else
       wf(2) = dcmplx( -hel*sqh )
       wf(3) = dcmplx( rZero , ifinal*sign(sqh,p(3)) )
    endif
  end subroutine v_ext1
  
  subroutine ThreeGluon(wf1,pwf1,wf2,pwf2,wf)
    ! Colour-ordered three-gluon interaction
    implicit none
    real(kind=8),dimension(4) :: wf1,wf2,wf
    real(kind=8),dimension(0:3) :: pwf1,pwf2
    real(kind=8),parameter :: prefact=1d0/sqrt(2d0)
    real(kind=8) :: TMP1,TMP2,TMP3
    TMP1 = (wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4))
    TMP2 = (wf1(1)*pwf2(0)-wf1(2)*pwf2(1)-wf1(3)*pwf2(2)-wf1(4)*pwf2(3))
    TMP3 = (wf2(1)*pwf1(0)-wf2(2)*pwf1(1)-wf2(3)*pwf1(2)-wf2(4)*pwf1(3))
    wf(1:4) = prefact*(TMP1*(pwf1(0:3)-pwf2(0:3))+2d0*(TMP2*wf2(1:4)-TMP3*wf1(1:4)))
  end subroutine ThreeGluon

  subroutine FourGluon(wf1,wf2,wf3,wf)
    ! Colour-ordered four-gluon interaction
    implicit none
    real(kind=8),dimension(4) :: wf1,wf2,wf3,wf
    real(kind=8),parameter :: prefact=0.5d0
    real(kind=8) :: TMP1,TMP2,TMP3
    TMP1 = (wf1(1)*wf2(1)-wf1(2)*wf2(2)-wf1(3)*wf2(3)-wf1(4)*wf2(4))
    TMP2 = (wf1(1)*wf3(1)-wf1(2)*wf3(2)-wf1(3)*wf3(3)-wf1(4)*wf3(4))
    TMP3 = (wf2(1)*wf3(1)-wf2(2)*wf3(2)-wf2(3)*wf3(3)-wf2(4)*wf3(4))
    wf(1:4) = prefact*(2d0*wf2(1:4)*TMP2-wf1(1:4)*TMP3-wf3(1:4)*TMP1)
  end subroutine FourGluon

  subroutine TwoGluontoTensor(wfg1,wfg2,wfT)
    ! This vertex includes the all factors such that the tensor "propagator"
    ! is simply the identity
    implicit none
    real(kind=8),dimension(4) :: wfg1,wfg2
    real(kind=8),dimension(6) :: wfT
    ! Since it is an anti-symmetric 4x4 tensor, take only the upper-right triangle.
    wfT(1)=(wfg1(1)*wfg2(2)-wfg1(2)*wfg2(1))
    wfT(2)=(wfg1(1)*wfg2(3)-wfg1(3)*wfg2(1))
    wfT(3)=(wfg1(1)*wfg2(4)-wfg1(4)*wfg2(1))
    wfT(4)=(wfg1(2)*wfg2(3)-wfg1(3)*wfg2(2))
    wfT(5)=(wfg1(2)*wfg2(4)-wfg1(4)*wfg2(2))
    wfT(6)=(wfg1(3)*wfg2(4)-wfg1(4)*wfg2(3))
!!$    do i=1,4
!!$       wfT(1:4,i)=(wfg1(1:4)*wfg2(i)-wfg2(1:4)*wfg1(i))
!!$    enddo
  end subroutine TwoGluontoTensor

  subroutine TensorGluontoGluon(wfT1,wfg2,wfg)
    implicit none
    real(kind=8),dimension(4) :: wfg2,wfg
    real(kind=8),dimension(6) :: wfT1
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1)=(wfT1(1)*wfg2(2)+wfT1(2)*wfg2(3)+wfT1(3)*wfg2(4))*prefact
    wfg(2)=(wfT1(1)*wfg2(1)+wfT1(4)*wfg2(3)+wfT1(5)*wfg2(4))*prefact
    wfg(3)=(wfT1(2)*wfg2(1)-wfT1(4)*wfg2(2)+wfT1(6)*wfg2(4))*prefact
    wfg(4)=(wfT1(3)*wfg2(1)-wfT1(5)*wfg2(2)-wfT1(6)*wfg2(3))*prefact
!!$    do i=1,4
!!$       wfg(i)=((wfT1(1,i)*wfg2(1)-wfT1(2,i)*wfg2(2)-wfT1(3,i)*wfg2(3)-wfT1(4,i)*wfg2(4))- &
!!$               (wfT1(i,1)*wfg2(1)-wfT1(i,2)*wfg2(2)-wfT1(i,3)*wfg2(3)-wfT1(i,4)*wfg2(4)))*0.25d0
!!$    enddo
  end subroutine TensorGluontoGluon

  subroutine GluonTensortoGluon(wfg1,wfT2,wfg)
    implicit none 
    real(kind=8),dimension(4) :: wfg1,wfg
    real(kind=8),dimension(6) :: wfT2
    real(kind=8),parameter :: prefact=0.5d0
    wfg(1)=(-wfg1(2)*wfT2(1)-wfg1(3)*wfT2(2)-wfg1(4)*wfT2(3))*prefact
    wfg(2)=(-wfg1(1)*wfT2(1)-wfg1(3)*wfT2(4)-wfg1(4)*wfT2(5))*prefact
    wfg(3)=(-wfg1(1)*wfT2(2)+wfg1(2)*wfT2(4)-wfg1(4)*wfT2(6))*prefact
    wfg(4)=(-wfg1(1)*wfT2(3)+wfg1(2)*wfT2(5)+wfg1(3)*wfT2(6))*prefact
!!$    do i=1,4
!!$       wfg(i)=-((wfg1(1)*wfT2(1,i)-wfg1(2)*wfT2(2,i)-wfg1(3)*wfT2(3,i)-wfg1(4)*wfT2(4,i))- &
!!$               (wfg1(1)*wfT2(i,1)-wfg1(2)*wfT2(i,2)-wfg1(3)*wfT2(i,3)-wfg1(4)*wfT2(i,4)))*0.25d0
!!$    enddo
  end subroutine GluonTensortoGluon

  subroutine setup_colmap_CSR(this,n,order)
    use color_algebra
    implicit none
    class(amplitude) :: this
    integer :: n,i,jperm,iperm,col_fac,imax,max_val
    integer,dimension(n),optional :: order
    logical,parameter :: color_flow=.true.
    integer,dimension(n) :: iper,jper
    integer,dimension(:),allocatable :: ic,ir
    if (present(order) .or. .not.color_flow) then
       write (*,*) 'CSR format only available when summing over all '// &
            'colour orders using color flow basis'
       stop 1
    endif
    allocate(this%col_value((n+1)/2))
    allocate(ic((n+1)/2))
    allocate(ir((n+1)/2))
    allocate(this%col_index(this%nperm**2,(n+1)/2))
    allocate(this%row_index(0:this%nperm,(n+1)/2))
    this%row_index(0,:)=0
    do i=1,(n+1)/2
       this%col_value(i)=3**(n-2*(i-1))
    enddo
    do i=1,n
       iper(i)=i
       jper(i)=i
    enddo
    ic=0
    ir=0
    do iperm=1,this%nperm
       do jperm=1,this%nperm
          if (jperm.ge.iperm) then
             call compute_color_factor(this%col_acc,n,iper,jper,col_fac,color_flow)
          endif
          call ipnext(jper,n-1)
          if (jperm.ge.iperm) then
             if (col_fac.eq.0) cycle
             do i=1,(n+1)/2
                if (this%col_value(i).eq.col_fac) exit
             enddo
             ic(i)=ic(i)+1
             ir(i)=ir(i)+1
             this%col_index(ic(i),i)=jperm
          endif
       enddo
       this%row_index(iperm,:)=ir(:)
       call ipnext(iper,n-1)
    enddo

    ! remove the largest one.
    imax=0
    max_val=0
    do i=1,(n+1)/2
       if (this%row_index(this%nperm,i).gt.max_val) then
          max_val=max(this%row_index(this%nperm,i),max_val)
          imax=i
       endif
    enddo
    if (all(this%row_index(this%nperm,:).ne.0)) then
       this%row_index(:,imax)=0
       this%col_value(:)=this%col_value(:)-this%col_value(imax)
    endif
  end subroutine setup_colmap_CSR

  subroutine setup_colmap(this,n,order)
    use color_algebra
    implicit none
    class(amplitude) :: this
    integer :: i,j,jperm,iperm,col_fac,n,n_col_one_row
    integer,dimension(n),optional :: order
    integer,dimension(n) :: iper,jper,permutation
    logical,parameter :: color_flow=.true.
    integer,dimension(:,:),allocatable :: temp2,permap
!!$    if (.not. (present(order) .or. (color_flow .and. this%col_acc.ge.2))) then
    if (.not. present(order) ) then
       if (.not.color_flow .and. this%col_acc.ge.2)  call Tr_allocate(n)
       ! allocate enough for one row in the colour matrix
       if (.not. allocated(this%colmap)) allocate(this%colmap(0:2,0:this%nperm))
       if (.not. allocated(permap)) allocate(permap(n,this%nperm))
       do i=1,n
          iper(i)=i
          jper(i)=i
       enddo
       ! create the one row of the colour matrix
       this%colmap(0,0)=0
       jperm=0
       do iperm=1,this%nperm
          call compute_color_factor(this%col_acc,n,iper,jper,col_fac,color_flow)
          if (col_fac.ne.0 .or. iperm.eq.1) then  ! include the iperm==1 case, since col_fac is zero for diagonal in fundamental basis at NLC
             this%colmap(0,0)=this%colmap(0,0)+1
             this%colmap(0,this%colmap(0,0))=col_fac
             this%colmap(1,this%colmap(0,0))=iperm
             permap(1:n,this%colmap(0,0))=iper(1:n)
             if (all(jper.eq.iper)) jperm=iperm
          endif
          call ipnext(iper,n-1)
       enddo
       if (jperm.eq.1) then
          this%colmap(2,1:this%colmap(0,0))=jperm
       else
          write (*,*) 'jperm not found',jperm
          write (*,*) 'requires improvements in the method'
          stop 1
       endif
       n_col_one_row=this%colmap(0,0)
       ! One row done.
       ! Now fill the rest by only checking the non-zero ones found and the
       ! relevant permutations.
       allocate(temp2(0:2,0:this%colmap(0,0)*this%nperm))
       temp2(0:2,0:this%colmap(0,0))=this%colmap(0:2,0:this%colmap(0,0))
       call move_alloc(temp2,this%colmap)
       do jperm=2,this%nperm
          call ipnext(jper,n-1)
          do j=1,n
             do i=1,n
                if (jper(i).eq.j) then
                   permutation(i)=j
                   exit
                endif
             enddo
          enddo
          do i=1,n_col_one_row
             call iperm_encode(n-1,permutation(permap(1:n-1,i)),iperm)
                this%colmap(0,0)=this%colmap(0,0)+1
                this%colmap(0,this%colmap(0,0))=this%colmap(0,i)
                this%colmap(1,this%colmap(0,0))=iperm
                this%colmap(2,this%colmap(0,0))=jperm
          enddo
       enddo
    else
       call old_setup_colmap(this,n,order,color_flow)
    endif
  end subroutine setup_colmap



  subroutine old_setup_colmap(this,n,order,color_flow)
    use color_algebra
    implicit none
    class(amplitude) :: this
    integer :: i,jperm,iperm,col_fac,n,acc
    integer,dimension(n),optional :: order
    integer,dimension(n) :: iper,jper
    character*(this%col_acc+1) :: buff
    logical :: color_flow
    real(kind=16) :: col_factor
    integer,dimension(:,:),allocatable :: temp2

    if (.not.color_flow .and. this%col_acc.ge.2)  call Tr_allocate(n)
    
    if (present(order)) then
       if (this%nperm.ne.1) then
          write (*,*) 'Fixed colour order, but not just 1 permutation',this%nperm
          stop 1
       endif
       if (.not. allocated(this%colmap)) allocate(this%colmap(0:2,0:1))
       col_fac=3**n ! leading colour factor
       this%colmap(0,0)=1
       this%colmap(0,1)=col_fac
       this%colmap(1,1)=1
       this%colmap(2,1)=1
    else
       ! allocate enough for one row in the colour matrix
       if (.not. allocated(this%colmap)) allocate(this%colmap(0:2,0:this%nperm))
       do i=1,n
          iper(i)=i
          jper(i)=i
       enddo
       this%colmap(0,0)=0
       do jperm=1,this%nperm
          do iperm=1,this%nperm
             if (iperm.ge.jperm) then
                ! include only leading and dominant subleading colour
                ! terms (to be set by col_acc in the init subroutine):
                if (color_flow) then
                   call color_flow_factor(n,jper,iper,col_fac)
                   if (col_fac.ge.n-2*this%col_acc) then
                      col_fac=3**col_fac
                   else
                      col_fac=0
                   endif
                else
                   if (this%col_acc.eq.0) then ! LC
                      if (iperm.eq.jperm) then
                         col_fac=3**n
                      else
                         col_fac=0
                      endif
                   elseif (this%col_acc.eq.1) then ! NLC
                      if (iperm.eq.jperm) then
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
                      do i=n,max(n-2*this%col_acc,0),-1 ! do not include any Nc
                                                        ! contributions with negative
                                                        ! powers, since they must cancel.
                         col_fac=col_fac+coef_nc(i,0)*3**i
                      enddo
                   endif
                endif
                if (col_fac.ne.0) then
                   this%colmap(0,0)=this%colmap(0,0)+1
                   if (iperm.ne.jperm) col_fac=col_fac *2
                   this%colmap(0,this%colmap(0,0))=col_fac
                   this%colmap(1,this%colmap(0,0))=iperm
                   this%colmap(2,this%colmap(0,0))=jperm
                endif
             endif
             call ipnext(iper,n-1)
          enddo
          if (jperm.eq.1) then
             ! allocate the full colmap
             allocate(temp2(0:2,0:this%colmap(0,0)*this%nperm))
             temp2(0:2,0:this%colmap(0,0))=this%colmap(0:2,0:this%colmap(0,0))
             call move_alloc(temp2,this%colmap)
          endif
          call ipnext(jper,n-1)
       enddo
       buff=''
       do i=1,this%col_acc
          buff(i:i)='N'
       enddo
       write (*,*) '- '//buff(1:this%col_acc)//'LC: Including ',this%colmap(0,0), &
            'colour terms out of ',this%nperm**2/2+this%nperm/2, &
            '(',this%colmap(0,0)/dble(this%nperm**2/2+this%nperm/2)*100d0,'%)'
    endif
  end subroutine old_setup_colmap


  
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
  
  subroutine setup_helmap(this,n)
    implicit none
    class(amplitude) :: this
    integer :: ih,iperm,ih1,ib,n
    integer,dimension(n) ::ipss
    if (this%sum_hel) then
       if (.not. allocated(this%helmap)) allocate(this%helmap(this%nperm,0:maskr(n)))
       do ih=0,maskr(n)
          do iperm=1,this%nperm
             ipss(1:n-1)=this%ips(1:n-1,iperm+this%istart(n-1)-1)
             ipss(n)=n
             ih1=0
             do ib=1,n
                call mvbits(ih,ipss(ib)-1,1,ih1,ib-1)
             enddo
             this%helmap(iperm,ih)=ih1
          enddo
       enddo
    else
       if (.not. allocated(this%helmap)) allocate(this%helmap(this%nperm,0:0))
       this%helmap(1:this%nperm,0)=0
    endif
  end subroutine setup_helmap

  integer function factorial(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    integer,dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*i
       enddo
       ifact_save=ifact
    endif
    factorial=factorial_save(ifact)
  end function factorial

  integer(kind=8) function factorial8(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    integer(kind=8),dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*i
       enddo
       ifact_save=ifact
    endif
    factorial8=factorial_save(ifact)
  end function factorial8

  real(kind=8) function factorial_dble(ifact)
    ! computes the factorial of 'ifact'.
    implicit none
    integer, value :: ifact
    integer :: i
    integer,save :: ifact_save=-1
    real(kind=8),dimension(:),allocatable,save :: factorial_save
    if (ifact.gt.ifact_save) then
       if (allocated(factorial_save)) deallocate(factorial_save)
       allocate(factorial_save(0:ifact))
       factorial_save(0)=1d0
       do i=1,ifact
          factorial_save(i)=factorial_save(i-1)*dble(i)
       enddo
       ifact_save=ifact
    endif
    factorial_dble=factorial_save(ifact)
  end function factorial_dble


  subroutine evaluate_order_v3(this,icol,next,p,order,hel,amp)
    ! simple algorithm to compute the all-gluon amplitude for a given
    ! helicity and color-order configuration
    implicit none
    class(col_amp) :: this
    integer :: next,i,ilength,ipart,ipart2,icol
    real(kind=8),dimension(0:3,1:next) :: p
    integer,dimension(1:next) :: hel,order
    real(kind=8) :: amp,propagator
    real(kind=8),dimension(1:4) :: wfout
    this%n=next
    allocate(this%order(1:this%n))
    this%order(1:this%n)=order(1:this%n)
    allocate(this%wf(1:4,this%n-1,this%n-1))
    allocate(this%pp(0:3,this%n-1,this%n))
    do ilength=1,this%n-1
       if (ilength.eq.1) then
          do i=1,this%n-ilength
             if (check_reuse(this,icol,ilength,i)) cycle
             if (this%order(i).le.2) then
                this%pp(0:3,1,i)=-p(0:3,this%order(i)) ! treat all momenta as outgoing
             else
                this%pp(0:3,1,i)=p(0:3,this%order(i))
             endif
             call v_ext(this%pp(0,1,i),hel(this%order(i)),1,this%wf(1,1,i))
          enddo
       else
          do i=1,this%n-ilength
             if (check_reuse(this,icol,ilength,i)) cycle
             this%wf(1:4,ilength,i)=0d0
             do ipart=1,ilength-1
                call ThreeGluon(this%wf(1,ipart,i),this%pp(0,ipart,i), &
                                this%wf(1,ilength-ipart,i+ipart),this%pp(0,ilength-ipart,i+ipart),wfout)
                this%wf(1:4,ilength,i)=this%wf(1:4,ilength,i)+wfout(1:4)
                do ipart2=1,ilength-ipart-1
                   call FourGluon(this%wf(1,ipart,i),this%wf(1,ipart2,i+ipart), &
                                  this%wf(1,ilength-ipart-ipart2,i+ipart+ipart2),wfout)
                   this%wf(1:4,ilength,i)=this%wf(1:4,ilength,i)+wfout(1:4)
                enddo
                if (ipart.eq.1) then
                   this%pp(0:3,ilength,i)=this%pp(0:3,ipart,i)+this%pp(0:3,ilength-ipart,i+ipart)
                endif
             enddo
             if (ilength.ne.this%n-1) then
                propagator=1d0/(this%pp(0,ilength,i)**2-this%pp(1,ilength,i)**2 &
                               -this%pp(2,ilength,i)**2-this%pp(3,ilength,i)**2)
                this%wf(1:4,ilength,i)=this%wf(1:4,ilength,i)*propagator
             endif
          enddo
       endif
    enddo
    if (this%order(this%n).le.2) then
       this%pp(0:3,1,this%n)=-p(0:3,this%order(this%n)) ! treat all momenta as outgoing
    else
       this%pp(0:3,1,this%n)=p(0:3,this%order(this%n))
    endif
    call v_ext(this%pp(0,1,this%n),hel(this%order(this%n)),1,wfout)
    amp=this%wf(1,this%n-1,1)*wfout(1)+this%wf(2,this%n-1,1)*wfout(2)+&
         this%wf(3,this%n-1,1)*wfout(3)+this%wf(4,this%n-1,1)*wfout(4)

  end subroutine evaluate_order_v3
  
  logical function check_reuse(this,icol,ilength,i)
    implicit none
    class(col_amp) :: this
    integer :: icol,ilength,i
    integer :: j,jj
    do j=1,icol-1
       do jj=1,this%n-ilength
          if (all(col_amp_list(j)%order(jj:jj+ilength-1).eq.this%order(i:i+ilength-1))) then
             check_reuse=.true.
             this%wf(1:4,ilength,i)=col_amp_list(j)%wf(1:4,ilength,jj)
             this%pp(0:3,ilength,i)=col_amp_list(j)%pp(0:3,ilength,jj)
             return
          endif
       enddo
    enddo
    check_reuse=.false.
  end function check_reuse
    
  subroutine evaluate_order_v2(n,p,order,hel,amp)
    ! simple algorithm to compute the all-gluon amplitude for a given
    ! helicity and color-order configuration
    implicit none
    integer :: n,i,ilength,ipart,ipart2
    real(kind=8),dimension(0:3,1:n) :: p
    integer,dimension(1:n) :: hel,order
    real(kind=8) :: amp,propagator
    real(kind=8),dimension(0:3,n-1,n) :: pp
    real(kind=8),dimension(1:4,n-1,n-1) :: wf
    real(kind=8),dimension(1:4) :: wfout
    do ilength=1,n-1
       if (ilength.eq.1) then
          do i=1,n-ilength
             if (order(i).le.2) then
                pp(0:3,1,i)=-p(0:3,order(i)) ! treat all momenta as outgoing
             else
                pp(0:3,1,i)=p(0:3,order(i))
             endif
             call v_ext(pp(0,1,i),hel(order(i)),1,wf(1,1,i))
          enddo
       else
          do i=1,n-ilength
             wf(1:4,ilength,i)=0d0
             do ipart=1,ilength-1
                call ThreeGluon(wf(1,ipart,i),pp(0,ipart,i), &
                                wf(1,ilength-ipart,i+ipart),pp(0,ilength-ipart,i+ipart),wfout)
                wf(1:4,ilength,i)=wf(1:4,ilength,i)+wfout(1:4)
                do ipart2=1,ilength-ipart-1
                   call FourGluon(wf(1,ipart,i),wf(1,ipart2,i+ipart), &
                                  wf(1,ilength-ipart-ipart2,i+ipart+ipart2),wfout)
                   wf(1:4,ilength,i)=wf(1:4,ilength,i)+wfout(1:4)
                enddo
                if (ipart.eq.1) then
                   pp(0:3,ilength,i)=pp(0:3,ipart,i)+pp(0:3,ilength-ipart,i+ipart)
                endif
             enddo
             if (ilength.ne.n-1) then
                propagator=1d0/(pp(0,ilength,i)**2-pp(1,ilength,i)**2 &
                               -pp(2,ilength,i)**2-pp(3,ilength,i)**2)
                wf(1:4,ilength,i)=wf(1:4,ilength,i)*propagator
             endif
          enddo
       endif
    enddo
    if (order(n).le.2) then
       pp(0:3,1,n)=-p(0:3,order(n)) ! treat all momenta as outgoing
    else
       pp(0:3,1,n)=p(0:3,order(n))
    endif
    call v_ext(pp(0,1,n),hel(order(n)),1,wfout)
    amp=wf(1,n-1,1)*wfout(1)+wf(2,n-1,1)*wfout(2)+wf(3,n-1,1)*wfout(3)+wf(4,n-1,1)*wfout(4)
  end subroutine evaluate_order_v2

    
  subroutine evaluate_order(this,n,p,order,hel,amp2)
    ! recursive algorithm to compute the all-gluon amplitude squared for
    ! a given helicity and color-order configuration
    use ext_wfs
    implicit none
    class(amplitude) :: this
    integer :: n,i
    real(kind=8),dimension(0:3,1:n) :: p
    integer,dimension(1:n) :: hel,order
    real(kind=8) :: amp2
    real(kind=8),dimension(0:3) :: pp
    real(kind=8),dimension(1:4) :: wf
    next=n
    if (.not.allocated(wf_ext)) allocate(wf_ext(4,next))
    if (.not.allocated(p_ext)) allocate(p_ext(0:3,next))
    do i=1,n
       if (i.le.2) then
          p_ext(0:3,i)=-p(0:3,i) ! treat all momenta as outgoing
       else
          p_ext(0:3,i)=p(0:3,i)
       endif
       call v_ext(p_ext(0,i),hel(i),1,wf_ext(1,i))
    enddo
    call eval_order(n-1,order(1:n-1), wf,pp)
    amp2=wf(1)*wf_ext(1,order(next))+wf(2)*wf_ext(2,order(next))+wf(3)*wf_ext(3,order(next))+wf(4)*wf_ext(4,order(next))
    amp2=amp2**2
  end subroutine evaluate_order

  recursive subroutine eval_order(len,order, wf,pp)
    use ext_wfs
    implicit none
    integer,intent(in) :: len
    integer,dimension(len),intent(in) :: order
    real(kind=8),dimension(0:3),intent(out) :: pp
    real(kind=8),dimension(4),intent(out) :: wf
    integer :: i,j
    real(kind=8),dimension(0:3) :: pp1,pp2,pp3
    real(kind=8),dimension(4) :: wf1,wf2,wf3,wfout
    real(kind=8) :: propagator
    if (len.eq.1) then
       wf(1:4)=wf_ext(1:4,order(1))
       pp(0:3)=p_ext(0:3,order(1))
       return
    endif
    wf(1:4)=0d0
    do i=1,len-1
       call eval_order(i,order(1:i),wf1,pp1)
       call eval_order(len-i,order(i+1:len),wf2,pp2)
       call ThreeGluon(wf1,pp1,wf2,pp2,wfout)
       wf(1:4)=wf(1:4)+wfout(1:4)
       if (i.eq.1) then
          pp(0:3)=pp1(0:3)+pp2(0:3)
       endif
       do j=i+1,len-1
          call eval_order(j-i,order(i+1:j),wf2,pp2)
          call eval_order(len-j,order(j+1:len),wf3,pp3)
          call FourGluon(wf1,wf2,wf3,wfout)
          wf(1:4)=wf(1:4)+wfout(1:4)
       enddo
    enddo
    if (len.lt.next-1) then
       propagator=1d0/(pp(0)**2-pp(1)**2-pp(2)**2-pp(3)**2)
       wf(1:4)=wf(1:4)*propagator
    endif
  end subroutine eval_order

  subroutine setup_imap_cache(this,decompose_4vert,next)
    implicit none
    class(amplitude_cache) :: this
    logical :: decompose_4vert
    integer :: n_cur,n_vert,nc,isize,next,n1,n2,n3,isplit,isplit2,ic1,ic2,ic3,max_cur,max_vert
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:,:),allocatable :: cur_part
    real(kind=4) :: tBefore,tAfter
    call cpu_time(tBefore)
    if (decompose_4vert) then
       write (*,*) 'setup imap with only 3-vertices...'
    else
       write (*,*) 'setup imap...'
    endif
    this%next=next
    allocate(this%n_cur_start(next-1))
    allocate(this%n_cur_end(next-1))
    allocate(this%n_vert_start(next-1))
    allocate(this%n_vert_end(next-1))
    call set_max_cur()
    call set_max_vert()
    allocate(current_dict(max_cur))
    allocate(cur_part(next-1,max_cur))
    allocate(this%cur_type(max_cur))
    allocate(this%cur_bin(max_cur))
    allocate(this%cur_n_vert(max_cur))
    this%cur_n_vert(1:max_cur)=0
    allocate(this%vert_type(max_vert))
    if (decompose_4vert) then
       allocate(this%vert_cur(2,max_vert))
    else
       allocate(this%vert_cur(3,max_vert))
    endif
    if (decompose_4vert) then
       ! need 3*(next-2), since we can combine gluon-gluon, tensor-gluon, and gluon-tensor to a gluon 
       allocate(this%cur_vertices(3*(next-2),max_cur)) 
       allocate(this%cur_vert_sign(3*(next-2),max_cur))
    else
       ! need (next-2) for 3-gluon and another (next-2)*(next-3)/2 for 4-gluon
       allocate(this%cur_vertices((next-2)+(next-2)*(next-3)/2,max_cur)) 
       allocate(this%cur_vert_sign((next-2)+(next-2)*(next-3)/2,max_cur))
    endif
    ! create a dictionary with all currents to be able to quickly find them in the list.
    call create_current_dict()
    n_cur=0  ! number of currents
    n_vert=0 ! number of vertices
    do isize=1,next-1
       call cpu_time(tAfter)
       write (*,*) '   isize',isize,tAfter-tBefore
       this%n_cur_start(isize)=n_cur+1
       this%n_vert_start(isize)=n_vert+1
       if (isize.eq.1) then
          ! external particles
          do nc=1,next-1
             n_cur=n_cur+1
             this%cur_type(n_cur)=0    ! gluon
             cur_part(1,n_cur)=nc
             this%cur_bin(n_cur)=ibset(0,nc-1) ! binary labeling for external particles included in this current
          enddo
       else
          ! try any combination of two (or three) previously computed currents
          ! that can give a current of size 'isize'
          do isplit=1,isize-1
             n1=isplit
             n2=isize-isplit
             do ic1=this%n_cur_start(n1),this%n_cur_end(n1)
                do ic2=this%n_cur_start(n2),this%n_cur_end(n2)
                   call add_if_allowed_threevertex()
                enddo
             enddo
             if (decompose_4vert .or. isize.le.2 .or. isplit.eq.isize-1) cycle
             do isplit2=isplit+1,isize-1
                n2=isplit2-isplit
                n3=isize-isplit2
                do ic1=this%n_cur_start(n1),this%n_cur_end(n1)
                   do ic2=this%n_cur_start(n2),this%n_cur_end(n2) 
                      do ic3=this%n_cur_start(n3),this%n_cur_end(n3)
                         call add_if_allowed_fourvertex()
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif
       this%n_cur_end(isize)=n_cur
       this%n_vert_end(isize)=n_vert
    enddo

    if (n_cur.ne.max_cur) then
       write (*,*) 'ERROR: amount of currents not equal to the expect amount',n_cur,max_cur
       stop 1
    endif
    if (n_vert.ne.max_vert) then
       write (*,*) 'ERROR: amount of vertices not equal to the expect amount',n_vert,max_vert
       stop 1
    endif
    
    call cpu_time(tAfter)
    write (*,*) '   isize',isize,tAfter-tBefore
    allocate(this%perm(1:next-1,1:(this%n_cur_end(next-1)-this%n_cur_start(next-1)+1)))
    this%perm(1:next-1,1:this%n_cur_end(next-1)-this%n_cur_start(next-1)+1)=&
         cur_part(1:next-1,this%n_cur_start(next-1):this%n_cur_end(next-1))
    write (*,*) '   total number of currents and vertices',n_cur,n_vert
    call cpu_time(tAfter)
    write (*,*) '... imap setup in ',tAfter-tBefore,'seconds'
  contains
    subroutine add_if_allowed_fourvertex()
      ! same as add_if_alllowed_threevertex, but for 4-vertices.
      implicit none
      ! only consider one ordering; the other will be obtained from symmetry:
      if (maxval(cur_part(1:n1,ic1)).ge.maxval(cur_part(1:n3,ic3))) return
      ! check that all particles are different in the three currents:
      if (popcnt(iparity([this%cur_bin(ic1),this%cur_bin(ic2),this%cur_bin(ic3)])).ne.isize) return
      ! check that types form a valid vertex. If so, add it to the list.
      if (this%cur_type(ic1).eq.0 .and. this%cur_type(ic2).eq.0 .and. this%cur_type(ic3).eq.0) then
         ! add a gluon-gluon-gluon to gluon vertex
         call add_valid_vertex(99)
         call add_all_4vert_to_currents()
      endif
    end subroutine add_if_allowed_fourvertex
    subroutine add_if_allowed_threevertex()
      ! check if we should consider the current combination, and if
      ! so, and the corresponding vertices to the list. Once the
      ! vertices are added, we need to check all the currents to which
      ! this vertex contributions and add it to all of them (using the
      ! 'add_all_3vert_to_currents()' subroutine).
      implicit none
      ! only consider one ordering; the other will be obtained from symmetry:
      if (maxval(cur_part(1:n1,ic1)).ge.maxval(cur_part(1:n2,ic2))) return
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(this%cur_bin(ic1),this%cur_bin(ic2))).ne.isize) return
      ! check that types form a valid vertex. If so, add it to the list.
      if (this%cur_type(ic1).eq.0 .and. this%cur_type(ic2).eq.0) then
         ! add a gluon-gluon to gluon vertex
         call add_valid_vertex(0)
         call add_all_3vert_to_currents()
         if (isize.ne.next-1 .and. decompose_4vert) then
            ! add a gluon-gluon to tensor vertex
            call add_valid_vertex(1)
            call add_all_3vert_to_currents()
         endif
      elseif (this%cur_type(ic1).eq.-1 .and. this%cur_type(ic2).eq.0) then
         ! add a tensor-gluon to gluon vertex
         call add_valid_vertex(2)
         call add_all_3vert_to_currents()
      elseif (this%cur_type(ic1).eq.0 .and. this%cur_type(ic2).eq.-1) then
         ! add a gluon-tensor to gluon vertex
         call add_valid_vertex(3)
         call add_all_3vert_to_currents()
      endif
    end subroutine add_if_allowed_threevertex
    subroutine add_valid_vertex(itype)
      ! We have a valid vertex, so add it to the list
      implicit none
      integer :: itype
      n_vert=n_vert+1
      this%vert_type(n_vert)=itype
      this%vert_cur(1,n_vert)=ic1
      this%vert_cur(2,n_vert)=ic2
      if (.not.decompose_4vert .and. itype.eq.99) this%vert_cur(3,n_vert)=ic3
    end subroutine add_valid_vertex
    subroutine add_all_4vert_to_currents()
      ! same as add_all_3vert_to_currents, but for 4-vertices
      implicit none
      logical :: vertex_sign
      integer,dimension(isize,16) :: ip
      integer :: i,cur_bin
      ! need to consider 16 possible permutations (1, 2, 4, or 8 permutations will be a valid order)
      ip(1:isize, 1)=[cur_part(1:n1   ,ic1),cur_part(1:n2   ,ic2),cur_part(1:n3   ,ic3)]
      ip(1:isize, 2)=[cur_part(1:n3   ,ic3),cur_part(1:n2   ,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize, 3)=[cur_part(n1:1:-1,ic1),cur_part(1:n2   ,ic2),cur_part(1:n3   ,ic3)]
      ip(1:isize, 4)=[cur_part(1:n3   ,ic3),cur_part(1:n2   ,ic2),cur_part(n1:1:-1,ic1)]
      ip(1:isize, 5)=[cur_part(1:n1   ,ic1),cur_part(n2:1:-1,ic2),cur_part(1:n3   ,ic3)]
      ip(1:isize, 6)=[cur_part(1:n3   ,ic3),cur_part(n2:1:-1,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize, 7)=[cur_part(1:n1   ,ic1),cur_part(1:n2   ,ic2),cur_part(n3:1:-1,ic3)]
      ip(1:isize, 8)=[cur_part(n3:1:-1,ic3),cur_part(1:n2   ,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize, 9)=[cur_part(n1:1:-1,ic1),cur_part(n2:1:-1,ic2),cur_part(1:n3   ,ic3)]
      ip(1:isize,10)=[cur_part(1:n3   ,ic3),cur_part(n2:1:-1,ic2),cur_part(n1:1:-1,ic1)]
      ip(1:isize,11)=[cur_part(n1:1:-1,ic1),cur_part(1:n2   ,ic2),cur_part(n3:1:-1,ic3)]
      ip(1:isize,12)=[cur_part(n3:1:-1,ic3),cur_part(1:n2   ,ic2),cur_part(n1:1:-1,ic1)]
      ip(1:isize,13)=[cur_part(1:n1   ,ic1),cur_part(n2:1:-1,ic2),cur_part(n3:1:-1,ic3)]
      ip(1:isize,14)=[cur_part(n3:1:-1,ic3),cur_part(n2:1:-1,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize,15)=[cur_part(n1:1:-1,ic1),cur_part(n2:1:-1,ic2),cur_part(n3:1:-1,ic3)]
      ip(1:isize,16)=[cur_part(n3:1:-1,ic3),cur_part(n2:1:-1,ic2),cur_part(n1:1:-1,ic1)]
      ! the binary label for the external particles included in this current
      cur_bin=this%cur_bin(ic1)+this%cur_bin(ic2)+this%cur_bin(ic3)
      do i=1,16
         if (n1.eq.1 .and. (i.eq.3 .or. i.eq.4 .or. i.eq.9 .or. i.eq.10 .or. i.eq.11 .or. i.eq.12 .or. i.eq.15 .or. i.eq.16)) cycle
         if (n2.eq.1 .and. (i.eq.5 .or. i.eq.6 .or. i.eq.9 .or. i.eq.10 .or. i.eq.13 .or. i.eq.14 .or. i.eq.15 .or. i.eq.16)) cycle
         if (n3.eq.1 .and. (i.eq.7 .or. i.eq.8 .or. i.eq.11 .or. i.eq.12 .or. i.eq.13 .or. i.eq.14 .or. i.eq.15 .or. i.eq.16)) cycle
         if (valid_current_order(ip(1:isize,i))) then
            if (i.eq.1 .or. i.eq.2 .or. &
                 (i.eq.3 .and. mod(n1,2).eq.1)    .or. (i.eq.4 .and. mod(n1,2).eq.1) .or. &
                 (i.eq.5 .and. mod(n2,2).eq.1)    .or. (i.eq.6 .and. mod(n2,2).eq.1) .or. &
                 (i.eq.7 .and. mod(n3,2).eq.1)    .or. (i.eq.8 .and. mod(n3,2).eq.1) .or. &
                 (i.eq.9 .and. mod(n1+n2,2).eq.0) .or. (i.eq.10.and. mod(n1+n2,2).eq.0) .or. &
                 (i.eq.11.and. mod(n1+n3,2).eq.0) .or. (i.eq.12.and. mod(n1+n3,2).eq.0) .or. &
                 (i.eq.13.and. mod(n2+n3,2).eq.0) .or. (i.eq.14.and. mod(n2+n3,2).eq.0) .or. &
                 (i.eq.15.and. mod(isize,2).eq.1) .or. (i.eq.16.and. mod(isize,2).eq.1) ) then
               vertex_sign=.false. ! no extra sign needed
            else
               vertex_sign=.true.  ! permutation requires a minus sign
            endif
            call add_one_to_currents(vertex_sign,cur_bin,ip(1:isize,i))
         endif
      enddo
    end subroutine add_all_4vert_to_currents
    subroutine add_all_3vert_to_currents()
      ! Given the vertex, we have to check all the permutations and
      ! add them to the corresponding currents. That is, if a vertex
      ! permutation contributes to a valid current, add that to that
      ! current with the 'add_one_to_currents()' subroutine. Also keep
      ! track of the sign: some permutations require a minus sign.
      implicit none
      logical :: vertex_sign
      integer,dimension(isize,8) :: ip
      integer :: i,cur_bin
      ! Need to consider the 8 possible permutations (1, 2 or 4 permutations will actually be a valid order)
      ip(1:isize,1)=[cur_part(1:n1   ,ic1),cur_part(1:n2   ,ic2)]
      ip(1:isize,2)=[cur_part(1:n2   ,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize,3)=[cur_part(n1:1:-1,ic1),cur_part(1:n2   ,ic2)]
      ip(1:isize,4)=[cur_part(1:n2   ,ic2),cur_part(n1:1:-1,ic1)]
      ip(1:isize,5)=[cur_part(1:n1   ,ic1),cur_part(n2:1:-1,ic2)]
      ip(1:isize,6)=[cur_part(n2:1:-1,ic2),cur_part(1:n1   ,ic1)]
      ip(1:isize,7)=[cur_part(n1:1:-1,ic1),cur_part(n2:1:-1,ic2)]
      ip(1:isize,8)=[cur_part(n2:1:-1,ic2),cur_part(n1:1:-1,ic1)]
      ! the binary label for the external particles included in this current
      cur_bin=this%cur_bin(ic1)+this%cur_bin(ic2)
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
            call add_one_to_currents(vertex_sign,cur_bin,ip(1:isize,i))
         endif
      enddo
    end subroutine add_all_3vert_to_currents
    logical function valid_current_order(ip)
      ! Checks that ip(1:isize) is an order for a current to be considered:
      ! the smallest number needs to come before the largest number in this
      ! list.
      implicit none
      integer :: i,maxi,mini,min_loc,max_loc
      integer,dimension(isize) :: ip
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
      else
         valid_current_order=.true.
      endif
    end function valid_current_order
    subroutine add_one_to_currents(vertex_sign,cur_bin,ip)
      ! The vertex contributions to the current 'ip'. Find this in the
      ! list of all currents and add this vertex to it.
      implicit none
      logical :: vertex_sign
      integer :: ic,cur_bin
      integer(kind=8) :: val
      integer,dimension(isize) :: ip
      if (this%vert_type(n_vert).ne.1) then
         ! gluon current
         call get_value(ip,0,val)
      else
         ! tensor current
         call get_value(ip,-1,val)
      endif
      call solve_dict(val,ic)
      if (this%cur_n_vert(ic).eq.0) then
         n_cur=n_cur+1
         cur_part(1:isize,ic)=ip(1:isize)
         this%cur_bin(ic)=cur_bin
         if (this%vert_type(n_vert).ne.1) then
            this%cur_type(ic)=0
         else
            this%cur_type(ic)=-1
         endif
      endif
      this%cur_n_vert(ic)=this%cur_n_vert(ic)+1
      this%cur_vertices(this%cur_n_vert(ic),ic)=n_vert
      this%cur_vert_sign(this%cur_n_vert(ic),ic)=vertex_sign
    end subroutine add_one_to_currents
    subroutine set_max_vert()
      ! computes the total number of needed vertices
      !
      ! For example, for next=6, we have for the 3-gluon vertices. Note that
      ! the denominators, i.e., symmetry factors due to inversion of the order
      ! in the currents, reduces the amount by a factor 2*2*2 (due to
      ! inversion of the combined current, and inversion of the two incoming
      ! currents separately, with the latter only applicable if the incoming
      ! currents contain more than 1 particle):
      ! - to compute the currents with 2 particles combined: (1+1)/2 ===> 5!/3! /2 = 10
      ! - to compute the currents with 3 particles combined: (1+2)/4+(2+1)/4 ===> 5!/2! /4 + 5!/2! /4 = 30
      ! - to compute the currents with 4 particles combined: (1+3)/4+(2+2)/8+(3+1)/4 ===> 5!/1! /4 + 5!/1! /8 + 5!/1! /4 = 75
      ! - to compute the currents with 5 particles combined: (1+4)/4+(2+3)/8+(3+2)/8+(4+1)/4 ===> 5!/0! + 5!/0! + 5!/0! = 90
      ! in total 205 3-gluon vertices.
      !
      ! To create the tensor particles, we need to double the amount of
      ! vertices we have to compute for th 2-4 particles combined currents:
      ! - 10+30+75 = 115 tensor creating currents.
      !
      ! To resolve the tensor particles we have more currents to choose from
      ! when combining to 3-5 particles. Note that a single particle currents
      ! cannot be a tensor currents and be careful not to include the
      ! contributions for which *both* incoming currents are tensor particles:
      ! - (1+2)/4+(2+1)/4 ===> 5!/2! /4 + 5!/2! /4 = 30
      ! - (1+3)/4+2*(2+2)/8+(3+1)/4 ===> 5!/1! /4 + 2* 5!/1! /8 + 5!/1! /4 = 90
      ! - (1+4)/4+2*(2+3)/8+2*(3+2)/8+(4+1)/4 ===> 5!/0! + 5!/0! + 5!/0! = 120
      ! resulting into 240 tensor resolving vertices.
      !
      ! Hence a total of 205+115+240=560 3-vertices need to be computed for 6
      ! gluon amplitudes.
      !
      ! If not decomposing the 4-gluon vertex into two three-vertices, a
      ! similar counting applies to the 4-gluon vertices, resulting in a total
      ! of 205+255=460 vertices to be computed.
      implicit none
      integer :: fact,fact2,iden,i,isplit,isplit2,itens
      real(kind=8) :: mv
      mv=0d0
      fact=factorial(next-1)
      do i=2,next-1
         fact2=fact/factorial(next-1-i)
         do isplit=1,i-1
            iden=2
            itens=1
            if (isplit.gt.1) iden=iden*2
            if (isplit.lt.i-1) iden=iden*2
            if (decompose_4vert) then
               if (i.ne.next-1) then
                  itens=itens+1
               endif
               if (i.ne.2) then
                  itens=itens+1
                  if (isplit.gt.1 .and. isplit.lt.i-1) itens=itens+1
               endif
            endif
            mv=mv+fact2/dble(iden)*itens
         enddo
      enddo
      max_vert=nint(mv)
      if (.not.decompose_4vert) then
         ! add the 4-gluon interactions
         mv=0d0
         do i=3,next-1
            fact2=fact/factorial(next-1-i)
            do isplit=1,i-2
               do isplit2=isplit+1,i-1
                  iden=2
                  if (isplit.gt.1) iden=iden*2
                  if (isplit2.gt.isplit+1) iden=iden*2
                  if (isplit2.lt.i-1) iden=iden*2
                  mv=mv+fact2/dble(iden)
               enddo
            enddo
         enddo
         max_vert=max_vert+nint(mv)
      endif
    end subroutine set_max_vert
    subroutine set_max_cur()
      ! Computes the total number of needed currents
      ! - Number of gluon currents:
      !   (next-1) + ( (next-1)*(next-2) + (next-1)*(next-2)*(next-3) + ... )/2
      ! - Number of tensor currents: 
      !   same as for the gluons except that the first and final terms are not present
      implicit none
      integer :: i,j,ifact
      max_cur=0
      do i=1,next-1
         ifact=next-1
         do j=1,i-1
            ifact=ifact*(next-1-j)
         enddo
         if (i.eq.1) then
            max_cur=max_cur+ifact
         elseif (i.lt.next-1) then
            if (decompose_4vert) then
               max_cur=max_cur+ifact
            else
               max_cur=max_cur+ifact/2
            endif
         else
            max_cur=max_cur+ifact/2
         endif
      enddo
    end subroutine set_max_cur
    subroutine create_current_dict()
      ! create an ordered dictionary that uniquely gives every current
      ! a label. This can be used to quickly find, (O(logN)), a
      ! current in the list of currents
      implicit none
      integer :: size,i,key
      integer(kind=8) :: val
      integer,dimension(:),allocatable :: ips_in,ips
      key=0
      size=1
      do isize=1,next-1
         size=size*(next-isize)
         allocate(ips_in(1:isize))
         do i=1,isize
            ips_in(i)=i
         enddo
         allocate(ips(1:isize))
         do i=1,size
            if (valid_current_order(ips_in)) then
               key=key+1
               call get_value(ips_in,0,val) ! add the gluon
               current_dict(key)=val
               if (isize.ne.1 .and. isize.ne.next-1 .and. decompose_4vert) then
                  key=key+1
                  call get_value(ips_in,-1,val) ! add the tensor
                  current_dict(key)=val
               endif
            endif
            call get_next_iperm(isize,ips_in,ips,next-1)
            ips_in=ips
         enddo
         deallocate(ips_in)
         deallocate(ips)
      enddo
    end subroutine create_current_dict
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
         val=val+int(ips(isize+1-j),kind=8)*int(next,kind=8)**int(j-1,kind=8)
      enddo
      ! Take the types into account (we have only 2 types (gluon and
      ! tensor), so multiply by two (and add one for the tensor))
      val=val*int(2,kind=8)
      if (itype.eq.-1) then
         val=val+int(1,kind=8)
      endif
    end subroutine get_value
    subroutine solve_dict(val,key)
      ! Given the value 'val', find the corresponding key in the
      ! 'current_dict' dictionary. Use a binary search
      ! algorithm. (This only works if the dictionary values are
      ! ordered, and all values only appear once).
      implicit none
      integer :: key,left,middle,right
      integer(kind=8) :: val
      left=1
      right=max_cur
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
      write (*,*) 'value not found in current dictionary',val
      stop 1
    end subroutine solve_dict
  end subroutine setup_imap_cache
  

  subroutine setup_imap(this,order)
    ! Initialises and fills the imap2 and imap3 arrays. Since we compute all
    ! permutations together, we must include --at all steps-- all possible
    ! ways to combine 2, 3, 4, ... external wavefunctions.  The loop over
    ! wavefunctions is as follows for 5 external particles (i.e., permutations
    ! over 4):
    !
    !   i      wavefunction                     imap2               imap3
    ! ----------------------------------------------------------------------------------
    !   1      A                                 -                    -
    !   2      B                                 -                    -
    !   3      C                                 -                    -
    !   4      D                                 -                    -
    !   5      AB=A+B                           1,2                   -
    !   6      AC=A+C                           1,2                   -
    !   7      AD=A+D                           1,2                   -
    !   8      BA=B+A                           2,3                   -
    !   9      BC=B+C                           2,3                   -
    !  10      BD=B+D                           2,3                   -
    !  11      CA=C+A                           3,4                   -
    !  12      CB=C+B                           3,4                   -
    !  13      CD=C+D                           3,4                   -
    !  14      DA=D+A                           3,4                   -
    !  15      DB=D+B                           3,4                   -
    !  16      DC=D+C                           3,4                   -
    !  17      ABC=A+BC + AB+C + A+B+C          1,9  / 5,3          1,2,3
    !  19      ABD=A+BD + AB+D + A+B+D          1,10 / 5,4          1,2,4
    !  20      ACB=A+CB + AC+B + A+C+B          1,12 / 6,2          1,3,2
    !  21      ACD=A+CD + AC+D + A+C+D          1,13 / 6,4          1,3,4
    !  22      ADB=A+DB + AD+B + A+D+B          1,15 / 7,2          1,4,2
    !  23      ADC=A+DC + AD+C + A+D+C          1,16 / 7,3          1,4,3
    !  24      BAC=B+AC + BA+C + B+A+C          1,6  / 8,3          2,1,3
    !  25      BAD=B+AD + BA+D + B+A+D          1,7  / 8,4          2,1,4
    !  26      BCA=B+CA + BC+A + B+C+A          1,11 / 9,1          2,3,1
    !  ...     ...                              ...                 ...
    !  39      DCA=D+CA + DC+A + D+C+A          4,11 / 16,1         4,3,1
    !  40      DCB=D+CB + DC+B + D+C+B          4,12 / 16,2         4,3,2
    !  41      ABCD=A+BCD + AB+CD + ABC+D +     1,27 / 5,13 / 17,4
    !               A+B+CD + A+BC+D + AB+C+D                        1,2,13 / 1,9,4 / 5,3,4
    !  42      ABDC=A+BDC + AB+DC + ABD+C +     1,28 / 5,16 / 19,3
    !               A+B+DC + A+BD+C + AB+D+C                        1,2,16 / 1,10,3 / 5,3,4
    !  ...     ...                              ...                 ...
    !  64      DCBA=D+CBA + DC+BA + DCB+A +     4,31 / 16,8 / 40,1
    !               D+C+BA + D+CB+A + DC+B+A                        4,3,8 / 4,12,1 / 16,2,1
    !
    !
    ! The 41st-64th wavefunctions can then be closed with the wavefunction for the
    ! 5th external particle (which will be at the i=0 spot) to form all the 24
    ! colour-ordered amplitudes.
    implicit none
    class(amplitude) :: this
    integer,dimension(this%next),optional :: order
    integer :: nnext
    integer,parameter :: zero=0
    integer :: i,j,isplit,s1,s2,jsplit,s3,ksplit,nb
    nnext=this%next-1
    call set_size(nnext,this%isize,this%nperm)
    if (.not. allocated(this%ncomb)) allocate(this%ncomb(0:this%isize)) 
    if (.not. allocated(this%nhel)) allocate(this%nhel(0:this%isize+1)) 
    if (.not. allocated(this%inverted)) allocate(this%inverted(0:this%isize))
    if (.not. allocated(this%nsplit2)) allocate(this%nsplit2(this%isize))
    if (.not. allocated(this%nsplit3)) allocate(this%nsplit3(this%isize))
    if (.not. allocated(this%istart)) allocate(this%istart(nnext))
    if (.not. allocated(this%ips)) allocate(this%ips(nnext,this%isize))
    if (.not. allocated(this%imap2)) allocate(this%imap2(2,nnext-1,this%isize))
    if (.not. allocated(this%imap3)) allocate(this%imap3(3,((nnext-1)*(nnext-2))/2,this%isize))
    ! number of ways the wavefunctions can be split in two:
    call set_nsplit2(nnext,this%isize,this%nperm,this%nsplit2)
    this%ncomb(0)=0
    if (this%sum_hel) then
       this%nhel(0)=2
    else
       this%nhel(0)=1
    endif
    do i=1,this%isize ! loop over wavefunction number
       ! number of particles that are combined:
       this%ncomb(i)=this%nsplit2(i)+1
       ! number of helicities/polarisations
       if (this%sum_hel) then
          this%nhel(i)=ibset(zero,this%ncomb(i))
       else
          this%nhel(i)=1
       endif
       ! number of ways it can be split in three:
       this%nsplit3(i)=max(0,this%nsplit2(i)*(this%nsplit2(i)-1)/2)
       if (this%ncomb(i).gt.this%ncomb(i-1)) then
          ! this%istart(#) tells where in the list of wavefunction we start combining
          ! # wavefunctions
          this%istart(this%ncomb(i))=i
          ! initialise this%ips(), will be updated by get_next_iperm() for increasing
          ! wavefunction number 'i'
          if (this%nperm.eq.1) then
             this%ips(1:this%ncomb(i),i)=order(1:this%ncomb(i))
          else
             do j=1,this%ncomb(i)
                this%ips(j,i)=j
             enddo
          endif
       else
          if (this%nperm.eq.1) then
             nb=i-this%istart(this%ncomb(i))
             this%ips(1:this%ncomb(i),i)=order(1+nb:this%ncomb(i)+nb)
          else
             call get_next_iperm(this%ncomb(i),this%ips(1,i-1),this%ips(1,i),nnext)
          endif
       endif
       ! Check if the corresponding inverted order already exits -- we can reuse it!
       this%inverted(i)=0
       if (this%nperm.gt.1) then
          do j=this%istart(this%ncomb(i)),i-1
             if(all(this%ips(this%ncomb(i):1:-1,i).eq.this%ips(1:this%ncomb(i),j))) then
                if (mod(this%ncomb(i),2).eq.0) then
                   this%inverted(i)=-j ! requires an additional minus sign in the wave-functions
                else
                   this%inverted(i)=j
                endif
                exit
             endif
          enddo
          if (this%inverted(i).ne.0)  cycle
       endif
       ksplit=0
       ! loop over the 2-splits (e.g. ABCD --> A+BCD , AB+CD , ABC+D):
       do isplit=1,this%nsplit2(i) 
          s1=isplit           ! length left of split
          s2=this%ncomb(i)-s1      ! length right of split
          call find_sub(this,i,0,s1,this%imap2(1,isplit,i))
          call find_sub(this,i,s1,s2,this%imap2(2,isplit,i))
          ! split the 2nd split once more to find the 3-splits
          ! (e.g. A+BCD --> A+B+CD , A+BC+D; AB+CD --> AB+C+D):
          do jsplit=1,this%nsplit2(i)-isplit
             ksplit=ksplit+1
             this%imap3(1,ksplit,i)=this%imap2(1,isplit,i) ! first split doesn't change
             s2=jsplit
             s3=this%ncomb(i)-(s1+s2)
             call find_sub(this,i,s1,s2,this%imap3(2,ksplit,i))
             call find_sub(this,i,s1+s2,s3,this%imap3(3,ksplit,i))
          enddo
       enddo
    enddo
    this%nhel(this%isize+1)=this%nhel(0)*this%nhel(this%isize)
  end subroutine setup_imap

  subroutine find_sub(this,ic,skip,s,isplits)
    ! By looping through the history, find which wavefunction number corresponds
    ! to the sub-list of permutations (for the current split).
    implicit none
    class(amplitude) :: this
    integer :: i,ic,s,isplits,b,e,skip
    b=skip+1
    e=skip+s
    do i=this%istart(s),this%istart(s+1)-1
       if (all(this%ips(1:s,i).eq.this%ips(b:e,ic))) then
          isplits=i
          exit
       endif
    enddo
  end subroutine find_sub


  subroutine set_size(ext,size,nperm)
    ! The total number of wavefunctions to compute:
    !    size = ext + ext*(ext-1) + ext*(ext-1)*(ext-2) + ...

    ! If nperm==1, then there is only a single permutation (instead of
    ! all), hence, the number of wavefunctions is much smaller:
    !    size = ext + ext-1 + ext-2 + ...
    implicit none
    integer :: ext,ifact,i,j,size,nperm
    if (nperm.eq.1) then
       size=0
       do i=1,ext
          size=size+(ext-i+1)
       enddo
    else
       size=0
       do i=1,ext
          ifact=ext
          do j=1,i-1
             ifact=ifact*(ext-j)
          enddo
          size=size+ifact
       enddo
    endif
  end subroutine set_size

  subroutine set_nsplit2(ext,isize,nperm,nsplit2)
    ! fills the array with the possible number of splits for
    ! wavefunction 'i'. This should be equal to the number of combined
    ! external wavefunctions minus one. However, since we don't know
    ! directly how many are combined at this stage, we have to infer
    ! it from the wavefunction number 'i'.
    implicit none
    integer :: ext,nperm,isize
    integer,dimension(isize) :: nsplit2
    integer :: i,j,size,ifact
    if (nperm.eq.1) then
       size=0
       do i=1,ext
          ifact=(ext-i+1)
          nsplit2(size+1:size+ifact)=i-1
          size=size+ifact
       enddo
    else
       size=0
       do i=1,ext
          ifact=ext
          do j=1,i-1
             ifact=ifact*(ext-j)
          enddo
          nsplit2(size+1:size+ifact)=i-1
          size=size+ifact
       enddo
    endif
  end subroutine set_nsplit2
    
  subroutine get_next_iperm(ip,ips_in,ips,n)
    ! Given a permutation ips_in, find the next one and return it through ips.
    ! For example for ip=3 (length of permutation list), n=4 (elements to be
    ! considered in the permutation) this gives
    !
    !    ips_in        ips
    !-------------------------
    !    1,2,3   -->   1,2,4
    !    1,2,4   -->   1,3,2
    !    1,3,2   -->   1,3,4
    !    1,3,4   -->   1,4,2
    !    1,4,2   -->   1,4,3
    !    1,4,3   -->   2,1,3
    !    2,1,3   -->   2,1,4
    !    2,1,4   -->   2,3,1
    !    2,3,1   -->   2,3,4
    !    2,3,4   -->   2,4,1
    !    2,4,1   -->   2,4,3
    !    2,4,3   -->   3,1,2
    !    3,1,2   -->   3,1,4
    !    3,1,4   -->   3,2,1
    !    3,2,1   -->   3,2,4
    !    3,2,4   -->   4,1,2
    !    4,1,2   -->   4,1,3
    !    4,1,3   -->   4,2,1
    !    4,2,1   -->   4,2,3
    !    4,2,3   -->   4,3,1
    !    4,3,1   -->   4,3,2
    !    4,3,2   -->   XXXXX
    !
    ! Note that when giving non-sensical inputs (e.g., the last one in the
    ! list above), the code either goes into an infinite loop, or returns some
    ! bogus result. There is no check on the consistency of the input.
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

     subroutine iperm_decode(iperm,n,ida)
       ! Given the number of elements 'n' and the permutation number
       ! 'iperm' (with 1<=iperm<=n!), returns the corresponding
       ! permutation 'ida(1:n)' in lexicographic order.
       implicit none
       integer,intent(in) :: iperm,n
       integer :: j,k
       integer,dimension(n),intent(out) :: ida
       logical,dimension(n) :: avail
       do j=1,n-1
          ida(j)=mod((iperm-1)/factorial(n-j),n-(j-1))+1
       enddo
       ida(n)=1
       avail(1:n)=.true.
       avail(ida(1))=.false.
       do j=2,n
          do k=1,n
             if (avail(k)) then
                if (ida(j).eq.1) then
                   ida(j)=k
                   avail(k)=.false.
                   exit
                endif
                ida(j)=ida(j)-1
             endif
          enddo
       enddo
     end subroutine iperm_decode

     subroutine iperm_encode(n,ida,iperm)
       ! Given the number of elements 'n' and the permutation
       ! 'ida(1:n)', returns the corresponding lexicographic
       ! permutation number 'iperm' (1<=iperm<=n!).
       implicit none
       integer,intent(in) :: n
       integer,dimension(n),intent(in) :: ida
       integer,intent(out) :: iperm
       integer :: j,k
       integer,dimension(n) :: idb
       idb=ida
       iperm=1
       do k=1,n-1
          do j=k+1,n
             if (idb(j).gt.idb(k)) idb(j)=idb(j)-1
          enddo
          iperm=iperm+(idb(k)-1)*factorial(n-k)
       enddo
     end subroutine iperm_encode

end module amplitude_mod
