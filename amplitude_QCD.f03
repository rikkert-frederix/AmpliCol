module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
  end type current
  type interaction
     integer :: type,singlet_move
     integer,dimension(2) :: currents
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
  end type interaction
  type amplitude_QCD
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: diff_col_vals,pp
     integer :: n_cur,n_vert,imode,nColOrd,n_qqbar,max_pp
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end,helmap,n_col_vals, &
          pp_bin_to_i,pp_i_to_bin
     integer,dimension(:,:),allocatable :: perm
     integer,dimension(:,:,:),allocatable :: row_index,col_index
   contains
     procedure :: init,evaluate,init_col2
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,part,order)
    use math_functions
    implicit none
    class(amplitude_QCD) :: this
    integer::n,imode
    integer,dimension(n)::part,order
    integer :: isize,nc,isplit,n1,n2,bin1,bin2,ic1,ic2,iv,i,max_cur,max_vert,max_key
    real(kind=4) :: tAfter,tBefore
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:),allocatable :: key_to_current

    if (imode.eq.1) then
       write (*,*) 'Initialising amplitude for:'
       write (*,*) '   - all polarisation/helicity configurations'
       write (*,*) '   - a single colour order'
    elseif (imode.eq.2) then
       write (*,*) 'Initialising amplitude for:'
       write (*,*) '   - a single polarisation/helicity configuration'
       write (*,*) '   - all colour orders'
    elseif (imode.eq.3) then
       write (*,*) 'Initialising amplitude for:'
       write (*,*) '   - a single polarisation/helicity configuration'
       write (*,*) '   - a single colour order'
    else
       write (*,*) 'ERROR unknown operation mode',imode
       stop 1
    endif
    this%imode=imode

    if (this%imode.eq.2) then
       call define_canonical_color_order()
    else
       this%nColOrd=1
    endif

    call check_input_consistency()

    call set_max_cur()
    call set_max_vert()

    if (this%imode.eq.2) then
       call cpu_time(tBefore)
       allocate(current_dict(max_cur)) 
       call create_current_dict()
       allocate(key_to_current(max_key)) 
       key_to_current(1:max_key)=0
       call cpu_time(tAfter)
       write (*,*) '   dictionary created ',tAfter-tBefore
       if (max_key.ne.max_cur) then
          write (*,*) 'Number of dictionary keys is expected to be identical to the maximum number of currents',&
               max_key,max_cur
          stop 1
       endif
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
             if (order(nc).le.2 .and. abs(part(order(nc))).le.6) then ! initial quark states
                this%current_list(this%n_cur)%type=anti_current(part(order(nc))) ! switch quark <--> anti-quark for initial states
             else
                this%current_list(this%n_cur)%type=part(order(nc))
             endif
             this%current_list(this%n_cur)%bin=ibset(0,order(nc)-1) ! give binary label
             if (this%imode.eq.1) then
                this%current_list(this%n_cur)%nhel=2 ! all possible helicities !!! MASSLESS ONLY
             elseif (this%imode.eq.2 .or. this%imode.eq.3) then
                this%current_list(this%n_cur)%nhel=1 ! only one helicity
             endif
             this%current_list(this%n_cur)%n_vert=0
          enddo
       else
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
       this%n_cur_end(isize)=this%n_cur
       if (isize.ge.2) this%n_vert_end(isize)=this%n_vert
    enddo

    call simple_consistency_checks()

    ! All done. But there could be currents that are not needed. Filter them out
    call filter_dead_trees()
    write (*,*) 'Total number of currents and vertices',this%n_cur,this%n_vert
    if (this%imode.eq.1) call create_helicity_map()
    if (this%imode.eq.2) call allocate_and_fill_colour_permutations()
    call setup_momentum_array()
  contains
    subroutine allocate_and_fill_colour_permutations()
      implicit none
      ! allocate and fill the colour orders
      if (this%n_qqbar.eq.0) then
         allocate(this%perm(1:n-1,1:this%nColOrd))
         do nc=1,this%nColOrd
            if ((.not.use_symmetry) .or. &
                 (use_symmetry .and. nc.le.this%nColOrd/2)) then
               this%perm(1:n-1,nc)=this%current_list(this%n_cur_start(n-1)-1+nc)%order(1:n-1)
            elseif (use_symmetry .and. nc.gt.this%nColOrd/2) then
               this%perm(1:n-1,nc)=this%current_list(this%n_cur_start(n-1)-1+nc-this%nColOrd/2)%order(n-1:1:-1)
            endif
         enddo
      elseif (this%n_qqbar.eq.1) then
         allocate(this%perm(1:n-2,1:this%nColOrd))
         do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
            this%perm(1:n-2,nc-this%n_cur_start(n-1)+1)=this%current_list(nc)%order(2:n-1)
         enddo
      endif
    end subroutine allocate_and_fill_colour_permutations
    subroutine setup_momentum_array()
      implicit none
      integer :: ic
      integer,dimension(1:maskr(n)) :: pp_i_to_bin
      allocate(this%pp_bin_to_i(1:maskr(n)))
      this%pp_bin_to_i(1:maskr(n))=0
      this%max_pp=0
      do ic=1,this%n_cur
         if (this%pp_bin_to_i(this%current_list(ic)%bin).eq.0) then
            this%max_pp=this%max_pp+1
            this%pp_bin_to_i(this%current_list(ic)%bin)=this%max_pp
            pp_i_to_bin(this%max_pp)=this%current_list(ic)%bin
         endif
      enddo
      allocate(this%pp(0:3,1:this%max_pp))
      allocate(this%pp_i_to_bin(this%max_pp))
      this%pp_i_to_bin(1:this%max_pp)=pp_i_to_bin(1:this%max_pp)
    end subroutine setup_momentum_array
    
    subroutine simple_consistency_checks()
      implicit none
      if (this%n_vert.gt.max_vert) then
         write (*,*) 'ERROR: too many interactions: max_vert not set correctly',max_vert,this%n_vert
         stop 1
      endif
      if (this%n_cur.gt.max_cur) then
         write (*,*) 'ERROR: too many currents: max_cur not set correctly',max_cur,this%n_cur
         stop 1
      endif
      if (((this%imode.eq.1 .or. this%imode.eq.3) .and. (this%nColOrd.ne.this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)) .or. &
           ((this%imode.eq.2) .and. (this%n_qqbar.ne.0) .and. (this%nColOrd.ne.this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)) .or. &
           ((this%imode.eq.2) .and. (this%n_qqbar.eq.0) .and. use_symmetry .and. &
           (this%nColOrd.ne.2*(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)))) &
           then
         write (*,*) 'The total number of colour orders to consider should be equal to the '// &
              'number of max-size currents (except for all-gluon and using symmetry)', &
              this%nColOrd,this%n_cur_start(n-1),this%n_cur_end(n-1),this%n_qqbar,use_symmetry
         stop 1
      endif
    end subroutine simple_consistency_checks
  
    subroutine define_canonical_color_order()
      ! canonical order: (q,glu,glu,glu,singlet,singlet,qbar,q,qbar)
      use math_functions
      implicit none
      integer :: iord,i
      integer :: nq,naq,nglu,nsing,iq,iaq,iglu,ising
      nq=0; naq=0 ; nglu=0 ; nsing=0
      do i=1,n
         if (part(i).eq.21) nglu=nglu+1
         if ((i.gt.2 .and. (part(i).ge.1 .and. part(i).le.6)) .or. &
             (i.le.2 .and. (part(i).le.-1 .and. part(i).ge.-6))) nq=nq+1
         if ((i.gt.2 .and. (part(i).le.-1 .and. part(i).ge.-6) ).or. &
             (i.le.2 .and. (part(i).ge.1 .and. part(i).le.6))) naq=naq+1
         if (part(i).eq.22) nsing=nsing+1
      enddo
      if (nq.ne.naq) then
         write (*,*) 'not the same number of quarks and anti-quarks',nq,naq
         stop 1
      endif
      if (nq.gt.2) then
         write (*,*) 'more than two quarks',nq
         stop 1
      endif
      if (nq+naq+nsing+nglu.ne.n) then
         write (*,*) 'particle types do not add up',nq,naq,nsing,nglu,':',n
         stop 1
      endif
      iq=0; iaq=0 ; iglu=0 ; ising=0
      do i=1,n
         if (part(i).eq.21) then
            iglu=iglu+1
            if (nq.ge.1) then
               order(iglu+1)=i
            else
               order(iglu)=i
            endif
         elseif((i.gt.2 .and. (part(i).ge.1 .and. part(i).le.6)) .or. &
              (i.le.2 .and. (part(i).le.-1 .and. part(i).ge.-6))) then
            iq=iq+1
            if (iq.eq.1) then
               order(1)=i
            else
               order(n-1)=i
            endif
         elseif ((i.gt.2 .and. (part(i).le.-1 .and. part(i).ge.-6) ).or. &
              (i.le.2 .and. (part(i).ge.1 .and. part(i).le.6))) then
            iaq=iaq+1
            if (iaq.eq.1) then
               order(n)=i
            else
               order(n-2)=i
            endif
         elseif (part(i).eq.22) then
            ising=ising+1
            if(nq.ne.0) then
               order(1+nglu+ising)=i
            else
               order(nglu+ising)=i
            endif
         endif
      enddo

      if (nq.eq.0) then
         this%nColOrd=factorial(nglu-1)
      elseif (nq.eq.1) then
         this%nColOrd=factorial(nglu)
      else
         write (*,*) 'Number of colour orders unknown',nq
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
            if (part(i).ne.21 .and. part(i).ne.22) then
               quark_flav(abs(part(i)))=quark_flav(abs(part(i)))-sign(1,part(i))
               if (part(i).lt.0) this%n_qqbar=this%n_qqbar+1
            endif
         else
            if (part(i).ne.21 .and. part(i).ne.22) then
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
      integer :: isize,j,ifact
      if (this%imode.eq.1 .or. this%imode.eq.3) then
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
      elseif(this%imode.eq.2) then
         if (this%n_qqbar.eq.0) then
            ! for increasing isize:
            ! - Number of gluon currents (remove the '/2' if use_symmetry=.false.):
            !   (next-1) + ( (next-1)*(next-2) + (next-1)*(next-2)*(next-3) + ... )/2
            ! - Number of tensor currents: 
            !   same as for the gluons except that the first and final terms are absent
            max_cur=0
            do isize=1,n-1
               ifact=n-1
               do j=1,isize-1
                  ifact=ifact*(n-1-j)
               enddo
               if (isize.eq.1) then
                  max_cur=max_cur+ifact
               elseif (isize.lt.n-1) then
                  if (use_symmetry) then
                     max_cur=max_cur+ifact
                  else 
                     max_cur=max_cur+ifact*2
                  endif
               else
                  if (use_symmetry) then
                     max_cur=max_cur+ifact/2
                  else
                     max_cur=max_cur+ifact
                  endif
               endif
            enddo
            max_cur=max_cur+1
         elseif(this%n_qqbar.eq.1) then
            ! for increasing isize:
            ! - Number of gluon currents (remove the '/2' if use_symmetry=.false.):
            !   (next-2) + ( (next-2)*(next-3) + (next-2)*(next-3)*(next-4) + ... )/2
            ! - Number of tensor currents: 
            !   same as for the gluons except that the first term is absent
            !   (final is absent as well, but we only know that after the dead
            !   tree-filtering)
            ! - Number of quark currents:
            !   1 + (next-1) + (next-1)*(next-2) + (next-1)*(next-2)*(next-3) + ...
            max_cur=0
            do isize=1,n-1
               ! gluons and tensors
               ifact=n-2
               do j=1,isize-1
                  ifact=ifact*(n-2-j)
               enddo
               if (isize.eq.1) then
                  max_cur=max_cur+ifact
               else
                  if (use_symmetry) then
                     max_cur=max_cur+ifact
                  else 
                     max_cur=max_cur+ifact*2
                  endif
               endif
               ! quarks
               ifact=1
               do j=1,isize-1
                  ifact=ifact*(n-1-j)
               enddo
               max_cur=max_cur+ifact
            enddo
            max_cur=max_cur+1
         endif
      endif
    end subroutine set_max_cur

    subroutine set_max_vert()
      ! rough upper bound on the maximum number of interactions
      implicit none
      integer :: isize,fact,iden,isplit,itens,fact2,next
      real(kind=8) :: mv
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         max_vert=0
         do isize=2,n-1
            if (isize.eq.2) then
               max_vert=max_vert+(n-isize)*3
            else
               max_vert=max_vert+isize*(n-isize)*3
            endif
         enddo
      elseif(this%imode.eq.2) then
         ! gluon and tensor vertices
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
         ! - to compute the currents with 5 particles combined: (1+4)/4+(2+3)/8+(3+2)/8+(4+1)/4 ===> 5!/0! /4 + 5!/0! /8 + 5!/0! /8 5!/0! /4 = 90
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
         if (this%n_qqbar.eq.0) then
            next=n-1
         elseif(this%n_qqbar.eq.1) then
            next=n-2
         endif
         mv=0d0
         fact=factorial(next)
         do isize=2,next
            fact2=fact/factorial(next-isize)
            do isplit=1,isize-1
               iden=2
               itens=1
               if (isplit.gt.1) iden=iden*2
               if (isplit.lt.isize-1) iden=iden*2
               if (isize.ne.n-1) then
                  itens=itens+1
               endif
               if (isize.ne.2) then
                  itens=itens+1
                  if (isplit.gt.1 .and. isplit.lt.isize-1) itens=itens+1
               endif
               if (.not.use_symmetry) iden=1
               mv=mv+fact2/dble(iden)*itens
            enddo
         enddo
         ! add the quark vertices
         if (this%n_qqbar.eq.1) then
            do isize=2,n-1
               do isplit=1,isize-1
                  iden=1
                  if (isplit.gt.1 .and. use_symmetry) iden=iden*2
                  mv=mv+fact/(factorial(n-isize-1))/dble(iden)
               enddo
            enddo
         endif
         max_vert=nint(mv)
      endif
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
      where_to_cur=0
      where_to_ver=0
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
      if (.not.valid_current_combination()) return
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
      elseif (this%current_list(ic1)%type.eq.22 .and. &
           (this%current_list(ic2)%type.ge.1 .and. this%current_list(ic2)%type.le.6)) then
         ! add a photon-quark to quark vertex
         call add_vertex(4,this%current_list(ic2)%type)
      elseif ((this%current_list(ic1)%type.ge.1 .and. this%current_list(ic1)%type.le.6) .and. &
           this%current_list(ic2)%type.eq.22) then
         ! add a quark-photon to quark vertex
         call add_vertex(6,this%current_list(ic1)%type)
      endif
    end subroutine add_if_allowed_threevertex


    logical function all_gluon_current(bin)
      ! returns .true. only if all external particles related to the binary
      ! label 'bin' are gluons
      implicit none
      integer :: bin,i,j
      all_gluon_current=.true.
      do i=1,n
         if (btest(bin,i-1) .and. part(i).ne.21) then
            all_gluon_current=.false.
            return
         endif
      enddo
    end function all_gluon_current
      
    logical function valid_current_combination()
      implicit none
      integer :: i,j,k
      logical :: gluon_current
      integer,dimension(isize) :: ip
      valid_current_combination=.false.
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(this%current_list(ic1)%bin,this%current_list(ic2)%bin)).ne.isize) return

      ! final particle should never be part of any combined currents: it will
      ! be used to close the amplitude instead
      if (n1.eq.1) then
         if (this%current_list(ic1)%order(n1).eq.order(n)) return
      endif
      if (n2.eq.1) then
         if (this%current_list(ic2)%order(n2).eq.order(n)) return
      endif

      ip(1:isize)=[this%current_list(ic1)%order(1:n1),this%current_list(ic2)%order(1:n2)]
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         ! check that current combination is compatible with the input colour
         ! order. Skip colour singlets
         do i=1,isize
            if (part(ip(i)).ne.22) exit
         enddo
         do j=1,n
            if (order(j).eq.ip(i)) exit
         enddo
         do k=1,isize-1
            if (i+k.gt.isize) exit            ! they are compatible
            if (part(ip(i+k)).eq.22) i=i+1
            if (part(order(j+k)).eq.22) j=j+1
            if (i+k.gt.isize) exit            ! they are compatible
            if (j+k.gt.n) return              ! incompatible: passed end of the order() array
            if (ip(i+k).ne.order(j+k)) return ! incompatible: order is different
         enddo
      endif

      ! If using symmetry and the current is a combination of all external
      ! gluons, take only one of the two possible orders
      gluon_current=all_gluon_current(this%current_list(ic1)%bin+this%current_list(ic2)%bin)
      if (use_symmetry .and. this%imode.eq.2 .and. gluon_current) then
         if (maxval(this%current_list(ic1)%order(1:n1)).ge.maxval(this%current_list(ic2)%order(1:n2))) return
      endif

!!$      write (*,*) any(this%current_list(ic1)%order(1:n1).eq.order(1)) &
!!$           ,this%current_list(ic1)%order(1).ne.order(1),any(this%current_list(ic2)%order(1:n2).eq.order(1)),&
!!$           gluon_current,this%current_list(ic1)%bin,this%current_list(ic2)%bin

      if (.not. gluon_current) then
         ! if quark is in there, it should be the very first particle
         if (any(this%current_list(ic1)%order(1:n1).eq.order(1)) .and. &
              this%current_list(ic1)%order(1).ne.order(1)) return
         if (any(this%current_list(ic2)%order(1:n2).eq.order(1))) then
            return
         endif
      endif
      valid_current_combination=.true.
    end function valid_current_combination
    
    subroutine add_vertex(itype,ctype)
      implicit none
      integer :: itype,ctype
      if (isize.eq.n-1 .and. ctype.ne.anti_current(this%current_list(n)%type)) return ! dead tree. Filter already here
      this%n_vert=this%n_vert+1
      this%interaction_list(this%n_vert)%type=itype
      this%interaction_list(this%n_vert)%currents(1)=ic1
      this%interaction_list(this%n_vert)%currents(2)=ic2
      this%interaction_list(this%n_vert)%singlet_move=0
      call add_all_currents(ctype)
    end subroutine add_vertex

    function combined_currents()
      implicit none
      integer,dimension(isize) :: combined_currents
      integer :: i
      do i=n1,1,-1
         if (part(this%current_list(ic1)%order(i)).eq.22) then
            this%interaction_list(this%n_vert)%singlet_move=this%interaction_list(this%n_vert)%singlet_move+1
         else
            exit
         endif
      enddo
      if (i.eq.n1) then
         combined_currents(1:isize)=[this%current_list(ic1)%order(1:n1),&
                                     this%current_list(ic2)%order(1:n2)]
      else
         combined_currents(1:isize)=[this%current_list(ic1)%order(1:i),&
                                     this%current_list(ic2)%order(1:n2),&
                                     this%current_list(ic1)%order(i+1:n1)]
      endif
    end function combined_currents

    subroutine add_all_currents(ctype)
      implicit none
      logical :: vertex_sign
      integer,dimension(isize,8) :: ip
      integer :: i,cur_bin,ctype
      if (.not.use_symmetry .or. this%imode.eq.1 .or. this%imode.eq.3) then
         cur_bin=this%current_list(ic1)%bin+this%current_list(ic2)%bin
         ip(1:isize,1)=combined_currents()
         call add_current(.false.,cur_bin,ip(1:isize,1),ctype)
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
      integer :: ctype,cur_bin,ic,key
      integer(kind=8) :: val
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         write (*,*) 'add_current',ctype,':',ip,':',this%current_list(ic1)%order(1:n1),':',this%current_list(ic2)%order(1:n2)
         ! Check if this interaction can be added to an existing current
         do ic=1,this%n_cur
            if (ctype.ne.this%current_list(ic)%type) cycle
            if (cur_bin.ne.this%current_list(ic)%bin) cycle
            if (any(this%current_list(ic)%order(1:isize).ne.ip(1:isize))) cycle
            this%current_list(ic)%n_vert=this%current_list(ic)%n_vert+1
            this%current_list(ic)%vertices(this%current_list(ic)%n_vert)=this%n_vert
            this%current_list(ic)%vertex_sign(this%current_list(ic)%n_vert)=vertex_sign
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
      elseif (this%imode.eq.2) then
         if (ctype.eq.21) then
            ! gluon current
            call get_value(ip,0,val)
         elseif (ctype.eq.-21) then
            ! tensor current
            call get_value(ip,-1,val)
         elseif (ctype.ge.1 .and. ctype.le.6) then
            ! quark current
            call get_value(ip,1,val)
         endif
         call solve_dict(val,key)
         ic=key_to_current(key)
         if (ic.eq.0) then
            ! initialise new current
            this%n_cur=this%n_cur+1
            key_to_current(key)=this%n_cur
            ic=this%n_cur
            allocate(this%current_list(ic)%order(isize))
            this%current_list(ic)%order(1:isize)=ip(1:isize)
            this%current_list(ic)%type=ctype
            this%current_list(ic)%bin=cur_bin
            this%current_list(ic)%nhel=this%current_list(ic1)%nhel*this%current_list(ic2)%nhel
            if (ctype.eq.21) then
               allocate(this%current_list(ic)%vertices(5*(isize-1)))
               allocate(this%current_list(ic)%vertex_sign(5*(isize-1)))
            elseif (ctype.eq.-21) then
               allocate(this%current_list(ic)%vertices(isize-1))
               allocate(this%current_list(ic)%vertex_sign(isize-1))
            else
               allocate(this%current_list(ic)%vertices(2*(isize-1)))
               allocate(this%current_list(ic)%vertex_sign(2*(isize-1)))
            endif
            this%current_list(ic)%n_vert=0
         endif
         ! add the vertex to the current
         this%current_list(ic)%n_vert=this%current_list(ic)%n_vert+1
         this%current_list(ic)%vertices(this%current_list(ic)%n_vert)=this%n_vert
         this%current_list(ic)%vertex_sign(this%current_list(ic)%n_vert)=vertex_sign
      endif
    end subroutine add_current
       
    subroutine create_current_dict()
      ! Create a dictionary that uniquely gives every current a label. This
      ! can be used to quickly find, (O(logN)), a current in the list of
      ! currents. Note that when we create the dictionary, we must make sure
      ! that the val's are created in ascending order, and that we add an
      ! element to the dictionary for all possible val's. Hence, better to
      ! create a larger dictionary than strictly needed.
      implicit none
      integer :: size,i,j,key
      integer(kind=8) :: val
      integer,dimension(:),allocatable :: ips_in,ips
      key=n  ! skip the external currents.
      size=n
      do isize=2,n-1
         size=size*(n-isize+1)
         allocate(ips_in(1:isize))
         allocate(ips(1:isize))
         ! simply add all possible permutations
         do i=1,size
            if (i.eq.1) then
               do j=1,isize
                  ips_in(j)=j
               enddo
            else
               call get_next_iperm(isize,ips_in,ips,n)
               ips_in=ips
            endif
            if ((.not.use_symmetry) .or. valid_current_order(ips_in)) then
               if (any(ips_in(1:isize).eq.order(n))) cycle
               if (this%n_qqbar.eq.0 .or. &
                    (this%n_qqbar.eq.1 .and. all(ips_in(1:isize).ne.order(1)))) then
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.1 .and. isize.ne.n-1) then
                     key=key+1
                     call get_value(ips_in,-1,val)
                     current_dict(key)=val
                  endif
               endif
               if (this%n_qqbar.eq.1 .and. ips_in(1).eq.order(1)) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif
            endif
         enddo
         deallocate(ips_in)
         deallocate(ips)
      enddo
      max_key=key
    end subroutine create_current_dict

    subroutine get_value(ips,itype,val)
      ! Give every current a unique value. This is based on the
      ! (external) particles that are part of the current as well as
      ! the current type.
      implicit none
      integer,dimension(isize) :: ips
      integer :: j,itype
      integer(kind=8) :: val
      if (isize.eq.1) then
         write (*,*) 'current_dict only setup for isize.ge.2',isize
         stop 1
      endif
      val=0
      ! Give a unique identifier based on the external
      ! particles. Simply convert the list to an integer with base
      ! equal to the number of external particles.
      do j=1,isize
         val=val+int(ips(isize+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
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
      integer :: key,left,middle,right
      integer(kind=8) :: val
      left=1
      right=max_key
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
    end subroutine solve_dict

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
  end subroutine init


  subroutine evaluate(this,n,p,hel)
    use FeynmanRules
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,hel
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih1,ih2,ih,ih_in,ip
    integer :: ifinal,ihm1
    logical :: ls,le

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
       if (this%imode.eq.1 .or. this%imode.eq.3) then
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

    call fill_momentum_array()
    
    do isize=1,n-1
       if (isize.eq.1) then
          ! fill the external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)   
             if (this%current_list(ic)%order(1).le.2) then
                ifinal=-1
             else
                ifinal=1
             endif

             do ih=1,this%current_list(ic)%nhel
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
                      call ext_gluon_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                           ih_in,ifinal,this%current_list(ic)%val_r(1:4,ih))
                   else
                      call ext_gluon_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                           ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                   endif
                elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
                   call ext_quark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                   call ext_antiquark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                elseif (this%current_list(ic)%type.eq.22) then
                   if (use_real_gluons) then
                      call ext_gluon_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                           ih_in,ifinal,this%current_list(ic)%val_r(1:4,ih))
                   else
                      call ext_gluon_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                           ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih))
                   endif
                else
                   write (*,*) 'External particle type unknown',ic,this%current_list(ic)%type,ih
                   stop 1
                endif
             enddo
          enddo
          cycle
       endif

       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          do ih2=1,this%current_list(this%interaction_list(iv)%currents(2))%nhel
             do ih1=1,this%current_list(this%interaction_list(iv)%currents(1))%nhel
                ih=(ih2-1)*this%current_list(this%interaction_list(iv)%currents(1))%nhel+ih1
                if (this%interaction_list(iv)%type.eq.0) then
                   if (use_real_gluons) then
                      call threeGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4,ih1),&
                           this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                           this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
                                                this%interaction_list(iv)%val_r(1:4,ih))
                   else
                      call threeGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                           this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                                      this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                           this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
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
                    ! MESSY CODE: IMPROVE
                    if (this%interaction_list(iv)%singlet_move.eq.1) then
                       ls=btest(ih-1,popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin)-1)
                       ihm1=ih-1
                       call mvbits(ihm1,popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin), &
                            popcnt(this%current_list(this%interaction_list(iv)%currents(2))%bin), &
                            ihm1,&
                            popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin)-1)
                       if (ls) then
                          ih=1+ibset(ihm1,popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin+&
                                                 this%current_list(this%interaction_list(iv)%currents(2))%bin)-1)
                       else
                          ih=1+ibclr(ihm1,popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin+&
                                                 this%current_list(this%interaction_list(iv)%currents(2))%bin)-1)
                       endif
                    elseif (this%interaction_list(iv)%singlet_move.gt.1) then
                       write (*,*) 'Cannot do more than one singlet move at once'
                       stop 1
                    endif
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
    subroutine fill_momentum_array
      implicit none
      integer :: ip,ibin,i
      do ip=1,this%max_pp
         ibin=this%pp_i_to_bin(ip)
         this%pp(0:3,ip)=0d0
         do i=1,2 ! treat incoming momenta as out-going
            if (btest(ibin,i-1)) this%pp(0:3,ip)=this%pp(0:3,ip)-p(0:3,i)
         enddo
         do i=3,n
            if (btest(ibin,i-1)) this%pp(0:3,ip)=this%pp(0:3,ip)+p(0:3,i)
         enddo
      enddo
    end subroutine fill_momentum_array
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
      elseif (this%imode.eq.3) then
         if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
            this%amps_r(1)=sum(this%current_list(this%n_cur)%val_r(1:4,1)*this%current_list(n)%val_r(1:4,2))
         else
            this%amps(1)=sum(this%current_list(this%n_cur)%val_c(1:4,1)*this%current_list(n)%val_c(1:4,1))
         endif
      endif
    end subroutine compute_amps_from_currents

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
         call GluonPropagator_real(this%current_list(ic)%val_r,this%current_list(ic)%nhel, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      else
         call GluonPropagator(this%current_list(ic)%val_c,this%current_list(ic)%nhel, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      endif
    end subroutine include_gluon_propagator

    subroutine include_quark_propagator()
      implicit none
      call QuarkPropagator(this%current_list(ic)%val_c,this%current_list(ic)%nhel, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
    end subroutine include_quark_propagator
  end subroutine evaluate

  subroutine init_col2(this,n,order,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer :: col_acc,n
    integer,dimension(n) :: order,iper,jper
    integer :: iperm,jperm,ival,iacc
    integer,dimension(1:3) :: n_vals
    integer,parameter :: max_vals=100
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(max_vals,1:3) :: diff_vals
    integer,dimension(:,:),allocatable :: ic,ir
    logical colour_flow
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
    do jperm=iperm,this%nColOrd
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
    allocate(this%col_index(1:this%nColOrd**2,1:maxval(n_vals(1:3)),1:3))
    allocate(this%row_index(0:this%nColOrd,1:maxval(n_vals(1:3)),1:3))
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
    do iperm=1,this%nColOrd
       if (this%n_qqbar.eq.0) then
          iper(1:n)=[this%perm(1:n-1,iperm),n]
       elseif (this%n_qqbar.eq.1) then
          iper(1:n)=[order(1),this%perm(1:n-2,iperm),order(n)]
       endif
       do jperm=iperm,this%nColOrd ! only include upper triangle (i.e., loop starts at iperm instead of 1)
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
!!$                  col_fac(2) = dble(3**(n-1) - (n-2) * 3**(n-3))
                  ! include the full expansion
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
                  col_fac(2)=dble(col_factor)
                  call Tr_deallocate
               else
                  call check_NLC_1qqbar(n,jper(2:n-1),iper(2:n-1),acc)
!!$                  col_fac(2)=dble(acc*(3)**(n-3))
                  ! include the full expansion
                  if (acc.ne.0) then
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
                     col_fac(2)=dble(col_factor)
                     call Tr_deallocate
                  endif
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
               col_fac(3)=dble(col_factor)
            endif
            call Tr_deallocate
         endif
      endif
    end subroutine compute_color_factor
  end subroutine init_col2
  
end module amplitude_QCD_mod
