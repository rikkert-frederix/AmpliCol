module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: use_symm_cm=.true.
  logical,parameter :: use_cm_dict=.true.
  type current
     integer :: type,bin,n_vert
     integer,dimension(:),allocatable :: vertices,order,spin
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8) :: mass,width
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
     integer,dimension(:),allocatable :: singlet_mv
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
  end type interaction
  type amplitude_QCD
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: pp
     real(kind=8),dimension(:,:,:),allocatable :: diff_col_vals
     integer :: n_cur,n_vert,imode,nColOrd,n_qqbar,max_pp,n_sing
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end, &
          pp_bin_to_i,pp_i_to_bin,quark_index,map_2qq_amps
     integer,dimension(:,:),allocatable :: n_col_vals,perm,col_index
     integer,dimension(:,:,:),allocatable :: spins,i_col_i
     integer,dimension(:,:,:,:),allocatable :: row_index
     logical :: same_flav
     logical,dimension(:),allocatable :: include_product
   contains
     procedure,public :: init,evaluate,init_col,filter_helicity
     procedure,private :: filter_dead_trees
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,part,spin,mass,width,order,it)
    use math_functions
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,imode,it
    integer,dimension(n)::part,order
    integer,dimension(0:3,n) :: spin
    real(kind=8),dimension(n) :: mass,width
    type(current),dimension(:),allocatable :: current_list_local
    type(interaction),dimension(:),allocatable :: interaction_list_local
    integer :: isize,nc,isplit,n1,n2,ic1,ic2,max_cur,max_vert,max_key,ispin
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

    call check_input_consistency()

    if (this%imode.eq.2) then
       call define_canonical_color_order()
    else
       this%nColOrd=1
    endif

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
    endif

    allocate(current_list_local(max_cur))
    allocate(interaction_list_local(max_vert))
    allocate(this%n_cur_start(n))
    allocate(this%n_cur_end(n))
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
             if (nc.eq.n) this%n_cur_start(n)=this%n_cur+1
             do ispin=1,spin(0,order(nc))
                call create_external_current(nc,spin(ispin,order(nc)))
             enddo
             if (nc.eq.n) this%n_cur_end(n)=this%n_cur
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
    call allocate_current_list_and_interaction_list()

    ! All done. But there could be currents that are not needed. Filter them out
    write (*,*) 'Total number of currents and vertices before filter',this%n_cur,this%n_vert
    call this%filter_dead_trees(n)
    write (*,*) 'Total number of currents and vertices',this%n_cur,this%n_vert
    if (this%imode.eq.1) call setup_spin_list()
    if (this%imode.eq.2) call allocate_and_fill_colour_permutations()
    call setup_momentum_array()

  contains

    subroutine setup_spin_list()
      implicit none
      integer :: ih1,ih2,ih,i
      allocate(this%spins(n,1,1:(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)))
      ! Note: this must be done in the same order as the amps() are computed in 'compute_amps_from_currents'
      do ih1=1,this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
         do ih2=1,this%n_cur_end(n)-this%n_cur_start(n)+1
            ih=(ih1-1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)+ih2
            do i=1,n
               if (i.lt.n) then
                  this%spins(order(i),1,ih)=this%current_list(this%n_cur_start(n-1)+ih1-1)%spin(i)
               elseif (i.eq.n) then
                  this%spins(order(i),1,ih)=this%current_list(this%n_cur_start(n  )+ih2-1)%spin(1)
               endif
            enddo
         enddo
      enddo
      allocate(this%include_product(1:(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)))
      this%include_product(:)=.true.
    end subroutine setup_spin_list
    
    subroutine create_external_current(nc,ispin)
      implicit none
      integer,intent(in) :: nc,ispin
      this%n_cur=this%n_cur+1
      allocate(current_list_local(this%n_cur)%order(isize))
      current_list_local(this%n_cur)%order(1)=order(nc)
      current_list_local(this%n_cur)%mass=mass(order(nc))
      current_list_local(this%n_cur)%width=width(order(nc))
      if (order(nc).le.2 .and. abs(part(order(nc))).le.6) then ! initial quark states
         current_list_local(this%n_cur)%type=anti_current(part(order(nc))) ! switch quark <--> anti-quark for initial states
      else
         current_list_local(this%n_cur)%type=part(order(nc))
      endif
      current_list_local(this%n_cur)%bin=ibset(0,order(nc)-1) ! give binary label
      allocate(current_list_local(this%n_cur)%spin(isize))
      current_list_local(this%n_cur)%spin=ispin
      current_list_local(this%n_cur)%n_vert=0
    end subroutine create_external_current
    
    subroutine allocate_and_fill_colour_permutations()
      implicit none
      integer :: ind,nc
      ! allocate and fill the colour orders in 'this%perm'. These are simply
      ! the orders of the elements in the 'this%current_list' (with size n-1)
      ! together with the final element). Exception: when there are colour
      ! singlets, they will not be part of the this%perm (while they are part
      ! of the elements in the this%current_list.
      allocate(this%perm(1:n-this%n_sing,1:this%nColOrd))
      if (this%n_cur_end(n).ne.this%n_cur_start(n)) then
         write (*,*) 'More than one element to close the current. Not possible when imode==2'
         write (*,*) this%n_cur_start
         write (*,*) this%n_cur_end
         stop 1
      elseif (this%current_list(this%n_cur_start(n))%type.eq.22) then
         write (*,*) 'Final current (that closes the amplitude) cannot be a colour singlet'
         write (*,*) this%current_list(this%n_cur_start(n))%type
         stop 1
      endif
      if (((.not.use_symmetry .or. this%n_qqbar.ne.0)  .and. this%nColOrd.ne.(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)) .or. &
           (use_symmetry .and. this%nColOrd.ne.2*(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1) .and. this%n_qqbar.eq.0)) then
         write (*,*) 'Number of expected colour orders not compatible with number of size n-1 currents'
         write (*,*) use_symmetry,this%n_qqbar,this%nColOrd
         write (*,*) this%n_cur_start
         write (*,*) this%n_cur_end
         stop 1
      endif

      if (this%n_qqbar.eq.0) then
         if (this%n_sing.ne.0) then
            write (*,*) 'For all-gluon processes, there should not be any colour singlets',this%n_sing
            stop 1
         endif
         do nc=1,this%nColOrd
            if ((.not.use_symmetry) .or. (use_symmetry .and. nc.le.this%nColOrd/2)) then
               this%perm(1:n,nc)=[this%current_list(this%n_cur_start(n-1)-1+nc)%order(1:n-1),&
                    this%current_list(this%n_cur_start(n))%order(1)]
            elseif (use_symmetry .and. nc.gt.this%nColOrd/2) then
               this%perm(1:n,nc)=[this%current_list(this%n_cur_start(n-1)-1+nc-this%nColOrd/2)%order(n-1:1:-1),&
                    this%current_list(this%n_cur_start(n))%order(1)]
            endif
         enddo
         return ! all gluons done: return
      endif

      ! The first particle in the order should be the quark. Double check that it is unique
      do nc=this%n_cur_start(1)+1,this%n_cur_end(1)
         if (this%current_list(nc)%order(1)  .eq. &
              this%current_list(this%n_cur_start(1))%order(1)) then
            write (*,*) 'First current is not unique. Not possible when imode==2'
            write (*,*) (this%current_list(ind)%order(1),ind=this%n_cur_start(1),this%n_cur_end(1))
            stop 1
         endif
      enddo
      ! First particle should not be a colour singlet
      if (this%current_list(this%n_cur_start(1))%type.eq.22) then
         write (*,*) 'First particle is a colour singlet. Not possible'
         stop 1
      endif
      
      if (this%n_qqbar.eq.1) then
         do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
            this%perm(1:n-this%n_sing,nc-this%n_cur_start(n-1)+1)=[this%current_list(this%n_cur_start(1))%order(1),&
                 this%current_list(nc)%order(2:n-1-this%n_sing),this%current_list(this%n_cur_start(n))%order(1)]
         enddo
      elseif (this%n_qqbar.eq.2) then
         call setup_map_2qq_amps()
         do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
            this%perm(1:n-this%n_sing,this%map_2qq_amps(nc-this%n_cur_start(n-1)+1)) = &
                 [this%current_list(nc)%order(1:n-1-this%n_sing),this%current_list(this%n_cur_start(n))%order(1)]
         enddo
      endif
    end subroutine allocate_and_fill_colour_permutations

    subroutine setup_map_2qq_amps
      use math_functions
      implicit none
      integer :: i,j,k,m,q,nc
      integer,dimension(1:n-4-this%n_sing) :: first,perm_out,perm_in,gluons
      integer,dimension(1:n-2-this%n_sing) :: ord
      integer,dimension(1:n-2-this%n_sing,this%nColOrd) :: buff
      
      allocate(this%map_2qq_amps(this%nColOrd))
      if (n-4-this%n_sing.gt.0) then

         ! first define 'buff', which will be the list of colour orders in canonical order
         k=1
         do i=1,n
            if (part(i).eq.21) then
               gluons(k)=i
               k=k+1
            endif
         enddo
         m = 1
         do i=0,n-4-this%n_sing ! loop over the number of gluons between the first quark and anti-quark in the order
            do j=1,n-4-this%n_sing
               first(j) = j
            enddo
            do k=1,factorial(n-4-this%n_sing) ! loop over all gluon permutations
               call get_next_iperm(n-4-this%n_sing,first,perm_out,n-4-this%n_sing)
               do q=1,n-4-this%n_sing
                  perm_in(q) = gluons(perm_out(q))
               enddo
               buff(1:i,m) = perm_in(1:i) ! gluons between the first quark and anti-quark
               buff(i+1:i+2,m)=0          ! the anti-quark and quark. 
               buff(i+3:n-2-this%n_sing,m) = perm_in(i+1:n-4-this%n_sing) ! gluons between the second quark and anti-quark
               first = perm_out
               m = m+1
            enddo
         enddo

         ! then, setup the map from the order of the colour orders as they
         ! appear in the amps to the canonical order as defined in 'buff'
         do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
            ord(1:n-2-this%n_sing) = this%current_list(nc)%order(2:n-1)
            do i=1,n-2
               if (is_quark(ord(i)).or.is_antiquark(ord(i))) then
                  ord(i)=0
               endif
            enddo
            do i=1,(n-3-this%n_sing)*factorial(n-4-this%n_sing)
               if (all(buff(1:n-2-this%n_sing,i).eq.ord(1:n-2-this%n_sing))) then
                  this%map_2qq_amps(nc-this%n_cur_start(n-1)+1) = i
                  exit
               endif
            enddo
            if (i.eq.(n-3-this%n_sing)*factorial(n-4-this%n_sing)+1) then
               write (*,*) 'check_2qq_order failed: order not found in list',i
               write (*,*) nc,':',this%current_list(nc)%order(1:n)
               stop 1
            endif
         enddo
      else
         this%map_2qq_amps(1)=1
      endif
    end subroutine setup_map_2qq_amps

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
    end subroutine simple_consistency_checks
  
    subroutine define_canonical_color_order()
      ! canonical order: (q,glu,glu,glu,singlet,singlet,qbar,q,qbar)
      use math_functions
      implicit none
      integer :: i,nq,naq,nglu,nsing,iq,iaq,iglu,ising

      nq=0; naq=0 ; nglu=0 ; nsing=0
      do i=1,n
         if (part(i).eq.21) nglu=nglu+1
         if (is_quark(i)) nq=nq+1
         if (is_antiquark(i)) naq=naq+1
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
      order= 0
      do i=1,n
         if (part(i).eq.21) then
            iglu=iglu+1
            if (nq.ge.1) then
               order(iglu+1)=i
            else
               order(iglu)=i
            endif
         elseif(is_quark(i)) then
            iq=iq+1
            if (iq.eq.1) then
               order(1)=i
            else
               order(n-1)=i
            endif
         elseif (is_antiquark(i)) then
            iaq=iaq+1
            !order(n)=i ! for one-qq
            if (it.eq.1)then 
               if (iaq.eq.1) then
                  order(n)=i
               else
                  order(n-2)=i
               endif
            elseif (it.eq.2) then
               if (iaq.eq.1) then
                  order(n-2)=i
               else
                  order(n)=i
               endif
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
      elseif (nq.eq.2) then
         this%nColOrd=factorial(nglu)*(nglu+1)
      else
         write (*,*) 'Number of colour orders unknown',nq
         stop 1
      endif
     end subroutine define_canonical_color_order

    subroutine check_input_consistency()
      implicit none
      integer :: i,j,k
! Setup the quark_index. Labels where the quarks and anti-quarks are in the
! process. Quarks are the odd entries (quark_index(1) and quark_index(3)),
! while the anti-quarks are the even entries. If quark flavours are different,
! quark_index(2) will be the anti-quark of quark_index(1) and quark_index(4)
! the anti-quark of quark_index(3).
      if (this%n_qqbar.ge.1) then
         allocate(this%quark_index(2*this%n_qqbar))
         this%quark_index=0 ! initialise all to zero
         k=1
         do i=1,n
            if(is_quark(i)) then
               ! found a quark
               this%quark_index(k)=i
               k=k+2
            endif
         enddo
         do i=1,n
            if (is_antiquark(i)) then
               ! found an anti-quark. Find the corresponding quark in the
               ! quark_index list.
               do k=1,2*this%n_qqbar-1,2
                  if (abs(part(this%quark_index(k))).eq.abs((part(i)))) then
                     ! if there are identical quarks, the 'k+1' label could
                     ! already have been filled. If that is the case, cycle to
                     ! the next label
                     if (this%quark_index(k+1).ne.0) cycle 
                     this%quark_index(k+1)=i
                     exit
                  endif
               enddo
            endif
         enddo
      endif
      
      if (this%n_qqbar.gt.3) then
         write (*,*) 'ERROR: code only working for 0, 1 or 2 qqbar pairs',this%n_qqbar
         write (*,*) part
         stop 1
      endif
      if (imode.eq.1.or.imode.eq.3) then
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
         if (.not.(is_quark(order(1)))) then
            write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
            write (*,*) order
            write (*,*) part
            stop 1
         endif
         if (.not.(is_antiquark(order(n)))) then
            write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
            write (*,*) order
            write (*,*) part
            stop 1
         endif
      endif

      if (this%n_qqbar.ge.2) then
         do i=2,n-1
            if (is_antiquark(order(i))) then
               ! next should be a quark
               if (.not.(is_quark(order(i+1)))) then
                  write (*,*) 'ERROR: in the colour order, after an initial state quark should come a final state quark'
                  write (*,*) order
                  write (*,*) part
                  stop 1
               endif
            endif
         enddo
         if (use_real_gluons) then
            write (*,*) 'Cannot use real gluons with two quark lines around'
            stop 1
         endif
      endif

      endif
    end subroutine check_input_consistency

    subroutine set_max_cur()
      ! rough upper bound for the maximum number of currents
      implicit none
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         max_cur=2*n*(n-1)*2**n
      else
         max_cur=factorial(n+1)*14
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
         if (this%n_sing.ge.1) max_vert=max_vert*this%n_sing
         max_vert=max_vert*2**(n-2) ! spins
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
         elseif (this%n_qqbar.eq.2) then
            next=n-4
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
         elseif (this%n_qqbar.eq.2) then
            do isize=2,n-1
               do isplit=1,isize-1
                  iden=1
                  if (isplit.gt.1 .and. use_symmetry) iden=iden*2
                  mv=mv+fact/(factorial(n-isize-1))/dble(iden)
               enddo
            enddo
            mv=100*mv ! TO CHANGE!
         endif
         max_vert=nint(mv)
      endif
      
    end subroutine set_max_vert


    subroutine allocate_current_list_and_interaction_list()
      ! allocate the minimum memory needed for the current_list and
      ! interaction_list to be able to perform the evaluate() procedure.
      implicit none
      integer :: isize,ic,iv
      allocate(this%current_list(1:this%n_cur))
      do isize=1,n-1
         do ic=this%n_cur_start(isize),this%n_cur_end(isize)
            this%current_list(ic)=current_list_local(ic)
            if(this%current_list(ic)%n_vert.gt.0) then
               deallocate(this%current_list(ic)%vertices)
               allocate(this%current_list(ic)%vertices(this%current_list(ic)%n_vert))
               this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)=&
                    current_list_local(ic)%vertices(1:current_list_local(ic)%n_vert)
               deallocate(this%current_list(ic)%vertex_sign)
               allocate(this%current_list(ic)%vertex_sign(this%current_list(ic)%n_vert))
               this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)=&
                    current_list_local(ic)%vertex_sign(1:current_list_local(ic)%n_vert)
            endif
         enddo
      enddo
      allocate(this%interaction_list(1:this%n_vert))
      do iv=1,this%n_vert
         this%interaction_list(iv)=interaction_list_local(iv)
         deallocate(this%interaction_list(iv)%singlet_mv)
         allocate(this%interaction_list(iv)%singlet_mv(0:interaction_list_local(iv)%singlet_mv(0)))
         this%interaction_list(iv)%singlet_mv(0:interaction_list_local(iv)%singlet_mv(0))=&
              interaction_list_local(iv)%singlet_mv(0:interaction_list_local(iv)%singlet_mv(0))
      enddo
    end subroutine allocate_current_list_and_interaction_list

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
      if (.not.valid_current_combination())  then
         return
      endif
      
      if (current_list_local(ic1)%type.eq.21 .and. current_list_local(ic2)%type.eq.21) then
         ! add the gluon-gluon to gluon vertex
         call add_vertex(0,21)
         ! add the gluon-gluon to tensor vertex
         call add_vertex(1,-21)

      elseif (current_list_local(ic1)%type.eq.-21 .and. current_list_local(ic2)%type.eq.21) then
         ! add a tensor-gluon to gluon vertex
         call add_vertex(2,21)

      elseif (current_list_local(ic1)%type.eq.21 .and. current_list_local(ic2)%type.eq.-21) then
         ! add a gluon-tensor to gluon vertex
         call add_vertex(3,21)

      elseif (current_list_local(ic1)%type.eq.21 .and. &
           (current_list_local(ic2)%type.ge.1 .and. current_list_local(ic2)%type.le.6)) then
         ! add a gluon-quark to quark vertex
         call add_vertex(4,current_list_local(ic2)%type)

      elseif (current_list_local(ic1)%type.eq.21 .and. &
           (current_list_local(ic2)%type.le.-1 .and. current_list_local(ic2)%type.ge.-6)) then
         ! add a gluon-antiquark to antiquark vertex
         call add_vertex(5,current_list_local(ic2)%type)

      elseif ((current_list_local(ic1)%type.ge.1 .and. current_list_local(ic1)%type.le.6) .and. &
           current_list_local(ic2)%type.eq.21) then
         ! add a quark-gluon to quark vertex
         call add_vertex(6,current_list_local(ic1)%type)

      elseif ((current_list_local(ic1)%type.le.-1 .and. current_list_local(ic1)%type.ge.-6) .and. &
           current_list_local(ic2)%type.eq.21) then
         ! add a antiquark-gluon to antiquark vertex
         call add_vertex(7,current_list_local(ic1)%type)

      elseif ((current_list_local(ic1)%type.ge.1 .and. current_list_local(ic1)%type.le.6) .and. &
           (current_list_local(ic2)%type.eq.anti_current(current_list_local(ic1)%type))) then
         ! add a quark-antiquark to gluon vertex
         call add_vertex(8,21)

      elseif ((current_list_local(ic1)%type.le.-1 .and. current_list_local(ic1)%type.ge.-6) .and. &
           (current_list_local(ic2)%type.eq.anti_current(current_list_local(ic1)%type))) then
         ! add a antiquark-quark to gluon vertex
         call add_vertex(9,21)

      elseif (current_list_local(ic1)%type.eq.22 .and. &
           (current_list_local(ic2)%type.ge.1 .and. current_list_local(ic2)%type.le.6)) then
         ! add a photon-quark to quark vertex
         call add_vertex(4,current_list_local(ic2)%type)

      elseif ((current_list_local(ic1)%type.ge.1 .and. current_list_local(ic1)%type.le.6) .and. &
           current_list_local(ic2)%type.eq.22) then
         ! add a quark-photon to quark vertex
         call add_vertex(6,current_list_local(ic1)%type)

      elseif (current_list_local(ic1)%type.eq.22 .and. &
           (current_list_local(ic2)%type.le.-1 .and. current_list_local(ic2)%type.ge.-6)) then
         ! add a photon-antiquark to quark vertex
         call add_vertex(5,current_list_local(ic2)%type)

      elseif ((current_list_local(ic1)%type.le.-1 .and. current_list_local(ic1)%type.ge.-6) .and. &
           current_list_local(ic2)%type.eq.22) then
         ! add a antiquark-photon to quark vertex
         call add_vertex(7,current_list_local(ic1)%type)
      endif
    end subroutine add_if_allowed_threevertex

    
    logical function valid_current_combination()
      ! Checks to see if the combination of currents ic1 and ic2 could be a
      ! valid combination. Checks to perform:
      ! 0. All particles must be different in the two currents & final
      !    particle should never be part of the combined currents (it will be
      !    used to close the currents)
      ! 1. For imode=1 or 3, we need to make sure that the colour order is
      !    compatible with the input colour order
      ! 2. For imode=2 and use_symmetry=.true. --> skip one of the two colour
      !    orders if all gluons in current
      ! 3. For colour singlets:
      !  --> if one (or two) of the two currents is (are) a singlet, only
      !      include one of the two orders
      !  --> order of singlets themselves should be ignored in this check: all
      !      must be included
      ! 5. If it is not a gluon current, and if the first particle (of the
      !    colour order) is in the current, this must come in the first
      !    place. Note that this is consistent with the 'it' parameter in
      !    define_canonical_order() for the two-quark-line case.
      ! 6. In general, the current must be of the format "q g..g qbar q g..g
      !    qbar" or any subset thereof.
      implicit none
      integer :: i,j,nc1,nc2
      logical :: gluon_current,colour_singlet1,colour_singlet2,found_quark,found_antiquark
      integer,dimension(isize) :: ip
      valid_current_combination=.false.
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(current_list_local(ic1)%bin,current_list_local(ic2)%bin)).ne.isize) return
      
      ! final particle should never be part of any combined currents: it will
      ! be used to close the amplitude instead
      if (n1.eq.1) then
         if (current_list_local(ic1)%order(n1).eq.order(n)) then
            return
          endif
      endif
      if (n2.eq.1) then
         if (current_list_local(ic2)%order(n2).eq.order(n)) then
            return
         endif
      endif

      colour_singlet1=all_singlet_current(current_list_local(ic1)%bin)
      colour_singlet2=all_singlet_current(current_list_local(ic2)%bin)
      ! If the first current is a singlet and the second is not, it is not a valid order
      if (colour_singlet1 .and. (.not.colour_singlet2)) return
      ! If both currents are colour singlets, only consider one of the two. 
      if (colour_singlet1 .and. colour_singlet2) then
         if (maxval(current_list_local(ic1)%order(1:n1)).ge.maxval(current_list_local(ic2)%order(1:n2))) then
            return
         else
            valid_current_combination=.true.
            return ! no need to check further: below is ony checks about the colours
         endif
      endif
      
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         ! check that current combination is compatible with the input colour
         ! order.
         do i=1,n1
            if (part(current_list_local(ic1)%order(i)).ge.22) exit
         enddo
         nc1=i-1
         do i=1,n2
            if (part(current_list_local(ic2)%order(i)).ge.22) exit
         enddo
         nc2=i-1
         ip(1:nc1+nc2)=[current_list_local(ic1)%order(1:nc1),current_list_local(ic2)%order(1:nc2)]
         do j=1,n
            if (order(j).eq.ip(1)) then
               do i=2,nc1+nc2
                  if (j-1+i.gt.n) return
                  if (order(j-1+i).ne.ip(i)) return
               enddo
               exit
            endif
         enddo
      endif

      ! If using symmetry and the current is a combination of all external
      ! gluons, take only one of the two possible orders
      gluon_current=all_gluon_current(current_list_local(ic1)%bin+current_list_local(ic2)%bin)
      if (use_symmetry .and. this%imode.eq.2 .and. gluon_current) then
         if (maxval(current_list_local(ic1)%order(1:n1)).ge.maxval(current_list_local(ic2)%order(1:n2))) return
      endif
      if (.not. gluon_current) then
         ! If the very first particle is in the current, it should be in the
         ! first position. Note that this must be a quark! This also means
         ! that for one *and* two-quark-line amplitudes, both the first and
         ! final particle in the orders are fixed. Hence, for the
         ! two-quark-line case only one of the two possibilities of the
         ! quarks/anti-quarks ordering is considered. This is consistent with
         ! the use of the 'it' variable in the 'define_canonical_color_order'.
         if (any(current_list_local(ic1)%order(1:n1).eq.order(1)) .and. &
              current_list_local(ic1)%order(1).ne.order(1)) then
            return
         endif
         if (any(current_list_local(ic2)%order(1:n2).eq.order(1))) then
            return
         endif

         ! for two quark lines. Should be of the form "q g..g qbar q g..g qbar" or any subset thereof. 
         ip(1:isize)=[current_list_local(ic1)%order(1:n1),current_list_local(ic2)%order(1:n2)]
         found_quark=.false.
         found_antiquark=.false.
         do i=1,isize
            if (is_quark(ip(i))) then
               ! found a quark.
               if (found_quark) then
                  ! no anti-quark between two quarks
                  return
               endif
               found_antiquark=.false.
               found_quark=.true.
            elseif (is_antiquark((ip(i)))) then
               ! found an anti-quark
               if (found_antiquark) then
                  ! no quark between two anti-quarks
                  return
               endif
               ! next one must be a quark:
               if (i.lt.isize) then
                  if (.not. (is_quark(ip(i+1)))) then
                     return
                  endif
               endif
               found_quark=.false.
               found_antiquark=.true.
            endif
         enddo
      endif

      valid_current_combination=.true.

    end function valid_current_combination
    
    subroutine add_vertex(itype,ctype)
      implicit none
      integer :: itype,ctype
      if (isize.eq.n-1 .and. ctype.ne.anti_current(current_list_local(this%n_cur_start(n))%type)) then 
        return ! dead tree. Filter already here
      endif
      this%n_vert=this%n_vert+1
      interaction_list_local(this%n_vert)%type=itype
      interaction_list_local(this%n_vert)%currents(1)=ic1
      interaction_list_local(this%n_vert)%currents(2)=ic2
      allocate(interaction_list_local(this%n_vert)%singlet_mv(0:isize))
      call add_all_currents(ctype)
    end subroutine add_vertex

    function combined_currents(n1,n2,ip1,ip2,singlet_mv)
      ! just concatenate the two colour orders, except if there is a colour
      ! singlet. Move the singlet to the end of the combined current order.
      implicit none
      integer,dimension(isize) :: combined_currents
      integer :: i,n1,n2,ipos,mv12,nc1,nc2,ns1,ns2
      integer,dimension(n1) :: ip1
      integer,dimension(n2) :: ip2
      integer,dimension(0:isize) :: singlet_mv
      
      do i=1,n1
         if (part(ip1(i)).ge.22) exit
      enddo
      nc1=i-1
      do i=1,n2
         if (part(ip2(i)).ge.22) exit
      enddo
      nc2=i-1

      combined_currents(1:nc1+nc2)=[ip1(1:nc1),ip2(1:nc2)]
      if (nc1.eq.n1) then
         ! No colour singlets or all colour singlets are in ip2
         singlet_mv(0)=0
         combined_currents(nc1+nc2+1:n1+n2)=ip2(nc2+1:n2)
         return
      elseif(nc2.eq.n2) then
         ! Some colour singlets in ip1, but no in ip2
         singlet_mv(0)=n1-nc1
         singlet_mv(1:singlet_mv(0))=nc1+1
         combined_currents(nc1+nc2+1:n1+n2)=ip1(nc1+1:n1)
         return
      else
         ! Some colour singlets in both ip1 and ip2
         singlet_mv(0)=0
         ns1=nc1+1
         ns2=nc2+1
         if (nc2.eq.0) then
            ! Special case: no coloured particles in ip2
            if (ip1(n1).lt.ip2(1)) then
               ! nothing to move
               combined_currents(1:n1+n2)=[ip1(1:n1),ip2(1:n2)]
               return
            endif
            do while (ip1(ns1).lt.ip2(1))
               ns1=ns1+1
            enddo
            combined_currents(nc1+1:ns1)=ip1(nc1+1:ns1)
         endif
         do while (ip2(ns2).lt.ip1(ns1))
            ns2=ns2+1
            if (ns2.gt.n2) exit
         enddo
         combined_currents(ns1+nc2:ns1+ns2-2)=ip2(nc2+1:ns2-1)
         do ipos=ns1+ns2-1,n1+n2
            if (ns1.gt.n1) then
               mv12=2
            elseif(ns2.gt.n2) then
               mv12=1
            elseif(ip1(ns1).lt.ip2(ns2)) then
               mv12=1
            else
               mv12=2
            endif
            singlet_mv(0)=singlet_mv(0)+1
            if (mv12.eq.1) then
               combined_currents(ipos)=ip1(ns1)
               singlet_mv(singlet_mv(0))=ns1 - (singlet_mv(0)-1)
               ns1=ns1+1
            elseif (mv12.eq.2) then
               combined_currents(ipos)=ip2(ns2)
               singlet_mv(singlet_mv(0))=n1+ns2+1 - (singlet_mv(0)-1)
               ns2=ns2+1
            endif
         enddo
      endif
    end function combined_currents

    function combined_spin(spin,singlet_mv)
      ! just concatenate the two spins, except if there is a colour
      ! singlet. Move the spin of the singlet to the end of the combined
      ! order.
      implicit none
      integer,dimension(isize) :: combined_spin,spin
      integer,dimension(0:isize) :: singlet_mv
      integer :: imv
      combined_spin(1:isize)=spin(1:isize)
      do imv=1,singlet_mv(0)
         combined_spin(1:isize)=[combined_spin(1:singlet_mv(imv)-1), &
                                 combined_spin(singlet_mv(imv)+1:isize-1),combined_spin(singlet_mv(imv))]
      enddo
    end function combined_spin
    
    subroutine add_all_currents(ctype)
      implicit none
      logical,dimension(8) :: vertex_sign
      integer,dimension(isize,8) :: ip,spin
      integer :: i,cur_bin,ctype,nperm
      integer,dimension(0:isize) :: singlet_mv
      if (.not.use_symmetry .or. this%imode.eq.1 .or. this%imode.eq.3) then
         cur_bin=current_list_local(ic1)%bin+current_list_local(ic2)%bin
         ip(1:isize,1)=combined_currents(n1,n2,current_list_local(ic1)%order(1:n1), &
              current_list_local(ic2)%order(1:n2),singlet_mv)
         spin(1:isize,1)=combined_spin([current_list_local(ic1)%spin(1:n1), &
              current_list_local(ic2)%spin(1:n2)],singlet_mv)
         interaction_list_local(this%n_vert)%singlet_mv(0:isize)=singlet_mv(0:isize)
         call add_current(.false.,cur_bin,ip(1:isize,1),ctype,spin(1:isize,1))
         return
      endif

      ! Need to consider all the possible permutations
      call check_all_permutations(nperm,ip,vertex_sign,spin)
      cur_bin=current_list_local(ic1)%bin+current_list_local(ic2)%bin
      do i=1,nperm
         call add_current(vertex_sign(i),cur_bin,ip(1:isize,i),ctype,spin(1:isize,i))
      enddo
    end subroutine add_all_currents

    subroutine check_all_permutations(nperm,ip,vertex_sign,spin)
      ! If a current only contains (external) gluons, we can use symmetry to
      ! relate them to eachother. This subroutine checks all permutations,
      ! and, if they give a valid current order, adds that current to the list
      ! that should be included.
      implicit none
      integer,intent(out) :: nperm
      integer,intent(out),dimension(isize,8) :: ip,spin
      logical,intent(out),dimension(8) :: vertex_sign
      integer,dimension(1:n1,2) :: ip1,spin1
      integer,dimension(1:n2,2) :: ip2,spin2
      logical :: ag1,ag2,iden
      integer,dimension(3) :: switch
      integer :: i,j,k
      integer,dimension(0:isize,8) :: singlet_mv
      switch(1:3)=1
      ag1=all_gluon_current(current_list_local(ic1)%bin)
      ag2=all_gluon_current(current_list_local(ic2)%bin)
      if (n1.ge.2 .and. ag1) switch(1)=2
      if (n2.ge.2 .and. ag2) switch(2)=2
      if (ag1 .and. ag2) switch(3)=2
      ip1(1:n1,1)=current_list_local(ic1)%order(1:n1)
      ip1(1:n1,2)=current_list_local(ic1)%order(n1:1:-1)
      ip2(1:n2,1)=current_list_local(ic2)%order(1:n2)
      ip2(1:n2,2)=current_list_local(ic2)%order(n2:1:-1)
      spin1(1:n1,1)=current_list_local(ic1)%spin(1:n1)
      spin1(1:n1,2)=current_list_local(ic1)%spin(n1:1:-1)
      spin2(1:n2,1)=current_list_local(ic2)%spin(1:n2)
      spin2(1:n2,2)=current_list_local(ic2)%spin(n2:1:-1)
      nperm=0
      do i=1,switch(1)
         do j=1,switch(2)
            do k=1,switch(3)
               nperm=nperm+1
               if (k.eq.1) then
                  ip(1:isize,nperm)=combined_currents(n1,n2,ip1(1:n1,i),ip2(1:n2,j),&
                       singlet_mv(0,nperm))
                  spin(1:isize,nperm)=combined_spin([spin1(1:n1,i),spin2(1:n2,j)],&
                       singlet_mv(0,nperm))
               else
                  ip(1:isize,nperm)=combined_currents(n2,n1,ip2(1:n2,j),ip1(1:n1,i),&
                       singlet_mv(0,nperm))
                  spin(1:isize,nperm)=combined_spin([spin2(1:n2,j),spin1(1:n1,i)],&
                       singlet_mv(0,nperm))
               endif
               vertex_sign(nperm)=(k.eq.2 .xor. (j.eq.2 .and. mod(n2,2).eq.0) .xor. (i.eq.2 .and. mod(n1,2).eq.0))
               if (.not.valid_current_order_excl_symmetry(ip(1:isize,nperm))) nperm=nperm-1
            enddo
         enddo
      enddo

      if (nperm.eq.0) then
         write (*,*) 'No valid permutations. This should not happen'
         write (*,*) switch,ag1,ag2
         write (*,*) current_list_local(ic1)%order(1:n1),quark_in_current(current_list_local(ic1)%order(1:n1),n1)
         write (*,*) current_list_local(ic2)%order(1:n2),quark_in_current(current_list_local(ic2)%order(1:n2),n2)
         stop 1
      endif
      
      iden=.true.
      if (all(singlet_mv(0,1:nperm).eq.singlet_mv(0,1))) then
         do i=2,nperm
            if (any(singlet_mv(1:singlet_mv(0,1),i).ne.singlet_mv(1:singlet_mv(0,1),1))) then
               iden=.false.
               exit
            endif
         enddo
      endif
      if (iden) then
         interaction_list_local(this%n_vert)%singlet_mv(0:singlet_mv(0,1))=singlet_mv(0:singlet_mv(0,1),1)
      else
         write (*,*) 'Singlet move not identical for all permutations',nperm
         do i=1,nperm
            write (*,*) nperm,':',singlet_mv(0,i),':',singlet_mv(1:singlet_mv(0,i),i)
         enddo
         stop 1
      endif
    end subroutine check_all_permutations
    
    logical function valid_current_order_excl_symmetry(ip)
      ! Checks that ip(1:isize) is an order for a current to be considered
      ! when use_symmetry=.true. --> the smallest number needs to come before
      ! the largest number in this list. 
      implicit none
      integer :: i,maxi,mini,min_loc,max_loc
      integer,dimension(isize) :: ip

      ! if there is a quark (or anti-quark) in the current, no symmetry can be
      ! used. Hence, this is a valid order
      if (popcnt(quark_in_current(ip,isize)).ge.1) then
         valid_current_order_excl_symmetry=.true.
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
         valid_current_order_excl_symmetry=.false.
         return
      endif
      valid_current_order_excl_symmetry=.true.
    end function valid_current_order_excl_symmetry


    subroutine add_current(vertex_sign,cur_bin,ip,ctype,spin)
      implicit none
      logical :: vertex_sign
      integer,dimension(isize) :: ip,spin   ! permutation and spin of the current
      integer :: ctype,cur_bin,ic,key
      integer(kind=8) :: val
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         ! Check if this interaction can be added to an existing current
         do ic=1,this%n_cur
            if (ctype.ne.current_list_local(ic)%type) cycle
            if (cur_bin.ne.current_list_local(ic)%bin) cycle
            if (any(current_list_local(ic)%order(1:isize).ne.ip(1:isize))) cycle
            if (any(current_list_local(ic)%spin(1:isize).ne.spin(1:isize))) cycle
            current_list_local(ic)%n_vert=current_list_local(ic)%n_vert+1
            current_list_local(ic)%vertices(current_list_local(ic)%n_vert)=this%n_vert
            current_list_local(ic)%vertex_sign(current_list_local(ic)%n_vert)=vertex_sign
            return
         enddo
         ! Need a new current
         this%n_cur=this%n_cur+1
         allocate(current_list_local(this%n_cur)%order(isize))
         current_list_local(this%n_cur)%order(1:isize)=ip(1:isize)
         current_list_local(this%n_cur)%type=ctype
         current_list_local(this%n_cur)%bin=cur_bin
         allocate(current_list_local(this%n_cur)%spin(isize))
         current_list_local(this%n_cur)%spin(1:isize)=spin(1:isize)
         if (current_list_local(ic1)%mass.eq.current_list_local(ic2)%mass)  then
            current_list_local(this%n_cur)%mass=0d0
         else
            current_list_local(this%n_cur)%mass=max(current_list_local(ic1)%mass,current_list_local(ic2)%mass)
         endif
         if (current_list_local(ic1)%width.eq.current_list_local(ic2)%width)  then
            current_list_local(this%n_cur)%width=0d0
         else
            current_list_local(this%n_cur)%width=max(current_list_local(ic1)%width,current_list_local(ic2)%width)
         endif
         if (ctype.eq.21) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         elseif (ctype.eq.-21) then
            allocate(current_list_local(this%n_cur)%vertices(isize-1))
            allocate(current_list_local(this%n_cur)%vertex_sign(isize-1))
         else
            allocate(current_list_local(this%n_cur)%vertices(2*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(2*(isize-1)))
         endif
         current_list_local(this%n_cur)%vertices(1)=this%n_vert
         current_list_local(this%n_cur)%vertex_sign(1)=vertex_sign
         current_list_local(this%n_cur)%n_vert=1
      elseif (this%imode.eq.2) then
         if (ctype.eq.21) then
            ! gluon current
            call get_value(ip,0,val)
         elseif (ctype.eq.-21) then
            ! tensor current
            call get_value(ip,-1,val)
         elseif (ctype.ge.1 .and. ctype.le.6) then
            ! quark current
            call get_value(ip,2*ctype-1,val)
         elseif (ctype.ge.-6 .and. ctype.le.-1) then
            ! anti-quark current
            call get_value(ip,2*abs(ctype),val)
         endif

         call solve_dict(val,key)
         ic=key_to_current(key)
         if (ic.eq.0) then
            ! initialise new current
            this%n_cur=this%n_cur+1
            key_to_current(key)=this%n_cur
            ic=this%n_cur
            allocate(current_list_local(ic)%order(isize))
            current_list_local(ic)%order(1:isize)=ip(1:isize)
            current_list_local(ic)%type=ctype
            current_list_local(ic)%bin=cur_bin
            allocate(current_list_local(ic)%spin(isize))
            current_list_local(ic)%spin(1:isize)=spin(1:isize)
            if (any(spin(1:isize).ne.-9)) then
               write (*,*) 'trying to combine currents with different spin: not possible',spin(1:isize)
               stop 1
            endif
            if (ctype.eq.21) then
               allocate(current_list_local(ic)%vertices(5*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(5*(isize-1)))
            elseif (ctype.eq.-21) then
               allocate(current_list_local(ic)%vertices(isize-1))
               allocate(current_list_local(ic)%vertex_sign(isize-1))
            else
               allocate(current_list_local(ic)%vertices(2*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(2*(isize-1)))
            endif
            current_list_local(ic)%n_vert=0
         endif
         ! add the vertex to the current
         current_list_local(ic)%n_vert=current_list_local(ic)%n_vert+1
         current_list_local(ic)%vertices(current_list_local(ic)%n_vert)=this%n_vert
         current_list_local(ic)%vertex_sign(current_list_local(ic)%n_vert)=vertex_sign
      endif
    end subroutine add_current
       
    subroutine create_current_dict()
      ! Create a dictionary that uniquely gives every current a label. This
      ! can be used to quickly find, (O(logN)), a current in the list of
      ! currents. Note that when we create the dictionary, we must make sure
      ! that the val's are created in ascending order, and that we add an
      ! element to the dictionary for all possible val's. Hence, better to
      ! create a larger dictionary than strictly needed.
      use math_functions
      implicit none
      integer :: size,i,j,key
      integer(kind=8) :: val,previous_val
      integer,dimension(:),allocatable :: ips_in,ips
      key=n  ! skip the external currents.
      size=n
      previous_val=0
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
            ! add the gluon:
            key=key+1
            call get_value(ips_in,0,val)
            if (val.le.previous_val) then
               write (*,*) 'inconsistent current dictionary #1',val,previous_val
               stop 1
            endif
            current_dict(key)=val
            ! add the tensor
            key=key+1
            call get_value(ips_in,-1,val)
            if (val.le.previous_val) then
               write (*,*) 'inconsistent current dictionary #2',val,previous_val
               stop 1
            endif
            current_dict(key)=val
            ! add the quarks and anti-quarks
            do j=1,6
               key=key+1
               call get_value(ips_in,2*j-1,val) ! quarks are the odd ones
               if (val.le.previous_val) then
                  write (*,*) 'inconsistent current dictionary #3',val,previous_val
                  stop 1
               endif
               current_dict(key)=val
               key=key+1
               call get_value(ips_in,2*j,val)   ! anti-quarks are the even ones
               if (val.le.previous_val) then
                  write (*,*) 'inconsistent current dictionary #4',val,previous_val
                  stop 1
               endif
               current_dict(key)=val
            enddo
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
      ! Take the types into account (we have only 14 types (gluon,
      ! tensor and 6 quarks and 6 anti-quarks):
      val=val*int(14,kind=8) ! gluon
      if (itype.eq.-1) then
         val=val+int(1,kind=8) ! tensor
      elseif (itype.ge.1) then
         val=val+int(itype+1,kind=8) ! quark or anti-quark
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

    logical function all_gluon_current(bin)
      ! returns .true. only if all external particles related to the binary
      ! label 'bin' are gluons
      implicit none
      integer :: bin,i
      all_gluon_current=.true.
      do i=1,n
         if (btest(bin,i-1) .and. part(i).ne.21) then
            all_gluon_current=.false.
            return
         endif
      enddo
    end function all_gluon_current
    logical function all_singlet_current(bin)
      ! returns .true. only if all external particles related to the binary
      ! label 'bin' are colour singlets
      implicit none
      integer :: bin,i
      all_singlet_current=.true.
      do i=1,n
         if (btest(bin,i-1) .and. abs(part(i)).lt.22) then
            all_singlet_current=.false.
            return
         endif
      enddo
    end function all_singlet_current
    logical function is_quark(io)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io
      if ( (io.le.2  .and. (part(io).le.-1 .and. part(io).ge.-6)) .or. &
           (io.gt.2  .and. (part(io).ge. 1 .and. part(io).le. 6))) then
         is_quark=.true.
      else
         is_quark=.false.
      endif
    end function is_quark
    logical function is_antiquark(io)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io
      if ( (io.le.2  .and. (part(io).ge. 1 .and. part(io).le. 6)) .or. &
           (io.gt.2  .and. (part(io).le.-1 .and. part(io).ge.-6))) then
         is_antiquark=.true.
      else
         is_antiquark=.false.
      endif
    end function is_antiquark
    integer function quark_in_current(ip,isize)
      ! binary function that checks which 'quark_index' is an external
      ! particle part of the current. (i.e., It sets the first bit if
      ! quark_index(1) is part of the current, the second bit if
      ! quark_index(2) is part of the current, etc.)
      implicit none
      integer :: isize,i
      integer,dimension(isize) :: ip
      quark_in_current=0
      do i=1,2*this%n_qqbar
         if (any(this%quark_index(i).eq.ip(1:isize))) quark_in_current=ibset(quark_in_current,i-1)
      enddo
    end function quark_in_current
  end subroutine init


  subroutine evaluate(this,n,p,mass,width,hel,part)
    use FeynmanRules
    implicit none
    class(amplitude_QCD) :: this
    integer :: n
    integer,dimension(n)::part,hel
    real(kind=8),dimension(n) :: mass,width
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih_in,ip,ifinal
    if (.not. allocated(this%current_list(1)%val_c) .and. .not.allocated(this%current_list(1)%val_r)) then
       do ic=1,this%n_cur
          if (this%current_list(ic)%type.eq.-21) then
             if (use_real_gluons) then
                allocate(this%current_list(ic)%val_r(1:6))
             else
                allocate(this%current_list(ic)%val_c(1:6))
             endif
          elseif (this%current_list(ic)%type.eq.21 .and. use_real_gluons) then
             allocate(this%current_list(ic)%val_r(1:4))
          else
             allocate(this%current_list(ic)%val_c(1:4))
          endif
       enddo
       do iv=1,this%n_vert
          if (this%interaction_list(iv)%type.eq.1) then
             if (use_real_gluons) then
                allocate(this%interaction_list(iv)%val_r(1:6))
             else
                allocate(this%interaction_list(iv)%val_c(1:6))
             endif
          elseif ((this%interaction_list(iv)%type.eq.0 .or. &
                   this%interaction_list(iv)%type.eq.2 .or. &
                   this%interaction_list(iv)%type.eq.3) .and. use_real_gluons  ) then
             allocate(this%interaction_list(iv)%val_r(1:4))
          else
             allocate(this%interaction_list(iv)%val_c(1:4))
          endif
       enddo
       if (this%imode.eq.1 .or. this%imode.eq.3) then
          if (use_real_gluons .and. this%n_qqbar.eq.0) then
             allocate(this%amps_r(1:(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)))
          else
             allocate(this%amps(1:(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)))
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

             if (this%current_list(ic)%spin(1).eq.-9) then
                ih_in=max(0,hel(this%current_list(ic)%order(1)))
             else
                ih_in=max(0,this%current_list(ic)%spin(1))
             endif
             if (this%current_list(ic)%type.eq.21) then
                if (use_real_gluons) then
                   call ext_gluon_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_r(1:4))
                else
                   call ext_gluon_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4))
                endif
             elseif (this%current_list(ic)%type.ge.1 .and. this%current_list(ic)%type.le.6 ) then
                call ext_quark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                call ext_antiquark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (this%current_list(ic)%type.eq.22) then
                if (use_real_gluons) then
                   call ext_gluon_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_r(1:4))
                else
                   call ext_gluon_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4))
                endif
             else
                write (*,*) 'External particle type unknown',ic,this%current_list(ic)%type,ih_in
                stop 1
             endif
          enddo
          cycle
       endif

       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          if (this%interaction_list(iv)%type.eq.0) then
             if (use_real_gluons) then
                call threeGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                     this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
                     this%interaction_list(iv)%val_r(1:4))
             else
                call threeGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                     this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                     this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
                     this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.1) then
             if (use_real_gluons) then
                call TwoGluonToTensor_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     this%interaction_list(iv)%val_r(1:6))
             else
                call TwoGluonToTensor(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                     this%interaction_list(iv)%val_c(1:6))
             endif

          elseif(this%interaction_list(iv)%type.eq.2) then
             if (use_real_gluons) then
                call TensorGluontoGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:6),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     this%interaction_list(iv)%val_r(1:4))
             else
                call TensorGluontoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:6),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                     this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.3) then
             if (use_real_gluons) then
                call GluonTensortoGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                             this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:6),&
                                             this%interaction_list(iv)%val_r(1:4))
             else
                call GluonTensortoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:6),&
                                        this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.4) then
             if (use_real_gluons) then
                call GluonQuarktoQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                            this%interaction_list(iv)%val_c(1:4))
             else
                call GluonQuarktoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                       this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.5) then
             if (use_real_gluons) then
                call GluonAquarktoAquark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                              this%interaction_list(iv)%val_c(1:4))
             else
                call GluonAquarktoAquark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:4))
             endif
          elseif(this%interaction_list(iv)%type.eq.6) then
             if (use_real_gluons) then
                call QuarkGluontoQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                            this%interaction_list(iv)%val_c(1:4))
             else
                call QuarkGluontoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                       this%interaction_list(iv)%val_c(1:4))
             endif
          elseif(this%interaction_list(iv)%type.eq.7) then
             if (use_real_gluons) then
                call AquarkGluontoAquark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                              this%interaction_list(iv)%val_c(1:4))
             else
                call AquarkGluontoAquark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:4))
             endif
                 
          elseif(this%interaction_list(iv)%type.eq.8) then
             call QuarkAquarktoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                     this%interaction_list(iv)%val_c(1:4))

          elseif(this%interaction_list(iv)%type.eq.9) then
             call AquarkQuarktoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                     this%interaction_list(iv)%val_c(1:4))

          else
             write (*,*) 'Unknown vertex type: not yet implemented',iv,this%interaction_list(iv)%type
             stop 1
          endif
       enddo

       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          if (this%current_list(ic)%type.eq.21) then
             call combine_interactions(4)
             ! a gluon current
             if (isize.ne.n-1)  then
                call include_gluon_propagator()
             endif
          elseif ((this%current_list(ic)%type.ge.1.and.this%current_list(ic)%type.le.6)) then
             ! a quark current
             call combine_interactions(4)
             if (isize.ne.n-1)  then
                call include_quark_propagator()
             endif
          elseif (this%current_list(ic)%type.eq.-21) then
             ! the non-propagating tensor current
             call combine_interactions(6)
          elseif ((this%current_list(ic)%type.le.-1.and.this%current_list(ic)%type.ge.-6)) then
             ! an anti-quark current
             call combine_interactions(4)
             if (isize.ne.n-1)  then
                call include_aquark_propagator()
             endif
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
      integer :: ih1,ih2,ih,ic,ihc
      if (this%imode.eq.1) then
         ! Note: this must be done in the same order as the this%spins() are setup in 'setup_spin_list()'
         ihc=0
         do ih1=1,this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
            do ih2=1,this%n_cur_end(n)-this%n_cur_start(n)+1
               ih=(ih1-1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)+ih2
               if (.not.this%include_product(ih)) cycle
               ihc=ihc+1
               if (use_real_gluons .and. this%current_list(ih2)%type.eq.21) then
                  this%amps_r(ihc)=sum(this%current_list(this%n_cur_start(n-1)+ih1-1)%val_r(1:4)*&
                                       this%current_list(this%n_cur_start(n  )+ih2-1)%val_r(1:4))
               else
                  this%amps(ihc)=sum(this%current_list(this%n_cur_start(n-1)+ih1-1)%val_c(1:4)*&
                                     this%current_list(this%n_cur_start(n  )+ih2-1)%val_c(1:4))
               endif
            enddo
         enddo

      elseif (this%imode.eq.2) then
         do ic=this%n_cur_start(n-1),this%n_cur_end(n-1)
            if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
               this%amps_r(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_r(1:4)*this%current_list(n)%val_r(1:4))
            else
               if (this%n_qqbar.eq.2) then
                  this%amps(this%map_2qq_amps(ic-this%n_cur_start(n-1)+1)) = &
                       sum(this%current_list(ic)%val_c(1:4)*this%current_list(n)%val_c(1:4))
               else        
                  this%amps(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_c(1:4)*this%current_list(n)%val_c(1:4))
               endif     
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
            this%amps_r(1)=sum(this%current_list(this%n_cur)%val_r(1:4)*this%current_list(n)%val_r(1:4))
         else
            this%amps(1)=sum(this%current_list(this%n_cur)%val_c(1:4)*this%current_list(n)%val_c(1:4))
         endif
      endif
    end subroutine compute_amps_from_currents

    subroutine combine_interactions(dim)
      implicit none
      integer :: dim,iv
      if (use_real_gluons .and. abs(this%current_list(ic)%type).eq.21) then
         this%current_list(ic)%val_r(1:dim)=0d0
         do iv=1,this%current_list(ic)%n_vert
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_r(1:dim)=&
                    this%current_list(ic)%val_r(1:dim)-this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim)
            else
               this%current_list(ic)%val_r(1:dim)=&
                    this%current_list(ic)%val_r(1:dim)+this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim)
            endif
         enddo

      else
         this%current_list(ic)%val_c(1:dim)=(0d0,0d0)
         do iv=1,this%current_list(ic)%n_vert
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_c(1:dim)=&
                    this%current_list(ic)%val_c(1:dim)-this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim)
            else
               this%current_list(ic)%val_c(1:dim)=&
                    this%current_list(ic)%val_c(1:dim)+this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim)
            endif
         enddo
      endif
    end subroutine combine_interactions
    subroutine include_gluon_propagator()
      implicit none
      if (use_real_gluons) then
         call GluonPropagator_real(this%current_list(ic)%val_r, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      else
         call GluonPropagator(this%current_list(ic)%val_c, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      endif
    end subroutine include_gluon_propagator

    subroutine include_quark_propagator()
      implicit none
      call QuarkPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), & 
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_quark_propagator

    subroutine include_aquark_propagator()
      implicit none
      call AquarkPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_aquark_propagator
  end subroutine evaluate

  subroutine init_col(this,n,part,it,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer,parameter :: max_vals=3000
    integer :: col_acc,n,iperm,jperm,ival,iacc,isum,i,j,gi,ui,uj,uj_upper,it,iperm_upper,&
         gi_iperm,key,max_keys,jperm_lower
    integer,dimension(n) :: iper,jper,part
    integer,dimension(:,:),allocatable :: n_vals,first_rows
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(:,:,:),allocatable :: diff_vals
    real(kind=8),dimension(:,:,:),allocatable :: col_vals
    integer,dimension(:,:,:),allocatable :: ic,ir,n_colour_elements
    integer(kind=8),dimension(:),allocatable :: perm_dict

    write (*,*) 'Initialising colour matrix ...'

    call create_perm_dict()

    if (this%n_qqbar.eq.0) then
       iperm_upper = 1 ! dummy, needed for 2qq
       gi_iperm = 1 ! dummy
       uj_upper = 1 ! dummy
       ui = 1 ! dummy
    elseif (this%n_qqbar.eq.1) then
       iperm_upper = 1 ! dummy, needed for 2qq
       gi_iperm = 1 ! dummy
       uj_upper = 1 ! dummy
       ui = 1 ! dummy
    elseif (this%n_qqbar.eq.2) then
       iperm_upper = (n-4)+1 ! number of gluon separations on two quark lines
       uj_upper = 2
       ui = it
    endif
    allocate(n_vals(1:3,iperm_upper))
    allocate(diff_vals(max_vals,1:3,iperm_upper))
    allocate(this%i_col_i(max_vals,1:3,iperm_upper))
    allocate(n_colour_elements(max_vals,1:3,iperm_upper))
    allocate(col_vals(1:3,max_keys,iperm_upper))
    allocate(first_rows(1:n,iperm_upper))

! first check a single row in the colour matrix to determine how many
! different colour factors there are. In the case of two quark lines, we need
! to consider multiple rows corresponding to the number of gluons between the
! first quark and anti-quark in the colour order.
    n_vals(1:3,:)=0
    col_vals(1:3,1:max_keys,:)=0d0
    do iperm=1,iperm_upper
       if (this%n_qqbar.eq.0 .or. this%n_qqbar.eq.1) then
          iper(1:n-this%n_sing)=this%perm(1:n-this%n_sing,iperm)
       elseif (this%n_qqbar.eq.2) then
          ! Loop through all orders and find the first 'iper' with 'iperm-1'
          ! gluons between the first quark and anti-quark in the colour order.
          do j=1,this%nColOrd
             iper(1:n-this%n_sing)=this%perm(1:n-this%n_sing,j)
             do i=2,n-2
                if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
                   gi = i - 2 ! number of gluons between the first quark and anti-quark in the colour order
                   exit
                endif
             enddo
             if (gi.eq.iperm-1) exit
          enddo
          gi_iperm = iperm
       endif
       if (use_cm_dict) first_rows(1:n-this%n_sing,gi_iperm) = iper(1:n-this%n_sing)
       
       do jperm=1,this%nColOrd 
          do uj=1,uj_upper
             jper(1:n-this%n_sing)=this%perm(1:n-this%n_sing,jperm)
             if (this%n_qqbar.eq.2 .and. uj.ne.ui) call get_other_quark_order(jper)
             key=solve_dict(get_value(jper(1:n)))
             call compute_color_factor(col_acc,n-this%n_sing,iper,jper,ui,uj,col_fac)
             if (use_symm_cm.and.this%n_qqbar.ne.2) then
                col_fac(1:3)=col_fac(1:3)*2d0 ! include a factor 2 for the off-diagonal terms
                if (iperm.eq.jperm.and.ui.eq.uj) col_fac(1:3)=col_fac(1:3)*0.5d0 
             endif
             do iacc=1,3
                if (col_fac(iacc).eq.0d0) cycle
                col_vals(iacc,key,iperm)=col_fac(iacc)
                do ival=1,n_vals(iacc,gi_iperm)
                   if (col_fac(iacc).eq.diff_vals(ival,iacc,gi_iperm)) then
                      n_colour_elements(ival,iacc,gi_iperm)=n_colour_elements(ival,iacc,gi_iperm)+1
                      exit
                   endif
                enddo
                if (ival.ge.max_vals) then
                   write (*,*) 'Too many different colour factors. Increase max_vals',&
                        ival,n_vals(1:3,gi_iperm),max_vals
                   stop 1
                elseif (ival.eq.n_vals(iacc,gi_iperm)+1) then
                   ! new colour factor
                   n_vals(iacc,gi_iperm)=ival
                   diff_vals(ival,iacc,gi_iperm)=col_fac(iacc)
                   n_colour_elements(ival,iacc,gi_iperm)=1
                endif
             enddo
          enddo
       enddo
       
       write (*,*) 'A single row in the colour matrix has',n_vals(1:3,gi_iperm),&
            ' different colour factors at LC, NLC and full colour, respectively'
    enddo

    ! determine i_col_i:
    isum=1
    do iacc=1,3
       do gi_iperm=1,iperm_upper
          do ival=1,n_vals(iacc,gi_iperm)
             this%i_col_i(ival,iacc,gi_iperm)=isum
             isum=isum+n_colour_elements(ival,iacc,gi_iperm)*this%nColOrd
          enddo
       enddo
    enddo

 ! Allocate the arrays now that we know their sizes
    allocate(ic(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    allocate(ir(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    allocate(this%col_index(1:isum,iperm_upper))
    allocate(this%row_index(0:this%nColOrd,1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper)) 
    this%row_index(0,1:maxval(n_vals(1:3,iperm_upper)),1:3,1:iperm_upper)=0
    this%col_index(1,1:iperm_upper)=0
    allocate(this%n_col_vals(1:3,iperm_upper))
    this%n_col_vals(1:3,1:iperm_upper)=n_vals(1:3,1:iperm_upper)
    allocate(this%diff_col_vals(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    do iacc=1,3
       this%diff_col_vals(1:n_vals(iacc,iperm_upper),iacc,1:iperm_upper)=diff_vals(1:n_vals(iacc,iperm_upper),iacc,1:iperm_upper)
    enddo

! Compute all the colour factors and fill the col_index and row_index arrays
    ic=0
    ir=0

    do iperm=1,this%nColOrd
       iper(1:n-this%n_sing)=this%perm(1:n-this%n_sing,iperm)
       if (this%n_qqbar.eq.2) then
          ! find out what channel it belongs to, find gi
          do i=2,n-2
             if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
                gi = i - 2
                exit
             endif
          enddo
          gi_iperm = gi + 1
       endif

       jperm_lower=1
       if (use_symm_cm.and.this%n_qqbar.ne.2) jperm_lower = iperm

       do jperm=jperm_lower,this%nColOrd
          do uj=1,uj_upper
             jper(1:n-this%n_sing)=this%perm(1:n-this%n_sing,jperm)
             if (this%n_qqbar.eq.2 .and. uj.ne.ui) call get_other_quark_order(jper)
             if (use_cm_dict) then
                ! GET color factors from permuting first row
                call get_col_fac(iper,jper,ui,uj,gi_iperm,col_fac)
             else
                ! COMPUTE color factors again
                call compute_color_factor(col_acc,n-this%n_sing,iper,jper,ui,uj,col_fac)
                if (use_symm_cm.and.this%n_qqbar.ne.2) then
                   col_fac(1:3)=col_fac(1:3)*2d0
                   if (iperm.eq.jperm.and.ui.eq.uj) col_fac(1:3)=col_fac(1:3)*0.5d0 ! include a factor 2 for the off-diagonal terms
                endif
             endif
             do iacc=1,3
                if (col_fac(iacc).eq.0d0) cycle
                do ival=1,n_vals(iacc,gi_iperm)
                   if (col_fac(iacc).eq.diff_vals(ival,iacc,gi_iperm)) exit
                enddo
                ic(ival,iacc,gi_iperm)=ic(ival,iacc,gi_iperm)+1
                ir(ival,iacc,gi_iperm)=ir(ival,iacc,gi_iperm)+1
                this%col_index(this%i_col_i(ival,iacc,gi_iperm)+ic(ival,iacc,gi_iperm),gi_iperm)=(uj-1)*this%nColOrd+jperm
             enddo
          enddo
       enddo
       do iacc=1,3
          this%row_index(iperm,1:n_vals(iacc,gi_iperm),iacc,gi_iperm)=ir(1:n_vals(iacc,gi_iperm),iacc,gi_iperm)
       enddo
    enddo

    write (*,*) '... colour matrix initialised'
  contains
   subroutine get_col_fac(iper,jper,ui,uj,gi_iperm,col_fac)
     implicit none
     integer,intent(in) :: gi_iperm,ui,uj
     integer,dimension(n),intent(in) :: iper,jper
     integer,dimension(n) :: col_new,row_first,row_per,col_per
     integer :: i,j,key
     real(kind=8),dimension(1:3),intent(out) :: col_fac
     
     ! First row
     row_first(1:n-this%n_sing)=first_rows(1:n-this%n_sing,gi_iperm)
     if (this%n_qqbar.eq.2 .and. uj.ne.ui) call get_other_quark_order(row_first)

     ! Row in consideration
     row_per(1:n-this%n_sing)=iper(1:n-this%n_sing)
     if (this%n_qqbar.eq.2 .and. uj.ne.ui) call get_other_quark_order(row_per)

     ! Column in consideration
     col_per(1:n-this%n_sing)=jper(1:n-this%n_sing)

     do i=1,n-this%n_sing
        do j=1,n-this%n_sing
           if (col_per(i) .eq. row_per(j)) exit
        enddo
        if (.not.(abs(part(col_per(i))).le.6.and.abs(part(col_per(i))).ge.1)) then
          col_new(i) = row_first(j)
        else 
          col_new(i) = col_per(i)
        endif
     enddo

     key=solve_dict(get_value(col_new(1:n)))
     
     col_fac(1:3)=col_vals(1:3,key,gi_iperm)
   end subroutine get_col_fac

   subroutine get_other_quark_order(jper)
     implicit none
     integer,dimension(n) :: jper,temp_part,jper_new
     integer :: i,aq1,aq2
     logical first
     temp_part=part
     do i=1,n
        if (abs(part(i)).ge.1.and.abs(part(i)).le.6) then
           if (i.le.2) temp_part(i)=-part(i)
        endif
     enddo
     first=.true.
     do i=1,n
        if (temp_part(jper(i)).le.-1..and.temp_part(jper(i)).ge.-6.and.first) then
           aq1=i
           first=.false.
        endif
        if (temp_part(jper(i)).le.-1.and.temp_part(jper(i)).ge.-6.and..not.first) then
           aq2=i
        endif
     enddo
     jper_new = jper
     jper_new(aq1)=jper(aq2)
     jper_new(aq2)=jper(aq1)
     jper=jper_new
   end subroutine get_other_quark_order

   integer function solve_dict(val)
     ! Given the value 'val', find the corresponding key in the 'perm_dict'
     ! dictionary. Use a binary search algorithm. (This only works if the
     ! dictionary values are ordered, and all values only appear once).
     implicit none
     integer :: left,middle,right
     integer(kind=8) :: val
     left=1
     right=max_keys
     do while (left.le.right)
        middle=(right+left)/2
        if (perm_dict(middle).eq.val) then
           solve_dict=middle
           return
        elseif(perm_dict(middle).gt.val) then
           right=middle-1
        else
           left=middle+1
        endif
     enddo
   end function solve_dict
    
    subroutine create_perm_dict()
      ! Create a dictionary that uniquely gives every colour permutation a
      ! label. This can be used to quickly find, (O(logN)), a permutation in
      ! the list of permutations. Note that when we create the dictionary, we
      ! must make sure that the val's are created in ascending order, and that
      ! we add an element to the dictionary for all possible val's. Hence,
      ! better to create a larger dictionary than strictly needed.
      use math_functions
      implicit none
      integer :: iperm,i
      integer(kind=8) :: val,previous_val
      integer,dimension(:),allocatable :: iper,iper_in
      if (this%n_sing.ne.0) then
         write (*,*) 'fix create_perm_dict when there are color singlets'
         stop 1
      endif
      allocate(iper(1:n))
      allocate(iper_in(1:n))
      max_keys=factorial(n)
      allocate(perm_dict(1:max_keys))
      do i=1,n
         iper(i)=i
      enddo
      previous_val=0
      do iperm=1,max_keys
         val=get_value(iper(1:n))
         if (val.le.previous_val) then
            write (*,*) 'In create_perm_dict need to get values in ascending order',val,previous_val
            stop 1
         else
            previous_val=val
         endif
         perm_dict(iperm)=val
         iper_in=iper
         call get_next_iperm(n,iper_in,iper,n)
      enddo
      deallocate(iper)
    end subroutine create_perm_dict

    integer(kind=8) function get_value(iper)
      ! Give a unique identifier based on the colour order. Simply convert the
      ! list to an integer with base equal to the number of elements in the
      ! order.
      implicit none
      integer :: j
      integer,dimension(1:n) :: iper
      get_value=0
      do j=1,n
         get_value=get_value+int(iper(n+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
      enddo
    end function get_value
    
    subroutine compute_color_factor(col_acc,n,iper,jper,ui,uj,col_fac)
      use color_algebra
      implicit none
      integer :: i,n,acc,col_acc,color_fac,k,ui,uj,gi,gj
      real(kind=8),dimension(1:3) :: col_fac
      integer,dimension(n) :: iper,jper
      integer,dimension(n-4) :: iper_glu,jper_glu,iper_ord,jper_ord
      real(kind=16) :: col_factor
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
         elseif (this%n_qqbar.eq.2) then
            if (all(iper.eq.jper)) then
               if (ui.eq.uj.and.ui.eq.1) then
                  col_fac(1)=dble(3**(n-2))
               elseif (ui.eq.uj.and.ui.eq.2.and..not.this%same_flav) then
                  col_fac(1)=dble(3**(n-4))
               elseif (ui.eq.uj.and.ui.eq.2.and.this%same_flav) then
                  col_fac(1)=dble(3**(n-2))
               endif
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
               ! include the full expansion
               call Tr_allocate(n)
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=1 ! that term is single string of matrices
               Tr(0,1,1)=2*(n-2)
               Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
               Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
               coef(1)=1d0
               coef_Nc(:,:)=0
               coef_Nc(0,1)=1
               call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
               col_fac(2)=dble(col_factor)
               call Tr_deallocate
            else
               call check_NLC_1qqbar(n,jper(2:n-1),iper(2:n-1),acc)
               col_fac(2)=dble(acc*(3)**(n-3))
               ! include the full expansion
               if (acc.ne.0) then
                  call Tr_allocate(n)
                  Tr(0,0,0)=1 ! one term
                  Tr(0,0,1)=1 ! that term is single string of matrices
                  Tr(0,1,1)=2*(n-2)
                  Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
                  Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
                  coef(1)=1d0
                  coef_Nc(:,:)=0
                  coef_Nc(0,1)=1
                  call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
                  col_fac(2)=dble(col_factor)
                  call Tr_deallocate
               endif
            endif
         elseif (this%n_qqbar.eq.2) then
            do i=1,n-1
               if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
                  if (i.ne.1) then
                     gi = i - 2
                     exit
                  endif
               endif
            enddo
            do i=1,n-1
               if ((abs(part(jper(i))).ge.1.and.abs(part(jper(i))).le.6)) then
                  if (i.ne.1) then
                     gj = i - 2
                     exit
                  endif
               endif
            enddo
            k=1
            do i=1,n
               if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) cycle
               iper_glu(k)=iper(i)
               k=k+1
            enddo
            k=1
            do i=1,n
               if ((abs(part(jper(i))).ge.1.and.abs(part(jper(i))).le.6)) cycle
               jper_glu(k)=jper(i)
               k=k+1
            enddo
            call Tr_allocate(n)
            call convert_gluon_string(n,iper_glu,jper_glu,iper_ord,jper_ord)
            call check_NLC_2qqbar(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
            if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
            if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-4)) ! NLC parts
            if (this%same_flav) then
               call check_NLC_2qqbar_SF(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
               if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
               if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-3)) ! NLC parts
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
            coef(1)=1d0
            coef_Nc(:,:)=0
            coef_Nc(0,1)=1
            ! compute the colour factor by simplifying the colour string
            call Tr_full_simplify(col_factor) 
            col_fac(3)=0d0
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
            coef(1)=1d0
            coef_Nc(:,:)=0
            coef_Nc(0,1)=1
            call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
            col_fac(3)=dble(col_factor)
         elseif (this%n_qqbar.eq.2) then
            if (ui.eq.uj.and.ui.eq.1) then
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=2
               Tr(0,1,1) = gi+gj  ! number of generators in first trace
               Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
               Tr(1:gi,1,1) = iper(2:2+gi-1)
               Tr(gi+1:gi+gj,1,1) = jper(2+gj-1:2:-1)
               Tr(1:n-4-gi,2,1) = iper(2+gi+2:n-1)
               Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-1:2+gj+2:-1)
               coef(1)=1d0
               call Tr_full_simplify(col_factor)
            elseif ((ui.eq.2.and.uj.eq.1) .or. (ui.eq.1.and.uj.eq.2)) then
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=1 ! a single trace
               Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
               Tr(1:gi,1,1) = iper(2:2+gi-1)
               Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-1:2+gj+2:-1)
               Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(2+gi+2:n-1)
               Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(2+gj-1:2:-1)
               if (.not.this%same_flav) then
                  coef(1)=-1d0/3d0
               else
                  coef(1)=-1d0
               endif
               call Tr_full_simplify(col_factor)
            elseif (ui.eq.uj.and.ui.eq.2) then
               Tr(0,0,0) = 1 ! one term
               Tr(0,0,1) = 2 ! product of two traces
               Tr(0,1,1) = gi+gj  ! number of generators in first trace
               Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
               Tr(1:gi,1,1) = iper(2:2+gi-1)
               Tr(gi+1:gi+gj,1,1) = jper(2+gj-1:2:-1)
               Tr(1:n-4-gi,2,1) = iper(2+gi+2:n-1)
               Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-1:2+gj+2:-1)
               if (.not.this%same_flav) then
                  coef(1)=1d0/9d0
               else
                  coef(1)=1d0
               endif
               call Tr_full_simplify(col_factor)
            endif
            col_fac(3)=dble(col_factor)
         endif
         call Tr_deallocate
      endif
    end subroutine compute_color_factor

    subroutine convert_gluon_string(next,iper,jper,iper_new,jper_new)
      implicit none
      integer :: next
      integer,dimension(next-4) :: iper,jper,iper_new,jper_new
      integer :: ik,jk,i,j
      integer,dimension(next-4) :: imax,jmax
      imax(1:next-4)=-1
      jmax(1:next-4)=-1
      do i=1,next-4
         do j=1,next-4
            if ((iper(j).gt.imax(i)).and..not.(any(iper(j).eq.imax(1:i-1)))) imax(i)=iper(j)
            if ((jper(j).gt.jmax(i)).and..not.(any(jper(j).eq.jmax(1:i-1)))) jmax(i)=jper(j)
         enddo
      enddo
      ik=next-4
      jk=next-4
      do i=1,next-4
         do j=1,next-4
            if (iper(j).eq.imax(i)) then
               iper_new(j) = next-4-i+1
            endif
            if (jper(j).eq.jmax(i)) then
               jper_new(j) = next-4-i+1
            endif
         enddo
      enddo
    end subroutine convert_gluon_string

  end subroutine init_col

  subroutine filter_helicity(this,n,nhel,include_hel)
    implicit none
    class(amplitude_qcd) :: this
    integer,intent(inout) :: nhel
    integer,intent(in) :: n
    integer,intent(inout),dimension(nhel) :: include_hel
    integer :: ih,ih1,ih2,ihc,nspin,ispin,ic,iv
    logical,dimension(:),allocatable :: include_current,include_product
    integer,dimension(:,:,:),allocatable :: tmp_spin

    ! deallocate a bunch
    do ic=1,this%n_cur
       if (allocated(this%current_list(ic)%val_c)) deallocate(this%current_list(ic)%val_c)
       if (allocated(this%current_list(ic)%val_r)) deallocate(this%current_list(ic)%val_r)
    enddo
    do iv=1,this%n_vert
       if (allocated(this%interaction_list(iv)%val_c)) deallocate(this%interaction_list(iv)%val_c)
       if (allocated(this%interaction_list(iv)%val_r)) deallocate(this%interaction_list(iv)%val_r)
    enddo
    if (allocated(this%amps)) deallocate(this%amps)
    if (allocated(this%amps_r)) deallocate(this%amps_r)
    
    allocate(include_current(this%n_cur))
    include_current(this%n_cur_start(n  ):this%n_cur_end(n  ))=.false.
    include_current(this%n_cur_start(n-1):this%n_cur_end(n-1))=.false.
    allocate(include_product(nhel))
    include_product(1:nhel)=.false.
    
    ! Note: this must be done in the same order as the amps() are computed in 'compute_amps_from_currents'
    nspin=0
    do ih1=1,this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
       do ih2=1,this%n_cur_end(n)-this%n_cur_start(n)+1
          ih=(ih1-1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)+ih2
          if (include_hel(ih).ge.1) then
             include_current(this%n_cur_start(n-1)+ih1-1)=.true.
             include_current(this%n_cur_start(n  )+ih2-1)=.true.
             include_product(ih)=.true.
             nspin=nspin+1
          endif
       enddo
    enddo

    allocate(tmp_spin(1:n,1:maxval(include_hel),nspin))
    nspin=nhel

    ihc=0
    nhel=0
    do ih1=1,this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
       if (.not.include_current(this%n_cur_start(n-1)+ih1-1)) cycle
       do ih2=1,this%n_cur_end(n)-this%n_cur_start(n)+1
          if (.not.include_current(this%n_cur_start(n)+ih2-1)) cycle
          ih=(ih1-1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)+ih2
          ihc=ihc+1
          include_product(ihc)=include_product(ih)
          if (include_product(ihc)) then
             nhel=nhel+1
             tmp_spin(1:n,1,nhel)=this%spins(1:n,1,ih)
             if (include_hel(ih).gt.1) then
                ic=1
                do ispin=ih+1,nspin
                   if (-include_hel(ispin).eq.ih) then
                      ic=ic+1
                      if (ic.gt.include_hel(ih)) then
                         write (*,*) 'inconsistent include_hel #1'
                         stop 1
                      endif
                      tmp_spin(1:n,ic,nhel)=this%spins(1:n,1,ispin)
                   endif
                enddo
             elseif (include_hel(ih).ne.1) then
                write (*,*) 'inconsistent include_hel #2'
                stop 1
             endif
             include_hel(nhel)=include_hel(ih)
          endif
       enddo
    enddo
    deallocate(this%spins)
    call move_alloc(tmp_spin,this%spins)

    deallocate(this%include_product)
    allocate(this%include_product(1:ihc))
    this%include_product(1:ihc)=include_product(1:ihc)
    
    call this%filter_dead_trees(n,include_current)
    write (*,*) 'Total number of currents and vertices after filtering helicities',this%n_cur,this%n_vert,ihc,nhel

  end subroutine filter_helicity

  
  subroutine filter_dead_trees(this,n,include_current)
    ! some currents can be removed, since the "tree" starting from some of
    ! the initial state particles might lead to a dead end with no possible
    ! interactions for that current and the remaining external
    ! particles. Hence, they do not need to be computed since they can not
    ! lead to a valid Feynman diagram. To filter them out, one starts at the
    ! end, and goes backwards through the list and see if there are any
    ! currents that were not needed (i.e., they are not the input to a
    ! vertex that is used anywhere).
    implicit none
    class(amplitude_qcd) :: this
    logical,dimension(:),allocatable :: is_needed_cur,is_needed_ver
    integer,dimension(:),allocatable :: where_to_cur,where_to_ver
    logical,dimension(*),optional :: include_current
    integer :: to_skip,isize,nc,iv,n
    allocate(is_needed_cur(this%n_cur))
    allocate(is_needed_ver(this%n_vert))
    allocate(where_to_cur(this%n_cur))
    allocate(where_to_ver(this%n_vert))
    ! assume nothing is needed
    is_needed_cur(:)=.false.
    is_needed_ver(:)=.false.
    where_to_cur=0
    where_to_ver=0
    if (.not.present(include_current)) then
       is_needed_cur(this%n_cur_start(n-1):this%n_cur_end(n-1))=.true.
       is_needed_cur(this%n_cur_start(n  ):this%n_cur_end(n  ))=.true.
    else
       is_needed_cur(this%n_cur_start(n-1):this%n_cur_end(n-1))=include_current(this%n_cur_start(n-1):this%n_cur_end(n-1))
       is_needed_cur(this%n_cur_start(n  ):this%n_cur_end(n  ))=include_current(this%n_cur_start(n  ):this%n_cur_end(n  ))
    endif
    ! loop through the list backward: if we got to the end, i.e.,
    ! nc.ge.this%n_cur_start(n-1), we know that it is a valid tree. This
    ! means that all inputs to that final current are also needed. By moving
    ! backwards through the list, we can filter out all the branches of the
    ! tree that are needed. (The nc.le.this%n_cur_end(1) is needed only to
    ! make sure that the current that corresponds to the n'th final state
    ! particle is marked as needed, since that current does not enter any of
    ! the trees: it is only needed to close the tree and get the amplitudes.)
    do nc=this%n_cur,1,-1
       if (is_needed_cur(nc)) then
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
    do isize=1,n
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
       if (isize.ge.2 .and. isize.le.n-1) then
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



end module amplitude_QCD_mod
