module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: color_flow=.false.
  type current
     integer :: type,bin,nhel,n_vert
     integer,dimension(:),allocatable :: vertices,order
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
     real(kind=8) :: mass,width
  end type current
  type interaction
     integer :: type
     integer,dimension(2) :: currents
     integer,dimension(:),allocatable :: singlet_mv
     complex(kind=8),dimension(:,:),allocatable :: val_c
     real(kind=8),dimension(:,:),allocatable :: val_r
  end type interaction
  type amplitude_QCD
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: pp
     real(kind=8),dimension(:,:,:),allocatable :: diff_col_vals
     integer :: n_cur,n_vert,imode,nColOrd,n_qqbar,max_pp,n_sing
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end,helmap, &
          pp_bin_to_i,pp_i_to_bin
     integer,dimension(:,:),allocatable ::  n_col_vals
     integer,dimension(:,:),allocatable :: perm
     integer,dimension(:),allocatable :: quark_index
     integer,dimension(:,:,:,:),allocatable :: row_index,col_index
     integer,dimension(:,:,:),allocatable :: u1_lin_comb
     integer,dimension(:,:),allocatable :: buff
     logical :: same_flav
   contains
     procedure :: init,evaluate,init_col2
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,orig_part,part,mass,width,order,it)
    use math_functions
    implicit none
    class(amplitude_QCD) :: this
    integer::n,imode
    integer,dimension(n)::part,orig_part,order
    real(kind=8),dimension(n) :: mass,width
    integer :: isize,nc,isplit,n1,n2,ic1,ic2,iv,i,max_cur,max_vert,max_key
    real(kind=4) :: tAfter,tBefore
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:),allocatable :: key_to_current
    integer :: it ! quark order type

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
       call define_canonical_color_order(it)
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
       if (max_key.ne.max_cur) then
          write (*,*) 'Number of dictionary keys is expected to be identical to the maximum number of currents',&
               max_key,max_cur
          !stop 1
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
             this%current_list(this%n_cur)%mass=mass(order(nc))
             this%current_list(this%n_cur)%width=width(order(nc))
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
    write (*,*) 'Total number of currents and vertices before filter',this%n_cur,this%n_vert
    call filter_dead_trees()
    write (*,*) 'Total number of currents and vertices',this%n_cur,this%n_vert
    if (this%imode.eq.1) call create_helicity_map()
    if (this%imode.eq.2) call allocate_and_fill_colour_permutations()
    call setup_momentum_array()
  contains
    subroutine allocate_and_fill_colour_permutations()
      implicit none
      integer :: ind
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
      elseif (this%n_qqbar.eq.2) then
         allocate(this%perm(1:n-1,1:this%nColOrd))
         call setup_buff_2qq()
         do nc=this%n_cur_start(n-1),this%n_cur_end(n-1)
            call check_2qq_order(nc,ind)
            !this%perm(1:n-1,nc-this%n_cur_start(n-1)+1)=this%current_list(nc)%order(1:n-1)
            this%perm(1:n-1,ind) = this%current_list(nc)%order(1:n-1)
         enddo
      endif

      do nc=1,this%n_cur_end(n-1)-this%n_cur_start(n-1)+1
      enddo
    end subroutine allocate_and_fill_colour_permutations

    subroutine setup_buff_2qq
      implicit none
      integer :: i,j,k,m,q
      integer,dimension(1:n-4-this%n_sing) :: first,perm_out,perm_in
      integer,dimension(1:n-4-this%n_sing) :: gluons

      if (n-4-this%n_sing.gt.0) then
      allocate(this%buff(1:n-2-this%n_sing,(n-2-this%n_sing)*factorial(n-4-this%n_sing))) ! for photons: change
      this%buff(:,:)=0

      k=1
      do i=1,n
         if (part(i).eq.21) then
                 gluons(k)=i
                 k=k+1
          endif
      enddo

      m = 1
      do i=1,n-2-this%n_sing-1
         do j=1,n-4-this%n_sing
           first(j) = j
         enddo
         do k=1,factorial(n-4-this%n_sing)
           call get_next_iperm(n-4-this%n_sing,first,perm_out,n-4-this%n_sing)
           do q=1,n-4-this%n_sing
              perm_in(q) = gluons(perm_out(q))
           enddo
           this%buff(1:i-1,m) = perm_in(1:i-1)
           this%buff(i+2:n-2-this%n_sing,m) = perm_in(i:n-4-this%n_sing)
           first = perm_out
           m = m+1
         enddo
      enddo

      endif
    end subroutine setup_buff_2qq

    subroutine check_2qq_order(nc,ind)
      implicit none
      integer :: nc,ind
      integer :: i
      integer,dimension(1:n-2-this%n_sing) :: ord

      if (n-4-this%n_sing.gt.0) then
      ord = this%current_list(nc)%order(2:n-1)
      do i=1,n-2
        if ((abs(part(ord(i))).ge.1.and.abs(part(ord(i))).le.6)) then
              ord(i)=0
        endif
      enddo

      do i=1,(n-2-this%n_sing)*factorial(n-4-this%n_sing)
          if (all(this%buff(:,i).eq.ord)) then
                  ind = i
                  return
          endif
      enddo
      
      else
      ind = 1 

      endif
    end subroutine check_2qq_order

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
         !stop 1
      endif
    end subroutine simple_consistency_checks
  
    subroutine define_canonical_color_order(it)
      ! canonical order: (q,glu,glu,glu,singlet,singlet,qbar,q,qbar)
      use math_functions
      implicit none
      integer :: i
      integer :: nq,naq,nglu,nsing,iq,iaq,iglu,ising
      integer,dimension(n) :: input
      integer :: it ! quark order type

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
      order= 0
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
      integer,dimension(6) :: quark_flav
      integer :: i,j,k,sgn
      integer,dimension(8) :: flav  ! fills the flavours of quarks it finds ( in abs)
      integer,dimension(n) :: temp_part

      this%n_qqbar=0
      this%n_sing=0
      quark_flav=0
      this%same_flav=.true.
      flav = 0
      k = 1
      do i=1,n
         if (i.le.2) then
            if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
               quark_flav(abs(orig_part(i)))=quark_flav(abs(orig_part(i)))-sign(1,part(i))
               flav(k) = abs(orig_part(i))
               k= k+1
               if (orig_part(i).lt.0) this%n_qqbar=this%n_qqbar+1
            endif
         else
            if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
               quark_flav(abs(orig_part(i)))=quark_flav(abs(orig_part(i)))+sign(1,orig_part(i))
               flav(k) = abs(orig_part(i))
               k= k+1
               if (orig_part(i).gt.0) this%n_qqbar=this%n_qqbar+1
            endif
         endif
         if (orig_part(i).ne.21 .and. abs(orig_part(i)).gt.6) this%n_sing=this%n_sing+1
      enddo

      if (any(flav(1:2*this%n_qqbar).ne.flav(1))) this%same_flav = .false.

      allocate(this%quark_index(2*this%n_qqbar))

      temp_part = part
      do i=1,n
       if ((i.le.2).and.(abs(temp_part(i)).le.6)) then
               temp_part(i) = -temp_part(i)
        endif
      enddo

      if (this%n_qqbar.ge.2) then
      k=1
      do i=1,n
         if((temp_part(i).ge.1).and.(temp_part(i).le.6)) then
           this%quark_index(k)=i
           k=k+2
         endif
      enddo
      k=2
      do i=1,n
         if ((temp_part(i).le.-1).and.(temp_part(i).ge.-6)) then
         if (temp_part(i).eq.-temp_part(this%quark_index(1))) then
           this%quark_index(k)=i
         elseif (temp_part(i).eq.-temp_part(this%quark_index(3))) then
           this%quark_index(k+2)=i
         endif
         endif
      enddo

      do i=1,n
       if (k.ge.2) then
       if ((temp_part(i).eq.-temp_part(this%quark_index(k-1)))) then
          this%quark_index(k)=i
          k=k+1
          exit
       endif
       else
          this%quark_index(k)=i
       endif
      enddo

      endif

      if (any(quark_flav(:).ne.0)) then
         write (*,*) 'ERROR: inconsistent quark flavours',part(1:n)
         write(*,*) quark_flav
         stop 1
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
            if (i.eq.1 .or. i.eq.n) cycle
            if (order(i).eq.2) then
               if (part(order(i)).ge.1 .and. part(order(i)).lt.6) then
                  ! next should be a quark
                  if (.not.(part(order(i+1)).ge.1 .and. part(order(i+1)).lt.6)) then
                     write (*,*) 'ERROR: in the colour order, after an initial state quark should come a final state quark'
                     write (*,*) order
                     write (*,*) part
                     stop 1
                  endif
               endif
            else
               if (part(order(i)).le.-1 .and. part(order(i)).gt.-6) then
                  ! next should be a quark
                  if (.not.(part(order(i+1)).ge.1 .and. part(order(i+1)).lt.6)) then
                     write (*,*) 'ERROR: in the colour order, after a final state anti-quark should come a quark'
                     write (*,*) order
                     write (*,*) part
                     stop 1
                  endif
               endif
            endif
         enddo
      endif

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
               if (this%n_qqbar.eq.2) then
                 if (n-isize.ge.4) then
                   max_cur = max_cur+(n-isize-3)*2+3 ! pure gluon currents
                   max_cur=max_cur+3 ! gq + qq currents
                 elseif (n-isize.ge.3) then
                   max_cur=max_cur+3
                 elseif (n-isize.ge.2) then
                   max_cur=max_cur+1
                 endif
               endif
            endif
         enddo
         max_cur=max_cur+1
         if (this%n_sing.ge.1) max_cur=max_cur*this%n_sing
         max_cur=max_cur*3 ! to REMOVE
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
            ! -Number of gluon currents (remove the '/2' if use_symmetry=.false.):
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

         elseif (this%n_qqbar.eq.2) then
            max_cur=0
            do isize=1,n-1
               ! gluons and tensors
               ifact=n-4 ! n-4 gluons in total -> if no singlets around!!
               do j=1,isize-1
                  ifact=ifact*(n-4-j)
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
               ! single quark currents
               ifact=1
               do j=1,isize-1
                  ifact=ifact*(n-4-j+1)
               enddo
               max_cur=max_cur+3*ifact ! factor 3 for 3 q/aq externals
               enddo
               ! double quark currents
               if (isize.ge.2) then
                  max_cur=max_cur+2*1 ! for q-aq of both flavors
                  ifact=1
                  do j=1,isize-1
                    ifact=ifact*(n-4-j+1)
                  enddo
                  !max_cur=max_cur+2*ifact ! for both flavors
               endif
               ! three quark currents
               if (isize.ge.3) then
                  max_cur=max_cur+2*1 
                  ifact=1
                  do j=1,isize-1
                    ifact=ifact*(n-4-j+1)
                  enddo
                  max_cur=max_cur+2*ifact ! for both flavors
               endif
               ! four quark currents
               if (isize.ge.4) then
                  max_cur=max_cur+1 ! only one 4-q
                  ifact=1
                  do j=1,isize-1
                    ifact=ifact*(n-4-j+1)
                  enddo
                  max_cur=max_cur+ifact 
               endif
               max_cur=max_cur+1
               max_cur=max_cur+50 ! TO CHANGE
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
         if (this%n_sing.ge.1) max_vert=max_vert*this%n_sing
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
            mv=10*mv ! TO CHANGE!
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
      if (.not.valid_current_combination())  then
        return
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

      elseif (this%current_list(ic1)%type.eq.22 .and. &
           (this%current_list(ic2)%type.ge.1 .and. this%current_list(ic2)%type.le.6)) then
         ! add a photon-quark to quark vertex
         call add_vertex(4,this%current_list(ic2)%type)

      elseif ((this%current_list(ic1)%type.ge.1 .and. this%current_list(ic1)%type.le.6) .and. &
           this%current_list(ic2)%type.eq.22) then
         ! add a quark-photon to quark vertex
         call add_vertex(6,this%current_list(ic1)%type)

      elseif (this%current_list(ic1)%type.eq.22 .and. &
           (this%current_list(ic2)%type.le.-1 .and. this%current_list(ic2)%type.ge.-6)) then
         ! add a photon-antiquark to quark vertex
         call add_vertex(5,this%current_list(ic2)%type)

      elseif ((this%current_list(ic1)%type.le.-1 .and. this%current_list(ic1)%type.ge.-6) .and. &
           this%current_list(ic2)%type.eq.22) then
         ! add a antiquark-photon to quark vertex
         call add_vertex(7,this%current_list(ic1)%type)
      endif
    end subroutine add_if_allowed_threevertex

    
    logical function valid_current_combination()
      ! Checks to see if the combination of currents ic1 and ic2 could be a
      ! valid combination. Checks to perform:
      ! 0. All particles must be different in the two currents & final
      !    particle should never be part of the combined currents (it will be
      !    used to close the currents)
      ! 1. For imode=1 or 3, we need to make sure that the colour order is compatible with the input colour order
      ! 2. For imode=2 and use_symmetry=.true. --> skip one of the two colour orders if all gluons in current
      ! 3. For colour singlets:
      !  --> if one (or two) of the two currents is (are) a singlet, only include one of the two orders
      !  --> order of singlets themselves should be ignored in this check: all must be included
      ! 4. For quark currents:
      ! --> Only the first of the two currents can be a quark current, and for
      !     that current the first particle in the colour order must be a
      !     quark.
      ! 5. For anti-quark currents ???
      implicit none
      integer :: i,j,k,nc1,nc2
      logical :: gluon_current,colour_singlet1,colour_singlet2
      integer,dimension(isize) :: ip
      valid_current_combination=.false.
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(this%current_list(ic1)%bin,this%current_list(ic2)%bin)).ne.isize) return
      
      ! final particle should never be part of any combined currents: it will
      ! be used to close the amplitude instead
      if (n1.eq.1) then
         if (this%current_list(ic1)%order(n1).eq.order(n)) then
                 return
          endif
      endif
      if (n2.eq.1) then
         if (this%current_list(ic2)%order(n2).eq.order(n)) then
                 return
         endif
      endif

      colour_singlet1=all_singlet_current(this%current_list(ic1)%bin)
      colour_singlet2=all_singlet_current(this%current_list(ic2)%bin)
      ! If the first current is a singlet and the second is not, it is not a valid order
      if (colour_singlet1 .and. (.not.colour_singlet2)) return
      ! If both currents are colour singlets, only consider one of the two. 
      if (colour_singlet1 .and. colour_singlet2) then
         if (maxval(this%current_list(ic1)%order(1:n1)).ge.maxval(this%current_list(ic2)%order(1:n2))) then
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
            if (part(this%current_list(ic1)%order(i)).ge.22) exit
         enddo
         nc1=i-1
         do i=1,n2
            if (part(this%current_list(ic2)%order(i)).ge.22) exit
         enddo
         nc2=i-1
         ip(1:nc1+nc2)=[this%current_list(ic1)%order(1:nc1),this%current_list(ic2)%order(1:nc2)]
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

      ! The second current cannot be a quark current
      !if ( this%current_list(ic2)%type.ge.1 .and. this%current_list(ic2)%type.le.6) return
      ! If the first current is a quark current, the first particle in the order must be the quark
      !if ((this%current_list(ic1)%type.ge.1 .and. this%current_list(ic1)%type.le.6) .and. &
      !     this%current_list(ic1)%order(1).ne.order(1)) return

      ! If using symmetry and the current is a combination of all external
      ! gluons, take only one of the two possible orders
      gluon_current=all_gluon_current(this%current_list(ic1)%bin+this%current_list(ic2)%bin)
      if (use_symmetry .and. this%imode.eq.2 .and. gluon_current) then
         if (maxval(this%current_list(ic1)%order(1:n1)).ge.maxval(this%current_list(ic2)%order(1:n2))) return
      endif
      if (.not. gluon_current) then
         ! if quark is in there, it should be the very first particle
         if (any(this%current_list(ic1)%order(1:n1).eq.order(1)) .and. &
              this%current_list(ic1)%order(1).ne.order(1)) then
              !stop 10
              return
         endif
         if (any(this%current_list(ic2)%order(1:n2).eq.order(1))) then
            !stop 11
            return
         endif
      endif
      valid_current_combination=.true.

    end function valid_current_combination
    
    subroutine add_vertex(itype,ctype)
      implicit none
      integer :: itype,ctype
      if (isize.eq.n-1 .and. ctype.ne.anti_current(this%current_list(n)%type)) then 
        return ! dead tree. Filter already here
      endif
      this%n_vert=this%n_vert+1
      this%interaction_list(this%n_vert)%type=itype
      this%interaction_list(this%n_vert)%currents(1)=ic1
      this%interaction_list(this%n_vert)%currents(2)=ic2
      allocate(this%interaction_list(this%n_vert)%singlet_mv(0:isize))
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

    subroutine add_all_currents(ctype)
      implicit none
      logical,dimension(8) :: vertex_sign
      integer,dimension(isize,8) :: ip
      integer :: i,cur_bin,ctype,nperm
      integer,dimension(0:isize) :: singlet_mv
      if (.not.use_symmetry .or. this%imode.eq.1 .or. this%imode.eq.3) then
         cur_bin=this%current_list(ic1)%bin+this%current_list(ic2)%bin
         ip(1:isize,1)=combined_currents(n1,n2,this%current_list(ic1)%order(1:n1), &
              this%current_list(ic2)%order(1:n2),singlet_mv)
         this%interaction_list(this%n_vert)%singlet_mv(0:isize)=singlet_mv(0:isize)
         call add_current(.false.,cur_bin,ip(1:isize,1),ctype)
         return
      endif

      ! Need to consider all the possible permutations
      call check_all_permutations(nperm,ip,vertex_sign)
      cur_bin=this%current_list(ic1)%bin+this%current_list(ic2)%bin
      do i=1,nperm
         call add_current(vertex_sign(i),cur_bin,ip(1:isize,i),ctype)
      enddo
    end subroutine add_all_currents

    subroutine check_all_permutations(nperm,ip,vertex_sign)
      ! If a current only contains (external) gluons, we can use symmetry to
      ! relate them to eachother. This subroutine checks all permutations,
      ! and, if they give a valid current order, adds that current to the list
      ! that should be included.
      implicit none
      integer,intent(out) :: nperm
      integer,intent(out),dimension(isize,8) :: ip
      logical,intent(out),dimension(8) :: vertex_sign
      integer,dimension(1:n1,2) :: ip1
      integer,dimension(1:n2,2) :: ip2
      logical :: ag1,ag2,iden
      integer,dimension(3) :: switch
      integer :: i,j,k
      integer,dimension(0:isize,8) :: singlet_mv
      switch(1:3)=1
      ag1=all_gluon_current(this%current_list(ic1)%bin)
      ag2=all_gluon_current(this%current_list(ic2)%bin)
      if (n1.ge.2 .and. ag1) switch(1)=2
      if (n2.ge.2 .and. ag2) switch(2)=2
      if (ag1 .and. ag2) switch(3)=2
      ip1(1:n1,1)=this%current_list(ic1)%order(1:n1)
      ip1(1:n1,2)=this%current_list(ic1)%order(n1:1:-1)
      ip2(1:n2,1)=this%current_list(ic2)%order(1:n2)
      ip2(1:n2,2)=this%current_list(ic2)%order(n2:1:-1)
      nperm=0
      do i=1,switch(1)
         do j=1,switch(2)
            do k=1,switch(3)
               nperm=nperm+1

               if (k.eq.1) then
                  ip(1:isize,nperm)=combined_currents(n1,n2,ip1(1:n1,i),ip2(1:n2,j),&
                       singlet_mv(0,nperm))
               else
                  ip(1:isize,nperm)=combined_currents(n2,n1,ip2(1:n2,j),ip1(1:n1,i),&
                       singlet_mv(0,nperm))
               endif
               vertex_sign(nperm)=(k.eq.2 .xor. (j.eq.2 .and. mod(n2,2).eq.0) .xor. (i.eq.2 .and. mod(n1,2).eq.0))
               if (.not.valid_current_order(ip(1:isize,nperm))) nperm=nperm-1
            enddo
         enddo
      enddo

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
         this%interaction_list(this%n_vert)%singlet_mv(0:singlet_mv(0,1))=singlet_mv(0:singlet_mv(0,1),1)
      else
         write (*,*) 'Singlet move not identical for all permutations',nperm
         do i=1,nperm
            write (*,*) nperm,':',singlet_mv(0,i),':',singlet_mv(1:singlet_mv(0,i),i)
         enddo
         stop 1
      endif
    end subroutine check_all_permutations
    
    logical function valid_current_order(ip)
      ! Checks that ip(1:isize) is an order for a current to be considered
      ! when use_symmetry=.true. --> the smallest number needs to come before
      ! the largest number in this list.
      implicit none
      integer :: i,j,maxi,mini,min_loc,max_loc,k,length
      integer,dimension(isize) :: ip,ip_q
      integer,dimension(2*this%n_qqbar) :: quarks
      
      if (this%n_qqbar.eq.1 .and. (any(ip(2:isize).eq.order(1)))) then
         ! if there is a quark, it needs to be in the first position
         valid_current_order=.false.
         !stop 10
         return
      endif
      if (this%n_qqbar.eq.1 .and. ip(1).eq.order(1)) then
         ! if there is a quark, and it is part of the current (it must be at
         ! position 1), then it is a valid order, since no symmetry can be
         ! used.
         valid_current_order=.true.
         return
      endif

      k=1
      do i=1,n
         if (abs(part(order(i))).ge.1.and.abs(part(order(i))).le.6) then
                 quarks(k)=order(i)
                 k=k+1
         endif
      enddo

      k=1
      ip_q=0
      if (this%n_qqbar.ge.2) then
         do i=1,isize
            if (any(quarks.eq.ip(i))) then
                 ip_q(k)=ip(i)
                 k=k+1
            endif
         enddo
         length=k-1
      endif


      if (this%n_qqbar.ge.2) then
         do i=1,4
           if (ip_q(1).eq.quarks(i)) then
              if (i+length-1.gt.4) then
                 valid_current_order=.false.
                 return
              endif

              if ((any(ip.eq.quarks(2))).and.(any(ip.eq.quarks(3)))) then ! gluons cannot be between the two internal quark pair in order
                  do j=1,isize
                     if (ip(j).eq.quarks(2)) then
                        if (j.eq.isize) then
                            valid_current_order=.false.
                            return
                        endif
                        if (ip(j+1).ne.quarks(3)) then
                                valid_current_order=.false.
                                return
                        endif
                     endif    
                  enddo
              endif
              if ((any(ip.eq.quarks(2))).and..not.(any(ip.eq.quarks(3)))) then ! gluons cannot be between the two internal quark pair in order
                 if ((ip(isize).ne.quarks(2))) then
                   valid_current_order=.false.
                   return
                 endif
              endif
              if ((any(ip.eq.quarks(3))).and..not.(any(ip.eq.quarks(2)))) then ! gluons cannot be between the two internal quark pair in order
                 if ((ip(1).ne.quarks(3))) then
                   valid_current_order=.false.
                   return
                 endif
              endif


              if (all(ip_q(1:length).eq.quarks(i:i+length-1))) then
                  valid_current_order=.true.
                  return
              else
                  valid_current_order=.false.
                  return
              endif
              
           endif
         enddo
      endif
     
      if (this%n_qqbar.eq.2) then
      if (quark_in_current(ip,isize).le.4) then ! one quark only-> OK
         valid_current_order=.true.
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
         return
      endif
      if (quark_in_current(ip,isize).ge.11.and.&
              quark_in_current(ip,isize).le.14) then ! three quarks -> OK
         valid_current_order=.true.
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
        return
      endif
      if ((quark_in_current(ip,isize).eq.5 &
              .or.quark_in_current(ip,isize).eq.10))  then
         valid_current_order=.true.
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
      return
      endif

      if (quark_in_current(ip,isize).eq.15) then ! four quarks -> OK
         valid_current_order=.true.
         return
      endif
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

    integer function quark_in_current(ip,isize)
      implicit none
      integer :: isize,i,j
      integer,dimension(isize) :: ip

      quark_in_current=0
      do i=1,isize
        ! q
        if(ip(i).eq.this%quark_index(1).and.(.not.any(ip(1:isize).eq.this%quark_index(2))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(3))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=1
            exit
        endif
        ! qx
        if(ip(i).eq.this%quark_index(2).and.(.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(3))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=2
            exit
        endif
        ! q'
        if(ip(i).eq.this%quark_index(3).and.(.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(2))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=3
            exit
        endif
        ! qx'
        if(ip(i).eq.this%quark_index(4).and.(.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(2))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(3)))) then
            quark_in_current=4
            exit
        endif
        ! q-qx
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(2)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(3))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=5
            exit
        endif
        ! q-q'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(3)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(2))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=6
            exit
        endif
        ! q-qx'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(4)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(2)))) then
            quark_in_current=7
            exit
        endif
        ! qx-q'
        if(ip(i).eq.this%quark_index(2).and.any(ip(1:isize).eq.this%quark_index(3)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=8
            exit
        endif
        ! qx-q'
        if(ip(i).eq.this%quark_index(2).and.any(ip(1:isize).eq.this%quark_index(4)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(3)))) then
            quark_in_current=9
            exit
        endif
        ! q'-qx'
        if(ip(i).eq.this%quark_index(3).and.any(ip(1:isize).eq.this%quark_index(4)).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(1))).and.&
               (.not.any(ip(1:isize).eq.this%quark_index(2)))) then
            quark_in_current=10
            exit
        endif
        ! q-qx-q'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(2)).and.&
               any(ip(1:isize).eq.this%quark_index(3)).and.&
              (.not.any(ip(1:isize).eq.this%quark_index(4)))) then
            quark_in_current=11
            exit
        endif
        ! q-qx-qx'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(2)).and.&
               any(ip(1:isize).eq.this%quark_index(4)).and.&
              (.not.any(ip(1:isize).eq.this%quark_index(3)))) then
            quark_in_current=12
            exit
        endif
        ! q-q'-qx'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(3)).and.&
               any(ip(1:isize).eq.this%quark_index(4)).and.&
              (.not.any(ip(1:isize).eq.this%quark_index(2)))) then
            quark_in_current=13
            exit
        endif
        ! qx-q'-qx'
        if(ip(i).eq.this%quark_index(2).and.any(ip(1:isize).eq.this%quark_index(3)).and.&
               any(ip(1:isize).eq.this%quark_index(4)).and.&
              (.not.any(ip(1:isize).eq.this%quark_index(1)))) then
            quark_in_current=14
            exit
        endif
        ! q-qx-q'-qx'
        if(ip(i).eq.this%quark_index(1).and.any(ip(1:isize).eq.this%quark_index(2)).and.&
               any(ip(1:isize).eq.this%quark_index(3)).and.&
               any(ip(1:isize).eq.this%quark_index(4))) then
            quark_in_current=15
            exit
        endif
      enddo

    end function quark_in_current


    subroutine add_current(vertex_sign,cur_bin,ip,ctype)
      implicit none
      logical :: vertex_sign
      integer,dimension(isize) :: ip ! permutation of the current
      integer :: ctype,cur_bin,ic,key
      integer(kind=8) :: val
      if (this%imode.eq.1 .or. this%imode.eq.3) then
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
         if (this%current_list(ic1)%mass.eq.this%current_list(ic2)%mass)  then
            this%current_list(this%n_cur)%mass=0d0
         else
            this%current_list(this%n_cur)%mass=max(this%current_list(ic1)%mass,this%current_list(ic2)%mass)
         endif
         if (this%current_list(ic1)%width.eq.this%current_list(ic2)%width)  then
            this%current_list(this%n_cur)%width=0d0
         else
            this%current_list(this%n_cur)%width=max(this%current_list(ic1)%width,this%current_list(ic2)%width)
         endif
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
         elseif (ctype.ge.-6 .and. ctype.le.-1) then
            ! anti-quark current
            call get_value(ip,2,val)
         endif

         !do ic=1,this%n_cur
         !   if (cur_bin.eq.this%current_list(ic)%bin) return
         !enddo
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
               if (any(ips_in(1:isize).eq.order(n))) cycle ! should not contain last closing particle 
               if (this%n_qqbar.eq.0 .or. &
                  (this%n_qqbar.eq.1.and.all(ips_in(1:isize).ne.order(1)))) then
                  
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.1 .and. isize.ne.n-1) then ! add the tensor
                     key=key+1
                     call get_value(ips_in,-1,val)
                     current_dict(key)=val
                  endif

                elseif (this%n_qqbar.eq.2) then
                     if ((this%n_qqbar.eq.2.and.(.not.any(ips_in(1:isize).eq.this%quark_index(1)))).and.&
                         (this%n_qqbar.eq.2.and.(.not.any(ips_in(1:isize).eq.this%quark_index(2)))).and.&
                         (this%n_qqbar.eq.2.and.(.not.any(ips_in(1:isize).eq.this%quark_index(3)))).and.&
                         (this%n_qqbar.eq.2.and.(.not.any(ips_in(1:isize).eq.this%quark_index(4))))) then
                        key=key+1
                        call get_value(ips_in,0,val) ! add the gluon
                        current_dict(key)=val
                        if (isize.ne.1 .and. isize.ne.n-1) then ! add the tensor
                           key=key+1
                           call get_value(ips_in,-1,val)
                           current_dict(key)=val
                        endif
                      endif
               endif

               if (this%n_qqbar.eq.1 .and. ips_in(1).eq.order(1)) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif

               if (this%n_qqbar.eq.2) then
               if (quark_in_current(ips_in,isize).eq.1) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif
               if (quark_in_current(ips_in,isize).eq.3) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif
              
               if (quark_in_current(ips_in,isize).eq.2) then
                  key=key+1
                  call get_value(ips_in,2,val) ! add a quark
                  current_dict(key)=val
               endif
               if (quark_in_current(ips_in,isize).eq.4) then
                  key=key+1
                  call get_value(ips_in,2,val) ! add a quark
                  current_dict(key)=val
               endif

               if (quark_in_current(ips_in,isize).eq.5) then
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.n-1) then ! add the tensor
                     key=key+1
                     call get_value(ips_in,-1,val)
                     current_dict(key)=val
                  endif
               endif
               if (quark_in_current(ips_in,isize).eq.10) then
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.n-1) then ! add the tensor
                     key=key+1
                     call get_value(ips_in,-1,val)
                     current_dict(key)=val
                  endif
               endif

               if (quark_in_current(ips_in,isize).eq.11) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif

               if (quark_in_current(ips_in,isize).eq.12) then
                  key=key+1
                  call get_value(ips_in,2,val) ! add a a-quark
                  current_dict(key)=val
               endif

               if (quark_in_current(ips_in,isize).eq.13) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a quark
                  current_dict(key)=val
               endif

               if (quark_in_current(ips_in,isize).eq.14) then
                  key=key+1
                  call get_value(ips_in,1,val) ! add a a-quark
                  current_dict(key)=val
               endif

               if (quark_in_current(ips_in,isize).eq.15) then
                  key=key+1
                  call get_value(ips_in,0,val) ! add the gluon
                  current_dict(key)=val
                  if (isize.ne.n-1) then ! add the tensor
                     key=key+1
                     call get_value(ips_in,-1,val)
                     current_dict(key)=val
                  endif
               endif
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
      ! Take the types into account (we have only 4 types (gluon,
      ! tensor and quark and anti-quark), so multiply by three (and add one for the
      ! tensor and two for quark))
      val=val*int(4,kind=8) ! gluon
      if (itype.eq.-1) then
         val=val+int(1,kind=8) ! tensor
      endif
      if (itype.eq.1) then
         val=val+int(2,kind=8) ! quark
      endif
      if (itype.eq.2) then
         val=val+int(3,kind=8) ! anti-quark
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
      if (this%n_qqbar.le.1) then
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

      elseif (this%n_qqbar.eq.2) then
        middle = 0
        do middle=1,max_key
         if (current_dict(middle).eq.val) then
            key=middle
            return
         endif
        enddo
        key=0
      endif
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
  end subroutine init


  subroutine evaluate(this,n,p,mass,width,hel,part)
    use FeynmanRules
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,hel
    integer,dimension(n)::part
    real(kind=8),dimension(n) :: mass,width
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih1,ih2,ih,ih_in,ip,imv
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
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih),this%current_list(ic)%mass)
                elseif (this%current_list(ic)%type.ge.-6 .and. this%current_list(ic)%type.le.-1 ) then
                        call ext_antiquark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4,ih),this%current_list(ic)%mass)
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

                 elseif(this%interaction_list(iv)%type.eq.4) then
                    if (use_real_gluons) then
                       call GluonQuarktoQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                                   this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                   this%interaction_list(iv)%val_c(1:4,ih))
                    else
                       call GluonQuarktoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))
                    endif

                 elseif(this%interaction_list(iv)%type.eq.5) then
                    call GluonAquarktoAquark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                           this%interaction_list(iv)%val_c(1:4,ih))

                 elseif(this%interaction_list(iv)%type.eq.6) then
                    ! makes sure the helicities in 'ih' are correctly set in case colour singlets are moved to the end
                    do imv=1,this%interaction_list(iv)%singlet_mv(0)
                       call move_ih(this%interaction_list(iv)%singlet_mv(imv),ih)
                    enddo
                    if (use_real_gluons) then
                       call QuarkGluontoQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                                   this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4,ih2),&
                                                   this%interaction_list(iv)%val_c(1:4,ih))
                    else
                       call QuarkGluontoQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))
                    endif
                 elseif(this%interaction_list(iv)%type.eq.7) then
                    ! makes sure the helicities in 'ih' are correctly set in case colour singlets are moved to the end
                    do imv=1,this%interaction_list(iv)%singlet_mv(0)
                       call move_ih(this%interaction_list(iv)%singlet_mv(imv),ih)
                    enddo
                    call AquarkGluontoAquark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                          this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                          this%interaction_list(iv)%val_c(1:4,ih))
                elseif(this%interaction_list(iv)%type.eq.8) then
                    ! makes sure the helicities in 'ih' are correctly set in case colour singlets are moved to the end
                    do imv=1,this%interaction_list(iv)%singlet_mv(0)
                       call move_ih(this%interaction_list(iv)%singlet_mv(imv),ih)
                    enddo

                   call QuarkAquarktoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))

                elseif(this%interaction_list(iv)%type.eq.9) then
                    ! makes sure the helicities in 'ih' are correctly set in case colour singlets are moved to the end
                    do imv=1,this%interaction_list(iv)%singlet_mv(0)
                       call move_ih(this%interaction_list(iv)%singlet_mv(imv),ih)
                    enddo
                    call AquarkQuarktoGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4,ih1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4,ih2),&
                                              this%interaction_list(iv)%val_c(1:4,ih))

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
    subroutine move_ih(imv,ih)
      ! If there is a colour singlet, we need to be careful, since it is
      ! always moved to the end of the color order. This might mess up the
      ! helicity assignment when summing over helicities. Hence, we need to
      ! move the helicity bit of the colour singlet to the end of 'ih' and
      ! move all the other bits one step towards the beginning.
      implicit none
      integer,intent(inout) :: ih
      integer,intent(in) :: imv
      integer :: ihm1
      integer :: ising
      ! when not summing over helicities (or when ih=1 (i.e. all bits are
      ! equal to zero)), there is nothing to do.
      if (ih.eq.1) return
      ihm1=ih-1
      ! ising will contain the helicity bit of the colour singlet that has
      ! been moved towards the end: this is the bit that needs to be put at
      ! the end of 'ih'
      call mvbits(ihm1,imv-1,1,ising,0)
      ! move all the other bits one step towards the beginning
      call mvbits(ihm1,imv, &
                       popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin+&
                              this%current_list(this%interaction_list(iv)%currents(2))%bin)-imv, &
                       ihm1,&
                       imv-1)
      ! place the colour singlet at the end
      call mvbits(ising,0,1,ihm1,popcnt(this%current_list(this%interaction_list(iv)%currents(1))%bin+&
                                        this%current_list(this%interaction_list(iv)%currents(2))%bin)-1)
      ih=ihm1+1
    end subroutine move_ih
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
         !this%pp(0:3,1)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
         !this%pp(0:3,4)=(/0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03/)
         !this%pp(0:3,2)=(/0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03/)
         !this%pp(0:3,3)=(/0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03/)

    end subroutine fill_momentum_array
    subroutine compute_amps_from_currents
      implicit none
      integer :: ind
      if (this%imode.eq.1) then
         do ih1=1,this%current_list(this%n_cur)%nhel ! helicities for combined current of particles 1 to n-1 (in the colour order)
            do ih2=1,this%current_list(n)%nhel       ! and for the current for particle n
               ih=(ih2-1)*this%current_list(this%n_cur)%nhel+ih1
               if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
                  this%amps_r(this%helmap(ih))=sum(this%current_list(this%n_cur)%val_r(1:4,ih1)*this%current_list(n)%val_r(1:4,ih2))
               else
                  this%amps(this%helmap(ih))=sum(this%current_list(this%n_cur)%val_c(1:4,ih1)*this%current_list(n)%val_c(1:4,ih2))
               !write(*,*) 'holahola0'
               !write(*,*) this%amps(this%helmap(ih))

               !write(*,*) this%current_list(n)%val_c(1:4,ih2)
               !write(*,*) this%current_list(this%n_cur)%val_c(1:4,ih1)
               endif
            enddo
         enddo

      elseif (this%imode.eq.2) then
         do ic=this%n_cur_start(n-1),this%n_cur_end(n-1)
            if (use_real_gluons .and. this%current_list(n)%type.eq.21) then
               this%amps_r(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_r(1:4,1)*this%current_list(n)%val_r(1:4,1))
            else
               if (this%n_qqbar.eq.2) then
                  call check_2qq_order(ic,ind)
                  this%amps(ind) = sum(this%current_list(ic)%val_c(1:4,1)*this%current_list(n)%val_c(1:4,1))
               else        
                  this%amps(ic-this%n_cur_start(n-1)+1)=sum(this%current_list(ic)%val_c(1:4,1)*this%current_list(n)%val_c(1:4,1))
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
            this%amps_r(1)=sum(this%current_list(this%n_cur)%val_r(1:4,1)*this%current_list(n)%val_r(1:4,2))
         else
            this%amps(1)=sum(this%current_list(this%n_cur)%val_c(1:4,1)*this%current_list(n)%val_c(1:4,1))
         endif
      endif
    end subroutine compute_amps_from_currents

    subroutine check_2qq_order(nc,ind)
      implicit none
      integer :: nc,ind
      integer :: i
      integer,dimension(1:n-2-this%n_sing) :: ord

      if (n-4-this%n_sing.gt.0) then
      ord = this%current_list(nc)%order(2:n-1)
      do i=1,n-2
        if ((abs(part(ord(i))).ge.1.and.abs(part(ord(i))).le.6)) then
              ord(i)=0
        endif
      enddo

      i = 1
      do 
          if (all(this%buff(:,i).eq.ord)) then
                  ind = i
                  return
          endif
          i = i+1
      enddo

      else
      ind = 1

      endif

    end subroutine check_2qq_order

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
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), & 
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_quark_propagator

    subroutine include_aquark_propagator()
      implicit none
      call AquarkPropagator(this%current_list(ic)%val_c,this%current_list(ic)%nhel, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_aquark_propagator
  end subroutine evaluate

  subroutine init_col2(this,n,part,order,it,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer :: col_acc,n
    integer,dimension(n) :: order,iper,jper,part
    integer :: iperm,jperm,ival,iacc
    !integer,dimension(1:3) :: n_vals
    integer,allocatable,dimension(:,:) :: n_vals
    integer,parameter :: max_vals=100
    real(kind=8),dimension(1:3) :: col_fac
    !real(kind=8),allocatable,dimension(:,:) :: col_fac
    !real(kind=8),dimension(max_vals,1:3) :: diff_vals
    real(kind=8),allocatable,dimension(:,:,:) :: diff_vals
    real(kind=8),allocatable,dimension(:,:) :: col_vals
    integer,dimension(:,:,:),allocatable :: ic,ir
    integer :: ri,rj,lim,y,t,maxterms_u1,i,j,gi
    integer :: ui,uj,uj_upper ! quark ordering type
    integer :: it ! quark order type (1 or 2)
    integer :: iperm_upper,gi_iperm  ! needed for 2qq
    integer :: val,max_val

    write (*,*) 'Initialising colour matrix ...'
    
    do i=0,n-1
        iper(i+1)=n-i
    enddo
    max_val=0
    do j=1,n
       max_val=max_val+int(iper(n+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
    enddo

    if (this%n_qqbar.eq.0) then
       allocate(n_vals(1:3,1))
       allocate(diff_vals(max_vals,1:3,1))
       allocate(col_vals(1:3,max_val))
       !allocate(col_fac(1:3,1))
       lim=0
       iperm_upper = 1
    elseif (this%n_qqbar.eq.1) then
       allocate(n_vals(1:3,1))
       allocate(diff_vals(max_vals,1:3,1))
       allocate(col_vals(1:3,max_val))
       !allocate(col_fac(1:3,1))
       lim=0
       iperm_upper = 1
       if (color_flow) then
          lim=1 ! for NLC only
          !lim=n-2 ! TV test: for FC
          lim = 0 ! if U(1) amps generated separately
          maxterms_u1 = 1d0
          do i=2,n-1
             maxterms_u1 = maxterms_u1 * (n-i)
          enddo
          allocate(this%u1_lin_comb(1:this%nColOrd,1:maxterms_u1,1:(n-2)))
          this%u1_lin_comb = 0
          call get_u1_lin_comb
       endif
    elseif (this%n_qqbar.eq.2) then
         lim = 0
         iperm_upper = (n-4)+1 ! number of gluon separations on two quark lines
         allocate(n_vals(1:3,iperm_upper))
         allocate(diff_vals(max_vals,1:3,iperm_upper))
         !allocate(col_fac(1:3,iperm_upper))
    endif

! first check a single row in the colour matrix to determine how many
! different colour factors there are
    n_vals(1:3,:)=0
    col_vals(1:3,1:max_val)=0d0
    do iperm=1,iperm_upper
    if (this%n_qqbar.eq.0) then
       iper(1:n)=[this%perm(1:n-1,iperm),n]
       gi_iperm = 1 ! dummy
    elseif (this%n_qqbar.eq.1) then
       iper(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,iperm),order(n)]
       gi_iperm = 1 ! dummy
    elseif (this%n_qqbar.eq.2) then
       do j=1,this%nColOrd
       iper(1:n-this%n_sing)=[this%perm(1:n-1-this%n_sing,j),order(n)] !
       do i=1,n-1
          if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
              if (i.ne.1) then
                 gi = i - 2
              exit
          endif
          endif
       enddo
       if (gi.eq.iperm-1) exit
       enddo
       gi_iperm = iperm
    endif

    if (this%n_qqbar.eq.2) then
        uj_upper = 2
    else
        uj_upper = 1
    endif

    do ri=0,lim ! number of U(1) gluons in amp
    do jperm=1,this%nColOrd ! TO CHANGE!
      do uj=1,uj_upper

       if (this%n_qqbar.eq.0) then
          jper(1:n)=[this%perm(1:n-1,jperm),n]
       elseif (this%n_qqbar.eq.1) then
          jper(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,jperm),order(n)]
       elseif (this%n_qqbar.eq.2) then
          jper(1:n-this%n_sing)=[this%perm(1:n-1-this%n_sing,jperm),order(n)]  ! last one is dummy
       endif

       val=0
       do j=1,n
          val=val+int(jper(n+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
       enddo

       do rj=0,lim
          if (this%n_qqbar.eq.2) then
              ui = it
          endif

          call compute_color_factor(col_acc,n-this%n_sing,iper,jper,ri,rj,ui,uj,col_fac,color_flow)
          col_fac(1:3)=col_fac(1:3)*2d0 ! include a factor 2 for the off-diagonal terms
          if (iperm.eq.jperm.and.ui.eq.uj) col_fac(1:3)=col_fac(1:3)*0.5d0 

          do iacc=1,3
            if (col_fac(iacc).eq.0d0) cycle
            col_vals(iacc,val)=col_fac(iacc)
             do ival=1,n_vals(iacc,gi_iperm)
               if (col_fac(iacc).eq.diff_vals(ival,iacc,gi_iperm)) exit
             enddo
             if (ival.ge.max_vals) then
               write (*,*) 'Too many different colour factors. Increase max_vals',&
                  ival,n_vals(1:3,gi_iperm),max_vals
               stop 1
             elseif (ival.eq.n_vals(iacc,gi_iperm)+1) then
               ! new colour factor
               n_vals(iacc,gi_iperm)=n_vals(iacc,gi_iperm)+1
               diff_vals(n_vals(iacc,gi_iperm),iacc,gi_iperm)=col_fac(iacc)
             endif
          enddo
       enddo
     enddo
    enddo
    enddo

    write (*,*) 'A single row in the colour matrix has',n_vals(1:3,gi_iperm),&
         ' different colour factors at LC, NLC and full colour, respectively'
    enddo


 ! Allocate the arrays now that we know their sizes
    allocate(ic(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    allocate(ir(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    if (color_flow) then
      allocate(this%col_index(1:((n-1)*this%nColOrd)**2,1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
      allocate(this%row_index(0:(n-1)*this%nColOrd,1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper)) ! for colour flow *(n-1)
    else
      allocate(this%col_index(1:(this%nColOrd)**2,1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
      allocate(this%row_index(0:this%nColOrd,1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper)) 
    endif
    this%row_index(0,1:maxval(n_vals(1:3,iperm_upper)),1:3,1:iperm_upper)=0
    this%col_index(1,1:maxval(n_vals(1:3,iperm_upper)),1:3,1:iperm_upper)=0
    allocate(this%n_col_vals(1:3,iperm_upper))
    this%n_col_vals(1:3,1:iperm_upper)=n_vals(1:3,1:iperm_upper)
    allocate(this%diff_col_vals(1:maxval(n_vals(1:3,iperm_upper)),1:3,iperm_upper))
    do iacc=1,3
       this%diff_col_vals(1:n_vals(iacc,iperm_upper),iacc,1:iperm_upper)=diff_vals(1:n_vals(iacc,iperm_upper),iacc,1:iperm_upper)
    enddo

! Compute all the colour factors and fill the col_index and row_index arrays
    ic=0
    ir=0

    do ri=0,lim
    do iperm=1,this%nColOrd
        
       if (this%n_qqbar.eq.0) then
          iper(1:n)=[this%perm(1:n-1,iperm),n]
       elseif (this%n_qqbar.eq.1) then
          iper(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,iperm),order(n)]
       elseif (this%n_qqbar.eq.2) then
          iper(1:n-this%n_sing)=[this%perm(1:n-1-this%n_sing,iperm),order(n)] 
          ! find out what channel it belongs to, find gi
          do i=1,n-1
          if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
              if (i.ne.1) then
                 gi = i - 2
              exit
              endif
          endif
          enddo
          gi_iperm = gi + 1
       endif

        ! for now, include all matrix ! TO CHANGE: put back off-diagonality
        do jperm=iperm,this%nColOrd ! only include upper triangle (i.e., loop starts at iperm instead of 1)
         do uj=1,uj_upper
          
          if (this%n_qqbar.eq.0) then
             jper(1:n)=[this%perm(1:n-1,jperm),n]
          elseif (this%n_qqbar.eq.1) then
             jper(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,jperm),order(n)]
          elseif (this%n_qqbar.eq.2) then
             ! all quark elements are dummy, only the gluon locations matter
             jper(1:n-this%n_sing)=[this%perm(1:n-1-this%n_sing,jperm),order(n)] 
          endif

          do rj=0,lim

            ! COMPUTE color factors again
            !call compute_color_factor(col_acc,n-this%n_sing,iper,jper,ri,rj,ui,uj,col_fac,color_flow)
            ! TO CHANGE: also here put back factor 2
            col_fac(1:3)=col_fac(1:3)*2d0
            if (iperm.eq.jperm.and.ui.eq.uj) col_fac(1:3)=col_fac(1:3)*0.5d0 ! include a factor 2 for the off-diagonal terms

            ! GET color factors from permuting first row
            call get_col_fac(col_fac)

            do iacc=1,3
              if (col_fac(iacc).eq.0d0) cycle
               do ival=1,n_vals(iacc,gi_iperm)
                 if (col_fac(iacc).eq.diff_vals(ival,iacc,gi_iperm)) then
                         exit
                 endif
               enddo 

               ic(ival,iacc,gi_iperm)=ic(ival,iacc,gi_iperm)+1
               ir(ival,iacc,gi_iperm)=ir(ival,iacc,gi_iperm)+1
               this%col_index(ic(ival,iacc,gi_iperm),ival,iacc,gi_iperm)=(rj*this%nColOrd)+((uj-1)*this%nColOrd)+jperm
            enddo

          enddo
          enddo

        enddo

        do iacc=1,3
          this%row_index((ri*this%nColOrd)+iperm,1:n_vals(iacc,gi_iperm),iacc,gi_iperm)=ir(1:n_vals(iacc,gi_iperm),iacc,gi_iperm)
        enddo
      enddo
    enddo
    write (*,*) '... colour matrix initialised'
  contains
   subroutine get_col_fac(col_fac)
     implicit none
     integer,dimension(n) :: col_new,row_first,row_per,col_per
     integer :: i,j,val
     real(kind=8),dimension(1:3) :: col_fac

     ! First row
     if (this%n_qqbar.eq.0) then
          row_first(1:n)=[this%perm(1:n-1,1),n]
     elseif (this%n_qqbar.eq.1) then
          row_first(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,1),order(n)]
     endif

     ! Row in consideration
     if (this%n_qqbar.eq.0) then
          row_per(1:n)=[this%perm(1:n-1,iperm),n]
     elseif (this%n_qqbar.eq.1) then
          row_per(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,iperm),order(n)]
     endif

     ! Column in consideration
     if (this%n_qqbar.eq.0) then
          col_per(1:n)=[this%perm(1:n-1,jperm),n]
     elseif (this%n_qqbar.eq.1) then
          col_per(1:n-this%n_sing)=[order(1),this%perm(1:n-2-this%n_sing,jperm),order(n)]
     endif

     do i=1,n
        do j=1,n
           if (col_per(i) .eq. row_per(j)) exit
        enddo
        col_new(i) = row_first(j)
     enddo

     val=0
     do j=1,n
          val=val+int(col_new(n+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
     enddo

     col_fac(1:3)=col_vals(1:3,val)
   end subroutine get_col_fac



    subroutine compute_color_factor(col_acc,n,iper,jper,ri,rj,ui,uj,col_fac,color_flow)
      use color_algebra
      implicit none
      integer :: i,n,acc,col_acc,color_fac,k
      real(kind=8),dimension(1:3) :: col_fac
      integer,dimension(n) :: iper,jper
      integer,dimension(n) :: iper_red, jper_red
      integer,dimension(n-4) :: iper_glu,jper_glu,iper_ord,jper_ord
      logical :: color_flow
      real(kind=16) :: col_factor
      integer :: ri,rj
      integer :: ui,uj ! for 2qq processes: type of quark ordering 
      integer :: gi,gj ! for 2qq: number of gluons on first quark color line

      if (color_flow) then
         if (this%n_qqbar.ne.0) then
            col_fac=0d0
            call convert_perm_color_flow(n,iper,jper,iper_red,jper_red)
            call color_flow_factor_1qqbar(n,iper_red,jper_red,ri,rj,color_fac)
            if (col_acc.ge.0) then
              if (ri.eq.rj.and.ri.eq.0.and.all(iper.eq.jper)) then
                    col_fac(1)=dble(3**color_fac) ! U(3)xU(3) LC
              endif
            endif

            if (col_acc.ge.1) then
            if (ri.eq.rj.and.ri.eq.0.and.all(iper.eq.jper)) then
                    col_fac(2)=dble(3**color_fac) ! U(3)xU(3) LC
            elseif (ri.eq.rj.and.ri.eq.0.and.color_fac.eq.(n-3)) then
                    col_fac(2)=dble(3**color_fac) ! U(3)xU(3) NLC
            ! comment out below...
            !elseif (((ri.eq.1.and.ri.ne.rj).or.(rj.eq.1.and.ri.ne.rj)).and.color_fac.eq.(n-3)) then
            !        col_fac(2)=-dble(3**color_fac)  ! U(1)xU(3) NLC
            ! ... and set this to (-)* (otherwise it is (+)*
            elseif (ri.eq.rj.and.ri.eq.1.and.color_fac.eq.(n-3)) then 
                    col_fac(2)=-dble(3**color_fac) ! U(1)xU(1) NLC
            endif
            endif
            if (col_acc.ge.2) then
               if (mod(ri+rj,2).eq.0) col_fac(3)=dble(3d0**color_fac)
               if (mod(ri+rj,2).eq.1) col_fac(3)=-dble(3d0**color_fac)
            endif
         else
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
                  coef(1)=(1d0,0d0)
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
                     coef(1)=(1d0,0d0)
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
                      coef(1)=(1d0,0d0)
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
                      coef(1)=(1d0,0d0)
                      if (.not.this%same_flav) coef(1)=((-1d0/3d0)**2)*(1d0,0d0)
                      call Tr_full_simplify(col_factor)
                  elseif (ui.eq.2.and.uj.eq.1) then
                      Tr(0,0,0)=1 ! one term
                      Tr(0,0,1)=1 ! a single trace
                      Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
                      Tr(1:gi,1,1) = iper(2:2+gi-1)
                      Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-1:2+gj+2:-1)
                      Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(2+gi+2:n-1)
                      Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(2+gj-1:2:-1)
                      coef(1)=(1d0,0d0)
                      if (.not.this%same_flav) coef(1)=((1d0/3d0))*(1d0,0d0)
                      coef(1)=coef(1)*(-1d0)
                      call Tr_full_simplify(col_factor)
                  elseif (ui.eq.1.and.uj.eq.2) then
                      Tr(0,0,0)=1 ! one term
                      Tr(0,0,1)=1 ! a single trace
                      Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
                      Tr(1:gi,1,1) = iper(2:2+gi-1)
                      Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-1:2+gj+2:-1)
                      Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(2+gi+2:n-1)
                      Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(2+gj-1:2:-1)
                      coef(1)=(1d0,0d0)
                      if (.not.this%same_flav) coef(1)=((1d0/3d0))*(1d0,0d0)
                      coef(1)=coef(1)*(-1d0)
                      call Tr_full_simplify(col_factor)
                  endif
                  col_fac(3)=dble(col_factor)
            endif
            call Tr_deallocate
         endif
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


    subroutine get_u1_lin_comb
      implicit none
      integer i,j,k,l,ii,m,top,add,num,j1,j2
      integer,dimension(:),allocatable :: u1_rem,u1_test
      integer,dimension(1:n-2-this%n_sing) :: perm,temp_perm

      do i=1,this%nColOrd
        do m=0,n-2 ! loop through all possible numbers of external U(1) gluons
          if(allocated(u1_rem)) deallocate(u1_rem)
          if(allocated(u1_test)) deallocate(u1_test)
          allocate(u1_rem(1:n-2-this%n_sing-m))
          allocate(u1_test(1:n-2-this%n_sing-m))
          ii=1
          perm(1:n-2-this%n_sing) = this%perm(1:n-2-this%n_sing,i)
          u1_rem(1:n-2-this%n_sing-m) = perm(1:n-2-this%n_sing-m)

          do j=1,this%nColOrd
             temp_perm(1:n-2-this%n_sing) = this%perm(1:n-2-this%n_sing,j)
             k=1
             do l=1,n-2-this%n_sing-m+1
               if (temp_perm(l).ne.perm(n-2-this%n_sing)) then
                  u1_test(k) = temp_perm(l)
                  k = k+1
               endif
             enddo

             if ((all(u1_test.eq.u1_rem))) then
              if (all(perm(n-2-this%n_sing-m+1:n-2-this%n_sing-1).eq.temp_perm(n-2-this%n_sing-m+2:n-2-this%n_sing)))  then
                 this%u1_lin_comb(i,ii,m) = j
                 ii = ii+1
              endif
             endif
          enddo
        enddo
      enddo
    end subroutine get_u1_lin_comb

    subroutine convert_perm_color_flow(n,iper,jper,iper_red,jper_red)
       implicit none
       integer :: n,i,j,maxval,minval,low,upp,min_x,max_x
       integer,dimension(n) :: iper,jper,iper_red,jper_red
       integer,dimension(2,n) :: dual,dual_red

       minval=100
       maxval=-1
       iper_red=0
       jper_red=0
       dual(1,1:n)=iper(1:n)
       dual(2,1:n)=jper(1:n)
       dual_red(1,:)=iper_red
       dual_red(2,:)=jper_red
       dual(1:2,1)=0
       dual(1:2,n)=0

       do j=1,2
         minval=100
         maxval=-1
         low=1
         upp=n-2
         do
         minval=100
         maxval=-1
         do i=2,n-1
           if (dual(j,i).lt.minval.and.dual(j,i).ne.0) then
               minval=dual(j,i)
               min_x = i
           endif
           if (dual(j,i).gt.maxval.and.dual(j,i).ne.0) then
              maxval=dual(j,i)
              max_x = i
           endif
         enddo
         dual_red(j,min_x)=low
         dual_red(j,max_x)=upp
         low=low+1
         upp=upp-1
         dual(j,min_x)=0
         dual(j,max_x)=0
         if (all(dual(j,:).eq.0)) exit
         enddo
       enddo

       iper_red(1:n) = dual_red(1,1:n)
       jper_red(1:n) = dual_red(2,1:n)
    end subroutine convert_perm_color_flow

  end subroutine init_col2
  
end module amplitude_QCD_mod
