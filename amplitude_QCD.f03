module amplitude_QCD_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: use_symm_cm=.true.
  logical,parameter :: use_cm_dict=.true.
  type :: current
     ! if adding variables here, also update the finalize_current and assign_current subroutines
     integer :: type,bin,n_vert
     integer(kind=16) :: iproc
     integer,dimension(:),allocatable :: vertices,order,spin,ext_type
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8) :: mass,width
   contains
     final :: finalize_current ! custom deallocation of current
  end type current
  type :: interaction
     ! if adding variables here, also update the finalize_interaction and assign_interaction subroutines
     integer :: type
     integer,dimension(2) :: currents
     integer,dimension(:),allocatable :: singlet_mv
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8),dimension(2) :: coupl
   contains
     final :: finalize_interaction ! custom deallocation of interaction
  end type interaction
  ! setting one current (or interaction) equal to another using equal
  ! sign, requires a special operation:
  interface assignment(=)
     module procedure assign_current
     module procedure assign_interaction
  end interface assignment(=)
  type :: amplitude_QCD
     ! if adding variables here, also update the finalize_amplitude_QCD subroutine
     integer :: n_cur,n_vert,imode,nColOrd,max_pp,n_amps,nprocs
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: pp,diff_col_vals
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end, &
          pp_bin_to_i,pp_i_to_bin,col_index,n_col_vals,iproc_start,n_sing,n_qqbar
     integer,dimension(:,:),allocatable :: perm,curr2amp,i_col_i,processes,&
          same_flavour_proc_map,same_flavour_sum
     integer,dimension(:,:,:),allocatable :: spins,row_index
     logical,dimension(:),allocatable :: include_amp,same_flav
   contains
     procedure,public :: init,evaluate,init_col,filter_helicity,write_init_amps_to_file,read_init_amps_from_file
     procedure,private :: filter_dead_trees
     final :: finalize_amplitude_QCD ! custom deallocation of amplitude_QCD
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,n_processes,part,spin,o,pm,read_file)
    use math_functions
    use particles
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    type(physics_model),intent(in) :: pm
    integer,intent(in) :: n,imode,n_processes
    integer,dimension(n,n_processes),intent(in) :: part,o
    integer,dimension(0:3,n),intent(in) :: spin
    integer,dimension(:,:),allocatable :: order
    type(current),dimension(:),allocatable :: current_list_local
    type(interaction),dimension(:),allocatable :: interaction_list_local
    integer :: isize,nc,isplit,n1,n2,ic1,ic2,max_cur,max_vert,max_key,ispin,iproc,jproc
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:,:),allocatable :: key_to_current
    logical :: read_file

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

    call check_input_consistency(part)

    if (this%imode.eq.2) then
       call define_canonical_color_order()
    else
       this%nColOrd=1
    endif

    call set_max_cur()
    call set_max_vert()
    
    if (this%imode.eq.2) then
       call create_current_dict()
       allocate(key_to_current(max_key,maskr(this%nprocs)))
       key_to_current(1:max_key,1:maskr(this%nprocs))=0
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
          do nc=n,1,-1 ! create first the currents that close the amplitude
             if (nc.eq.n) this%n_cur_start(n)=this%n_cur+1
             do iproc=1,this%nprocs
                ! Do not include the same-flavour processes in the creation of
                ! the currents; they will be reconstructed from the
                ! corresponding two different-flavour amplitudes
                if (this%same_flav(iproc)) cycle
                do ispin=1, spin(0,order(nc,iproc))
                   call create_external_current(nc,iproc,spin(ispin,order(nc,iproc)),&
                        this%processes(order(nc,iproc),iproc),order(nc,iproc))
                enddo
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

    call allocate_and_fill_currents_to_amps_map()
    
    call allocate_current_list_and_interaction_list()

    if (this%imode.eq.1) call allocate_and_fill_spins()
    call allocate_and_fill_colour_permutations()
    call allocate_and_fill_momentum_array()

    ! All done. But there could be currents that are not needed. Filter them out
    write (*,*) 'Total number of currents and vertices before filter',this%n_cur,this%n_vert
    call this%filter_dead_trees(n)
    write (*,*) 'Total number of currents and vertices',this%n_cur,this%n_vert

    call deallocate_unneeded()

  contains

    subroutine create_external_current(nc,iproc,ispin,ipart,iorder)
      implicit none
      integer,intent(in) :: nc,ispin,ipart,iorder,iproc
      integer :: i,ic
      do ic=1,this%n_cur
         if (current_list_local(ic)%order(1).ne.iorder) cycle
         if (iorder.le.2 .and. abs(ipart).le.6) then
            if (current_list_local(ic)%type.ne.anti_current(ipart)) cycle
         else
            if (current_list_local(ic)%type.ne.ipart) cycle
         endif
         if (current_list_local(ic)%bin.ne.ibset(0,iorder-1)) cycle
         if (current_list_local(ic)%spin(1).ne.ispin) cycle
         ! existing current.
         current_list_local(ic)%iproc=ibset(current_list_local(ic)%iproc,iproc-1)
         return
      enddo
      ! new external current
      this%n_cur=this%n_cur+1
      allocate(current_list_local(this%n_cur)%order(isize))
      current_list_local(this%n_cur)%order(1)=iorder
      if (iorder.le.2 .and. abs(ipart).le.6) then ! initial quark states
         current_list_local(this%n_cur)%type=anti_current(ipart) ! switch quark <--> anti-quark for initial states
      else
         current_list_local(this%n_cur)%type=ipart
      endif
      current_list_local(this%n_cur)%mass=pm%get_mass(current_list_local(this%n_cur)%type)
      current_list_local(this%n_cur)%width=pm%get_width(current_list_local(this%n_cur)%type)
      allocate(current_list_local(this%n_cur)%ext_type(isize))
      current_list_local(this%n_cur)%ext_type(1)=current_list_local(this%n_cur)%type
      current_list_local(this%n_cur)%bin=ibset(0,iorder-1) ! give binary label
      allocate(current_list_local(this%n_cur)%spin(isize))
      current_list_local(this%n_cur)%spin(1)=ispin
      current_list_local(this%n_cur)%n_vert=0
      current_list_local(this%n_cur)%iproc=ibset(int(0,kind=16),iproc-1)
    end subroutine create_external_current
    
    subroutine allocate_and_fill_currents_to_amps_map()
      ! The 'curr2amp(1:2,iamp)' variable lists which two currents (one of
      ! size n-1 and one of size 1) result in the amplitude 'iamp'. This
      ! subroutine also sets the include_amp(iamp) to .true. for all
      ! amplitudes
      implicit none
      integer :: icur,jcur,iamp,jamp,i,j,iproc,n_diff_flavour_amps
      integer(kind=16) :: proc
      integer,dimension(:,:),allocatable :: curr2amp,same_flavour_sum
      integer,dimension(n) :: iord,jord,ispn,jspn
      integer,dimension(n,this%nprocs) :: procs
      do j=1,this%nprocs
         do i=1,n
            if (order(i,j).le.2) then
               procs(i,j)=anti_current(this%processes(order(i,j),j))
            else
               procs(i,j)=this%processes(order(i,j),j)
            endif
         enddo
      enddo
      allocate(curr2amp(1:2,1:(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*(this%n_cur_end(n)-this%n_cur_start(n)+1)))
      allocate(this%iproc_start(1:this%nprocs+1))
      this%n_amps=0
      do iproc=1,this%nprocs
         this%iproc_start(iproc)=this%n_amps+1
         if (this%same_flav(iproc)) cycle
         do icur=this%n_cur_start(n-1),this%n_cur_end(n-1)
            do jcur=this%n_cur_start(n),this%n_cur_end(n)
               if ( current_list_local(icur)%type .ne. anti_current(current_list_local(jcur)%type) ) cycle
               proc=iand(current_list_local(icur)%iproc,current_list_local(jcur)%iproc)
               if (popcnt(proc).eq.0) then
                  ! combination of icur and jcur does not contribute to any of the processes
                  cycle
               elseif (popcnt(proc).ne.1) then
                  write (*,*) 'A given amplitude should only contribute to one process'
                  write (*,*) proc,icur,jcur
                  write (*,'(a,i40,B64)') 'cur-i',current_list_local(icur)%iproc,current_list_local(icur)%iproc
                  write (*,'(a,i40,B64)') 'cur-j',current_list_local(jcur)%iproc,current_list_local(jcur)%iproc
                  write (*,*) current_list_local(icur)%ext_type(1:n-1),'   , ',current_list_local(jcur)%ext_type(1)
                  stop 1
               elseif (.not. btest(proc,iproc-1)) then
                  ! one process, but it is not equal to process 'iproc'
                  cycle
               endif
               this%n_amps=this%n_amps+1
               curr2amp(1,this%n_amps)=icur
               curr2amp(2,this%n_amps)=jcur
            enddo
         enddo
      enddo
      n_diff_flavour_amps=this%n_amps
      allocate(same_flavour_sum(2*this%n_amps,2))
      same_flavour_sum=-1
      do iproc=1,this%nprocs
         if (.not. this%same_flav(iproc)) cycle
         this%iproc_start(iproc)=this%n_amps+1
         do iamp=this%iproc_start(this%same_flavour_proc_map(iproc,1)),this%iproc_start(this%same_flavour_proc_map(iproc,1)+1)-1
            iord=[current_list_local(curr2amp(1,iamp))%order(1:n-1),current_list_local(curr2amp(2,iamp))%order(1)]
            ispn=[current_list_local(curr2amp(1,iamp))%spin(1:n-1) ,current_list_local(curr2amp(2,iamp))%spin(1) ]
            do jamp=this%iproc_start(this%same_flavour_proc_map(iproc,2)),this%iproc_start(this%same_flavour_proc_map(iproc,2)+1)-1
               ! if they have the same colour order and spins, they need to be added together
               jord=[current_list_local(curr2amp(1,jamp))%order(1:n-1),current_list_local(curr2amp(2,jamp))%order(1)]
               jspn=[current_list_local(curr2amp(1,jamp))%spin(1:n-1) ,current_list_local(curr2amp(2,jamp))%spin(1) ]
               ! check that the two quarks are in similar order (the anti-quarks might be different order). 
               if (iord(1).ne.jord(1)) then
                  ! different order of the quarks, switch the two colour strings:
                  do i=1,n
                     if (jord(i).eq.iord(1)) then
                        jord(1:n)=[jord(i:n),jord(1:i-1)]
                        jspn(1:n)=[jspn(i:n),jspn(1:i-1)]
                        exit
                     endif
                  enddo
               endif
               if (any(jord(1:n).ne.iord(1:n)) .or. any(jspn(1:n).ne.ispn(1:n))) cycle
               this%n_amps=this%n_amps+1
               if ( abs(this%processes(iord(1),this%same_flavour_proc_map(iproc,1))).eq. &
                    abs(this%processes(iord(n),this%same_flavour_proc_map(iproc,1)))) then
                  same_flavour_sum(this%n_amps,1)=iamp ! in iamp, the different-flavour quarks are connected
                  same_flavour_sum(this%n_amps,2)=jamp ! in jamp, the different-flavour quarks are not connected
               else
                  same_flavour_sum(this%n_amps,1)=jamp ! in jamp, the different-flavour quarks are connected    
                  same_flavour_sum(this%n_amps,2)=iamp ! in iamp, the different-flavour quarks are not connected
               endif
               exit
            enddo
            if (jamp.gt.this%n_amps) then
               write (*,*) 'permutation not found',jamp,this%n_amps
               stop 1
            endif
         enddo
         if (any(same_flavour_sum(this%iproc_start(iproc):this%n_amps,1:2).eq.-1)) then
            write (*,*) same_flavour_sum(this%iproc_start(iproc):this%n_amps,1)
            write (*,*) same_flavour_sum(this%iproc_start(iproc):this%n_amps,2)
            write (*,*) 'not all permutations mapped'
            stop 1
         endif
      enddo

      if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. this%imode.eq.2) then
         allocate(this%curr2amp(1:2,1:2*this%n_amps))
         this%curr2amp(1:2,1:this%n_amps)=curr2amp(1:2,1:this%n_amps)
         this%curr2amp(1:2,this%n_amps+1:2*this%n_amps)=curr2amp(1:2,1:this%n_amps)
         this%n_amps=this%n_amps*2
      else
         allocate(this%curr2amp(1:2,1:n_diff_flavour_amps))
         this%curr2amp(1:2,1:n_diff_flavour_amps)=curr2amp(1:2,1:n_diff_flavour_amps)
      endif
      if (this%imode.eq.3 .and. this%n_amps.ne.1) then
         write (*,*) 'For this%imode==3, there should only be one amplitude',this%n_amps
         write (*,*) this%n_cur_start
         write (*,*) this%n_cur_end
         stop 1
      endif

      allocate(this%same_flavour_sum(this%n_amps,2))
      this%same_flavour_sum(1:this%n_amps,1:2)=same_flavour_sum(1:this%n_amps,1:2)
      this%iproc_start(this%nprocs+1)=this%n_amps+1

      allocate(this%include_amp(1:this%n_amps))
      this%include_amp(:)=.true.
    end subroutine allocate_and_fill_currents_to_amps_map

    subroutine allocate_and_fill_spins()
      implicit none
      integer :: iamp,i,iproc
      allocate(this%spins(n,1,1:this%n_amps))
      do iproc=1,this%nprocs
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            if (.not.this%same_flav(iproc)) then
               do i=1,n
                  if (i.lt.n) then
                     this%spins(this%current_list(this%curr2amp(1,iamp))%order(i),1,iamp)= &
                          this%current_list(this%curr2amp(1,iamp))%spin(i)
                  elseif (i.eq.n) then
                     this%spins(this%current_list(this%curr2amp(2,iamp))%order(1),1,iamp)= &
                          this%current_list(this%curr2amp(2,iamp))%spin(1)
                  endif
               enddo
            else
               if (any(this%spins(1:n,1,this%same_flavour_sum(iamp,1)).ne.this%spins(1:n,1,this%same_flavour_sum(iamp,2)))) then
                  write (*,*) 'ERROR in spin mapping'
                  stop 1
               endif
               this%spins(1:n,1,iamp)=this%spins(1:n,1,this%same_flavour_sum(iamp,1))
            endif
         enddo
      enddo
    end subroutine allocate_and_fill_spins
    
    subroutine allocate_and_fill_colour_permutations()
      implicit none
      integer :: iamp,i,iproc,iamp_to_compare
      ! allocate and fill the colour orders in 'this%perm'. These are simply
      ! the orders of the elements in the 'this%current_list' (with size n-1)
      ! together with the final element). Exception: when there are colour
      ! singlets, they will not be part of the this%perm (while they are part
      ! of the elements in the this%current_list.
      allocate(this%perm(1:n-this%n_sing(1),1:this%n_amps))
      if (is_singlet(this%current_list(this%n_cur_start(n))%type)) then
         write (*,*) 'Final current (that closes the amplitude) cannot be a colour singlet'
         write (*,*) this%current_list(this%n_cur_start(n))%type
         stop 1
      endif
      do iproc=1,this%nprocs
         if (this%imode.ne.2) then
            iamp_to_compare=this%iproc_start(iproc)
         else
            iamp_to_compare=1
         endif
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            if (.not.this%same_flav(iproc)) then
               this%perm(1:n-this%n_sing(1),iamp)=[this%current_list(this%curr2amp(1,iamp))%order(1:n-1-this%n_sing(1)),&
                                                   this%current_list(this%curr2amp(2,iamp))%order(1)]
               if (use_symmetry .and. this%n_qqbar(iproc).eq.0 .and. iamp.gt.this%n_amps/2) then
                  this%perm(1:n-this%n_sing(1),iamp)=[this%current_list(this%curr2amp(1,iamp))%order(n-1-this%n_sing(1):1:-1),&
                                                      this%current_list(this%curr2amp(2,iamp))%order(1)]
               endif
               if (this%n_qqbar(iproc).eq.2) then
                  ! make sure that orders of the quarks is fixed among all
                  ! perm's. (The order of the anti-quarks may vary). Without
                  ! this, the computation of the colour factor will be incorrect.
                  if (iamp.eq.iamp_to_compare) cycle
                  if (this%perm(1,iamp).ne.this%perm(1,iamp_to_compare)) then
                     ! different order, switch the two colour strings:
                     do i=1,n-this%n_sing(1)
                        if (this%perm(i,iamp).eq.this%perm(1,iamp_to_compare)) then
                           this%perm(1:n-this%n_sing(1),iamp)=[this%perm(i:n-this%n_sing(1),iamp),this%perm(1:i-1,iamp)]
                           exit
                        endif
                     enddo
                  endif
               endif
            else
               this%perm(1:n-this%n_sing(1),iamp)=this%perm(1:n-this%n_sing(1),this%same_flavour_sum(iamp,1))
            endif
         enddo
      enddo
    end subroutine allocate_and_fill_colour_permutations

    subroutine allocate_and_fill_momentum_array()
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
    end subroutine allocate_and_fill_momentum_array

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
      if (this%n_cur_start(n-1).gt.this%n_cur_end(n-1)) then
         write (*,*) 'ERROR: no valid matrix elements found: check your process definition'
         stop 1
      endif
    end subroutine simple_consistency_checks

    subroutine define_canonical_color_order()
      ! canonical order: (q,glu,glu,glu,qbar,q,singlet,singlet,qbar)
      use math_functions
      implicit none
      integer :: i,nq,naq,nglu,nsing,iq,iaq,iglu,ising
      do iproc=1,this%nprocs
         nq=0; naq=0 ; nglu=0 ; nsing=0
         do i=1,n
            if (is_gluon(this%processes(i,iproc))) nglu=nglu+1
            if (is_quark_from_order(i,iproc)) nq=nq+1
            if (is_antiquark_from_order(i,iproc)) naq=naq+1
            if (is_singlet(this%processes(i,iproc))) nsing=nsing+1
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
         order(1:n,iproc)= 0
         do i=1,n
            if (is_gluon(this%processes(i,iproc))) then
               iglu=iglu+1
               if (nq.ge.1) then
                  order(iglu+1,iproc)=i
               else
                  order(iglu,iproc)=i
               endif
            elseif(is_quark_from_order(i,iproc)) then
               iq=iq+1
               if (iq.eq.1) then
                  order(1,iproc)=i
               else
                  order(n-1-nsing,iproc)=i
               endif
            elseif (is_antiquark_from_order(i,iproc)) then
               iaq=iaq+1
               if (iaq.eq.1) then
                  order(n,iproc)=i
               else
                  order(n-2-nsing,iproc)=i
               endif
            elseif (is_singlet(this%processes(i,iproc))) then
               ising=ising+1
               if(nq.ne.0) then
                  order(n-ising,iproc)=i
               else
                  order(nglu+ising,iproc)=i
               endif
            endif
         enddo
      enddo
      if (nq.eq.0) then
         this%nColOrd=factorial(nglu-1)
      elseif (nq.eq.1) then
         this%nColOrd=factorial(nglu)
      elseif (nq.eq.2) then
         this%nColOrd=factorial(nglu)*(nglu+1)*2
      else
         write (*,*) 'Number of colour orders unknown',nq
         stop 1
      endif
    end subroutine define_canonical_color_order

    integer function number_of_quark_lines(process,is_same_flavour_process)
      implicit none
      integer,dimension(n),intent(in) :: process
      logical,intent(out) :: is_same_flavour_process
      integer :: i,iflav
      is_same_flavour_process=.true.
      number_of_quark_lines=0
      iflav=0
      do i=1,n
         if (is_quark(process(i)) .or. is_antiquark(process(i))) then
            number_of_quark_lines=number_of_quark_lines+1
            if (iflav.eq.0) iflav=abs(process(i))
            if (abs(process(i)).ne.iflav) is_same_flavour_process=.false.
         endif
      enddo
      number_of_quark_lines=number_of_quark_lines/2
      if (number_of_quark_lines.lt.2) is_same_flavour_process=.false.
    end function number_of_quark_lines
    
    subroutine check_input_consistency(part)
      implicit none
      integer :: i,j,k,iflav,iproc,jproc,idum,ichan
      integer,dimension(n,n_processes),intent(in) :: part
      integer,dimension(n,2) :: part_sf
      integer,dimension(n) :: jord
      logical :: sf
      if (this%imode.eq.2) then
         if (n_processes.ne.1) then
            write (*,*) 'There should only be one process when doing imode=2'
            stop 1
         endif
         idum=number_of_quark_lines(part(1,1),sf)
         if (sf) then
            ! add the two different-flavour processes that make up the same-flavour process
            this%nprocs=3
            allocate(this%processes(n,this%nprocs))
            call define_symm_2qq(part(1,1),part_sf,1)
            this%processes(1:n,1)=part_sf(1:n,1)
            call define_symm_2qq(part(1,1),part_sf,2)
            this%processes(1:n,2)=part_sf(1:n,1)
            this%processes(1:n,3)=part(1:n,1)
            allocate(order(1:n,this%nprocs))
            order(1:n,1)=o(1:n,1)
            order(1:n,2)=o(1:n,1)
            order(1:n,3)=o(1:n,1)
         else
            this%nprocs=n_processes
            allocate(this%processes(n,this%nprocs))
            this%processes(1:n,1:this%nprocs)=part(1:n,1:this%nprocs)
            allocate(order(1:n,this%nprocs))
            order(1:n,1:this%nprocs)=o(1:n,1:this%nprocs)
         endif
      else
         idum=number_of_quark_lines(part(1,1),sf)
         if (sf .and. .not.read_file) then
            ! add the two different-flavour processes that make up the same-flavour process
            this%nprocs=3
            allocate(this%processes(n,this%nprocs))
            call define_symm_2qq(part(1,1),part_sf,1)
            this%processes(1:n,1)=part_sf(1:n,1)
            call define_symm_2qq(part(1,1),part_sf,2)
            this%processes(1:n,2)=part_sf(1:n,1)
            this%processes(1:n,3)=part(1:n,1)
            allocate(order(1:n,this%nprocs))
            order(1:n,1)=o(1:n,1)
            order(1:n,2)=o(1:n,1)
            order(1:n,3)=o(1:n,1)
         else
            this%nprocs=n_processes
            allocate(this%processes(n,this%nprocs))
            this%processes(1:n,1:this%nprocs)=part(1:n,1:this%nprocs)
            allocate(order(1:n,this%nprocs))
            order(1:n,1:this%nprocs)=o(1:n,1:this%nprocs)
         endif
      endif
      
      allocate(this%n_sing(1:this%nprocs))
      allocate(this%n_qqbar(1:this%nprocs))
      allocate(this%same_flav(1:this%nprocs))
      allocate(this%same_flavour_proc_map(1:this%nprocs,2))
      do iproc=1,this%nprocs
         this%n_qqbar(iproc)=number_of_quark_lines(this%processes(1,iproc),this%same_flav(iproc))
         this%same_flavour_proc_map(iproc,1:2)=0
         this%n_sing(iproc)=0
         do i=1,n
            if (is_singlet(this%processes(i,iproc))) this%n_sing(iproc)=this%n_sing(iproc)+1
         enddo
         if (iproc.gt.1) then
            if (this%n_qqbar(iproc-1).gt.this%n_qqbar(iproc)) then
               write (*,*) 'ERROR: processes not correctly ordered in the list.'
               write (*,*) 'Need to be in increasing number of quark lines.'
               stop 1
            endif
            if (this%same_flav(iproc-1) .and. (.not.this%same_flav(iproc))) then
               write (*,*) 'ERROR: processes not correctly ordered in the list.'
               write (*,*) 'Need first different-flavour and then same-flavour processes.'
               stop 1
            endif
         endif
         if (this%n_qqbar(iproc).gt.3) then
            write (*,*) 'ERROR: code only working for 0, 1 or 2 qqbar pairs',this%n_qqbar(iproc),iproc
            write (*,*) this%processes(1:n,iproc)
            stop 1
         endif
         if (imode.eq.1.or.imode.eq.3) then
            if (any(order(:,iproc).gt.n) .or. any(order(:,iproc).lt.1)) then
               write (*,*) 'ERROR: inconsistent colour order. An element is too large or too small',order(1:n,iproc),iproc
               stop 1
            endif
            do i=1,n-1
               do j=i+1,n
                  if (order(i,iproc).eq.order(j,iproc)) then
                     write (*,*) 'ERROR: inconsistent colour order. An element appears twice',order(1:n,iproc),iproc
                     stop 1
                  endif
               enddo
            enddo
            if (this%n_qqbar(iproc).gt.0) then
               if (.not.(is_quark_from_order(order(1,iproc),iproc))) then
                  write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
                  write (*,*) iproc
                  write (*,*) order(1:n,iproc)
                  write (*,*) this%processes(1:n,iproc)
                  stop 1
               endif
               if (.not.(is_antiquark_from_order(order(n,iproc),iproc))) then
                  write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
                  write (*,*) iproc
                  write (*,*) order(1:n,iproc)
                  write (*,*) this%processes(1:n,iproc)
                  stop 1
               endif
            endif
            if (this%n_qqbar(iproc).ge.2) then
               do i=2,n-1
                  if (is_antiquark_from_order(order(i,iproc),iproc)) then
                     ! next should be a quark
                     if (.not.(is_quark_from_order(order(i+1,iproc),iproc))) then
                        write (*,*) 'ERROR: in the colour order, after an initial state quark should come a final state quark'
                        write (*,*) iproc
                        write (*,*) order(1:n,iproc)
                        write (*,*) this%processes(1:n,iproc)
                        stop 1
                     endif
                  endif
               enddo
               if (use_real_gluons) then
                  write (*,*) 'ERROR: cannot use real gluons with two quark lines around'
                  stop 1
               endif
            endif
         endif
         if (this%same_flav(iproc)) then
            ! check that the corresponding different-flavour processes are included.
            if (this%n_qqbar(iproc).lt.2) then
               write (*,*) 'ERROR: Same-flavour process, but there are less than two quark lines'
               write (*,*) iproc,':',this%processes(1:n,iproc)
               stop 1
            endif
            do ichan=1,2
               call define_symm_2qq(this%processes(1,iproc),part_sf,ichan)
               do jproc=1,iproc-1
                  if (this%n_qqbar(jproc).ne.this%n_qqbar(iproc) .or. this%same_flav(jproc)) cycle
                  jord(1:n)=order(1:n,jproc)
                  if (jord(1).ne.order(1,iproc)) then
                     ! different order of the quarks, switch the two colour strings:
                     do i=1,n
                        if (jord(i).eq.order(1,iproc)) then
                           jord(1:n)=[jord(i:n),jord(1:i-1)]
                           exit
                        endif
                     enddo
                  endif
                  if (any(jord(1:n).ne.order(1:n,iproc))) cycle
                  if (all(this%processes(1:n,jproc).eq.part_sf(1:n,1))) then
                     this%same_flavour_proc_map(iproc,ichan)=jproc
                  elseif (all(this%processes(1:n,jproc).eq.part_sf(1:n,2))) then
                     this%same_flavour_proc_map(iproc,ichan)=jproc
                  endif
               enddo
            enddo
            if (any(this%same_flavour_proc_map(iproc,1:2).eq.0)) then
               write (*,*) 'Same flavour process not found',iproc
               write (*,*) part_sf(1:n,1)
               write (*,*) part_sf(1:n,2)
               stop 1
            endif
         endif
      enddo
    end subroutine check_input_consistency

    subroutine define_symm_2qq(part_in,part_out,chan)
      implicit none
      integer :: chan
      integer,dimension(n),intent(in) :: part_in
      integer,dimension(n,2),intent(out) :: part_out
      integer :: i,iq,ia,i_same,n_sing,add_or_subtract
      integer,dimension(2,2) :: connection
      n_sing=0
      do i=1,n
         if (is_singlet(part_in(i))) n_sing=n_sing+1
      enddo
      if (n_sing.eq.0) then
         i_same=1
      else
         i_same=2
      endif
      do i=1,n
         if (is_quark(part_in(i)).or.is_antiquark(part_in(i))) then
            if (abs(part_in(i)).gt.2) then
               add_or_subtract=-1
            else
               add_or_subtract=1
            endif
            exit
         endif
      enddo
      part_out(1:n,1)=part_in(1:n)
      part_out(1:n,2)=part_in(1:n)
      iq=0
      ia=0
      do i=1,n
         if (i.le.2 .and. is_quark(part_in(i))) then
            ia=ia+1
            connection(2,ia)=i
         elseif(i.le.2 .and. is_antiquark(part_in(i))) then
            iq=iq+1
            connection(1,iq)=i
         elseif (i.gt.2 .and. is_quark(part_in(i))) then
            iq=iq+1
            connection(1,iq)=i
         elseif (i.gt.2 .and. is_antiquark(part_in(i))) then
            ia=ia+1
            connection(2,ia)=i
         endif
      enddo
      if (chan.eq.1) then
         ! change the 2nd quark and an anti-quark in the process
         part_out(connection(1,2),1)=sign(abs(part_in(connection(1,2)))+add_or_subtract*i_same,part_in(connection(1,2)))
         part_out(connection(2,2),1)=sign(abs(part_in(connection(2,2)))+add_or_subtract*i_same,part_in(connection(2,2)))
         ! change the 1st quark and an anti-quark in the process
         part_out(connection(1,1),2)=sign(abs(part_in(connection(1,1)))+add_or_subtract*i_same,part_in(connection(1,1)))
         part_out(connection(2,1),2)=sign(abs(part_in(connection(2,1)))+add_or_subtract*i_same,part_in(connection(2,1)))

      elseif(chan.eq.2) then
!!$         if (abs(part_in(connection(1,1))).lt.4) then
            ! change the mixed quark and an anti-quark in the process; leave the
            ! first (anti-)quark unchanged.
               part_out(connection(1,2),1)=sign(abs(part_in(connection(1,2)))+&
                       add_or_subtract*i_same,part_in(connection(1,2)))
               part_out(connection(2,1),1)=sign(abs(part_in(connection(2,1)))+&
                       add_or_subtract*i_same,part_in(connection(2,1)))
               part_out(connection(1,1),2)=sign(abs(part_in(connection(1,1)))+&
                       add_or_subtract*i_same,part_in(connection(1,1)))
               part_out(connection(2,2),2)=sign(abs(part_in(connection(2,2)))+&
                       add_or_subtract*i_same,part_in(connection(2,2)))
!!$         else
!!$            ! change the mixed quark and an anti-quark in the process; leave the
!!$            ! second (anti-)quark unchanged.
!!$               part_out(connection(1,1),1)=sign(abs(part_in(connection(1,1)))+&
!!$                       add_or_subtract*i_same,part_in(connection(1,1)))
!!$               part_out(connection(2,2),1)=sign(abs(part_in(connection(2,2)))+&
!!$                       add_or_subtract*i_same,part_in(connection(2,2)))
!!$               part_out(connection(1,2),2)=sign(abs(part_in(connection(1,2)))+&
!!$                       add_or_subtract*i_same,part_in(connection(1,2)))
!!$               part_out(connection(2,1),2)=sign(abs(part_in(connection(2,1)))+&
!!$                       add_or_subtract*i_same,part_in(connection(2,1)))
!!$         endif
      endif
    end subroutine define_symm_2qq

    
    subroutine set_max_cur()
      ! rough upper bound for the maximum number of currents
      implicit none
      max_cur=1024
    end subroutine set_max_cur

    subroutine set_max_vert()
      ! rough upper bound on the maximum number of interactions
      implicit none
      max_vert=1024
    end subroutine set_max_vert

    subroutine increase_max_cur()
      implicit none
      integer :: new_max_cur,ic
      type(current),dimension(:),allocatable :: tmp
      new_max_cur=2*max_cur
      allocate(tmp(new_max_cur))
      do ic=1,max_cur
         allocate(tmp(ic)%vertices(size(current_list_local(ic)%vertices)))
         allocate(tmp(ic)%vertex_sign(size(current_list_local(ic)%vertex_sign)))
         tmp(ic)=current_list_local(ic)
      enddo
      do ic=1,max_cur
         call finalize_current(current_list_local(ic))
      enddo
      deallocate(current_list_local)
      allocate(current_list_local(new_max_cur))
      do ic=1,max_cur
         allocate(current_list_local(ic)%vertices(size(tmp(ic)%vertices)))
         allocate(current_list_local(ic)%vertex_sign(size(tmp(ic)%vertex_sign)))
         current_list_local(ic)=tmp(ic)
      enddo
      do ic=1,max_cur
         call finalize_current(tmp(ic))
      enddo
      deallocate(tmp)
      max_cur=new_max_cur
    end subroutine increase_max_cur
    
    subroutine increase_max_vert()
      implicit none
      integer :: new_max_vert,iv
      type(interaction),dimension(:),allocatable :: tmp
      new_max_vert=2*max_vert
      allocate(tmp(new_max_vert))
      ! copy old list into tmp
      do iv=1,max_vert
         if (allocated(interaction_list_local(iv)%singlet_mv)) &
              allocate(tmp(iv)%singlet_mv(0:size(interaction_list_local(iv)%singlet_mv)-1))
         tmp(iv)=interaction_list_local(iv)
      enddo
      ! empty old list
      do iv=1,max_vert
         call finalize_interaction(interaction_list_local(iv))
      enddo
      deallocate(interaction_list_local)
      ! allocate new list
      allocate(interaction_list_local(new_max_vert))
      ! copy tmp into new list
      do iv=1,max_vert
         if (allocated(tmp(iv)%singlet_mv)) &
              allocate(interaction_list_local(iv)%singlet_mv(0:size(tmp(iv)%singlet_mv)-1))
         interaction_list_local(iv)=tmp(iv)
      enddo
      ! empty tmp
      do iv=1,max_vert
         call finalize_interaction(tmp(iv))
      enddo
      deallocate(tmp)
      max_vert=new_max_vert
    end subroutine increase_max_vert
      
    subroutine allocate_current_list_and_interaction_list()
      ! allocate the minimum memory needed for the current_list and
      ! interaction_list to be able to perform the evaluate() procedure.
      implicit none
      integer :: isize,ic,iv
      allocate(this%current_list(1:this%n_cur))
      do isize=1,n-1
         do ic=this%n_cur_start(isize),this%n_cur_end(isize)
            this%current_list(ic)=current_list_local(ic) ! use non-custom 'assignment'
         enddo
      enddo
      ! make sure to deallocate currents consistently
      do ic=1,size(current_list_local)
         call finalize_current(current_list_local(ic))
      enddo
      deallocate(current_list_local)
      allocate(this%interaction_list(1:this%n_vert))
      do iv=1,this%n_vert
         this%interaction_list(iv)=interaction_list_local(iv) ! use non-custom 'assignment'
      enddo
      ! make sure to deallocate interactions consistently
      do iv=1,size(interaction_list_local)
         call finalize_interaction(interaction_list_local(iv))
      enddo
      deallocate(interaction_list_local)
    end subroutine allocate_current_list_and_interaction_list

    subroutine add_if_allowed_threevertex()
      ! check if we should consider the current combination, and if
      ! so, and the corresponding vertices to the list.
      implicit none
      integer :: i
      if (.not.valid_current_combination())  then
         return
      endif
      do i=1,pm%nint
         if ( current_list_local(ic1)%type.eq.pm%vertex_list(i)%particles(1) .and. &
              current_list_local(ic2)%type.eq.pm%vertex_list(i)%particles(2) ) then
            call add_vertex(pm%vertex_list(i)%type, &
                            pm%vertex_list(i)%particles(3), &
                            pm%vertex_list(i)%coupl)
         endif
      enddo
    end subroutine add_if_allowed_threevertex

    
    logical function valid_current_combination()
      ! Checks to see if the combination of currents ic1 and ic2 could be a
      ! valid combination. Checks to perform:
      ! 0. All particles must be different in the two currents & final
      !    particle should never be part of the combined currents (it will be
      !    used to close the currents), and both currents must be able to
      !    contribute to the same process.
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
      logical :: gluon_current,colour_singlet1,colour_singlet2,found_quark,found_antiquark,not_valid
      integer,dimension(isize) :: ip,et
      valid_current_combination=.false.
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(current_list_local(ic1)%bin,current_list_local(ic2)%bin)).ne.isize) return
      ! final particle should never be part of any combined currents: it will
      ! be used to close the amplitude instead
      if (n1.eq.1) then
         if (all(current_list_local(ic1)%order(n1).eq.order(n,1:this%nprocs))) then
            return
          endif
      endif
      if (n2.eq.1) then
         if (all(current_list_local(ic2)%order(n2).eq.order(n,1:this%nprocs))) then
            return
         endif
      endif
      ! check that both currents can contribute to the same process
      if (iand(current_list_local(ic1)%iproc,current_list_local(ic2)%iproc).eq.0) return
      ! Check for colour singlets:
      colour_singlet1=all_singlet_current(current_list_local(ic1),n1)
      colour_singlet2=all_singlet_current(current_list_local(ic2),n2)
      ! If the first current is a singlet and the second is not, it is not a valid order
      if (colour_singlet1 .and. (.not.colour_singlet2)) return
      ! If both currents are colour singlets, only consider one of the two. 
      if (colour_singlet1 .and. colour_singlet2) then
         if (maxval(current_list_local(ic1)%order(1:n1)).ge.maxval(current_list_local(ic2)%order(1:n2))) then
            return
         else
            valid_current_combination=.true.
            return ! no need to check further: below are only checks about the colours
         endif
      endif
      
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         ! check that current combination is compatible with the input colour
         ! order. First, find where the singlets are, since they do not matter
         ! for the colour order
         do i=1,n1
            if (is_singlet(current_list_local(ic1)%ext_type(i))) exit
         enddo
         nc1=i-1
         do i=1,n2
            if (is_singlet(current_list_local(ic2)%ext_type(i))) exit
         enddo
         nc2=i-1
         ip(1:nc1+nc2)=[current_list_local(ic1)%order(1:nc1),current_list_local(ic2)%order(1:nc2)]
         ! do they actual checking:
         not_valid=.true.
         do_iproc: do iproc=1,this%nprocs
            ! check that both currents contribute to the iproc process:
            if (.not. btest(iand(current_list_local(ic1)%iproc,current_list_local(ic2)%iproc),iproc-1)) cycle
            ! check that the final particle is not part of the combined
            ! current (it will be used to close the amplitude instead):
            if (btest(current_list_local(ic1)%bin+current_list_local(ic2)%bin,order(n,iproc)-1)) cycle
            ! Check if they are compatible with the colour order of the iproc:
            do_j: do j=1,n
               if (order(j,iproc).eq.ip(1)) then
                  do i=2,nc1+nc2
                     if (j-1+i.gt.n) exit do_j
                     if (order(j-1+i,iproc).ne.ip(i)) exit do_j
                  enddo
                  not_valid=.false. ! it's compatible with the input colour order of iproc
                  exit do_iproc
               endif
            enddo do_j
         enddo do_iproc
         if (not_valid) return ! not compatible with any of the iprocs
      endif

      ! If using symmetry and the current is a combination of all external
      ! gluons, take only one of the two possible orders
      gluon_current=all_gluon_current(current_list_local(ic1),n1).and.all_gluon_current(current_list_local(ic2),n2)
      if (use_symmetry .and. this%imode.eq.2 .and. gluon_current) then
         if (maxval(current_list_local(ic1)%order(1:n1)).ge.maxval(current_list_local(ic2)%order(1:n2))) return
      endif
      if (.not. gluon_current) then
         ! for two quark lines. Should be of the form "q g..g qbar q g..g qbar" or any subset thereof. 
         et(1:isize)=[current_list_local(ic1)%ext_type(1:n1),current_list_local(ic2)%ext_type(1:n2)]
         found_quark=.false.
         found_antiquark=.false.
         do i=1,isize
            if (is_quark(et(i))) then
               ! found a quark.
               if (found_quark) then
                  ! no anti-quark between two quarks
                  return
               endif
               found_antiquark=.false.
               found_quark=.true.
            elseif (is_antiquark(et(i))) then
               ! found an anti-quark
               if (found_antiquark) then
                  ! no quark between two anti-quarks
                  return
               endif
               ! next one must be a quark:
               j=i
               do while (j.lt.isize)
                  if (is_singlet(et(j+1))) then
                     j=j+1
                  elseif (.not.(is_quark(et(j+1)))) then
                     return
                  else
                     exit
                  endif
               enddo
               found_quark=.false.
               found_antiquark=.true.
            endif
         enddo
         ! if the current is of length n-1, the first should be a quark
         if (isize.eq.n-1 .and. .not.is_quark(et(1))) return
      endif
      ! Got all the way to the end. This must be a valid current combination
      valid_current_combination=.true.
    end function valid_current_combination
    
    subroutine add_vertex(itype,ctype,coupl)
      implicit none
      integer :: itype,ctype,ic
      real(kind=8),dimension(2) :: coupl
      if (isize.eq.n-1) then
         do ic=this%n_cur_start(n),this%n_cur_end(n)
            if (ctype.eq.anti_current(current_list_local(ic)%type)) exit
         enddo
         if (ic.eq.this%n_cur_end(n)+1) return ! dead tree. Filter already here
      endif
      this%n_vert=this%n_vert+1
      if (this%n_vert.gt.max_vert) call increase_max_vert()
      interaction_list_local(this%n_vert)%type=itype
      interaction_list_local(this%n_vert)%currents(1)=ic1
      interaction_list_local(this%n_vert)%currents(2)=ic2
      interaction_list_local(this%n_vert)%coupl=coupl
      allocate(interaction_list_local(this%n_vert)%singlet_mv(0:isize))
      call add_all_currents(ctype)
    end subroutine add_vertex

    function combine_lists(current,singlet_mv)
      ! just concatenate the two currents, except if there is a colour
      ! singlet. Move the label of the singlet to the end of the combined
      ! order.
      implicit none
      integer,dimension(isize) :: combine_lists,current
      integer,dimension(0:isize),intent(in) :: singlet_mv
      integer :: imv
      combine_lists(1:isize)=current(1:isize)
      do imv=1,singlet_mv(0)
         combine_lists(1:isize)=[combine_lists(1:singlet_mv(imv)-1), &
                                 combine_lists(singlet_mv(imv)+1:isize),combine_lists(singlet_mv(imv))]
      enddo
    end function combine_lists

    type(current) function combine_currents(ic1,ic2,ctype,singlet_mv,invert)
      ! combine the currents corresponding to ic1 and ic2 into a new current
      ! of type 'ctype'. This also sets up 'singlet_mv' that determines how to
      ! move the colour singlets to the correct position. If the first (or
      ! second) bit of 'invert' is set to 1, the colour order of ic1 (or ic2)
      ! is reversed before the two currents are combined.
      implicit none
      integer,intent(in) :: ic1,ic2,ctype
      integer,intent(in) :: invert
      integer,dimension(0:isize),intent(out) :: singlet_mv
      integer :: i,n1,n2,ipos,mv12,nc1,nc2,ns1,ns2
      integer,dimension(isize) :: ord
      integer,dimension(:),allocatable :: ord1,spin1,et1,ord2,spin2,et2
      allocate(combine_currents%order(1:isize))
      allocate(combine_currents%spin(1:isize))
      allocate(combine_currents%ext_type(1:isize))
      combine_currents%type=ctype
      combine_currents%bin=current_list_local(ic1)%bin+current_list_local(ic2)%bin
      combine_currents%iproc=iand(current_list_local(ic1)%iproc,current_list_local(ic2)%iproc)
      n1=popcnt(current_list_local(ic1)%bin)
      n2=popcnt(current_list_local(ic2)%bin)
      allocate(ord1(n1))
      allocate(spin1(n1))
      allocate(et1(n1))
      allocate(ord2(n2))
      allocate(spin2(n2))
      allocate(et2(n2))
      if (btest(invert,0)) then
         ord1(1:n1)=current_list_local(ic1)%order(n1:1:-1)
         spin1(1:n1)=current_list_local(ic1)%spin(n1:1:-1)
         et1(1:n1)=current_list_local(ic1)%ext_type(n1:1:-1)
      else
         ord1(1:n1)=current_list_local(ic1)%order(1:n1)
         spin1(1:n1)=current_list_local(ic1)%spin(1:n1)
         et1(1:n1)=current_list_local(ic1)%ext_type(1:n1)
      endif
      if (btest(invert,1)) then
         ord2=current_list_local(ic2)%order(n2:1:-1)
         spin2=current_list_local(ic2)%spin(n2:1:-1)
         et2=current_list_local(ic2)%ext_type(n2:1:-1)
      else
         ord2=current_list_local(ic2)%order(1:n2)
         spin2=current_list_local(ic2)%spin(1:n2)
         et2=current_list_local(ic2)%ext_type(1:n2)
      endif

      ! nc1 and nc2 lists where the coloured particles end in currents ic1 and
      ! ic2, respectively:
      do i=1,n1
         if (is_singlet(et1(i))) exit
      enddo
      nc1=i-1
      do i=1,n2
         if (is_singlet(et2(i))) exit
      enddo
      nc2=i-1

      ! The order of the coloured particles can be concatinated:
      ord(1:nc1+nc2)=[ord1(1:nc1),ord2(1:nc2)]
      ! Setup the singlet_mv and put the colour singlets in the right order in
      ! the combined ord():
      if (nc1.eq.n1) then
         ! No colour singlets or all colour singlets are in ic2
         singlet_mv(0)=0
         ord(nc1+nc2+1:n1+n2)=ord2(nc2+1:n2)
      elseif(nc2.eq.n2) then
         ! Some colour singlets in ic1, but no in ic2
         singlet_mv(0)=n1-nc1
         singlet_mv(1:singlet_mv(0))=nc1+1
         ord(nc1+nc2+1:n1+n2)=ord1(nc1+1:n1)
      else
         ! Some colour singlets in both ic1 and ic2
         singlet_mv(0)=0
         ns1=nc1+1
         ns2=nc2+1
         if (nc2.eq.0) then
            ! Special case: no coloured particles in ic2
            if (ord1(n1).lt.ord2(1)) then
               ! nothing to move
               ord(1:n1+n2)=[ord1(1:n1),ord2(1:n2)]
               combine_currents%order(1:isize)=ord(1:isize)
               combine_currents%spin(1:isize)=combine_lists([spin1(1:n1),spin2(1:n2)],singlet_mv)
               combine_currents%ext_type(1:isize)=combine_lists([et1(1:n1),et2(1:n2)],singlet_mv)
               return
            endif
            do while (ord1(ns1).lt.ord2(1))
               ns1=ns1+1
            enddo
            ord(nc1+1:ns1)=ord1(nc1+1:ns1)
         endif
         do while (ord2(ns2).lt.ord1(ns1))
            ns2=ns2+1
            if (ns2.gt.n2) exit
         enddo
         ord(ns1+nc2:ns1+ns2-2)=ord2(nc2+1:ns2-1)
         do ipos=ns1+ns2-1,n1+n2
            if (ns1.gt.n1) then
               mv12=2
            elseif(ns2.gt.n2) then
               mv12=1
            elseif(ord1(ns1).lt.ord2(ns2)) then
               mv12=1
            else
               mv12=2
            endif
            singlet_mv(0)=singlet_mv(0)+1
            if (mv12.eq.1) then
               ord(ipos)=ord1(ns1)
               singlet_mv(singlet_mv(0))=ns1 - (singlet_mv(0)-1)
               ns1=ns1+1
            elseif (mv12.eq.2) then
               ord(ipos)=ord2(ns2)
               singlet_mv(singlet_mv(0))=n1+ns2 - (singlet_mv(0)-1)
               ns2=ns2+1
            endif
         enddo
      endif
      combine_currents%order(1:isize)=ord(1:isize)
      combine_currents%spin(1:isize)=combine_lists([spin1(1:n1),spin2(1:n2)],singlet_mv)
      combine_currents%ext_type(1:isize)=combine_lists([et1(1:n1),et2(1:n2)],singlet_mv)
    end function combine_currents
    
    subroutine add_all_currents(ctype)
      ! combine currents ic1 and ic2 and add them to the list of currents to
      ! compute. If use_symmetry=.true., need to consider all possible
      ! permutations allowed under the symmetry (at most 8 if both ic1 and ic2
      ! are made up of only gluons).
      implicit none
      logical,dimension(8) :: vertex_sign
      integer :: i,ctype,nperm
      integer,dimension(0:isize) :: singlet_mv
      type(current),dimension(8) :: new_currents
      if (.not.use_symmetry .or. this%imode.eq.1 .or. this%imode.eq.3) then
         new_currents(1)=combine_currents(ic1,ic2,ctype,singlet_mv,0)
         interaction_list_local(this%n_vert)%singlet_mv(0:isize)=singlet_mv(0:isize)
         call add_current(.false.,new_currents(1))
         return
      endif
      ! Need to consider all the possible permutations
      call check_all_permutations(ctype,nperm,new_currents,vertex_sign)
      do i=1,nperm
         call add_current(vertex_sign(i),new_currents(i))
      enddo
    end subroutine add_all_currents

    subroutine check_all_permutations(ctype,nperm,new_currents,vertex_sign)
      ! If a current only contains (external) gluons, we can use symmetry to
      ! relate them to eachother. This subroutine checks all permutations,
      ! and, if they give a valid current order, adds that current to the list
      ! that should be included.
      implicit none
      integer,intent(in) :: ctype
      integer,intent(out) :: nperm
      type(current),dimension(8),intent(out) :: new_currents
      logical,intent(out),dimension(8) :: vertex_sign
      logical :: ag1,ag2,iden
      integer,dimension(3) :: switch
      integer :: i,j,k,invert
      integer,dimension(0:isize,8) :: singlet_mv
      switch(1:3)=1
      ag1=all_gluon_current(current_list_local(ic1),n1)
      ag2=all_gluon_current(current_list_local(ic2),n2)
      if (n1.ge.2 .and. ag1) switch(1)=2
      if (n2.ge.2 .and. ag2) switch(2)=2
      if (ag1 .and. ag2) switch(3)=2
      nperm=0
      do i=1,switch(1)
         do j=1,switch(2)
            do k=1,switch(3)
               invert=0
               if ((i.eq.2 .and. k.eq.1) .or. (j.eq.2 .and. k.eq.2)) invert=ibset(invert,0)
               if ((i.eq.2 .and. k.eq.2) .or. (j.eq.2 .and. k.eq.1)) invert=ibset(invert,1)
               nperm=nperm+1
               if (k.eq.1) then
                  new_currents(nperm)=combine_currents(ic1,ic2,ctype,singlet_mv(0,nperm),invert)
               else
                  new_currents(nperm)=combine_currents(ic2,ic1,ctype,singlet_mv(0,nperm),invert)
               endif
               vertex_sign(nperm)=(k.eq.2 .xor. (j.eq.2 .and. mod(n2,2).eq.0) .xor. (i.eq.2 .and. mod(n1,2).eq.0))
               if (.not.valid_current_order_excl_symmetry(new_currents(nperm))) nperm=nperm-1
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
         interaction_list_local(this%n_vert)%singlet_mv(0:singlet_mv(0,1))=singlet_mv(0:singlet_mv(0,1),1)
      else
         write (*,*) 'Singlet move not identical for all permutations',nperm
         do i=1,nperm
            write (*,*) nperm,':',singlet_mv(0,i),':',singlet_mv(1:singlet_mv(0,i),i)
         enddo
         stop 1
      endif
    end subroutine check_all_permutations
    
    logical function valid_current_order_excl_symmetry(new_current)
      ! Checks that ip(1:isize) is an order for a current to be considered
      ! when use_symmetry=.true. --> the smallest number needs to come before
      ! the largest number in this list. 
      implicit none
      type(current),intent(in) :: new_current
      integer :: i,maxi,mini,min_loc,max_loc
      ! If it is not an all-gluon current, no symmetry can be used to reduce
      ! the number of currents to compute and therefore we have to include
      ! this current combination
      if (.not.all_gluon_current(new_current,isize)) then
         valid_current_order_excl_symmetry=.true.
         return
      endif
      ! This is an all-gluon current. Here we take only one single
      ! order. Define it such that smallest label comes before the
      ! biggest. This must be compatible with what orders are skipped in
      ! 'add_if_allowed_threevertex()'.
      maxi=0
      mini=100
      do i=1,isize
         if (new_current%order(i).gt.maxi) then
            maxi=new_current%order(i)
            max_loc=i
         endif
         if (new_current%order(i).lt.mini) then
            mini=new_current%order(i)
            min_loc=i
         endif
      enddo
      if (min_loc.gt.max_loc) then
         valid_current_order_excl_symmetry=.false.
         return
      endif
      valid_current_order_excl_symmetry=.true.
    end function valid_current_order_excl_symmetry


    subroutine add_current(vertex_sign,new_current)
      ! Adds the 'new_current' to the list of currents to consider. If this
      ! current has the same order (and type etc.) of an existing current, we
      ! can add this current to that existing one. If not, add it to the end
      ! of the list.
      implicit none
      type(current),intent(in) :: new_current
      logical,intent(in) :: vertex_sign
      integer :: ic,key
      integer(kind=8) :: val
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         ! Check if this interaction can be added to an existing current
         do ic=1,this%n_cur
            if (new_current%type.ne.current_list_local(ic)%type) cycle
            if (new_current%bin.ne.current_list_local(ic)%bin) cycle
            if (any(new_current%order(1:isize).ne.current_list_local(ic)%order(1:isize))) cycle
            if (any(new_current%spin(1:isize).ne.current_list_local(ic)%spin(1:isize))) cycle
            if (any(new_current%ext_type(1:isize).ne.current_list_local(ic)%ext_type(1:isize))) cycle
            current_list_local(ic)%n_vert=current_list_local(ic)%n_vert+1
            current_list_local(ic)%vertices(current_list_local(ic)%n_vert)=this%n_vert
            current_list_local(ic)%vertex_sign(current_list_local(ic)%n_vert)=vertex_sign
            return
         enddo
         ! Need a new current
         this%n_cur=this%n_cur+1
         if (this%n_cur.gt.max_cur) call increase_max_cur()
         current_list_local(this%n_cur)=new_current
         current_list_local(this%n_cur)%mass=pm%get_mass(new_current%type)
         current_list_local(this%n_cur)%width=pm%get_width(new_current%type)
         if (is_gluon(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         elseif (is_tensor(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(isize-1))
            allocate(current_list_local(this%n_cur)%vertex_sign(isize-1))
         elseif (is_massiveboson(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         else
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         endif
         current_list_local(this%n_cur)%vertices(1)=this%n_vert
         current_list_local(this%n_cur)%vertex_sign(1)=vertex_sign
         current_list_local(this%n_cur)%n_vert=1
      elseif (this%imode.eq.2) then
         call get_value(new_current%order,new_current%type,val)
         call solve_dict(val,key)
         ic=key_to_current(key,new_current%iproc)
         if (ic.eq.0) then
            ! initialise new current
            this%n_cur=this%n_cur+1
            if (this%n_cur.gt.max_cur) call increase_max_cur()
            key_to_current(key,new_current%iproc)=this%n_cur
            ic=this%n_cur
            current_list_local(ic)=new_current
            current_list_local(ic)%mass=pm%get_mass(new_current%type)
            current_list_local(ic)%width=pm%get_width(new_current%type)
            if (any(current_list_local(ic)%spin(1:isize).ne.-9)) then
               write (*,*) 'trying to combine currents with different spin: not possible',&
                    current_list_local(ic)%spin(1:isize)
               stop 1
            endif
            if (is_gluon(new_current%type)) then
               allocate(current_list_local(ic)%vertices(5*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(5*(isize-1)))
            elseif (is_tensor(new_current%type)) then
               allocate(current_list_local(ic)%vertices(isize-1))
               allocate(current_list_local(ic)%vertex_sign(isize-1))
            elseif (is_massiveboson(new_current%type)) then
               allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
               allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
            else
               allocate(current_list_local(ic)%vertices(5*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(5*(isize-1)))
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
      size=n
      max_key=n
      do isize=2,n-1
         size=size*(n-isize+1)
         do i=1,size
            max_key=max_key+14
         enddo
      enddo
      allocate(current_dict(max_key)) 
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
            do j=1,pm%npart
               key=key+1
               call get_value(ips_in,pm%particle_list(j)%type,val)
               if (val.le.previous_val) then
                  write (*,*) 'inconsistent current dictionary #1',val,previous_val
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
      integer :: j,itype,i
      integer(kind=8) :: val,offset
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
      ! Take the types into account (don't worry about particle
      ! vs. anti-particle, since there should be no confusion given
      ! the (external) particles that are part of the current).
      do j=1,pm%npart
         if (itype.eq.pm%particle_list(j)%type .or. itype.eq.pm%particle_list(j)%anti_type) then
            offset=int(j-1,kind=8)
         endif
      enddo
      val=val*int(pm%npart,kind=8) + offset
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

    integer function anti_current(ctype)
      implicit none
      integer :: ctype
      anti_current=pm%get_antipart(ctype)
    end function anti_current
    logical function all_gluon_current(curr,len)
      ! returns .true. only if all external particles are gluons
      implicit none
      type(current),intent(in) :: curr
      integer,intent(in) :: len
      if (any(curr%ext_type(1:len).ne.21)) then
         all_gluon_current=.false.
      else
         all_gluon_current=.true.
      endif
    end function all_gluon_current
    logical function all_singlet_current(curr,len)
      ! returns .true. only if all external particles are colour singlets
      implicit none
      type(current),intent(in) :: curr
      integer,intent(in) :: len
      integer :: i
      all_singlet_current=.true.
      do i=1,len
         if (.not.is_singlet(curr%ext_type(i))) then
            all_singlet_current=.false.
            return
         endif
      enddo
    end function all_singlet_current
    logical function is_quark_from_order(io,iproc)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io,iproc
      if ( (io.le.2  .and. is_antiquark(this%processes(io,iproc))) .or. &
           (io.gt.2  .and. is_quark(this%processes(io,iproc)))) then
         is_quark_from_order=.true.
      else
         is_quark_from_order=.false.
      endif
    end function is_quark_from_order
    logical function is_antiquark_from_order(io,iproc)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io,iproc
      if ( (io.le.2  .and. is_quark(this%processes(io,iproc))) .or. &
           (io.gt.2  .and. is_antiquark(this%processes(io,iproc)))) then
         is_antiquark_from_order=.true.
      else
         is_antiquark_from_order=.false.
      endif
    end function is_antiquark_from_order
    logical function quark_in_current(curr,len)
      implicit none
      type(current),intent(in) :: curr
      integer,intent(in) :: len
      integer :: i
      do i=1,len
         if (is_quark(curr%ext_type(i)).or.is_antiquark(curr%ext_type(i))) then
            quark_in_current=.true.
            return
         endif
      enddo
      quark_in_current=.false.
    end function quark_in_current

    subroutine deallocate_unneeded()
      ! Deallocate some of the current_list information that will not be used
      ! later. This saves a little memory.
      implicit none
      integer :: ic
      do isize=1,n
         do ic=this%n_cur_start(isize),this%n_cur_end(isize)
            if (allocated(this%current_list(ic)%ext_type)) deallocate(this%current_list(ic)%ext_type)
            if (isize.ne.1 .and. isize.ne.n) then
               if (allocated(this%current_list(ic)%spin)) deallocate(this%current_list(ic)%spin)
               if (allocated(this%current_list(ic)%order)) deallocate(this%current_list(ic)%order)
            endif
         enddo
      enddo
      deallocate(order)
      if (this%imode.eq.2) then
         deallocate(current_dict)
         deallocate(key_to_current)
      endif
    end subroutine deallocate_unneeded
    
  end subroutine init

  subroutine write_init_amps_to_file(this,n,iunit)
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,iunit,ic,iv,isize,iamp,iproc
    write (iunit) this%n_cur,this%n_vert,this%imode,this%nColOrd,this%max_pp,this%n_amps,this%nprocs
    write (iunit) this%n_cur_start(1:n)
    write (iunit) this%n_cur_end(1:n)
    write (iunit) this%n_vert_start(2:n-1)
    write (iunit) this%n_vert_end(2:n-1)
    ! current_list
    do isize=1,n-1
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          write (iunit) this%current_list(ic)%type,this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%iproc,this%current_list(ic)%mass,this%current_list(ic)%width
          write (iunit) this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)
          write (iunit) this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)
          if (isize.eq.1 .or. isize.eq.n)  write (iunit) this%current_list(ic)%order(1),this%current_list(ic)%spin(1)
       enddo
    enddo
    ! interaction_list
    do iv=1,this%n_vert
       write (iunit) this%interaction_list(iv)%type,this%interaction_list(iv)%currents(1:2),&
            this%interaction_list(iv)%coupl(1:2)
       if (allocated(this%interaction_list(iv)%singlet_mv)) then
          write (iunit) this%interaction_list(iv)%singlet_mv(0:this%interaction_list(iv)%singlet_mv(0))
       else
          write (iunit) 0
       endif
    enddo
    ! momenta array
    write (iunit) this%pp_bin_to_i(1:maskr(n))
    write (iunit) this%pp_i_to_bin(1:this%max_pp)
    ! process specific information
    do iproc=1,this%nprocs
       write (iunit) this%iproc_start(iproc),this%same_flav(iproc),this%same_flavour_proc_map(iproc,1:2),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       write (iunit) this%processes(1:n,iproc)
    enddo
    write(iunit) this%iproc_start(this%nprocs+1)
    ! amp specific information
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          write (iunit) this%include_amp(iamp),this%same_flavour_sum(iamp,1:2)
          write (iunit) this%spins(1:n,1,iamp)
          write (iunit) this%perm(1:n-this%n_sing(1),iamp)
          if (.not.this%same_flav(iproc)) write (iunit) this%curr2amp(1:2,iamp)
       enddo
    enddo
  end subroutine write_init_amps_to_file

  subroutine read_init_amps_from_file(this,n,iunit)
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,iunit,ic,iv,isize,iamp,iproc,itmp
    call deallocate_all()
    read (iunit) this%n_cur,this%n_vert,this%imode,this%nColOrd,this%max_pp,this%n_amps,this%nprocs
    allocate(this%n_cur_start(1:n))
    allocate(this%n_cur_end(1:n))
    read (iunit) this%n_cur_start(1:n)
    read (iunit) this%n_cur_end(1:n)
    allocate(this%n_vert_start(2:n-1))
    allocate(this%n_vert_end(2:n-1))
    read (iunit) this%n_vert_start(2:n-1)
    read (iunit) this%n_vert_end(2:n-1)
    ! current_list
    allocate(this%current_list(this%n_cur))
    do isize=1,n-1
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          read (iunit) this%current_list(ic)%type,this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%iproc,this%current_list(ic)%mass,this%current_list(ic)%width
          allocate(this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert))
          allocate(this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert))
          read (iunit) this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)
          read (iunit) this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)
          if (isize.eq.1 .or. isize.eq.n) then
             allocate(this%current_list(ic)%order(1))
             allocate(this%current_list(ic)%spin(1))
             read (iunit) this%current_list(ic)%order(1),this%current_list(ic)%spin(1)
          endif
       enddo
    enddo
    ! interaction_list
    allocate(this%interaction_list(1:this%n_vert))
    do iv=1,this%n_vert
       read (iunit) this%interaction_list(iv)%type,this%interaction_list(iv)%currents(1:2),this%interaction_list(iv)%coupl(1:2),itmp
       if (itmp.gt.0) then
          allocate(this%interaction_list(iv)%singlet_mv(0:itmp))
          this%interaction_list(iv)%singlet_mv(0)=itmp
          read (iunit) this%interaction_list(iv)%singlet_mv(1:itmp)
       endif
    enddo
    ! momenta array
    allocate(this%pp_bin_to_i(1:maskr(n)))
    allocate(this%pp_i_to_bin(1:this%max_pp))
    allocate(this%pp(0:3,1:this%max_pp))
    read (iunit) this%pp_bin_to_i(1:maskr(n))
    read (iunit) this%pp_i_to_bin(1:this%max_pp)
    ! process specific information
    allocate(this%iproc_start(1:this%nprocs+1))
    allocate(this%same_flav(1:this%nprocs))
    allocate(this%same_flavour_proc_map(1:this%nprocs,1:2))
    allocate(this%n_qqbar(1:this%nprocs))
    allocate(this%n_sing(1:this%nprocs))
    allocate(this%processes(1:n,1:this%nprocs))
    do iproc=1,this%nprocs
       read (iunit) this%iproc_start(iproc),this%same_flav(iproc),this%same_flavour_proc_map(iproc,1:2),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       read (iunit) this%processes(1:n,iproc)
    enddo
    read(iunit) this%iproc_start(this%nprocs+1)
    ! amp specific information
    allocate(this%include_amp(1:this%n_amps))
    allocate(this%same_flavour_sum(1:this%n_amps,1:2))
    allocate(this%spins(1:n,1,1:this%n_amps))
    allocate(this%perm(1:n-this%n_sing(1),1:this%n_amps))
    do iproc=1,this%nprocs
       if (this%same_flav(iproc)) exit
    enddo
    allocate(this%curr2amp(1:2,1:this%iproc_start(iproc)-1))
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          read (iunit) this%include_amp(iamp),this%same_flavour_sum(iamp,1:2)
          read (iunit) this%spins(1:n,1,iamp)
          read (iunit) this%perm(1:n-this%n_sing(1),iamp)
          if (.not.this%same_flav(iproc)) read (iunit) this%curr2amp(1:2,iamp)
       enddo
    enddo
  contains
    subroutine deallocate_all()
      implicit none
      integer :: i
      if (allocated(this%n_cur_start)) deallocate(this%n_cur_start)
      if (allocated(this%n_cur_end)) deallocate(this%n_cur_end)
      if (allocated(this%n_vert_start)) deallocate(this%n_vert_start)
      if (allocated(this%n_vert_end)) deallocate(this%n_vert_end)
      if (allocated(this%current_list)) then
         do i=1,size(this%current_list)
            call finalize_current(this%current_list(i))
         enddo
         deallocate(this%current_list)
      endif
      if (allocated(this%interaction_list)) then
         do i=1,size(this%interaction_list)
            call finalize_interaction(this%interaction_list(i))
         enddo
         deallocate(this%interaction_list)
      endif
      if (allocated(this%pp)) deallocate(this%pp)
      if (allocated(this%pp_bin_to_i)) deallocate(this%pp_bin_to_i)
      if (allocated(this%pp_i_to_bin)) deallocate(this%pp_i_to_bin)
      if (allocated(this%iproc_start)) deallocate(this%iproc_start)
      if (allocated(this%same_flav)) deallocate(this%same_flav)
      if (allocated(this%same_flavour_proc_map)) deallocate(this%same_flavour_proc_map)
      if (allocated(this%n_qqbar)) deallocate(this%n_qqbar)
      if (allocated(this%n_sing)) deallocate(this%n_sing)
      if (allocated(this%processes)) deallocate(this%processes)
      if (allocated(this%include_amp)) deallocate(this%include_amp)
      if (allocated(this%same_flavour_sum)) deallocate(this%same_flavour_sum)
      if (allocated(this%spins)) deallocate(this%spins)
      if (allocated(this%perm)) deallocate(this%perm)
      if (allocated(this%curr2amp)) deallocate(this%curr2amp)
    end subroutine deallocate_all
  end subroutine read_init_amps_from_file
  
  subroutine evaluate(this,n,p,hel,read_file,pm)
    use FeynmanRules
    use particles
    implicit none
    class(amplitude_QCD) :: this
    type(physics_model),intent(in) :: pm
    integer :: n
    integer,dimension(n)::hel
    real(kind=8),dimension(0:3,n) :: p
    integer :: ic,iv,isize,ih_in,ip,ifinal,dim
    logical :: read_file 
    if (.not. allocated(this%current_list(1)%val_c) .and. .not.allocated(this%current_list(1)%val_r)) then
       do ic=1,this%n_cur
          if (use_real_gluons .and. &
               (is_gluon(this%current_list(ic)%type) .or. is_tensor6(this%current_list(ic)%type))) then
             dim=pm%get_dim(this%current_list(ic)%type)
             allocate(this%current_list(ic)%val_r(1:dim))
          else
             dim=pm%get_dim(this%current_list(ic)%type)
             allocate(this%current_list(ic)%val_c(1:dim))
          endif
       enddo
       do iv=1,this%n_vert
          if (use_real_gluons .and. &
               (this%interaction_list(iv)%type.ge.0 .and. this%interaction_list(iv)%type.le.3)) then
             dim=pm%get_inter_dim(this%interaction_list(iv)%type)
             allocate(this%interaction_list(iv)%val_r(1:dim))
          else
             dim=pm%get_inter_dim(this%interaction_list(iv)%type)
             allocate(this%interaction_list(iv)%val_c(1:dim))
          endif
       enddo
       if (use_real_gluons .and. this%n_qqbar(1).eq.0) then
          allocate(this%amps_r(1:this%n_amps))
       else
          allocate(this%amps(1:this%n_amps))
       endif
    endif
    
    call fill_momentum_array()

    do isize=1,n-1
       if (isize.eq.1) then
          ! fill the external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)
             ifinal=1
             if (this%current_list(ic)%spin(1).eq.-9) then
                ih_in=0
             else
                ih_in=this%current_list(ic)%spin(1)
             endif
             if (is_gluon(this%current_list(ic)%type) .or. is_photon(this%current_list(ic)%type)) then
                if (use_real_gluons) then
                   call ext_gluon_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_r(1:4))
                else
                   call ext_gluon_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4))
                endif
             elseif (is_quark(this%current_list(ic)%type)) then
                call ext_quark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (is_antiquark(this%current_list(ic)%type)) then
                call ext_antiquark(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (is_massiveboson(this%current_list(ic)%type)) then
                call ext_gluon_mass(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (is_higgs(this%current_list(ic)%type)) then
                call ext_scalar(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ifinal,this%current_list(ic)%val_c(1))
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

          elseif(this%interaction_list(iv)%type.eq.10) then
             if (use_real_gluons) then
                call QuarkGluontoQuark_coupl_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                                  this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                                  this%interaction_list(iv)%val_c(1:4),&
                                                  this%interaction_list(iv)%coupl(1:2))
             else
                call QuarkGluontoQuark_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                             this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                             this%interaction_list(iv)%val_c(1:4),&
                                             this%interaction_list(iv)%coupl(1:2))
             endif
          elseif(this%interaction_list(iv)%type.eq.11) then
             if (use_real_gluons) then
                call AquarkGluontoAquark_coupl_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                                    this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                                    this%interaction_list(iv)%val_c(1:4),&
                                                    this%interaction_list(iv)%coupl(1:2))
             else
                call AquarkGluontoAquark_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                               this%interaction_list(iv)%val_c(1:4),&
                                               this%interaction_list(iv)%coupl(1:2))
             endif
          elseif (this%interaction_list(iv)%type.eq.12) then
             call threeGluon_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                       this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                       this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
                       this%interaction_list(iv)%val_c(1:4),&
                       this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.13) then
             call TwoGluonToTensor_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:6),&
                                         this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.14) then
             call TensorGluontoGluon_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:6),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                           this%interaction_list(iv)%val_c(1:4),&
                                           this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.15) then
             call GluonTensortoGluon_coupl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:6),&
                                           this%interaction_list(iv)%val_c(1:4),&
                                           this%interaction_list(iv)%coupl(1:2))

          else
             write (*,*) 'Unknown vertex type: not yet implemented',iv,this%interaction_list(iv)%type
             stop 1
          endif
       enddo

       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          if (is_gluon(this%current_list(ic)%type) .or. is_photon(this%current_list(ic)%type)) then
             call combine_interactions(4)
             ! a gluon current
             if (isize.ne.n-1)  then
                call include_gluon_propagator()
             endif
          elseif (is_quark(this%current_list(ic)%type)) then
             ! a quark current
             call combine_interactions(4)
             if (isize.ne.n-1)  then
                call include_quark_propagator()
             endif
          elseif (is_tensor6(this%current_list(ic)%type)) then
             ! the non-propagating tensor current
             call combine_interactions(6)
          elseif (is_antiquark(this%current_list(ic)%type)) then
             ! an anti-quark current
             call combine_interactions(4)
             if (isize.ne.n-1)  then
                call include_aquark_propagator()
             endif
          elseif (is_massiveboson(this%current_list(ic)%type)) then
             call combine_interactions(4)
             ! a massive vector boson current
             if (isize.ne.n-1)  then
                call include_gluon_propagator_mass()
             endif
          elseif (is_higgs(this%current_list(ic)%type)) then
             ! a scalar current
             call combine_interactions(1)
             if (isize.ne.n-1)  then
                call include_scalar_propagator()
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
      integer :: iamp,ih1,ih2,ih,ic,ihc,iproc
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         if (use_real_gluons .and. all(this%n_qqbar(1:this%nprocs).eq.0)) then
            do iamp=1,this%n_amps
               this%amps_r(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_r(1:4)* &
                                     this%current_list(this%curr2amp(2,iamp))%val_r(1:4))
            enddo
         else
            if (.not.read_file) then
               do iproc=1,this%nprocs
               do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (.not.this%same_flav(iproc)) then
                     this%amps(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_c(1:4)* &
                                         this%current_list(this%curr2amp(2,iamp))%val_c(1:4))
                  else
                     ! same-flavour amps are build from two different-flavour amps
                     if (this%same_flavour_sum(iamp,1).gt.0 .and. this%same_flavour_sum(iamp,2).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))+ &
                                        this%amps(this%same_flavour_sum(iamp,2))/3d0
                     elseif (this%same_flavour_sum(iamp,1).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))
                     elseif (this%same_flavour_sum(iamp,2).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,2))/3d0
                     else
                        write (*,*) 'At least one should contribute'
                        stop 1
                     endif
                  endif
               enddo
               enddo
            else   
             do iproc=1,this%nprocs
                do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (.not.this%same_flav(iproc)) then
                     this%amps(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_c(1:4)* &
                                         this%current_list(this%curr2amp(2,iamp))%val_c(1:4))
                  else
                     ! same-flavour amps are build from two different-flavour amps
                     if (this%same_flavour_sum(iamp,1).gt.0 .and. this%same_flavour_sum(iamp,2).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))+ &
                                        this%amps(this%same_flavour_sum(iamp,2))/3d0
                     elseif (this%same_flavour_sum(iamp,1).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))
                     elseif (this%same_flavour_sum(iamp,2).gt.0) then
                        this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,2))/3d0
                     else
                        write (*,*) 'At least one should contribute'
                        stop 1
                     endif
                  endif
               enddo
            enddo
            endif
         endif
      elseif(this%imode.eq.2) then
         if (use_real_gluons .and. this%n_qqbar(1).eq.0) then
            do iamp=1,this%n_amps
               if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. iamp.gt.this%n_amps/2 .and. mod(n,2).eq.1) then
                  this%amps_r(iamp)=-sum(this%current_list(this%curr2amp(1,iamp))%val_r(1:4)* &
                       this%current_list(this%curr2amp(2,iamp))%val_r(1:4))
               else
                  this%amps_r(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_r(1:4)* &
                       this%current_list(this%curr2amp(2,iamp))%val_r(1:4))
               endif
            enddo
         else
            do iproc=1,this%nprocs
               do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. iamp.gt.this%n_amps/2 .and. mod(n,2).eq.1) then
                     if (iproc.ne.1) then
                        write (*,*) 'this can only be one process at the time'
                        stop 1
                     endif
                     this%amps(iamp)=-sum(this%current_list(this%curr2amp(1,iamp))%val_c(1:4)* &
                                          this%current_list(this%curr2amp(2,iamp))%val_c(1:4))
                  else
                     if (.not.this%same_flav(iproc)) then
                        this%amps(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_c(1:4)* &
                                            this%current_list(this%curr2amp(2,iamp))%val_c(1:4))
                     else
                        ! same-flavour amps are build from two different-flavour amps
                        if (this%same_flavour_sum(iamp,1).gt.0 .and. this%same_flavour_sum(iamp,2).gt.0) then
                           this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))+ &
                                           this%amps(this%same_flavour_sum(iamp,2))/3d0
                        elseif (this%same_flavour_sum(iamp,1).gt.0) then
                           this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,1))
                        elseif (this%same_flavour_sum(iamp,2).gt.0) then
                           this%amps(iamp)=this%amps(this%same_flavour_sum(iamp,2))/3d0
                        else
                           write (*,*) 'At least one should contribute'
                           stop 1
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
      endif
    end subroutine compute_amps_from_currents

    subroutine combine_interactions(dim)
      implicit none
      integer :: dim,iv,i
      if (use_real_gluons .and. (is_gluon(this%current_list(ic)%type).or.is_tensor_g(this%current_list(ic)%type))) then
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

    subroutine include_gluon_propagator_mass()
      implicit none
      call GluonPropagator_mass(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%mass,this%current_list(ic)%width)
    end subroutine include_gluon_propagator_mass

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
    subroutine include_scalar_propagator()
      implicit none
      call ScalarPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_scalar_propagator
  end subroutine evaluate

  subroutine init_col(this,n,col_acc)
    use color_algebra
    use math_functions
    implicit none
    class(amplitude_qcd) :: this
    integer,parameter :: max_vals=10000
    integer :: col_acc,n,iperm,jperm,ival,iacc,isum,i,j,gi,gj,ui,uj,uj_upper,iperm_upper,&
         gi_iperm,key,max_keys,jperm_lower,ui_upper,n_unique_rows,irow,iunique,iproc,ioff,nOrd
    integer,dimension(n) :: iper,jper,part
    integer,dimension(:),allocatable :: n_vals
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(:,:),allocatable :: diff_vals
    real(kind=8),dimension(:,:,:),allocatable :: col_vals
    integer,dimension(:,:),allocatable :: ic,ir,n_colour_elements,unique_rows
    integer(kind=8),dimension(:),allocatable :: perm_dict

    write (*,*) 'Initialising colour matrix ...'
    if (this%nprocs.eq.1) then
       iproc=1
    elseif(this%nprocs.eq.3) then
       iproc=3
    else
       write (*,*) 'computation of color factor only for a single process at the time',this%nprocs
    endif
    part(1:n)=this%processes(1:n,1)
    ioff=this%iproc_start(iproc)-1
    nOrd=n-this%n_sing(iproc)

    allocate(n_vals(1:3))
    allocate(diff_vals(max_vals,1:3))
    allocate(this%i_col_i(max_vals,1:3))
    allocate(n_colour_elements(max_vals,1:3))
    
    if (this%n_qqbar(iproc).eq.0 .or. this%n_qqbar(iproc).eq.1) then
       n_unique_rows=1 ! all rows are similar
    elseif (this%n_qqbar(iproc).eq.2) then
       n_unique_rows=(nOrd-4)+1 ! number of gluon separations among the two quark lines
       n_unique_rows=n_unique_rows*2 ! two ways of combining quarks with anti-quarks
    else
       write (*,*) 'Inconsistent number of quark pairs',this%n_qqbar(iproc)
       stop 1
    endif
    if (use_cm_dict) then
       call create_perm_dict()
       allocate(col_vals(1:3,max_keys,n_unique_rows))
       allocate(unique_rows(1:nOrd,n_unique_rows))
    endif

! first check the unique rows in the colour matrix to determine how many
! different colour factors there are. This also sets up the library for the
! colour factors if use_cm_dict=.true.
    n_vals(1:3)=0
    if (use_cm_dict) col_vals(1:3,1:max_keys,:)=0d0
    do iunique=1,n_unique_rows
       if (this%n_qqbar(iproc).eq.0 .or. this%n_qqbar(iproc).eq.1) then
          ! just take the first row: all rows are similar
          irow=1
       elseif (this%n_qqbar(iproc).eq.2) then
          ! determine how the quarks are connected and how many gluons are on
          ! each quark line and make sure that it is not similar to one
          ! already considered
          call get_unique_row(iunique,irow,gi,ui)
       endif
       iper(1:nOrd)=this%perm(1:nOrd,ioff+irow)
       if (use_cm_dict) unique_rows(1:nOrd,iunique) = iper(1:nOrd)
       ! loop over the columns
       do jperm=1,this%nColOrd
          jper(1:nOrd)=this%perm(1:nOrd,ioff+jperm)
          if (this%n_qqbar(iproc).eq.2) then
             ! determine how the quarks are connected
             call determine_gi_ui(jper,gj,uj)
          endif
          call compute_color_factor(col_acc,nOrd,iper,jper,ui,uj,gi,gj,col_fac)
          if (use_symm_cm) then
             ! include a factor 2 for the off-diagonal terms
             if (irow.ne.jperm) col_fac(1:3)=col_fac(1:3)*2d0 
          endif
          do iacc=1,3
             if (col_fac(iacc).eq.0d0) cycle
             if (use_cm_dict) then
                key=solve_dict(get_value(nOrd,jper(1:nOrd)))
                col_vals(iacc,key,iunique)=col_fac(iacc)
             endif
             do ival=1,n_vals(iacc)
                if (col_fac(iacc).eq.diff_vals(ival,iacc)) then
                   n_colour_elements(ival,iacc)=n_colour_elements(ival,iacc)+1
                   exit
                endif
             enddo
             if (ival.ge.max_vals) then
                write (*,*) 'Too many different colour factors. Increase max_vals',&
                     ival,n_vals(1:3),max_vals
                stop 1
             elseif (ival.eq.n_vals(iacc)+1) then
                ! new colour factor
                n_vals(iacc)=ival
                diff_vals(ival,iacc)=col_fac(iacc)
                n_colour_elements(ival,iacc)=1
             endif
          enddo
       enddo
    enddo
    write (*,*) 'A single row in the colour matrix has',n_vals(1:3),&
         ' different colour factors at LC, NLC and full colour, respectively'
    
    ! determine i_col_i:
    isum=1
    do iacc=1,3
       do ival=1,n_vals(iacc)
          this%i_col_i(ival,iacc)=isum
          isum=isum+n_colour_elements(ival,iacc)*this%nColOrd
       enddo
    enddo

 ! Allocate the arrays now that we know their sizes
    allocate(ic(1:maxval(n_vals(1:3)),1:3))
    allocate(ir(1:maxval(n_vals(1:3)),1:3))
    allocate(this%col_index(1:isum))
    allocate(this%row_index(0:this%nColOrd,1:maxval(n_vals(1:3)),1:3)) 
    this%row_index(0,1:maxval(n_vals(1:3)),1:3)=0
    this%col_index(1)=0
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
       iper(1:nOrd)=this%perm(1:nOrd,ioff+iperm)
       if (this%n_qqbar(iproc).eq.2)  call determine_gi_ui(iper,gi,ui)
       if (use_symm_cm) then
          jperm_lower=iperm
       else
          jperm_lower=1
       endif
       do jperm=jperm_lower,this%nColOrd
          jper(1:nOrd)=this%perm(1:nOrd,ioff+jperm)
          if (this%n_qqbar(iproc).eq.2)  call determine_gi_ui(jper,gj,uj)
          if (use_cm_dict) then
             ! GET color factors from permuting first row
             call get_col_fac(iper,jper,ui,uj,gi,gj,col_fac)
          else
             ! COMPUTE color factors again
             call compute_color_factor(col_acc,nOrd,iper,jper,ui,uj,gi,gj,col_fac)
             if (use_symm_cm) then
                ! include a factor 2 for the off-diagonal terms
                if (iperm.ne.jperm) col_fac(1:3)=col_fac(1:3)*2d0
             endif
          endif
          do iacc=1,3
             if (col_fac(iacc).eq.0d0) cycle
             do ival=1,n_vals(iacc)
                if (col_fac(iacc).eq.diff_vals(ival,iacc)) exit
             enddo
             if (ival.gt.n_vals(iacc)) then
                write (*,*) 'colour factor not encountered during setup',ival,n_vals(iacc)
                stop 1
             endif
             ic(ival,iacc)=ic(ival,iacc)+1
             ir(ival,iacc)=ir(ival,iacc)+1
             this%col_index(this%i_col_i(ival,iacc)+ic(ival,iacc))=jperm
          enddo
       enddo
       do iacc=1,3
          this%row_index(iperm,1:n_vals(iacc),iacc)=ir(1:n_vals(iacc),iacc)
       enddo
    enddo

    write (*,*) '... colour matrix initialised'
  contains
    subroutine get_unique_row(iunique,irow,gi,ui)
      ! get a new row in the colour matrix corresponding to 'iunique'
      implicit none
      integer,intent(in) :: iunique
      integer,intent(out) :: irow,gi,ui
      do irow=1,this%nColOrd
         call determine_gi_ui(this%perm(1,ioff+irow),gi,ui)
         if (ui.eq.1 .and. gi.eq.iunique-1) return
         if (ui.eq.2 .and. gi.eq.iunique-1-((nOrd-4)+1)) return
      enddo
    end subroutine get_unique_row

    subroutine determine_gi(iper,gi)
      ! determine how many gluons are on the first quark line
      implicit none
      integer :: gi,i
      integer,dimension(n) :: iper
      do i=2,n-2
         if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
            gi=i-2
            return
         endif
      enddo
    end subroutine determine_gi

    subroutine determine_ui(iper,ui)
      implicit none
      integer :: ui
      integer,dimension(n) :: iper
      if (abs(part(iper(1))).eq.abs(part(iper(n-this%n_sing(1))))) then
         ui=1
      else
         ui=2
      endif
    end subroutine determine_ui
    
    subroutine determine_gi_ui(iper,gi,ui)
      implicit none
      integer :: ui,gi
      integer,dimension(n) :: iper
      call determine_gi(iper,gi)
      call determine_ui(iper,ui)
    end subroutine determine_gi_ui
    
   subroutine get_col_fac(iper,jper,ui,uj,gi,gj,col_fac)
     implicit none
     integer,intent(in) :: gi,gj,ui,uj
     integer,dimension(n),intent(in) :: iper,jper
     integer,dimension(n) :: col_new,row_first,row_per,col_per
     integer :: i,j,key,iunique
     real(kind=8),dimension(1:3),intent(out) :: col_fac
     if (this%n_qqbar(iproc).eq.2) then
        iunique=(gi+1)+(ui-1)*((nOrd-4)+1)
     else
        iunique=1
     endif
     ! First row
     row_first(1:nOrd)=unique_rows(1:nOrd,iunique)
     ! Row in consideration
     row_per(1:nOrd)=iper(1:nOrd)
     ! Column in consideration
     col_per(1:nOrd)=jper(1:nOrd)
     ! Map one into the other
     do i=1,nOrd
        do j=1,nOrd
           if (col_per(i) .eq. row_per(j)) exit
        enddo
        if (.not.(abs(part(col_per(i))).le.6.and.abs(part(col_per(i))).ge.1)) then
          col_new(i) = row_first(j)
        else 
          col_new(i) = col_per(i)
        endif
     enddo
     key=solve_dict(get_value(nOrd,col_new(1:nOrd)))
     col_fac(1:3)=col_vals(1:3,key,iunique)
   end subroutine get_col_fac

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
      allocate(iper(1:nOrd))
      allocate(iper_in(1:nOrd))
      max_keys=factorial(n)/factorial(n-nOrd)
      allocate(perm_dict(1:max_keys))
      do i=1,nOrd
         iper(i)=i
      enddo
      previous_val=0
      do iperm=1,max_keys
         val=get_value(nOrd,iper(1:nOrd))
         if (val.le.previous_val) then
            write (*,*) 'In create_perm_dict need to get values in ascending order',val,previous_val
            write (*,*) iper(1:nOrd)
            stop 1
         else
            previous_val=val
         endif
         perm_dict(iperm)=val
         iper_in=iper
         call get_next_iperm(nOrd,iper_in,iper,n)
      enddo
      deallocate(iper)
    end subroutine create_perm_dict

    integer(kind=8) function get_value(nOrd,iper)
      ! Give a unique identifier based on the colour order. Simply convert the
      ! list to an integer with base equal to the number of elements in the
      ! order.
      implicit none
      integer :: j,nOrd
      integer,dimension(1:nOrd) :: iper
      get_value=0
      do j=1,nOrd
         get_value=get_value+int(iper(nOrd+1-j),kind=8)*int(n+1,kind=8)**int(j-1,kind=8)
      enddo
    end function get_value
    
    subroutine compute_color_factor(col_acc,n,iper,jper,ui,uj,gi,gj,col_fac)
      use color_algebra
      implicit none
      integer :: i,n,acc,col_acc,k,ui,uj,gi,gj
      real(kind=8),dimension(1:3) :: col_fac
      integer,dimension(n) :: iper,jper
      integer,dimension(n-4) :: iper_glu,jper_glu,iper_ord,jper_ord
      real(kind=16) :: col_factor
      col_fac(1:3)=0d0
      if (col_acc.ge.0) then ! LC
         if (this%n_qqbar(iproc).eq.0) then
            if (all(iper.eq.jper)) then
               col_fac(1)=dble(3**n)
            endif
         elseif (this%n_qqbar(iproc).eq.1) then
            if (all(iper.eq.jper)) then
               col_fac(1)=dble(3**(n-1))
            endif
         elseif (this%n_qqbar(iproc).eq.2) then
            if (all(iper.eq.jper)) then
               if (ui.eq.1 .and. uj.eq.1) then
                  col_fac(1)=dble(3**(n-2))
               elseif (ui.eq.2 .and. uj.eq.2 .and. .not.this%same_flav(iproc)) then
                  col_fac(1)=dble(3**(n-4))
               elseif (ui.eq.2 .and. uj.eq.2 .and. this%same_flav(iproc)) then
                  col_fac(1)=dble(3**(n-2))
               endif
            endif
         endif
      endif
      if (col_acc.ge.1) then ! NLC
         if (this%n_qqbar(iproc).eq.0) then
            if (all(iper.eq.jper)) then
               col_fac(2) = dble(3**n - n * 3**(n-2))
            else
               call check_NLC(n,jper,iper,acc)
               col_fac(2)=dble(acc*3**(n-2))
            endif
         elseif (this%n_qqbar(iproc).eq.1) then
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
         elseif (this%n_qqbar(iproc).eq.2) then
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
            call convert_gluon_string(n,iper_glu,jper_glu,iper_ord,jper_ord)
            if (.not.this%same_flav(iproc)) then
               call check_NLC_2qqbar(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
               if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
               if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-4)) ! NLC parts
            else
               call check_NLC_2qqbar_SF(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
               if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
               if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-3)) ! NLC parts
            endif
            ! include the full expansion
            if (acc.ne.0) then
               call Tr_allocate(n)
               if (ui.eq.uj) then
                  Tr(0,0,0)=1 ! one term
                  Tr(0,0,1)=2
                  Tr(0,1,1) = gi+gj  ! number of generators in first trace
                  Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
                  Tr(1:gi,1,1) = iper(2:2+gi-1)
                  Tr(gi+1:gi+gj,1,1) = jper(2+gj-1:2:-1)
                  Tr(1:n-4-gi,2,1) = iper(2+gi+2:n-1)
                  Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-1:2+gj+2:-1)
                  if (ui.eq.2 .and. uj.eq.2 .and. .not.this%same_flav(iproc)) then
                     coef(1)=1d0/9d0
                  else
                     coef(1)=1d0
                  endif
                  call Tr_full_simplify(col_factor)
               elseif ((ui.eq.2 .and. uj.eq.1) .or. (ui.eq.1 .and. uj.eq.2)) then
                  Tr(0,0,0)=1 ! one term
                  Tr(0,0,1)=1 ! a single trace
                  Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
                  Tr(1:gi,1,1) = iper(2:2+gi-1)
                  Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-1:2+gj+2:-1)
                  Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(2+gi+2:n-1)
                  Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(2+gj-1:2:-1)
                  if (.not.this%same_flav(iproc)) then
                     coef(1)=-1d0/3d0
                  else
                     coef(1)=-1d0
                  endif
                  call Tr_full_simplify(col_factor)
               else
                  write (*,*) 'Inconsistent ui and uj',ui,uj
                  stop 1
               endif
               col_fac(2)=dble(col_factor)
               call Tr_deallocate
            endif
         endif
      endif
      if (col_acc.ge.2) then
         call Tr_allocate(n)
         if (this%n_qqbar(iproc).eq.0) then
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
         elseif (this%n_qqbar(iproc).eq.1) then
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
         elseif (this%n_qqbar(iproc).eq.2) then
            if (ui.eq.uj) then
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=2
               Tr(0,1,1) = gi+gj  ! number of generators in first trace
               Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
               Tr(1:gi,1,1) = iper(2:2+gi-1)
               Tr(gi+1:gi+gj,1,1) = jper(2+gj-1:2:-1)
               Tr(1:n-4-gi,2,1) = iper(2+gi+2:n-1)
               Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-1:2+gj+2:-1)
               if (ui.eq.2 .and. uj.eq.2 .and. .not.this%same_flav(iproc)) then
                  coef(1)=1d0/9d0
               else
                  coef(1)=1d0
               endif
               call Tr_full_simplify(col_factor)
            elseif ((ui.eq.2 .and. uj.eq.1) .or. (ui.eq.1 .and. uj.eq.2)) then
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=1 ! a single trace
               Tr(0,1,1) = 2*(n-4) ! all gluon generators appear in the single trace
               Tr(1:gi,1,1) = iper(2:2+gi-1)
               Tr(gi+1:gi+(n-4-gj),1,1) = jper(n-1:2+gj+2:-1)
               Tr(gi+(n-4-gj)+1:2*(n-4)-gj,1,1) = iper(2+gi+2:n-1)
               Tr(2*(n-4)-gj+1:2*(n-4),1,1) = jper(2+gj-1:2:-1)
               if (.not.this%same_flav(iproc)) then
                  coef(1)=-1d0/3d0
               else
                  coef(1)=-1d0
               endif
               call Tr_full_simplify(col_factor)
            else
               write (*,*) 'Inconsistent ui and uj',ui,uj
               stop 1
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
    integer,intent(in) :: n
    integer,intent(inout) :: nhel
    integer,intent(inout),dimension(nhel) :: include_hel
    integer :: nspin,ispin,ic,iv,iamp
    logical,dimension(:),allocatable :: include_current
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

    this%include_amp(1:this%n_amps)=.false.
    
    allocate(tmp_spin(1:n,1:maxval(include_hel),nhel))

    nspin=0
    do iamp=1,nhel
       if (include_hel(iamp).ge.1) then
          this%include_amp(iamp)=.true.
          if (this%same_flavour_sum(iamp,1).le.0) then
             include_current(this%curr2amp(1,iamp))=.true.
             include_current(this%curr2amp(2,iamp))=.true.
          else
             ! same-flavour amplitude.
             if (all(include_hel(this%same_flavour_sum(iamp,1:2)).eq.0)) then
                write (*,*) 'inconsistency in helicity filter for same-flavour process'
                stop 1
             endif
             if (include_hel(this%same_flavour_sum(iamp,1)).lt.0) &
                  this%same_flavour_sum(iamp,1)=-include_hel(this%same_flavour_sum(iamp,1))
             if (include_hel(this%same_flavour_sum(iamp,2)).lt.0) &
                  this%same_flavour_sum(iamp,2)=-include_hel(this%same_flavour_sum(iamp,2))
          endif
          nspin=nspin+1
          tmp_spin(1:n,1,nspin)=this%spins(1:n,1,iamp)
          ic=1
          do ispin=iamp+1,nhel
             if (-include_hel(ispin).eq.iamp) then
                ic=ic+1
                tmp_spin(1:n,ic,nspin)=this%spins(1:n,1,ispin)
             endif
          enddo
          include_hel(nspin)=include_hel(iamp)
       endif
    enddo

    deallocate(this%spins)
    allocate(this%spins(1:n,1:maxval(include_hel),nspin))
    this%spins(1:n,1:maxval(include_hel),1:nspin)=tmp_spin(1:n,1:maxval(include_hel),1:nspin)

    call this%filter_dead_trees(n,include_current)
    
    nhel=this%n_amps
    write (*,*) 'Total number of currents and vertices after filtering helicities',this%n_cur,this%n_vert,this%n_amps
    deallocate(this%include_amp)

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
    integer,dimension(:),allocatable :: where_to_cur,where_to_ver,where_to_amp
    logical,dimension(*),optional :: include_current
    integer :: to_skip,isize,nc,iv,n,iamp,i,iproc
    allocate(is_needed_cur(this%n_cur))
    allocate(is_needed_ver(this%n_vert))
    allocate(where_to_cur(this%n_cur))
    allocate(where_to_ver(this%n_vert))
    allocate(where_to_amp(0:this%n_amps))
    ! assume nothing is needed
    is_needed_cur(:)=.false.
    is_needed_ver(:)=.false.
    where_to_cur=0
    where_to_ver=0
    where_to_amp=0
    if (.not.present(include_current)) then
       is_needed_cur(this%n_cur_start(n-1):this%n_cur_end(n-1))=.true.
       is_needed_cur(this%n_cur_start(n  ):this%n_cur_end(n  ))=.true.
    else
       is_needed_cur(this%n_cur_start(n-1):this%n_cur_end(n-1))=include_current(this%n_cur_start(n-1):this%n_cur_end(n-1))
       is_needed_cur(this%n_cur_start(n  ):this%n_cur_end(n  ))=include_current(this%n_cur_start(n  ):this%n_cur_end(n  ))
    endif
    ! Since currents are created from previous ones, we should go backwards
    ! throught the list. If we encounter a 'is_needed_cur=.true.', it means
    ! that the (two, or more) currents that were combined to created that
    ! current, are also 'needed'. This determines all the currents (and
    ! vertices) that need to be kept.
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
    to_skip=0
    do iamp=1,this%n_amps
       if (.not. this%include_amp(iamp)) then
          to_skip=to_skip+1
          cycle
       endif
       where_to_amp(iamp)=iamp-to_skip
    enddo
    ! do the actual shifting of the currents in the list
    do nc=1,this%n_cur
       if (.not.is_needed_cur(nc)) cycle
       if (where_to_cur(nc).ne.nc) then
          if (allocated(this%current_list(where_to_cur(nc))%vertices)) &
               deallocate(this%current_list(where_to_cur(nc))%vertices)
          if (allocated(this%current_list(where_to_cur(nc))%vertex_sign)) &
               deallocate(this%current_list(where_to_cur(nc))%vertex_sign)
          this%current_list(where_to_cur(nc))=this%current_list(nc)
       endif
       do iv=1,this%current_list(where_to_cur(nc))%n_vert
          this%current_list(where_to_cur(nc))%vertices(iv)= &
               where_to_ver(this%current_list(where_to_cur(nc))%vertices(iv))
       enddo
    enddo
    ! do the actual shifting of the interactions in the list
    do iv=1,this%n_vert
       if (.not.is_needed_ver(iv)) cycle
       if (where_to_ver(iv).ne.iv) then
          if (allocated(this%interaction_list(where_to_ver(iv))%singlet_mv)) &
               deallocate(this%interaction_list(where_to_ver(iv))%singlet_mv)
          this%interaction_list(where_to_ver(iv))=this%interaction_list(iv)
       endif
       this%interaction_list(where_to_ver(iv))%currents(1:2)= &
            where_to_cur(this%interaction_list(where_to_ver(iv))%currents(1:2))
    enddo
    ! do the actual shifting of the amplitudes in the list
    do iamp=1,this%n_amps
       if (.not.this%include_amp(iamp)) cycle
       if (this%same_flavour_sum(iamp,1).le.0) then
          this%curr2amp(1,where_to_amp(iamp))=where_to_cur(this%curr2amp(1,iamp))
          this%curr2amp(2,where_to_amp(iamp))=where_to_cur(this%curr2amp(2,iamp))
       else
          this%same_flavour_sum(where_to_amp(iamp),1)=where_to_amp(this%same_flavour_sum(iamp,1))
          this%same_flavour_sum(where_to_amp(iamp),2)=where_to_amp(this%same_flavour_sum(iamp,2))
       endif
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
    do iproc=2,this%nprocs
       do iamp=this%iproc_start(iproc),this%n_amps
          if (where_to_amp(iamp).ne.0) then
             this%iproc_start(iproc)=where_to_amp(iamp)
             exit
          endif
       enddo
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
    do iamp=this%n_amps,1,-1
       if (where_to_amp(iamp).ne.0) then
          this%n_amps=where_to_amp(iamp)
          exit
       endif
    enddo
    this%iproc_start(this%nprocs+1)=this%n_amps+1
    deallocate(is_needed_ver)
    deallocate(is_needed_cur)
    deallocate(where_to_ver)
    deallocate(where_to_cur)
  end subroutine filter_dead_trees

  subroutine assign_interaction(lhs,rhs)
    ! sets non-custom 'lhs' = 'rhs' for interactions
    use particles
    implicit none
    type(interaction),intent(inout) :: lhs
    type(interaction),intent(in) :: rhs
    integer :: val_size
    lhs%type=rhs%type
    lhs%currents(1:2)=rhs%currents(1:2)
    lhs%coupl(1:2)=rhs%coupl(1:2)
    if (allocated(rhs%singlet_mv)) then
       if (rhs%singlet_mv(0).gt.0) then
          if (.not.allocated(lhs%singlet_mv)) allocate(lhs%singlet_mv(0:rhs%singlet_mv(0)))
          lhs%singlet_mv(0:rhs%singlet_mv(0))=rhs%singlet_mv(0:rhs%singlet_mv(0))
       elseif (allocated(lhs%singlet_mv)) then
          lhs%singlet_mv(0)=0
       endif
    endif
    if (allocated(lhs%val_c)) deallocate(lhs%val_c)
    if (allocated(rhs%val_c)) then
       if(is_tensor6(lhs%type)) then
          val_size=6
       else
          val_size=4
       endif
       allocate(lhs%val_c(1:val_size))
       lhs%val_c(1:val_size)=rhs%val_c(1:val_size)
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       if(is_tensor6(lhs%type)) then
          val_size=6
       else
          val_size=4
       endif
       allocate(lhs%val_r(1:val_size))
       lhs%val_r(1:val_size)=rhs%val_r(1:val_size)
    endif
  end subroutine assign_interaction
  
  subroutine assign_current(lhs,rhs)
    ! sets non-custom 'lhs' = 'rhs' for currents
    use particles
    implicit none
    type(current),intent(inout) :: lhs
    type(current),intent(in) :: rhs
    integer :: isize,val_size
    lhs%type=rhs%type
    lhs%bin=rhs%bin
    isize=popcnt(lhs%bin)
    lhs%n_vert=rhs%n_vert
    lhs%iproc=rhs%iproc
    lhs%mass=rhs%mass
    lhs%width=rhs%width
    if (allocated(rhs%vertices) .and. rhs%n_vert.gt.0) then
       if (.not.allocated(lhs%vertices)) allocate(lhs%vertices(1:lhs%n_vert))
       lhs%vertices(1:lhs%n_vert)=rhs%vertices(1:lhs%n_vert)
    endif
    if (allocated(lhs%order)) deallocate(lhs%order)
    if (allocated(rhs%order) .and. isize.gt.0) then
       allocate(lhs%order(1:isize))
       lhs%order(1:isize)=rhs%order(1:isize)
    endif
    if (allocated(lhs%spin)) deallocate(lhs%spin)
    if (allocated(rhs%spin) .and. isize.gt.0) then
       allocate(lhs%spin(1:isize))
       lhs%spin(1:isize)=rhs%spin(1:isize)
    endif
    if (allocated(lhs%ext_type)) deallocate(lhs%ext_type)
    if (allocated(rhs%ext_type) .and. isize.gt.0) then
       allocate(lhs%ext_type(1:isize))
       lhs%ext_type(1:isize)=rhs%ext_type(1:isize)
    endif
    if (allocated(rhs%vertex_sign) .and. rhs%n_vert.gt.0) then
       if (.not.allocated(lhs%vertex_sign)) allocate(lhs%vertex_sign(1:lhs%n_vert))
       lhs%vertex_sign(1:lhs%n_vert)=rhs%vertex_sign(1:lhs%n_vert)
    endif
    if (allocated(lhs%val_c)) deallocate(lhs%val_c)
    if (allocated(rhs%val_c)) then
       if(is_tensor6(lhs%type)) then
          val_size=6
       else
          val_size=4
       endif
       allocate(lhs%val_c(1:val_size))
       lhs%val_c(1:val_size)=rhs%val_c(1:val_size)
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       if(is_tensor6(lhs%type)) then
          val_size=6
       else
          val_size=4
       endif
       allocate(lhs%val_r(1:val_size))
       lhs%val_r(1:val_size)=rhs%val_r(1:val_size)
    endif
  end subroutine assign_current
  subroutine finalize_amplitude_QCD(amp)
    type(amplitude_QCD),intent(inout) :: amp
    integer :: i
    if (allocated(amp%current_list)) then
       do i=1,size(amp%current_list)
          call finalize_current(amp%current_list(i))
       enddo
       deallocate(amp%current_list)
    endif
    if (allocated(amp%interaction_list)) then
       do i=1,size(amp%interaction_list)
          call finalize_interaction(amp%interaction_list(i))
       enddo
       deallocate(amp%interaction_list)
    endif
    if (allocated(amp%amps)) deallocate(amp%amps)
    if (allocated(amp%amps_r)) deallocate(amp%amps_r)
    if (allocated(amp%pp)) deallocate(amp%pp)
    if (allocated(amp%diff_col_vals)) deallocate(amp%diff_col_vals)
    if (allocated(amp%n_cur_start)) deallocate(amp%n_cur_start)
    if (allocated(amp%n_cur_end)) deallocate(amp%n_cur_end)
    if (allocated(amp%n_vert_start)) deallocate(amp%n_vert_start)
    if (allocated(amp%n_vert_end)) deallocate(amp%n_vert_end)
    if (allocated(amp%pp_bin_to_i)) deallocate(amp%pp_bin_to_i)
    if (allocated(amp%pp_i_to_bin)) deallocate(amp%pp_i_to_bin)
    if (allocated(amp%col_index)) deallocate(amp%col_index)
    if (allocated(amp%n_col_vals)) deallocate(amp%n_col_vals)
    if (allocated(amp%iproc_start)) deallocate(amp%iproc_start)
    if (allocated(amp%n_sing)) deallocate(amp%n_sing)
    if (allocated(amp%n_qqbar)) deallocate(amp%n_qqbar)
    if (allocated(amp%perm)) deallocate(amp%perm)
    if (allocated(amp%curr2amp)) deallocate(amp%curr2amp)
    if (allocated(amp%i_col_i)) deallocate(amp%i_col_i)
    if (allocated(amp%processes)) deallocate(amp%processes)
    if (allocated(amp%same_flavour_proc_map)) deallocate(amp%same_flavour_proc_map)
    if (allocated(amp%same_flavour_sum)) deallocate(amp%same_flavour_sum)
    if (allocated(amp%spins)) deallocate(amp%spins)
    if (allocated(amp%row_index)) deallocate(amp%row_index)
    if (allocated(amp%include_amp)) deallocate(amp%include_amp)
    if (allocated(amp%same_flav)) deallocate(amp%same_flav)
  end subroutine finalize_amplitude_QCD
  subroutine finalize_interaction(vert)
    type(interaction),intent(inout) :: vert
    if (allocated(vert%singlet_mv)) deallocate(vert%singlet_mv)
    if (allocated(vert%val_c)) deallocate(vert%val_c)
    if (allocated(vert%val_r)) deallocate(vert%val_r)
  end subroutine finalize_interaction
  subroutine finalize_current(cur)
    type(current),intent(inout) :: cur
    if (allocated(cur%vertices)) deallocate(cur%vertices)
    if (allocated(cur%order)) deallocate(cur%order)
    if (allocated(cur%spin)) deallocate(cur%spin)
    if (allocated(cur%ext_type)) deallocate(cur%ext_type)
    if (allocated(cur%vertex_sign)) deallocate(cur%vertex_sign)
    if (allocated(cur%val_c)) deallocate(cur%val_c)
    if (allocated(cur%val_r)) deallocate(cur%val_r)
  end subroutine finalize_current
end module amplitude_QCD_mod
