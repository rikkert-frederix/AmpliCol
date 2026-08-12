module amplitude_QCD_mod
  use bitset_mod
  implicit none
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: use_symm_cm=.true.
  logical,parameter :: use_cm_dict=.true.
  integer(kind=8),parameter :: max_three_line_color_orders=5000_8
  type :: current
     ! if adding variables here, also update the finalize_current and assign_current subroutines
     integer :: type,bin,n_vert,chirality,n_gs=0,n_ew=0,open_quark_leg=0
     integer,dimension(3) :: ew_pairs=0,u1_pairs=0,gluon_pairs=0
     integer,dimension(3) :: u1_links=0,gluon_links=0
     type(bitset) :: iproc
     integer(kind=16) :: ext_cur
     integer,dimension(:),allocatable :: vertices,order,spin,ext_type
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8) :: mass,width
     integer,dimension(:),allocatable :: fermi_list
   contains
     final :: finalize_current ! custom deallocation of current
  end type current
  type :: interaction
     ! if adding variables here, also update the finalize_interaction and assign_interaction subroutines
     integer :: type,chirality,n_gs=0,n_ew=0
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
     integer :: n_cur,n_vert,imode,nColOrd,max_pp,n_amps,nprocs,n_sectors=0
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     complex(kind=8),dimension(:,:),allocatable :: amps_by_order
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: pp,diff_col_vals
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end, &
          pp_bin_to_i,pp_i_to_bin,col_index,n_col_vals,iproc_start,n_sing,n_qqbar
     integer,dimension(:,:),allocatable :: perm,curr2amp,i_col_i,processes,&
          same_flavour_sum,same_flavour_sum_operation,three_line_partner_curr2amp
     integer,dimension(:,:),allocatable :: sector_powers
     integer,dimension(:,:),allocatable :: sector_sign
     integer,dimension(:,:),allocatable :: sector_term_start,sector_term_curr2amp
     integer,dimension(:),allocatable :: sector_term_sign
     integer,dimension(:,:,:),allocatable :: sector_curr2amp,sector_three_line_partner_curr2amp
     integer,dimension(:,:,:),allocatable :: spins,row_index
     logical,dimension(:),allocatable :: include_amp,same_flav
     logical,dimension(:,:),allocatable :: sector_present
     logical,dimension(:,:),allocatable :: sector_retained
     logical :: lib_created=.false.,sectors_pruned_empty=.false.,sectors_pruned=.false.
   contains
     procedure,public :: init,evaluate,init_col,filter_helicity,write_init_amps_to_file,read_init_amps_from_file &
          ,create_library,optimise_evaluation,sector_index,prune_coupling_sectors
     procedure,private :: filter_dead_trees
     final :: finalize_amplitude_QCD ! custom deallocation of amplitude_QCD
  end type amplitude_QCD
contains
  subroutine init(this,imode,n,n_processes,part,spin,o,pm,max_as_ew_order)
    use math_functions
    use particles
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    type(physics_model),intent(in) :: pm
    integer,intent(in) :: n,imode,n_processes
    integer,dimension(n,n_processes),intent(in) :: part,o
    integer,dimension(0:3,n),intent(in) :: spin
    integer,optional,intent(in) :: max_as_ew_order
    integer,dimension(:,:,:),allocatable :: order
    type(current),dimension(:),allocatable :: current_list_local
    type(interaction),dimension(:),allocatable :: interaction_list_local
    integer :: isize,nc,isplit,n1,n2,ic1,ic2,max_cur,max_vert,max_key,ispin,iproc
    integer :: n_external_quarks,n_external_gluons,expected_ew_order
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:,:),allocatable :: key_to_current
    logical :: compact_max_as_build,compact_qcd_backbone,compact_inputs_valid
    integer :: compact_ew_order

    compact_max_as_build=present(max_as_ew_order)
    compact_qcd_backbone=.false.
    compact_ew_order=-1
    if (compact_max_as_build) compact_ew_order=max_as_ew_order
    if (compact_max_as_build .and. (imode.ne.1 .or. compact_ew_order.lt.0)) then
       write (*,*) 'The maximum-aS construction fast path requires imode=1 and a non-negative EW order'
       stop 1
    endif
    this%sectors_pruned_empty=.false.
    this%sectors_pruned=.false.
    if (allocated(this%sector_retained)) deallocate(this%sector_retained)
    if (imode.eq.1) then
       write (99,*) 'Initialising amplitude for:'
       write (99,*) '   - all polarisation/helicity configurations'
       write (99,*) '   - a single colour order'
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
    if (compact_max_as_build) then
       compact_inputs_valid=.not.any(this%n_qqbar.lt.0) .and. &
            .not.any(this%n_qqbar.gt.3) .and. &
            .not.any(.not.((abs(this%processes).ge.1 .and. &
            abs(this%processes).le.6) .or. this%processes.eq.21 .or. &
            this%processes.eq.22 .or. this%processes.eq.23 .or. &
            abs(this%processes).eq.24 .or. this%processes.eq.25 .or. &
            (abs(this%processes).ge.11 .and. abs(this%processes).le.16)))
       do iproc=1,this%nprocs
          n_external_quarks=count(abs(this%processes(:,iproc)).ge.1 .and. &
               abs(this%processes(:,iproc)).le.6)
          n_external_gluons=count(this%processes(:,iproc).eq.21)
          if (n_external_quarks.ne.2*this%n_qqbar(iproc)) &
               compact_inputs_valid=.false.
          if (this%n_qqbar(iproc).eq.0) then
             if (n_external_gluons.eq.n) then
                expected_ew_order=0
             else
                compact_inputs_valid=.false.
                expected_ew_order=-1
             endif
          else
             expected_ew_order=count(this%processes(:,iproc).eq.22 .or. &
                  this%processes(:,iproc).eq.23 .or. &
                  abs(this%processes(:,iproc)).eq.24 .or. &
                  this%processes(:,iproc).eq.25 .or. &
                  (abs(this%processes(:,iproc)).ge.11 .and. &
                  abs(this%processes(:,iproc)).le.16))
          endif
          if (expected_ew_order.ne.compact_ew_order) compact_inputs_valid=.false.
       enddo
       if ((any(this%n_qqbar.le.1) .and. any(this%n_qqbar.ge.2)) .or. &
            .not.compact_inputs_valid) then
          write (*,*) 'The maximum-aS construction fast path received an incompatible process group'
          stop 1
       endif
       compact_qcd_backbone=all(this%n_qqbar.ge.2)
    endif

    if (this%imode.eq.2) then
       call define_canonical_color_order()
    else
       this%nColOrd=1
    endif

    call set_max_cur()
    call set_max_vert()
    
    if (this%imode.eq.2) then
       if (this%nprocs.ne.1) then
          write (*,*) 'For imode==2 there should only be a single process at a time'
          stop 1
       endif
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
                do ispin=1, spin(0,order(nc,1,iproc))
                   call create_external_current(iproc,spin(ispin,order(nc,1,iproc)),&
                        this%processes(order(nc,1,iproc),iproc),order(nc,1,iproc))
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
    call remap_two_line_singlet_exchange_sectors()
    if (this%imode.ne.2 .and. .not.compact_max_as_build) &
         call filter_fixed_three_line_sector_terms()
    if (this%imode.eq.2 .and. this%n_qqbar(1).eq.3) then
       call group_three_line_colour_flows()
    endif
    call allocate_and_fill_momentum_array()

    ! All done. But there could be currents that are not needed. Filter them out
    write (99,*) 'Total number of currents, vertices and amplitudes before filter',this%n_cur,this%n_vert,this%n_amps
    call this%filter_dead_trees(n)
    write (99,*) 'Total number of currents, vertices and amplitudes after filter',this%n_cur,this%n_vert,this%n_amps

    call deallocate_unneeded()
  contains

    subroutine create_external_current(iproc,ispin,ipart,iorder)
      implicit none
      integer,intent(in) :: ispin,ipart,iorder,iproc
      integer :: ic,ctype,ichir

      if (iorder.le.2 .and. abs(ipart).le.6) then
         ctype=anti_current(ipart)
      else
         ctype=ipart
      endif
      ichir=external_chirality(ctype,iorder,ispin)
      do ic=1,this%n_cur
         if (current_list_local(ic)%order(1).ne.iorder) cycle
         if (current_list_local(ic)%type.ne.ctype) cycle
         if (current_list_local(ic)%chirality.ne.ichir) cycle
         if (current_list_local(ic)%bin.ne.ibset(0,iorder-1)) cycle
         if (current_list_local(ic)%spin(1).ne.ispin) cycle
         ! existing current.
         call current_list_local(ic)%iproc%set_bit(iproc)
         return
      enddo
      ! new external current
      this%n_cur=this%n_cur+1
      allocate(current_list_local(this%n_cur)%order(isize))
      current_list_local(this%n_cur)%order(1)=iorder
      current_list_local(this%n_cur)%type=ctype
      current_list_local(this%n_cur)%chirality=ichir
      current_list_local(this%n_cur)%mass=pm%get_mass(current_list_local(this%n_cur)%type)
      current_list_local(this%n_cur)%width=pm%get_width(current_list_local(this%n_cur)%type)
      allocate(current_list_local(this%n_cur)%ext_type(isize))
      current_list_local(this%n_cur)%ext_type(1)=current_list_local(this%n_cur)%type
      current_list_local(this%n_cur)%bin=ibset(0,iorder-1) ! give binary label
      allocate(current_list_local(this%n_cur)%spin(isize))
      current_list_local(this%n_cur)%spin(1)=ispin
      current_list_local(this%n_cur)%n_vert=0
      current_list_local(this%n_cur)%n_gs=0
      current_list_local(this%n_cur)%n_ew=0
      current_list_local(this%n_cur)%open_quark_leg=0
      current_list_local(this%n_cur)%ew_pairs=0
      current_list_local(this%n_cur)%u1_pairs=0
      current_list_local(this%n_cur)%gluon_pairs=0
      current_list_local(this%n_cur)%u1_links=0
      current_list_local(this%n_cur)%gluon_links=0
      if (pm%is_quark(ctype).or.pm%is_antiquark(ctype)) &
           current_list_local(this%n_cur)%open_quark_leg=iorder
      call current_list_local(this%n_cur)%iproc%init(this%nprocs)
      call current_list_local(this%n_cur)%iproc%set_bit(iproc)
      current_list_local(this%n_cur)%ext_cur=ibset(int(0,kind=16),this%n_cur-1)      
      ! create the lepton-related properties
!!$      allocate(current_list_local(this%n_cur)%fermi_list(1+nl))
!!$      current_list_local(this%n_cur)%fermi_list=0
!!$      if (pm%is_lepton(current_list_local(this%n_cur)%type)) then
!!$           current_list_local(this%n_cur)%fermi_list(1)=1
!!$           current_list_local(this%n_cur)%fermi_list(2)=ext_from_cur(this%n_cur)
!!$      elseif (pm%is_antilepton(current_list_local(this%n_cur)%type)) then
!!$           current_list_local(this%n_cur)%fermi_list(1)=1
!!$           current_list_local(this%n_cur)%fermi_list(3)=-ext_from_cur(this%n_cur)
!!$      endif
    end subroutine create_external_current
    
    subroutine allocate_and_fill_currents_to_amps_map()
      ! The 'curr2amp(1:2,iamp)' variable lists which two currents (one of
      ! size n-1 and one of size 1) result in the amplitude 'iamp'. This
      ! subroutine also sets the include_amp(iamp) to .true. for all
      ! amplitudes
      implicit none
      integer :: icur,jcur,i,j,iproc,iraw,iphys,isector,nraw,maxraw
      integer :: nphys_base,nterms,offset,slot
      type(bitset) :: proc
      integer,dimension(:,:),allocatable :: raw_curr2amp,raw_power,counts,cursor
      integer,dimension(:),allocatable :: raw_iproc,raw_phys,physical_representative
      integer,dimension(n,this%nprocs) :: procs
      do j=1,this%nprocs
         do i=1,n
            if (order(i,1,j).le.2) then
               procs(i,j)=anti_current(this%processes(order(i,1,j),j))
            else
               procs(i,j)=this%processes(order(i,1,j),j)
            endif
         enddo
      enddo
      maxraw=(this%n_cur_end(n-1)-this%n_cur_start(n-1)+1)*&
           (this%n_cur_end(n)-this%n_cur_start(n)+1)*this%nprocs
      allocate(raw_curr2amp(2,maxraw),raw_power(2,maxraw),raw_iproc(maxraw),raw_phys(maxraw))
      allocate(physical_representative(maxraw))
      raw_curr2amp=0
      raw_power=0
      raw_iproc=0
      raw_phys=0
      physical_representative=0
      allocate(this%iproc_start(1:this%nprocs+1))
      nraw=0
      this%n_amps=0
      do iproc=1,this%nprocs
         this%iproc_start(iproc)=this%n_amps+1
         if (this%same_flav(iproc)) cycle
         do icur=this%n_cur_start(n-1),this%n_cur_end(n-1)
            do jcur=this%n_cur_start(n),this%n_cur_end(n)
               if ( current_list_local(icur)%type .ne. anti_current(current_list_local(jcur)%type) ) cycle
               if (current_list_local(icur)%chirality.ne.0 .and. current_list_local(jcur)%chirality.ne.0) then
                  if (current_list_local(jcur)%chirality.ne.-current_list_local(icur)%chirality) cycle
               endif
               if (iand(current_list_local(icur)%bin,current_list_local(jcur)%bin).ne.0) cycle
               proc=current_list_local(icur)%iproc.and.current_list_local(jcur)%iproc
               if (proc%count_bits().eq.0) then
                  ! combination of icur and jcur does not contribute to any of the processes
                  cycle
               elseif (proc%count_bits().ne.1) then
                  write (*,*) 'A given amplitude should only contribute to one process',proc%count_bits()
                  do i=1,this%nprocs
                     if (proc%test_bit(i)) write (*,*) i
                  enddo
                  write (*,*) current_list_local(icur)%ext_type(1:n-1),'   , ',current_list_local(jcur)%ext_type(1)
                  stop 1
               elseif (.not.proc%test_bit(iproc)) then
                  ! one process, but it is not equal to process 'iproc'
                  cycle
               endif
               nraw=nraw+1
               raw_curr2amp(:,nraw)=[icur,jcur]
               raw_power(:,nraw)=[current_list_local(icur)%n_gs+current_list_local(jcur)%n_gs,&
                    current_list_local(icur)%n_ew+current_list_local(jcur)%n_ew]
               raw_iproc(nraw)=iproc
               if (sum(raw_power(:,nraw)).ne.n-2) then
                  write (*,*) 'Inconsistent terminal coupling order',raw_power(:,nraw),n-2
                  write (*,*) this%processes(:,iproc)
                  stop 1
               endif
               do iphys=this%iproc_start(iproc),this%n_amps
                  if (same_physical_closure(raw_curr2amp(:,nraw),&
                       raw_curr2amp(:,physical_representative(iphys)),iproc)) exit
               enddo
               if (iphys.gt.this%n_amps) then
                  this%n_amps=this%n_amps+1
                  physical_representative(this%n_amps)=nraw
                  iphys=this%n_amps
               endif
               raw_phys(nraw)=iphys
            enddo
         enddo
      enddo

      call build_sector_list(nraw,raw_power)
      nphys_base=this%n_amps
      if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. this%imode.eq.2) then
         this%n_amps=this%n_amps*2
      endif
      allocate(this%curr2amp(2,this%n_amps))
      allocate(this%sector_curr2amp(2,this%n_amps,this%n_sectors))
      allocate(this%sector_present(this%n_amps,this%n_sectors))
      allocate(this%sector_sign(this%n_amps,this%n_sectors))
      allocate(counts(this%n_amps,this%n_sectors))
      this%curr2amp=0
      this%sector_curr2amp=0
      this%sector_present=.false.
      this%sector_sign=1
      counts=0
      do iraw=1,nraw
         isector=find_sector(raw_power(:,iraw))
         iphys=raw_phys(iraw)
         counts(iphys,isector)=counts(iphys,isector)+1
      enddo
      if (this%n_amps.gt.nphys_base) then
         counts(nphys_base+1:this%n_amps,:)=counts(1:nphys_base,:)
      endif
      nterms=sum(counts)
      allocate(this%sector_term_start(0:this%n_amps,this%n_sectors))
      allocate(this%sector_term_curr2amp(2,nterms))
      allocate(this%sector_term_sign(nterms))
      allocate(cursor(this%n_amps,this%n_sectors))
      offset=0
      do isector=1,this%n_sectors
         this%sector_term_start(0,isector)=offset
         do iphys=1,this%n_amps
            offset=offset+counts(iphys,isector)
            this%sector_term_start(iphys,isector)=offset
            cursor(iphys,isector)=offset-counts(iphys,isector)
         enddo
      enddo
      do iraw=1,nraw
         isector=find_sector(raw_power(:,iraw))
         iphys=raw_phys(iraw)
         cursor(iphys,isector)=cursor(iphys,isector)+1
         slot=cursor(iphys,isector)
         this%sector_term_curr2amp(:,slot)=raw_curr2amp(:,iraw)
         this%sector_term_sign(slot)=1
         if (this%n_amps.gt.nphys_base) then
            iphys=iphys+nphys_base
            cursor(iphys,isector)=cursor(iphys,isector)+1
            slot=cursor(iphys,isector)
            this%sector_term_curr2amp(:,slot)=raw_curr2amp(:,iraw)
            this%sector_term_sign(slot)=1
         endif
      enddo
      do iphys=1,this%n_amps
         do isector=1,this%n_sectors
            if (counts(iphys,isector).eq.0) cycle
            slot=this%sector_term_start(iphys-1,isector)+1
            this%sector_curr2amp(:,iphys,isector)=this%sector_term_curr2amp(:,slot)
            this%sector_sign(iphys,isector)=this%sector_term_sign(slot)
            this%sector_present(iphys,isector)=.true.
            if (all(this%curr2amp(:,iphys).eq.0)) &
                 this%curr2amp(:,iphys)=this%sector_curr2amp(:,iphys,isector)
         enddo
      enddo
      if (this%imode.eq.3 .and. this%n_amps.ne.1) then
         write (*,*) 'For this%imode==3, there should only be one amplitude',this%n_amps
         write (*,*) this%n_cur_start
         write (*,*) this%n_cur_end
         stop 1
      endif

      allocate(this%same_flavour_sum(this%n_amps,2))
      allocate(this%same_flavour_sum_operation(this%n_amps,2))
      this%same_flavour_sum=-1
      this%iproc_start(this%nprocs+1)=this%n_amps+1

      allocate(this%include_amp(1:this%n_amps))
      this%include_amp(:)=.true.
      deallocate(raw_curr2amp,raw_power,raw_iproc,raw_phys,physical_representative,&
           counts,cursor)
    end subroutine allocate_and_fill_currents_to_amps_map

    logical function same_physical_closure(left,right,iproc)
        implicit none
        integer,dimension(2),intent(in) :: left,right
        integer,intent(in) :: iproc
        integer :: il,ir
        il=left(1)
        ir=right(1)
        same_physical_closure=.false.
        if (this%imode.ne.2 .and. this%n_qqbar(iproc).ge.2) then
           same_physical_closure=same_terminal_helicity(left,right)
           return
        endif
        if (left(2).ne.right(2)) return
        if (current_list_local(il)%type.ne.current_list_local(ir)%type) return
        if (current_list_local(il)%chirality.ne.current_list_local(ir)%chirality) return
        if (current_list_local(il)%bin.ne.current_list_local(ir)%bin) return
        if (current_list_local(il)%ext_cur.ne.current_list_local(ir)%ext_cur) return
        ! Coupling order and colour-singlet closure history select terms of
        ! one physical ordered coefficient; they are not extra amplitudes.
        ! The sparse sector-term map retains those distinct roots.
        if (any(current_list_local(il)%order(1:n-1).ne.&
             current_list_local(ir)%order(1:n-1))) return
        same_physical_closure=.true.
      end function same_physical_closure

      logical function same_terminal_helicity(left,right)
        implicit none
        integer,dimension(2),intent(in) :: left,right
        integer :: side,pos,leg
        integer,dimension(n,2) :: helicity
        integer,dimension(2,2) :: roots
        roots(:,1)=left
        roots(:,2)=right
        helicity=999
        do side=1,2
           do pos=1,n-1
              leg=current_list_local(roots(1,side))%order(pos)
              helicity(leg,side)=current_list_local(roots(1,side))%spin(pos)
           enddo
           leg=current_list_local(roots(2,side))%order(1)
           helicity(leg,side)=current_list_local(roots(2,side))%spin(1)
        enddo
        same_terminal_helicity=all(helicity(:,1).eq.helicity(:,2))
      end function same_terminal_helicity

      subroutine build_sector_list(nentries,powers)
        implicit none
        integer,intent(in) :: nentries
        integer,dimension(:,:),intent(in) :: powers
        integer :: entry,sector,left,right
        integer,dimension(2) :: candidate
        integer,dimension(:,:),allocatable :: unique_powers
        allocate(unique_powers(2,max(1,nentries)))
        this%n_sectors=0
        do entry=1,nentries
           candidate=powers(:,entry)
           do sector=1,this%n_sectors
              if (all(unique_powers(:,sector).eq.candidate)) exit
           enddo
           if (sector.le.this%n_sectors) cycle
           this%n_sectors=this%n_sectors+1
           unique_powers(:,this%n_sectors)=candidate
        enddo
        ! Stable, deterministic order: leading QCD first, then lower EW.
        do left=1,this%n_sectors-1
           do right=left+1,this%n_sectors
              if (unique_powers(1,right).gt.unique_powers(1,left) .or. &
                   (unique_powers(1,right).eq.unique_powers(1,left) .and. &
                    unique_powers(2,right).lt.unique_powers(2,left))) then
                 candidate=unique_powers(:,left)
                 unique_powers(:,left)=unique_powers(:,right)
                 unique_powers(:,right)=candidate
              endif
           enddo
        enddo
        allocate(this%sector_powers(2,this%n_sectors))
        this%sector_powers=unique_powers(:,1:this%n_sectors)
        deallocate(unique_powers)
      end subroutine build_sector_list

      integer function find_sector(power)
        implicit none
        integer,dimension(2),intent(in) :: power
        do find_sector=1,this%n_sectors
           if (all(this%sector_powers(:,find_sector).eq.power)) return
        enddo
        write (*,*) 'Internal error: unknown amplitude coupling sector',power
        stop 1
      end function find_sector

    subroutine allocate_and_fill_spins()
      implicit none
      integer :: iamp,i,iproc
      allocate(this%spins(n,1,1:this%n_amps))
      do iproc=1,this%nprocs
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            do i=this%n_cur_start(1),this%n_cur_end(1)
               if (btest(this%current_list(this%curr2amp(1,iamp))%ext_cur,i-1)) then
                  this%spins(this%current_list(i)%order(1),1,iamp)=&
                       this%current_list(i)%spin(1)
               endif
            enddo
            this%spins(order(n,1,iproc),1,iamp)=this&
                 %current_list(this%curr2amp(2,iamp))%spin(1)
         enddo
      enddo
    end subroutine allocate_and_fill_spins
    
    subroutine allocate_and_fill_colour_permutations()
      implicit none
      integer :: iamp,i,iproc,iamp_to_compare,pos
      ! allocate and fill the colour orders in 'this%perm'. These are simply
      ! the orders of the elements in the 'this%current_list' (with size n-1)
      ! together with the final element). Exception: when there are colour
      ! singlets, they will not be part of the this%perm (while they are part
      ! of the elements in the this%current_list.
      allocate(this%perm(1:n-this%n_sing(1),1:this%n_amps))
      if (pm%is_singlet(this%current_list(this%n_cur_start(n))%type)) then
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
            if (this%imode.ne.2 .and. this%n_qqbar(iproc).ge.2) then
               pos=0
               do i=1,n
                  if (pm%is_singlet(this%processes(order(i,1,iproc),iproc))) cycle
                  pos=pos+1
                  this%perm(pos,iamp)=order(i,1,iproc)
               enddo
               cycle
            endif
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
         enddo
      enddo
    end subroutine allocate_and_fill_colour_permutations

    subroutine remap_two_line_singlet_exchange_sectors()
      ! A colour-singlet exchange connects each quark to the antiquark on the
      ! same fermion line.  The recursive planar order labels that coefficient
      ! by the complementary open-string flow used for a gluon exchange.  Move
      ! it to the physical delta-flow row and include the fermion-order sign.
      implicit none
      integer :: iproc,isector,source,target,first_amp,last_amp,nord,iterm,&
           new_term,nterms,term_sign,offset
      integer,dimension(:),allocatable :: transformed
      integer,dimension(:),allocatable :: old_term_sign,term_target,term_sector
      integer,dimension(:,:),allocatable :: old_term_start,old_term_curr2amp,&
           counts,cursor

      if (this%imode.ne.2) then
         ! Compact QCD-backbone amplitudes have already excluded the
         ! colour-singlet quark-line closure whose endpoint exchange this
         ! routine resolves.  They use only the prescribed physical planar
         ! order, so there is nothing left to remap or filter here.
         if (.not.compact_max_as_build) call filter_fixed_two_line_sector_terms()
         return
      endif
      if (all(this%n_qqbar(1:this%nprocs).ne.2)) return
      allocate(old_term_start(0:this%n_amps,this%n_sectors))
      allocate(old_term_curr2amp(2,size(this%sector_term_sign)))
      allocate(old_term_sign(size(this%sector_term_sign)))
      old_term_start=this%sector_term_start
      old_term_curr2amp=this%sector_term_curr2amp
      old_term_sign=this%sector_term_sign
      nterms=size(old_term_sign)
      allocate(term_target(max(1,nterms)),term_sector(max(1,nterms)))
      allocate(counts(this%n_amps,this%n_sectors))
      counts=0
      new_term=0
      do iproc=1,this%nprocs
         if (this%n_qqbar(iproc).ne.2) cycle
         nord=n-this%n_sing(iproc)
         allocate(transformed(nord))
         first_amp=this%iproc_start(iproc)
         last_amp=this%iproc_start(iproc+1)-1
         do isector=1,this%n_sectors
            do source=first_amp,last_amp
               do iterm=old_term_start(source-1,isector)+1,&
                    old_term_start(source,isector)
                  call transform_quark_line_flow(this%perm(:,source),&
                       this%current_list(old_term_curr2amp(1,iterm))%ew_pairs,&
                       iproc,transformed,term_sign)
                  do target=first_amp,last_amp
                     if (all(this%perm(:,target).eq.transformed)) exit
                  enddo
                  if (target.gt.last_amp) then
                     write (*,*) 'Cannot find colour-singlet exchange partner flow',&
                          this%perm(:,source),transformed
                     stop 1
                  endif
                  new_term=new_term+1
                  term_target(new_term)=target
                  term_sector(new_term)=isector
                  old_term_sign(iterm)=old_term_sign(iterm)*term_sign
                  counts(target,isector)=counts(target,isector)+1
               enddo
            enddo
         enddo
         deallocate(transformed)
      enddo
      if (new_term.ne.nterms) then
         write (*,*) 'Internal error while remapping two-line sector terms',&
              new_term,nterms
         stop 1
      endif
      deallocate(this%sector_term_start,this%sector_term_curr2amp,&
           this%sector_term_sign)
      allocate(this%sector_term_start(0:this%n_amps,this%n_sectors))
      allocate(this%sector_term_curr2amp(2,nterms))
      allocate(this%sector_term_sign(nterms))
      allocate(cursor(this%n_amps,this%n_sectors))
      offset=0
      do isector=1,this%n_sectors
         this%sector_term_start(0,isector)=offset
         do target=1,this%n_amps
            offset=offset+counts(target,isector)
            this%sector_term_start(target,isector)=offset
            cursor(target,isector)=offset-counts(target,isector)
         enddo
      enddo
      new_term=0
      do isector=1,this%n_sectors
         do source=1,this%n_amps
            do iterm=old_term_start(source-1,isector)+1,&
                 old_term_start(source,isector)
               new_term=new_term+1
               target=term_target(new_term)
               cursor(target,isector)=cursor(target,isector)+1
               offset=cursor(target,isector)
               this%sector_term_curr2amp(:,offset)=old_term_curr2amp(:,iterm)
               this%sector_term_sign(offset)=old_term_sign(iterm)
            enddo
         enddo
      enddo
      this%sector_curr2amp=0
      this%sector_present=.false.
      this%sector_sign=1
      this%curr2amp=0
      do target=1,this%n_amps
         do isector=1,this%n_sectors
            if (counts(target,isector).eq.0) cycle
            iterm=this%sector_term_start(target-1,isector)+1
            this%sector_curr2amp(:,target,isector)=&
                 this%sector_term_curr2amp(:,iterm)
            this%sector_sign(target,isector)=this%sector_term_sign(iterm)
            this%sector_present(target,isector)=.true.
            if (all(this%curr2amp(:,target).eq.0)) &
                 this%curr2amp(:,target)=this%sector_curr2amp(:,target,isector)
         enddo
      enddo
      deallocate(old_term_start,old_term_curr2amp,old_term_sign,term_target,&
           term_sector,counts,cursor)
    end subroutine remap_two_line_singlet_exchange_sectors

    subroutine filter_fixed_two_line_sector_terms()
      ! A fixed physical colour flow is assembled from two raw planar
      ! representations.  Keep, for each helicity and coupling sector, only
      ! the roots whose colour-singlet endpoint transformation lands on the
      ! prescribed flow.  This aligns QCD and EW sectors before squaring.
      implicit none
      integer :: source,isector,iterm,iproc,nord,nold,nkeep,offset,slot
      integer :: root,term_sign
      integer,dimension(:),allocatable :: keep_source,keep_sector,keep_sign
      integer,dimension(:,:),allocatable :: keep_curr2amp,counts,cursor
      integer,dimension(:),allocatable :: raw_word,transformed

      if (all(this%n_qqbar.ne.2)) return
      nold=size(this%sector_term_sign)
      allocate(keep_source(max(1,nold)),keep_sector(max(1,nold)),&
           keep_sign(max(1,nold)),keep_curr2amp(2,max(1,nold)))
      allocate(counts(this%n_amps,this%n_sectors))
      counts=0
      nkeep=0
      iproc=1
      do source=1,this%n_amps
         do while (source.ge.this%iproc_start(iproc+1) .and. iproc.lt.this%nprocs)
            iproc=iproc+1
         enddo
         do isector=1,this%n_sectors
            do iterm=this%sector_term_start(source-1,isector)+1,&
                 this%sector_term_start(source,isector)
               if (this%n_qqbar(iproc).eq.2) then
                  nord=n-this%n_sing(iproc)
                  allocate(raw_word(nord),transformed(nord))
                  root=this%sector_term_curr2amp(1,iterm)
                  raw_word(1:nord-1)=this%current_list(root)%order(1:nord-1)
                  raw_word(nord)=this%current_list(&
                       this%sector_term_curr2amp(2,iterm))%order(1)
                  call transform_quark_line_flow(raw_word,&
                       this%current_list(root)%ew_pairs,iproc,transformed,term_sign)
                  if (.not.same_two_line_flow(transformed,this%perm(:,source),iproc)) then
                     deallocate(raw_word,transformed)
                     cycle
                  endif
                  deallocate(raw_word,transformed)
               else
                  term_sign=1
               endif
               nkeep=nkeep+1
               keep_source(nkeep)=source
               keep_sector(nkeep)=isector
               keep_sign(nkeep)=this%sector_term_sign(iterm)*term_sign
               keep_curr2amp(:,nkeep)=this%sector_term_curr2amp(:,iterm)
               counts(source,isector)=counts(source,isector)+1
            enddo
         enddo
      enddo

      deallocate(this%sector_term_start,this%sector_term_curr2amp,&
           this%sector_term_sign)
      allocate(this%sector_term_start(0:this%n_amps,this%n_sectors))
      allocate(this%sector_term_curr2amp(2,nkeep))
      allocate(this%sector_term_sign(nkeep))
      allocate(cursor(this%n_amps,this%n_sectors))
      offset=0
      do isector=1,this%n_sectors
         this%sector_term_start(0,isector)=offset
         do source=1,this%n_amps
            offset=offset+counts(source,isector)
            this%sector_term_start(source,isector)=offset
            cursor(source,isector)=offset-counts(source,isector)
         enddo
      enddo
      do iterm=1,nkeep
         source=keep_source(iterm)
         isector=keep_sector(iterm)
         cursor(source,isector)=cursor(source,isector)+1
         slot=cursor(source,isector)
         this%sector_term_curr2amp(:,slot)=keep_curr2amp(:,iterm)
         this%sector_term_sign(slot)=keep_sign(iterm)
      enddo
      this%curr2amp=0
      this%sector_curr2amp=0
      this%sector_present=.false.
      this%sector_sign=1
      do source=1,this%n_amps
         do isector=1,this%n_sectors
            if (counts(source,isector).eq.0) cycle
            iterm=this%sector_term_start(source-1,isector)+1
            this%sector_curr2amp(:,source,isector)=this%sector_term_curr2amp(:,iterm)
            this%sector_sign(source,isector)=this%sector_term_sign(iterm)
            this%sector_present(source,isector)=.true.
            if (all(this%curr2amp(:,source).eq.0)) &
                 this%curr2amp(:,source)=this%sector_curr2amp(:,source,isector)
         enddo
      enddo
      deallocate(keep_source,keep_sector,keep_sign,keep_curr2amp,counts,cursor)
    end subroutine filter_fixed_two_line_sector_terms

    subroutine filter_fixed_three_line_sector_terms()
      ! The recursive construction admits all terminal-preserving raw
      ! three-string representations.  Transform each sparse root using its
      ! actual singlet/Fierz closure history and retain only terms belonging
      ! to the one physical flow prescribed by this fixed-order amplitude.
      implicit none
      integer :: source,isector,iterm,iproc,nord,nold,nkeep,offset,slot
      integer :: root,closing,term_sign
      integer,dimension(:),allocatable :: keep_source,keep_sector,keep_sign
      integer,dimension(:,:),allocatable :: keep_curr2amp,counts,cursor
      integer,dimension(:),allocatable :: raw_word,transformed

      if (all(this%n_qqbar.ne.3)) return
      nold=size(this%sector_term_sign)
      allocate(keep_source(max(1,nold)),keep_sector(max(1,nold)),&
           keep_sign(max(1,nold)),keep_curr2amp(2,max(1,nold)))
      allocate(counts(this%n_amps,this%n_sectors))
      counts=0
      nkeep=0
      do iproc=1,this%nprocs
         do source=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            do isector=1,this%n_sectors
               do iterm=this%sector_term_start(source-1,isector)+1,&
                    this%sector_term_start(source,isector)
                  term_sign=1
                  if (this%n_qqbar(iproc).eq.3) then
                     nord=n-this%n_sing(iproc)
                     allocate(raw_word(nord),transformed(nord))
                     root=this%sector_term_curr2amp(1,iterm)
                     closing=this%sector_term_curr2amp(2,iterm)
                     raw_word(1:nord-1)=&
                          this%current_list(root)%order(1:nord-1)
                     raw_word(nord)=this%current_list(closing)%order(1)
                     if (any(this%current_list(root)%ew_pairs.ne.0) .and. &
                          any(this%current_list(root)%gluon_links.ne.0)) then
                        call construct_three_line_flow(raw_word,&
                             this%current_list(root)%ew_pairs,&
                             this%current_list(root)%gluon_pairs,&
                             this%current_list(root)%gluon_links,&
                             this%current_list(root)%open_quark_leg,&
                             this%current_list(closing)%open_quark_leg,&
                             iproc,transformed,term_sign)
                     else
                        call transform_quark_line_flow(raw_word,&
                             this%current_list(root)%ew_pairs,iproc,&
                             transformed,term_sign)
                     endif
                     if (.not.same_three_line_flow(transformed,&
                          this%perm(:,source),iproc)) then
                        deallocate(raw_word,transformed)
                        cycle
                     endif
                     deallocate(raw_word,transformed)
                  endif
                  nkeep=nkeep+1
                  keep_source(nkeep)=source
                  keep_sector(nkeep)=isector
                  keep_sign(nkeep)=this%sector_term_sign(iterm)*term_sign
                  keep_curr2amp(:,nkeep)=this%sector_term_curr2amp(:,iterm)
                  counts(source,isector)=counts(source,isector)+1
               enddo
            enddo
         enddo
      enddo

      deallocate(this%sector_term_start,this%sector_term_curr2amp,&
           this%sector_term_sign)
      allocate(this%sector_term_start(0:this%n_amps,this%n_sectors))
      allocate(this%sector_term_curr2amp(2,nkeep))
      allocate(this%sector_term_sign(nkeep))
      allocate(cursor(this%n_amps,this%n_sectors))
      offset=0
      do isector=1,this%n_sectors
         this%sector_term_start(0,isector)=offset
         do source=1,this%n_amps
            offset=offset+counts(source,isector)
            this%sector_term_start(source,isector)=offset
            cursor(source,isector)=offset-counts(source,isector)
         enddo
      enddo
      do iterm=1,nkeep
         source=keep_source(iterm)
         isector=keep_sector(iterm)
         cursor(source,isector)=cursor(source,isector)+1
         slot=cursor(source,isector)
         this%sector_term_curr2amp(:,slot)=keep_curr2amp(:,iterm)
         this%sector_term_sign(slot)=keep_sign(iterm)
      enddo
      this%curr2amp=0
      this%sector_curr2amp=0
      this%sector_present=.false.
      this%sector_sign=1
      do source=1,this%n_amps
         do isector=1,this%n_sectors
            if (counts(source,isector).eq.0) cycle
            iterm=this%sector_term_start(source-1,isector)+1
            this%sector_curr2amp(:,source,isector)=&
                 this%sector_term_curr2amp(:,iterm)
            this%sector_sign(source,isector)=this%sector_term_sign(iterm)
            this%sector_present(source,isector)=.true.
            if (all(this%curr2amp(:,source).eq.0)) &
                 this%curr2amp(:,source)=this%sector_curr2amp(:,source,isector)
         enddo
      enddo
      deallocate(keep_source,keep_sector,keep_sign,keep_curr2amp,counts,cursor)
    end subroutine filter_fixed_three_line_sector_terms

    logical function same_two_line_flow(left,right,iproc)
      implicit none
      integer,dimension(:),intent(in) :: left,right
      integer,intent(in) :: iproc
      integer :: line,candidate,other
      integer,dimension(2) :: left_len,right_len
      integer,dimension(n,2) :: left_lines,right_lines
      call split_two_line_flow(left,iproc,left_len,left_lines)
      call split_two_line_flow(right,iproc,right_len,right_lines)
      same_two_line_flow=.true.
      do line=1,2
         other=0
         do candidate=1,2
            if (right_lines(1,candidate).eq.left_lines(1,line)) then
               other=candidate
               exit
            endif
         enddo
         if (other.eq.0 .or. left_len(line).ne.right_len(other)) then
            same_two_line_flow=.false.
            return
         endif
         if (any(left_lines(1:left_len(line),line).ne.&
              right_lines(1:right_len(other),other))) then
            same_two_line_flow=.false.
            return
         endif
      enddo
    end function same_two_line_flow

    subroutine split_two_line_flow(word,iproc,line_len,lines)
      implicit none
      integer,dimension(:),intent(in) :: word
      integer,intent(in) :: iproc
      integer,dimension(2),intent(out) :: line_len
      integer,dimension(n,2),intent(out) :: lines
      integer :: pos,line
      line_len=0
      lines=0
      line=0
      do pos=1,size(word)
         if (is_quark_from_order(word(pos),iproc)) line=line+1
         if (line.lt.1 .or. line.gt.2) then
            write (*,*) 'Malformed two-line colour flow',word
            stop 1
         endif
         line_len(line)=line_len(line)+1
         lines(line_len(line),line)=word(pos)
      enddo
      if (line.ne.2) then
         write (*,*) 'Incomplete two-line colour flow',word
         stop 1
      endif
    end subroutine split_two_line_flow

    subroutine group_three_line_colour_flows()
      ! Collapse the alternative recursive closures onto the physical
      ! three-open-string basis.  Coupling sectors can contain several roots
      ! for one physical row; keep those roots in a compact sparse term map.
      implicit none
      integer :: old_amp,flow,flow_count,old_namps,isector,target,nterms,iterm,&
           old_term
      integer :: offset,root,term_sign
      integer,dimension(:),allocatable :: representative,flow_of_old,term_target,term_sector
      integer,dimension(:),allocatable :: term_signs,cursor
      integer,dimension(:,:),allocatable :: old_curr2amp,old_perm,&
           old_same_flavour_sum,old_same_flavour_sum_operation,old_sector_sign
      integer,dimension(:,:),allocatable :: old_term_start,old_term_curr2amp
      integer,dimension(:),allocatable :: old_term_sign
      integer,dimension(:,:),allocatable :: term_curr2amp,counts
      integer,dimension(:,:,:),allocatable :: old_sector_curr2amp
      logical,dimension(:),allocatable :: old_include_amp
      logical,dimension(:,:),allocatable :: old_sector_present
      integer,dimension(size(this%perm,1)) :: transformed

      old_namps=this%n_amps
      allocate(representative(old_namps),flow_of_old(old_namps))
      representative=0
      flow_of_old=0
      flow_count=0
      do old_amp=1,old_namps
         do flow=1,flow_count
            if (same_three_line_flow(this%perm(:,old_amp),&
                 this%perm(:,representative(flow)))) exit
         enddo
         if (flow.gt.flow_count) then
            flow_count=flow_count+1
            representative(flow_count)=old_amp
         endif
         flow_of_old(old_amp)=flow
      enddo
      if (flow_count.ne.this%nColOrd) then
         write (*,*) 'Unexpected three-line colour-flow closure multiplicity',&
              flow_count,this%nColOrd,old_namps
         stop 1
      endif

      allocate(old_curr2amp(2,old_namps))
      allocate(old_perm(size(this%perm,1),old_namps))
      allocate(old_include_amp(old_namps))
      allocate(old_same_flavour_sum(old_namps,2))
      allocate(old_same_flavour_sum_operation(old_namps,2))
      allocate(old_sector_curr2amp(2,old_namps,this%n_sectors))
      allocate(old_sector_present(old_namps,this%n_sectors))
      allocate(old_sector_sign(old_namps,this%n_sectors))
      allocate(old_term_start(0:old_namps,this%n_sectors))
      allocate(old_term_curr2amp(2,size(this%sector_term_sign)))
      allocate(old_term_sign(size(this%sector_term_sign)))
      old_curr2amp=this%curr2amp
      old_perm=this%perm
      old_include_amp=this%include_amp
      old_same_flavour_sum=this%same_flavour_sum
      old_same_flavour_sum_operation=this%same_flavour_sum_operation
      old_sector_curr2amp=this%sector_curr2amp
      old_sector_present=this%sector_present
      old_sector_sign=this%sector_sign
      old_term_start=this%sector_term_start
      old_term_curr2amp=this%sector_term_curr2amp
      old_term_sign=this%sector_term_sign

      nterms=size(old_term_sign)
      allocate(term_target(max(1,nterms)),term_sector(max(1,nterms)))
      allocate(term_signs(max(1,nterms)),term_curr2amp(2,max(1,nterms)))
      allocate(counts(this%nColOrd,this%n_sectors))
      counts=0
      iterm=0
      do old_amp=1,old_namps
         do isector=1,this%n_sectors
            do old_term=old_term_start(old_amp-1,isector)+1,&
                 old_term_start(old_amp,isector)
               root=old_term_curr2amp(1,old_term)
               if (any(this%current_list(root)%ew_pairs.ne.0) .and. &
                    any(this%current_list(root)%gluon_links.ne.0)) then
                  call construct_three_line_flow(old_perm(:,old_amp),&
                       this%current_list(root)%ew_pairs,&
                       this%current_list(root)%gluon_pairs,&
                       this%current_list(root)%gluon_links,&
                       this%current_list(root)%open_quark_leg,&
                       this%current_list(old_term_curr2amp(2,old_term))%open_quark_leg,&
                       1,transformed,term_sign)
               else
                  call transform_quark_line_flow(old_perm(:,old_amp),&
                       this%current_list(root)%ew_pairs,1,transformed,term_sign)
               endif
               do target=1,this%nColOrd
                  if (same_three_line_flow(transformed,&
                       old_perm(:,representative(target)))) exit
               enddo
               if (target.gt.this%nColOrd) then
                  write (*,*) 'Cannot map coupling-sector root to a three-line colour flow',&
                       old_perm(:,old_amp),transformed,this%current_list(root)%ew_pairs
                  stop 1
               endif
               iterm=iterm+1
               term_target(iterm)=target
               term_sector(iterm)=isector
               term_signs(iterm)=term_sign*old_term_sign(old_term)
               term_curr2amp(:,iterm)=old_term_curr2amp(:,old_term)
               counts(target,isector)=counts(target,isector)+1
            enddo
         enddo
      enddo

      deallocate(this%curr2amp,this%perm,this%include_amp,this%same_flavour_sum,&
           this%same_flavour_sum_operation,this%sector_curr2amp,this%sector_present,&
           this%sector_sign,this%sector_term_start,this%sector_term_curr2amp,&
           this%sector_term_sign)
      allocate(this%curr2amp(2,this%nColOrd))
      allocate(this%perm(size(old_perm,1),this%nColOrd))
      allocate(this%include_amp(this%nColOrd))
      allocate(this%same_flavour_sum(this%nColOrd,2))
      allocate(this%same_flavour_sum_operation(this%nColOrd,2))
      allocate(this%sector_curr2amp(2,this%nColOrd,this%n_sectors))
      allocate(this%sector_present(this%nColOrd,this%n_sectors))
      allocate(this%sector_sign(this%nColOrd,this%n_sectors))
      allocate(this%sector_term_start(0:this%nColOrd,this%n_sectors))
      allocate(this%sector_term_curr2amp(2,nterms))
      allocate(this%sector_term_sign(nterms))
      this%curr2amp=0
      this%sector_curr2amp=0
      this%sector_present=.false.
      this%sector_sign=1
      offset=0
      do isector=1,this%n_sectors
         this%sector_term_start(0,isector)=offset
         do flow=1,this%nColOrd
            offset=offset+counts(flow,isector)
            this%sector_term_start(flow,isector)=offset
         enddo
      enddo
      allocate(cursor(this%nColOrd*this%n_sectors))
      do isector=1,this%n_sectors
         do flow=1,this%nColOrd
            cursor((isector-1)*this%nColOrd+flow)=&
                 this%sector_term_start(flow-1,isector)
         enddo
      enddo
      do iterm=1,nterms
         target=term_target(iterm)
         isector=term_sector(iterm)
         flow=(isector-1)*this%nColOrd+target
         cursor(flow)=cursor(flow)+1
         this%sector_term_curr2amp(:,cursor(flow))=term_curr2amp(:,iterm)
         this%sector_term_sign(cursor(flow))=term_signs(iterm)
      enddo
      do flow=1,this%nColOrd
         this%perm(:,flow)=old_perm(:,representative(flow))
         this%include_amp(flow)=old_include_amp(representative(flow))
         this%same_flavour_sum(flow,:)=&
              old_same_flavour_sum(representative(flow),:)
         this%same_flavour_sum_operation(flow,:)=&
              old_same_flavour_sum_operation(representative(flow),:)
         do isector=1,this%n_sectors
            if (counts(flow,isector).eq.0) cycle
            iterm=this%sector_term_start(flow-1,isector)+1
            this%sector_curr2amp(:,flow,isector)=this%sector_term_curr2amp(:,iterm)
            this%sector_sign(flow,isector)=this%sector_term_sign(iterm)
            this%sector_present(flow,isector)=.true.
            if (all(this%curr2amp(:,flow).eq.0)) &
                 this%curr2amp(:,flow)=this%sector_curr2amp(:,flow,isector)
         enddo
      enddo
      this%n_amps=this%nColOrd
      this%iproc_start(this%nprocs+1)=this%n_amps+1
      deallocate(representative,flow_of_old,term_target,term_sector,term_signs,&
           term_curr2amp,counts,cursor,old_curr2amp,old_perm,old_include_amp,&
           old_same_flavour_sum,old_same_flavour_sum_operation,old_sector_curr2amp,&
           old_sector_present,old_sector_sign,old_term_start,old_term_curr2amp,&
           old_term_sign)
    end subroutine group_three_line_colour_flows

    subroutine construct_three_line_flow(word,ew_pairs,gluon_pairs,gluon_links,&
         root_leg,closing_leg,iproc,transformed,flow_sign)
      ! Reconstruct a three-quark-line colour tensor from its actual fermion
      ! closures.  A type-9 q-qbar closure supplies the U(N) part of the
      ! Fierz identity and therefore crosses the closed line with the line to
      ! which that gluon is attached.  The recursive planar word alone is not
      ! sufficient once a colour-singlet exchange has changed the ordering.
      implicit none
      integer,dimension(:),intent(in) :: word
      integer,dimension(3),intent(in) :: ew_pairs,gluon_pairs,gluon_links
      integer,intent(in) :: root_leg,closing_leg,iproc
      integer,dimension(size(word)),intent(out) :: transformed
      integer,intent(out) :: flow_sign
      integer :: i,pos,line,nlines,pair_code,first,second,qleg,aleg,&
           attach,attach_q,temp
      integer,dimension(3) :: qlabel,aq_for_q

      transformed=word
      qlabel=0
      aq_for_q=0
      nlines=0
      do pos=1,size(word)
         if (is_quark_from_order(word(pos),iproc)) then
            nlines=nlines+1
            if (nlines.gt.3) then
               write (*,*) 'Too many quark lines while reconstructing colour flow',word
               stop 1
            endif
            qlabel(nlines)=word(pos)
         endif
      enddo
      if (nlines.ne.3) then
         write (*,*) 'Incomplete quark lines while reconstructing colour flow',word
         stop 1
      endif

      do i=1,3
         if (ew_pairs(i).ne.0) call add_three_line_baseline_pair(&
              abs(ew_pairs(i)),iproc,qlabel,aq_for_q,word)
         if (gluon_pairs(i).ne.0) call add_three_line_baseline_pair(&
              gluon_pairs(i),iproc,qlabel,aq_for_q,word)
      enddo
      if (root_leg.ne.0 .and. closing_leg.ne.0) then
         pair_code=min(root_leg,closing_leg)*(n+1)+max(root_leg,closing_leg)
         call add_three_line_baseline_pair(pair_code,iproc,qlabel,aq_for_q,word)
      endif
      if (any(aq_for_q.eq.0)) then
         write (*,*) 'Incomplete fermion-pair history for three-line colour flow',&
              word,ew_pairs,gluon_pairs,root_leg,closing_leg,aq_for_q
         stop 1
      endif

      do i=1,3
         if (gluon_links(i).eq.0) cycle
         pair_code=gluon_links(i)/(n+1)
         attach=mod(gluon_links(i),n+1)
         call decode_three_line_pair(pair_code,iproc,qleg,aleg)
         if (is_quark_from_order(attach,iproc)) then
            attach_q=attach
         else
            attach_q=0
            do line=1,3
               if (aq_for_q(line).eq.attach) then
                  attach_q=qlabel(line)
                  exit
               endif
            enddo
         endif
         if (attach_q.eq.0) then
            write (*,*) 'Cannot locate gluon attachment line',word,gluon_links,attach
            stop 1
         endif
         qleg=find_three_line_label(qleg,qlabel)
         attach_q=find_three_line_label(attach_q,qlabel)
         if (qleg.eq.0 .or. attach_q.eq.0) then
            write (*,*) 'Cannot locate gluon Fierz endpoints',word,gluon_links
            stop 1
         endif
         temp=aq_for_q(qleg)
         aq_for_q(qleg)=aq_for_q(attach_q)
         aq_for_q(attach_q)=temp
      enddo

      line=0
      do pos=1,size(transformed)
         if (is_quark_from_order(transformed(pos),iproc)) then
            line=find_three_line_label(transformed(pos),qlabel)
         elseif (is_antiquark_from_order(transformed(pos),iproc)) then
            if (line.eq.0) then
               write (*,*) 'Antiquark before quark in reconstructed colour flow',word
               stop 1
            endif
            transformed(pos)=aq_for_q(line)
         endif
      enddo
      flow_sign=1
      do i=1,3
         if (ew_pairs(i).ne.0) flow_sign=-flow_sign
      enddo

    end subroutine construct_three_line_flow

    subroutine decode_three_line_pair(code,iproc,qout,aout)
      implicit none
      integer,intent(in) :: code,iproc
      integer,intent(out) :: qout,aout
      integer :: left,right
      left=code/(n+1)
      right=mod(code,n+1)
      if (is_quark_from_order(left,iproc)) then
         qout=left
         aout=right
      elseif (is_quark_from_order(right,iproc)) then
         qout=right
         aout=left
      else
         write (*,*) 'Closure does not contain a quark and antiquark',code,left,right
         stop 1
      endif
    end subroutine decode_three_line_pair

    integer function find_three_line_label(label,labels)
      implicit none
      integer,intent(in) :: label
      integer,dimension(3),intent(in) :: labels
      do find_three_line_label=1,3
         if (labels(find_three_line_label).eq.label) return
      enddo
      find_three_line_label=0
    end function find_three_line_label

    subroutine add_three_line_baseline_pair(code,iproc,qlabel,aq_for_q,word)
      implicit none
      integer,intent(in) :: code,iproc
      integer,dimension(3),intent(in) :: qlabel
      integer,dimension(3),intent(inout) :: aq_for_q
      integer,dimension(:),intent(in) :: word
      integer :: qout,aout,qline
      call decode_three_line_pair(code,iproc,qout,aout)
      qline=find_three_line_label(qout,qlabel)
      if (qline.eq.0) then
         write (*,*) 'Cannot locate baseline quark endpoint',word,code
         stop 1
      endif
      if (aq_for_q(qline).ne.0 .and. aq_for_q(qline).ne.aout) then
         write (*,*) 'Conflicting fermion-pair history',word,code,aq_for_q(qline)
         stop 1
      endif
      aq_for_q(qline)=aout
    end subroutine add_three_line_baseline_pair

    subroutine transform_quark_line_flow(word,pairs,iproc,transformed,flow_sign)
      implicit none
      integer,dimension(:),intent(in) :: word
      integer,dimension(3),intent(in) :: pairs
      integer,intent(in) :: iproc
      integer,dimension(size(word)),intent(out) :: transformed
      integer,intent(out) :: flow_sign
      integer :: ipair,first,second,qleg,aleg,line,pos,nlines,other,temp
      integer,dimension(3) :: qlabel,aqpos

      transformed=word
      flow_sign=1
      do ipair=1,3
         if (pairs(ipair).eq.0) cycle
         ! Closing a colour-singlet fermion pair changes the canonical
         ! Grassmann ordering once, irrespective of whether its antiquark
         ! endpoint was already on the target string.
         flow_sign=-flow_sign
         first=abs(pairs(ipair))/(n+1)
         second=mod(abs(pairs(ipair)),n+1)
         if (is_quark_from_order(first,iproc)) then
            qleg=first
            aleg=second
         elseif (is_quark_from_order(second,iproc)) then
            qleg=second
            aleg=first
         else
            write (*,*) 'EW closure does not contain a quark and antiquark',first,second
            stop 1
         endif
         qlabel=0
         aqpos=0
         nlines=0
         do pos=1,size(transformed)
            if (is_quark_from_order(transformed(pos),iproc)) then
               if (nlines.gt.0) aqpos(nlines)=pos-1
               nlines=nlines+1
               qlabel(nlines)=transformed(pos)
            endif
         enddo
         aqpos(nlines)=size(transformed)
         line=0
         other=0
         do pos=1,nlines
            if (qlabel(pos).eq.qleg) line=pos
            if (transformed(aqpos(pos)).eq.aleg) other=pos
         enddo
         if (line.eq.0 .or. other.eq.0) then
            write (*,*) 'Cannot locate EW quark-pair endpoints in colour flow',&
                 transformed,qleg,aleg
            stop 1
         endif
         if (line.ne.other) then
            temp=transformed(aqpos(line))
            transformed(aqpos(line))=transformed(aqpos(other))
            transformed(aqpos(other))=temp
         endif
      enddo
    end subroutine transform_quark_line_flow

    logical function same_three_line_flow(left,right,iproc_in)
      implicit none
      integer,dimension(:),intent(in) :: left,right
      integer,intent(in),optional :: iproc_in
      integer :: line,other,candidate,flow_iproc
      integer,dimension(3) :: left_len,right_len
      integer,dimension(n,3) :: left_lines,right_lines

      flow_iproc=1
      if (present(iproc_in)) flow_iproc=iproc_in
      call split_three_line_flow(left,flow_iproc,left_len,left_lines)
      call split_three_line_flow(right,flow_iproc,right_len,right_lines)
      same_three_line_flow=.true.
      do line=1,3
         other=0
         do candidate=1,3
            if (right_lines(1,candidate).eq.left_lines(1,line)) then
               other=candidate
               exit
            endif
         enddo
         if (other.eq.0) then
            same_three_line_flow=.false.
            return
         endif
         if (left_len(line).ne.right_len(other)) then
            same_three_line_flow=.false.
            return
         endif
         if (any(left_lines(1:left_len(line),line).ne.&
              right_lines(1:right_len(other),other))) then
            same_three_line_flow=.false.
            return
         endif
      enddo
    end function same_three_line_flow

    subroutine split_three_line_flow(word,iproc,line_len,lines)
      implicit none
      integer,dimension(:),intent(in) :: word
      integer,intent(in) :: iproc
      integer,dimension(3),intent(out) :: line_len
      integer,dimension(n,3),intent(out) :: lines
      integer :: pos,label,line,nord

      line_len=0
      lines=0
      line=0
      nord=size(word)
      do pos=1,nord
         label=word(pos)
         if (is_quark_from_order(label,iproc)) then
            line=line+1
            if (line.gt.3) then
               write (*,*) 'Too many strings in three-line colour flow',word
               stop 1
            endif
         endif
         if (line.eq.0) then
            write (*,*) 'Three-line colour flow does not begin with a quark',word
            stop 1
         endif
         line_len(line)=line_len(line)+1
         lines(line_len(line),line)=label
      enddo
      if (line.ne.3) then
         write (*,*) 'Incomplete three-line colour flow',word
         stop 1
      endif
    end subroutine split_three_line_flow

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
!!$      if (this%nprocs.gt.128) then
!!$         write (*,*) 'ERROR: too many processes. Not compatible with the "integer(kind=16) :: iproc" variable'
!!$         stop 1
!!$      endif
    end subroutine simple_consistency_checks

    subroutine define_canonical_color_order()
      ! canonical order: (q,glu,glu,glu,qbar,q,singlet,singlet,qbar)
      use math_functions
      implicit none
      integer :: i,nq,naq,nglu,nsing,iq,iaq,iglu,ising
      integer(kind=8) :: color_orders_64
      do iproc=1,this%nprocs
         nq=0; naq=0 ; nglu=0 ; nsing=0
         do i=1,n
            if (pm%is_colour_flow_vector(this%processes(i,iproc))) nglu=nglu+1
            if (is_quark_from_order(i,iproc)) nq=nq+1
            if (is_antiquark_from_order(i,iproc)) naq=naq+1
            if (pm%is_singlet(this%processes(i,iproc))) nsing=nsing+1
         enddo
         if (nq.ne.naq) then
            write (*,*) 'not the same number of quarks and anti-quarks',nq,naq
            stop 1
         endif
         if (nq.gt.3) then
            write (*,*) 'more than three quarks',nq
            stop 1
         endif
         if (nq+naq+nsing+nglu.ne.n) then
            write (*,*) 'particle types do not add up',nq,naq,nsing,nglu,':',n
            stop 1
         endif
         if (nq.eq.3) then
            ! The three-line order was canonicalised while checking the input.
            ! Rebuilding it with the two-line layout below would drop a line.
            cycle
         endif
         iq=0; iaq=0 ; iglu=0 ; ising=0
         order(1:n,1,iproc)= 0
         do i=1,n
            if (pm%is_colour_flow_vector(this%processes(i,iproc))) then
               iglu=iglu+1
               if (nq.ge.1) then
                  order(iglu+1,1,iproc)=i
               else
                  order(iglu,1,iproc)=i
               endif
            elseif(is_quark_from_order(i,iproc)) then
               iq=iq+1
               if (iq.eq.1) then
                  order(1,1,iproc)=i
               else
                  order(n-1-nsing,1,iproc)=i
               endif
            elseif (is_antiquark_from_order(i,iproc)) then
               iaq=iaq+1
               if (iaq.eq.1) then
                  order(n,1,iproc)=i
               else
                  order(n-2-nsing,1,iproc)=i
               endif
            elseif (pm%is_singlet(this%processes(i,iproc))) then
               ising=ising+1
               if(nq.ne.0) then
                  order(n-ising,1,iproc)=i
               else
                  order(nglu+ising,1,iproc)=i
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
      elseif (nq.eq.3) then
         ! Distribute the gluons over three ordered strings and connect the
         ! three quarks to the three antiquarks in all 3! ways.
         color_orders_64=3_8*int(nglu+1,kind=8)*int(nglu+2,kind=8)
         do i=2,nglu
            if (color_orders_64.gt.max_three_line_color_orders/int(i,kind=8)) then
               color_orders_64=max_three_line_color_orders+1_8
               exit
            endif
            color_orders_64=color_orders_64*int(i,kind=8)
         enddo
         if (color_orders_64.gt.max_three_line_color_orders) then
            write (*,*) 'Three-line colour basis exceeds supported size',&
                 color_orders_64,max_three_line_color_orders
            stop 1
         endif
         this%nColOrd=int(color_orders_64)
      else
         write (*,*) 'Number of colour orders unknown',nq
         stop 1
      endif
    end subroutine define_canonical_color_order

    integer function number_of_quark_lines(process)
      implicit none
      integer,dimension(n),intent(in) :: process
      integer :: i
      number_of_quark_lines=0
      do i=1,n
         if (pm%is_quark(process(i)) .or. pm%is_antiquark(process(i))) then
            number_of_quark_lines=number_of_quark_lines+1
         endif
      enddo
      number_of_quark_lines=number_of_quark_lines/2
    end function number_of_quark_lines
    
    subroutine check_input_consistency(part)
      implicit none
      integer :: i,j,iproc
      integer,dimension(n,n_processes),intent(in) :: part
      if (this%imode.eq.2) then
         if (n_processes.ne.1) then
            write (*,*) 'There should only be one process when doing imode=2'
            stop 1
         endif
      endif
      this%nprocs=n_processes
      allocate(this%processes(n,this%nprocs))
      allocate(this%n_qqbar(1:this%nprocs))
      this%processes(1:n,1:this%nprocs)=part(1:n,1:this%nprocs)
      do iproc=2,this%nprocs
         do i=1,n
            if (pm%is_lepton_any(this%processes(i,iproc)) .neqv. &
                 pm%is_lepton_any(this%processes(i,1))) then
               write (*,*) 'Lepton positions must agree within a process group',i,iproc
               stop 1
            endif
         enddo
      enddo
      do iproc=1,this%nprocs
         this%n_qqbar(iproc)=number_of_quark_lines(this%processes(1,iproc))
      enddo
      if (all(this%n_qqbar.le.2)) then
         if (this%imode.ne.2 .and. any(this%n_qqbar.eq.2) .and. &
              .not.compact_max_as_build) then
            allocate(order(1:n,2,this%nprocs))
         else
            allocate(order(1:n,1,this%nprocs))
         endif
      elseif (all(this%n_qqbar.le.3)) then
         if (this%imode.ne.2 .and. .not.compact_max_as_build) then
            ! For a fixed three-line coefficient, colour-singlet closures can
            ! map any endpoint pairing into the requested physical flow.  The
            ! twelve entries are the six endpoint permutations, with both
            ! orders of the two non-terminal strings.  They remain internal
            ! recursive representations, not additional physical amplitudes.
            allocate(order(1:n,12,this%nprocs))
         else
            allocate(order(1:n,2,this%nprocs))
         endif
      else
         write (*,*) 'Cannot allocated all the needed orders'
         stop 1
      endif
      order(1:n,1,1:this%nprocs)=o(1:n,1:this%nprocs)
      do j=2,size(order,2)
         order(1:n,j,1:this%nprocs)=order(1:n,1,1:this%nprocs)
      enddo
      allocate(this%n_sing(1:this%nprocs))
      allocate(this%same_flav(1:this%nprocs))
      do iproc=1,this%nprocs
         if (this%n_qqbar(iproc).eq.3) then
            if (this%imode.eq.2) then
               call canonicalize_three_line_order(iproc)
               call fill_alternative_quark_order(iproc)
            elseif (compact_max_as_build) then
               call fill_alternative_quark_order(iproc)
            else
               call fill_fixed_three_line_orders(iproc)
            endif
         elseif (this%n_qqbar(iproc).eq.2 .and. this%imode.ne.2 .and. &
              .not.compact_max_as_build) then
            call fill_alternative_two_line_order(iproc)
         endif
         this%same_flav(iproc)=.false. ! This will be updated once the numerical check using 'find_same_flavour' is done
         this%n_sing(iproc)=0
         do i=1,n
            if (pm%is_singlet(this%processes(i,iproc))) this%n_sing(iproc)=this%n_sing(iproc)+1
         enddo
         if (iproc.gt.1) then
!!$            if (this%n_qqbar(iproc-1).gt.this%n_qqbar(iproc)) then
!!$               write (*,*) 'ERROR: processes not correctly ordered in the list.'
!!$               write (*,*) 'Need to be in increasing number of quark lines.',iproc
!!$               do i=1,iproc
!!$                  write (*,*) this%processes(:,i),':',this%n_qqbar(i)
!!$               enddo
!!$               stop 1
!!$            endif
!!$            if (this%same_flav(iproc-1) .and. (.not.this%same_flav(iproc))) then
!!$               write (*,*) 'ERROR: processes not correctly ordered in the list.'
!!$               write (*,*) 'Need first different-flavour and then same-flavour processes.'
!!$               stop 1
!!$            endif
         endif
         if (this%n_qqbar(iproc).gt.3) then
            write (*,*) 'ERROR: code only working for 0, 1, 2 or 3 qqbar pairs',this%n_qqbar(iproc),iproc
            write (*,*) this%processes(1:n,iproc)
            stop 1
         endif
         if (imode.eq.1.or.imode.eq.3) then
            if (any(order(:,1,iproc).gt.n) .or. any(order(:,1,iproc).lt.1)) then
               write (*,*) 'ERROR: inconsistent colour order. An element is too large or too small',order(1:n,1,iproc),iproc
               stop 1
            endif
            do i=1,n-1
               do j=i+1,n
                  if (order(i,1,iproc).eq.order(j,1,iproc)) then
                     write (*,*) 'ERROR: inconsistent colour order. An element appears twice',order(1:n,1,iproc),iproc
                     stop 1
                  endif
               enddo
            enddo
            if (this%n_qqbar(iproc).gt.0) then
               if (.not.(is_quark_from_order(order(1,1,iproc),iproc))) then
                  write (*,*) 'ERROR: first particle in order is not a final state quark (or initial state anti-quark)'
                  write (*,*) iproc
                  write (*,*) order(1:n,1,iproc)
                  write (*,*) this%processes(1:n,iproc)
                  stop 1
               endif
               if (.not.(is_antiquark_from_order(order(n,1,iproc),iproc))) then
                  write (*,*) 'ERROR: final particle in order is not a final state anti-quark (or initial state quark)'
                  write (*,*) iproc
                  write (*,*) order(1:n,1,iproc)
                  write (*,*) this%processes(1:n,iproc)
                  stop 1
               endif
            endif
            if (this%n_qqbar(iproc).ge.2) then
               do i=2,n-1
                  if (is_antiquark_from_order(order(i,1,iproc),iproc)) then
                     ! next should be a quark
                     if (.not.(is_quark_from_order(order(i+1,1,iproc),iproc))) then
                        write (*,*) 'ERROR: in the colour order, after an initial state quark should come a final state quark'
                        write (*,*) iproc
                        write (*,*) order(1:n,1,iproc)
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
      enddo
    end subroutine check_input_consistency

    subroutine canonicalize_three_line_order(iproc)
      implicit none
      integer,intent(in) :: iproc
      integer :: i,pos,nq,naq,nglu,nsing
      integer,dimension(3) :: q,aq
      integer,dimension(n) :: gluons,singlets

      q=0
      aq=0
      gluons=0
      singlets=0
      nq=0
      naq=0
      nglu=0
      nsing=0
      do i=1,n
         if (is_quark_from_order(i,iproc)) then
            nq=nq+1
            if (nq.le.3) q(nq)=i
         elseif (is_antiquark_from_order(i,iproc)) then
            naq=naq+1
            if (naq.le.3) aq(naq)=i
         elseif (pm%is_colour_flow_vector(this%processes(i,iproc))) then
            nglu=nglu+1
            gluons(nglu)=i
         elseif (pm%is_singlet(this%processes(i,iproc))) then
            nsing=nsing+1
            singlets(nsing)=i
         else
            write (*,*) 'Unknown particle in three-line colour order',&
                 this%processes(i,iproc)
            stop 1
         endif
      enddo
      if (nq.ne.3 .or. naq.ne.3) then
         write (*,*) 'Three-line colour order has inconsistent endpoints',nq,naq
         stop 1
      endif

      pos=1
      order(pos,1,iproc)=q(1)
      do i=1,nglu
         pos=pos+1
         order(pos,1,iproc)=gluons(i)
      enddo
      pos=pos+1
      order(pos,1,iproc)=aq(1)
      pos=pos+1
      order(pos,1,iproc)=q(2)
      pos=pos+1
      order(pos,1,iproc)=aq(2)
      pos=pos+1
      order(pos,1,iproc)=q(3)
      do i=1,nsing
         pos=pos+1
         order(pos,1,iproc)=singlets(i)
      enddo
      pos=pos+1
      order(pos,1,iproc)=aq(3)
      if (pos.ne.n) then
         write (*,*) 'Three-line canonical colour order has wrong size',pos,n
         stop 1
      endif
    end subroutine canonicalize_three_line_order

    subroutine fill_alternative_quark_order(iproc)
      implicit none
      integer,intent(in) :: iproc
      integer :: i,ncolored
      integer,dimension(3) :: q,aq
      integer,dimension(n) :: colored_positions,colored_order,alternative
      q=0
      aq=0
      ncolored=0
      colored_positions=0
      colored_order=0
      alternative=0
      do i=1,n
         if (pm%is_singlet(this%processes(order(i,1,iproc),iproc))) cycle
         ncolored=ncolored+1
         colored_positions(ncolored)=i
         colored_order(ncolored)=order(i,1,iproc)
         if (is_quark_from_order(order(i,1,iproc),iproc)) then
            if (q(1).eq.0) then
               q(1)=ncolored
            elseif (q(2).eq.0) then
               q(2)=ncolored
            elseif (q(3).eq.0) then
               q(3)=ncolored
            endif
         endif
         if (is_antiquark_from_order(order(i,1,iproc),iproc)) then
            if (aq(1).eq.0) then
               aq(1)=ncolored
            elseif (aq(2).eq.0) then
               aq(2)=ncolored
            elseif (aq(3).eq.0) then
               aq(3)=ncolored
            endif
         endif
      enddo
      if (aq(1).ne.q(2)-1) then
         write (*,*) 'Second quark should come right after first anti-quark in colour order'
         write (*,*) q,aq,ncolored,colored_order(1:ncolored)
         stop 1
      endif
      if (aq(2).ne.q(3)-1) then
         write (*,*) 'Third quark should come right after second anti-quark in colour order'
         write (*,*) q,aq,ncolored,colored_order(1:ncolored)
         stop 1
      endif
      if (aq(3).ne.ncolored) then
         write (*,*) 'there are more coloured particles after final anti-quark'
         stop 1
      endif
      alternative(1:ncolored)=[colored_order(q(2):aq(2)),&
           colored_order(q(1):aq(1)),colored_order(q(3):aq(3))]
      order(:,2,iproc)=order(:,1,iproc)
      do i=1,ncolored
         order(colored_positions(i),2,iproc)=alternative(i)
      enddo
    end subroutine fill_alternative_quark_order

    subroutine fill_fixed_three_line_orders(iproc)
      ! Enumerate every raw three-string representation which can transform
      ! into one prescribed physical flow while keeping the closing external
      ! antiquark fixed.  For each endpoint permutation, the two strings not
      ! ending on the closing leg may occur in either recursive order.
      implicit none
      integer,intent(in) :: iproc
      integer,parameter,dimension(3,6) :: endpoint_permutations=reshape(&
           [1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[3,6])
      integer :: i,ncolored,nq,naq,iperm,iorient,iorder,line,last_line,&
           nremaining,pos,prefix_pos
      integer,dimension(3) :: q,aq,remaining,line_order
      integer,dimension(n) :: colored_positions,colored_order,alternative

      ncolored=0
      nq=0
      naq=0
      q=0
      aq=0
      colored_positions=0
      colored_order=0
      do i=1,n
         if (pm%is_singlet(this%processes(order(i,1,iproc),iproc))) cycle
         ncolored=ncolored+1
         colored_positions(ncolored)=i
         colored_order(ncolored)=order(i,1,iproc)
         if (is_quark_from_order(colored_order(ncolored),iproc)) then
            nq=nq+1
            if (nq.le.3) q(nq)=ncolored
         elseif (is_antiquark_from_order(colored_order(ncolored),iproc)) then
            naq=naq+1
            if (naq.le.3) aq(naq)=ncolored
         endif
      enddo
      if (nq.ne.3 .or. naq.ne.3 .or. q(1).ne.1 .or. &
           aq(1).ne.q(2)-1 .or. aq(2).ne.q(3)-1 .or. &
           aq(3).ne.ncolored) then
         write (*,*) 'Malformed three-line fixed colour order',order(:,1,iproc)
         stop 1
      endif

      do iperm=1,6
         last_line=0
         remaining=0
         nremaining=0
         do line=1,3
            if (endpoint_permutations(line,iperm).eq.3) then
               last_line=line
            else
               nremaining=nremaining+1
               remaining(nremaining)=line
            endif
         enddo
         if (last_line.eq.0 .or. nremaining.ne.2) then
            write (*,*) 'Invalid internal three-line endpoint permutation',iperm
            stop 1
         endif
         do iorient=1,2
            if (iorient.eq.1) then
               line_order=[remaining(1),remaining(2),last_line]
            else
               line_order=[remaining(2),remaining(1),last_line]
            endif
            alternative=0
            pos=0
            do i=1,3
               line=line_order(i)
               do prefix_pos=q(line),aq(line)-1
                  pos=pos+1
                  alternative(pos)=colored_order(prefix_pos)
               enddo
               pos=pos+1
               alternative(pos)=colored_order(&
                    aq(endpoint_permutations(line,iperm)))
            enddo
            if (pos.ne.ncolored .or. alternative(ncolored).ne.&
                 colored_order(aq(3))) then
               write (*,*) 'Cannot construct fixed three-line raw order',&
                    iperm,iorient,alternative(1:pos)
               stop 1
            endif
            iorder=2*(iperm-1)+iorient
            order(:,iorder,iproc)=order(:,1,iproc)
            do i=1,ncolored
               order(colored_positions(i),iorder,iproc)=alternative(i)
            enddo
         enddo
      enddo
    end subroutine fill_fixed_three_line_orders

    subroutine fill_alternative_two_line_order(iproc)
      ! Add the complementary two-string representation without changing the
      ! prescribed terminal leg.  This lets a fixed-order object import the
      ! colour-singlet coefficient which is generated in the opposite raw
      ! planar ordering.
      implicit none
      integer,intent(in) :: iproc
      integer :: i,ncolored,nq,naq,first_q,second_q,first_aq,second_aq,&
           first_len,second_len,pos
      integer,dimension(n) :: colored_positions,colored_order,alternative

      ncolored=0
      nq=0
      naq=0
      first_q=0
      second_q=0
      first_aq=0
      second_aq=0
      colored_positions=0
      colored_order=0
      alternative=0
      do i=1,n
         if (pm%is_singlet(this%processes(order(i,1,iproc),iproc))) cycle
         ncolored=ncolored+1
         colored_positions(ncolored)=i
         colored_order(ncolored)=order(i,1,iproc)
         if (is_quark_from_order(colored_order(ncolored),iproc)) then
            nq=nq+1
            if (nq.eq.1) first_q=ncolored
            if (nq.eq.2) second_q=ncolored
         elseif (is_antiquark_from_order(colored_order(ncolored),iproc)) then
            naq=naq+1
            if (naq.eq.1) first_aq=ncolored
            if (naq.eq.2) second_aq=ncolored
         endif
      enddo
      if (nq.ne.2 .or. naq.ne.2 .or. first_q.ne.1 .or. &
           first_aq.ne.second_q-1 .or. second_aq.ne.ncolored) then
         write (*,*) 'Malformed two-line fixed colour order',order(:,1,iproc)
         stop 1
      endif
      first_len=first_aq-first_q
      second_len=second_aq-second_q
      pos=0
      ! Prefix of the second string, followed by the first antiquark.
      alternative(1:second_len)=colored_order(second_q:second_aq-1)
      pos=second_len
      pos=pos+1
      alternative(pos)=colored_order(first_aq)
      ! Prefix of the first string, followed by the original terminal
      ! antiquark.  Hence order(n,2,iproc)==order(n,1,iproc).
      alternative(pos+1:pos+first_len)=colored_order(first_q:first_aq-1)
      pos=pos+first_len
      pos=pos+1
      alternative(pos)=colored_order(second_aq)
      if (pos.ne.ncolored) then
         write (*,*) 'Cannot construct complementary two-line order',&
              colored_order(1:ncolored),alternative(1:pos)
         stop 1
      endif
      order(:,2,iproc)=order(:,1,iproc)
      do i=1,ncolored
         order(colored_positions(i),2,iproc)=alternative(i)
      enddo
    end subroutine fill_alternative_two_line_order
    
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
         if (allocated(current_list_local(ic)%vertices)) &
              allocate(tmp(ic)%vertices(size(current_list_local(ic)%vertices)))
         if (allocated(current_list_local(ic)%vertex_sign)) &
              allocate(tmp(ic)%vertex_sign(size(current_list_local(ic)%vertex_sign)))
         tmp(ic)=current_list_local(ic)
      enddo
      do ic=1,max_cur
         call finalize_current(current_list_local(ic))
      enddo
      deallocate(current_list_local)
      allocate(current_list_local(new_max_cur))
      do ic=1,max_cur
         if (allocated(tmp(ic)%vertices)) &
              allocate(current_list_local(ic)%vertices(size(tmp(ic)%vertices)))
         if (allocated(tmp(ic)%vertex_sign)) &
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
      real(kind=4) :: sgn
      integer :: ichir
      if (.not.valid_current_combination())  then
         return
      endif
      do i=1,pm%nint
         if ( current_list_local(ic1)%type.eq.pm%vertex_list(i)%particles(1) .and. &
              current_list_local(ic2)%type.eq.pm%vertex_list(i)%particles(2) ) then
            ! Colour-singlet triple-vector vertices are unordered physical
            ! vertices.  With several quark lines the two planar current
            ! representations otherwise insert the same VVV diagram twice
            ! (once for each ordering of its two child currents).
            if (all(this%n_qqbar(1:this%nprocs).eq.3) .and. &
                 pm%vertex_list(i)%type.eq.12 .and. &
                 current_list_local(ic1)%bin.gt.current_list_local(ic2)%bin) cycle
            ichir=vertex_result_chirality(pm%vertex_list(i)%type, &
                 pm%vertex_list(i)%particles(3),pm%vertex_list(i)%coupl)
            if (ichir.eq.-99) cycle
            sgn=1d0
              call add_vertex(pm%vertex_list(i)%type, &
                            pm%vertex_list(i)%particles(3), &
                            sgn*pm%vertex_list(i)%coupl,ichir,&
                            pm%vertex_list(i)%n_gs,pm%vertex_list(i)%n_ew)
         endif
      enddo
    end subroutine add_if_allowed_threevertex

    integer function ext_from_cur(ic)
      implicit none
      integer :: ic,ncur,nc,irproc,ispin
      ncur=0
      do nc=n,1,-1
         do iproc=1,this%nprocs
            if (this%same_flav(iproc)) cycle
            do ispin=1, spin(0,order(nc,1,iproc))
               ncur=ncur+1
               if (ncur.eq.ic) then
                       goto 13
               endif
            enddo
         enddo
      enddo
13   ext_from_cur=order(nc,1,iproc)
    end function ext_from_cur
    
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
      integer :: i,j,nc1,nc2,c
      logical :: gluon_current,colour_singlet1,colour_singlet2,found_quark,found_antiquark,valid
      integer,dimension(isize) :: ip,et
      type(bitset) :: iproc_combined
      valid_current_combination=.false.
      ! check that all particles are different in the two currents:
      if (popcnt(ieor(current_list_local(ic1)%bin,current_list_local(ic2)%bin)).ne.isize) return
      ! final particle should never be part of any combined currents: it will
      ! be used to close the amplitude instead
      if (n1.eq.1) then
         if (all(current_list_local(ic1)%order(n1).eq.order(n,1,1:this%nprocs))) then
            return
          endif
      endif
      if (n2.eq.1) then
         if (all(current_list_local(ic2)%order(n2).eq.order(n,1,1:this%nprocs))) then
            return
         endif
      endif
      ! check that both currents can contribute to the same process
      iproc_combined=current_list_local(ic1)%iproc.and.current_list_local(ic2)%iproc
      if (iproc_combined%count_bits().eq.0) return
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
            if (pm%is_singlet(current_list_local(ic1)%ext_type(i))) exit
         enddo
         nc1=i-1
         do i=1,n2
            if (pm%is_singlet(current_list_local(ic2)%ext_type(i))) exit
         enddo
         nc2=i-1
         ip(1:nc1+nc2)=[current_list_local(ic1)%order(1:nc1),current_list_local(ic2)%order(1:nc2)]
         ! do they actual checking:
         valid=.false.
         do_iproc: do iproc=1,this%nprocs
            ! check that both currents contribute to the iproc process:
            if (.not. iproc_combined%test_bit(iproc)) cycle
            ! check that the final particle is not part of the combined
            ! current (it will be used to close the amplitude instead):
            if (btest(current_list_local(ic1)%bin+current_list_local(ic2)%bin,order(n,1,iproc)-1)) cycle
            do_c: do c=1,size(order,2) ! Fixed multi-line coefficients can need complementary raw orders.
               ! Check if they are compatible with the colour order of the iproc:
               do_j: do j=1,n
                  if (order(j,c,iproc).eq.ip(1)) then
                     do i=2,nc1+nc2
                        if (j-1+i.gt.n) exit do_j
                        if (order(j-1+i,c,iproc).ne.ip(i)) exit do_j
                     enddo
                     valid=.true. ! it's compatible with the input colour order of iproc
                     exit do_iproc
                  endif
               enddo do_j
            enddo do_c
         enddo do_iproc
         if (.not.valid) return ! not compatible with any of the iprocs
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
            if (pm%is_quark(et(i))) then
               ! found a quark.
               if (found_quark) then
                  ! no anti-quark between two quarks
                  return
               endif
               found_antiquark=.false.
               found_quark=.true.
            elseif (pm%is_antiquark(et(i))) then
               ! found an anti-quark
               if (found_antiquark) then
                  ! no quark between two anti-quarks
                  return
               endif
               ! next one must be a quark:
               j=i
               do while (j.lt.isize)
                  if (pm%is_singlet(et(j+1))) then
                     j=j+1
                  elseif (.not.(pm%is_quark(et(j+1)))) then
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
         if (isize.eq.n-1 .and. .not.pm%is_quark(et(1))) return
      endif
      ! Got all the way to the end. This must be a valid current combination
      valid_current_combination=.true.
    end function valid_current_combination
    
    subroutine add_vertex(itype,ctype,coupl,ichir,n_gs,n_ew)
      implicit none
      integer :: itype,ctype,ichir,ic,n_gs,n_ew
      real(kind=8),dimension(2) :: coupl
      if (compact_max_as_build) then
         if (current_list_local(ic1)%n_ew+current_list_local(ic2)%n_ew+n_ew.gt.&
              compact_ew_order) return
         ! With at least two quark lines, the compact recursion represents a
         ! fixed QCD colour backbone with electroweak particles radiated from
         ! its fermion lines.  Crossed q-qbar-to-vector closures replace a QCD
         ! link and therefore cannot enter this leading-aS sector.  Keep them
         ! for zero/one-line amplitudes, where they can be the physical
         ! annihilation current (for example Drell--Yan production).
         if (compact_qcd_backbone .and. (itype.eq.21 .or. itype.eq.22) .and. &
              current_list_local(ic1)%open_quark_leg.ne.0 .and. &
              current_list_local(ic2)%open_quark_leg.ne.0) return
      endif
      if (isize.eq.n-1) then
         do ic=this%n_cur_start(n),this%n_cur_end(n)
            if (ctype.eq.anti_current(current_list_local(ic)%type)) exit
         enddo
         if (ic.eq.this%n_cur_end(n)+1) return ! dead tree. Filter already here
      endif
      this%n_vert=this%n_vert+1
      if (this%n_vert.gt.max_vert) call increase_max_vert()
      interaction_list_local(this%n_vert)%type=itype
      interaction_list_local(this%n_vert)%chirality=ichir
      interaction_list_local(this%n_vert)%n_gs=&
           current_list_local(ic1)%n_gs+current_list_local(ic2)%n_gs+n_gs
      interaction_list_local(this%n_vert)%n_ew=&
           current_list_local(ic1)%n_ew+current_list_local(ic2)%n_ew+n_ew
      interaction_list_local(this%n_vert)%currents(1)=ic1
      interaction_list_local(this%n_vert)%currents(2)=ic2
      interaction_list_local(this%n_vert)%coupl=coupl
      allocate(interaction_list_local(this%n_vert)%singlet_mv(0:isize))
      call add_all_currents(ctype,ichir)
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

    type(current) function combine_currents(ic1,ic2,ctype,ichir,singlet_mv,invert)
      ! combine the currents corresponding to ic1 and ic2 into a new current
      ! of type 'ctype'. This also sets up 'singlet_mv' that determines how to
      ! move the colour singlets to the correct position. If the first (or
      ! second) bit of 'invert' is set to 1, the colour order of ic1 (or ic2)
      ! is reversed before the two currents are combined.
      implicit none
      integer,intent(in) :: ic1,ic2,ctype,ichir
      integer,intent(in) :: invert
      integer,dimension(0:isize),intent(out) :: singlet_mv
      integer :: i,j,n1,n2,ipos,mv12,nc1,nc2,ns1,ns2,n_mv12_1
      integer :: n_pairs,pair_code,tmp_pair
      integer,dimension(isize) :: ord
      integer,dimension(:),allocatable :: ord1,spin1,et1,ord2,spin2,et2
      allocate(combine_currents%order(1:isize))
      allocate(combine_currents%spin(1:isize))
      allocate(combine_currents%ext_type(1:isize))
      combine_currents%type=ctype
      combine_currents%chirality=ichir
      combine_currents%n_gs=interaction_list_local(this%n_vert)%n_gs
      combine_currents%n_ew=interaction_list_local(this%n_vert)%n_ew
      combine_currents%open_quark_leg=0
      combine_currents%ew_pairs=0
      combine_currents%u1_pairs=0
      combine_currents%gluon_pairs=0
      combine_currents%u1_links=0
      combine_currents%gluon_links=0
      if (pm%is_quark(ctype).or.pm%is_antiquark(ctype)) then
         if (current_list_local(ic1)%open_quark_leg.ne.0) then
            combine_currents%open_quark_leg=current_list_local(ic1)%open_quark_leg
         else
            combine_currents%open_quark_leg=current_list_local(ic2)%open_quark_leg
         endif
      endif
      n_pairs=0
      do j=1,3
         if (current_list_local(ic1)%ew_pairs(j).eq.0) cycle
         n_pairs=n_pairs+1
         combine_currents%ew_pairs(n_pairs)=current_list_local(ic1)%ew_pairs(j)
      enddo
      do j=1,3
         if (current_list_local(ic2)%ew_pairs(j).eq.0) cycle
         if (any(combine_currents%ew_pairs.eq.current_list_local(ic2)%ew_pairs(j))) cycle
         n_pairs=n_pairs+1
         if (n_pairs.gt.3) then
            write (*,*) 'Too many colour-singlet quark-line closures'
            stop 1
         endif
         combine_currents%ew_pairs(n_pairs)=current_list_local(ic2)%ew_pairs(j)
      enddo
      if ((interaction_list_local(this%n_vert)%type.eq.21 .or. &
           interaction_list_local(this%n_vert)%type.eq.22) .and. &
           current_list_local(ic1)%open_quark_leg.ne.0 .and. &
           current_list_local(ic2)%open_quark_leg.ne.0) then
         pair_code=min(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)*(n+1)+&
              max(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)
         ! Preserve which spinor-chain orientation closed the fermion pair.
         ! The colour endpoints are the same, but the two orientations carry
         ! different fermion-order conventions and must not be merged before
         ! the physical-flow signs are assigned.
         if (interaction_list_local(this%n_vert)%type.eq.22) &
              pair_code=-pair_code
         if (.not.any(combine_currents%ew_pairs.eq.pair_code)) then
            n_pairs=n_pairs+1
            if (n_pairs.gt.3) then
               write (*,*) 'Too many colour-singlet quark-line closures'
               stop 1
            endif
            combine_currents%ew_pairs(n_pairs)=pair_code
         endif
      endif
      do i=1,n_pairs-1
         do j=i+1,n_pairs
            if (combine_currents%ew_pairs(j).lt.combine_currents%ew_pairs(i)) then
               tmp_pair=combine_currents%ew_pairs(i)
               combine_currents%ew_pairs(i)=combine_currents%ew_pairs(j)
               combine_currents%ew_pairs(j)=tmp_pair
            endif
         enddo
      enddo
      n_pairs=0
      do j=1,3
         if (current_list_local(ic1)%gluon_links(j).eq.0) cycle
         n_pairs=n_pairs+1
         combine_currents%gluon_links(n_pairs)=current_list_local(ic1)%gluon_links(j)
      enddo
      do j=1,3
         if (current_list_local(ic2)%gluon_links(j).eq.0) cycle
         if (any(combine_currents%gluon_links.eq.current_list_local(ic2)%gluon_links(j))) cycle
         n_pairs=n_pairs+1
         if (n_pairs.gt.3) then
            write (*,*) 'Too many gluon quark-line links'
            stop 1
         endif
         combine_currents%gluon_links(n_pairs)=current_list_local(ic2)%gluon_links(j)
      enddo
      if (interaction_list_local(this%n_vert)%type.ge.4 .and. &
           interaction_list_local(this%n_vert)%type.le.7 .and. &
           (pm%is_quark(ctype).or.pm%is_antiquark(ctype))) then
         if (pm%is_colour_flow_vector(current_list_local(ic1)%type)) then
            do j=1,3
               if (current_list_local(ic1)%gluon_pairs(j).eq.0) cycle
               if (current_list_local(ic2)%open_quark_leg.eq.0) cycle
               pair_code=current_list_local(ic1)%gluon_pairs(j)*(n+1)+&
                    current_list_local(ic2)%open_quark_leg
               if (any(combine_currents%gluon_links.eq.pair_code)) cycle
               n_pairs=n_pairs+1
               if (n_pairs.gt.3) then
                  write (*,*) 'Too many gluon quark-line links'
                  stop 1
               endif
               combine_currents%gluon_links(n_pairs)=pair_code
            enddo
         elseif (pm%is_colour_flow_vector(current_list_local(ic2)%type)) then
            do j=1,3
               if (current_list_local(ic2)%gluon_pairs(j).eq.0) cycle
               if (current_list_local(ic1)%open_quark_leg.eq.0) cycle
               pair_code=current_list_local(ic2)%gluon_pairs(j)*(n+1)+&
                    current_list_local(ic1)%open_quark_leg
               if (any(combine_currents%gluon_links.eq.pair_code)) cycle
               n_pairs=n_pairs+1
               if (n_pairs.gt.3) then
                  write (*,*) 'Too many gluon quark-line links'
                  stop 1
               endif
               combine_currents%gluon_links(n_pairs)=pair_code
            enddo
         endif
      endif
      do i=1,n_pairs-1
         do j=i+1,n_pairs
            if (combine_currents%gluon_links(j).lt.combine_currents%gluon_links(i)) then
               tmp_pair=combine_currents%gluon_links(i)
               combine_currents%gluon_links(i)=combine_currents%gluon_links(j)
               combine_currents%gluon_links(j)=tmp_pair
            endif
         enddo
      enddo
      n_pairs=0
      do j=1,3
         if (current_list_local(ic1)%gluon_pairs(j).eq.0) cycle
         n_pairs=n_pairs+1
         combine_currents%gluon_pairs(n_pairs)=current_list_local(ic1)%gluon_pairs(j)
      enddo
      do j=1,3
         if (current_list_local(ic2)%gluon_pairs(j).eq.0) cycle
         if (any(combine_currents%gluon_pairs.eq.current_list_local(ic2)%gluon_pairs(j))) cycle
         n_pairs=n_pairs+1
         if (n_pairs.gt.3) then
            write (*,*) 'Too many gluon quark-line closures'
            stop 1
         endif
         combine_currents%gluon_pairs(n_pairs)=current_list_local(ic2)%gluon_pairs(j)
      enddo
      if (interaction_list_local(this%n_vert)%type.eq.9 .and. &
           current_list_local(ic1)%open_quark_leg.ne.0 .and. &
           current_list_local(ic2)%open_quark_leg.ne.0) then
         pair_code=min(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)*(n+1)+&
              max(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)
         if (.not.any(combine_currents%gluon_pairs.eq.pair_code)) then
            n_pairs=n_pairs+1
            if (n_pairs.gt.3) then
               write (*,*) 'Too many gluon quark-line closures'
               stop 1
            endif
            combine_currents%gluon_pairs(n_pairs)=pair_code
         endif
      endif
      do i=1,n_pairs-1
         do j=i+1,n_pairs
            if (combine_currents%gluon_pairs(j).lt.combine_currents%gluon_pairs(i)) then
               tmp_pair=combine_currents%gluon_pairs(i)
               combine_currents%gluon_pairs(i)=combine_currents%gluon_pairs(j)
               combine_currents%gluon_pairs(j)=tmp_pair
            endif
         enddo
      enddo
      n_pairs=0
      do j=1,3
         if (current_list_local(ic1)%u1_links(j).eq.0) cycle
         n_pairs=n_pairs+1
         combine_currents%u1_links(n_pairs)=current_list_local(ic1)%u1_links(j)
      enddo
      do j=1,3
         if (current_list_local(ic2)%u1_links(j).eq.0) cycle
         if (any(combine_currents%u1_links.eq.current_list_local(ic2)%u1_links(j))) cycle
         n_pairs=n_pairs+1
         if (n_pairs.gt.3) then
            write (*,*) 'Too many auxiliary-U(1) quark-line links'
            stop 1
         endif
         combine_currents%u1_links(n_pairs)=current_list_local(ic2)%u1_links(j)
      enddo
      if (interaction_list_local(this%n_vert)%type.ge.4 .and. &
           interaction_list_local(this%n_vert)%type.le.7 .and. &
           (pm%is_quark(ctype).or.pm%is_antiquark(ctype))) then
         if (current_list_local(ic1)%type.eq.99) then
            do j=1,3
               if (current_list_local(ic1)%u1_pairs(j).eq.0) cycle
               if (current_list_local(ic2)%open_quark_leg.eq.0) cycle
               pair_code=current_list_local(ic1)%u1_pairs(j)*(n+1)+&
                    current_list_local(ic2)%open_quark_leg
               if (any(combine_currents%u1_links.eq.pair_code)) cycle
               n_pairs=n_pairs+1
               if (n_pairs.gt.3) then
                  write (*,*) 'Too many auxiliary-U(1) quark-line links'
                  stop 1
               endif
               combine_currents%u1_links(n_pairs)=pair_code
            enddo
         elseif (current_list_local(ic2)%type.eq.99) then
            do j=1,3
               if (current_list_local(ic2)%u1_pairs(j).eq.0) cycle
               if (current_list_local(ic1)%open_quark_leg.eq.0) cycle
               pair_code=current_list_local(ic2)%u1_pairs(j)*(n+1)+&
                    current_list_local(ic1)%open_quark_leg
               if (any(combine_currents%u1_links.eq.pair_code)) cycle
               n_pairs=n_pairs+1
               if (n_pairs.gt.3) then
                  write (*,*) 'Too many auxiliary-U(1) quark-line links'
                  stop 1
               endif
               combine_currents%u1_links(n_pairs)=pair_code
            enddo
         endif
      endif
      do i=1,n_pairs-1
         do j=i+1,n_pairs
            if (combine_currents%u1_links(j).lt.combine_currents%u1_links(i)) then
               tmp_pair=combine_currents%u1_links(i)
               combine_currents%u1_links(i)=combine_currents%u1_links(j)
               combine_currents%u1_links(j)=tmp_pair
            endif
         enddo
      enddo
      n_pairs=0
      do j=1,3
         if (current_list_local(ic1)%u1_pairs(j).eq.0) cycle
         n_pairs=n_pairs+1
         combine_currents%u1_pairs(n_pairs)=current_list_local(ic1)%u1_pairs(j)
      enddo
      do j=1,3
         if (current_list_local(ic2)%u1_pairs(j).eq.0) cycle
         if (any(combine_currents%u1_pairs.eq.current_list_local(ic2)%u1_pairs(j))) cycle
         n_pairs=n_pairs+1
         if (n_pairs.gt.3) then
            write (*,*) 'Too many auxiliary-U(1) quark-line closures'
            stop 1
         endif
         combine_currents%u1_pairs(n_pairs)=current_list_local(ic2)%u1_pairs(j)
      enddo
      if (interaction_list_local(this%n_vert)%type.eq.8 .and. &
           current_list_local(ic1)%open_quark_leg.ne.0 .and. &
           current_list_local(ic2)%open_quark_leg.ne.0) then
         pair_code=min(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)*(n+1)+&
              max(current_list_local(ic1)%open_quark_leg,&
              current_list_local(ic2)%open_quark_leg)
         if (.not.any(combine_currents%u1_pairs.eq.pair_code)) then
            n_pairs=n_pairs+1
            if (n_pairs.gt.3) then
               write (*,*) 'Too many auxiliary-U(1) quark-line closures'
               stop 1
            endif
            combine_currents%u1_pairs(n_pairs)=pair_code
         endif
      endif
      do i=1,n_pairs-1
         do j=i+1,n_pairs
            if (combine_currents%u1_pairs(j).lt.combine_currents%u1_pairs(i)) then
               tmp_pair=combine_currents%u1_pairs(i)
               combine_currents%u1_pairs(i)=combine_currents%u1_pairs(j)
               combine_currents%u1_pairs(j)=tmp_pair
            endif
         enddo
      enddo
      combine_currents%bin=current_list_local(ic1)%bin+current_list_local(ic2)%bin
      combine_currents%iproc=current_list_local(ic1)%iproc.and.current_list_local(ic2)%iproc
      combine_currents%ext_cur=current_list_local(ic1)%ext_cur+current_list_local(ic2)%ext_cur
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
         if (pm%is_singlet(et1(i))) exit
      enddo
      nc1=i-1
      do i=1,n2
         if (pm%is_singlet(et2(i))) exit
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
         n_mv12_1=0
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
               singlet_mv(singlet_mv(0))=ns1 - n_mv12_1
               n_mv12_1=n_mv12_1+1
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

    subroutine add_all_currents(ctype,ichir)
      ! combine currents ic1 and ic2 and add them to the list of currents to
      ! compute. If use_symmetry=.true., need to consider all possible
      ! permutations allowed under the symmetry (at most 8 if both ic1 and ic2
      ! are made up of only gluons).
      implicit none
      logical,dimension(8) :: vertex_sign
      logical :: lepton_sign
      integer :: i,ctype,ichir,nperm
      integer,dimension(0:isize) :: singlet_mv
      type(current),dimension(8) :: new_currents
      ! External spinors are Grassmann odd.  The numerical wavefunctions used
      ! below commute, so restore the sign needed to put the leptons contained
      ! in the two child currents into one fixed (external-leg) order.
      if (interaction_list_local(this%n_vert)%type.eq.22) then
         ! The antilepton-lepton rule evaluates its spinor chain in the
         ! opposite order; use that same order for the Grassmann factors.
         lepton_sign=lepton_reordering_is_odd(ic2,ic1)
      else
         lepton_sign=lepton_reordering_is_odd(ic1,ic2)
      endif
      if (.not.use_symmetry .or. this%imode.eq.1 .or. this%imode.eq.3) then
         new_currents(1)=combine_currents(ic1,ic2,ctype,ichir,singlet_mv,0)
         interaction_list_local(this%n_vert)%singlet_mv(0:isize)=singlet_mv(0:isize)
         call add_current(lepton_sign,new_currents(1))
         return
      endif
      ! Need to consider all the possible permutations
      call check_all_permutations(ctype,ichir,nperm,new_currents,vertex_sign)
      do i=1,nperm
         call add_current(vertex_sign(i) .neqv. lepton_sign,new_currents(i))
      enddo
    end subroutine add_all_currents

    logical function lepton_reordering_is_odd(current1,current2)
      ! Return the parity of the permutation which changes
      !
      !   (leptons in current1), (leptons in current2)
      !
      ! into ascending external-leg order.  Applying this at every recursive
      ! merge gives all diagrams the same external-fermion convention, including
      ! diagrams with different pairings of identical leptons.
      implicit none
      integer,intent(in) :: current1,current2
      integer :: i,j,ncross
      ncross=0
      do i=1,n
         if (.not.btest(current_list_local(current1)%bin,i-1)) cycle
         if (.not.pm%is_lepton_any(this%processes(i,1))) cycle
         do j=1,i-1
            if (.not.btest(current_list_local(current2)%bin,j-1)) cycle
            if (pm%is_lepton_any(this%processes(j,1))) ncross=ncross+1
         enddo
      enddo
      lepton_reordering_is_odd=mod(ncross,2).eq.1
    end function lepton_reordering_is_odd

    subroutine check_all_permutations(ctype,ichir,nperm,new_currents,vertex_sign)
      ! If a current only contains (external) gluons, we can use symmetry to
      ! relate them to eachother. This subroutine checks all permutations,
      ! and, if they give a valid current order, adds that current to the list
      ! that should be included.
      implicit none
      integer,intent(in) :: ctype,ichir
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
                  new_currents(nperm)=combine_currents(ic1,ic2,ctype,ichir,singlet_mv(0,nperm),invert)
              else
                  new_currents(nperm)=combine_currents(ic2,ic1,ctype,ichir,singlet_mv(0,nperm),invert)
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
      logical :: sgn,three_line_current
      integer :: i,ic,key
      integer(kind=8) :: val
      integer :: lep1,lep2
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         three_line_current=.false.
         do i=1,this%nprocs
            if (new_current%iproc%test_bit(i) .and. this%n_qqbar(i).eq.3) &
                 three_line_current=.true.
         enddo
         ! Check if this interaction can be added to an existing current
         do ic=1,this%n_cur
            if (new_current%type.ne.current_list_local(ic)%type) cycle
            if (new_current%chirality.ne.current_list_local(ic)%chirality) cycle
            if (new_current%n_gs.ne.current_list_local(ic)%n_gs) cycle
            if (new_current%n_ew.ne.current_list_local(ic)%n_ew) cycle
            if (.not.compact_max_as_build) then
               if (new_current%open_quark_leg.ne.current_list_local(ic)%open_quark_leg) cycle
               if (any(new_current%ew_pairs.ne.current_list_local(ic)%ew_pairs)) cycle
               if (any(new_current%u1_pairs.ne.current_list_local(ic)%u1_pairs)) cycle
               if (any(new_current%gluon_pairs.ne.current_list_local(ic)%gluon_pairs)) cycle
               if (any(new_current%u1_links.ne.current_list_local(ic)%u1_links)) cycle
               if (any(new_current%gluon_links.ne.current_list_local(ic)%gluon_links)) cycle
            endif
            if (new_current%bin.ne.current_list_local(ic)%bin) cycle
            if (new_current%ext_cur.ne.current_list_local(ic)%ext_cur) cycle
            if (three_line_current .and. .not.compact_max_as_build) then
               if (any(new_current%order.ne.current_list_local(ic)%order)) cycle
            endif
            call append_current_vertex(ic,this%n_vert,vertex_sign)
            return
         enddo
         ! Need a new current
         this%n_cur=this%n_cur+1
         if (this%n_cur.gt.max_cur) call increase_max_cur()
         current_list_local(this%n_cur)=new_current
         current_list_local(this%n_cur)%mass=pm%get_mass(new_current%type)
         current_list_local(this%n_cur)%width=pm%get_width(new_current%type)
         
         if (pm%is_colour_flow_vector(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         elseif (pm%is_auxiliary_tensor(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(2*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(2*(isize-1)))
         elseif (pm%is_massive_vector(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
         else
            allocate(current_list_local(this%n_cur)%vertices(8*(isize-1)))
            allocate(current_list_local(this%n_cur)%vertex_sign(8*(isize-1)))
         endif
         current_list_local(this%n_cur)%vertices(1)=this%n_vert
         current_list_local(this%n_cur)%vertex_sign(1)=vertex_sign
         current_list_local(this%n_cur)%n_vert=1
      elseif (this%imode.eq.2) then
         call get_value(new_current%order,new_current%type,val)
         call solve_dict(val,key)
         ic=key_to_current(key,new_current%iproc%bitset_to_integer())
         ! Coupling order is part of current identity.  The colour-order
         ! dictionary predates this axis, so search the (usually very short)
         ! equal-key family when its representative has another order.
         if (ic.ne.0) then
            if (current_list_local(ic)%n_gs.ne.new_current%n_gs .or. &
                 current_list_local(ic)%n_ew.ne.new_current%n_ew .or. &
                 current_list_local(ic)%chirality.ne.new_current%chirality .or. &
                 current_list_local(ic)%open_quark_leg.ne.new_current%open_quark_leg .or. &
                 any(current_list_local(ic)%ew_pairs.ne.new_current%ew_pairs) .or. &
                 any(current_list_local(ic)%u1_pairs.ne.new_current%u1_pairs) .or. &
                 any(current_list_local(ic)%gluon_pairs.ne.new_current%gluon_pairs) .or. &
                 any(current_list_local(ic)%u1_links.ne.new_current%u1_links) .or. &
                 any(current_list_local(ic)%gluon_links.ne.new_current%gluon_links)) then
               ic=0
               do i=1,this%n_cur
                  if (current_list_local(i)%type.ne.new_current%type) cycle
                  if (current_list_local(i)%bin.ne.new_current%bin) cycle
                  if (current_list_local(i)%n_gs.ne.new_current%n_gs) cycle
                 if (current_list_local(i)%n_ew.ne.new_current%n_ew) cycle
                 if (current_list_local(i)%chirality.ne.new_current%chirality) cycle
                  if (current_list_local(i)%open_quark_leg.ne.new_current%open_quark_leg) cycle
                  if (any(current_list_local(i)%ew_pairs.ne.new_current%ew_pairs)) cycle
                  if (any(current_list_local(i)%u1_pairs.ne.new_current%u1_pairs)) cycle
                  if (any(current_list_local(i)%gluon_pairs.ne.new_current%gluon_pairs)) cycle
                  if (any(current_list_local(i)%u1_links.ne.new_current%u1_links)) cycle
                  if (any(current_list_local(i)%gluon_links.ne.new_current%gluon_links)) cycle
                  if (any(current_list_local(i)%order(1:isize).ne.new_current%order(1:isize))) cycle
                  if (current_list_local(i)%iproc%bitset_to_integer().ne.&
                       new_current%iproc%bitset_to_integer()) cycle
                  ic=i
                  exit
               enddo
            endif
         endif
         if (ic.eq.0) then
            ! initialise new current
            this%n_cur=this%n_cur+1
            if (this%n_cur.gt.max_cur) call increase_max_cur()
            if (key_to_current(key,new_current%iproc%bitset_to_integer()).eq.0) &
                 key_to_current(key,new_current%iproc%bitset_to_integer())=this%n_cur
            ic=this%n_cur
            current_list_local(ic)=new_current
            current_list_local(ic)%mass=pm%get_mass(new_current%type)
            current_list_local(ic)%width=pm%get_width(new_current%type)

            if (any(current_list_local(ic)%spin(1:isize).ne.-9)) then
               write (*,*) 'trying to combine currents with different spin: not possible',&
                    current_list_local(ic)%spin(1:isize)
               stop 1
            endif
            if (pm%is_colour_flow_vector(new_current%type)) then
               allocate(current_list_local(ic)%vertices(5*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(5*(isize-1)))
            elseif (pm%is_auxiliary_tensor(new_current%type)) then
               allocate(current_list_local(ic)%vertices(isize-1))
               allocate(current_list_local(ic)%vertex_sign(isize-1))
            elseif (pm%is_massive_vector(new_current%type)) then
               allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)))
               allocate(current_list_local(this%n_cur)%vertex_sign(5*(isize-1)))
            else
               allocate(current_list_local(ic)%vertices(5*(isize-1)))
               allocate(current_list_local(ic)%vertex_sign(5*(isize-1)))
            endif
            current_list_local(ic)%n_vert=0
         endif
         ! add the vertex to the current
         call append_current_vertex(ic,this%n_vert,vertex_sign)
      endif
    end subroutine add_current

    subroutine append_current_vertex(ic,iv,vertex_sign)
      implicit none
      integer,intent(in) :: ic,iv
      logical,intent(in) :: vertex_sign
      integer :: old_capacity,new_capacity
      integer,dimension(:),allocatable :: vertices
      logical,dimension(:),allocatable :: vertex_signs

      if (.not.allocated(current_list_local(ic)%vertices)) then
         allocate(current_list_local(ic)%vertices(1))
         allocate(current_list_local(ic)%vertex_sign(1))
      elseif (current_list_local(ic)%n_vert.eq.size(current_list_local(ic)%vertices)) then
         old_capacity=size(current_list_local(ic)%vertices)
         new_capacity=max(1,2*old_capacity)
         allocate(vertices(new_capacity),vertex_signs(new_capacity))
         vertices(1:current_list_local(ic)%n_vert)=&
              current_list_local(ic)%vertices(1:current_list_local(ic)%n_vert)
         vertex_signs(1:current_list_local(ic)%n_vert)=&
              current_list_local(ic)%vertex_sign(1:current_list_local(ic)%n_vert)
         call move_alloc(vertices,current_list_local(ic)%vertices)
         call move_alloc(vertex_signs,current_list_local(ic)%vertex_sign)
      endif
      current_list_local(ic)%n_vert=current_list_local(ic)%n_vert+1
      current_list_local(ic)%vertices(current_list_local(ic)%n_vert)=iv
      current_list_local(ic)%vertex_sign(current_list_local(ic)%n_vert)=vertex_sign
    end subroutine append_current_vertex

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
            max_key=max_key+pm%npart
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
      integer :: j,itype
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

    integer function external_chirality(ctype,iorder,ispin)
      implicit none
      integer,intent(in) :: ctype,iorder,ispin
      external_chirality=0
      if (this%imode.eq.2) return
      if (.not.pm%is_chiral_eligible(ctype)) return
      if (abs(ispin).ne.1) return
      if (iorder.le.2) then
         external_chirality=-ispin
      else
         external_chirality=ispin
      endif
    end function external_chirality

    integer function vertex_result_chirality(itype,ctype,coupl)
      implicit none
      integer,intent(in) :: itype,ctype
      real(kind=8),dimension(2),intent(in) :: coupl
      integer :: ichir
      vertex_result_chirality=0
      if (this%imode.eq.2) return
      if (.not.pm%is_chiral_eligible(ctype)) return
      ichir=fermion_input_chirality()
      if (ichir.eq.0) return
      if (.not.weyl_vertex_allowed(itype,ctype,coupl,ichir)) then
         vertex_result_chirality=-99
         return
      endif
      vertex_result_chirality=ichir
    end function vertex_result_chirality

    integer function fermion_input_chirality()
      implicit none
      fermion_input_chirality=0
      if (pm%is_fermion(current_list_local(ic1)%type)) then
         fermion_input_chirality=current_list_local(ic1)%chirality
         if (fermion_input_chirality.ne.0) return
      endif
      if (pm%is_fermion(current_list_local(ic2)%type)) then
         fermion_input_chirality=current_list_local(ic2)%chirality
      endif
    end function fermion_input_chirality

    logical function weyl_vertex_allowed(itype,ctype,coupl,ichir)
      implicit none
      integer,intent(in) :: itype,ctype,ichir
      real(kind=8),dimension(2),intent(in) :: coupl
      integer :: icoupl
      weyl_vertex_allowed=.true.
      if (itype.eq.16) then
         weyl_vertex_allowed=.false.
         return
      endif
      if (itype.eq.10 .or. itype.eq.11 .or. itype.eq.23 .or. itype.eq.24) then
         icoupl=fermion_coupling_index(ctype,ichir)
         if (coupl(icoupl).eq.0d0) weyl_vertex_allowed=.false.
      endif
    end function weyl_vertex_allowed

    integer function fermion_coupling_index(ctype,ichir)
      implicit none
      integer,intent(in) :: ctype,ichir
      if (ctype.gt.0) then
         if (ichir.eq.-1) then
            fermion_coupling_index=1
         else
            fermion_coupling_index=2
         endif
      else
         if (ichir.eq.1) then
            fermion_coupling_index=1
         else
            fermion_coupling_index=2
         endif
      endif
    end function fermion_coupling_index

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
         if (.not.pm%is_singlet(curr%ext_type(i))) then
            all_singlet_current=.false.
            return
         endif
      enddo
    end function all_singlet_current
    logical function is_quark_from_order(io,iproc)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io,iproc
      if ( (io.le.2  .and. pm%is_antiquark(this%processes(io,iproc))) .or. &
           (io.gt.2  .and. pm%is_quark(this%processes(io,iproc)))) then
         is_quark_from_order=.true.
      else
         is_quark_from_order=.false.
      endif
    end function is_quark_from_order
    logical function is_antiquark_from_order(io,iproc)
      ! 'io' should be a label in the colour order
      implicit none
      integer :: io,iproc
      if ( (io.le.2  .and. pm%is_quark(this%processes(io,iproc))) .or. &
           (io.gt.2  .and. pm%is_antiquark(this%processes(io,iproc)))) then
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
         if (pm%is_quark(curr%ext_type(i)).or.pm%is_antiquark(curr%ext_type(i))) then
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
    integer :: n,iunit,ic,iv,isize,iamp,iproc,nterms,itmp
    logical :: include_this_amp
    write (iunit) this%n_cur,this%n_vert,this%imode,this%nColOrd,this%max_pp,this%n_amps,this%nprocs,this%n_sectors
    write (iunit) this%sector_powers(:,1:this%n_sectors)
    write (iunit) this%n_cur_start(1:n)
    write (iunit) this%n_cur_end(1:n)
    write (iunit) this%n_vert_start(2:n-1)
    write (iunit) this%n_vert_end(2:n-1)
    ! current_list
    do isize=1,n-1
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          call this%current_list(ic)%iproc%bitset_write_unformatted(iunit)
          write (iunit) this%current_list(ic)%type,this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%chirality,this%current_list(ic)%n_gs,this%current_list(ic)%n_ew, &
               this%current_list(ic)%open_quark_leg,this%current_list(ic)%ew_pairs, &
               this%current_list(ic)%u1_pairs, &
               this%current_list(ic)%gluon_pairs, &
               this%current_list(ic)%u1_links,this%current_list(ic)%gluon_links, &
               this%current_list(ic)%mass,this%current_list(ic)%width
          write (iunit) this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)
          write (iunit) this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)
          if (isize.eq.1 .or. isize.eq.n)  write (iunit) this%current_list(ic)%order(1),this%current_list(ic)%spin(1)
       enddo
    enddo
    ! interaction_list
    do iv=1,this%n_vert
       itmp=0
       if (allocated(this%interaction_list(iv)%singlet_mv)) &
            itmp=this%interaction_list(iv)%singlet_mv(0)
       write (iunit) this%interaction_list(iv)%type,this%interaction_list(iv)%chirality,&
            this%interaction_list(iv)%n_gs,this%interaction_list(iv)%n_ew,&
            this%interaction_list(iv)%currents(1:2),&
            this%interaction_list(iv)%coupl(1:2),itmp
       if (itmp.gt.0) write (iunit) this%interaction_list(iv)%singlet_mv(1:itmp)
    enddo
    ! momenta array
    write (iunit) this%pp_bin_to_i(1:maskr(n))
    write (iunit) this%pp_i_to_bin(1:this%max_pp)
    ! process specific information
    do iproc=1,this%nprocs
       write (iunit) this%iproc_start(iproc),this%same_flav(iproc),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       write (iunit) this%processes(1:n,iproc)
    enddo
    write(iunit) this%iproc_start(this%nprocs+1)
    ! amp specific information
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          include_this_amp=.true.
          if (allocated(this%include_amp)) include_this_amp=this%include_amp(iamp)
          write (iunit) include_this_amp,this%same_flavour_sum(iamp,1:2),this%same_flavour_sum_operation(iamp,1:2)
          if (allocated(this%spins)) then
             write (iunit) this%spins(1:n,1,iamp)
          else
             write (iunit) [(0,iv=1,n)]
          endif
          write (iunit) this%perm(1:n-this%n_sing(1),iamp)
          if (.not.this%same_flav(iproc)) write (iunit) this%curr2amp(1:2,iamp)
          write (iunit) this%sector_present(iamp,1:this%n_sectors)
          write (iunit) this%sector_sign(iamp,1:this%n_sectors)
          write (iunit) this%sector_curr2amp(:,iamp,1:this%n_sectors)
          if (allocated(this%sector_three_line_partner_curr2amp)) then
             write (iunit) this%sector_three_line_partner_curr2amp(:,iamp,1:this%n_sectors)
          else
             write (iunit) reshape([(0,iv=1,2*this%n_sectors)],[2,this%n_sectors])
          endif
       enddo
    enddo
    nterms=size(this%sector_term_sign)
    write (iunit) nterms
    write (iunit) this%sector_term_start(0:this%n_amps,1:this%n_sectors)
    if (nterms.gt.0) then
       write (iunit) this%sector_term_curr2amp(:,1:nterms)
       write (iunit) this%sector_term_sign(1:nterms)
    endif
    ! Pruned amplitudes can be cached after the selector has been resolved.
    ! Persist the gate as well as the compact sparse roots so a later
    ! same-flavour reconstruction cannot resurrect a removed sector.
    write (iunit) this%sectors_pruned,this%sectors_pruned_empty
    if (this%sectors_pruned) then
       if (.not.allocated(this%sector_retained)) then
          write (*,*) 'Missing retained-sector gate while writing amplitude cache'
          stop 1
       endif
       write (iunit) this%sector_retained(1:this%n_amps,1:this%n_sectors)
    endif
  end subroutine write_init_amps_to_file

  subroutine read_init_amps_from_file(this,n,iunit)
    implicit none
    class(amplitude_QCD) :: this
    integer :: n,iunit,ic,iv,isize,iamp,iproc,itmp,nterms
    this%sectors_pruned_empty=.false.
    this%sectors_pruned=.false.
    call deallocate_all()
    read (iunit) this%n_cur,this%n_vert,this%imode,this%nColOrd,this%max_pp,this%n_amps,this%nprocs,this%n_sectors
    allocate(this%sector_powers(2,this%n_sectors))
    read (iunit) this%sector_powers
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
          call this%current_list(ic)%iproc%bitset_read_unformatted(iunit)
          read (iunit) this%current_list(ic)%type,this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%chirality,this%current_list(ic)%n_gs,this%current_list(ic)%n_ew,&
               this%current_list(ic)%open_quark_leg,this%current_list(ic)%ew_pairs,&
               this%current_list(ic)%u1_pairs,&
               this%current_list(ic)%gluon_pairs,&
               this%current_list(ic)%u1_links,this%current_list(ic)%gluon_links,&
               this%current_list(ic)%mass,this%current_list(ic)%width
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
       read (iunit) this%interaction_list(iv)%type,this%interaction_list(iv)%chirality,&
            this%interaction_list(iv)%n_gs,this%interaction_list(iv)%n_ew,&
            this%interaction_list(iv)%currents(1:2),this%interaction_list(iv)%coupl(1:2),itmp
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
    allocate(this%n_qqbar(1:this%nprocs))
    allocate(this%n_sing(1:this%nprocs))
    allocate(this%processes(1:n,1:this%nprocs))
    do iproc=1,this%nprocs
       read (iunit) this%iproc_start(iproc),this%same_flav(iproc),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       read (iunit) this%processes(1:n,iproc)
    enddo
    read(iunit) this%iproc_start(this%nprocs+1)
    ! amp specific information
    allocate(this%include_amp(1:this%n_amps))
    allocate(this%same_flavour_sum(1:this%n_amps,1:2))
    allocate(this%same_flavour_sum_operation(1:this%n_amps,1:2))
    allocate(this%spins(1:n,1,1:this%n_amps))
    allocate(this%perm(1:n-this%n_sing(1),1:this%n_amps))
    do iproc=1,this%nprocs
       if (this%same_flav(iproc)) exit
    enddo
    allocate(this%curr2amp(1:2,1:this%iproc_start(iproc)-1))
    allocate(this%sector_present(this%n_amps,this%n_sectors))
    allocate(this%sector_sign(this%n_amps,this%n_sectors))
    allocate(this%sector_curr2amp(2,this%n_amps,this%n_sectors))
    allocate(this%sector_three_line_partner_curr2amp(2,this%n_amps,this%n_sectors))
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          read (iunit) this%include_amp(iamp),this%same_flavour_sum(iamp,1:2),this%same_flavour_sum_operation(iamp,1:2)
          read (iunit) this%spins(1:n,1,iamp)
          read (iunit) this%perm(1:n-this%n_sing(1),iamp)
          if (.not.this%same_flav(iproc)) read (iunit) this%curr2amp(1:2,iamp)
          read (iunit) this%sector_present(iamp,1:this%n_sectors)
          read (iunit) this%sector_sign(iamp,1:this%n_sectors)
          read (iunit) this%sector_curr2amp(:,iamp,1:this%n_sectors)
          read (iunit) this%sector_three_line_partner_curr2amp(:,iamp,1:this%n_sectors)
       enddo
    enddo
    read (iunit) nterms
    allocate(this%sector_term_start(0:this%n_amps,this%n_sectors))
    allocate(this%sector_term_curr2amp(2,nterms))
    allocate(this%sector_term_sign(nterms))
    read (iunit) this%sector_term_start
    if (nterms.gt.0) then
       read (iunit) this%sector_term_curr2amp
       read (iunit) this%sector_term_sign
    endif
    read (iunit) this%sectors_pruned,this%sectors_pruned_empty
    if (this%sectors_pruned) then
       allocate(this%sector_retained(this%n_amps,this%n_sectors))
       read (iunit) this%sector_retained
    endif
    if (this%sectors_pruned_empty .neqv. &
         (nterms.eq.0 .and. .not.any(this%sector_present))) then
       write (*,*) 'Inconsistent empty coupling-sector state in amplitude cache'
       stop 1
    endif
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
      if (allocated(this%n_qqbar)) deallocate(this%n_qqbar)
      if (allocated(this%n_sing)) deallocate(this%n_sing)
      if (allocated(this%processes)) deallocate(this%processes)
      if (allocated(this%include_amp)) deallocate(this%include_amp)
      if (allocated(this%same_flavour_sum)) deallocate(this%same_flavour_sum)
      if (allocated(this%same_flavour_sum_operation)) deallocate(this%same_flavour_sum_operation)
      if (allocated(this%spins)) deallocate(this%spins)
      if (allocated(this%perm)) deallocate(this%perm)
      if (allocated(this%curr2amp)) deallocate(this%curr2amp)
      if (allocated(this%sector_powers)) deallocate(this%sector_powers)
      if (allocated(this%sector_curr2amp)) deallocate(this%sector_curr2amp)
      if (allocated(this%sector_three_line_partner_curr2amp)) &
           deallocate(this%sector_three_line_partner_curr2amp)
      if (allocated(this%sector_present)) deallocate(this%sector_present)
      if (allocated(this%sector_retained)) deallocate(this%sector_retained)
      if (allocated(this%sector_sign)) deallocate(this%sector_sign)
      if (allocated(this%sector_term_start)) deallocate(this%sector_term_start)
      if (allocated(this%sector_term_curr2amp)) deallocate(this%sector_term_curr2amp)
      if (allocated(this%sector_term_sign)) deallocate(this%sector_term_sign)
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
    integer :: ic,iv,isize,ih_in,ifinal,dim
    logical :: read_file
    if (this%sectors_pruned_empty) then
       if (allocated(this%amps)) deallocate(this%amps)
       if (allocated(this%amps_by_order)) deallocate(this%amps_by_order)
       if (allocated(this%amps_r)) deallocate(this%amps_r)
       allocate(this%amps(this%n_amps))
       allocate(this%amps_by_order(this%n_amps,this%n_sectors))
       this%amps=(0d0,0d0)
       this%amps_by_order=(0d0,0d0)
       if (use_real_gluons .and. all(this%n_qqbar(1:this%nprocs).eq.0)) then
          allocate(this%amps_r(this%n_amps))
          this%amps_r=0d0
       endif
       return
    endif
    if (.not. allocated(this%current_list(1)%val_c) .and. .not.allocated(this%current_list(1)%val_r)) then
       do ic=1,this%n_cur
          if (use_real_gluons .and. &
               (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. &
               pm%is_auxiliary_tensor(this%current_list(ic)%type))) then
             dim=current_dim(ic)
             allocate(this%current_list(ic)%val_r(1:dim))
          else
             dim=current_dim(ic)
             allocate(this%current_list(ic)%val_c(1:dim))
          endif
       enddo
       do iv=1,this%n_vert
          if (use_real_gluons .and. &
               (this%interaction_list(iv)%type.ge.0 .and. this%interaction_list(iv)%type.le.3)) then
             dim=interaction_dim(iv)
             allocate(this%interaction_list(iv)%val_r(1:dim))
          else
             dim=interaction_dim(iv)
             allocate(this%interaction_list(iv)%val_c(1:dim))
          endif
       enddo
       if (.not.allocated(this%amps)) allocate(this%amps(1:this%n_amps))
       if (.not.allocated(this%amps_by_order)) &
            allocate(this%amps_by_order(1:this%n_amps,1:this%n_sectors))
       if (use_real_gluons .and. this%n_qqbar(1).eq.0) then
          if (.not. allocated(this%amps_r)) allocate(this%amps_r(1:this%n_amps))
       endif
    endif
    
    call fill_momentum_array()

    do isize=1,n-1
       if (isize.eq.1) then
          ! fill the external wave_functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize)
             ifinal=1
             if (this%current_list(ic)%spin(1).eq.-9) then
                ih_in=hel(this%current_list(ic)%order(1))
             else
                ih_in=this%current_list(ic)%spin(1)
             endif
             if (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. pm%is_photon(this%current_list(ic)%type)) then
                if (use_real_gluons) then
                   call ext_massless_vector_real(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_r(1:4))
                else
                   call ext_massless_vector_cmplx(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4))
                endif
             elseif (pm%is_quark(this%current_list(ic)%type) .or. &
                  pm%is_lepton(this%current_list(ic)%type)) then
                if (this%current_list(ic)%chirality.ne.0) then
                   call ext_fermion_outflow_weyl(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:2),this%current_list(ic)%chirality)
                else
                   call ext_fermion_outflow(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
                endif
             elseif (pm%is_antiquark(this%current_list(ic)%type) .or. &
                  pm%is_antilepton(this%current_list(ic)%type)) then
                if (this%current_list(ic)%chirality.ne.0) then
                   call ext_fermion_inflow_weyl(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:2),this%current_list(ic)%chirality)
                else
                   call ext_fermion_inflow(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
                endif
             elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
                call ext_massive_vector(this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (pm%is_higgs(this%current_list(ic)%type)) then
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
                call TwoGluonToAuxTensor_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     this%interaction_list(iv)%val_r(1:6))
             else
                call TwoGluonToAuxTensor(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                     this%interaction_list(iv)%val_c(1:6))
             endif

          elseif(this%interaction_list(iv)%type.eq.2) then
             if (use_real_gluons) then
                call AuxTensorGluonToGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:6),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     this%interaction_list(iv)%val_r(1:4))
             else
                call AuxTensorGluonToGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:6),&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                     this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.3) then
             if (use_real_gluons) then
                call GluonAuxTensorToGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                             this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:6),&
                                             this%interaction_list(iv)%val_r(1:4))
             else
                call GluonAuxTensorToGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:6),&
                                        this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.4) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                if (use_real_gluons) then
                   write (*,*) 'Weyl fermions with real gluons are not implemented'
                   stop 1
                endif
                call ColourFlowVectorQuarkToQuark_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                            this%interaction_list(iv)%val_c(1),&
                                            this%interaction_list(iv)%chirality)
             elseif (use_real_gluons) then
                call ColourFlowVectorQuarkToQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                            this%interaction_list(iv)%val_c(1:4))
             else
                call ColourFlowVectorQuarkToQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                       this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.5) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                if (use_real_gluons) then
                   write (*,*) 'Weyl fermions with real gluons are not implemented'
                   stop 1
                endif
                call ColourFlowVectorAntiquarkToAntiquark_weyl(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                              this%interaction_list(iv)%val_c(1),&
                                              this%interaction_list(iv)%chirality)
             elseif (use_real_gluons) then
                call ColourFlowVectorAntiquarkToAntiquark_real(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                              this%interaction_list(iv)%val_c(1:4))
             else
                 call ColourFlowVectorAntiquarkToAntiquark(&
                      this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:4))
             endif
          elseif(this%interaction_list(iv)%type.eq.6) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                if (use_real_gluons) then
                   write (*,*) 'Weyl fermions with real gluons are not implemented'
                   stop 1
                endif
                call QuarkColourFlowVectorToQuark_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                            this%interaction_list(iv)%val_c(1),&
                                            this%interaction_list(iv)%chirality)
             elseif (use_real_gluons) then
                call QuarkColourFlowVectorToQuark_real(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                            this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                            this%interaction_list(iv)%val_c(1:4))
             else
                call QuarkColourFlowVectorToQuark(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                       this%interaction_list(iv)%val_c(1:4))
             endif
          elseif(this%interaction_list(iv)%type.eq.7) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                if (use_real_gluons) then
                   write (*,*) 'Weyl fermions with real gluons are not implemented'
                   stop 1
                endif
                call AntiquarkColourFlowVectorToAntiquark_weyl(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                              this%interaction_list(iv)%val_c(1),&
                                              this%interaction_list(iv)%chirality)
             elseif (use_real_gluons) then
                call AntiquarkColourFlowVectorToAntiquark_real(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                              this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                                              this%interaction_list(iv)%val_c(1:4))
             else
                call AntiquarkColourFlowVectorToAntiquark(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:4))
             endif
                 
          elseif(this%interaction_list(iv)%type.eq.8) then
             if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0 .or. &
                 this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                call QuarkAntiquarkToColourFlowU1Vector_weyl(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                             this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                             this%interaction_list(iv)%val_c(1:4),&
                                             this%interaction_list(iv)%coupl(1:2),&
                                             this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                             this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call QuarkAntiquarkToColourFlowU1Vector(&
                     this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                        this%interaction_list(iv)%val_c(1:4),&
                                        this%interaction_list(iv)%coupl(1:2))
             endif

          elseif(this%interaction_list(iv)%type.eq.9) then
             if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0 .or. &
                 this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                call AntiquarkQuarkToGluon_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                             this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                             this%interaction_list(iv)%val_c(1:4),&
                                             this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                             this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call AntiquarkQuarkToGluon(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                        this%interaction_list(iv)%val_c(1:4))
             endif

          elseif(this%interaction_list(iv)%type.eq.10) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                call FermionVectorToFermion_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                                  this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                                  this%interaction_list(iv)%val_c(1),&
                                                  this%interaction_list(iv)%coupl(1:2),&
                                                  this%interaction_list(iv)%chirality)
             else
                if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0) then
                   call FermionVectorToFermion_mixed(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2),&
                        this%current_list(this%interaction_list(iv)%currents(1))%chirality)
                else
                   call FermionVectorToFermion(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2))
                endif
             endif
          elseif(this%interaction_list(iv)%type.eq.11) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                call AntifermionVectorToAntifermion_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                                    this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                                    this%interaction_list(iv)%val_c(1),&
                                                    this%interaction_list(iv)%coupl(1:2),&
                                                    this%interaction_list(iv)%chirality)
             else
                if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0) then
                   call AntifermionVectorToAntifermion_mixed(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2),&
                        this%current_list(this%interaction_list(iv)%currents(1))%chirality)
                else
                   call AntifermionVectorToAntifermion(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2))
                endif
             endif
          elseif (this%interaction_list(iv)%type.eq.12) then
             call VectorVectorToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                       this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)),&
                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                       this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)),&
                       this%interaction_list(iv)%val_c(1:4),&
                       this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.13) then
             call VectorVectorToAuxTensor(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                         this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                         this%interaction_list(iv)%val_c(1:6),&
                                         this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.14) then
             call AuxTensorVectorToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:6),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                           this%interaction_list(iv)%val_c(1:4),&
                                           this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.15) then
             call VectorAuxTensorToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                           this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:6),&
                                           this%interaction_list(iv)%val_c(1:4),&
                                           this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.16) then
             call FermionScalarToFermion(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                     this%interaction_list(iv)%val_c(1:4),&
                                     this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.17) then
             call VectorVectorToScalar(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                     this%interaction_list(iv)%val_c(1),&
                                     this%interaction_list(iv)%coupl(1:2))

          elseif(this%interaction_list(iv)%type.eq.18) then
             call ScalarVectorToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                                     this%interaction_list(iv)%val_c(1:4),&
                                     this%interaction_list(iv)%coupl)

          elseif(this%interaction_list(iv)%type.eq.19) then
             call VectorScalarToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                                     this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                     this%interaction_list(iv)%val_c(1:4),&
                                     this%interaction_list(iv)%coupl)

          elseif(this%interaction_list(iv)%type.eq.20) then
             call ScalarScalarToScalar(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                       this%interaction_list(iv)%val_c(1),&
                                       this%interaction_list(iv)%coupl(1:2))


          elseif(this%interaction_list(iv)%type.eq.21) then
             if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0 .or. &
                 this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                call FermionAntifermionToVector_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                               this%interaction_list(iv)%val_c(1:4),&
                                               this%interaction_list(iv)%coupl(1:2),&
                                               this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                               this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call FermionAntifermionToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                          this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                          this%interaction_list(iv)%val_c(1:4),&
                                          this%interaction_list(iv)%coupl(1:2))
             endif

          elseif(this%interaction_list(iv)%type.eq.22) then
             if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0 .or. &
                 this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                call AntifermionFermionToVector_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                               this%interaction_list(iv)%val_c(1:4),&
                                               this%interaction_list(iv)%coupl(1:2),&
                                               this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                               this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call AntifermionFermionToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                          this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                          this%interaction_list(iv)%val_c(1:4),&
                                          this%interaction_list(iv)%coupl(1:2))
             endif

          elseif(this%interaction_list(iv)%type.eq.23) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                call VectorFermionToFermion_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                                  this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                                  this%interaction_list(iv)%val_c(1),&
                                                  this%interaction_list(iv)%coupl(1:2),&
                                                  this%interaction_list(iv)%chirality)
             else
                if (this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                   call VectorFermionToFermion_mixed(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2),&
                        this%current_list(this%interaction_list(iv)%currents(2))%chirality)
                else
                   call VectorFermionToFermion(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2))
                endif
             endif
          elseif(this%interaction_list(iv)%type.eq.24) then
             if (this%interaction_list(iv)%chirality.ne.0) then
                call VectorAntifermionToAntifermion_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                                    this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                                    this%interaction_list(iv)%val_c(1),&
                                                    this%interaction_list(iv)%coupl(1:2),&
                                                    this%interaction_list(iv)%chirality)
             else
                if (this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                   call VectorAntifermionToAntifermion_mixed(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2),&
                        this%current_list(this%interaction_list(iv)%currents(2))%chirality)
                else
                   call VectorAntifermionToAntifermion(&
                        this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                        this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                        this%interaction_list(iv)%val_c(1:4),&
                        this%interaction_list(iv)%coupl(1:2))
                endif
             endif

          else
             write (*,*) 'Unknown vertex type: not yet implemented',iv,this%interaction_list(iv)%type
             stop 1
          endif
       enddo

       ! compute the currents by combining the interactions
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          if (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. pm%is_photon(this%current_list(ic)%type)) then
             call combine_interactions(4)
             ! a massless-vector current
             if (isize.ne.n-1)  then
                call include_massless_vector_propagator()
             endif
          elseif (pm%is_quark(this%current_list(ic)%type).or.pm%is_lepton(this%current_list(ic)%type)) then
             ! a fermion current
             if (this%current_list(ic)%chirality.ne.0) then
                call combine_interactions(2)
             else
                call combine_interactions(4)
             endif
             if (isize.ne.n-1)  then
                if (this%current_list(ic)%chirality.ne.0) then
                   call include_fermion_propagator_weyl()
                else
                   call include_fermion_propagator()
                endif
             endif
          elseif (pm%is_auxiliary_tensor(this%current_list(ic)%type)) then
             ! the non-propagating tensor current
             call combine_interactions(6)
          elseif (pm%is_antiquark(this%current_list(ic)%type).or.pm%is_antilepton(this%current_list(ic)%type)) then
             ! an antifermion current
             if (this%current_list(ic)%chirality.ne.0) then
                call combine_interactions(2)
             else
                call combine_interactions(4)
             endif
             if (isize.ne.n-1)  then
                if (this%current_list(ic)%chirality.ne.0) then
                   call include_antifermion_propagator_weyl()
                else
                   call include_antifermion_propagator()
                endif
             endif
          elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
             call combine_interactions(4)
             ! a massive vector boson current
             if (isize.ne.n-1)  then
                call include_massive_vector_propagator()
             endif
          elseif (pm%is_higgs(this%current_list(ic)%type)) then
             ! a scalar current
             call combine_interactions(1)
             if (isize.ne.n-1)  then
                call include_scalar_propagator()
             endif
          elseif (pm%is_auxiliary_scalar(this%current_list(ic)%type)) then
             ! a non-propagating scalar current
             call combine_interactions(1)

          else
             write (*,*) 'Unknown current type',ic,this%current_list(ic)%type
             stop 1
          endif
       enddo
    enddo

    call compute_amps_from_currents

  contains

    integer function current_dim(icur)
      implicit none
      integer,intent(in) :: icur
      if (this%current_list(icur)%chirality.ne.0) then
         current_dim=2
      else
         current_dim=pm%get_dim(this%current_list(icur)%type)
      endif
    end function current_dim

    integer function interaction_dim(ivert)
      implicit none
      integer,intent(in) :: ivert
      if (this%interaction_list(ivert)%chirality.ne.0) then
         interaction_dim=2
      else
         interaction_dim=pm%get_inter_dim(this%interaction_list(ivert)%type)
      endif
    end function interaction_dim

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
      integer :: iamp,iproc,idau,isector,iterm
      this%amps_by_order=(0d0,0d0)
      do iproc=1,this%nprocs
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            if (this%same_flav(iproc)) cycle
            do isector=1,this%n_sectors
               if (.not.this%sector_present(iamp,isector)) cycle
               if (allocated(this%sector_term_start)) then
                  do iterm=this%sector_term_start(iamp-1,isector)+1,&
                       this%sector_term_start(iamp,isector)
                     this%amps_by_order(iamp,isector)=this%amps_by_order(iamp,isector)+&
                          dble(this%sector_term_sign(iterm))*contract_current_pair(&
                          this%sector_term_curr2amp(1,iterm),&
                          this%sector_term_curr2amp(2,iterm))
                  enddo
               else
                  this%amps_by_order(iamp,isector)=contract_sector(iamp,isector)
                  if (allocated(this%sector_three_line_partner_curr2amp)) then
                     if (this%sector_three_line_partner_curr2amp(1,iamp,isector).ne.0) then
                        this%amps_by_order(iamp,isector)=this%amps_by_order(iamp,isector)+&
                             contract_current_pair(&
                             this%sector_three_line_partner_curr2amp(1,iamp,isector),&
                             this%sector_three_line_partner_curr2amp(2,iamp,isector))
                     endif
                  endif
                  this%amps_by_order(iamp,isector)=dble(this%sector_sign(iamp,isector))*&
                       this%amps_by_order(iamp,isector)
               endif
               if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. &
                    iamp.gt.this%n_amps/2 .and. mod(n,2).eq.1) &
                    this%amps_by_order(iamp,isector)=-this%amps_by_order(iamp,isector)
            enddo
         enddo
      enddo
      ! Same-flavour amplitudes are linear combinations of physical amplitudes;
      ! apply exactly the same operation independently in every sector.
      do iproc=1,this%nprocs
         if (.not.this%same_flav(iproc)) cycle
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            do isector=1,this%n_sectors
               if (this%sectors_pruned) then
                  if (.not.this%sector_retained(iamp,isector)) then
                     this%sector_present(iamp,isector)=.false.
                     cycle
                  endif
               endif
               this%sector_present(iamp,isector)=.false.
               do idau=1,2
                  if (this%same_flavour_sum(iamp,idau).gt.0) then
                     this%amps_by_order(iamp,isector)=this%amps_by_order(iamp,isector)+&
                          apply_sector_operation(iamp,idau,isector)
                     if (this%sector_present(this%same_flavour_sum(iamp,idau),isector)) &
                          this%sector_present(iamp,isector)=.true.
                  endif
               enddo
            enddo
         enddo
      enddo
      this%amps=sum(this%amps_by_order,dim=2)
      if (use_real_gluons .and. all(this%n_qqbar(1:this%nprocs).eq.0)) this%amps_r=dble(this%amps)
    end subroutine compute_amps_from_currents

    complex(kind=8) function contract_sector(iamp,isector)
      implicit none
      integer,intent(in) :: iamp,isector
      integer :: ic1,ic2
      ic1=this%sector_curr2amp(1,iamp,isector)
      ic2=this%sector_curr2amp(2,iamp,isector)
      if (use_real_gluons .and. all(this%n_qqbar(1:this%nprocs).eq.0)) then
         contract_sector=cmplx(sum(this%current_list(ic1)%val_r(1:4)*&
              this%current_list(ic2)%val_r(1:4)),0d0,kind=8)
      else
         contract_sector=contract_current_pair(ic1,ic2)
      endif
    end function contract_sector

    complex(kind=8) function contract_currents(iamp)
      implicit none
      integer,intent(in) :: iamp
      integer :: ic1,ic2
      ic1=this%curr2amp(1,iamp)
      ic2=this%curr2amp(2,iamp)
      contract_currents=contract_current_pair(ic1,ic2)
    end function contract_currents

    complex(kind=8) function contract_current_pair(ic1,ic2)
      implicit none
      integer,intent(in) :: ic1,ic2
      contract_current_pair=ContractFermionCurrents(&
           this%current_list(ic1)%val_c,this%current_list(ic1)%chirality,&
           this%current_list(ic2)%val_c,this%current_list(ic2)%chirality)
    end function contract_current_pair

    complex (kind=8) function apply_operation(iamp,idau)
      implicit none
      integer,intent(in) :: iamp,idau
      apply_operation=this%amps(this%same_flavour_sum(iamp,idau))
      if (btest(this%same_flavour_sum_operation(iamp,idau),2)) apply_operation=cmplx(aimag(apply_operation),dble(apply_operation))
      if (btest(this%same_flavour_sum_operation(iamp,idau),0)) apply_operation=-conjg(apply_operation)
      if (btest(this%same_flavour_sum_operation(iamp,idau),1)) apply_operation=conjg(apply_operation)
    end function apply_operation

    complex(kind=8) function apply_sector_operation(iamp,idau,isector)
      implicit none
      integer,intent(in) :: iamp,idau,isector
      apply_sector_operation=this%amps_by_order(this%same_flavour_sum(iamp,idau),isector)
      if (btest(this%same_flavour_sum_operation(iamp,idau),2)) &
           apply_sector_operation=cmplx(aimag(apply_sector_operation),dble(apply_sector_operation))
      if (btest(this%same_flavour_sum_operation(iamp,idau),0)) &
           apply_sector_operation=-conjg(apply_sector_operation)
      if (btest(this%same_flavour_sum_operation(iamp,idau),1)) &
           apply_sector_operation=conjg(apply_sector_operation)
    end function apply_sector_operation
    
    subroutine combine_interactions(dim)
      implicit none
      integer :: dim,iv
      if (use_real_gluons .and. &
           (pm%is_colour_flow_vector(this%current_list(ic)%type).or. &
           pm%is_gluon_aux_tensor(this%current_list(ic)%type))) then
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
    subroutine include_massless_vector_propagator()
      implicit none
      if (use_real_gluons) then
         call MasslessVectorPropagator_real(this%current_list(ic)%val_r, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      else
         call MasslessVectorPropagator(this%current_list(ic)%val_c, &
              this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)))
      endif
    end subroutine include_massless_vector_propagator

    subroutine include_massive_vector_propagator()
      implicit none
      call MassiveVectorPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%mass,this%current_list(ic)%width)
    end subroutine include_massive_vector_propagator

    subroutine include_fermion_propagator()
      implicit none
      call FermionPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), & 
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_fermion_propagator

    subroutine include_fermion_propagator_weyl()
      implicit none
      call FermionPropagator_weyl(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)), &
           this%current_list(ic)%chirality)
    end subroutine include_fermion_propagator_weyl

    subroutine include_antifermion_propagator()
      implicit none
      call AntifermionPropagator(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_antifermion_propagator
    subroutine include_antifermion_propagator_weyl()
      implicit none
      call AntifermionPropagator_weyl(this%current_list(ic)%val_c, &
           this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin)),&
           this%current_list(ic)%chirality)
    end subroutine include_antifermion_propagator_weyl
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
    integer :: col_acc,n,iperm,jperm,ival,iacc,isum,gi,gj,ui,uj,key &
         ,max_keys,jperm_lower,n_unique_rows,irow,iunique,iproc,ioff &
         ,nOrd
    integer,dimension(n) :: iper,jper,part
    integer,dimension(:),allocatable :: n_vals
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(:,:),allocatable :: diff_vals
    real(kind=8),dimension(:,:,:),allocatable :: col_vals
    integer,dimension(:,:),allocatable :: ic,ir,n_colour_elements,unique_rows
    integer(kind=8),dimension(:),allocatable :: perm_dict

    write (99,*) 'Initialising colour matrix ...'
    if (this%nprocs.eq.1) then
       iproc=1
    elseif(this%nprocs.eq.3) then
       iproc=3
    else
       write (*,*) 'computation of color factor only for a single process at the time',this%nprocs
       stop 1
    endif
    part(1:n)=this%processes(1:n,1)
    ioff=this%iproc_start(iproc)-1
    nOrd=n-this%n_sing(iproc)

    if (this%n_qqbar(iproc).eq.3) then
       if (col_acc.eq.0) then
          call init_three_quark_line_lc_diagonal()
       else
          call init_three_quark_line_colour_matrix()
       endif
       write (99,*) '... colour matrix initialised'
       return
    endif

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
    write (99,*) 'A single row in the colour matrix has',n_vals(1:3),&
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
             call get_col_fac(iper,jper,ui,gi,col_fac)
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

    write (99,*) '... colour matrix initialised'
  contains
    subroutine init_three_quark_line_lc_diagonal()
      ! At leading colour the physical three-line basis is diagonal with one
      ! common coefficient. Store only those rows rather than materialising
      ! the quadratic NLC/full-colour Gram matrix.
      implicit none
      integer :: row,nrows
      integer,dimension(:),allocatable :: flow_sign
      integer,dimension(:,:),allocatable :: endpoints,ngluons
      integer,dimension(:,:,:),allocatable :: gluons
      real(kind=8),dimension(3) :: three_line_col_fac

      nrows=this%nColOrd
      if (this%n_amps.ne.nrows) then
         write (*,*) 'Three-line colour basis does not match generated amplitudes',&
              this%n_amps,nrows
         stop 1
      endif
      allocate(endpoints(3,nrows))
      allocate(ngluons(3,nrows))
      allocate(gluons(nOrd,3,nrows))
      allocate(flow_sign(nrows))
      call build_three_line_flow_metadata(endpoints,ngluons,gluons,flow_sign)
      call compute_three_line_color_factor(1,1,endpoints,ngluons,gluons,&
           flow_sign,three_line_col_fac)
      if (three_line_col_fac(1).eq.0d0) then
         write (*,*) 'Three-line leading-colour diagonal is zero'
         stop 1
      endif

      allocate(this%n_col_vals(3))
      this%n_col_vals=(/1,0,0/)
      allocate(this%diff_col_vals(1,3))
      this%diff_col_vals=0d0
      this%diff_col_vals(1,1)=three_line_col_fac(1)
      allocate(this%i_col_i(1,3))
      this%i_col_i=0
      this%i_col_i(1,1)=1
      allocate(this%row_index(0:nrows,1,3))
      this%row_index=0
      allocate(this%col_index(nrows+1))
      this%col_index=0
      do row=1,nrows
         this%col_index(row+1)=row
         this%row_index(row,1,1)=row
      enddo
    end subroutine init_three_quark_line_lc_diagonal

    subroutine init_three_quark_line_colour_matrix()
      ! Sew each pair of open fundamental-string flows into closed traces and
      ! let color_algebra perform the exact SU(Nc) reduction.
      implicit none
      integer :: row,col,iacc_local,ival_local,nrows,max_pairs,max_nvals,total_entries
      integer :: offset
      integer,dimension(:),allocatable :: flow_sign,nvalues
      integer,dimension(:,:),allocatable :: endpoints,ngluons,counts,cursor
      integer,dimension(:,:,:),allocatable :: gluons
      real(kind=8),dimension(:,:,:),allocatable :: factors
      real(kind=8),dimension(:,:),allocatable :: values
      real(kind=8),dimension(3) :: three_line_col_fac

      nrows=this%nColOrd
      if (int(nrows,kind=8).gt.max_three_line_color_orders) then
         write (*,*) 'Three-line colour matrix exceeds supported size',&
              nrows,max_three_line_color_orders
         stop 1
      endif
      if (this%n_amps.ne.nrows) then
         write (*,*) 'Three-line colour basis does not match generated amplitudes',&
              this%n_amps,nrows
         stop 1
      endif

      allocate(endpoints(3,nrows))
      allocate(ngluons(3,nrows))
      allocate(gluons(nOrd,3,nrows))
      allocate(flow_sign(nrows))
      call build_three_line_flow_metadata(endpoints,ngluons,gluons,flow_sign)

      allocate(factors(nrows,nrows,3))
      factors=0d0
      do row=1,nrows
         if (use_symm_cm) then
            col=row
         else
            col=1
         endif
         do while (col.le.nrows)
            call compute_three_line_color_factor(row,col,endpoints,ngluons,&
                 gluons,flow_sign,three_line_col_fac)
            if (use_symm_cm .and. row.ne.col) &
                 three_line_col_fac=2d0*three_line_col_fac
            factors(row,col,1:3)=three_line_col_fac(1:3)
            col=col+1
         enddo
      enddo

      if (use_symm_cm) then
         max_pairs=nrows*(nrows+1)/2
      else
         max_pairs=nrows*nrows
      endif
      allocate(nvalues(3))
      allocate(values(max_pairs,3))
      allocate(counts(max_pairs,3))
      nvalues=0
      values=0d0
      counts=0
      do row=1,nrows
         if (use_symm_cm) then
            col=row
         else
            col=1
         endif
         do while (col.le.nrows)
            do iacc_local=1,3
               if (factors(row,col,iacc_local).eq.0d0) cycle
               do ival_local=1,nvalues(iacc_local)
                  if (factors(row,col,iacc_local).eq.values(ival_local,iacc_local)) exit
               enddo
               if (ival_local.eq.nvalues(iacc_local)+1) then
                  nvalues(iacc_local)=ival_local
                  values(ival_local,iacc_local)=factors(row,col,iacc_local)
               endif
               counts(ival_local,iacc_local)=counts(ival_local,iacc_local)+1
            enddo
            col=col+1
         enddo
      enddo

      max_nvals=maxval(nvalues)
      allocate(this%n_col_vals(3))
      this%n_col_vals=nvalues
      allocate(this%diff_col_vals(max_nvals,3))
      allocate(this%i_col_i(max_nvals,3))
      allocate(this%row_index(0:nrows,max_nvals,3))
      this%diff_col_vals=0d0
      this%i_col_i=0
      this%row_index=0
      do iacc_local=1,3
         this%diff_col_vals(1:nvalues(iacc_local),iacc_local)=&
              values(1:nvalues(iacc_local),iacc_local)
      enddo

      total_entries=1+sum(counts)
      allocate(this%col_index(total_entries))
      this%col_index=0
      offset=1
      do iacc_local=1,3
         do ival_local=1,nvalues(iacc_local)
            this%i_col_i(ival_local,iacc_local)=offset
            offset=offset+counts(ival_local,iacc_local)
         enddo
      enddo

      allocate(cursor(max_nvals,3))
      cursor=0
      do row=1,nrows
         if (use_symm_cm) then
            col=row
         else
            col=1
         endif
         do while (col.le.nrows)
            do iacc_local=1,3
               if (factors(row,col,iacc_local).eq.0d0) cycle
               do ival_local=1,nvalues(iacc_local)
                  if (factors(row,col,iacc_local).eq.values(ival_local,iacc_local)) exit
               enddo
               cursor(ival_local,iacc_local)=cursor(ival_local,iacc_local)+1
               this%col_index(this%i_col_i(ival_local,iacc_local)+&
                    cursor(ival_local,iacc_local))=col
            enddo
            col=col+1
         enddo
         do iacc_local=1,3
            this%row_index(row,1:nvalues(iacc_local),iacc_local)=&
                 cursor(1:nvalues(iacc_local),iacc_local)
         enddo
      enddo
    end subroutine init_three_quark_line_colour_matrix

    subroutine build_three_line_flow_metadata(endpoints,ngluons,gluons,flow_sign)
      implicit none
      integer,dimension(3,this%nColOrd),intent(out) :: endpoints,ngluons
      integer,dimension(nOrd,3,this%nColOrd),intent(out) :: gluons
      integer,dimension(this%nColOrd),intent(out) :: flow_sign
      integer :: row,line,qrank,arank,i,j,inversions
      integer,dimension(3) :: ref_q,ref_aq,row_q,row_aq,row_ngluons
      integer,dimension(nOrd,3) :: row_gluons

      endpoints=0
      ngluons=0
      gluons=0
      call parse_three_line_word(this%perm(1:nOrd,ioff+1),ref_q,ref_aq,&
           row_ngluons,row_gluons)
      do row=1,this%nColOrd
         call parse_three_line_word(this%perm(1:nOrd,ioff+row),row_q,row_aq,&
              row_ngluons,row_gluons)
         do line=1,3
            qrank=0
            arank=0
            do i=1,3
               if (row_q(line).eq.ref_q(i)) qrank=i
               if (row_aq(line).eq.ref_aq(i)) arank=i
            enddo
            if (qrank.eq.0 .or. arank.eq.0) then
               write (*,*) 'Inconsistent external labels in three-line colour flow',row
               stop 1
            endif
            endpoints(qrank,row)=arank
            ngluons(qrank,row)=row_ngluons(line)
            if (row_ngluons(line).gt.0) then
               gluons(1:row_ngluons(line),qrank,row)=&
                    row_gluons(1:row_ngluons(line),line)
            endif
         enddo
         inversions=0
         do i=1,2
            do j=i+1,3
               if (endpoints(i,row).gt.endpoints(j,row)) inversions=inversions+1
            enddo
         enddo
         if (mod(inversions,2).eq.0) then
            flow_sign(row)=1
         else
            flow_sign(row)=-1
         endif
      enddo
    end subroutine build_three_line_flow_metadata

    subroutine parse_three_line_word(word,q_labels,aq_labels,line_ngluons,line_gluons)
      implicit none
      integer,dimension(nOrd),intent(in) :: word
      integer,dimension(3),intent(out) :: q_labels,aq_labels,line_ngluons
      integer,dimension(nOrd,3),intent(out) :: line_gluons
      integer :: pos,label,line

      q_labels=0
      aq_labels=0
      line_ngluons=0
      line_gluons=0
      line=0
      do pos=1,nOrd
         label=word(pos)
         if (label_is_quark(label)) then
            if (line.ne.0) then
               if (aq_labels(line).eq.0) then
                  write (*,*) 'Adjacent quarks in three-line colour word',word
                  stop 1
               endif
            endif
            line=line+1
            if (line.gt.3) then
               write (*,*) 'Too many quark strings in three-line colour word',word
               stop 1
            endif
            q_labels(line)=label
         elseif (label_is_antiquark(label)) then
            if (line.eq.0) then
               write (*,*) 'Unmatched antiquark in three-line colour word',word
               stop 1
            endif
            if (q_labels(line).eq.0 .or. aq_labels(line).ne.0) then
               write (*,*) 'Unmatched antiquark in three-line colour word',word
               stop 1
            endif
            aq_labels(line)=label
         elseif (abs(part(label)).eq.21) then
            if (line.eq.0) then
               write (*,*) 'Gluon outside an open quark string',word
               stop 1
            endif
            if (aq_labels(line).ne.0) then
               write (*,*) 'Gluon outside an open quark string',word
               stop 1
            endif
            line_ngluons(line)=line_ngluons(line)+1
            line_gluons(line_ngluons(line),line)=label
         else
            write (*,*) 'Unknown coloured label in three-line colour word',label,word
            stop 1
         endif
      enddo
      if (line.ne.3 .or. any(q_labels.eq.0) .or. any(aq_labels.eq.0)) then
         write (*,*) 'Incomplete three-line colour word',word
         stop 1
      endif
    end subroutine parse_three_line_word

    logical function label_is_quark(label)
      implicit none
      integer,intent(in) :: label
      label_is_quark=(label.le.2 .and. part(label).le.-1 .and. part(label).ge.-6) .or. &
           (label.gt.2 .and. part(label).ge.1 .and. part(label).le.6)
    end function label_is_quark

    logical function label_is_antiquark(label)
      implicit none
      integer,intent(in) :: label
      label_is_antiquark=(label.le.2 .and. part(label).ge.1 .and. part(label).le.6) .or. &
           (label.gt.2 .and. part(label).le.-1 .and. part(label).ge.-6)
    end function label_is_antiquark

    subroutine compute_three_line_color_factor(row,col,endpoints,ngluons,&
         gluons,flow_sign,three_line_col_fac)
      implicit none
      integer,intent(in) :: row,col
      integer,dimension(3,this%nColOrd),intent(in) :: endpoints,ngluons
      integer,dimension(nOrd,3,this%nColOrd),intent(in) :: gluons
      integer,dimension(this%nColOrd),intent(in) :: flow_sign
      real(kind=8),dimension(3),intent(out) :: three_line_col_fac
      integer :: q,qnext,anti,start,line,nloops,pos,power,leading_power,max_power,sgn
      integer,dimension(3) :: trace_len
      integer,dimension(2*nOrd,3) :: trace_word
      logical,dimension(3) :: visited
      real(kind=16) :: exact

      trace_len=0
      trace_word=0
      visited=.false.
      nloops=0
      do start=1,3
         if (visited(start)) cycle
         nloops=nloops+1
         q=start
         do
            if (visited(q)) then
               if (q.eq.start) exit
               write (*,*) 'Malformed three-line colour-flow cycle',row,col
               stop 1
            endif
            visited(q)=.true.
            do pos=1,ngluons(q,row)
               trace_len(nloops)=trace_len(nloops)+1
               trace_word(trace_len(nloops),nloops)=gluons(pos,q,row)
            enddo
            anti=endpoints(q,row)
            qnext=0
            do line=1,3
               if (endpoints(line,col).eq.anti) then
                  qnext=line
                  exit
               endif
            enddo
            if (qnext.eq.0) then
               write (*,*) 'Cannot sew three-line colour-flow endpoints',row,col,anti
               stop 1
            endif
            do pos=ngluons(qnext,col),1,-1
               trace_len(nloops)=trace_len(nloops)+1
               trace_word(trace_len(nloops),nloops)=gluons(pos,qnext,col)
            enddo
            q=qnext
         enddo
      enddo

      call Tr_allocate(nOrd)
      Tr=0
      coef=0d0
      coef_Nc=0
      Tr(0,0,0)=1
      Tr(0,0,1)=nloops
      do line=1,nloops
         Tr(0,line,1)=trace_len(line)
         if (trace_len(line).gt.0) then
            Tr(1:trace_len(line),line,1)=trace_word(1:trace_len(line),line)
         endif
      enddo
      sgn=flow_sign(row)*flow_sign(col)
      coef(1)=dble(sgn)
      coef_Nc(0,1)=sgn
      call Tr_full_simplify(exact)

      leading_power=nOrd-3
      three_line_col_fac=0d0
      if (leading_power.ge.-nOrd .and. leading_power.le.nOrd) then
         three_line_col_fac(1)=dble(coef_Nc(leading_power,0))*3d0**leading_power
      endif
      max_power=-nOrd-1
      do power=nOrd,-nOrd,-1
         if (coef_Nc(power,0).ne.0) then
            max_power=power
            exit
         endif
      enddo
      if (max_power.ge.leading_power-2) three_line_col_fac(2)=dble(exact)
      three_line_col_fac(3)=dble(exact)
      call Tr_deallocate
    end subroutine compute_three_line_color_factor

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
      if (irow.gt.this%nColOrd) then
         write (*,*) 'Could not determine ui and gi correctly'
         write (*,*) this%nColOrd,ui,gi,iunique,nOrd
         write (*,*) part(1:n)
         write (*,*) part(iper(1:n))
         stop 1
      endif
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
      ! check if quarks are connected in the same (or opposite way) as
      ! compared to the very first permutation
      if ((iper(1).eq.this%perm(1,ioff+1) .and. iper(nOrd).eq.this%perm(nOrd,ioff+1)) .or. &
           (iper(1).ne.this%perm(1,ioff+1) .and. iper(nOrd).ne.this%perm(nOrd,ioff+1))) then
         ui=1
      else
         ui=2
      endif
!!$      if (abs(part(iper(1))).eq.abs(part(iper(n-this%n_sing(1))))) then
!!$         ui=1
!!$      else
!!$         ui=2
!!$      endif
    end subroutine determine_ui
    
    subroutine determine_gi_ui(iper,gi,ui)
      implicit none
      integer :: ui,gi
      integer,dimension(n) :: iper
      call determine_gi(iper,gi)
      call determine_ui(iper,ui)
    end subroutine determine_gi_ui
    
   subroutine get_col_fac(iper,jper,ui,gi,col_fac)
     implicit none
     integer,intent(in) :: gi,ui
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
                  col_fac(1)=dble(3**(n-4))  * 9d0 ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
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
!!$            if (.not.this%same_flav(iproc)) then
!!$               call check_NLC_2qqbar(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
!!$               if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
!!$               if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-4)) ! NLC parts
!!$               write (*,*) 'DF',col_fac(2),acc
!!$            else
               call check_NLC_2qqbar_SF(n,iper_ord,jper_ord,gi,gj,ui,uj,acc)
               if (acc.eq.99) col_fac(2)=dble((3)**(n-2))-dble((n-4)*(3)**(n-4)) ! LC interfence
               if (acc.le.1) col_fac(2)=dble(acc*(3)**(n-3)) ! NLC parts
!!$               write (*,*) 'SF',col_fac(2),acc
!!$            endif
            ! include the full expansion
            if (acc.ne.0) then
               call Tr_allocate(n)
               if (ui.eq.uj) then
                  Tr(0,0,0)=1 ! one term
                  Tr(0,0,1)=2 ! two traces
                  Tr(0,1,1) = gi+gj  ! number of generators in first trace
                  Tr(0,2,1) = 2*(n-4)-(gi+gj)  ! number of generators in second trace
                  Tr(1:gi,1,1) = iper(2:2+gi-1)
                  Tr(gi+1:gi+gj,1,1) = jper(2+gj-1:2:-1)
                  Tr(1:n-4-gi,2,1) = iper(2+gi+2:n-1)
                  Tr(n-4-gi+1:2*(n-4)-(gi+gj),2,1) = jper(n-1:2+gj+2:-1)
                  if (ui.eq.2 .and. uj.eq.2 .and. .not.this%same_flav(iproc)) then
                     coef(1)=1d0/9d0 *9d0  ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
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
                     coef(1)=-1d0/3d0 *3d0  ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
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
                  coef(1)=1d0/9d0  *9d0  ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
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
                  coef(1)=-1d0/3d0   * 3d0  ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
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

  subroutine optimise_evaluation(this,n)
    !
    ! Checks all computed currents and checks if some are equal. If
    ! equal, do not recompute, rather re-use already computed values
    !
    ! POTENTIAL OTHER OPTIMISATIONS:
    !
    ! 1. REMOVE INTERACTIONS THAT YIELD ZERO RESULT
    ! 2. INCLUDE MULTIPLICATIVE COUPLING CONSTANT AT LATER STAGE
    ! 3. WEYL SPINORS & SEPARATE VERTEX ROUTINES FOR LEFT AND RIGHT-HANDED INTERACTIONS?
    !      
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    integer :: isize,ic1,ic2,iv1,iv2,i,n_vert,n,ic,iv
    integer,dimension(:,:),allocatable :: map_cur,map_vert
    integer,dimension(n-1) :: identical_curr,identical_vert
    real(kind=8),parameter :: tiny=1d-10
    logical,dimension(:),allocatable :: include_cur,include_vert,reordered_interactions
    type(current),dimension(:),allocatable :: current_list_local
    type(interaction),dimension(:),allocatable :: interaction_list_local
    integer,dimension(:),allocatable :: interactions_map
    if (this%sectors_pruned_empty) then
       write (99,*) 'Skipping amplitude optimisation: no selected coupling sector is present'
       return
    endif
    allocate(include_cur(1:this%n_cur))
    allocate(include_vert(1:this%n_vert))
    include_cur=.true.
    include_vert=.true.
    allocate(map_cur(0:this%n_cur,1:2))
    allocate(map_vert(0:this%n_vert,1:2))
    map_cur(0,1)=0
    map_vert(0,1)=0
    identical_curr=0
    identical_vert=0
    do isize=1,n-1
       do ic1=this%n_cur_start(isize),this%n_cur_end(isize)-1
          if (.not.include_cur(ic1)) cycle
          if (sum(abs(this%current_list(ic1)%val_c)).eq.0d0) then
             map_cur(0,1)=map_cur(0,1)+1
             map_cur(map_cur(0,1),1)=ic1
!!$               map_cur(map_cur(0,1),2)=0
             map_cur(map_cur(0,1),2)=ic1
             include_cur(ic1)=.false.
             cycle
          endif
          do ic2=ic1+1,this%n_cur_end(isize)
             if (.not.include_cur(ic2)) cycle
             if (this%current_list(ic1)%n_gs.ne.this%current_list(ic2)%n_gs) cycle
             if (this%current_list(ic1)%n_ew.ne.this%current_list(ic2)%n_ew) cycle
             if (size(this%current_list(ic1)%val_c).ne.size(this%current_list(ic2)%val_c)) cycle
             if (this%current_list(ic1)%type.ne.this%current_list(ic2)%type) cycle
             if (this%current_list(ic1)%chirality.ne.this%current_list(ic2)%chirality) cycle
             if (this%current_list(ic1)%bin.ne.this%current_list(ic2)%bin) cycle
             if ( sum(abs(this%current_list(ic1)%val_c-this%current_list(ic2)%val_c))/ &
                  sum(abs(this%current_list(ic1)%val_c)+abs(this%current_list(ic2)%val_c)).lt.tiny) then
                map_cur(0,1)=map_cur(0,1)+1
                map_cur(map_cur(0,1),1)=ic2
                map_cur(map_cur(0,1),2)=ic1
                identical_curr(isize)=identical_curr(isize)+1
                include_cur(ic2)=.false.
             endif
          enddo
       enddo
    enddo
    do i=1,map_cur(0,1)
       do iv1=1,this%n_vert
          if (this%interaction_list(iv1)%currents(1).eq.map_cur(i,1)) then
             this%interaction_list(iv1)%currents(1)=map_cur(i,2)
          endif
          if (this%interaction_list(iv1)%currents(2).eq.map_cur(i,1)) then
             this%interaction_list(iv1)%currents(2)=map_cur(i,2)
          endif
       enddo
    enddo
    do isize=2,n-1
       do iv1=this%n_vert_start(isize),this%n_vert_end(isize)-1
          if (.not.include_vert(iv1)) cycle
          if (sum(abs(this%interaction_list(iv1)%val_c)).eq.0d0) then
             map_vert(0,1)=map_vert(0,1)+1
             map_vert(map_vert(0,1),1)=iv1
!!$               map_vert(map_vert(0,1),2)=0
             map_vert(map_vert(0,1),2)=iv1
             include_vert(iv1)=.false.
             cycle
          endif
          do iv2=iv1+1,this%n_vert_end(isize)
             if (.not.include_vert(iv2)) cycle
             if (this%interaction_list(iv1)%n_gs.ne.this%interaction_list(iv2)%n_gs) cycle
             if (this%interaction_list(iv1)%n_ew.ne.this%interaction_list(iv2)%n_ew) cycle
             if (size(this%interaction_list(iv1)%val_c).ne.size(this%interaction_list(iv2)%val_c)) cycle
             if ( sum(abs(this%interaction_list(iv1)%val_c-this%interaction_list(iv2)%val_c))/ &
                  sum(abs(this%interaction_list(iv1)%val_c)+abs(this%interaction_list(iv2)%val_c)).lt.tiny) then
                map_vert(0,1)=map_vert(0,1)+1
                map_vert(map_vert(0,1),1)=iv2
                map_vert(map_vert(0,1),2)=iv1
                identical_vert(isize)=identical_vert(isize)+1
                include_vert(iv2)=.false.
             endif
          enddo
       enddo
    enddo
    do i=1,map_vert(0,1)
       do ic1=1,this%n_cur
          do iv1=1,this%current_list(ic1)%n_vert
             if (this%current_list(ic1)%vertices(iv1).eq.map_vert(i,1)) then
                this%current_list(ic1)%vertices(iv1)=map_vert(i,2)
             endif
          enddo
       enddo
    enddo

    if (allocated(this%include_amp)) deallocate(this%include_amp)
    allocate(this%include_amp(1:this%n_amps))
    this%include_amp=.true.
    call this%filter_dead_trees(n)
    deallocate(this%include_amp)
    do ic=1,this%n_cur
       if (allocated(this%current_list(ic)%val_c)) deallocate(this%current_list(ic)%val_c)
       if (allocated(this%current_list(ic)%val_r)) deallocate(this%current_list(ic)%val_r)
    enddo
    do iv=1,this%n_vert
       if (allocated(this%interaction_list(iv)%val_c)) deallocate(this%interaction_list(iv)%val_c)
       if (allocated(this%interaction_list(iv)%val_r)) deallocate(this%interaction_list(iv)%val_r)
    enddo
    write (99,*) 'Total number of currents, vertices and amplitudes after optimisation',this%n_cur,this%n_vert,this%n_amps
  end subroutine optimise_evaluation

  subroutine create_library(this,n,hel,igroup,iint,pm,p)
    use particles
    implicit none
    class(amplitude_QCD) :: this
    type(physics_model),intent(in) :: pm
    integer,intent(in) :: n,igroup,iint
    real(kind=8),dimension(0:3,n),intent(in) :: p
    integer,parameter :: iunit=14
    integer,dimension(n),intent(in)::hel
    character(len=512) :: line,tmp,idx
    integer :: ip,ibin,i,isize,ih_in,ifinal,ic,iv,iamp,iproc,itype,j,ii,jj,idau,vkey,isector,iterm
    integer :: max_current_vertices
    integer :: chir1,chir2,chiri
    integer,dimension(0:24,0:8) :: icount
    integer,dimension(:,:),allocatable :: icount_type
    integer,dimension(:,:),allocatable :: curs
    integer,dimension(:),allocatable :: pp
    real(kind=8),dimension(:),allocatable :: m,w
    integer,dimension(this%n_vert,0:24,0:8) :: cur1,cur2,int1,pp1,pp2
    real(kind=8),dimension(2,this%n_vert,0:24,0:8) :: coupl
    max_current_vertices=max(1,maxval(this%current_list(:)%n_vert))
    allocate(icount_type(max_current_vertices,11))
    write(tmp,*) igroup
    write(line,*) iint
    line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
    open(file=line,unit=iunit,form='unformatted',access='stream',status='unknown')
    write(iunit) p
    write(iunit) this%amps
    write(iunit) this%amps_by_order
    close(iunit)
    
    write(tmp,*) igroup
    write(line,*) iint
    line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.f03'
    open(file=line,unit=iunit,status='unknown')
    write(line,*) iint
    write(iunit,'(a)') 'module amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib'
    write(iunit,'(2x,a)') 'use FeynmanRules'
    write(iunit,'(2x,a)') 'implicit none'
    write(iunit,'(2x,a)') 'private'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(2x,a)') 'public :: evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))
    write(iunit,'(2x,a)') 'public :: evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_by_order'
    write(iunit,'(2x,a)') 'contains'
    write(iunit,'(2x,a)') 'subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'(p,amps)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) n
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(tmp,*) this%n_amps
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//'),intent(out) :: amps'
    write(line,*) this%n_sectors
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//','//&
         trim(adjustl(line))//') :: amps_by_order'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(4x,a)') 'call evaluate_amp'//trim(adjustl(tmp))//'_'//&
         trim(adjustl(line))//'_by_order(p,amps_by_order)'
    write(iunit,'(4x,a)') 'amps=sum(amps_by_order,dim=2)'
    write(iunit,'(2x,a)') 'end subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))
    write(iunit,'(a)') ''
    write(iunit,'(2x,a)') 'subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//&
         trim(adjustl(line))//'_by_order(p,amps_by_order)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) n
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(tmp,*) this%n_amps
    write(line,*) this%n_sectors
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//','//&
         trim(adjustl(line))//'),intent(out) :: amps_by_order'
    write(tmp,*) this%max_pp
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//') :: pp'
    write(tmp,*) this%n_cur
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//') :: val_c'
    write(tmp,*) this%n_vert
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//') :: int_c'
    if (this%sectors_pruned_empty) then
       ! A locally valid subprocess/order can contain none of the globally
       ! selected sectors.  Its generated routine must be an exact-zero fast
       ! path too: do not build even external currents for it.
       write(iunit,'(4x,a)') 'amps_by_order=(0d0,0d0)'
       write(iunit,'(4x,a)') 'return'
    else
       write(iunit,'(4x,a)') 'call fill_momentum_array(p,pp)'
       write(iunit,'(4x,a)') 'call compute_external_currents(pp,val_c)'
       do isize=2,n-1
          write(tmp,*) isize
          write(iunit,'(4x,a)') 'call compute_vertices'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
          write(iunit,'(4x,a)') 'call compute_currents'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
       enddo
       write(iunit,'(4x,a)') 'call compute_amps_by_order(amps_by_order,val_c)'
    endif
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(2x,a)') 'end subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_by_order'
    write(iunit,'(a)') ''
    write(iunit,'(2x,a)') 'subroutine fill_momentum_array(p,pp)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) n
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(tmp,*) this%max_pp
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(out) :: pp'
    ! fill_momentum_array
    do ip=1,this%max_pp
       ibin=this%pp_i_to_bin(ip)
       write(tmp,*) ip
       line='pp(0:3,'//trim(adjustl(tmp))//')='
       do i=1,n
          write(tmp,*) i
          if (btest(ibin,i-1) .and. i.le.2) &
               line=trim(adjustl(line))//'-p(0:3,'//trim(adjustl(tmp))//')'
          if (btest(ibin,i-1) .and. i.ge.3) &
               line=trim(adjustl(line))//'+p(0:3,'//trim(adjustl(tmp))//')'
       enddo
       write (iunit,'(4x,a)') trim(adjustl(line))
    enddo
    write(iunit,'(2x,a)') 'end subroutine fill_momentum_array'
    write(iunit,'(a)') ''
    do isize=1,n-1
       if (isize.eq.1) then
          write(iunit,'(2x,a)') 'subroutine compute_external_currents(pp,val_c)'
          write(iunit,'(4x,a)') 'implicit none'
          write(tmp,*) this%max_pp
          write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: pp'
          write(tmp,*) this%n_cur
          write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(out) :: val_c'
          ! external wave-functions
          do ic=this%n_cur_start(isize),this%n_cur_end(isize) 
             ifinal=1
             if (this%current_list(ic)%spin(1).eq.-9) then
                ih_in=hel(this%current_list(ic)%order(1))
             else
                ih_in=this%current_list(ic)%spin(1)
             endif
             if (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. pm%is_photon(this%current_list(ic)%type)) then
                write(tmp,*) this%pp_bin_to_i(this%current_list(ic)%bin)
                line='call ext_massless_vector_cmplx(pp(0,'//trim(adjustl(tmp))//'),'
                write(tmp,*) ih_in
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ifinal
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ic
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//'))'
             elseif (pm%is_quark(this%current_list(ic)%type).or. &
                  pm%is_lepton(this%current_list(ic)%type)) then
                write(tmp,*) this%pp_bin_to_i(this%current_list(ic)%bin)
                if (this%current_list(ic)%chirality.ne.0) then
                   line='call ext_fermion_outflow_weyl(pp(0,'//trim(adjustl(tmp))//'),'
                else
                   line='call ext_fermion_outflow(pp(0,'//trim(adjustl(tmp))//'),'
                endif
                write(tmp,*) ih_in
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ifinal
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ic
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//')'
                if (this%current_list(ic)%chirality.ne.0) then
                   write(tmp,*) this%current_list(ic)%chirality
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))//')'
                else
                   line=trim(adjustl(line))//','
                   write(tmp,'(d20.12)') this%current_list(ic)%mass
                   line=trim(adjustl(line))//trim(adjustl(tmp))
                   line=trim(adjustl(line))//')'
                endif
             elseif (pm%is_antiquark(this%current_list(ic)%type).or. &
                  pm%is_antilepton(this%current_list(ic)%type)) then
                write(tmp,*) this%pp_bin_to_i(this%current_list(ic)%bin)
                if (this%current_list(ic)%chirality.ne.0) then
                   line='call ext_fermion_inflow_weyl(pp(0,'//trim(adjustl(tmp))//'),'
                else
                   line='call ext_fermion_inflow(pp(0,'//trim(adjustl(tmp))//'),'
                endif
                write(tmp,*) ih_in
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ifinal
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ic
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//')'
                if (this%current_list(ic)%chirality.ne.0) then
                   write(tmp,*) this%current_list(ic)%chirality
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))//')'
                else
                   line=trim(adjustl(line))//','
                   write(tmp,'(d20.12)') this%current_list(ic)%mass
                   line=trim(adjustl(line))//trim(adjustl(tmp))
                   line=trim(adjustl(line))//')'
                endif
             elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
                write(tmp,*) this%pp_bin_to_i(this%current_list(ic)%bin)
                line='call ext_massive_vector(pp(0,'//trim(adjustl(tmp))//'),'
                write(tmp,*) ih_in
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ifinal
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ic
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//'),'
                write(tmp,'(d20.12)') this%current_list(ic)%mass
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
             elseif (pm%is_higgs(this%current_list(ic)%type)) then
                write(tmp,*) this%pp_bin_to_i(this%current_list(ic)%bin)
                line='call ext_scalar(pp(0,'//trim(adjustl(tmp))//'),'
                write(tmp,*) ifinal
                line=trim(adjustl(line))//trim(adjustl(tmp))//','
                write(tmp,*) ic
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//'))'
             else
                write (*,*) 'External particle type unknown',ic,this%current_list(ic)%type,ih_in
                stop 1
             endif
             write(iunit,'(4x,a)') trim(adjustl(line))
          enddo
          write(iunit,'(2x,a)') 'end subroutine compute_external_currents'
          write(iunit,'(a)') ''
          cycle
       endif

       if (use_real_gluons) then
          write (*,*) 'create library not implemented for use_real_gluons'
          stop 1
       endif
       
       ! interactions
       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       write(tmp,*) isize
       write(iunit,'(2x,a)') 'subroutine compute_vertices'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
       write(iunit,'(4x,a)') 'implicit none'
       write(tmp,*) this%max_pp
       write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: pp'
       write(tmp,*) this%n_cur
       write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: val_c'
       write(tmp,*) this%n_vert
       write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(inout) :: int_c'

       icount(0:24,0:8)=0
       do itype=0,24 ! vertex type
          do iv=this%n_vert_start(isize),this%n_vert_end(isize)
             if (this%interaction_list(iv)%type.eq.itype) then
                chir1=this%current_list(this%interaction_list(iv)%currents(1))%chirality
                chir2=this%current_list(this%interaction_list(iv)%currents(2))%chirality
                chiri=this%interaction_list(iv)%chirality
                vkey=4
                if (itype.eq.4 .or. itype.eq.5 .or. itype.eq.6 .or. itype.eq.7) then
                   vkey=chiri+4
                elseif (itype.eq.10 .or. itype.eq.11 .or. itype.eq.23 .or. itype.eq.24) then
                   if (chiri.ne.0) then
                      ! A compact result uses the Weyl rule and is keyed by
                      ! its result chirality (vkey 3 or 5).
                      vkey=chiri+4
                   elseif (itype.eq.10 .or. itype.eq.11) then
                      ! A full result may still have a compact first child.
                      ! Keep the child chirality in the generated vertex key:
                      ! -1,0,+1 map to vkey 0,4,8.
                      vkey=4*(chir1+1)
                   else
                      ! Types 23 and 24 have the fermion as their second child.
                      vkey=4*(chir2+1)
                   endif
                elseif (itype.eq.8 .or. itype.eq.9 .or. itype.eq.21 .or. itype.eq.22) then
                   vkey=(chir1+1)*3+(chir2+1)
                endif
                icount(itype,vkey)=icount(itype,vkey)+1
                cur1(icount(itype,vkey),itype,vkey)=this%interaction_list(iv)%currents(1)
                cur2(icount(itype,vkey),itype,vkey)=this%interaction_list(iv)%currents(2)
                pp1(icount(itype,vkey),itype,vkey)=this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin)
                pp2(icount(itype,vkey),itype,vkey)=this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin)
                int1(icount(itype,vkey),itype,vkey)=iv
                coupl(1:2,icount(itype,vkey),itype,vkey)=this%interaction_list(iv)%coupl(1:2)
             endif
          enddo
       enddo
       do itype=0,24
          do vkey=0,8
             if (icount(itype,vkey).eq.0) cycle
             write(tmp,*) isize
             line='call vertex_type'//trim(adjustl(tmp))//'_'
             write(tmp,*) itype
             line=trim(adjustl(line))//trim(adjustl(tmp))
             if (vkey.ne.4) then
                write(tmp,*) vkey
                line=trim(adjustl(line))//'_v'//trim(adjustl(tmp))
             endif
             line=trim(adjustl(line))//'(pp,val_c,int_c)'
             write(iunit,'(4x,a)') trim(adjustl(line))
          enddo
       enddo
       write(tmp,*) isize
       write(iunit,'(2x,a)') 'end subroutine compute_vertices'//trim(adjustl(tmp))
       write(iunit,'(a)') ''

       do itype=0,24
          do vkey=0,8
          if (icount(itype,vkey).eq.0) cycle
          write(tmp,*) isize
          line='subroutine vertex_type'//trim(adjustl(tmp))//'_'
          write(tmp,*) itype
          line=trim(adjustl(line))//trim(adjustl(tmp))
          if (vkey.ne.4) then
             write(tmp,*) vkey
             line=trim(adjustl(line))//'_v'//trim(adjustl(tmp))
          endif
          line=trim(adjustl(line))//'(pp,val_c,int_c)'
          write(iunit,'(2x,a)') trim(adjustl(line))
          write(iunit,'(4x,a)') 'implicit none'
          write(tmp,*) this%max_pp
          write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: pp'
          write(tmp,*) this%n_cur
          write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: val_c'
          write(tmp,*) this%n_vert
          write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(inout) :: int_c'
          write(iunit,'(4x,a)') 'integer :: i'
          write(tmp,*) icount(itype,vkey)
          write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: cur1=[ &'
          line=''
          do i=1,icount(itype,vkey)
             write(tmp,*) cur1(i,itype,vkey)
             if (i.eq.1) then
                line=trim(adjustl(tmp))
             else
                line=trim(adjustl(line))//','//trim(adjustl(tmp))
             endif
             if (mod(i,12).eq.0 .and. i.ne.icount(itype,vkey)) then
                line=trim(adjustl(line))//' &'
                write(iunit,'(6x,a)') trim(adjustl(line))
                line=''
             endif
          enddo
          write(iunit,'(6x,a)') trim(adjustl(line))//']'
          write(tmp,*) icount(itype,vkey)
          write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: cur2=[ &'
          line=''
          do i=1,icount(itype,vkey)
             write(tmp,*) cur2(i,itype,vkey)
             if (i.eq.1) then
                line=trim(adjustl(tmp))
             else
                line=trim(adjustl(line))//','//trim(adjustl(tmp))
             endif
             if (mod(i,12).eq.0 .and. i.ne.icount(itype,vkey)) then
                line=trim(adjustl(line))//' &'
                write(iunit,'(6x,a)') trim(adjustl(line))
                line=''
             endif
          enddo
          write(iunit,'(6x,a)') trim(adjustl(line))//']'
          write(tmp,*) icount(itype,vkey)
          write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: int1=[ &'
          line=''
          do i=1,icount(itype,vkey)
             write(tmp,*) int1(i,itype,vkey)
             if (i.eq.1) then
                line=trim(adjustl(tmp))
             else
                line=trim(adjustl(line))//','//trim(adjustl(tmp))
             endif
             if (mod(i,12).eq.0 .and. i.ne.icount(itype,vkey)) then
                line=trim(adjustl(line))//' &'
                write(iunit,'(6x,a)') trim(adjustl(line))
                line=''
             endif
          enddo
          write(iunit,'(6x,a)') trim(adjustl(line))//']'
          if (itype.eq.0 .or. itype.eq.12) then
             write(tmp,*) icount(itype,vkey)
             write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: pp1=[ &'
             line=''
             do i=1,icount(itype,vkey)
                write(tmp,*) pp1(i,itype,vkey)
                if (i.eq.1) then
                   line=trim(adjustl(tmp))
                else
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))
                endif
                if (mod(i,12).eq.0 .and. i.ne.icount(itype,vkey)) then
                   line=trim(adjustl(line))//' &'
                   write(iunit,'(6x,a)') trim(adjustl(line))
                   line=''
                endif
             enddo
             write(iunit,'(6x,a)') trim(adjustl(line))//']'
             write(tmp,*) icount(itype,vkey)
             write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: pp2=[ &'
             line=''
             do i=1,icount(itype,vkey)
                write(tmp,*) pp2(i,itype,vkey)
                if (i.eq.1) then
                   line=trim(adjustl(tmp))
                else
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))
                endif
                if (mod(i,12).eq.0 .and. i.ne.icount(itype,vkey)) then
                   line=trim(adjustl(line))//' &'
                   write(iunit,'(6x,a)') trim(adjustl(line))
                   line=''
                endif
             enddo
             write(iunit,'(6x,a)') trim(adjustl(line))//']'
          endif
          if (itype.eq.8 .or. itype.ge.10) then
             write(tmp,*) icount(itype,vkey)*2
             write(iunit,'(4x,a)') 'real(kind=8),parameter,dimension('//trim(adjustl(tmp))//') :: coupl=[ &'
             line=''
             do i=1,icount(itype,vkey)
                write(tmp,'(D24.16)') coupl(1,i,itype,vkey)
                if (i.eq.1) then
                   line=trim(adjustl(tmp))
                   write(tmp,'(D24.16)') coupl(2,i,itype,vkey)
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))
                else
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))
                   write(tmp,'(D24.16)') coupl(2,i,itype,vkey)
                   line=trim(adjustl(line))//','//trim(adjustl(tmp))
                endif
                if (mod(i,2).eq.0 .and. i.ne.icount(itype,vkey)) then
                   line=trim(adjustl(line))//' &'
                   write(iunit,'(6x,a)') trim(adjustl(line))
                   line=''
                endif
             enddo
             write(iunit,'(6x,a)') trim(adjustl(line))//']'
          endif

       
          write(tmp,*) icount(itype,vkey)
          write(iunit,'(4x,a)')'do i=1,'//trim(adjustl(tmp))
          if (itype.eq.0) then
             line='call threeGluon(val_c(1,cur1(i)),pp(0,pp1(i)),val_c(1,cur2(i)),pp(0,pp2(i)),int_c(1,int1(i)))'
          elseif(itype.eq.1) then
             line='call TwoGluonToAuxTensor(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
          elseif(itype.eq.2) then
             line='call AuxTensorGluonToGluon(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
          elseif(itype.eq.3) then
             line='call GluonAuxTensorToGluon(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
          elseif(itype.eq.4) then
             if (vkey.eq.4) then
                line='call ColourFlowVectorQuarkToQuark(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
             else
                write(tmp,*) vkey-4
                line='call ColourFlowVectorQuarkToQuark_weyl('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'// &
                     trim(adjustl(tmp))//')'
             endif
          elseif(itype.eq.5) then
             if (vkey.eq.4) then
                line='call ColourFlowVectorAntiquarkToAntiquark('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
             else
                write(tmp,*) vkey-4
                line='call ColourFlowVectorAntiquarkToAntiquark_weyl('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'// &
                     trim(adjustl(tmp))//')'
             endif
          elseif(itype.eq.6) then
             if (vkey.eq.4) then
                line='call QuarkColourFlowVectorToQuark(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
             else
                write(tmp,*) vkey-4
                line='call QuarkColourFlowVectorToQuark_weyl('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'// &
                     trim(adjustl(tmp))//')'
             endif
          elseif(itype.eq.7) then
             if (vkey.eq.4) then
                line='call AntiquarkColourFlowVectorToAntiquark('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
             else
                write(tmp,*) vkey-4
                line='call AntiquarkColourFlowVectorToAntiquark_weyl('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'// &
                     trim(adjustl(tmp))//')'
             endif
          elseif(itype.eq.8) then
             if (vkey.eq.4) then
                write(iunit,'(6x,a)') 'call QuarkAntiquarkToColourFlowU1Vector('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             else
                write(iunit,'(6x,a)') 'call QuarkAntiquarkToColourFlowU1Vector_weyl('// &
                     'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/3-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//','
                write(tmp,*) mod(vkey,3)-1
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             endif
             line=''
          elseif(itype.eq.9) then
             if (vkey.eq.4) then
                line='call AntiquarkQuarkToGluon(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)))'
             else
                write(iunit,'(6x,a)') 'call AntiquarkQuarkToGluon_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/3-1
                line=trim(adjustl(tmp))//','
                write(tmp,*) mod(vkey,3)-1
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
                line=''
             endif
          elseif(itype.eq.10) then
             if (vkey.eq.3 .or. vkey.eq.5) then
                write(iunit,'(6x,a)') 'call FermionVectorToFermion_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey-4
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             elseif (vkey.eq.0 .or. vkey.eq.8) then
                write(iunit,'(6x,a)') 'call FermionVectorToFermion_mixed(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/4-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             else
                write(iunit,'(6x,a)') 'call FermionVectorToFermion(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             endif
             line=''
          elseif(itype.eq.11) then
             if (vkey.eq.3 .or. vkey.eq.5) then
                write(iunit,'(6x,a)') 'call AntifermionVectorToAntifermion_weyl( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey-4
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             elseif (vkey.eq.0 .or. vkey.eq.8) then
                write(iunit,'(6x,a)') 'call AntifermionVectorToAntifermion_mixed( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/4-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             else
                write(iunit,'(6x,a)') 'call AntifermionVectorToAntifermion( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             endif
             line=''
          elseif(itype.eq.12) then
             write(iunit,'(6x,a)') 'call VectorVectorToVector( &'
             write(iunit,'(8x,a)') 'val_c(1,cur1(i)),pp(0,pp1(i)),val_c(1,cur2(i)),pp(0,pp2(i)), &'
             write(iunit,'(8x,a)') 'int_c(1,int1(i)),[coupl(2*i-1),coupl(2*i)])'
             line=''
          elseif(itype.eq.13) then
             line='call VectorVectorToAuxTensor(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.14) then
             line='call AuxTensorVectorToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.15) then
             line='call VectorAuxTensorToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.16) then
             line='call FermionScalarToFermion(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.17) then
             line='call VectorVectorToScalar(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.18) then
             line='call ScalarVectorToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.19) then
             line='call VectorScalarToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.20) then
             line='call ScalarScalarToScalar(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)),'//&
                  '[coupl(2*i-1),coupl(2*i)])'
          elseif(itype.eq.21) then
             if (vkey.eq.4) then
                write(iunit,'(6x,a)') 'call FermionAntifermionToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             else
                write(iunit,'(6x,a)') 'call FermionAntifermionToVector_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/3-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//','
                write(tmp,*) mod(vkey,3)-1
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             endif
             line=''
          elseif(itype.eq.22) then
             if (vkey.eq.4) then
                write(iunit,'(6x,a)') 'call AntifermionFermionToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             else
                write(iunit,'(6x,a)') 'call AntifermionFermionToVector_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/3-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//','
                write(tmp,*) mod(vkey,3)-1
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             endif
             line=''
          elseif(itype.eq.23) then
             if (vkey.eq.3 .or. vkey.eq.5) then
                write(iunit,'(6x,a)') 'call VectorFermionToFermion_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey-4
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             elseif (vkey.eq.0 .or. vkey.eq.8) then
                write(iunit,'(6x,a)') 'call VectorFermionToFermion_mixed(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/4-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             else
                write(iunit,'(6x,a)') 'call VectorFermionToFermion(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             endif
             line=''
          elseif(itype.eq.24) then
             if (vkey.eq.3 .or. vkey.eq.5) then
                write(iunit,'(6x,a)') 'call VectorAntifermionToAntifermion_weyl( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey-4
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             elseif (vkey.eq.0 .or. vkey.eq.8) then
                write(iunit,'(6x,a)') 'call VectorAntifermionToAntifermion_mixed( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/4-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             else
                write(iunit,'(6x,a)') 'call VectorAntifermionToAntifermion( &'
                write(iunit,'(8x,a)') 'val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             endif
             line=''
          endif
          if (len_trim(line).gt.0) write(iunit,'(6x,a)')trim(adjustl(line))
          write(iunit,'(4x,a)')'enddo'

          write(tmp,*) isize
          line='end subroutine vertex_type'//trim(adjustl(tmp))//'_'
          write(tmp,*) itype
          line=trim(adjustl(line))//trim(adjustl(tmp))
          if (vkey.ne.4) then
             write(tmp,*) vkey
             line=trim(adjustl(line))//'_v'//trim(adjustl(tmp))
          endif
          write(iunit,'(2x,a)') trim(adjustl(line))
          write(iunit,'(a)') ''

          enddo
       enddo

       write(tmp,*) isize
       write(iunit,'(2x,a)') 'subroutine compute_currents'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
       write(iunit,'(4x,a)') 'implicit none'
       write(tmp,*) this%max_pp
       write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: pp'
       write(tmp,*) this%n_cur
       write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(inout) :: val_c'
       write(tmp,*) this%n_vert
       write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: int_c'

       icount_type=0
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          if (pm%is_colour_flow_vector(this%current_list(ic)%type).or. &
               pm%is_photon(this%current_list(ic)%type)) then
             itype=1
          elseif (pm%is_quark(this%current_list(ic)%type).or. &
               pm%is_lepton(this%current_list(ic)%type)) then
             if (this%current_list(ic)%chirality.ne.0) then
                if (this%current_list(ic)%chirality.eq.1) then
                   itype=8
                else
                   itype=9
                endif
             else
                itype=2
             endif
          elseif (pm%is_antiquark(this%current_list(ic)%type).or. &
               pm%is_antilepton(this%current_list(ic)%type)) then
             if (this%current_list(ic)%chirality.ne.0) then
                if (this%current_list(ic)%chirality.eq.1) then
                   itype=10
                else
                   itype=11
                endif
             else
                itype=3
             endif
          elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
             itype=4
          elseif (pm%is_higgs(this%current_list(ic)%type)) then
             itype=5
          elseif (pm%is_auxiliary_tensor(this%current_list(ic)%type)) then
             itype=6
          elseif (pm%is_auxiliary_scalar(this%current_list(ic)%type)) then
             itype=7
          else
             write (*,*) 'not found:',this%current_list(ic)%type
             stop 1
          endif
          icount_type(this%current_list(ic)%n_vert,itype)=icount_type(this%current_list(ic)%n_vert,itype)+1
       enddo

       do i=1,max_current_vertices
          do j=1,11
             if (icount_type(i,j).eq.0) cycle
             write(tmp,*) isize
             line='call combine_currents_'//trim(adjustl(tmp))
             write(tmp,*) i
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))
             write(tmp,*) j
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
             write(iunit,'(4x,a)') trim(adjustl(line))
          enddo
       enddo
       write(tmp,*) isize
       write(iunit,'(2x,a)') 'end subroutine compute_currents'//trim(adjustl(tmp))
       write(iunit,'(a)') ''

       
       do i=1,max_current_vertices
          do j=1,11
             if (icount_type(i,j).eq.0) cycle

             allocate(curs(0:i,icount_type(i,j)))
             allocate(pp(icount_type(i,j)))
             allocate(m(icount_type(i,j)))
             allocate(w(icount_type(i,j)))
             curs=0
             ii=0
             do ic=this%n_cur_start(isize),this%n_cur_end(isize)
                if (pm%is_colour_flow_vector(this%current_list(ic)%type).or. &
                     pm%is_photon(this%current_list(ic)%type)) then
                   itype=1
                elseif (pm%is_quark(this%current_list(ic)%type).or.&
                     pm%is_lepton(this%current_list(ic)%type)) then
                   if (this%current_list(ic)%chirality.ne.0) then
                      if (this%current_list(ic)%chirality.eq.1) then
                         itype=8
                      else
                         itype=9
                      endif
                   else
                      itype=2
                   endif
                elseif (pm%is_antiquark(this%current_list(ic)%type).or. &
                     pm%is_antilepton(this%current_list(ic)%type)) then
                   if (this%current_list(ic)%chirality.ne.0) then
                      if (this%current_list(ic)%chirality.eq.1) then
                         itype=10
                      else
                         itype=11
                      endif
                   else
                      itype=3
                   endif
                elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
                   itype=4
                elseif (pm%is_higgs(this%current_list(ic)%type)) then
                   itype=5
                elseif (pm%is_auxiliary_tensor(this%current_list(ic)%type)) then
                   itype=6
                elseif (pm%is_auxiliary_scalar(this%current_list(ic)%type)) then
                   itype=7
                else
                   write (*,*) 'not found',this%current_list(ic)%type
                   stop 1
                endif
                if (itype.ne.j) cycle
                if (this%current_list(ic)%n_vert.ne.i) cycle
                ii=ii+1
                ! A negative label tells the generated library to subtract
                ! this interaction when assembling the current.
                curs(1:i,ii)=merge(-this%current_list(ic)%vertices(1:i),&
                     this%current_list(ic)%vertices(1:i),&
                     this%current_list(ic)%vertex_sign(1:i))
                curs(0,ii)=ic
                pp(ii)=this%pp_bin_to_i(this%current_list(ic)%bin)
                m(ii)=this%current_list(ic)%mass
                w(ii)=this%current_list(ic)%width
             enddo
             write(tmp,*) isize
             line='subroutine combine_currents_'//trim(adjustl(tmp))
             write(tmp,*) i
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))
             write(tmp,*) j
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
             write(iunit,'(2x,a)') trim(adjustl(line))
             write(iunit,'(4x,a)') 'implicit none'
             write(tmp,*) this%max_pp
             write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: pp'
             write(tmp,*) this%n_cur
             write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(inout) :: val_c'
             write(tmp,*) this%n_vert
             write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: int_c'
             write(iunit,'(4x,a)') 'integer :: i'
             write(tmp,*) i
             line='integer,parameter,dimension(0:'//trim(adjustl(tmp))//','
             write(tmp,*) icount_type(i,j)
             line=trim(adjustl(line))//trim(adjustl(tmp))//') :: int1=reshape([ &'
             write(iunit,'(4x,a)') trim(adjustl(line))
             line=''
             do ii=1,icount_type(i,j)
                do jj=0,i
                   write(tmp,*) curs(jj,ii)
                   if (ii.eq.1.and.jj.eq.0) then
                      line=trim(adjustl(tmp))
                   else
                      line=trim(adjustl(line))//','//trim(adjustl(tmp))
                   endif
                   if (mod(jj+1+(ii-1)*(i+1),12).eq.0 .and. .not.(ii.eq.icount_type(i,j) .and. jj.eq.i)) then
                      line=trim(adjustl(line))//' &'
                      write(iunit,'(6x,a)') trim(adjustl(line))
                      line=''
                   endif
                enddo
             enddo
             write(tmp,*) i+1
             line=trim(adjustl(line))//'], shape=['//trim(adjustl(tmp))//','
             write(tmp,*) icount_type(i,j)
             line=trim(adjustl(line))//trim(adjustl(tmp))//'])'
             write(iunit,'(6x,a)') trim(adjustl(line))

             if (j.ne.6 .and. j.ne.7 .and. isize.ne.n-1) then
                write(tmp,*) icount_type(i,j)
                write(iunit,'(4x,a)') 'integer,parameter,dimension('//trim(adjustl(tmp))//') :: pp1=[ &'
                line=''
                do ii=1,icount_type(i,j)
                   write(tmp,*) pp(ii)
                   if (ii.eq.1) then
                      line=trim(adjustl(tmp))
                   else
                      line=trim(adjustl(line))//','//trim(adjustl(tmp))
                   endif
                   if (mod(ii,12).eq.0 .and. ii.ne.icount_type(i,j)) then
                      line=trim(adjustl(line))//' &'
                      write(iunit,'(6x,a)') trim(adjustl(line))
                      line=''
                   endif
                enddo
                write(iunit,'(6x,a)') trim(adjustl(line))//']'
             endif

             if (isize.ne.n-1 .and. (j.ge.2 .and. j.le.5)) then
                write(tmp,*) icount_type(i,j)
                write(iunit,'(4x,a)') 'real(kind=8),parameter,dimension('//trim(adjustl(tmp))//') :: m=[ &'
                line=''
                do ii=1,icount_type(i,j)
                   write(tmp,'(D24.16)') m(ii)
                   if (ii.eq.1) then
                      line=trim(adjustl(tmp))
                   else
                      line=trim(adjustl(line))//','//trim(adjustl(tmp))
                   endif
                   if (mod(ii,5).eq.0 .and. ii.ne.icount_type(i,j)) then
                      line=trim(adjustl(line))//' &'
                      write(iunit,'(6x,a)') trim(adjustl(line))
                      line=''
                   endif
                enddo
                write(iunit,'(6x,a)') trim(adjustl(line))//']'
                write(tmp,*) icount_type(i,j)
                write(iunit,'(4x,a)') 'real(kind=8),parameter,dimension('//trim(adjustl(tmp))//') :: w=[ &'
                line=''
                do ii=1,icount_type(i,j)
                   write(tmp,'(D24.16)') w(ii)
                   if (ii.eq.1) then
                      line=trim(adjustl(tmp))
                   else
                      line=trim(adjustl(line))//','//trim(adjustl(tmp))
                   endif
                   if (mod(ii,5).eq.0 .and. ii.ne.icount_type(i,j)) then
                      line=trim(adjustl(line))//' &'
                      write(iunit,'(6x,a)') trim(adjustl(line))
                      line=''
                   endif
                enddo
                write(iunit,'(6x,a)') trim(adjustl(line))//']'
             endif

             write(tmp,*) icount_type(i,j)
             write(iunit,'(4x,a)') 'do i=1,'//trim(adjustl(tmp))
             write(tmp,*) i
             if (j.eq.6) then
                write(iunit,'(6x,a)') 'val_c(1:6,int1(0,i))=sum(int_c(1:6,abs(int1(1:'//trim(adjustl(tmp))//&
                     ',i)))*spread(sign(1,int1(1:'//trim(adjustl(tmp))//',i)),1,6),dim=2)'
             elseif (j.eq.5 .or. j.eq.7) then
                write(iunit,'(6x,a)') 'val_c(1,int1(0,i))=sum(int_c(1,abs(int1(1:'//trim(adjustl(tmp))//&
                     ',i)))*sign(1,int1(1:'//trim(adjustl(tmp))//',i)))'
             elseif (j.ge.8 .and. j.le.11) then
                write(iunit,'(6x,a)') 'val_c(1:2,int1(0,i))=sum(int_c(1:2,abs(int1(1:'//trim(adjustl(tmp))//&
                     ',i)))*spread(sign(1,int1(1:'//trim(adjustl(tmp))//',i)),1,2),dim=2)'
             else
                write(iunit,'(6x,a)') 'val_c(1:4,int1(0,i))=sum(int_c(1:4,abs(int1(1:'//trim(adjustl(tmp))//&
                     ',i)))*spread(sign(1,int1(1:'//trim(adjustl(tmp))//',i)),1,4),dim=2)'
             endif
             if (j.eq.1 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call MasslessVectorPropagator(val_c(1,int1(0,i)),pp(0,pp1(i)))'
             elseif(j.eq.2 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call FermionPropagator(val_c(1,int1(0,i)),pp(0,pp1(i)),m(i),w(i))'
             elseif(j.eq.3 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call AntifermionPropagator(val_c(1,int1(0,i)),pp(0,pp1(i)),m(i),w(i))'
             elseif(j.eq.4 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call MassiveVectorPropagator(val_c(1,int1(0,i)),pp(0,pp1(i)),m(i),w(i))'
             elseif(j.eq.5 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call ScalarPropagator(val_c(1,int1(0,i)),pp(0,pp1(i)),m(i),w(i))'
             elseif(j.eq.8 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call FermionPropagator_weyl(val_c(1,int1(0,i)),pp(0,pp1(i)),1)'
             elseif(j.eq.9 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call FermionPropagator_weyl(val_c(1,int1(0,i)),pp(0,pp1(i)),-1)'
             elseif(j.eq.10 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call AntifermionPropagator_weyl(val_c(1,int1(0,i)),pp(0,pp1(i)),1)'
             elseif(j.eq.11 .and. isize.ne.n-1) then
                write(iunit,'(6x,a)') 'call AntifermionPropagator_weyl(val_c(1,int1(0,i)),pp(0,pp1(i)),-1)'
             endif
             write(iunit,'(4x,a)') 'enddo'
             write(tmp,*) isize
             line='end subroutine combine_currents_'//trim(adjustl(tmp))
             write(tmp,*) i
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))
             write(tmp,*) j
             line=trim(adjustl(line))//'_'//trim(adjustl(tmp))
             write(iunit,'(2x,a)') trim(adjustl(line))
             write(iunit,'(a)') ''
             deallocate(curs)
             deallocate(pp)
             deallocate(m)
             deallocate(w)
          enddo
       enddo
    enddo

    write(iunit,'(2x,a)') 'subroutine compute_amps_by_order(amps_by_order,val_c)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) this%n_amps
    write(line,*) this%n_sectors
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//','//&
         trim(adjustl(line))//'),intent(out) :: amps_by_order'
    write(tmp,*) this%n_cur
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: val_c'
    write(iunit,'(4x,a)') 'amps_by_order=(0d0,0d0)'
    if (this%sectors_pruned_empty) write(iunit,'(4x,a)') 'return'
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          if (this%same_flav(iproc)) cycle
          do isector=1,this%n_sectors
             if (.not.this%sector_present(iamp,isector)) cycle
             do iterm=this%sector_term_start(iamp-1,isector)+1,&
                  this%sector_term_start(iamp,isector)
                write(tmp,*) iamp
                write(idx,*) isector
                line='amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')='//&
                     'amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')'
                if (this%sector_term_sign(iterm).eq.1) then
                   line=trim(adjustl(line))//'+'
                elseif (this%sector_term_sign(iterm).eq.-1) then
                   line=trim(adjustl(line))//'-'
                else
                   write(tmp,*) this%sector_term_sign(iterm)
                   line=trim(adjustl(line))//'+'//trim(adjustl(tmp))//'.0d0*'
                endif
                line=trim(adjustl(line))//'ContractFermionCurrents('
                write(tmp,*) this%sector_term_curr2amp(1,iterm)
                line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//'),'
                write(tmp,*) this%current_list(this%sector_term_curr2amp(1,iterm))%chirality
                line=trim(adjustl(line))//trim(adjustl(tmp))//',val_c(1,'
                write(tmp,*) this%sector_term_curr2amp(2,iterm)
                line=trim(adjustl(line))//trim(adjustl(tmp))//'),'
                write(tmp,*) this%current_list(this%sector_term_curr2amp(2,iterm))%chirality
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(4x,a)') trim(adjustl(line))
             enddo
             if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. &
                  iamp.gt.this%n_amps/2 .and. mod(n,2).eq.1) then
                write(tmp,*) iamp
                write(idx,*) isector
                write(iunit,'(4x,a)') 'amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')=-'//&
                     'amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')'
             endif
          enddo
       enddo
    enddo
    do iproc=1,this%nprocs
       if (.not.this%same_flav(iproc)) cycle
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          do isector=1,this%n_sectors
             if (this%sectors_pruned) then
                if (.not.this%sector_retained(iamp,isector)) cycle
             endif
             write(tmp,*) iamp
             write(idx,*) isector
             line='amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')='
             do idau=1,2
                if (this%same_flavour_sum(iamp,idau).le.0) cycle
                write(tmp,*) this%same_flavour_sum(iamp,idau)
                if (this%same_flavour_sum_operation(iamp,idau).eq.0) then
                   line=trim(adjustl(line))//'+amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.1) then
                   line=trim(adjustl(line))//'-conjg(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//'))'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.2) then
                   line=trim(adjustl(line))//'+conjg(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//'))'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.3) then
                   line=trim(adjustl(line))//'-amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.4) then
                   line=trim(adjustl(line))//'+cmplx(aimag(amps_by_order('//trim(adjustl(tmp))//','//&
                        trim(adjustl(idx))//')),dble(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')))'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.5) then
                   line=trim(adjustl(line))//'+cmplx(-aimag(amps_by_order('//trim(adjustl(tmp))//','//&
                        trim(adjustl(idx))//')),dble(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')))'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.6) then
                   line=trim(adjustl(line))//'+cmplx(aimag(amps_by_order('//trim(adjustl(tmp))//','//&
                        trim(adjustl(idx))//')),-dble(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')))'
                elseif (this%same_flavour_sum_operation(iamp,idau).eq.7) then
                   line=trim(adjustl(line))//'+cmplx(-aimag(amps_by_order('//trim(adjustl(tmp))//','//&
                        trim(adjustl(idx))//')),-dble(amps_by_order('//trim(adjustl(tmp))//','//trim(adjustl(idx))//')))'
                else
                   write (*,*) 'ERROR: unknown sector operation in creating library',&
                        this%same_flavour_sum_operation(iamp,idau)
                   stop 1
                endif
             enddo
             write(iunit,'(4x,a)') trim(adjustl(line))
          enddo
       enddo
    enddo
    write(iunit,'(2x,a)') 'end subroutine compute_amps_by_order'
    write(iunit,'(a)') ''

    write(iunit,'(2x,a)') 'subroutine compute_amps(amps,val_c)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) this%n_amps
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//'),intent(out) :: amps'
    write(tmp,*) this%n_cur
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: val_c'
    write(tmp,*) this%n_amps
    write(line,*) this%n_sectors
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//','//&
         trim(adjustl(line))//') :: amps_by_order'
    write(iunit,'(4x,a)') 'call compute_amps_by_order(amps_by_order,val_c)'
    write(iunit,'(4x,a)') 'amps=sum(amps_by_order,dim=2)'
    write(iunit,'(2x,a)') 'end subroutine compute_amps'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(a)') 'end module amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib'
    close(iunit)
    deallocate(icount_type)
  end subroutine create_library



  subroutine filter_helicity(this,n,nhel,include_hel)
    implicit none
    class(amplitude_qcd) :: this
    integer,intent(in) :: n
    integer,intent(inout) :: nhel
    integer,intent(inout),dimension(nhel) :: include_hel
    integer :: nspin,ispin,ic,iv,iamp,isector,iterm
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
    if (allocated(this%amps_by_order)) deallocate(this%amps_by_order)
    if (allocated(this%amps_r)) deallocate(this%amps_r)
    
    allocate(include_current(this%n_cur))
    include_current=.false.
    include_current(this%n_cur_start(n  ):this%n_cur_end(n  ))=.false.
    include_current(this%n_cur_start(n-1):this%n_cur_end(n-1))=.false.

    this%include_amp(1:this%n_amps)=.false.
    
    allocate(tmp_spin(1:n,1:maxval(include_hel),nhel))

    nspin=0
    do iamp=1,nhel
       if (include_hel(iamp).ge.1) then
          this%include_amp(iamp)=.true.
          if (this%same_flavour_sum(iamp,1).le.0) then
             do isector=1,this%n_sectors
                if (.not.this%sector_present(iamp,isector)) cycle
                if (allocated(this%sector_term_start)) then
                   do iterm=this%sector_term_start(iamp-1,isector)+1,&
                        this%sector_term_start(iamp,isector)
                      include_current(this%sector_term_curr2amp(1,iterm))=.true.
                      include_current(this%sector_term_curr2amp(2,iterm))=.true.
                   enddo
                else
                   include_current(this%sector_curr2amp(1,iamp,isector))=.true.
                   include_current(this%sector_curr2amp(2,iamp,isector))=.true.
                endif
                if (allocated(this%sector_three_line_partner_curr2amp)) then
                   if (this%sector_three_line_partner_curr2amp(1,iamp,isector).ne.0) then
                      include_current(this%sector_three_line_partner_curr2amp(1,iamp,isector))=.true.
                      include_current(this%sector_three_line_partner_curr2amp(2,iamp,isector))=.true.
                   endif
                endif
             enddo
          else
             ! same-flavour amplitude.
             if (all(include_hel(this%same_flavour_sum(iamp,1:2)).eq.0)) then
                write (*,*) 'inconsistency in helicity filter for same-flavour process #1'
                write (*,*) iamp,this%same_flavour_sum(iamp,1:2)
                stop 1
             endif
             if (include_hel(this%same_flavour_sum(iamp,1)).lt.0) then
                this%same_flavour_sum(iamp,1)=-include_hel(this%same_flavour_sum(iamp,1))
                if (.not. this%include_amp(this%same_flavour_sum(iamp,1))) then
                   write (*,*) 'inconsistency in helicity filter for same-flavour process #2'
                   stop 1
                endif
             endif
             if (include_hel(this%same_flavour_sum(iamp,2)).lt.0) then
                this%same_flavour_sum(iamp,2)=-include_hel(this%same_flavour_sum(iamp,2))
                if (.not. this%include_amp(this%same_flavour_sum(iamp,2))) then
                   write (*,*) 'inconsistency in helicity filter for same-flavour process #3'
                   stop 1
                endif
             endif
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
    write (99,*) 'Total number of currents, vertices and amplitudes after filtering helicities',this%n_cur,this%n_vert,this%n_amps
    deallocate(this%include_amp)

  end subroutine filter_helicity

  subroutine prune_coupling_sectors(this,allowed_sector_pairs)
    ! Remove amplitude-level coupling sectors which cannot contribute to an
    ! already-resolved squared-order selection.  The selector is deliberately
    ! passed as a sector-pair mask so this core module does not depend on the
    ! input-language/coupling_orders module.
    !
    ! For fixed-colour imode=1 a sector is retained separately for every
    ! physical helicity amplitude iff it has a partner in that same slot.  For
    ! full-colour imode=2, the colour Gram matrix can couple different flow
    ! slots, so use the conservative global endpoint union: every present
    ! (flow,sector) is retained if its sector has an allowed partner in any
    ! other present flow.  The mask need not be symmetric; either orientation
    ! is accepted.  Interference-only endpoints and shared currents therefore
    ! remain live.  The public sector and physical-amplitude axes remain
    ! stable; pruned amp-sector entries have sector_present=.false. and empty
    ! sparse-terminal ranges.
    implicit none
    class(amplitude_qcd),intent(inout) :: this
    logical,dimension(:,:),intent(in) :: allowed_sector_pairs
    logical,dimension(:,:),allocatable :: retain_sector,effective_present
    logical,dimension(:),allocatable :: include_current,global_present
    integer,dimension(:,:),allocatable :: new_term_start,new_term_curr2amp
    integer,dimension(:),allocatable :: new_term_sign
    integer :: n,old_nterms,new_nterms,iamp,isector,jsector,iterm,new_term,&
         iproc,idau,daughter
    logical :: changed,had_include_amp

    this%sectors_pruned_empty=.false.

    if (size(allowed_sector_pairs,1).ne.this%n_sectors .or. &
         size(allowed_sector_pairs,2).ne.this%n_sectors) then
       write (*,*) 'Coupling-sector pair mask has incompatible shape',&
            shape(allowed_sector_pairs),this%n_sectors
       stop 1
    endif
    if (this%imode.ne.1 .and. this%imode.ne.2) then
       write (*,*) 'Coupling-sector pruning only supports imode=1 or imode=2 amplitudes',&
            this%imode
       stop 1
    endif
    if (this%imode.eq.2 .and. .not.allocated(this%col_index)) then
       write (*,*) 'Full-colour coupling-sector pruning must be applied after init_col'
       stop 1
    endif
    if (this%sectors_pruned) then
       write (*,*) 'Coupling-sector pruning may only be applied once to an amplitude'
       stop 1
    endif
    if (.not.allocated(this%sector_present) .or. &
         .not.allocated(this%sector_term_start) .or. &
         .not.allocated(this%sector_term_curr2amp) .or. &
         .not.allocated(this%sector_term_sign)) then
       write (*,*) 'Coupling-sector pruning requires initialized sparse sector metadata'
       stop 1
    endif
    if (this%nprocs.lt.1 .or. size(this%processes,1).lt.2) then
       write (*,*) 'Coupling-sector pruning requires an initialized amplitude'
       stop 1
    endif
    n=size(this%processes,1)

    allocate(retain_sector(this%n_amps,this%n_sectors))
    allocate(effective_present(this%n_amps,this%n_sectors))
    effective_present=this%sector_present
    ! Same-flavour amplitudes are derived linear combinations and have no
    ! terminal roots of their own.  Resolve their effective sector support
    ! from the daughter graph before applying the pair selector.
    do
       changed=.false.
       do iproc=1,this%nprocs
          if (.not.this%same_flav(iproc)) cycle
          do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
             do idau=1,2
                daughter=this%same_flavour_sum(iamp,idau)
                if (daughter.le.0 .or. daughter.gt.this%n_amps) cycle
                do isector=1,this%n_sectors
                   if (effective_present(iamp,isector) .or. &
                        .not.effective_present(daughter,isector)) cycle
                   effective_present(iamp,isector)=.true.
                   changed=.true.
                enddo
             enddo
          enddo
       enddo
       if (.not.changed) exit
    enddo
    retain_sector=.false.
    if (this%imode.eq.1) then
       do iamp=1,this%n_amps
          do isector=1,this%n_sectors
             if (.not.effective_present(iamp,isector)) cycle
             do jsector=1,this%n_sectors
                if (.not.effective_present(iamp,jsector)) cycle
                if (allowed_sector_pairs(isector,jsector) .or. &
                     allowed_sector_pairs(jsector,isector)) then
                   retain_sector(iamp,isector)=.true.
                   exit
                endif
             enddo
          enddo
       enddo
    else
       allocate(global_present(this%n_sectors))
       do isector=1,this%n_sectors
          global_present(isector)=any(effective_present(:,isector))
       enddo
       do iamp=1,this%n_amps
          do isector=1,this%n_sectors
             if (.not.effective_present(iamp,isector)) cycle
             do jsector=1,this%n_sectors
                if (.not.global_present(jsector)) cycle
                if (allowed_sector_pairs(isector,jsector) .or. &
                     allowed_sector_pairs(jsector,isector)) then
                   retain_sector(iamp,isector)=.true.
                   exit
                endif
             enddo
          enddo
       enddo
       deallocate(global_present)
    endif
    ! A retained sector of a derived amplitude requires the same sector from
    ! each daughter that contributes it.  Close that dependency graph before
    ! compacting terminal roots so same-flavour interference remains exact.
    do
       changed=.false.
       do iproc=1,this%nprocs
          if (.not.this%same_flav(iproc)) cycle
          do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
             do isector=1,this%n_sectors
                if (.not.retain_sector(iamp,isector)) cycle
                do idau=1,2
                   daughter=this%same_flavour_sum(iamp,idau)
                   if (daughter.le.0 .or. daughter.gt.this%n_amps) cycle
                   if (.not.effective_present(daughter,isector)) cycle
                   if (.not.retain_sector(daughter,isector)) then
                      retain_sector(daughter,isector)=.true.
                      changed=.true.
                   endif
                enddo
             enddo
          enddo
       enddo
       if (.not.changed) exit
    enddo

    old_nterms=size(this%sector_term_sign)
    new_nterms=0
    do isector=1,this%n_sectors
       do iamp=1,this%n_amps
          if (.not.retain_sector(iamp,isector)) cycle
          new_nterms=new_nterms+this%sector_term_start(iamp,isector)-&
               this%sector_term_start(iamp-1,isector)
       enddo
    enddo
    allocate(new_term_start(0:this%n_amps,this%n_sectors))
    allocate(new_term_curr2amp(2,new_nterms))
    allocate(new_term_sign(new_nterms))
    new_term=0
    do isector=1,this%n_sectors
       new_term_start(0,isector)=new_term
       do iamp=1,this%n_amps
          if (retain_sector(iamp,isector)) then
             do iterm=this%sector_term_start(iamp-1,isector)+1,&
                  this%sector_term_start(iamp,isector)
                new_term=new_term+1
                new_term_curr2amp(:,new_term)=this%sector_term_curr2amp(:,iterm)
                new_term_sign(new_term)=this%sector_term_sign(iterm)
             enddo
          endif
          new_term_start(iamp,isector)=new_term
       enddo
    enddo
    if (new_term.ne.new_nterms) then
       write (*,*) 'Internal error compacting coupling-sector terminal terms',&
            new_term,new_nterms,old_nterms
       stop 1
    endif
    call move_alloc(new_term_start,this%sector_term_start)
    call move_alloc(new_term_curr2amp,this%sector_term_curr2amp)
    call move_alloc(new_term_sign,this%sector_term_sign)

    this%sector_present=retain_sector
    if (allocated(this%sector_retained)) deallocate(this%sector_retained)
    allocate(this%sector_retained(this%n_amps,this%n_sectors))
    this%sector_retained=retain_sector
    this%sectors_pruned=.true.
    do iamp=1,this%n_amps
       do isector=1,this%n_sectors
          if (retain_sector(iamp,isector)) cycle
          this%sector_curr2amp(:,iamp,isector)=0
          this%sector_sign(iamp,isector)=1
          if (allocated(this%sector_three_line_partner_curr2amp)) &
               this%sector_three_line_partner_curr2amp(:,iamp,isector)=0
       enddo
    enddo
    ! The legacy coherent root is only a compatibility representative.  Point
    ! it at the first retained sparse term, or clear it for a locally empty
    ! physical amplitude.  This does not alter n_amps or any public indexing.
    this%curr2amp=0
    do iamp=1,this%n_amps
       do isector=1,this%n_sectors
          if (.not.retain_sector(iamp,isector)) cycle
          if (this%sector_term_start(iamp,isector).eq.&
               this%sector_term_start(iamp-1,isector)) cycle
          iterm=this%sector_term_start(iamp-1,isector)+1
          this%sector_curr2amp(:,iamp,isector)=&
               this%sector_term_curr2amp(:,iterm)
          this%sector_sign(iamp,isector)=this%sector_term_sign(iterm)
          if (iamp.le.size(this%curr2amp,2)) then
             if (all(this%curr2amp(:,iamp).eq.0)) &
                  this%curr2amp(:,iamp)=this%sector_term_curr2amp(:,iterm)
          endif
       enddo
    enddo

    ! filter_dead_trees normally keeps every terminal current.  Supplying this
    ! explicit mask makes the recursion start only at retained sparse roots,
    ! while shared lower currents/interactions are discovered and preserved by
    ! its backwards traversal.
    allocate(include_current(this%n_cur))
    include_current=.false.
    do iterm=1,size(this%sector_term_sign)
       include_current(this%sector_term_curr2amp(1,iterm))=.true.
       include_current(this%sector_term_curr2amp(2,iterm))=.true.
    enddo
    had_include_amp=allocated(this%include_amp)
    if (had_include_amp) then
       this%include_amp=.true.
    else
       allocate(this%include_amp(this%n_amps))
       this%include_amp=.true.
    endif

    ! Evaluation work arrays refer to the old current/interaction arrays and
    ! dimensions.  Recreate them lazily on the next evaluate call.
    do iterm=1,this%n_cur
       if (allocated(this%current_list(iterm)%val_c)) &
            deallocate(this%current_list(iterm)%val_c)
       if (allocated(this%current_list(iterm)%val_r)) &
            deallocate(this%current_list(iterm)%val_r)
    enddo
    do iterm=1,this%n_vert
       if (allocated(this%interaction_list(iterm)%val_c)) &
            deallocate(this%interaction_list(iterm)%val_c)
       if (allocated(this%interaction_list(iterm)%val_r)) &
            deallocate(this%interaction_list(iterm)%val_r)
    enddo
    if (allocated(this%amps)) deallocate(this%amps)
    if (allocated(this%amps_by_order)) deallocate(this%amps_by_order)
    if (allocated(this%amps_r)) deallocate(this%amps_r)

    if (new_nterms.eq.0) then
       call retain_external_currents_only()
       this%sectors_pruned_empty=.true.
    else
       call this%filter_dead_trees(n,include_current)
    endif
    ! Unlike helicity filtering, squared-order pruning never owns or compacts
    ! the physical amplitude axis.  Keep include_amp only if it already
    ! existed on entry; otherwise restore the ordinary post-init state.
    if (.not.had_include_amp .and. allocated(this%include_amp)) &
         deallocate(this%include_amp)
    deallocate(include_current,retain_sector,effective_present)
  contains
    subroutine retain_external_currents_only()
      implicit none
      type(current),dimension(:),allocatable :: old_currents
      integer :: old_ncur,old_nvert,ic,iv,new_cur

      old_ncur=this%n_cur
      old_nvert=this%n_vert
      allocate(old_currents(old_ncur))
      do ic=1,old_ncur
         old_currents(ic)=this%current_list(ic)
      enddo
      new_cur=0
      do ic=this%n_cur_start(1),this%n_cur_end(1)
         new_cur=new_cur+1
         if (new_cur.ne.ic) this%current_list(new_cur)=old_currents(ic)
         this%current_list(new_cur)%n_vert=0
         if (allocated(this%current_list(new_cur)%vertices)) &
              deallocate(this%current_list(new_cur)%vertices)
         if (allocated(this%current_list(new_cur)%vertex_sign)) &
              deallocate(this%current_list(new_cur)%vertex_sign)
      enddo
      do ic=1,old_ncur
         call finalize_current(old_currents(ic))
      enddo
      deallocate(old_currents)
      do ic=new_cur+1,old_ncur
         call finalize_current(this%current_list(ic))
      enddo
      do iv=1,old_nvert
         call finalize_interaction(this%interaction_list(iv))
      enddo
      this%n_cur=new_cur
      this%n_vert=0
      this%n_cur_start(1)=1
      this%n_cur_end(1)=new_cur
      do ic=2,n
         this%n_cur_start(ic)=new_cur+1
         this%n_cur_end(ic)=new_cur
      enddo
      do ic=2,n-1
         this%n_vert_start(ic)=1
         this%n_vert_end(ic)=0
      enddo
      this%curr2amp=0
      this%sector_curr2amp=0
      this%sector_sign=1
      if (allocated(this%three_line_partner_curr2amp)) &
           this%three_line_partner_curr2amp=0
      if (allocated(this%sector_three_line_partner_curr2amp)) &
           this%sector_three_line_partner_curr2amp=0
    end subroutine retain_external_currents_only
  end subroutine prune_coupling_sectors

  
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
    integer,dimension(:),allocatable :: new_sector_term_sign
    integer,dimension(:,:),allocatable :: old_sector_term_start,&
         new_sector_term_start,new_sector_term_curr2amp
    integer,dimension(:,:),allocatable :: old_sector_term_curr2amp
    integer,dimension(:),allocatable :: old_sector_term_sign
    logical,dimension(*),optional :: include_current
    integer :: to_skip,isize,nc,iv,n,iamp,iproc,i,isector,iterm,&
         old_term,nterms_new,new_amp,old_cur_start,old_cur_end,&
         old_vert_start,old_vert_end
    logical :: found_range
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
       is_needed_cur(:)=include_current(1:this%n_cur)
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
       if (iamp.le.size(this%curr2amp,2) .and. &
            where_to_amp(iamp).le.size(this%curr2amp,2)) then
          do i=1,2
             if (this%curr2amp(i,iamp).ne.0) then
                this%curr2amp(i,where_to_amp(iamp))=where_to_cur(this%curr2amp(i,iamp))
             endif
          enddo
       endif
       do isector=1,this%n_sectors
          this%sector_present(where_to_amp(iamp),isector)=this%sector_present(iamp,isector)
          if (allocated(this%sector_retained)) &
               this%sector_retained(where_to_amp(iamp),isector)=&
               this%sector_retained(iamp,isector)
          this%sector_sign(where_to_amp(iamp),isector)=this%sector_sign(iamp,isector)
          do i=1,2
             if (this%sector_curr2amp(i,iamp,isector).ne.0) then
                this%sector_curr2amp(i,where_to_amp(iamp),isector)=&
                     where_to_cur(this%sector_curr2amp(i,iamp,isector))
             else
                this%sector_curr2amp(i,where_to_amp(iamp),isector)=0
             endif
          enddo
       enddo
    enddo
    if (allocated(this%three_line_partner_curr2amp)) then
       do iamp=1,this%n_amps
          if (.not.this%include_amp(iamp)) cycle
          do i=1,2
             if (this%three_line_partner_curr2amp(i,iamp).ne.0) then
                this%three_line_partner_curr2amp(i,where_to_amp(iamp))=&
                     where_to_cur(this%three_line_partner_curr2amp(i,iamp))
             else
                this%three_line_partner_curr2amp(i,where_to_amp(iamp))=0
             endif
          enddo
       enddo
    endif
    if (allocated(this%sector_three_line_partner_curr2amp)) then
       do iamp=1,this%n_amps
          if (.not.this%include_amp(iamp)) cycle
          do isector=1,this%n_sectors
             do i=1,2
                if (this%sector_three_line_partner_curr2amp(i,iamp,isector).ne.0) then
                   this%sector_three_line_partner_curr2amp(i,where_to_amp(iamp),isector)=&
                        where_to_cur(this%sector_three_line_partner_curr2amp(i,iamp,isector))
                else
                   this%sector_three_line_partner_curr2amp(i,where_to_amp(iamp),isector)=0
                endif
             enddo
          enddo
       enddo
    endif
    if (allocated(this%sector_term_start)) then
       allocate(old_sector_term_start(0:this%n_amps,this%n_sectors))
       allocate(old_sector_term_curr2amp(2,size(this%sector_term_sign)))
       allocate(old_sector_term_sign(size(this%sector_term_sign)))
       old_sector_term_start=this%sector_term_start
       old_sector_term_curr2amp=this%sector_term_curr2amp
       old_sector_term_sign=this%sector_term_sign
       nterms_new=0
       do isector=1,this%n_sectors
          do iamp=1,this%n_amps
             if (.not.this%include_amp(iamp)) cycle
             nterms_new=nterms_new+old_sector_term_start(iamp,isector)-&
                  old_sector_term_start(iamp-1,isector)
          enddo
       enddo
       allocate(new_sector_term_start(0:count(this%include_amp),this%n_sectors))
       allocate(new_sector_term_curr2amp(2,nterms_new))
       allocate(new_sector_term_sign(nterms_new))
       iterm=0
       do isector=1,this%n_sectors
          new_sector_term_start(0,isector)=iterm
          new_amp=0
          do iamp=1,this%n_amps
             if (.not.this%include_amp(iamp)) cycle
             new_amp=new_amp+1
             do old_term=old_sector_term_start(iamp-1,isector)+1,&
                  old_sector_term_start(iamp,isector)
                iterm=iterm+1
                new_sector_term_curr2amp(:,iterm)=where_to_cur(&
                     old_sector_term_curr2amp(:,old_term))
                new_sector_term_sign(iterm)=old_sector_term_sign(old_term)
             enddo
             new_sector_term_start(new_amp,isector)=iterm
          enddo
       enddo
       call move_alloc(new_sector_term_start,this%sector_term_start)
       call move_alloc(new_sector_term_curr2amp,this%sector_term_curr2amp)
       call move_alloc(new_sector_term_sign,this%sector_term_sign)
       deallocate(old_sector_term_start,old_sector_term_curr2amp,&
            old_sector_term_sign)
    endif
    ! Keep legacy representative roots synchronized with the compact sparse
    ! terminal map.  A physical amplitude is intentionally allowed to have no
    ! retained sector after squared-order pruning.
    this%curr2amp=0
    do new_amp=1,count(this%include_amp)
       do isector=1,this%n_sectors
          if (.not.this%sector_present(new_amp,isector)) cycle
          if (this%sector_term_start(new_amp,isector).eq.&
               this%sector_term_start(new_amp-1,isector)) cycle
          old_term=this%sector_term_start(new_amp-1,isector)+1
          this%sector_curr2amp(:,new_amp,isector)=&
               this%sector_term_curr2amp(:,old_term)
          this%sector_sign(new_amp,isector)=this%sector_term_sign(old_term)
          if (new_amp.le.size(this%curr2amp,2)) &
               this%curr2amp(:,new_amp)=this%sector_term_curr2amp(:,old_term)
          exit
       enddo
    enddo
    do iamp=1,this%n_amps
       if (.not.this%include_amp(iamp)) cycle
       do i=1,2
          if (this%same_flavour_sum(iamp,i).eq.0) then
             this%same_flavour_sum(where_to_amp(iamp),i)=0
             this%same_flavour_sum_operation(where_to_amp(iamp),i)=0
          elseif (this%same_flavour_sum(iamp,i).gt.0) then
             this%same_flavour_sum(where_to_amp(iamp),i)=where_to_amp(this%same_flavour_sum(iamp,i))
             this%same_flavour_sum_operation(where_to_amp(iamp),i)=this%same_flavour_sum_operation(iamp,i)
          endif
       enddo
    enddo
    ! and also the shifting of the auxiliary arrays and variables
    do isize=1,n
       old_cur_start=this%n_cur_start(isize)
       old_cur_end=this%n_cur_end(isize)
       found_range=.false.
       do nc=old_cur_start,old_cur_end
          if (where_to_cur(nc).ne.0) then
             this%n_cur_start(isize)=where_to_cur(nc)
             found_range=.true.
             exit
          endif
       enddo
       if (found_range) then
          do nc=old_cur_end,old_cur_start,-1
             if (where_to_cur(nc).ne.0) then
                this%n_cur_end(isize)=where_to_cur(nc)
                exit
             endif
          enddo
       elseif (isize.eq.1) then
          this%n_cur_start(isize)=1
          this%n_cur_end(isize)=0
       elseif (isize.eq.n) then
          ! The size-n range is a sentinel for external closing currents and
          ! is not ordered after the composite-current blocks.
          this%n_cur_start(isize)=count(is_needed_cur)+1
          this%n_cur_end(isize)=count(is_needed_cur)
       else
          this%n_cur_start(isize)=this%n_cur_end(isize-1)+1
          this%n_cur_end(isize)=this%n_cur_end(isize-1)
       endif
       if (isize.ge.2 .and. isize.le.n-1) then
          old_vert_start=this%n_vert_start(isize)
          old_vert_end=this%n_vert_end(isize)
          found_range=.false.
          do iv=old_vert_start,old_vert_end
             if (where_to_ver(iv).ne.0) then
                this%n_vert_start(isize)=where_to_ver(iv)
                found_range=.true.
                exit
             endif
          enddo
          if (found_range) then
             do iv=old_vert_end,old_vert_start,-1
                if (where_to_ver(iv).ne.0) then
                   this%n_vert_end(isize)=where_to_ver(iv)
                   exit
                endif
             enddo
          elseif (isize.eq.2) then
             this%n_vert_start(isize)=1
             this%n_vert_end(isize)=0
          else
             this%n_vert_start(isize)=this%n_vert_end(isize-1)+1
             this%n_vert_end(isize)=this%n_vert_end(isize-1)
          endif
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
    do iamp=1,this%n_amps
       do i=1,2
          if (this%same_flavour_sum(iamp,i).gt.this%n_amps) then
             this%same_flavour_sum(iamp,i)=0
          endif
       enddo
    enddo
    do iproc=this%nprocs+1,2,-1
       if (this%iproc_start(iproc).gt.this%n_amps+1) then
          this%iproc_start(iproc)=this%n_amps+1
       endif
    enddo
    deallocate(is_needed_ver)
    deallocate(is_needed_cur)
    deallocate(where_to_ver)
    deallocate(where_to_cur)
    deallocate(where_to_amp)
  end subroutine filter_dead_trees

  integer function sector_index(this,n_gs,n_ew)
    implicit none
    class(amplitude_QCD),intent(in) :: this
    integer,intent(in) :: n_gs,n_ew
    integer :: isector
    sector_index=0
    do isector=1,this%n_sectors
       if (this%sector_powers(1,isector).eq.n_gs .and. &
            this%sector_powers(2,isector).eq.n_ew) then
          sector_index=isector
          return
       endif
    enddo
  end function sector_index

  subroutine assign_interaction(lhs,rhs)
    ! sets non-custom 'lhs' = 'rhs' for interactions
    use particles
    implicit none
    type(interaction),intent(inout) :: lhs
    type(interaction),intent(in) :: rhs
    integer :: val_size
    lhs%type=rhs%type
    lhs%chirality=rhs%chirality
    lhs%n_gs=rhs%n_gs
    lhs%n_ew=rhs%n_ew
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
       val_size=size(rhs%val_c)
       allocate(lhs%val_c(1:val_size))
       lhs%val_c(1:val_size)=rhs%val_c(1:val_size)
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       val_size=size(rhs%val_r)
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
    integer :: isize,val_size, lsize
    lhs%type=rhs%type
    lhs%bin=rhs%bin
    lhs%chirality=rhs%chirality
    lhs%n_gs=rhs%n_gs
    lhs%n_ew=rhs%n_ew
    lhs%open_quark_leg=rhs%open_quark_leg
    lhs%ew_pairs=rhs%ew_pairs
    lhs%u1_pairs=rhs%u1_pairs
    lhs%gluon_pairs=rhs%gluon_pairs
    lhs%u1_links=rhs%u1_links
    lhs%gluon_links=rhs%gluon_links
    isize=popcnt(lhs%bin)
    lhs%n_vert=rhs%n_vert
    if (allocated(lhs%iproc%bits)) deallocate(lhs%iproc%bits)
    if (allocated(rhs%iproc%bits)) then
       call lhs%iproc%init(rhs%iproc%n_bits)
       lhs%iproc%bits=rhs%iproc%bits
    endif
    lhs%mass=rhs%mass
    lhs%width=rhs%width
    lhs%ext_cur=rhs%ext_cur
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
       val_size=size(rhs%val_c)
       allocate(lhs%val_c(1:val_size))
       lhs%val_c(1:val_size)=rhs%val_c(1:val_size)
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       val_size=size(rhs%val_r)
       allocate(lhs%val_r(1:val_size))
       lhs%val_r(1:val_size)=rhs%val_r(1:val_size)
    endif
    if (allocated(lhs%fermi_list)) deallocate(lhs%fermi_list)
    if (allocated(rhs%fermi_list)) then
       lsize=size(rhs%fermi_list)
       allocate(lhs%fermi_list(1:lsize))
       lhs%fermi_list(1:lsize)=rhs%fermi_list(1:lsize)
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
    if (allocated(amp%amps_by_order)) deallocate(amp%amps_by_order)
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
    if (allocated(amp%sector_powers)) deallocate(amp%sector_powers)
    if (allocated(amp%sector_curr2amp)) deallocate(amp%sector_curr2amp)
    if (allocated(amp%sector_three_line_partner_curr2amp)) &
         deallocate(amp%sector_three_line_partner_curr2amp)
    if (allocated(amp%sector_present)) deallocate(amp%sector_present)
    if (allocated(amp%sector_retained)) deallocate(amp%sector_retained)
    if (allocated(amp%sector_sign)) deallocate(amp%sector_sign)
    if (allocated(amp%sector_term_start)) deallocate(amp%sector_term_start)
    if (allocated(amp%sector_term_curr2amp)) deallocate(amp%sector_term_curr2amp)
    if (allocated(amp%sector_term_sign)) deallocate(amp%sector_term_sign)
    if (allocated(amp%three_line_partner_curr2amp)) &
         deallocate(amp%three_line_partner_curr2amp)
    if (allocated(amp%i_col_i)) deallocate(amp%i_col_i)
    if (allocated(amp%processes)) deallocate(amp%processes)
    if (allocated(amp%same_flavour_sum)) deallocate(amp%same_flavour_sum)
    if (allocated(amp%same_flavour_sum_operation)) deallocate(amp%same_flavour_sum_operation)
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
    if (allocated(cur%iproc%bits)) deallocate(cur%iproc%bits)
    if (allocated(cur%vertices)) deallocate(cur%vertices)
    if (allocated(cur%order)) deallocate(cur%order)
    if (allocated(cur%spin)) deallocate(cur%spin)
    if (allocated(cur%ext_type)) deallocate(cur%ext_type)
    if (allocated(cur%vertex_sign)) deallocate(cur%vertex_sign)
    if (allocated(cur%val_c)) deallocate(cur%val_c)
    if (allocated(cur%val_r)) deallocate(cur%val_r)
    if (allocated(cur%fermi_list)) deallocate(cur%fermi_list)
  end subroutine finalize_current

end module amplitude_QCD_mod
