module amplitude_QCD_mod
  use bitset_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: complex_value_is_finite,complex_amplitude_is_safe,amplitude_value_limit,&
       valid_external_colour_order,is_coloured_external_particle,&
       valid_external_particle,valid_current_particle
  logical,parameter :: use_symmetry=.true.
  logical,parameter :: use_real_gluons=.false.
  logical,parameter :: use_symm_cm=.true.
  logical,parameter :: use_cm_dict=.true.
  real(kind=8),parameter :: helicity_zero_tolerance=1d-24
  real(kind=8),parameter :: amplitude_value_limit=0.25d0*huge(1d0)**0.25d0
  integer(kind=8),parameter :: max_three_line_color_orders=5000_8
  integer,parameter,public :: max_amplitude_external_particles=20
  integer,parameter :: max_current_dictionary_entries=20000000
  integer,parameter :: max_amplitude_current_records=20000000
  integer,parameter :: max_amplitude_interaction_records=20000000
  integer(kind=8),parameter :: max_amplitude_workspace_bytes=2147483648_8
  integer(kind=8),parameter :: max_amplitude_optimisation_comparisons=100000000_8
  integer(kind=8),parameter :: max_colour_matrix_elements=100000000_8
  type :: current
     ! if adding variables here, also update the finalize_current and assign_current subroutines
     integer :: type=0,bin=0,n_vert=0,chirality=0
     type(bitset) :: iproc
     integer(kind=16) :: ext_cur=0_16
     integer,dimension(:),allocatable :: vertices,order,spin,ext_type
     logical,dimension(:),allocatable :: vertex_sign
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8) :: mass=0d0,width=0d0
     integer,dimension(:),allocatable :: fermi_list
   contains
     final :: finalize_current ! custom deallocation of current
  end type current
  type :: interaction
     ! if adding variables here, also update the finalize_interaction and assign_interaction subroutines
     integer :: type=0,chirality=0
     integer,dimension(2) :: currents=0
     integer,dimension(:),allocatable :: singlet_mv
     complex(kind=8),dimension(:),allocatable :: val_c
     real(kind=8),dimension(:),allocatable :: val_r
     real(kind=8),dimension(2) :: coupl=0d0
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
     integer :: n_cur=0,n_vert=0,imode=0,nColOrd=0,max_pp=0,n_amps=0,nprocs=0
     type(current),dimension(:),allocatable :: current_list
     type(interaction),dimension(:),allocatable :: interaction_list
     complex(kind=8),dimension(:),allocatable :: amps
     real(kind=8),dimension(:),allocatable :: amps_r
     real(kind=8),dimension(:,:),allocatable :: pp,diff_col_vals
     complex(kind=8),dimension(:,:,:),allocatable :: optimisation_current_samples
     integer :: optimisation_sample_count=0
     integer,dimension(:),allocatable :: n_cur_start,n_cur_end,n_vert_start,n_vert_end, &
          pp_bin_to_i,pp_i_to_bin,col_index,n_col_vals,iproc_start,n_sing,n_qqbar
     integer,dimension(:,:),allocatable :: perm,curr2amp,i_col_i,processes,&
          same_flavour_sum,same_flavour_sum_operation,three_line_partner_curr2amp
     integer,dimension(:,:,:),allocatable :: spins,row_index
     logical,dimension(:),allocatable :: include_amp,same_flav
     logical :: lib_created=.false.
     logical :: evaluation_workspace_ready=.false.
   contains
     procedure,public :: init,evaluate,init_col,filter_helicity,write_init_amps_to_file,read_init_amps_from_file &
          ,create_library,record_optimisation_sample,clear_optimisation_samples,optimise_evaluation
     procedure,private :: filter_dead_trees
     final :: finalize_amplitude_QCD ! custom deallocation of amplitude_QCD
  end type amplitude_QCD
contains
  pure logical function is_coloured_external_particle(particle) result(coloured)
    implicit none
    integer,intent(in) :: particle

    ! Avoid ABS here: ABS(HUGE-negative) is not representable for a default
    ! integer supplied by a corrupt process or serialized-amplitude file.
    coloured=(particle.ge.1 .and. particle.le.6) .or. &
         (particle.le.-1 .and. particle.ge.-6) .or. particle.eq.21
  end function is_coloured_external_particle

  pure logical function valid_external_particle(particle) result(valid)
    implicit none
    integer,intent(in) :: particle

    valid=is_coloured_external_particle(particle) .or. &
         (particle.ge.11 .and. particle.le.16) .or. &
         (particle.le.-11 .and. particle.ge.-16) .or. &
         particle.eq.22 .or. particle.eq.23 .or. &
         particle.eq.24 .or. particle.eq.-24 .or. particle.eq.25
  end function valid_external_particle

  pure logical function valid_current_particle(particle) result(valid)
    implicit none
    integer,intent(in) :: particle

    valid=valid_external_particle(particle) .or. particle.eq.99 .or. &
         particle.eq.-21 .or. particle.eq.-23 .or. &
         particle.eq.26 .or. particle.eq.-26 .or. &
         particle.eq.125 .or. particle.eq.126 .or. particle.eq.127
  end function valid_current_particle

  logical function valid_external_colour_order(labels,process) result(valid)
    ! A colour order must contain every coloured external label exactly once
    ! and no singlet label.  Keep this check independent of short-circuit
    ! evaluation so malformed serialized or in-memory metadata cannot be used
    ! as an array subscript by the colour-factor helpers below.
    implicit none
    integer,dimension(:),intent(in) :: labels,process
    integer :: label
    logical :: coloured

    valid=.false.
    if (size(labels).gt.size(process)) return
    if (any(labels.lt.1) .or. any(labels.gt.size(process))) return
    do label=1,size(process)
       coloured=is_coloured_external_particle(process(label))
       if (coloured) then
          if (count(labels.eq.label).ne.1) return
       elseif (any(labels.eq.label)) then
          return
       endif
    enddo
    valid=.true.
  end function valid_external_colour_order

  subroutine reset_amplitude_QCD(this)
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    integer :: i

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
    if (allocated(this%amps)) deallocate(this%amps)
    if (allocated(this%amps_r)) deallocate(this%amps_r)
    if (allocated(this%pp)) deallocate(this%pp)
    if (allocated(this%diff_col_vals)) deallocate(this%diff_col_vals)
    if (allocated(this%optimisation_current_samples)) &
         deallocate(this%optimisation_current_samples)
    if (allocated(this%n_cur_start)) deallocate(this%n_cur_start)
    if (allocated(this%n_cur_end)) deallocate(this%n_cur_end)
    if (allocated(this%n_vert_start)) deallocate(this%n_vert_start)
    if (allocated(this%n_vert_end)) deallocate(this%n_vert_end)
    if (allocated(this%pp_bin_to_i)) deallocate(this%pp_bin_to_i)
    if (allocated(this%pp_i_to_bin)) deallocate(this%pp_i_to_bin)
    if (allocated(this%col_index)) deallocate(this%col_index)
    if (allocated(this%n_col_vals)) deallocate(this%n_col_vals)
    if (allocated(this%iproc_start)) deallocate(this%iproc_start)
    if (allocated(this%n_sing)) deallocate(this%n_sing)
    if (allocated(this%n_qqbar)) deallocate(this%n_qqbar)
    if (allocated(this%perm)) deallocate(this%perm)
    if (allocated(this%curr2amp)) deallocate(this%curr2amp)
    if (allocated(this%three_line_partner_curr2amp)) &
         deallocate(this%three_line_partner_curr2amp)
    if (allocated(this%i_col_i)) deallocate(this%i_col_i)
    if (allocated(this%processes)) deallocate(this%processes)
    if (allocated(this%same_flavour_sum)) deallocate(this%same_flavour_sum)
    if (allocated(this%same_flavour_sum_operation)) &
         deallocate(this%same_flavour_sum_operation)
    if (allocated(this%spins)) deallocate(this%spins)
    if (allocated(this%row_index)) deallocate(this%row_index)
    if (allocated(this%include_amp)) deallocate(this%include_amp)
    if (allocated(this%same_flav)) deallocate(this%same_flav)
    this%n_cur=0
    this%n_vert=0
    this%imode=0
    this%nColOrd=0
    this%max_pp=0
    this%n_amps=0
    this%nprocs=0
    this%optimisation_sample_count=0
    this%lib_created=.false.
    this%evaluation_workspace_ready=.false.
  end subroutine reset_amplitude_QCD

  subroutine init(this,imode,n,n_processes,part,spin,o,pm,valid)
    use math_functions
    use particles
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    type(physics_model),intent(in) :: pm
    integer,intent(in) :: n,imode,n_processes
    integer,dimension(n,n_processes),intent(in) :: part,o
    integer,dimension(0:3,n),intent(in) :: spin
    logical,intent(out),optional :: valid
    integer,dimension(:,:,:),allocatable :: order
    type(current),dimension(:),allocatable :: current_list_local
    type(interaction),dimension(:),allocatable :: interaction_list_local
    integer :: isize,nc,isplit,n1,n2,ic1,ic2,max_cur,max_vert,max_key,ispin,iproc,&
         allocation_status
    integer(kind=8),dimension(:),allocatable :: current_dict
    integer,dimension(:,:),allocatable :: key_to_current
    character(len=256) :: allocation_message

    if (present(valid)) valid=.true.
    if (n.lt.3 .or. n.gt.max_amplitude_external_particles) then
       write (*,*) 'Unsupported number of external particles in amplitude:',&
            n,max_amplitude_external_particles
       stop 1
    endif
    if (n_processes.lt.1) then
       write (*,*) 'An amplitude requires at least one subprocess'
       stop 1
    endif
    do iproc=1,n
       if (spin(0,iproc).lt.1 .or. spin(0,iproc).gt.3) then
          write (*,*) 'Invalid number of external spin states:',iproc,spin(0,iproc)
          stop 1
       endif
       if (any((spin(1:spin(0,iproc),iproc).lt.-1 .or. &
            spin(1:spin(0,iproc),iproc).gt.1) .and. &
            spin(1:spin(0,iproc),iproc).ne.-9)) then
          write (*,*) 'Invalid external spin label:',iproc,spin(1:spin(0,iproc),iproc)
          stop 1
       endif
    enddo
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
    call reset_amplitude_QCD(this)
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
       if (this%nprocs.ne.1) then
          write (*,*) 'For imode==2 there should only be a single process at a time'
          stop 1
       endif
       call create_current_dict()
       allocate(key_to_current(max_key,maskr(this%nprocs)),stat=allocation_status,&
            errmsg=allocation_message)
       call require_amplitude_allocation('current lookup table')
       key_to_current(1:max_key,1:maskr(this%nprocs))=0
    endif

    allocate(current_list_local(max_cur),interaction_list_local(max_vert),&
         this%n_cur_start(n),this%n_cur_end(n),this%n_vert_start(2:n-1),&
         this%n_vert_end(2:n-1),stat=allocation_status,errmsg=allocation_message)
    call require_amplitude_allocation('initial current and interaction workspaces')
   
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
    if (present(valid)) then
       if (.not.valid) return
    endif

    call allocate_and_fill_currents_to_amps_map()
    if (this%n_amps.eq.0) then
       if (present(valid)) then
          valid=.false.
          return
       endif
       write (*,*) 'ERROR: no valid matrix elements found: check your process definition'
       stop 1
    endif
    
    call allocate_current_list_and_interaction_list()

    if (this%imode.eq.1) call allocate_and_fill_spins()
    call allocate_and_fill_colour_permutations()
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
      if (this%n_cur.gt.int(bit_size(int(0,kind=16)),kind=kind(this%n_cur))) then
         write (*,*) 'Too many distinct external currents for the amplitude bit mask:',this%n_cur
         stop 1
      endif
      allocate(current_list_local(this%n_cur)%order(isize),&
           current_list_local(this%n_cur)%ext_type(isize),&
           current_list_local(this%n_cur)%spin(isize),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('external-current metadata')
      current_list_local(this%n_cur)%order(1)=iorder
      current_list_local(this%n_cur)%type=ctype
      current_list_local(this%n_cur)%chirality=ichir
      current_list_local(this%n_cur)%mass=pm%get_mass(current_list_local(this%n_cur)%type)
      current_list_local(this%n_cur)%width=pm%get_width(current_list_local(this%n_cur)%type)
      current_list_local(this%n_cur)%ext_type(1)=current_list_local(this%n_cur)%type
      current_list_local(this%n_cur)%bin=ibset(0,iorder-1) ! give binary label
      current_list_local(this%n_cur)%spin(1)=ispin
      current_list_local(this%n_cur)%n_vert=0
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
      integer :: icur,jcur,i,j,iproc,capacity,new_capacity
      type(bitset) :: proc
      integer,dimension(:,:),allocatable :: curr2amp,tmp_curr2amp
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
      capacity=min(1024,max_amplitude_current_records)
      allocate(curr2amp(1:2,1:capacity),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('current-to-amplitude builder')
      curr2amp=0
      allocate(this%iproc_start(1:this%nprocs+1),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('subprocess amplitude offsets')
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
               if (this%n_amps.ge.max_amplitude_current_records) then
                  write (*,*) 'Amplitude count exceeds supported size:',&
                       max_amplitude_current_records
                  stop 1
               endif
               if (this%n_amps.eq.capacity) then
                  new_capacity=min(max_amplitude_current_records,2*capacity)
                  if (new_capacity.le.capacity) then
                     write (*,*) 'Could not grow current-to-amplitude builder:',capacity
                     stop 1
                  endif
                  allocate(tmp_curr2amp(1:2,1:new_capacity),stat=allocation_status,&
                       errmsg=allocation_message)
                  call require_amplitude_allocation('grown current-to-amplitude builder')
                  tmp_curr2amp=0
                  tmp_curr2amp(:,1:this%n_amps)=curr2amp(:,1:this%n_amps)
                  call move_alloc(tmp_curr2amp,curr2amp)
                  capacity=new_capacity
               endif
               this%n_amps=this%n_amps+1
               curr2amp(1,this%n_amps)=icur
               curr2amp(2,this%n_amps)=jcur
            enddo
         enddo
      enddo

      if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. this%imode.eq.2) then
         if (this%n_amps.gt.max_amplitude_current_records/2) then
            write (*,*) 'Symmetry-expanded amplitude count exceeds supported size:',&
                 this%n_amps,max_amplitude_current_records
            stop 1
         endif
         allocate(this%curr2amp(1:2,1:2*this%n_amps),stat=allocation_status,&
              errmsg=allocation_message)
         call require_amplitude_allocation('symmetry-expanded amplitude map')
         this%curr2amp(1:2,1:this%n_amps)=curr2amp(1:2,1:this%n_amps)
         this%curr2amp(1:2,this%n_amps+1:2*this%n_amps)=curr2amp(1:2,1:this%n_amps)
         this%n_amps=this%n_amps*2
      else
         allocate(this%curr2amp(1:2,1:this%n_amps),stat=allocation_status,&
              errmsg=allocation_message)
         call require_amplitude_allocation('amplitude map')
         this%curr2amp(1:2,1:this%n_amps)=curr2amp(1:2,1:this%n_amps)
      endif
      if (this%imode.eq.3 .and. this%n_amps.ne.1) then
         write (*,*) 'For this%imode==3, there should only be one amplitude',this%n_amps
         write (*,*) this%n_cur_start
         write (*,*) this%n_cur_end
         stop 1
      endif

      allocate(this%same_flavour_sum(this%n_amps,2),&
           this%same_flavour_sum_operation(this%n_amps,2),&
           this%include_amp(1:this%n_amps),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('per-amplitude metadata')
      this%same_flavour_sum=-1
      this%same_flavour_sum_operation=0
      this%iproc_start(this%nprocs+1)=this%n_amps+1

      this%include_amp(:)=.true.
    end subroutine allocate_and_fill_currents_to_amps_map

    subroutine allocate_and_fill_spins()
      implicit none
      integer :: iamp,i,iproc
      allocate(this%spins(n,1,1:this%n_amps),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('amplitude helicities')
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
      integer :: iamp,i,iproc,iamp_to_compare
      ! allocate and fill the colour orders in 'this%perm'. These are simply
      ! the orders of the elements in the 'this%current_list' (with size n-1)
      ! together with the final element). Exception: when there are colour
      ! singlets, they will not be part of the this%perm (while they are part
      ! of the elements in the this%current_list.
      allocate(this%perm(1:n-this%n_sing(1),1:this%n_amps),&
           stat=allocation_status,errmsg=allocation_message)
      call require_amplitude_allocation('amplitude colour permutations')
      if (this%n_sing(1).eq.n) then
         ! A process with no coloured external particles has one empty colour
         ! order.  Its currents still close normally, but there is no label to
         ! append to the zero-extent permutation table.
         return
      endif
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

    subroutine group_three_line_colour_flows()
      ! The alternative three-line canonical order is needed while building
      ! currents, but it can close the same physical open-string flow more
      ! than once. Keep one amplitude for each set of q[g...]qbar strings and
      ! retain the second closure so that both contributions are evaluated.
      implicit none
      integer :: old_amp,flow,flow_count
      integer,dimension(:),allocatable :: representative,partner
      integer,dimension(:,:),allocatable :: old_curr2amp,old_perm,&
           old_same_flavour_sum,old_same_flavour_sum_operation
      logical,dimension(:),allocatable :: old_include_amp

      allocate(representative(this%n_amps),partner(this%n_amps),&
           old_perm(size(this%perm,1),this%n_amps),&
           old_curr2amp(2,this%n_amps),old_include_amp(this%n_amps),&
           old_same_flavour_sum(this%n_amps,2),&
           old_same_flavour_sum_operation(this%n_amps,2),&
           this%three_line_partner_curr2amp(2,this%nColOrd),&
           stat=allocation_status,errmsg=allocation_message)
      ! Keep this check at the allocation site.  Besides retaining a precise
      ! diagnostic, it makes the successful descriptor postcondition visible
      ! to compilers before the local allocatables are copied below.
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate three-line colour-flow grouping: ',&
              trim(allocation_message)
         stop 1
      endif
      if (.not.allocated(partner)) then
         write (*,*) 'Three-line colour-flow partner workspace was not allocated'
         stop 1
      endif
      ! Keep a local snapshot while comparing array sections.  In particular,
      ! avoid repeatedly passing a section of an allocatable component through
      ! nested internal procedures; this also preserves the ungrouped ordering
      ! needed below.
      old_perm=this%perm
      representative=0
      partner=0
      flow_count=0
      do old_amp=1,this%n_amps
         flow=1
         do flow=1,flow_count
            if (representative(flow).lt.1 .or. representative(flow).gt.this%n_amps) then
               write (*,*) 'Invalid three-line colour-flow representative',flow,representative(flow)
               stop 1
            endif
            if (same_three_line_flow(old_perm(:,old_amp),&
                 old_perm(:,representative(flow)))) exit
         enddo
         if (flow.gt.flow_count) then
            flow_count=flow_count+1
            representative(flow_count)=old_amp
         elseif (partner(flow).eq.0) then
            partner(flow)=old_amp
         else
            write (*,*) 'More than two closures for a three-line colour flow',flow
            stop 1
         endif
      enddo
      if (flow_count.ne.this%nColOrd) then
         write (*,*) 'Unexpected three-line colour-flow closure multiplicity',&
              flow_count,this%nColOrd,this%n_amps
         stop 1
      endif

      old_curr2amp=this%curr2amp
      old_include_amp=this%include_amp
      old_same_flavour_sum=this%same_flavour_sum
      old_same_flavour_sum_operation=this%same_flavour_sum_operation
      this%three_line_partner_curr2amp=0
      do flow=1,this%nColOrd
         this%curr2amp(:,flow)=old_curr2amp(:,representative(flow))
         this%perm(:,flow)=old_perm(:,representative(flow))
         this%include_amp(flow)=old_include_amp(representative(flow))
         this%same_flavour_sum(flow,:)=&
              old_same_flavour_sum(representative(flow),:)
         this%same_flavour_sum_operation(flow,:)=&
              old_same_flavour_sum_operation(representative(flow),:)
         if (partner(flow).ne.0) then
            this%three_line_partner_curr2amp(:,flow)=old_curr2amp(:,partner(flow))
         endif
      enddo
      this%n_amps=this%nColOrd
      this%iproc_start(this%nprocs+1)=this%n_amps+1
    end subroutine group_three_line_colour_flows

    logical function same_three_line_flow(left,right)
      implicit none
      integer,dimension(:),intent(in) :: left,right
      integer :: line,other,candidate
      integer,dimension(3) :: left_len,right_len
      integer,dimension(n,3) :: left_lines,right_lines

      call split_three_line_flow(left,left_len,left_lines)
      call split_three_line_flow(right,right_len,right_lines)
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

    subroutine split_three_line_flow(word,line_len,lines)
      implicit none
      integer,dimension(:),intent(in) :: word
      integer,dimension(3),intent(out) :: line_len
      integer,dimension(n,3),intent(out) :: lines
      integer :: pos,label,line,nord

      line_len=0
      lines=0
      line=0
      nord=size(word)
      do pos=1,nord
         label=word(pos)
         if (is_quark_from_order(label,1)) then
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
      integer :: ic,momentum_bins
      integer,dimension(:),allocatable :: pp_i_to_bin
      momentum_bins=maskr(n)
      if (int(momentum_bins,kind=8)*8_8.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Amplitude momentum lookup exceeds supported workspace:',&
              momentum_bins,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(this%pp_bin_to_i(1:momentum_bins),&
           pp_i_to_bin(1:momentum_bins),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('amplitude momentum lookup')
      if (.not.allocated(pp_i_to_bin)) then
         write (*,*) 'Amplitude momentum lookup workspace was not allocated'
         stop 1
      endif
      this%pp_bin_to_i=0
      pp_i_to_bin=0
      this%max_pp=0
      do ic=1,this%n_cur
         if (this%pp_bin_to_i(this%current_list(ic)%bin).eq.0) then
            this%max_pp=this%max_pp+1
            this%pp_bin_to_i(this%current_list(ic)%bin)=this%max_pp
            pp_i_to_bin(this%max_pp)=this%current_list(ic)%bin
         endif
      enddo
      allocate(this%pp(0:3,1:this%max_pp),&
           this%pp_i_to_bin(this%max_pp),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('amplitude momentum storage')
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
         if (present(valid)) then
            valid=.false.
            return
         endif
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
         ! A completely colour-singlet process has one empty colour order.
         this%nColOrd=factorial(max(nglu-1,0))
      elseif (nq.eq.1) then
         this%nColOrd=factorial(nglu)
      elseif (nq.eq.2) then
         color_orders_64=checked_multiply8(factorial8(nglu),int(nglu,kind=8)+1_8,&
              'two-line colour-order count')
         color_orders_64=checked_multiply8(color_orders_64,2_8,'two-line colour-order count')
         if (color_orders_64.gt.int(huge(this%nColOrd),kind=8)) then
            write (*,*) 'Two-line colour basis exceeds supported integer size',color_orders_64
            stop 1
         endif
         this%nColOrd=int(color_orders_64)
      elseif (nq.eq.3) then
         ! Distribute the gluons over three ordered strings and connect the
         ! three quarks to the three antiquarks in all 3! ways.
         color_orders_64=3_8
         if (int(nglu,kind=8)+1_8.gt.max_three_line_color_orders/color_orders_64) then
            color_orders_64=max_three_line_color_orders+1_8
         else
            color_orders_64=color_orders_64*(int(nglu,kind=8)+1_8)
         endif
         if (color_orders_64.le.max_three_line_color_orders) then
            if (int(nglu,kind=8)+2_8.gt.max_three_line_color_orders/color_orders_64) then
               color_orders_64=max_three_line_color_orders+1_8
            else
               color_orders_64=color_orders_64*(int(nglu,kind=8)+2_8)
            endif
         endif
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
      if (this%nprocs.gt.max_bitset_bits .or. &
           int(n,kind=8)*int(this%nprocs,kind=8)*4_8.gt.&
           max_amplitude_workspace_bytes) then
         write (*,*) 'Process group exceeds supported amplitude workspace:',&
              n,this%nprocs,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(this%processes(n,this%nprocs),this%n_qqbar(1:this%nprocs),&
           stat=allocation_status,errmsg=allocation_message)
      call require_amplitude_allocation('amplitude process metadata')
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
         allocate(order(1:n,1,this%nprocs),stat=allocation_status,&
              errmsg=allocation_message)
      elseif (all(this%n_qqbar.le.3)) then
         allocate(order(1:n,2,this%nprocs),stat=allocation_status,&
              errmsg=allocation_message)
      else
         write (*,*) 'Cannot allocated all the needed orders'
         stop 1
      endif
      call require_amplitude_allocation('amplitude construction orders')
      order(1:n,1,1:this%nprocs)=o(1:n,1:this%nprocs)
      allocate(this%n_sing(1:this%nprocs),this%same_flav(1:this%nprocs),&
           stat=allocation_status,errmsg=allocation_message)
      call require_amplitude_allocation('amplitude subprocess classifications')
      do iproc=1,this%nprocs
         if (this%n_qqbar(iproc).eq.3) then
            if (this%imode.eq.2) call canonicalize_three_line_order(iproc)
            call fill_alternative_quark_order(iproc)
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
      integer :: new_max_cur,old_n_cur,ic,allocation_status
      integer(kind=8) :: descriptor_bytes,peak_descriptor_bytes
      type(current),dimension(:),allocatable :: tmp
      character(len=256) :: allocation_message
      old_n_cur=min(this%n_cur-1,max_cur)
      if (max_cur.ge.max_amplitude_current_records) then
         write (*,*) 'Amplitude current workspace exceeds supported size:',&
              max_amplitude_current_records
         stop 1
      endif
      new_max_cur=min(max_amplitude_current_records,2*max_cur)
      descriptor_bytes=max(1_8,int(storage_size(current_list_local),kind=8)/8_8)
      peak_descriptor_bytes=2_8*int(new_max_cur,kind=8)*descriptor_bytes
      if (peak_descriptor_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Amplitude current descriptors exceed supported workspace:',&
              new_max_cur,peak_descriptor_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(tmp(new_max_cur),stat=allocation_status,errmsg=allocation_message)
      if (allocation_status.ne.0) then
         write (*,*) 'Could not grow amplitude current workspace: ',trim(allocation_message)
         stop 1
      endif
      do ic=1,old_n_cur
         tmp(ic)=current_list_local(ic)
      enddo
      do ic=1,old_n_cur
         call finalize_current(current_list_local(ic))
      enddo
      deallocate(current_list_local)
      call move_alloc(tmp,current_list_local)
      max_cur=new_max_cur
    end subroutine increase_max_cur
    
    subroutine increase_max_vert()
      implicit none
      integer :: new_max_vert,iv,allocation_status
      integer(kind=8) :: descriptor_bytes,peak_descriptor_bytes
      type(interaction),dimension(:),allocatable :: tmp
      character(len=256) :: allocation_message
      if (max_vert.ge.max_amplitude_interaction_records) then
         write (*,*) 'Amplitude interaction workspace exceeds supported size:',&
              max_amplitude_interaction_records
         stop 1
      endif
      new_max_vert=min(max_amplitude_interaction_records,2*max_vert)
      descriptor_bytes=max(1_8,int(storage_size(interaction_list_local),kind=8)/8_8)
      peak_descriptor_bytes=2_8*int(new_max_vert,kind=8)*descriptor_bytes
      if (peak_descriptor_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Amplitude interaction descriptors exceed supported workspace:',&
              new_max_vert,peak_descriptor_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(tmp(new_max_vert),stat=allocation_status,errmsg=allocation_message)
      if (allocation_status.ne.0) then
         write (*,*) 'Could not grow amplitude interaction workspace: ',trim(allocation_message)
         stop 1
      endif
      ! copy old list into tmp
      do iv=1,max_vert
         tmp(iv)=interaction_list_local(iv)
      enddo
      ! empty old list
      do iv=1,max_vert
         call finalize_interaction(interaction_list_local(iv))
      enddo
      deallocate(interaction_list_local)
      call move_alloc(tmp,interaction_list_local)
      max_vert=new_max_vert
    end subroutine increase_max_vert
      
    subroutine allocate_current_list_and_interaction_list()
      ! allocate the minimum memory needed for the current_list and
      ! interaction_list to be able to perform the evaluate() procedure.
      implicit none
      integer :: isize,ic,iv
      integer(kind=8) :: descriptor_bytes
      descriptor_bytes=int(this%n_cur,kind=8)*&
           max(1_8,int(storage_size(current_list_local),kind=8)/8_8)+&
           int(this%n_vert,kind=8)*&
           max(1_8,int(storage_size(interaction_list_local),kind=8)/8_8)
      if (descriptor_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Final amplitude descriptors exceed supported workspace:',&
              descriptor_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(this%current_list(1:this%n_cur),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('final amplitude currents')
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
      allocate(this%interaction_list(1:this%n_vert),stat=allocation_status,&
           errmsg=allocation_message)
      call require_amplitude_allocation('final amplitude interactions')
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
      real(kind=8) :: sgn
      integer :: ichir
      if (.not.valid_current_combination())  then
         return
      endif
      do i=1,pm%nint
         if ( current_list_local(ic1)%type.eq.pm%vertex_list(i)%particles(1) .and. &
              current_list_local(ic2)%type.eq.pm%vertex_list(i)%particles(2) ) then
            ichir=vertex_result_chirality(pm%vertex_list(i)%type, &
                 pm%vertex_list(i)%particles(3),pm%vertex_list(i)%coupl)
            if (ichir.eq.-99) cycle
            sgn=1d0
              call add_vertex(pm%vertex_list(i)%type, &
                            pm%vertex_list(i)%particles(3), &
                            sgn*pm%vertex_list(i)%coupl,ichir)
         endif
      enddo
    end subroutine add_if_allowed_threevertex

    integer function ext_from_cur(ic)
      implicit none
      integer :: ic,ncur,nc,ispin
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
            do_c: do c=1,2 ! check both orders for the quark-antiquark blocks when n_qqbar=3
               if (this%n_qqbar(iproc).ne.3 .and. c.eq.2) cycle
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
    
    subroutine add_vertex(itype,ctype,coupl,ichir)
      implicit none
      integer :: itype,ctype,ichir,ic
      real(kind=8),dimension(2) :: coupl
      if (isize.eq.n-1) then
         do ic=this%n_cur_start(n),this%n_cur_end(n)
            if (ctype.eq.anti_current(current_list_local(ic)%type)) exit
         enddo
         if (ic.eq.this%n_cur_end(n)+1) return ! dead tree. Filter already here
      endif
      if (this%n_vert.ge.max_amplitude_interaction_records) then
         write (*,*) 'Amplitude contains too many interactions:',&
              max_amplitude_interaction_records
         stop 1
      endif
      this%n_vert=this%n_vert+1
      if (this%n_vert.gt.max_vert) call increase_max_vert()
      interaction_list_local(this%n_vert)%type=itype
      interaction_list_local(this%n_vert)%chirality=ichir
      interaction_list_local(this%n_vert)%currents(1)=ic1
      interaction_list_local(this%n_vert)%currents(2)=ic2
      interaction_list_local(this%n_vert)%coupl=coupl
      allocate(interaction_list_local(this%n_vert)%singlet_mv(0:isize),&
           stat=allocation_status,errmsg=allocation_message)
      call require_amplitude_allocation('interaction singlet map')
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
      integer :: i,n1,n2,ipos,mv12,nc1,nc2,ns1,ns2,n_mv12_1
      integer,dimension(isize) :: ord
      integer,dimension(:),allocatable :: ord1,spin1,et1,ord2,spin2,et2
      n1=popcnt(current_list_local(ic1)%bin)
      n2=popcnt(current_list_local(ic2)%bin)
      allocate(combine_currents%order(1:isize),&
           combine_currents%spin(1:isize),&
           combine_currents%ext_type(1:isize),ord1(n1),spin1(n1),et1(n1),&
           ord2(n2),spin2(n2),et2(n2),stat=allocation_status,&
           errmsg=allocation_message)
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate combined-current metadata: ',&
              trim(allocation_message)
         stop 1
      endif
      combine_currents%type=ctype
      combine_currents%chirality=ichir
      combine_currents%n_vert=0
      combine_currents%mass=0d0
      combine_currents%width=0d0
      combine_currents%bin=current_list_local(ic1)%bin+current_list_local(ic2)%bin
      combine_currents%iproc=current_list_local(ic1)%iproc.and.current_list_local(ic2)%iproc
      combine_currents%ext_cur=current_list_local(ic1)%ext_cur+current_list_local(ic2)%ext_cur
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
      integer :: min_loc,max_loc
      if (isize.lt.1) then
         write (*,*) 'Cannot test the symmetry order of an empty current'
         stop 1
      endif
      if (.not.allocated(new_current%order) .or. &
           .not.allocated(new_current%ext_type)) then
         write (*,*) 'Cannot test the symmetry order of an incomplete current'
         stop 1
      endif
      if (size(new_current%order).lt.isize .or. &
           size(new_current%ext_type).lt.isize) then
         write (*,*) 'Current metadata is shorter than its symmetry order:',&
              size(new_current%order),size(new_current%ext_type),isize
         stop 1
      endif
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
      max_loc=maxloc(new_current%order(1:isize),dim=1)
      min_loc=minloc(new_current%order(1:isize),dim=1)
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
            if (new_current%chirality.ne.current_list_local(ic)%chirality) cycle
            if (new_current%bin.ne.current_list_local(ic)%bin) cycle
            if (new_current%ext_cur.ne.current_list_local(ic)%ext_cur) cycle
            call append_current_vertex(ic,this%n_vert,vertex_sign)
            return
         enddo
         ! Need a new current
         if (this%n_cur.ge.max_amplitude_current_records) then
            write (*,*) 'Amplitude contains too many currents:',max_amplitude_current_records
            stop 1
         endif
         this%n_cur=this%n_cur+1
         if (this%n_cur.gt.max_cur) call increase_max_cur()
         current_list_local(this%n_cur)=new_current
         current_list_local(this%n_cur)%mass=pm%get_mass(new_current%type)
         current_list_local(this%n_cur)%width=pm%get_width(new_current%type)
         
         if (pm%is_colour_flow_vector(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)),&
                 current_list_local(this%n_cur)%vertex_sign(5*(isize-1)),&
                 stat=allocation_status,errmsg=allocation_message)
         elseif (pm%is_auxiliary_tensor(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(2*(isize-1)),&
                 current_list_local(this%n_cur)%vertex_sign(2*(isize-1)),&
                 stat=allocation_status,errmsg=allocation_message)
         elseif (pm%is_massive_vector(new_current%type)) then
            allocate(current_list_local(this%n_cur)%vertices(5*(isize-1)),&
                 current_list_local(this%n_cur)%vertex_sign(5*(isize-1)),&
                 stat=allocation_status,errmsg=allocation_message)
         else
            allocate(current_list_local(this%n_cur)%vertices(8*(isize-1)),&
                 current_list_local(this%n_cur)%vertex_sign(8*(isize-1)),&
                 stat=allocation_status,errmsg=allocation_message)
         endif
         call require_amplitude_allocation('current interaction references')
         current_list_local(this%n_cur)%vertices(1)=this%n_vert
         current_list_local(this%n_cur)%vertex_sign(1)=vertex_sign
         current_list_local(this%n_cur)%n_vert=1
      elseif (this%imode.eq.2) then
         call get_value(new_current%order,new_current%type,val)
         call solve_dict(val,key)
         ic=key_to_current(key,new_current%iproc%bitset_to_integer())
         if (ic.eq.0) then
            ! initialise new current
            if (this%n_cur.ge.max_amplitude_current_records) then
               write (*,*) 'Amplitude contains too many currents:',max_amplitude_current_records
               stop 1
            endif
            this%n_cur=this%n_cur+1
            if (this%n_cur.gt.max_cur) call increase_max_cur()
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
               allocate(current_list_local(ic)%vertices(5*(isize-1)),&
                    current_list_local(ic)%vertex_sign(5*(isize-1)),&
                    stat=allocation_status,errmsg=allocation_message)
            elseif (pm%is_auxiliary_tensor(new_current%type)) then
               allocate(current_list_local(ic)%vertices(isize-1),&
                    current_list_local(ic)%vertex_sign(isize-1),&
                    stat=allocation_status,errmsg=allocation_message)
            elseif (pm%is_massive_vector(new_current%type)) then
               allocate(current_list_local(ic)%vertices(5*(isize-1)),&
                    current_list_local(ic)%vertex_sign(5*(isize-1)),&
                    stat=allocation_status,errmsg=allocation_message)
            else
               allocate(current_list_local(ic)%vertices(5*(isize-1)),&
                    current_list_local(ic)%vertex_sign(5*(isize-1)),&
                    stat=allocation_status,errmsg=allocation_message)
            endif
            call require_amplitude_allocation('current interaction references')
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
      integer :: old_capacity,new_capacity,allocation_status
      integer,dimension(:),allocatable :: vertices
      logical,dimension(:),allocatable :: vertex_signs
      character(len=256) :: allocation_message

      if (.not.allocated(current_list_local(ic)%vertices)) then
         allocate(current_list_local(ic)%vertices(1),&
              current_list_local(ic)%vertex_sign(1),stat=allocation_status,&
              errmsg=allocation_message)
         if (allocation_status.ne.0) then
            write (*,*) 'Could not allocate current interaction list: ',&
                 trim(allocation_message)
            stop 1
         endif
      elseif (current_list_local(ic)%n_vert.eq.size(current_list_local(ic)%vertices)) then
         old_capacity=size(current_list_local(ic)%vertices)
         if (old_capacity.ge.max_amplitude_interaction_records) then
            write (*,*) 'Current contains too many interaction references:',ic,old_capacity
            stop 1
         endif
         new_capacity=min(max_amplitude_interaction_records,max(1,2*old_capacity))
         allocate(vertices(new_capacity),vertex_signs(new_capacity),stat=allocation_status,&
              errmsg=allocation_message)
         if (allocation_status.ne.0) then
            write (*,*) 'Could not grow current interaction list: ',trim(allocation_message)
            stop 1
         endif
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
      integer :: size,i,j,key,istat,level_entries
      integer(kind=8) :: val,previous_val
      integer,dimension(:),allocatable :: ips_in,ips
      character(len=256) :: allocation_message
      size=n
      max_key=n
      do isize=2,n-1
         size=checked_multiply(size,n-isize+1,'current-dictionary permutation count')
         level_entries=checked_multiply(size,pm%npart,'current-dictionary level size')
         max_key=checked_add(max_key,level_entries,'current-dictionary size')
      enddo
      if (max_key.gt.max_current_dictionary_entries) then
         write (*,*) 'Current dictionary exceeds supported workspace:',max_key,&
              max_current_dictionary_entries
         stop 1
      endif
      allocate(current_dict(max_key),stat=istat,errmsg=allocation_message)
      if (istat.ne.0) then
         write (*,*) 'Could not allocate current dictionary:',trim(allocation_message)
         stop 1
      endif
      key=n  ! skip the external currents.
      size=n
      previous_val=0
      do isize=2,n-1
         size=checked_multiply(size,n-isize+1,'current-dictionary permutation count')
         allocate(ips_in(1:isize),ips(1:isize),stat=istat,&
              errmsg=allocation_message)
         if (istat.ne.0) then
            write (*,*) 'Could not allocate current-dictionary permutation: ',&
                 trim(allocation_message)
            stop 1
         endif
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
               key=checked_add(key,1,'current-dictionary fill position')
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
      if (key.ne.max_key) then
         write (*,*) 'Current dictionary was not filled consistently:',key,max_key
         stop 1
      endif
      max_key=key
    end subroutine create_current_dict

    subroutine get_value(ips,itype,val)
      ! Give every current a unique value. This is based on the
      ! (external) particles that are part of the current as well as
      ! the current type.
      implicit none
      integer,dimension(isize) :: ips
      integer :: j,itype
      integer(kind=8) :: val,offset,base,digit
      if (isize.eq.1) then
         write (*,*) 'current_dict only setup for isize.ge.2',isize
         stop 1
      endif
      if (any(ips.lt.1) .or. any(ips.gt.n)) then
         write (*,*) 'Current dictionary received an out-of-range external label:',ips
         stop 1
      endif
      do j=1,isize-1
         if (any(ips(j+1:isize).eq.ips(j))) then
            write (*,*) 'Current dictionary received duplicate external labels:',ips
            stop 1
         endif
      enddo
      val=0_8
      base=int(n,kind=8)+1_8
      ! Give a unique identifier based on the external
      ! particles. Simply convert the list to an integer with base
      ! equal to the number of external particles.
      do j=1,isize
         digit=int(ips(j),kind=8)
         if (val.gt.(huge(val)-digit)/base) then
            write (*,*) 'Current dictionary key exceeds 64-bit integer range:',ips
            stop 1
         endif
         val=val*base+digit
      enddo
      ! Take the types into account (don't worry about particle
      ! vs. anti-particle, since there should be no confusion given
      ! the (external) particles that are part of the current).
      offset=-1_8
      do j=1,pm%npart
         if (itype.eq.pm%particle_list(j)%type .or. itype.eq.pm%particle_list(j)%anti_type) then
            offset=int(j-1,kind=8)
            exit
         endif
      enddo
      if (offset.lt.0_8 .or. pm%npart.lt.1) then
         write (*,*) 'Unknown current type in current dictionary:',itype
         stop 1
      endif
      if (val.gt.(huge(val)-offset)/int(pm%npart,kind=8)) then
         write (*,*) 'Typed current dictionary key exceeds 64-bit integer range:',ips,itype
         stop 1
      endif
      val=val*int(pm%npart,kind=8)+offset
    end subroutine get_value

    subroutine solve_dict(val,key)
      ! Given the value 'val', find the corresponding key in the
      ! 'current_dict' dictionary. Use a binary search
      ! algorithm. (This only works if the dictionary values are
      ! ordered, and all values only appear once).
      implicit none
      integer :: key,left,middle,right
      integer(kind=8) :: val
      key=0
      left=1
      right=max_key
      do while (left.le.right)
         middle=left+(right-left)/2
         if (current_dict(middle).eq.val) then
            key=middle
            return
         elseif(current_dict(middle).gt.val) then
            right=middle-1
         else
            left=middle+1
         endif
      enddo
      write(*,*) 'ERROR: current is missing from dictionary:',val
      stop 1
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
      if (ispin.ne.-1 .and. ispin.ne.1) return
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

    subroutine require_amplitude_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate ',trim(label),': ',trim(allocation_message)
         stop 1
      endif
    end subroutine require_amplitude_allocation

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
    class(amplitude_QCD),intent(in) :: this
    integer,intent(in) :: n,iunit
    integer :: ic,iv,isize,iamp,iproc,ios
    integer :: process_buffer(n),spin_buffer(n),permutation_buffer(n)
    logical :: unit_open
    character(len=20) :: unit_access,unit_form
    character(len=256) :: io_message

    inquire(unit=iunit,opened=unit_open,access=unit_access,form=unit_form,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect serialized-amplitude output unit:',&
            ios,trim(io_message)
       stop 1
    endif
    if (.not.unit_open .or. trim(unit_access).ne.'STREAM' .or. &
         trim(unit_form).ne.'UNFORMATTED') then
       write (*,*) 'Serialized amplitudes require an open unformatted stream unit:',&
            iunit,trim(unit_access),trim(unit_form)
       stop 1
    endif
    call validate_serialized_write_state()

    write (iunit,iostat=ios,iomsg=io_message) this%n_cur,this%n_vert,this%imode,&
         this%nColOrd,this%max_pp,this%n_amps,this%nprocs
    call require_serialized_write('amplitude header')
    write (iunit,iostat=ios,iomsg=io_message) this%n_cur_start(1:n)
    call require_serialized_write('current-range starts')
    write (iunit,iostat=ios,iomsg=io_message) this%n_cur_end(1:n)
    call require_serialized_write('current-range ends')
    write (iunit,iostat=ios,iomsg=io_message) this%n_vert_start(2:n-1)
    call require_serialized_write('interaction-range starts')
    write (iunit,iostat=ios,iomsg=io_message) this%n_vert_end(2:n-1)
    call require_serialized_write('interaction-range ends')
    ! current_list
    do isize=1,n-1
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          call this%current_list(ic)%iproc%bitset_write_unformatted(iunit,ios,io_message)
          call require_serialized_write('current process mask')
          write (iunit,iostat=ios,iomsg=io_message) this%current_list(ic)%type,&
               this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%chirality, &
               this%current_list(ic)%mass,this%current_list(ic)%width
          call require_serialized_write('current metadata')
          if (this%current_list(ic)%n_vert.gt.0) then
             write (iunit,iostat=ios,iomsg=io_message) &
                  this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)
             call require_serialized_write('current vertex indices')
             write (iunit,iostat=ios,iomsg=io_message) &
                  this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)
             call require_serialized_write('current vertex signs')
          endif
          if (isize.eq.1) then
             write (iunit,iostat=ios,iomsg=io_message) &
                  this%current_list(ic)%order(1),this%current_list(ic)%spin(1)
             call require_serialized_write('external-current labels')
          endif
       enddo
    enddo
    ! interaction_list
    do iv=1,this%n_vert
       write (iunit,iostat=ios,iomsg=io_message) this%interaction_list(iv)%type,&
            this%interaction_list(iv)%chirality,&
            this%interaction_list(iv)%currents(1:2),&
            this%interaction_list(iv)%coupl(1:2)
       call require_serialized_write('interaction metadata')
       if (allocated(this%interaction_list(iv)%singlet_mv)) then
          write (iunit,iostat=ios,iomsg=io_message) &
               this%interaction_list(iv)%singlet_mv(&
               0:this%interaction_list(iv)%singlet_mv(0))
       else
          write (iunit,iostat=ios,iomsg=io_message) 0
       endif
       call require_serialized_write('interaction singlet map')
    enddo
    ! momenta array
    write (iunit,iostat=ios,iomsg=io_message) this%pp_bin_to_i(1:maskr(n))
    call require_serialized_write('binary-to-momentum map')
    write (iunit,iostat=ios,iomsg=io_message) this%pp_i_to_bin(1:this%max_pp)
    call require_serialized_write('momentum-to-binary map')
    ! process specific information
    do iproc=1,this%nprocs
       write (iunit,iostat=ios,iomsg=io_message) this%iproc_start(iproc),&
            this%same_flav(iproc),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       call require_serialized_write('process counters')
       ! A local contiguous snapshot avoids a gfortran bounds-check defect for
       ! rank-reduced allocatable components of polymorphic objects in an I/O
       ! list.  It also makes the serialized column layout explicit.
       process_buffer=this%processes(1:n,iproc)
       write (iunit,iostat=ios,iomsg=io_message) process_buffer
       call require_serialized_write('process particle identities')
    enddo
    write(iunit,iostat=ios,iomsg=io_message) this%iproc_start(this%nprocs+1)
    call require_serialized_write('final subprocess offset')
    ! amp specific information
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          write (iunit,iostat=ios,iomsg=io_message) this%include_amp(iamp),&
               this%same_flavour_sum(iamp,1),this%same_flavour_sum(iamp,2),&
               this%same_flavour_sum_operation(iamp,1),&
               this%same_flavour_sum_operation(iamp,2)
          call require_serialized_write('amplitude inclusion metadata')
          spin_buffer=this%spins(1:n,1,iamp)
          write (iunit,iostat=ios,iomsg=io_message) spin_buffer
          call require_serialized_write('amplitude helicities')
          permutation_buffer(1:n-this%n_sing(1))=&
               this%perm(1:n-this%n_sing(1),iamp)
          write (iunit,iostat=ios,iomsg=io_message) &
               permutation_buffer(1:n-this%n_sing(1))
          call require_serialized_write('amplitude permutation')
          if (.not.this%same_flav(iproc)) then
             write (iunit,iostat=ios,iomsg=io_message) &
                  this%curr2amp(1,iamp),this%curr2amp(2,iamp)
             call require_serialized_write('current-to-amplitude map')
          endif
       enddo
    enddo
  contains
    subroutine validate_serialized_write_state()
      implicit none
      integer :: current_index,interaction_index,process_index,amplitude_index,&
           label_index,momentum_index,perm_extent,momentum_bins,n_quark,n_antiquark,&
           n_singlet
      integer :: local_process(n),local_permutation(n)
      logical :: first_is_lepton,current_is_lepton

      momentum_bins=maskr(n)
      if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
           this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
           this%n_vert.lt.0 .or. &
           this%n_vert.gt.max_amplitude_interaction_records .or. &
           this%imode.lt.1 .or. this%imode.gt.3 .or. &
           this%nColOrd.lt.1 .or. &
           this%nColOrd.gt.max_current_dictionary_entries .or. &
           this%max_pp.lt.1 .or. this%max_pp.gt.momentum_bins .or. &
           this%n_amps.lt.1 .or. &
           this%n_amps.gt.max_amplitude_current_records .or. &
           this%nprocs.lt.1 .or. this%nprocs.gt.max_bitset_bits) then
         call invalid_serialized_write_state('invalid amplitude dimensions')
      endif
      if (this%imode.eq.2 .and. this%nprocs.ne.1) &
           call invalid_serialized_write_state(&
           'mode 2 requires exactly one subprocess')

      if (.not.allocated(this%n_cur_start) .or. &
           .not.allocated(this%n_cur_end) .or. &
           .not.allocated(this%n_vert_start) .or. &
           .not.allocated(this%n_vert_end)) &
           call invalid_serialized_write_state('missing current/interaction ranges')
      if (lbound(this%n_cur_start,1).ne.1 .or. &
           ubound(this%n_cur_start,1).ne.n .or. &
           lbound(this%n_cur_end,1).ne.1 .or. &
           ubound(this%n_cur_end,1).ne.n .or. &
           lbound(this%n_vert_start,1).ne.2 .or. &
           ubound(this%n_vert_start,1).ne.n-1 .or. &
           lbound(this%n_vert_end,1).ne.2 .or. &
           ubound(this%n_vert_end,1).ne.n-1) &
           call invalid_serialized_write_state(&
           'incompatible current/interaction range bounds')
      call validate_serialized_write_ranges()

      if (.not.allocated(this%current_list)) &
           call invalid_serialized_write_state('missing current descriptors')
      if (lbound(this%current_list,1).ne.1 .or. &
           ubound(this%current_list,1).lt.this%n_cur) &
           call invalid_serialized_write_state('short current descriptor array')
      if (this%n_vert.gt.0) then
         if (.not.allocated(this%interaction_list)) &
              call invalid_serialized_write_state('missing interaction descriptors')
         if (lbound(this%interaction_list,1).ne.1 .or. &
              ubound(this%interaction_list,1).lt.this%n_vert) &
              call invalid_serialized_write_state(&
              'short interaction descriptor array')
      endif

      do current_index=1,this%n_cur
         if (this%current_list(current_index)%iproc%n_bits.ne.this%nprocs) &
              call invalid_serialized_write_state(&
              'current has an incompatible process mask')
         if (.not.valid_current_particle(&
              this%current_list(current_index)%type)) &
              call invalid_serialized_write_state(&
              'current has an unsupported particle type')
         if (this%current_list(current_index)%bin.lt.1 .or. &
              this%current_list(current_index)%bin.gt.momentum_bins .or. &
              this%current_list(current_index)%chirality.lt.-1 .or. &
              this%current_list(current_index)%chirality.gt.1 .or. &
              this%current_list(current_index)%n_vert.lt.0 .or. &
              this%current_list(current_index)%n_vert.gt.this%n_vert) &
              call invalid_serialized_write_state('invalid current metadata')
         if (.not.ieee_is_finite(this%current_list(current_index)%mass) .or. &
              .not.ieee_is_finite(this%current_list(current_index)%width)) &
              call invalid_serialized_write_state(&
              'non-finite current mass or width')
         if (this%current_list(current_index)%mass.lt.0d0 .or. &
              this%current_list(current_index)%width.lt.0d0 .or. &
              this%current_list(current_index)%mass.gt.amplitude_value_limit .or. &
              this%current_list(current_index)%width.gt.amplitude_value_limit) &
              call invalid_serialized_write_state('unsafe current mass or width')
         if (this%current_list(current_index)%n_vert.gt.0) then
            if (.not.allocated(this%current_list(current_index)%vertices) .or. &
                 .not.allocated(this%current_list(current_index)%vertex_sign)) &
                 call invalid_serialized_write_state(&
                 'missing current interaction entries')
            if (lbound(this%current_list(current_index)%vertices,1).gt.1 .or. &
                 ubound(this%current_list(current_index)%vertices,1).lt.&
                 this%current_list(current_index)%n_vert .or. &
                 lbound(this%current_list(current_index)%vertex_sign,1).gt.1 .or. &
                 ubound(this%current_list(current_index)%vertex_sign,1).lt.&
                 this%current_list(current_index)%n_vert) &
                 call invalid_serialized_write_state(&
                 'short current interaction entries')
            if (any(this%current_list(current_index)%vertices(&
                 1:this%current_list(current_index)%n_vert).lt.1) .or. &
                 any(this%current_list(current_index)%vertices(&
                 1:this%current_list(current_index)%n_vert).gt.this%n_vert)) &
                 call invalid_serialized_write_state(&
                 'out-of-range current interaction entry')
         endif
         if (current_index.ge.this%n_cur_start(1) .and. &
              current_index.le.this%n_cur_end(1)) then
            if (.not.allocated(this%current_list(current_index)%order) .or. &
                 .not.allocated(this%current_list(current_index)%spin)) &
                 call invalid_serialized_write_state(&
                 'missing external-current labels')
            if (lbound(this%current_list(current_index)%order,1).gt.1 .or. &
                 ubound(this%current_list(current_index)%order,1).lt.1 .or. &
                 lbound(this%current_list(current_index)%spin,1).gt.1 .or. &
                 ubound(this%current_list(current_index)%spin,1).lt.1) &
                 call invalid_serialized_write_state(&
                 'empty external-current labels')
            if (this%current_list(current_index)%order(1).lt.1 .or. &
                 this%current_list(current_index)%order(1).gt.n .or. &
                 this%current_list(current_index)%spin(1).lt.-9 .or. &
                 this%current_list(current_index)%spin(1).gt.1 .or. &
                 (this%current_list(current_index)%spin(1).lt.-1 .and. &
                 this%current_list(current_index)%spin(1).ne.-9)) &
                 call invalid_serialized_write_state(&
                 'invalid external-current labels')
         endif
      enddo

      do interaction_index=1,this%n_vert
         if (this%interaction_list(interaction_index)%type.lt.0 .or. &
              this%interaction_list(interaction_index)%type.gt.24 .or. &
              this%interaction_list(interaction_index)%chirality.lt.-1 .or. &
              this%interaction_list(interaction_index)%chirality.gt.1 .or. &
              any(this%interaction_list(interaction_index)%currents.lt.1) .or. &
              any(this%interaction_list(interaction_index)%currents.gt.this%n_cur)) &
              call invalid_serialized_write_state('invalid interaction metadata')
         if (.not.all(ieee_is_finite(&
              this%interaction_list(interaction_index)%coupl))) &
              call invalid_serialized_write_state('non-finite interaction coupling')
         if (any(abs(this%interaction_list(interaction_index)%coupl).gt.&
              amplitude_value_limit)) &
              call invalid_serialized_write_state('unsafe interaction coupling')
         if (allocated(this%interaction_list(interaction_index)%singlet_mv)) then
            if (lbound(this%interaction_list(interaction_index)%singlet_mv,1).ne.0 .or. &
                 ubound(this%interaction_list(interaction_index)%singlet_mv,1).lt.0) &
                 call invalid_serialized_write_state(&
                 'malformed interaction singlet map')
            if (this%interaction_list(interaction_index)%singlet_mv(0).lt.0 .or. &
                 this%interaction_list(interaction_index)%singlet_mv(0).gt.n .or. &
                 ubound(this%interaction_list(interaction_index)%singlet_mv,1).lt.&
                 this%interaction_list(interaction_index)%singlet_mv(0)) &
                 call invalid_serialized_write_state(&
                 'invalid interaction singlet-map extent')
            if (this%interaction_list(interaction_index)%singlet_mv(0).gt.0) then
               if (any(this%interaction_list(interaction_index)%singlet_mv(&
                    1:this%interaction_list(interaction_index)%singlet_mv(0)).lt.1) .or. &
                    any(this%interaction_list(interaction_index)%singlet_mv(&
                    1:this%interaction_list(interaction_index)%singlet_mv(0)).gt.n)) &
                    call invalid_serialized_write_state(&
                    'out-of-range interaction singlet map')
            endif
         endif
      enddo

      if (.not.allocated(this%pp_bin_to_i) .or. &
           .not.allocated(this%pp_i_to_bin)) &
           call invalid_serialized_write_state('missing momentum maps')
      if (lbound(this%pp_bin_to_i,1).ne.1 .or. &
           ubound(this%pp_bin_to_i,1).ne.momentum_bins .or. &
           lbound(this%pp_i_to_bin,1).ne.1 .or. &
           ubound(this%pp_i_to_bin,1).ne.this%max_pp) &
           call invalid_serialized_write_state('incompatible momentum maps')
      if (any(this%pp_bin_to_i.lt.0) .or. &
           any(this%pp_bin_to_i.gt.this%max_pp) .or. &
           any(this%pp_i_to_bin.lt.1) .or. &
           any(this%pp_i_to_bin.gt.momentum_bins)) &
           call invalid_serialized_write_state('out-of-range momentum map')
      do momentum_index=1,this%max_pp
         if (this%pp_bin_to_i(this%pp_i_to_bin(momentum_index)).ne.&
              momentum_index) call invalid_serialized_write_state(&
              'inconsistent momentum maps')
      enddo
      do current_index=1,this%n_cur
         if (this%pp_bin_to_i(this%current_list(current_index)%bin).lt.1) &
              call invalid_serialized_write_state('unmapped current momentum')
      enddo

      if (.not.allocated(this%n_sing) .or. &
           .not.allocated(this%processes) .or. &
           .not.allocated(this%iproc_start) .or. &
           .not.allocated(this%same_flav) .or. &
           .not.allocated(this%n_qqbar)) &
           call invalid_serialized_write_state('missing process metadata')
      if (lbound(this%n_sing,1).ne.1 .or. &
           ubound(this%n_sing,1).ne.this%nprocs .or. &
           lbound(this%n_qqbar,1).ne.1 .or. &
           ubound(this%n_qqbar,1).ne.this%nprocs .or. &
           lbound(this%same_flav,1).ne.1 .or. &
           ubound(this%same_flav,1).ne.this%nprocs .or. &
           lbound(this%iproc_start,1).ne.1 .or. &
           ubound(this%iproc_start,1).ne.this%nprocs+1 .or. &
           lbound(this%processes,1).ne.1 .or. &
           ubound(this%processes,1).ne.n .or. &
           lbound(this%processes,2).ne.1 .or. &
           ubound(this%processes,2).ne.this%nprocs) &
           call invalid_serialized_write_state(&
           'incompatible process metadata bounds')
      if (this%iproc_start(1).ne.1 .or. &
           this%iproc_start(this%nprocs+1).ne.this%n_amps+1 .or. &
           any(this%iproc_start(2:this%nprocs+1).lt.&
           this%iproc_start(1:this%nprocs))) &
           call invalid_serialized_write_state('invalid subprocess offsets')

      do process_index=1,this%nprocs
         n_quark=0
         n_antiquark=0
         n_singlet=0
         do label_index=1,n
            if (.not.valid_external_particle(&
                 this%processes(label_index,process_index))) &
                 call invalid_serialized_write_state(&
                 'unsupported external particle')
            if (label_index.le.2) then
               if (this%processes(label_index,process_index).ge.1 .and. &
                    this%processes(label_index,process_index).le.6) &
                    n_antiquark=n_antiquark+1
               if (this%processes(label_index,process_index).le.-1 .and. &
                    this%processes(label_index,process_index).ge.-6) &
                    n_quark=n_quark+1
            else
               if (this%processes(label_index,process_index).ge.1 .and. &
                    this%processes(label_index,process_index).le.6) &
                    n_quark=n_quark+1
               if (this%processes(label_index,process_index).le.-1 .and. &
                    this%processes(label_index,process_index).ge.-6) &
                    n_antiquark=n_antiquark+1
            endif
            if (.not.is_coloured_external_particle(&
                 this%processes(label_index,process_index))) &
                 n_singlet=n_singlet+1
            if (process_index.gt.1) then
               first_is_lepton=(this%processes(label_index,1).ge.11 .and. &
                    this%processes(label_index,1).le.16) .or. &
                    (this%processes(label_index,1).le.-11 .and. &
                    this%processes(label_index,1).ge.-16)
               current_is_lepton=&
                    (this%processes(label_index,process_index).ge.11 .and. &
                    this%processes(label_index,process_index).le.16) .or. &
                    (this%processes(label_index,process_index).le.-11 .and. &
                    this%processes(label_index,process_index).ge.-16)
               if (first_is_lepton.neqv.current_is_lepton) &
                    call invalid_serialized_write_state(&
                    'inconsistent lepton positions')
            endif
         enddo
         if (n_quark.ne.n_antiquark .or. &
              n_quark.ne.this%n_qqbar(process_index) .or. &
              this%n_qqbar(process_index).lt.0 .or. &
              this%n_qqbar(process_index).gt.3) &
              call invalid_serialized_write_state(&
              'inconsistent subprocess quark-line count')
         if (n_singlet.ne.this%n_sing(process_index)) &
              call invalid_serialized_write_state(&
              'inconsistent subprocess singlet count')
      enddo
      if (any(this%n_sing.ne.this%n_sing(1))) &
           call invalid_serialized_write_state(&
           'subprocesses disagree on colour-singlet count')

      if (.not.allocated(this%include_amp) .or. &
           .not.allocated(this%same_flavour_sum) .or. &
           .not.allocated(this%same_flavour_sum_operation) .or. &
           .not.allocated(this%spins) .or. .not.allocated(this%perm) .or. &
           .not.allocated(this%curr2amp)) &
           call invalid_serialized_write_state('missing per-amplitude metadata')
      perm_extent=n-this%n_sing(1)
      if (lbound(this%include_amp,1).ne.1 .or. &
           ubound(this%include_amp,1).lt.this%n_amps .or. &
           lbound(this%same_flavour_sum,1).ne.1 .or. &
           ubound(this%same_flavour_sum,1).lt.this%n_amps .or. &
           lbound(this%same_flavour_sum,2).ne.1 .or. &
           ubound(this%same_flavour_sum,2).ne.2 .or. &
           lbound(this%same_flavour_sum_operation,1).ne.1 .or. &
           ubound(this%same_flavour_sum_operation,1).lt.this%n_amps .or. &
           lbound(this%same_flavour_sum_operation,2).ne.1 .or. &
           ubound(this%same_flavour_sum_operation,2).ne.2 .or. &
           lbound(this%spins,1).ne.1 .or. ubound(this%spins,1).ne.n .or. &
           lbound(this%spins,2).ne.1 .or. ubound(this%spins,2).ne.1 .or. &
           lbound(this%spins,3).ne.1 .or. &
           ubound(this%spins,3).lt.this%n_amps .or. &
           lbound(this%perm,1).ne.1 .or. &
           ubound(this%perm,1).ne.perm_extent .or. &
           lbound(this%perm,2).ne.1 .or. ubound(this%perm,2).lt.this%n_amps .or. &
           lbound(this%curr2amp,1).ne.1 .or. &
           ubound(this%curr2amp,1).ne.2 .or. &
           lbound(this%curr2amp,2).ne.1 .or. &
           ubound(this%curr2amp,2).lt.this%n_amps) &
           call invalid_serialized_write_state(&
           'incompatible per-amplitude metadata bounds')
      if (any(this%spins(1:n,1,1:this%n_amps).lt.-1) .or. &
           any(this%spins(1:n,1,1:this%n_amps).gt.1)) &
           call invalid_serialized_write_state('invalid amplitude helicity')
      if (any(this%same_flavour_sum(1:this%n_amps,1:2).lt.-1) .or. &
           any(this%same_flavour_sum(1:this%n_amps,1:2).gt.this%n_amps)) &
           call invalid_serialized_write_state('invalid same-flavour map')
      do amplitude_index=1,this%n_amps
         if ((this%same_flavour_sum(amplitude_index,1).gt.0) .neqv. &
              (this%same_flavour_sum(amplitude_index,2).gt.0)) &
              call invalid_serialized_write_state(&
              'incomplete same-flavour map')
         do label_index=1,2
            if (this%same_flavour_sum(amplitude_index,label_index).gt.0) then
               if (this%same_flavour_sum_operation(amplitude_index,label_index).lt.0 .or. &
                    this%same_flavour_sum_operation(amplitude_index,label_index).gt.7) &
                    call invalid_serialized_write_state(&
                    'invalid same-flavour operation')
            elseif (this%same_flavour_sum_operation(amplitude_index,label_index).ne.0) then
               call invalid_serialized_write_state(&
                    'orphan same-flavour operation')
            endif
         enddo
      enddo

      do process_index=1,this%nprocs
         local_process=this%processes(1:n,process_index)
         do amplitude_index=this%iproc_start(process_index),&
              this%iproc_start(process_index+1)-1
            if (perm_extent.gt.0) &
                 local_permutation(1:perm_extent)=&
                 this%perm(1:perm_extent,amplitude_index)
            if (.not.valid_external_colour_order(&
                 local_permutation(1:perm_extent),local_process)) &
                 call invalid_serialized_write_state(&
                 'invalid amplitude colour permutation')
            if (.not.this%same_flav(process_index)) then
               if (any(this%curr2amp(1:2,amplitude_index).lt.1) .or. &
                    any(this%curr2amp(1:2,amplitude_index).gt.this%n_cur)) &
                    call invalid_serialized_write_state(&
                    'invalid current-to-amplitude map')
            endif
         enddo
      enddo
    end subroutine validate_serialized_write_state

    subroutine validate_serialized_write_ranges()
      implicit none
      integer :: range_index

      if (this%n_cur_start(1).ne.1 .or. &
           this%n_cur_end(n-1).ne.this%n_cur) &
           call invalid_serialized_write_state('incomplete current ranges')
      do range_index=1,n-1
         if (this%n_cur_start(range_index).lt.1 .or. &
              this%n_cur_start(range_index).gt.this%n_cur+1 .or. &
              this%n_cur_end(range_index).lt.0 .or. &
              this%n_cur_end(range_index).gt.this%n_cur .or. &
              this%n_cur_start(range_index).gt.&
              this%n_cur_end(range_index)+1) &
              call invalid_serialized_write_state('out-of-range current range')
         if (range_index.gt.1) then
            if (this%n_cur_start(range_index).ne.&
                 this%n_cur_end(range_index-1)+1) &
                 call invalid_serialized_write_state(&
                 'non-contiguous current ranges')
         endif
      enddo
      if (this%n_cur_start(n).lt.this%n_cur_start(1) .or. &
           this%n_cur_end(n).gt.this%n_cur_end(1) .or. &
           this%n_cur_start(n).gt.this%n_cur_end(n)+1) &
           call invalid_serialized_write_state('invalid closing-current range')
      if (this%n_vert.eq.0) then
         if (any(this%n_vert_start(2:n-1).ne.1) .or. &
              any(this%n_vert_end(2:n-1).ne.0)) &
              call invalid_serialized_write_state(&
              'invalid empty interaction ranges')
      else
         if (this%n_vert_start(2).ne.1 .or. &
              this%n_vert_end(n-1).ne.this%n_vert) &
              call invalid_serialized_write_state('incomplete interaction ranges')
         do range_index=2,n-1
            if (this%n_vert_start(range_index).lt.1 .or. &
                 this%n_vert_start(range_index).gt.this%n_vert+1 .or. &
                 this%n_vert_end(range_index).lt.0 .or. &
                 this%n_vert_end(range_index).gt.this%n_vert .or. &
                 this%n_vert_start(range_index).gt.&
                 this%n_vert_end(range_index)+1) &
                 call invalid_serialized_write_state(&
                 'out-of-range interaction range')
            if (range_index.gt.2) then
               if (this%n_vert_start(range_index).ne.&
                    this%n_vert_end(range_index-1)+1) &
                    call invalid_serialized_write_state(&
                    'non-contiguous interaction ranges')
            endif
         enddo
      endif
    end subroutine validate_serialized_write_ranges

    subroutine invalid_serialized_write_state(message)
      implicit none
      character(len=*),intent(in) :: message

      write (*,*) 'Invalid amplitude serialization state: ',trim(message)
      stop 1
    end subroutine invalid_serialized_write_state

    subroutine require_serialized_write(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Could not serialize amplitude while writing ',trim(label),&
              ios,trim(io_message)
         stop 1
      endif
    end subroutine require_serialized_write
  end subroutine write_init_amps_to_file

  subroutine read_init_amps_from_file(this,n,iunit)
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    integer,intent(in) :: n,iunit
    integer :: ic,iv,isize,iamp,iproc,itmp,ios,allocation_status,perm_extent,&
         momentum_bins,label_index
    integer :: process_buffer(n),spin_buffer(n),permutation_buffer(n)
    integer(kind=8) :: workspace_bytes
    integer(kind=8),parameter :: max_serialized_amplitude_workspace_bytes=2147483648_8
    logical :: unit_open
    character(len=20) :: unit_access,unit_form
    character(len=256) :: io_message,allocation_message

    if (n.lt.3 .or. n.gt.max_amplitude_external_particles) then
       write (*,*) 'Invalid external multiplicity in serialized amplitude:',n
       stop 1
    endif
    inquire(unit=iunit,opened=unit_open,access=unit_access,form=unit_form,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect serialized-amplitude input unit:',&
            ios,trim(io_message)
       stop 1
    endif
    if (.not.unit_open .or. trim(unit_access).ne.'STREAM' .or. &
         trim(unit_form).ne.'UNFORMATTED') then
       write (*,*) 'Serialized amplitudes require an open unformatted stream unit:',&
            iunit,trim(unit_access),trim(unit_form)
       stop 1
    endif
    call reset_amplitude_QCD(this)
    workspace_bytes=0_8
    read (iunit,iostat=ios,iomsg=io_message) this%n_cur,this%n_vert,this%imode,&
         this%nColOrd,this%max_pp,this%n_amps,this%nprocs
    call require_serialized_read('amplitude header')
    momentum_bins=maskr(n)
    if (this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
         this%n_vert.lt.0 .or. this%n_vert.gt.max_amplitude_interaction_records .or. &
         this%imode.lt.1 .or. this%imode.gt.3 .or. &
         this%nColOrd.lt.1 .or. this%nColOrd.gt.max_current_dictionary_entries .or. &
         this%max_pp.lt.1 .or. this%max_pp.gt.momentum_bins .or. &
         this%n_amps.lt.1 .or. this%n_amps.gt.max_amplitude_current_records .or. &
         this%nprocs.lt.1 .or. this%nprocs.gt.max_bitset_bits) then
       write (*,*) 'Invalid dimensions in serialized amplitude:',this%n_cur,this%n_vert,&
            this%imode,this%nColOrd,this%max_pp,this%n_amps,this%nprocs
       stop 1
    endif
    if (this%imode.eq.2 .and. this%nprocs.ne.1) then
       write (*,*) 'Serialized mode-2 amplitude has multiple subprocesses:',&
            this%nprocs
       stop 1
    endif
    call reserve_serialized_workspace(4_8*int(n,kind=8)*4_8,&
         'amplitude range tables')
    allocate(this%n_cur_start(1:n),this%n_cur_end(1:n),&
         this%n_vert_start(2:n-1),this%n_vert_end(2:n-1),&
         stat=allocation_status,errmsg=allocation_message)
    call require_serialized_allocation('amplitude range tables')
    read (iunit,iostat=ios,iomsg=io_message) this%n_cur_start(1:n)
    call require_serialized_read('current-range starts')
    read (iunit,iostat=ios,iomsg=io_message) this%n_cur_end(1:n)
    call require_serialized_read('current-range ends')
    read (iunit,iostat=ios,iomsg=io_message) this%n_vert_start(2:n-1)
    call require_serialized_read('interaction-range starts')
    read (iunit,iostat=ios,iomsg=io_message) this%n_vert_end(2:n-1)
    call require_serialized_read('interaction-range ends')
    call validate_serialized_ranges()
    ! current_list
    call reserve_serialized_workspace(256_8*int(this%n_cur,kind=8),&
         'serialized current descriptors')
    allocate(this%current_list(this%n_cur),stat=allocation_status,&
         errmsg=allocation_message)
    call require_serialized_allocation('serialized current descriptors')
    do isize=1,n-1
       do ic=this%n_cur_start(isize),this%n_cur_end(isize)
          call this%current_list(ic)%iproc%bitset_read_unformatted(iunit,ios,io_message)
          call require_serialized_read('current process mask')
          if (this%current_list(ic)%iproc%n_bits.ne.this%nprocs) then
             write (*,*) 'Serialized current has an incompatible process mask:',&
                  ic,this%current_list(ic)%iproc%n_bits,this%nprocs
             stop 1
          endif
          read (iunit,iostat=ios,iomsg=io_message) this%current_list(ic)%type,&
               this%current_list(ic)%bin,this%current_list(ic)%n_vert, &
               this%current_list(ic)%chirality,this%current_list(ic)%mass,this%current_list(ic)%width
          call require_serialized_read('current metadata')
          if (.not.valid_current_particle(this%current_list(ic)%type) .or. &
               this%current_list(ic)%bin.lt.1 .or. &
               this%current_list(ic)%bin.gt.momentum_bins .or. &
               this%current_list(ic)%n_vert.lt.0 .or. &
               this%current_list(ic)%n_vert.gt.this%n_vert .or. &
               this%current_list(ic)%chirality.lt.-1 .or. &
               this%current_list(ic)%chirality.gt.1) then
             write (*,*) 'Invalid serialized current metadata:',ic,&
                  this%current_list(ic)%type,this%current_list(ic)%bin,&
                  this%current_list(ic)%n_vert,this%current_list(ic)%chirality
             stop 1
          endif
          if (.not.ieee_is_finite(this%current_list(ic)%mass) .or. &
               .not.ieee_is_finite(this%current_list(ic)%width)) then
             write (*,*) 'Non-finite mass or width in serialized current:',ic
             stop 1
          endif
          if (this%current_list(ic)%mass.gt.amplitude_value_limit .or. &
               this%current_list(ic)%width.gt.amplitude_value_limit .or. &
               this%current_list(ic)%mass.lt.0d0 .or. &
               this%current_list(ic)%width.lt.0d0) then
             write (*,*) 'Unsafe mass or width in serialized current:',ic,&
                  this%current_list(ic)%mass,this%current_list(ic)%width
             stop 1
          endif
          call reserve_serialized_workspace(&
               5_8*int(this%current_list(ic)%n_vert,kind=8),&
               'serialized current vertices')
          allocate(this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert),&
               this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert),&
               stat=allocation_status,errmsg=allocation_message)
          call require_serialized_allocation('serialized current vertices')
          read (iunit,iostat=ios,iomsg=io_message) &
               this%current_list(ic)%vertices(1:this%current_list(ic)%n_vert)
          call require_serialized_read('current vertex indices')
          read (iunit,iostat=ios,iomsg=io_message) &
               this%current_list(ic)%vertex_sign(1:this%current_list(ic)%n_vert)
          call require_serialized_read('current vertex signs')
          if (this%current_list(ic)%n_vert.gt.0) then
             if (any(this%current_list(ic)%vertices.lt.1) .or. &
                  any(this%current_list(ic)%vertices.gt.this%n_vert)) then
                write (*,*) 'Out-of-range interaction index in serialized current:',ic
                stop 1
             endif
          endif
          if (isize.eq.1) then
             call reserve_serialized_workspace(8_8,'external-current labels')
             allocate(this%current_list(ic)%order(1),this%current_list(ic)%spin(1),&
                  stat=allocation_status,errmsg=allocation_message)
             call require_serialized_allocation('external-current labels')
             read (iunit,iostat=ios,iomsg=io_message) this%current_list(ic)%order(1),&
                  this%current_list(ic)%spin(1)
             call require_serialized_read('external-current labels')
             if (this%current_list(ic)%order(1).lt.1 .or. &
                  this%current_list(ic)%order(1).gt.n .or. &
                  this%current_list(ic)%spin(1).lt.-9 .or. &
                  this%current_list(ic)%spin(1).gt.1 .or. &
                  (this%current_list(ic)%spin(1).lt.-1 .and. &
                  this%current_list(ic)%spin(1).ne.-9)) then
                write (*,*) 'Invalid external-current labels in serialized amplitude:',ic
                stop 1
             endif
          endif
       enddo
    enddo
    ! interaction_list
    call reserve_serialized_workspace(128_8*int(this%n_vert,kind=8),&
         'serialized interaction descriptors')
    allocate(this%interaction_list(1:this%n_vert),stat=allocation_status,&
         errmsg=allocation_message)
    call require_serialized_allocation('serialized interaction descriptors')
    do iv=1,this%n_vert
       read (iunit,iostat=ios,iomsg=io_message) this%interaction_list(iv)%type,&
            this%interaction_list(iv)%chirality,&
            this%interaction_list(iv)%currents(1:2),this%interaction_list(iv)%coupl(1:2),itmp
       call require_serialized_read('interaction metadata')
       if (this%interaction_list(iv)%type.lt.0 .or. &
            this%interaction_list(iv)%type.gt.24 .or. &
            this%interaction_list(iv)%chirality.lt.-1 .or. &
            this%interaction_list(iv)%chirality.gt.1 .or. &
            any(this%interaction_list(iv)%currents.lt.1) .or. &
            any(this%interaction_list(iv)%currents.gt.this%n_cur) .or. &
            itmp.lt.0 .or. itmp.gt.n) then
          write (*,*) 'Invalid serialized interaction metadata:',iv,&
               this%interaction_list(iv)%currents,itmp
          stop 1
       endif
       if (.not.all(ieee_is_finite(this%interaction_list(iv)%coupl))) then
          write (*,*) 'Non-finite coupling in serialized interaction:',iv
          stop 1
       endif
       if (any(abs(this%interaction_list(iv)%coupl).gt.amplitude_value_limit)) then
          write (*,*) 'Unsafe coupling in serialized interaction:',iv,&
               this%interaction_list(iv)%coupl
          stop 1
       endif
       if (itmp.gt.0) then
          call reserve_serialized_workspace(4_8*int(itmp+1,kind=8),&
               'interaction singlet map')
          allocate(this%interaction_list(iv)%singlet_mv(0:itmp),&
               stat=allocation_status,errmsg=allocation_message)
          call require_serialized_allocation('interaction singlet map')
          this%interaction_list(iv)%singlet_mv(0)=itmp
          read (iunit,iostat=ios,iomsg=io_message) &
               this%interaction_list(iv)%singlet_mv(1:itmp)
          call require_serialized_read('interaction singlet map')
          if (any(this%interaction_list(iv)%singlet_mv(1:itmp).lt.1) .or. &
               any(this%interaction_list(iv)%singlet_mv(1:itmp).gt.n)) then
             write (*,*) 'Invalid interaction singlet map in serialized amplitude:',iv
             stop 1
          endif
       endif
    enddo
    ! momenta array
    call reserve_serialized_workspace(4_8*int(momentum_bins,kind=8)+&
         36_8*int(this%max_pp,kind=8),'serialized momentum maps')
    allocate(this%pp_bin_to_i(1:momentum_bins),this%pp_i_to_bin(1:this%max_pp),&
         this%pp(0:3,1:this%max_pp),stat=allocation_status,&
         errmsg=allocation_message)
    call require_serialized_allocation('serialized momentum maps')
    read (iunit,iostat=ios,iomsg=io_message) this%pp_bin_to_i(1:momentum_bins)
    call require_serialized_read('binary-to-momentum map')
    read (iunit,iostat=ios,iomsg=io_message) this%pp_i_to_bin(1:this%max_pp)
    call require_serialized_read('momentum-to-binary map')
    if (any(this%pp_bin_to_i.lt.0) .or. any(this%pp_bin_to_i.gt.this%max_pp) .or. &
         any(this%pp_i_to_bin.lt.1) .or. any(this%pp_i_to_bin.gt.momentum_bins)) then
       write (*,*) 'Out-of-range momentum map in serialized amplitude'
       stop 1
    endif
    do iproc=1,this%max_pp
       if (this%pp_bin_to_i(this%pp_i_to_bin(iproc)).ne.iproc) then
          write (*,*) 'Inconsistent momentum map in serialized amplitude:',iproc
          stop 1
       endif
    enddo
    this%pp=0d0
    ! process specific information
    call reserve_serialized_workspace(17_8*int(this%nprocs,kind=8)+4_8+&
         4_8*int(n,kind=8)*int(this%nprocs,kind=8),&
         'serialized process metadata')
    allocate(this%iproc_start(1:this%nprocs+1),this%same_flav(1:this%nprocs),&
         this%n_qqbar(1:this%nprocs),this%n_sing(1:this%nprocs),&
         this%processes(1:n,1:this%nprocs),stat=allocation_status,&
         errmsg=allocation_message)
    call require_serialized_allocation('serialized process metadata')
    do iproc=1,this%nprocs
       read (iunit,iostat=ios,iomsg=io_message) this%iproc_start(iproc),&
            this%same_flav(iproc),&
            this%n_qqbar(iproc),this%n_sing(iproc)
       call require_serialized_read('process counters')
       if (this%iproc_start(iproc).lt.1 .or. &
            this%iproc_start(iproc).gt.this%n_amps+1 .or. &
            this%n_qqbar(iproc).lt.0 .or. this%n_qqbar(iproc).gt.3 .or. &
            this%n_sing(iproc).lt.0 .or. this%n_sing(iproc).gt.n) then
          write (*,*) 'Invalid process counters in serialized amplitude:',iproc,&
               this%iproc_start(iproc),this%n_qqbar(iproc),this%n_sing(iproc)
          stop 1
       endif
       read (iunit,iostat=ios,iomsg=io_message) process_buffer
       call require_serialized_read('process particle identities')
       do label_index=1,n
          if (.not.valid_external_particle(process_buffer(label_index))) then
             write (*,*) 'Unsupported external particle in serialized amplitude:',&
                  iproc,label_index,process_buffer(label_index)
             stop 1
          endif
       enddo
       this%processes(1:n,iproc)=process_buffer
    enddo
    read(iunit,iostat=ios,iomsg=io_message) this%iproc_start(this%nprocs+1)
    call require_serialized_read('final subprocess offset')
    if (this%iproc_start(1).ne.1 .or. &
         this%iproc_start(this%nprocs+1).ne.this%n_amps+1 .or. &
         any(this%iproc_start(2:this%nprocs+1).lt.&
         this%iproc_start(1:this%nprocs))) then
       write (*,*) 'Invalid subprocess offsets in serialized amplitude:',this%iproc_start
       stop 1
    endif
    if (any(this%n_sing.ne.this%n_sing(1))) then
       write (*,*) 'Serialized subprocesses disagree on colour-singlet count'
       stop 1
    endif
    call validate_serialized_processes()
    ! amp specific information
    perm_extent=n-this%n_sing(1)
    call reserve_serialized_workspace(&
         33_8*int(this%n_amps,kind=8)+&
         4_8*int(n+perm_extent,kind=8)*int(this%n_amps,kind=8),&
         'serialized amplitude metadata')
    allocate(this%include_amp(1:this%n_amps),&
         this%same_flavour_sum(1:this%n_amps,1:2),&
         this%same_flavour_sum_operation(1:this%n_amps,1:2),&
         this%spins(1:n,1,1:this%n_amps),&
         this%perm(1:perm_extent,1:this%n_amps),&
         this%curr2amp(1:2,1:this%n_amps),stat=allocation_status,&
         errmsg=allocation_message)
    call require_serialized_allocation('serialized amplitude metadata')
    this%curr2amp=0
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          read (iunit,iostat=ios,iomsg=io_message) this%include_amp(iamp),&
               this%same_flavour_sum(iamp,1),this%same_flavour_sum(iamp,2),&
               this%same_flavour_sum_operation(iamp,1),&
               this%same_flavour_sum_operation(iamp,2)
          call require_serialized_read('amplitude inclusion metadata')
          read (iunit,iostat=ios,iomsg=io_message) spin_buffer
          call require_serialized_read('amplitude helicities')
          this%spins(1:n,1,iamp)=spin_buffer
          if (any(this%spins(1:n,1,iamp).lt.-1) .or. &
               any(this%spins(1:n,1,iamp).gt.1)) then
             write (*,*) 'Invalid helicity in serialized amplitude:',iamp
             stop 1
          endif
          read (iunit,iostat=ios,iomsg=io_message) &
               permutation_buffer(1:perm_extent)
          call require_serialized_read('amplitude permutation')
          this%perm(1:perm_extent,iamp)=permutation_buffer(1:perm_extent)
          ! Use local contiguous buffers here as well: bounds-enabled gfortran
          ! mishandles rank-reduced allocatable components of a polymorphic
          ! object in some assumed-shape argument descriptors.
          process_buffer=this%processes(1:n,iproc)
          if (.not.valid_external_colour_order(&
               permutation_buffer(1:perm_extent),process_buffer)) then
             write (*,*) 'Invalid colour permutation in serialized amplitude:',iamp,&
                  permutation_buffer(1:perm_extent)
             stop 1
          endif
          if (.not.this%same_flav(iproc)) then
             read (iunit,iostat=ios,iomsg=io_message) &
                  this%curr2amp(1,iamp),this%curr2amp(2,iamp)
             call require_serialized_read('current-to-amplitude map')
             if (any(this%curr2amp(1:2,iamp).lt.1) .or. &
                  any(this%curr2amp(1:2,iamp).gt.this%n_cur)) then
                write (*,*) 'Invalid current-to-amplitude map in serialized amplitude:',iamp
                stop 1
             endif
          endif
       enddo
    enddo
    if (any(this%same_flavour_sum.lt.-1) .or. &
         any(this%same_flavour_sum.gt.this%n_amps)) then
       write (*,*) 'Invalid same-flavour map in serialized amplitude'
       stop 1
    endif
    do iamp=1,this%n_amps
       if ((this%same_flavour_sum(iamp,1).gt.0) .neqv. &
            (this%same_flavour_sum(iamp,2).gt.0)) then
          write (*,*) 'Incomplete same-flavour map in serialized amplitude:',iamp,&
               this%same_flavour_sum(iamp,1:2)
          stop 1
       endif
    enddo
    where (this%same_flavour_sum.le.0)
       this%same_flavour_sum_operation=0
    endwhere
    if (any((this%same_flavour_sum.gt.0) .and. &
         (this%same_flavour_sum_operation.lt.0 .or. &
         this%same_flavour_sum_operation.gt.7))) then
       write (*,*) 'Invalid same-flavour operation in serialized amplitude'
       stop 1
    endif
  contains
    subroutine validate_serialized_processes()
      implicit none
      integer :: process_index,particle_index,n_quark,n_antiquark,n_singlet
      logical :: first_is_lepton,current_is_lepton

      do process_index=1,this%nprocs
         n_quark=0
         n_antiquark=0
         n_singlet=0
         do particle_index=1,n
            if (particle_index.le.2) then
               if (this%processes(particle_index,process_index).ge.1 .and. &
                    this%processes(particle_index,process_index).le.6) &
                    n_antiquark=n_antiquark+1
               if (this%processes(particle_index,process_index).le.-1 .and. &
                    this%processes(particle_index,process_index).ge.-6) &
                    n_quark=n_quark+1
            else
               if (this%processes(particle_index,process_index).ge.1 .and. &
                    this%processes(particle_index,process_index).le.6) &
                    n_quark=n_quark+1
               if (this%processes(particle_index,process_index).le.-1 .and. &
                    this%processes(particle_index,process_index).ge.-6) &
                    n_antiquark=n_antiquark+1
            endif
            if (.not.is_coloured_external_particle(&
                 this%processes(particle_index,process_index))) &
                 n_singlet=n_singlet+1
            if (process_index.gt.1) then
               first_is_lepton=&
                    (this%processes(particle_index,1).ge.11 .and. &
                    this%processes(particle_index,1).le.16) .or. &
                    (this%processes(particle_index,1).le.-11 .and. &
                    this%processes(particle_index,1).ge.-16)
               current_is_lepton=&
                    (this%processes(particle_index,process_index).ge.11 .and. &
                    this%processes(particle_index,process_index).le.16) .or. &
                    (this%processes(particle_index,process_index).le.-11 .and. &
                    this%processes(particle_index,process_index).ge.-16)
               if (first_is_lepton.neqv.current_is_lepton) then
                  write (*,*) 'Serialized subprocesses disagree on lepton positions:',&
                       particle_index,process_index
                  stop 1
               endif
            endif
         enddo
         if (n_quark.ne.n_antiquark .or. &
              n_quark.ne.this%n_qqbar(process_index)) then
            write (*,*) 'Inconsistent quark-line count in serialized amplitude:',&
                 process_index,n_quark,n_antiquark,this%n_qqbar(process_index)
            stop 1
         endif
         if (n_singlet.ne.this%n_sing(process_index)) then
            write (*,*) 'Inconsistent singlet count in serialized amplitude:',&
                 process_index,n_singlet,this%n_sing(process_index)
            stop 1
         endif
      enddo
    end subroutine validate_serialized_processes

    subroutine require_serialized_read(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Truncated/corrupt serialized amplitude while reading ',&
              trim(label),ios,trim(io_message)
         stop 1
      endif
    end subroutine require_serialized_read

    subroutine require_serialized_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate ',trim(label),': ',trim(allocation_message)
         stop 1
      endif
    end subroutine require_serialized_allocation

    subroutine reserve_serialized_workspace(bytes,label)
      implicit none
      integer(kind=8),intent(in) :: bytes
      character(len=*),intent(in) :: label
      if (bytes.lt.0_8 .or. bytes.gt.max_serialized_amplitude_workspace_bytes .or. &
           workspace_bytes.gt.max_serialized_amplitude_workspace_bytes-bytes) then
         write (*,*) 'Serialized amplitude exceeds the supported workspace for ',&
              trim(label),workspace_bytes,bytes,max_serialized_amplitude_workspace_bytes
         stop 1
      endif
      workspace_bytes=workspace_bytes+bytes
    end subroutine reserve_serialized_workspace

    subroutine validate_serialized_ranges()
      implicit none
      integer :: i
      if (this%n_cur_start(1).ne.1 .or. this%n_cur_end(n-1).ne.this%n_cur) then
         write (*,*) 'Serialized amplitude has incomplete current ranges'
         stop 1
      endif
      do i=1,n-1
         if (this%n_cur_start(i).lt.1 .or. this%n_cur_start(i).gt.this%n_cur+1 .or. &
              this%n_cur_end(i).lt.0 .or. this%n_cur_end(i).gt.this%n_cur .or. &
              this%n_cur_start(i).gt.this%n_cur_end(i)+1) then
            write (*,*) 'Invalid serialized current range:',i,&
                 this%n_cur_start(i),this%n_cur_end(i)
            stop 1
         endif
         if (i.gt.1) then
            if (this%n_cur_start(i).ne.this%n_cur_end(i-1)+1) then
               write (*,*) 'Non-contiguous serialized current ranges:',i
               stop 1
            endif
         endif
      enddo
      if (this%n_cur_start(n).lt.this%n_cur_start(1) .or. &
           this%n_cur_end(n).gt.this%n_cur_end(1) .or. &
           this%n_cur_start(n).gt.this%n_cur_end(n)+1) then
         write (*,*) 'Invalid closing-current range in serialized amplitude'
         stop 1
      endif
      if (this%n_vert.eq.0) then
         if (any(this%n_vert_start(2:n-1).ne.1) .or. &
              any(this%n_vert_end(2:n-1).ne.0)) then
            write (*,*) 'Invalid empty interaction ranges in serialized amplitude'
            stop 1
         endif
      else
         if (this%n_vert_start(2).ne.1 .or. &
              this%n_vert_end(n-1).ne.this%n_vert) then
            write (*,*) 'Serialized amplitude has incomplete interaction ranges'
            stop 1
         endif
         do i=2,n-1
            if (this%n_vert_start(i).lt.1 .or. &
                 this%n_vert_start(i).gt.this%n_vert+1 .or. &
                 this%n_vert_end(i).lt.0 .or. this%n_vert_end(i).gt.this%n_vert .or. &
                 this%n_vert_start(i).gt.this%n_vert_end(i)+1) then
               write (*,*) 'Invalid serialized interaction range:',i,&
                    this%n_vert_start(i),this%n_vert_end(i)
               stop 1
            endif
            if (i.gt.2) then
               if (this%n_vert_start(i).ne.this%n_vert_end(i-1)+1) then
                  write (*,*) 'Non-contiguous serialized interaction ranges:',i
                  stop 1
               endif
            endif
         enddo
      endif
    end subroutine validate_serialized_ranges

  end subroutine read_init_amps_from_file
  
  subroutine evaluate(this,n,p,hel,read_file,pm,status)
    use FeynmanRules
    use particles
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    type(physics_model),intent(in) :: pm
    integer,intent(in) :: n
    integer,dimension(n),intent(in) :: hel
    real(kind=8),dimension(0:3,n),intent(in) :: p
    integer :: ic,iv,isize,ih_in,ifinal,dim,allocation_status
    integer(kind=8) :: evaluation_workspace_bytes
    real(kind=8),dimension(0:3) :: pp_loc,pp1_loc,pp2_loc
    complex(kind=8),dimension(6) :: wf1_loc,wf2_loc
    character(len=256) :: allocation_message
    logical :: use_real_amplitudes
    logical,intent(in) :: read_file
    integer,intent(out),optional :: status
    call validate_evaluation_header()
    call ensure_evaluation_workspace()
    ! This argument is retained for source compatibility.  File-backed and
    ! freshly constructed amplitudes deliberately use the same evaluator.
    if (read_file) continue
    if (present(status)) status=0
    call reset_feynman_numerical_status()
    call validate_feynman_momenta(p)
    if (.not.feynman_numerical_status_ok()) then
       if (allocated(this%amps)) this%amps=(0d0,0d0)
       if (allocated(this%amps_r)) this%amps_r=0d0
       if (present(status)) status=-20
       return
    endif

    call fill_momentum_array()
    call validate_feynman_momenta(this%pp)
    if (.not.feynman_numerical_status_ok()) then
       if (allocated(this%amps)) this%amps=(0d0,0d0)
       if (allocated(this%amps_r)) this%amps_r=0d0
       if (present(status)) status=-20
       return
    endif

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
             pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
             if (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. pm%is_photon(this%current_list(ic)%type)) then
                if (use_real_gluons) then
                   call ext_massless_vector_real(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_r(1:4))
                else
                   call ext_massless_vector_cmplx(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4))
                endif
             elseif (pm%is_quark(this%current_list(ic)%type) .or. &
                  pm%is_lepton(this%current_list(ic)%type)) then
                if (this%current_list(ic)%chirality.ne.0) then
                   call ext_fermion_outflow_weyl(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:2),this%current_list(ic)%chirality)
                else
                   call ext_fermion_outflow(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
                endif
             elseif (pm%is_antiquark(this%current_list(ic)%type) .or. &
                  pm%is_antilepton(this%current_list(ic)%type)) then
                if (this%current_list(ic)%chirality.ne.0) then
                   call ext_fermion_inflow_weyl(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:2),this%current_list(ic)%chirality)
                else
                   call ext_fermion_inflow(pp_loc, &
                        ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
                endif
             elseif (pm%is_massive_vector(this%current_list(ic)%type)) then
                call ext_massive_vector(pp_loc, &
                     ih_in,ifinal,this%current_list(ic)%val_c(1:4),this%current_list(ic)%mass)
             elseif (pm%is_higgs(this%current_list(ic)%type)) then
                call ext_scalar(pp_loc, &
                     ifinal,this%current_list(ic)%val_c(1))
             else
                write (*,*) 'External particle type unknown',ic,this%current_list(ic)%type,ih_in
                stop 1
             endif
          enddo
          if (.not.feynman_numerical_status_ok()) exit
          cycle
       endif
       ! loop over the vertices required to create all the currents with isize
       ! number of external particles combined
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          if (this%interaction_list(iv)%type.eq.0) then
             pp1_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin))
             pp2_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin))
             wf1_loc=(0d0,0d0)
             wf2_loc=(0d0,0d0)
             wf1_loc(1:size(this%current_list(this%interaction_list(iv)%currents(1))%val_c))=&
                  this%current_list(this%interaction_list(iv)%currents(1))%val_c
             wf2_loc(1:size(this%current_list(this%interaction_list(iv)%currents(2))%val_c))=&
                  this%current_list(this%interaction_list(iv)%currents(2))%val_c
             if (use_real_gluons) then
                call threeGluon_real(this%current_list(this%interaction_list(iv)%currents(1))%val_r(1:4),&
                     pp1_loc,&
                     this%current_list(this%interaction_list(iv)%currents(2))%val_r(1:4),&
                     pp2_loc,&
                     this%interaction_list(iv)%val_r(1:4))
             else
                call threeGluon(wf1_loc(1:4),&
                     pp1_loc,&
                     wf2_loc(1:4),&
                     pp2_loc,&
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
             pp1_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(1))%bin))
             pp2_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(this%interaction_list(iv)%currents(2))%bin))
             call VectorVectorToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1:4),&
                       pp1_loc,&
                       this%current_list(this%interaction_list(iv)%currents(2))%val_c(1:4),&
                       pp2_loc,&
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
                call LeptonAntileptonToVector_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                               this%interaction_list(iv)%val_c(1:4),&
                                               this%interaction_list(iv)%coupl(1:2),&
                                               this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                               this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call LeptonAntileptonToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                          this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                          this%interaction_list(iv)%val_c(1:4),&
                                          this%interaction_list(iv)%coupl(1:2))
             endif

          elseif(this%interaction_list(iv)%type.eq.22) then
             if (this%current_list(this%interaction_list(iv)%currents(1))%chirality.ne.0 .or. &
                 this%current_list(this%interaction_list(iv)%currents(2))%chirality.ne.0) then
                call AntileptonLeptonToVector_weyl(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
                                               this%current_list(this%interaction_list(iv)%currents(2))%val_c(1),&
                                               this%interaction_list(iv)%val_c(1:4),&
                                               this%interaction_list(iv)%coupl(1:2),&
                                               this%current_list(this%interaction_list(iv)%currents(1))%chirality,&
                                               this%current_list(this%interaction_list(iv)%currents(2))%chirality)
             else
                call  AntileptonLeptonToVector(this%current_list(this%interaction_list(iv)%currents(1))%val_c(1),&
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
          if (allocated(this%interaction_list(iv)%val_c)) &
               call validate_complex_wavefunction(this%interaction_list(iv)%val_c)
          if (allocated(this%interaction_list(iv)%val_r)) &
               call validate_real_wavefunction(this%interaction_list(iv)%val_r)
          if (.not.feynman_numerical_status_ok()) exit
       enddo
       if (.not.feynman_numerical_status_ok()) exit

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
          if (.not.feynman_numerical_status_ok()) exit
       enddo
       if (.not.feynman_numerical_status_ok()) exit
    enddo

    if (feynman_numerical_status_ok()) call compute_amps_from_currents
    if (allocated(this%amps)) then
       if (.not.all(complex_amplitude_is_safe(this%amps))) &
            call invalidate_feynman_point()
    endif
    if (allocated(this%amps_r)) then
       if (.not.all(ieee_is_finite(this%amps_r))) then
          call invalidate_feynman_point()
       else if (any(abs(this%amps_r).gt.amplitude_value_limit)) then
          call invalidate_feynman_point()
       endif
    endif
    if (.not.feynman_numerical_status_ok()) then
       if (allocated(this%amps)) this%amps=(0d0,0d0)
       if (allocated(this%amps_r)) this%amps_r=0d0
       if (present(status)) status=-20
    endif

  contains

    subroutine validate_evaluation_header()
      implicit none
      if (n.lt.3 .or. n.gt.max_amplitude_external_particles) then
         write (*,*) 'Invalid external multiplicity for amplitude evaluation:',n
         stop 1
      endif
      if (any(hel.lt.-1) .or. any(hel.gt.1)) then
         write (*,*) 'Invalid helicity for amplitude evaluation:',hel
         stop 1
      endif
      if (this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
           this%n_vert.lt.0 .or. this%n_vert.gt.max_amplitude_interaction_records .or. &
           this%n_amps.lt.1 .or. this%n_amps.gt.max_amplitude_current_records .or. &
           this%nprocs.lt.1 .or. this%imode.lt.1 .or. this%imode.gt.3 .or. &
           this%max_pp.lt.1) then
         write (*,*) 'Amplitude is not initialized for evaluation:',this%n_cur,&
              this%n_vert,this%n_amps,this%nprocs,this%imode,this%max_pp
         stop 1
      endif
      if (.not.allocated(this%current_list)) &
           call invalid_evaluation_state('missing current list')
      if (size(this%current_list).lt.this%n_cur) &
           call invalid_evaluation_state('inconsistent current-list extent')
      if (.not.allocated(this%interaction_list)) &
           call invalid_evaluation_state('missing interaction list')
      if (size(this%interaction_list).lt.this%n_vert) &
           call invalid_evaluation_state('inconsistent interaction-list extent')
      if (.not.allocated(this%processes)) &
           call invalid_evaluation_state('missing subprocess metadata')
      if (size(this%processes,1).ne.n .or. &
           size(this%processes,2).ne.this%nprocs) &
           call invalid_evaluation_state('incompatible subprocess dimensions')
      if (.not.allocated(this%n_cur_start) .or. &
           .not.allocated(this%n_cur_end)) &
           call invalid_evaluation_state('missing current ranges')
      if (size(this%n_cur_start).ne.n .or. size(this%n_cur_end).ne.n) &
           call invalid_evaluation_state('incompatible current ranges')
      if (.not.allocated(this%n_vert_start) .or. &
           .not.allocated(this%n_vert_end)) &
           call invalid_evaluation_state('missing interaction ranges')
      if (size(this%n_vert_start).ne.n-2 .or. &
           size(this%n_vert_end).ne.n-2) &
           call invalid_evaluation_state('incompatible interaction ranges')
      if (.not.allocated(this%pp) .or. &
           .not.allocated(this%pp_bin_to_i) .or. &
           .not.allocated(this%pp_i_to_bin)) &
           call invalid_evaluation_state('missing momentum workspace')
      if (size(this%pp,1).ne.4 .or. size(this%pp,2).ne.this%max_pp .or. &
           size(this%pp_i_to_bin).ne.this%max_pp .or. &
           size(this%pp_bin_to_i).ne.maskr(n)) &
           call invalid_evaluation_state('incompatible momentum workspace')
      if (.not.allocated(this%curr2amp)) &
           call invalid_evaluation_state('missing current-to-amplitude map')
      if (size(this%curr2amp,1).ne.2 .or. &
           size(this%curr2amp,2).lt.this%n_amps) &
           call invalid_evaluation_state('incompatible current-to-amplitude map')
      if (.not.allocated(this%iproc_start) .or. &
           .not.allocated(this%same_flav) .or. &
           .not.allocated(this%n_qqbar)) &
           call invalid_evaluation_state('missing subprocess evaluation metadata')
      if (size(this%iproc_start).ne.this%nprocs+1 .or. &
           size(this%same_flav).ne.this%nprocs .or. &
           size(this%n_qqbar).ne.this%nprocs) &
           call invalid_evaluation_state('incompatible subprocess evaluation metadata')
      if (.not.allocated(this%same_flavour_sum) .or. &
           .not.allocated(this%same_flavour_sum_operation)) &
           call invalid_evaluation_state('missing same-flavour amplitude metadata')
      if (size(this%same_flavour_sum,1).lt.this%n_amps .or. &
           size(this%same_flavour_sum,2).ne.2 .or. &
           size(this%same_flavour_sum_operation,1).lt.this%n_amps .or. &
           size(this%same_flavour_sum_operation,2).ne.2) &
           call invalid_evaluation_state('incompatible same-flavour amplitude metadata')
    end subroutine validate_evaluation_header

    subroutine ensure_evaluation_workspace()
      implicit none

      use_real_amplitudes=.false.
      if (use_real_gluons) then
         if (this%imode.eq.2) then
            use_real_amplitudes=this%n_qqbar(1).eq.0
         else
            use_real_amplitudes=all(this%n_qqbar.eq.0)
         endif
      endif
      if (this%evaluation_workspace_ready) then
         if (use_real_amplitudes) then
            if (.not.allocated(this%amps_r)) &
                 call invalid_evaluation_state('missing real amplitude workspace')
            if (size(this%amps_r).ne.this%n_amps) &
                 call invalid_evaluation_state('incompatible real amplitude workspace')
         else
            if (.not.allocated(this%amps)) &
                 call invalid_evaluation_state('missing complex amplitude workspace')
            if (size(this%amps).ne.this%n_amps) &
                 call invalid_evaluation_state('incompatible complex amplitude workspace')
         endif
         return
      endif

      call clear_evaluation_workspace()
      call validate_evaluation_metadata()
      evaluation_workspace_bytes=0_8
      do ic=1,this%n_cur
         if (this%current_list(ic)%chirality.lt.-1 .or. &
              this%current_list(ic)%chirality.gt.1) then
            write (*,*) 'Invalid current chirality for evaluation:',ic,&
                 this%current_list(ic)%chirality
            stop 1
         endif
         dim=current_dim(ic)
         call validate_wavefunction_dimension('current',ic,dim)
         if (use_real_gluons .and. &
              (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. &
              pm%is_auxiliary_tensor(this%current_list(ic)%type))) then
            call reserve_evaluation_workspace(8_8*int(dim,kind=8),'current values')
         else
            call reserve_evaluation_workspace(16_8*int(dim,kind=8),'current values')
         endif
      enddo
      do iv=1,this%n_vert
         if (this%interaction_list(iv)%type.lt.0 .or. &
              this%interaction_list(iv)%type.gt.24 .or. &
              this%interaction_list(iv)%chirality.lt.-1 .or. &
              this%interaction_list(iv)%chirality.gt.1 .or. &
              any(this%interaction_list(iv)%currents.lt.1) .or. &
              any(this%interaction_list(iv)%currents.gt.this%n_cur)) then
            write (*,*) 'Invalid interaction metadata for evaluation:',iv,&
                 this%interaction_list(iv)%type,this%interaction_list(iv)%chirality,&
                 this%interaction_list(iv)%currents
            stop 1
         endif
         dim=interaction_dim(iv)
         call validate_wavefunction_dimension('interaction',iv,dim)
         if (use_real_gluons .and. this%interaction_list(iv)%type.ge.0 .and. &
              this%interaction_list(iv)%type.le.3) then
            call reserve_evaluation_workspace(8_8*int(dim,kind=8),'interaction values')
         else
            call reserve_evaluation_workspace(16_8*int(dim,kind=8),'interaction values')
         endif
      enddo
      if (any(this%curr2amp(:,1:this%n_amps).lt.0) .or. &
           any(this%curr2amp(:,1:this%n_amps).gt.this%n_cur)) &
           call invalid_evaluation_state('out-of-range current-to-amplitude map')
      if (any(this%same_flavour_sum(1:this%n_amps,:).lt.-1) .or. &
           any(this%same_flavour_sum(1:this%n_amps,:).gt.this%n_amps)) &
           call invalid_evaluation_state('out-of-range same-flavour amplitude map')
      if (use_real_amplitudes) then
         call reserve_evaluation_workspace(8_8*int(this%n_amps,kind=8),&
              'real amplitudes')
      else
         call reserve_evaluation_workspace(16_8*int(this%n_amps,kind=8),&
              'complex amplitudes')
      endif

      do ic=1,this%n_cur
         dim=current_dim(ic)
         if (use_real_gluons .and. &
              (pm%is_colour_flow_vector(this%current_list(ic)%type) .or. &
              pm%is_auxiliary_tensor(this%current_list(ic)%type))) then
            allocate(this%current_list(ic)%val_r(1:dim),stat=allocation_status,&
                 errmsg=allocation_message)
         else
            allocate(this%current_list(ic)%val_c(1:dim),stat=allocation_status,&
                 errmsg=allocation_message)
         endif
         call require_evaluation_allocation('current values')
      enddo
      do iv=1,this%n_vert
         dim=interaction_dim(iv)
         if (use_real_gluons .and. this%interaction_list(iv)%type.ge.0 .and. &
              this%interaction_list(iv)%type.le.3) then
            allocate(this%interaction_list(iv)%val_r(1:dim),stat=allocation_status,&
                 errmsg=allocation_message)
         else
            allocate(this%interaction_list(iv)%val_c(1:dim),stat=allocation_status,&
                 errmsg=allocation_message)
         endif
         call require_evaluation_allocation('interaction values')
      enddo
      if (use_real_amplitudes) then
         allocate(this%amps_r(1:this%n_amps),stat=allocation_status,&
              errmsg=allocation_message)
         call require_evaluation_allocation('real amplitudes')
         this%amps_r=0d0
      else
         allocate(this%amps(1:this%n_amps),stat=allocation_status,&
              errmsg=allocation_message)
         call require_evaluation_allocation('complex amplitudes')
         this%amps=(0d0,0d0)
      endif
      this%evaluation_workspace_ready=.true.
    end subroutine ensure_evaluation_workspace

    subroutine validate_evaluation_metadata()
      implicit none
      integer :: i,iproc,iamp

      if (this%n_cur_start(1).ne.1 .or. &
           this%n_cur_end(n-1).ne.this%n_cur) &
           call invalid_evaluation_state('incomplete current ranges')
      do i=1,n-1
         if (this%n_cur_start(i).lt.1 .or. &
              this%n_cur_start(i).gt.this%n_cur+1 .or. &
              this%n_cur_end(i).lt.0 .or. &
              this%n_cur_end(i).gt.this%n_cur .or. &
              this%n_cur_start(i).gt.this%n_cur_end(i)+1) &
              call invalid_evaluation_state('out-of-range current range')
         if (i.gt.1) then
            if (this%n_cur_start(i).ne.this%n_cur_end(i-1)+1) &
                 call invalid_evaluation_state('non-contiguous current ranges')
         endif
      enddo
      if (this%n_cur_start(n).lt.this%n_cur_start(1) .or. &
           this%n_cur_end(n).gt.this%n_cur_end(1) .or. &
           this%n_cur_start(n).gt.this%n_cur_end(n)+1) &
           call invalid_evaluation_state('invalid closing-current range')
      if (this%n_vert.eq.0) then
         if (any(this%n_vert_start.ne.1) .or. any(this%n_vert_end.ne.0)) &
              call invalid_evaluation_state('invalid empty interaction ranges')
      else
         if (this%n_vert_start(2).ne.1 .or. &
              this%n_vert_end(n-1).ne.this%n_vert) &
              call invalid_evaluation_state('incomplete interaction ranges')
         do i=2,n-1
            if (this%n_vert_start(i).lt.1 .or. &
                 this%n_vert_start(i).gt.this%n_vert+1 .or. &
                 this%n_vert_end(i).lt.0 .or. &
                 this%n_vert_end(i).gt.this%n_vert .or. &
                 this%n_vert_start(i).gt.this%n_vert_end(i)+1) &
                 call invalid_evaluation_state('out-of-range interaction range')
            if (i.gt.2) then
               if (this%n_vert_start(i).ne.this%n_vert_end(i-1)+1) &
                    call invalid_evaluation_state('non-contiguous interaction ranges')
            endif
         enddo
      endif

      do i=1,this%max_pp
         if (this%pp_i_to_bin(i).lt.1 .or. &
              this%pp_i_to_bin(i).gt.size(this%pp_bin_to_i)) &
              call invalid_evaluation_state('out-of-range momentum reverse map')
         if (this%pp_bin_to_i(this%pp_i_to_bin(i)).ne.i) &
              call invalid_evaluation_state('inconsistent momentum maps')
      enddo
      do i=1,this%n_cur
         if (this%current_list(i)%bin.lt.1 .or. &
              this%current_list(i)%bin.gt.size(this%pp_bin_to_i)) &
              call invalid_evaluation_state('out-of-range current momentum label')
         if (this%pp_bin_to_i(this%current_list(i)%bin).lt.1 .or. &
              this%pp_bin_to_i(this%current_list(i)%bin).gt.this%max_pp) &
              call invalid_evaluation_state('unmapped current momentum label')
         if (this%current_list(i)%n_vert.lt.0 .or. &
              this%current_list(i)%n_vert.gt.this%n_vert) &
              call invalid_evaluation_state('invalid current interaction count')
         if (this%current_list(i)%n_vert.gt.0) then
            if (.not.allocated(this%current_list(i)%vertices) .or. &
                 .not.allocated(this%current_list(i)%vertex_sign)) &
                 call invalid_evaluation_state('missing current interaction list')
            if (size(this%current_list(i)%vertices).lt.&
                 this%current_list(i)%n_vert .or. &
                 size(this%current_list(i)%vertex_sign).lt.&
                 this%current_list(i)%n_vert) &
                 call invalid_evaluation_state('short current interaction list')
            if (any(this%current_list(i)%vertices(1:&
                 this%current_list(i)%n_vert).lt.1) .or. &
                 any(this%current_list(i)%vertices(1:&
                 this%current_list(i)%n_vert).gt.this%n_vert)) &
                 call invalid_evaluation_state('out-of-range current interaction')
         endif
      enddo
      do i=this%n_cur_start(1),this%n_cur_end(1)
         if (.not.allocated(this%current_list(i)%order) .or. &
              .not.allocated(this%current_list(i)%spin)) &
              call invalid_evaluation_state('missing external-current labels')
         if (size(this%current_list(i)%order).lt.1 .or. &
              size(this%current_list(i)%spin).lt.1) &
              call invalid_evaluation_state('empty external-current labels')
         if (this%current_list(i)%order(1).lt.1 .or. &
              this%current_list(i)%order(1).gt.n .or. &
              this%current_list(i)%spin(1).lt.-9 .or. &
              this%current_list(i)%spin(1).gt.1 .or. &
              (this%current_list(i)%spin(1).lt.-1 .and. &
              this%current_list(i)%spin(1).ne.-9)) &
              call invalid_evaluation_state('invalid external-current labels')
      enddo

      if (this%iproc_start(1).ne.1 .or. &
           this%iproc_start(this%nprocs+1).ne.this%n_amps+1) &
           call invalid_evaluation_state('incomplete subprocess offsets')
      if (any(this%iproc_start(2:this%nprocs+1).lt.&
           this%iproc_start(1:this%nprocs))) &
           call invalid_evaluation_state('non-monotonic subprocess offsets')
      do iproc=1,this%nprocs
         if (this%iproc_start(iproc).lt.1 .or. &
              this%iproc_start(iproc+1).gt.this%n_amps+1) &
              call invalid_evaluation_state('out-of-range subprocess offsets')
         if (this%same_flav(iproc)) cycle
         do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
            if (any(this%curr2amp(:,iamp).lt.1)) &
                 call invalid_evaluation_state('missing current-to-amplitude entry')
         enddo
      enddo
      if (any(this%same_flavour_sum_operation(1:this%n_amps,:).lt.0) .or. &
           any(this%same_flavour_sum_operation(1:this%n_amps,:).gt.7)) &
           call invalid_evaluation_state('invalid same-flavour operation')
    end subroutine validate_evaluation_metadata

    subroutine require_evaluation_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (allocation_status.ne.0) then
         call clear_evaluation_workspace()
         write (*,*) 'Could not allocate amplitude-evaluation ',trim(label),': ',&
              trim(allocation_message)
         stop 1
      endif
    end subroutine require_evaluation_allocation

    subroutine reserve_evaluation_workspace(bytes,label)
      implicit none
      integer(kind=8),intent(in) :: bytes
      character(len=*),intent(in) :: label
      if (bytes.lt.0_8 .or. bytes.gt.max_amplitude_workspace_bytes .or. &
           evaluation_workspace_bytes.gt.max_amplitude_workspace_bytes-bytes) then
         write (*,*) 'Amplitude evaluation exceeds the supported workspace for ',&
              trim(label),evaluation_workspace_bytes,bytes,&
              max_amplitude_workspace_bytes
         stop 1
      endif
      evaluation_workspace_bytes=evaluation_workspace_bytes+bytes
    end subroutine reserve_evaluation_workspace

    subroutine clear_evaluation_workspace()
      implicit none
      integer :: i
      if (allocated(this%current_list)) then
         do i=1,size(this%current_list)
            if (allocated(this%current_list(i)%val_c)) &
                 deallocate(this%current_list(i)%val_c)
            if (allocated(this%current_list(i)%val_r)) &
                 deallocate(this%current_list(i)%val_r)
         enddo
      endif
      if (allocated(this%interaction_list)) then
         do i=1,size(this%interaction_list)
            if (allocated(this%interaction_list(i)%val_c)) &
                 deallocate(this%interaction_list(i)%val_c)
            if (allocated(this%interaction_list(i)%val_r)) &
                 deallocate(this%interaction_list(i)%val_r)
         enddo
      endif
      if (allocated(this%amps)) deallocate(this%amps)
      if (allocated(this%amps_r)) deallocate(this%amps_r)
      this%evaluation_workspace_ready=.false.
    end subroutine clear_evaluation_workspace

    subroutine validate_wavefunction_dimension(label,index,dimension)
      implicit none
      character(len=*),intent(in) :: label
      integer,intent(in) :: index,dimension
      if (dimension.lt.1 .or. dimension.gt.6) then
         write (*,*) 'Invalid ',trim(label),' wavefunction dimension:',index,dimension
         stop 1
      endif
    end subroutine validate_wavefunction_dimension

    subroutine invalid_evaluation_state(message)
      implicit none
      character(len=*),intent(in) :: message
      write (*,*) 'Invalid amplitude evaluation state: ',trim(message)
      stop 1
    end subroutine invalid_evaluation_state

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
      integer :: iamp,iproc,idau
      if (this%imode.eq.1 .or. this%imode.eq.3) then
         if (use_real_gluons .and. all(this%n_qqbar(1:this%nprocs).eq.0)) then
            do iamp=1,this%n_amps
               this%amps_r(iamp)=sum(this%current_list(this%curr2amp(1,iamp))%val_r(1:4)* &
                                     this%current_list(this%curr2amp(2,iamp))%val_r(1:4))
            enddo
         else
            do iproc=1,this%nprocs
               do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (.not.this%same_flav(iproc)) then
                     this%amps(iamp)=contract_currents(iamp)
                  endif
               enddo
            enddo
            do iproc=1,this%nprocs
               do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (this%same_flav(iproc)) then
                     ! same-flavour amps are build from two different-flavour amps
                     this%amps(iamp)=(0d0,0d0)
                     do idau=1,2
                        if (this%same_flavour_sum(iamp,idau).gt.0) then
                           this%amps(iamp)=this%amps(iamp)+apply_operation(iamp,idau)
                        endif
                     enddo
                  endif
               enddo
            enddo
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
                     this%amps(iamp)=-contract_currents(iamp)
                  else
                     if (.not.this%same_flav(iproc)) then
                        this%amps(iamp)=contract_currents(iamp)
                        ! Fortran does not require short-circuit evaluation of
                        ! .and.; keep the allocation test separate from the
                        ! array reference (zero-sized partner maps are valid).
                        if (allocated(this%three_line_partner_curr2amp)) then
                           if (iamp.le.size(this%three_line_partner_curr2amp,2)) then
                              if (this%three_line_partner_curr2amp(1,iamp).ne.0) then
                                 this%amps(iamp)=this%amps(iamp)+contract_current_pair(&
                                      this%three_line_partner_curr2amp(1,iamp),&
                                      this%three_line_partner_curr2amp(2,iamp))
                              endif
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
            do iproc=1,this%nprocs
               do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
                  if (use_symmetry .and. this%n_qqbar(1).eq.0 .and. iamp.gt.this%n_amps/2 .and. mod(n,2).eq.1) then
                     cycle
                  else
                     if (this%same_flav(iproc)) then
                        ! same-flavour amps are build from two different-flavour amps
                        this%amps(iamp)=(0d0,0d0)
                        do idau=1,2
                           if (this%same_flavour_sum(iamp,idau).gt.0) then
                              this%amps(iamp)=this%amps(iamp)+apply_operation(iamp,idau)
                           endif
                        enddo
                     endif
                  endif
               enddo
            enddo
         endif
      endif
    end subroutine compute_amps_from_currents

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
      if (btest(this%same_flavour_sum_operation(iamp,idau),2)) &
           apply_operation=cmplx(aimag(apply_operation),dble(apply_operation),kind=8)
      if (btest(this%same_flavour_sum_operation(iamp,idau),0)) apply_operation=-conjg(apply_operation)
      if (btest(this%same_flavour_sum_operation(iamp,idau),1)) apply_operation=conjg(apply_operation)
    end function apply_operation
    
    subroutine combine_interactions(dim)
      implicit none
      integer :: dim,iv
      if (use_real_gluons .and. &
           (pm%is_colour_flow_vector(this%current_list(ic)%type).or. &
           pm%is_gluon_aux_tensor(this%current_list(ic)%type))) then
         this%current_list(ic)%val_r(1:dim)=0d0
         do iv=1,this%current_list(ic)%n_vert
            call validate_real_wavefunction(&
                 this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim))
            if (.not.feynman_numerical_status_ok()) exit
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_r(1:dim)=&
                    this%current_list(ic)%val_r(1:dim)-this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim)
            else
               this%current_list(ic)%val_r(1:dim)=&
                    this%current_list(ic)%val_r(1:dim)+this%interaction_list(this%current_list(ic)%vertices(iv))%val_r(1:dim)
            endif
            call validate_real_wavefunction(this%current_list(ic)%val_r(1:dim))
            if (.not.feynman_numerical_status_ok()) exit
         enddo

      else
         this%current_list(ic)%val_c(1:dim)=(0d0,0d0)
         do iv=1,this%current_list(ic)%n_vert
            call validate_complex_wavefunction(&
                 this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim))
            if (.not.feynman_numerical_status_ok()) exit
            if (this%current_list(ic)%vertex_sign(iv))then
               this%current_list(ic)%val_c(1:dim)=&
                    this%current_list(ic)%val_c(1:dim)-this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim)
            else
               this%current_list(ic)%val_c(1:dim)=&
                    this%current_list(ic)%val_c(1:dim)+this%interaction_list(this%current_list(ic)%vertices(iv))%val_c(1:dim)
            endif
            call validate_complex_wavefunction(this%current_list(ic)%val_c(1:dim))
            if (.not.feynman_numerical_status_ok()) exit
         enddo
      endif
    end subroutine combine_interactions
    subroutine include_massless_vector_propagator()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      if (use_real_gluons) then
         call MasslessVectorPropagator_real(this%current_list(ic)%val_r, &
              pp_loc)
      else
         call MasslessVectorPropagator(this%current_list(ic)%val_c, &
              pp_loc)
      endif
    end subroutine include_massless_vector_propagator

    subroutine include_massive_vector_propagator()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call MassiveVectorPropagator(this%current_list(ic)%val_c, &
           pp_loc,&
           this%current_list(ic)%mass,this%current_list(ic)%width)
    end subroutine include_massive_vector_propagator

    subroutine include_fermion_propagator()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call FermionPropagator(this%current_list(ic)%val_c, &
           pp_loc, &
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_fermion_propagator

    subroutine include_fermion_propagator_weyl()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call FermionPropagator_weyl(this%current_list(ic)%val_c, &
           pp_loc, &
           this%current_list(ic)%chirality)
    end subroutine include_fermion_propagator_weyl

    subroutine include_antifermion_propagator()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call AntifermionPropagator(this%current_list(ic)%val_c, &
           pp_loc,&
           this%current_list(ic)%mass,&
           this%current_list(ic)%width)
    end subroutine include_antifermion_propagator
    subroutine include_antifermion_propagator_weyl()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call AntifermionPropagator_weyl(this%current_list(ic)%val_c, &
           pp_loc,&
           this%current_list(ic)%chirality)
    end subroutine include_antifermion_propagator_weyl
    subroutine include_scalar_propagator()
      implicit none
      real(kind=8),dimension(0:3) :: pp_loc
      pp_loc=this%pp(0:3,this%pp_bin_to_i(this%current_list(ic)%bin))
      call ScalarPropagator(this%current_list(ic)%val_c, &
           pp_loc, &
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
         ,nOrd,allocation_status,max_nvals,max_row_classes,row_class,candidate_row
    integer(kind=8) :: workspace_bytes
    integer,dimension(n) :: iper,jper,part
    integer,dimension(:),allocatable :: n_vals
    integer,dimension(:,:),allocatable :: perm_local
    real(kind=8),dimension(1:3) :: col_fac
    real(kind=8),dimension(:,:),allocatable :: diff_vals
    real(kind=8),dimension(:,:,:),allocatable :: col_vals
    integer,dimension(:,:),allocatable :: ic,ir,n_colour_elements,unique_rows
    integer(kind=8),dimension(:),allocatable :: perm_dict
    integer,dimension(2*max_amplitude_external_particles) :: representative_rows,&
         representative_gi,representative_ui,row_class_to_unique,class_representative
    character(len=256) :: allocation_message

    write (99,*) 'Initialising colour matrix ...'
    if (col_acc.lt.0 .or. col_acc.gt.max_amplitude_external_particles) then
       write (*,*) 'Invalid colour-accuracy depth:',col_acc,&
            max_amplitude_external_particles
       stop 1
    endif
    if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
         this%nColOrd.lt.1 .or. this%n_amps.lt.1 .or. this%nprocs.lt.1) then
       write (*,*) 'Invalid amplitude dimensions while initialising colour matrix:',&
            n,this%nColOrd,this%n_amps,this%nprocs
       stop 1
    endif
    if (.not.allocated(this%processes) .or. .not.allocated(this%n_sing) .or. &
         .not.allocated(this%n_qqbar) .or. .not.allocated(this%perm) .or. &
         .not.allocated(this%iproc_start)) then
       write (*,*) 'Amplitude metadata is incomplete while initialising colour matrix'
       stop 1
    endif
    if (size(this%processes,1).ne.n .or. size(this%processes,2).ne.this%nprocs .or. &
         size(this%n_sing).ne.this%nprocs .or. size(this%n_qqbar).ne.this%nprocs .or. &
         size(this%iproc_start).lt.this%nprocs+1) then
       write (*,*) 'Amplitude metadata has incompatible colour-matrix dimensions'
       stop 1
    endif
    if (int(this%nColOrd,kind=8)*int(this%nColOrd,kind=8).gt.&
         max_colour_matrix_elements) then
       write (*,*) 'Colour matrix exceeds the supported comparison count:',&
            this%nColOrd,max_colour_matrix_elements
       stop 1
    endif
    if (allocated(this%i_col_i)) deallocate(this%i_col_i)
    if (allocated(this%col_index)) deallocate(this%col_index)
    if (allocated(this%row_index)) deallocate(this%row_index)
    if (allocated(this%n_col_vals)) deallocate(this%n_col_vals)
    if (allocated(this%diff_col_vals)) deallocate(this%diff_col_vals)
    if (this%nprocs.eq.1) then
       iproc=1
    elseif(this%nprocs.eq.3) then
       iproc=3
    else
       write (*,*) 'computation of color factor only for a single process at the time',this%nprocs
       stop 1
    endif
    if (this%n_sing(iproc).lt.0 .or. this%n_sing(iproc).gt.n .or. &
         this%n_qqbar(iproc).lt.0 .or. this%n_qqbar(iproc).gt.3) then
       write (*,*) 'Invalid subprocess classification while initialising colour matrix:',&
            iproc,this%n_sing(iproc),this%n_qqbar(iproc)
       stop 1
    endif
    part(1:n)=this%processes(1:n,iproc)
    ioff=this%iproc_start(iproc)-1
    nOrd=n-this%n_sing(iproc)
    if (nOrd.eq.0) then
       if (this%nColOrd.ne.1 .or. this%n_qqbar(iproc).ne.0 .or. &
            size(this%perm,1).ne.0 .or. ioff.lt.0 .or. &
            ioff+1.gt.size(this%perm,2) .or. ioff+1.gt.this%n_amps) then
          write (*,*) 'Invalid empty colour order while initialising colour matrix:',&
               this%nColOrd,this%n_qqbar(iproc),shape(this%perm),ioff,this%n_amps
          stop 1
       endif
       workspace_bytes=3_8*8_8+3_8*4_8+6_8*4_8+4_8*4_8
       allocate(this%n_col_vals(3),this%diff_col_vals(1,3),&
            this%i_col_i(1,3),this%row_index(0:1,1,3),&
            this%col_index(4),stat=allocation_status,errmsg=allocation_message)
       call require_colour_allocation('empty colour matrix')
       this%n_col_vals=1
       this%diff_col_vals=1d0
       this%i_col_i(1,:)=[1,2,3]
       this%row_index=0
       this%row_index(1,1,:)=1
       this%col_index=0
       this%col_index(2:4)=1
       write (99,*) '... empty colour matrix initialised'
       return
    endif
    if (nOrd.lt.1 .or. nOrd.gt.size(this%perm,1) .or. ioff.lt.0 .or. &
         ioff+this%nColOrd.gt.size(this%perm,2) .or. &
         ioff+this%nColOrd.gt.this%n_amps) then
       write (*,*) 'Invalid colour-order range while initialising colour matrix:',&
            nOrd,ioff,this%nColOrd,size(this%perm,1),size(this%perm,2),this%n_amps
       stop 1
    endif
    workspace_bytes=int(nOrd,kind=8)*int(this%nColOrd,kind=8)*4_8
    if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Colour-order workspace exceeds supported size:',&
            workspace_bytes,max_amplitude_workspace_bytes
       stop 1
    endif
    allocate(perm_local(nOrd,this%nColOrd),stat=allocation_status,&
         errmsg=allocation_message)
    call require_colour_allocation('colour permutations')
    do iperm=1,this%nColOrd
       do jperm=1,nOrd
          perm_local(jperm,iperm)=this%perm(jperm,ioff+iperm)
       enddo
       if (.not.valid_external_colour_order(perm_local(:,iperm),part)) then
          write (*,*) 'Invalid external colour order while initialising colour matrix:',&
               iperm,perm_local(:,iperm),part
          stop 1
       endif
    enddo

    if (this%n_qqbar(iproc).eq.3) then
       if (col_acc.eq.0) then
          call init_three_quark_line_lc_diagonal()
       else
          call init_three_quark_line_colour_matrix()
       endif
       write (99,*) '... colour matrix initialised'
       return
    endif

    allocate(n_vals(1:3),diff_vals(max_vals,1:3),&
         this%i_col_i(max_vals,1:3),n_colour_elements(max_vals,1:3),&
         stat=allocation_status,errmsg=allocation_message)
    call require_colour_allocation('colour-factor classification')
    n_vals=0
    diff_vals=0d0
    this%i_col_i=0
    n_colour_elements=0
    
    representative_rows=0
    representative_gi=0
    representative_ui=0
    row_class_to_unique=0
    class_representative=0
    if (this%n_qqbar(iproc).eq.0 .or. this%n_qqbar(iproc).eq.1) then
       n_unique_rows=1 ! all rows are similar
    elseif (this%n_qqbar(iproc).eq.2) then
       max_row_classes=2*(nOrd-3)
       if (max_row_classes.lt.2 .or. max_row_classes.gt.size(class_representative)) then
          write (*,*) 'Invalid two-quark-line colour-row class count:',&
               nOrd,max_row_classes
          stop 1
       endif
       do candidate_row=1,this%nColOrd
          call determine_gi_ui(perm_local(1:nOrd,candidate_row),gi,ui)
          if (gi.lt.0 .or. gi.gt.nOrd-4 .or. ui.lt.1 .or. ui.gt.2) then
             write (*,*) 'Invalid two-quark-line colour-row class:',&
                  candidate_row,gi,ui,nOrd
             stop 1
          endif
          row_class=(gi+1)+(ui-1)*(nOrd-3)
          if (class_representative(row_class).eq.0) &
               class_representative(row_class)=candidate_row
       enddo
       n_unique_rows=0
       do row_class=1,max_row_classes
          if (class_representative(row_class).eq.0) cycle
          n_unique_rows=n_unique_rows+1
          representative_rows(n_unique_rows)=class_representative(row_class)
          representative_gi(n_unique_rows)=mod(row_class-1,nOrd-3)
          representative_ui(n_unique_rows)=(row_class-1)/(nOrd-3)+1
          row_class_to_unique(row_class)=n_unique_rows
       enddo
       if (n_unique_rows.lt.1) then
          write (*,*) 'Two-quark-line colour matrix contains no row classes'
          stop 1
       endif
    else
       write (*,*) 'Inconsistent number of quark pairs',this%n_qqbar(iproc)
       stop 1
    endif
    if (use_cm_dict) then
       call create_perm_dict()
       workspace_bytes=int(max_keys,kind=8)*int(n_unique_rows,kind=8)*24_8+&
            int(nOrd,kind=8)*int(n_unique_rows,kind=8)*4_8
       if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
          write (*,*) 'Colour-factor lookup exceeds supported workspace:',&
               workspace_bytes,max_amplitude_workspace_bytes
          stop 1
       endif
       allocate(col_vals(1:3,max_keys,n_unique_rows),&
            unique_rows(1:nOrd,n_unique_rows),stat=allocation_status,&
            errmsg=allocation_message)
       call require_colour_allocation('colour-factor lookup')
    endif

! first check the unique rows in the colour matrix to determine how many
! different colour factors there are. This also sets up the library for the
! colour factors if use_cm_dict=.true.
    if (use_cm_dict) col_vals(1:3,1:max_keys,:)=0d0
    do iunique=1,n_unique_rows
       irow=0
       gi=0
       ui=0
       select case(this%n_qqbar(iproc))
       case(0,1)
          ! just take the first row: all rows are similar
          irow=1
       case(2)
          ! determine how the quarks are connected and how many gluons are on
          ! each quark line and make sure that it is not similar to one
          ! already considered
          call get_unique_row(iunique,irow,gi,ui)
       case default
          write (*,*) 'Invalid quark-line count while selecting a colour row:',&
               this%n_qqbar(iproc)
          stop 1
       end select
       if (irow.lt.1 .or. irow.gt.this%nColOrd) then
          write (*,*) 'Invalid representative colour row:',iunique,irow,this%nColOrd
          stop 1
       endif
       iper(1:nOrd)=perm_local(1:nOrd,irow)
       if (use_cm_dict) unique_rows(1:nOrd,iunique) = iper(1:nOrd)
       ! loop over the columns
       do jperm=1,this%nColOrd
          jper(1:nOrd)=perm_local(1:nOrd,jperm)
          gj=0
          uj=0
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
          isum=checked_add(isum,checked_multiply(n_colour_elements(ival,iacc),&
               this%nColOrd,'colour-index row count'),'colour-index size')
       enddo
    enddo

 ! Allocate the arrays now that we know their sizes
    max_nvals=maxval(n_vals(1:3))
    if (max_nvals.lt.1) then
       write (*,*) 'Colour matrix contains no non-zero colour factors'
       stop 1
    endif
    workspace_bytes=int(max_nvals,kind=8)*3_8*8_8+int(isum,kind=8)*4_8+&
         int(this%nColOrd+1,kind=8)*int(max_nvals,kind=8)*3_8*4_8+&
         int(max_nvals,kind=8)*3_8*8_8
    if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Compressed colour matrix exceeds supported workspace:',&
            workspace_bytes,max_amplitude_workspace_bytes
       stop 1
    endif
    allocate(ic(1:max_nvals,1:3),ir(1:max_nvals,1:3),&
         this%col_index(1:isum),&
         this%row_index(0:this%nColOrd,1:max_nvals,1:3),&
         this%n_col_vals(1:3),this%diff_col_vals(1:max_nvals,1:3),&
         stat=allocation_status,errmsg=allocation_message)
    call require_colour_allocation('compressed colour matrix')
    if (.not.allocated(ir)) then
       write (*,*) 'Compressed colour-row workspace was not allocated'
       stop 1
    endif
    this%row_index=0
    this%col_index(1)=0
    this%n_col_vals(1:3)=n_vals(1:3)
    this%diff_col_vals=0d0
    do iacc=1,3
       this%diff_col_vals(1:n_vals(iacc),iacc)=diff_vals(1:n_vals(iacc),iacc)
    enddo

! Compute all the colour factors and fill the col_index and row_index arrays
    ic=0
    ir=0
    do iperm=1,this%nColOrd
       iper(1:nOrd)=perm_local(1:nOrd,iperm)
       if (this%n_qqbar(iproc).eq.2)  call determine_gi_ui(iper,gi,ui)
       if (use_symm_cm) then
          jperm_lower=iperm
       else
          jperm_lower=1
       endif
       do jperm=jperm_lower,this%nColOrd
          jper(1:nOrd)=perm_local(1:nOrd,jperm)
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
    subroutine require_colour_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate ',trim(label),': ',&
              trim(allocation_message)
         stop 1
      endif
    end subroutine require_colour_allocation

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
      workspace_bytes=int(nrows,kind=8)*(28_8+int(nOrd,kind=8)*12_8)
      if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Three-line colour metadata exceeds supported workspace:',&
              workspace_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(endpoints(3,nrows),ngluons(3,nrows),&
           gluons(nOrd,3,nrows),flow_sign(nrows),stat=allocation_status,&
           errmsg=allocation_message)
      call require_colour_allocation('three-line colour metadata')
      call build_three_line_flow_metadata(endpoints,ngluons,gluons,flow_sign)
      call compute_three_line_color_factor(1,1,endpoints,ngluons,gluons,&
           flow_sign,three_line_col_fac)
      if (three_line_col_fac(1).eq.0d0) then
         write (*,*) 'Three-line leading-colour diagonal is zero'
         stop 1
      endif

      allocate(this%n_col_vals(3),this%diff_col_vals(1,3),&
           this%i_col_i(1,3),this%row_index(0:nrows,1,3),&
           this%col_index(nrows+1),stat=allocation_status,&
           errmsg=allocation_message)
      call require_colour_allocation('three-line leading-colour matrix')
      this%n_col_vals=(/1,0,0/)
      this%diff_col_vals=0d0
      this%diff_col_vals(1,1)=three_line_col_fac(1)
      this%i_col_i=0
      this%i_col_i(1,1)=1
      this%row_index=0
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

      if (use_symm_cm) then
         max_pairs=nrows*(nrows+1)/2
      else
         max_pairs=nrows*nrows
      endif
      workspace_bytes=int(nrows,kind=8)*(28_8+int(nOrd,kind=8)*12_8)+&
           int(nrows,kind=8)*int(nrows,kind=8)*24_8+&
           int(max_pairs,kind=8)*36_8
      if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Three-line colour matrix workspace exceeds supported size:',&
              workspace_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(endpoints(3,nrows),ngluons(3,nrows),&
           gluons(nOrd,3,nrows),flow_sign(nrows),&
           factors(nrows,nrows,3),nvalues(3),values(max_pairs,3),&
           counts(max_pairs,3),stat=allocation_status,&
           errmsg=allocation_message)
      call require_colour_allocation('three-line full-colour workspace')
      if (.not.allocated(factors)) then
         write (*,*) 'Three-line colour-factor workspace was not allocated'
         stop 1
      endif
      if (.not.allocated(values)) then
         write (*,*) 'Three-line colour-value workspace was not allocated'
         stop 1
      endif
      if (.not.allocated(counts)) then
         write (*,*) 'Three-line colour-count workspace was not allocated'
         stop 1
      endif
      call build_three_line_flow_metadata(endpoints,ngluons,gluons,flow_sign)

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
                  if (ival_local.gt.max_pairs) then
                     write (*,*) 'Too many distinct three-line colour factors:',&
                          iacc_local,ival_local,max_pairs
                     stop 1
                  endif
                  nvalues(iacc_local)=ival_local
                  values(ival_local,iacc_local)=factors(row,col,iacc_local)
               endif
               counts(ival_local,iacc_local)=counts(ival_local,iacc_local)+1
            enddo
            col=col+1
         enddo
      enddo

      max_nvals=maxval(nvalues)
      if (max_nvals.lt.1) then
         write (*,*) 'Three-line colour matrix contains no non-zero factors'
         stop 1
      endif
      total_entries=1+sum(counts)
      workspace_bytes=workspace_bytes+int(max_nvals,kind=8)*36_8+&
           int(nrows+1,kind=8)*int(max_nvals,kind=8)*12_8+&
           int(total_entries,kind=8)*4_8
      if (workspace_bytes.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Compressed three-line colour matrix exceeds supported workspace:',&
              workspace_bytes,max_amplitude_workspace_bytes
         stop 1
      endif
      allocate(this%n_col_vals(3),this%diff_col_vals(max_nvals,3),&
           this%i_col_i(max_nvals,3),&
           this%row_index(0:nrows,max_nvals,3),&
           this%col_index(total_entries),cursor(max_nvals,3),&
           stat=allocation_status,errmsg=allocation_message)
      call require_colour_allocation('compressed three-line colour matrix')
      if (.not.allocated(cursor)) then
         write (*,*) 'Compressed three-line colour cursor was not allocated'
         stop 1
      endif
      this%n_col_vals=nvalues
      this%diff_col_vals=0d0
      this%i_col_i=0
      this%row_index=0
      do iacc_local=1,3
         this%diff_col_vals(1:nvalues(iacc_local),iacc_local)=&
              values(1:nvalues(iacc_local),iacc_local)
      enddo

      this%col_index=0
      offset=1
      do iacc_local=1,3
         do ival_local=1,nvalues(iacc_local)
            this%i_col_i(ival_local,iacc_local)=offset
            offset=offset+counts(ival_local,iacc_local)
         enddo
      enddo

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
               if (ival_local.gt.nvalues(iacc_local)) then
                  write (*,*) 'Three-line colour factor changed during compression:',&
                       row,col,iacc_local,factors(row,col,iacc_local)
                  stop 1
               endif
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
      call parse_three_line_word(perm_local(1:nOrd,1),ref_q,ref_aq,&
           row_ngluons,row_gluons)
      do row=1,this%nColOrd
         call parse_three_line_word(perm_local(1:nOrd,row),row_q,row_aq,&
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
      label_is_quark=.false.
      if (label.lt.1 .or. label.gt.n) return
      if (label.le.2) then
         label_is_quark=part(label).le.-1 .and. part(label).ge.-6
      else
         label_is_quark=part(label).ge.1 .and. part(label).le.6
      endif
    end function label_is_quark

    logical function label_is_antiquark(label)
      implicit none
      integer,intent(in) :: label
      label_is_antiquark=.false.
      if (label.lt.1 .or. label.gt.n) return
      if (label.le.2) then
         label_is_antiquark=part(label).ge.1 .and. part(label).le.6
      else
         label_is_antiquark=part(label).le.-1 .and. part(label).ge.-6
      endif
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
      coef=0.0_16
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
      coef(1)=real(sgn,kind=16)
      coef_Nc(0,1)=int(sgn,kind=8)
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
      if (iunique.lt.1 .or. iunique.gt.n_unique_rows) then
         write (*,*) 'Invalid representative colour-row index:',&
              iunique,n_unique_rows
         stop 1
      endif
      irow=representative_rows(iunique)
      gi=representative_gi(iunique)
      ui=representative_ui(iunique)
      if (irow.lt.1 .or. irow.gt.this%nColOrd) then
         write (*,*) 'Missing representative colour row:',&
              iunique,irow,this%nColOrd,gi,ui
         stop 1
      endif
    end subroutine get_unique_row

    subroutine determine_gi(iper,gi)
      ! determine how many gluons are on the first quark line
      implicit none
      integer,intent(out) :: gi
      integer :: i
      integer,dimension(:),intent(in) :: iper
      gi=-1
      do i=2,nOrd-2
         if ((abs(part(iper(i))).ge.1.and.abs(part(iper(i))).le.6)) then
            gi=i-2
            return
         endif
      enddo
      write (*,*) 'Could not locate the second quark string in colour order:',iper
      stop 1
    end subroutine determine_gi

    subroutine determine_ui(iper,ui)
      implicit none
      integer,intent(out) :: ui
      integer,dimension(:),intent(in) :: iper
      ! check if quarks are connected in the same (or opposite way) as
      ! compared to the very first permutation
      if ((iper(1).eq.perm_local(1,1) .and. iper(nOrd).eq.perm_local(nOrd,1)) .or. &
           (iper(1).ne.perm_local(1,1) .and. iper(nOrd).ne.perm_local(nOrd,1))) then
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
      integer,intent(out) :: ui,gi
      integer,dimension(:),intent(in) :: iper
      if (size(iper).lt.nOrd) then
         write (*,*) 'Colour order is shorter than the classified coloured multiplicity:',&
              size(iper),nOrd
         stop 1
      endif
      call determine_gi(iper,gi)
      call determine_ui(iper,ui)
    end subroutine determine_gi_ui
    
   subroutine get_col_fac(iper,jper,ui,gi,col_fac)
     implicit none
     integer,intent(in) :: gi,ui
     integer,dimension(n),intent(in) :: iper,jper
     integer,dimension(n) :: col_new,row_first,row_per,col_per
     integer :: i,j,key,iunique,row_class
     real(kind=8),dimension(1:3),intent(out) :: col_fac
     if (this%n_qqbar(iproc).eq.2) then
        if (gi.lt.0 .or. gi.gt.nOrd-4 .or. ui.lt.1 .or. ui.gt.2) then
           write (*,*) 'Invalid two-quark-line colour lookup class:',gi,ui,nOrd
           stop 1
        endif
        row_class=(gi+1)+(ui-1)*(nOrd-3)
        iunique=row_class_to_unique(row_class)
        if (iunique.lt.1 .or. iunique.gt.n_unique_rows) then
           write (*,*) 'Colour row class is absent from the compressed lookup:',&
                gi,ui,row_class
           stop 1
        endif
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
        if (j.gt.nOrd) then
          write (*,*) 'Inconsistent colour-order label sets:',row_per(1:nOrd),&
               col_per(1:nOrd)
          stop 1
        endif
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
     solve_dict=0
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
     write(*,*) 'ERROR: colour permutation is missing from dictionary:',val
     stop 1
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
      max_keys=permutation_count(n,nOrd)
      if (max_keys.lt.1 .or. max_keys.gt.max_current_dictionary_entries .or. &
           int(max_keys,kind=8)*8_8.gt.max_amplitude_workspace_bytes) then
         write (*,*) 'Colour-permutation dictionary exceeds supported size:',&
              max_keys,max_current_dictionary_entries
         stop 1
      endif
      allocate(iper(1:nOrd),iper_in(1:nOrd),perm_dict(1:max_keys),&
           stat=allocation_status,errmsg=allocation_message)
      call require_colour_allocation('colour-permutation dictionary')
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
         if (iperm.lt.max_keys) then
            iper_in=iper
            call get_next_iperm(nOrd,iper_in,iper,n)
         endif
      enddo
      deallocate(iper)
    end subroutine create_perm_dict

    integer(kind=8) function get_value(nOrd,iper)
      ! Give a unique identifier based on the colour order. Simply convert the
      ! list to an integer with base equal to the number of elements in the
      ! order.
      implicit none
      integer :: j,nOrd
      integer(kind=8) :: base,digit
      integer,dimension(1:nOrd) :: iper
      if (n.lt.1 .or. nOrd.lt.0 .or. nOrd.gt.n) then
         write (*,*) 'ERROR: invalid colour-permutation key dimensions:',n,nOrd
         stop 1
      endif
      if (any(iper.lt.1) .or. any(iper.gt.n)) then
         write (*,*) 'ERROR: invalid entry in colour permutation:',iper
         stop 1
      endif
      base=int(n,kind=8)+1_8
      get_value=0_8
      do j=1,nOrd
         digit=int(iper(j),kind=8)
         if (get_value.gt.(huge(get_value)-digit)/base) then
            write (*,*) 'ERROR: colour-permutation dictionary key overflows 64-bit integer:',n,nOrd
            stop 1
         endif
         get_value=get_value*base+digit
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
               col_fac(1)=3d0**n
            endif
         elseif (this%n_qqbar(iproc).eq.1) then
            if (all(iper.eq.jper)) then
               col_fac(1)=3d0**(n-1)
            endif
         elseif (this%n_qqbar(iproc).eq.2) then
            if (all(iper.eq.jper)) then
               if (ui.eq.1 .and. uj.eq.1) then
                  col_fac(1)=3d0**(n-2)
               elseif (ui.eq.2 .and. uj.eq.2 .and. .not.this%same_flav(iproc)) then
                  col_fac(1)=3d0**(n-4) * 9d0 ! compensate for the already included factor 1/3 in qqbar->g Feynman rule
               elseif (ui.eq.2 .and. uj.eq.2 .and. this%same_flav(iproc)) then
                  col_fac(1)=3d0**(n-2)
               endif
            endif
         endif
      endif
      if (col_acc.ge.1) then ! NLC
         if (this%n_qqbar(iproc).eq.0) then
            if (all(iper.eq.jper)) then
               col_fac(2) = 3d0**n-dble(n)*3d0**(n-2)
            else
               call check_NLC(n,jper,iper,acc)
               col_fac(2)=dble(acc)*3d0**(n-2)
            endif
         elseif (this%n_qqbar(iproc).eq.1) then
            if (all(iper.eq.jper)) then
               col_fac(2) = 3d0**(n-1)-dble(n-2)*3d0**(n-3)
               ! include the full expansion
               call Tr_allocate(n)
               Tr(0,0,0)=1 ! one term
               Tr(0,0,1)=1 ! that term is single string of matrices
               Tr(0,1,1)=2*(n-2)
               Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
               Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
               coef(1)=1.0_16
               coef_Nc(:,:)=0
               coef_Nc(0,1)=1
               call Tr_full_simplify(col_factor) ! compute the colour factor by simplifying the product of traces
               col_fac(2)=dble(col_factor)
               call Tr_deallocate
            else
               call check_NLC_1qqbar(n,jper(2:n-1),iper(2:n-1),acc)
               col_fac(2)=dble(acc)*3d0**(n-3)
               ! include the full expansion
               if (acc.ne.0) then
                  call Tr_allocate(n)
                  Tr(0,0,0)=1 ! one term
                  Tr(0,0,1)=1 ! that term is single string of matrices
                  Tr(0,1,1)=2*(n-2)
                  Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
                  Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
                  coef(1)=1.0_16
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
               if (acc.eq.99) col_fac(2)=3d0**(n-2)-dble(n-4)*3d0**(n-4) ! LC interfence
               if (acc.le.1) col_fac(2)=dble(acc)*3d0**(n-3) ! NLC parts
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
                     coef(1)=1.0_16/9.0_16*9.0_16 ! compensate for the included 1/3 in qqbar->g
                  else
                     coef(1)=1.0_16
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
                     coef(1)=-1.0_16/3.0_16*3.0_16 ! compensate for the included 1/3 in qqbar->g
                  else
                     coef(1)=-1.0_16
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
            coef(1)=1.0_16
            coef_Nc(:,:)=0
            coef_Nc(0,1)=1
            ! compute the colour factor by simplifying the colour string
            call Tr_full_simplify(col_factor) 
            col_fac(3)=0d0
            do i=n,max(n-2*col_acc,0),-1  ! do not include any Nc
                                          ! contributions with negative
                                          ! powers, since they must cancel.
               if (i.ge.0) then
                  col_fac(3)=col_fac(3)+dble(coef_nc(i,0))*3d0**i
               else
                  col_fac(3)=col_fac(3)+dble(coef_nc(i,0))*3d0**i
               endif
            enddo
         elseif (this%n_qqbar(iproc).eq.1) then
            Tr(0,0,0)=1 ! one term
            Tr(0,0,1)=1 ! that term is single string of matrices
            Tr(0,1,1)=2*(n-2)
            Tr(1:n-2,1,1)=iper(2:n-1) ! the order of the matrices in each term
            Tr(n-1:2*(n-2),1,1)=jper(n-1:2:-1)
            coef(1)=1.0_16
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
                  coef(1)=1.0_16/9.0_16*9.0_16 ! compensate for the included 1/3 in qqbar->g
               else
                  coef(1)=1.0_16
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
                  coef(1)=-1.0_16/3.0_16*3.0_16 ! compensate for the included 1/3 in qqbar->g
               else
                  coef(1)=-1.0_16
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
      integer,intent(in) :: next
      integer,dimension(next-4),intent(in) :: iper,jper
      integer,dimension(next-4),intent(out) :: iper_new,jper_new
      integer :: i,j
      integer,dimension(next-4) :: imax,jmax
      if (next.lt.4) then
         write (*,*) 'Invalid gluon-string conversion size:',next
         stop 1
      endif
      iper_new=0
      jper_new=0
      if (next.eq.4) return
      do i=1,next-4
         if (count(iper.eq.iper(i)).ne.1 .or. &
              count(jper.eq.jper(i)).ne.1 .or. count(jper.eq.iper(i)).ne.1) then
            write (*,*) 'Invalid gluon-string permutations:',iper,jper
            stop 1
         endif
      enddo
      imax(1:next-4)=-1
      jmax(1:next-4)=-1
      do i=1,next-4
         do j=1,next-4
            if ((iper(j).gt.imax(i)).and..not.(any(iper(j).eq.imax(1:i-1)))) imax(i)=iper(j)
            if ((jper(j).gt.jmax(i)).and..not.(any(jper(j).eq.jmax(1:i-1)))) jmax(i)=jper(j)
         enddo
      enddo
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

  subroutine record_optimisation_sample(this,sample_index,nsamples)
    ! Retain several independent evaluations before identifying reusable
    ! currents.  A single phase-space point is not enough: distinct helicity
    ! currents can agree accidentally at special kinematics.
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    integer,intent(in) :: sample_index,nsamples
    integer :: ic,max_dim,value_dim,allocation_status
    integer(kind=8) :: sample_elements
    character(len=256) :: allocation_message

    if (sample_index.lt.1 .or. sample_index.gt.nsamples) then
       write (*,*) 'Invalid amplitude-optimisation sample index',sample_index,nsamples
       stop 1
    endif
    if (.not.this%evaluation_workspace_ready) then
       write (*,*) 'Cannot sample an amplitude before a successful evaluation'
       stop 1
    endif
    if (sample_index.gt.this%optimisation_sample_count+1) then
       write (*,*) 'Amplitude-optimisation samples must be recorded contiguously:',&
            sample_index,this%optimisation_sample_count
       stop 1
    endif
    if (.not.allocated(this%optimisation_current_samples)) then
       max_dim=1
       do ic=1,this%n_cur
          if (allocated(this%current_list(ic)%val_c)) &
               max_dim=max(max_dim,size(this%current_list(ic)%val_c))
          if (allocated(this%current_list(ic)%val_r)) &
               max_dim=max(max_dim,size(this%current_list(ic)%val_r))
       enddo
       if (max_dim.lt.1 .or. max_dim.gt.6 .or. this%n_cur.lt.1 .or. &
            this%n_cur.gt.max_amplitude_current_records .or. nsamples.lt.1) then
          write (*,*) 'Invalid amplitude-optimisation sample dimensions:',&
               max_dim,this%n_cur,nsamples
          stop 1
       endif
       sample_elements=int(max_dim,kind=8)*int(this%n_cur,kind=8)*&
            int(nsamples,kind=8)
       if (sample_elements.gt.max_amplitude_workspace_bytes/16_8) then
          write (*,*) 'Amplitude-optimisation samples exceed supported workspace:',&
               max_dim,this%n_cur,nsamples,max_amplitude_workspace_bytes
          stop 1
       endif
       allocate(this%optimisation_current_samples(max_dim,this%n_cur,nsamples),&
            stat=allocation_status,errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate amplitude-optimisation samples: ',&
               trim(allocation_message)
          stop 1
       endif
       this%optimisation_current_samples=(0d0,0d0)
    elseif (size(this%optimisation_current_samples,1).lt.1 .or. &
         size(this%optimisation_current_samples,1).gt.6 .or. &
         size(this%optimisation_current_samples,2).ne.this%n_cur .or. &
         size(this%optimisation_current_samples,3).ne.nsamples) then
       write (*,*) 'Amplitude-optimisation sample storage is inconsistent'
       stop 1
    endif

    this%optimisation_current_samples(:,:,sample_index)=(0d0,0d0)
    do ic=1,this%n_cur
       if (allocated(this%current_list(ic)%val_c)) then
          value_dim=size(this%current_list(ic)%val_c)
          this%optimisation_current_samples(1:value_dim,ic,sample_index)=&
               this%current_list(ic)%val_c
       elseif (allocated(this%current_list(ic)%val_r)) then
          value_dim=size(this%current_list(ic)%val_r)
          this%optimisation_current_samples(1:value_dim,ic,sample_index)=&
               cmplx(this%current_list(ic)%val_r,0d0,kind=8)
       else
          write (*,*) 'Cannot sample an unevaluated amplitude current',ic
          stop 1
       endif
    enddo
    if (.not.all(complex_value_is_finite(&
         this%optimisation_current_samples(:,:,sample_index)))) then
       write (*,*) 'Non-finite current encountered during amplitude optimisation',&
            sample_index
       stop 1
    endif
    this%optimisation_sample_count=max(this%optimisation_sample_count,sample_index)
  end subroutine record_optimisation_sample

  subroutine clear_optimisation_samples(this)
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    if (allocated(this%optimisation_current_samples)) &
         deallocate(this%optimisation_current_samples)
    this%optimisation_sample_count=0
  end subroutine clear_optimisation_samples

  subroutine optimise_evaluation(this,n)
    ! Re-use currents only when their metadata is compatible and their
    ! values agree at every recorded warm-up point.  Interactions are merged
    ! only when they are structurally identical after current remapping.
    implicit none
    class(amplitude_QCD),intent(inout) :: this
    integer,intent(in) :: n
    integer :: isize,ic1,ic2,iv1,iv2,i,ic,iv,allocation_status
    integer(kind=8) :: workspace_bytes,comparison_count
    integer,dimension(:,:),allocatable :: map_cur,map_vert
    logical,dimension(:),allocatable :: include_cur,include_vert
    character(len=256) :: allocation_message

    if (.not.allocated(this%optimisation_current_samples) .or. &
         this%optimisation_sample_count.lt.2) then
       write (*,*) 'Amplitude-current optimisation requires at least two warm-up samples'
       stop 1
    endif
    if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
         this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
         this%n_vert.lt.0 .or. this%n_vert.gt.max_amplitude_interaction_records .or. &
         size(this%optimisation_current_samples,1).lt.1 .or. &
         size(this%optimisation_current_samples,1).gt.6 .or. &
         size(this%optimisation_current_samples,2).ne.this%n_cur .or. &
         this%optimisation_sample_count.gt.&
         size(this%optimisation_current_samples,3)) then
       write (*,*) 'Invalid amplitude-optimisation state:',n,this%n_cur,&
            this%n_vert,this%optimisation_sample_count,&
            shape(this%optimisation_current_samples)
       stop 1
    endif
    if (.not.allocated(this%n_cur_start) .or. .not.allocated(this%n_cur_end) .or. &
         .not.allocated(this%n_vert_start) .or. .not.allocated(this%n_vert_end) .or. &
         .not.allocated(this%curr2amp)) then
       write (*,*) 'Amplitude-optimisation metadata is incomplete'
       stop 1
    endif
    if (size(this%n_cur_start).ne.n .or. size(this%n_cur_end).ne.n .or. &
         size(this%n_vert_start).ne.n-2 .or. size(this%n_vert_end).ne.n-2 .or. &
         size(this%curr2amp,1).ne.2 .or. size(this%curr2amp,2).lt.this%n_amps) then
       write (*,*) 'Amplitude-optimisation metadata has incompatible dimensions'
       stop 1
    endif
    workspace_bytes=16_8*size(this%optimisation_current_samples,kind=8)+&
         int(this%n_cur+this%n_vert,kind=8)*&
         max(1_8,int(storage_size(.true.),kind=8)/8_8)+&
         8_8*(int(this%n_cur,kind=8)+int(this%n_vert,kind=8)+2_8)
    if (workspace_bytes.lt.0_8 .or. &
         workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Amplitude optimisation exceeds the supported workspace:',&
            workspace_bytes,max_amplitude_workspace_bytes
       stop 1
    endif
    allocate(include_cur(1:this%n_cur),include_vert(1:this%n_vert),&
         map_cur(0:this%n_cur,1:2),map_vert(0:this%n_vert,1:2),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate amplitude-optimisation maps: ',&
            trim(allocation_message)
       stop 1
    endif
    include_cur=.true.
    include_vert=.true.
    map_cur=0
    map_vert=0
    map_cur(0,1)=0
    map_vert(0,1)=0
    comparison_count=0_8

    do isize=1,n-1
       do ic1=this%n_cur_start(isize),this%n_cur_end(isize)-1
          if (.not.include_cur(ic1)) cycle
          do ic2=ic1+1,this%n_cur_end(isize)
             if (.not.include_cur(ic2)) cycle
             call count_optimisation_comparison('current')
             if (.not.current_metadata_matches(ic1,ic2)) cycle
             if (.not.current_samples_match(ic1,ic2)) cycle
             map_cur(0,1)=map_cur(0,1)+1
             map_cur(map_cur(0,1),1)=ic2
             map_cur(map_cur(0,1),2)=ic1
             include_cur(ic2)=.false.
          enddo
       enddo
    enddo
    do i=1,map_cur(0,1)
       do iv1=1,this%n_vert
          if (this%interaction_list(iv1)%currents(1).eq.map_cur(i,1)) &
               this%interaction_list(iv1)%currents(1)=map_cur(i,2)
          if (this%interaction_list(iv1)%currents(2).eq.map_cur(i,1)) &
               this%interaction_list(iv1)%currents(2)=map_cur(i,2)
       enddo
       where (this%curr2amp.eq.map_cur(i,1)) this%curr2amp=map_cur(i,2)
       if (allocated(this%three_line_partner_curr2amp)) then
          where (this%three_line_partner_curr2amp.eq.map_cur(i,1)) &
               this%three_line_partner_curr2amp=map_cur(i,2)
       endif
    enddo

    do isize=2,n-1
       do iv1=this%n_vert_start(isize),this%n_vert_end(isize)-1
          if (.not.include_vert(iv1)) cycle
          do iv2=iv1+1,this%n_vert_end(isize)
             if (.not.include_vert(iv2)) cycle
             call count_optimisation_comparison('interaction')
             if (.not.interactions_match(iv1,iv2)) cycle
             map_vert(0,1)=map_vert(0,1)+1
             map_vert(map_vert(0,1),1)=iv2
             map_vert(map_vert(0,1),2)=iv1
             include_vert(iv2)=.false.
          enddo
       enddo
    enddo
    do i=1,map_vert(0,1)
       do ic1=1,this%n_cur
          do iv1=1,this%current_list(ic1)%n_vert
             if (this%current_list(ic1)%vertices(iv1).eq.map_vert(i,1)) &
                  this%current_list(ic1)%vertices(iv1)=map_vert(i,2)
          enddo
       enddo
    enddo

    if (.not.allocated(this%include_amp)) then
       allocate(this%include_amp(1:this%n_amps),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate amplitude-inclusion map: ',&
               trim(allocation_message)
          stop 1
       endif
    elseif (size(this%include_amp).lt.this%n_amps) then
       write (*,*) 'Amplitude-inclusion map is too short during optimisation'
       stop 1
    endif
    this%include_amp(1:this%n_amps)=.true.
    ! filter_dead_trees conservatively keeps every terminal current unless
    ! told which ones still close an amplitude.  Supplying the remapped
    ! closing-current set is what makes terminal-current sharing effective.
    include_cur=.false.
    do i=1,this%n_amps
       do ic=1,2
          if (this%curr2amp(ic,i).ne.0) include_cur(this%curr2amp(ic,i))=.true.
          if (allocated(this%three_line_partner_curr2amp)) then
             if (i.le.size(this%three_line_partner_curr2amp,2)) then
                if (this%three_line_partner_curr2amp(ic,i).ne.0) &
                     include_cur(this%three_line_partner_curr2amp(ic,i))=.true.
             endif
          endif
       enddo
    enddo
    deallocate(map_cur,map_vert,include_vert)
    call this%filter_dead_trees(n,include_cur)
    do ic=1,this%n_cur
       if (allocated(this%current_list(ic)%val_c)) deallocate(this%current_list(ic)%val_c)
       if (allocated(this%current_list(ic)%val_r)) deallocate(this%current_list(ic)%val_r)
    enddo
    do iv=1,this%n_vert
       if (allocated(this%interaction_list(iv)%val_c)) deallocate(this%interaction_list(iv)%val_c)
       if (allocated(this%interaction_list(iv)%val_r)) deallocate(this%interaction_list(iv)%val_r)
    enddo
    this%evaluation_workspace_ready=.false.
    call this%clear_optimisation_samples()
    write (99,*) 'Total number of currents, vertices and amplitudes after optimisation',&
         this%n_cur,this%n_vert,this%n_amps

  contains

    subroutine count_optimisation_comparison(label)
      character(len=*),intent(in) :: label
      if (comparison_count.ge.max_amplitude_optimisation_comparisons) then
         write (*,*) 'Amplitude ',trim(label),&
              ' optimisation exceeds the supported comparison budget:',&
              max_amplitude_optimisation_comparisons
         stop 1
      endif
      comparison_count=comparison_count+1_8
    end subroutine count_optimisation_comparison

    logical function current_metadata_matches(first,second)
      integer,intent(in) :: first,second
      current_metadata_matches=.false.
      if (this%current_list(first)%type.ne.this%current_list(second)%type) return
      if (this%current_list(first)%bin.ne.this%current_list(second)%bin) return
      if (this%current_list(first)%chirality.ne.this%current_list(second)%chirality) return
      if (this%current_list(first)%mass.ne.this%current_list(second)%mass) return
      if (this%current_list(first)%width.ne.this%current_list(second)%width) return
      if (allocated(this%current_list(first)%val_c).neqv.&
           allocated(this%current_list(second)%val_c)) return
      if (allocated(this%current_list(first)%val_r).neqv.&
           allocated(this%current_list(second)%val_r)) return
      if (allocated(this%current_list(first)%val_c)) then
         if (size(this%current_list(first)%val_c).ne.&
              size(this%current_list(second)%val_c)) return
      endif
      if (allocated(this%current_list(first)%val_r)) then
         if (size(this%current_list(first)%val_r).ne.&
              size(this%current_list(second)%val_r)) return
      endif
      current_metadata_matches=.true.
    end function current_metadata_matches

    logical function current_samples_match(first,second)
      integer,intent(in) :: first,second
      integer :: isample
      real(kind=8) :: difference,scale,normalized_scale
      real(kind=8),parameter :: equality_tolerance=1d-11
      logical :: compared
      current_samples_match=.false.
      compared=.false.
      do isample=1,this%optimisation_sample_count
         if (.not.all(complex_value_is_finite(&
              this%optimisation_current_samples(:,first,isample)))) return
         if (.not.all(complex_value_is_finite(&
              this%optimisation_current_samples(:,second,isample)))) return
         scale=max(&
              maxval(abs(real(this%optimisation_current_samples(:,first,isample),kind=8))),&
              maxval(abs(aimag(this%optimisation_current_samples(:,first,isample)))),&
              maxval(abs(real(this%optimisation_current_samples(:,second,isample),kind=8))),&
              maxval(abs(aimag(this%optimisation_current_samples(:,second,isample)))))
         if (scale.le.tiny(1d0)) cycle
         compared=.true.
         normalized_scale=&
              sum(abs(this%optimisation_current_samples(:,first,isample)/scale))+&
              sum(abs(this%optimisation_current_samples(:,second,isample)/scale))
         difference=sum(abs(&
              this%optimisation_current_samples(:,first,isample)/scale-&
              this%optimisation_current_samples(:,second,isample)/scale))
         if (difference.gt.equality_tolerance*normalized_scale) return
      enddo
      current_samples_match=compared
    end function current_samples_match

    logical function interactions_match(first,second)
      integer,intent(in) :: first,second
      interactions_match=.false.
      if (this%interaction_list(first)%type.ne.this%interaction_list(second)%type) return
      if (this%interaction_list(first)%chirality.ne.&
           this%interaction_list(second)%chirality) return
      if (any(this%interaction_list(first)%currents.ne.&
           this%interaction_list(second)%currents)) return
      if (any(this%interaction_list(first)%coupl.ne.&
           this%interaction_list(second)%coupl)) return
      if (allocated(this%interaction_list(first)%singlet_mv).neqv.&
           allocated(this%interaction_list(second)%singlet_mv)) return
      if (allocated(this%interaction_list(first)%singlet_mv)) then
         if (size(this%interaction_list(first)%singlet_mv).ne.&
              size(this%interaction_list(second)%singlet_mv)) return
         if (any(this%interaction_list(first)%singlet_mv.ne.&
              this%interaction_list(second)%singlet_mv)) return
      endif
      interactions_match=.true.
    end function interactions_match

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
    character(len=512) :: line
    character(len=64) :: tmp
    character(len=256) :: io_message,allocation_message
    integer :: ip,ibin,i,isize,ih_in,ifinal,ic,iv,iamp,iproc,itype,j,ii,jj,idau,vkey,idim,ios
    integer :: max_current_vertices,max_interaction_bucket,allocation_status
    integer(kind=8) :: generator_workspace_bytes
    integer,dimension(0:24,0:8) :: icount
    integer,dimension(:,:),allocatable :: icount_type
    integer,dimension(:,:),allocatable :: curs
    integer,dimension(:),allocatable :: pp
    real(kind=8),dimension(:),allocatable :: m,w
    integer,dimension(:,:,:),allocatable :: cur1,cur2,int1,pp1,pp2
    real(kind=8),dimension(:,:,:,:),allocatable :: coupl
    if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
         this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
         this%n_vert.lt.0 .or. this%n_vert.gt.max_amplitude_interaction_records .or. &
         this%n_amps.lt.1) then
       write (*,*) 'Invalid dimensions while creating amplitude library:',&
            n,this%n_cur,this%n_vert,this%n_amps
       stop 1
    endif
    if (.not.allocated(this%current_list) .or. &
         .not.allocated(this%interaction_list) .or. .not.allocated(this%amps) .or. &
         .not.allocated(this%n_cur_start) .or. .not.allocated(this%n_cur_end) .or. &
         .not.allocated(this%n_vert_start) .or. .not.allocated(this%n_vert_end) .or. &
         .not.allocated(this%pp_bin_to_i) .or. .not.allocated(this%pp_i_to_bin)) then
       write (*,*) 'Amplitude library requested before the amplitude was fully evaluated'
       stop 1
    endif
    if (size(this%current_list).lt.this%n_cur .or. &
         size(this%interaction_list).lt.this%n_vert .or. &
         size(this%amps).lt.this%n_amps) then
       write (*,*) 'Amplitude-library storage is inconsistent with its dimensions'
       stop 1
    endif
    if (any(this%current_list(1:this%n_cur)%n_vert.lt.0) .or. &
         any(this%current_list(1:this%n_cur)%n_vert.gt.this%n_vert)) then
       write (*,*) 'Invalid current interaction count while creating amplitude library'
       stop 1
    endif
    max_current_vertices=max(1,maxval(this%current_list(1:this%n_cur)%n_vert))
    generator_workspace_bytes=4_8*int(max_current_vertices,kind=8)*11_8
    if (generator_workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Amplitude-library current counts exceed supported workspace:',&
            max_current_vertices,max_amplitude_workspace_bytes
       stop 1
    endif
    allocate(icount_type(max_current_vertices,11),stat=allocation_status,&
         errmsg=allocation_message)
    call require_generator_allocation('amplitude-library current counts')
    write(tmp,*) igroup
    write(line,*) iint
    line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
    open(file=line,unit=iunit,form='unformatted',access='stream',status='replace',&
         action='write',iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not create amplitude-library reference data: ',trim(line),&
            trim(io_message)
       stop 1
    endif
    write(iunit,iostat=ios,iomsg=io_message) p
    call require_generator_io('amplitude-library reference momenta')
    write(iunit,iostat=ios,iomsg=io_message) this%amps
    call require_generator_io('amplitude-library reference amplitudes')
    close(iunit,iostat=ios,iomsg=io_message)
    call require_generator_io('amplitude-library reference data')
    
    write(tmp,*) igroup
    write(line,*) iint
    line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.f03'
    open(file=line,unit=iunit,status='replace',action='write',iostat=ios,&
         iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not create amplitude-library source: ',trim(line),&
            trim(io_message)
       stop 1
    endif
    write(line,*) iint
    write(iunit,'(a)') 'module amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib'
    write(iunit,'(2x,a)') 'use FeynmanRules'
    write(iunit,'(2x,a)') 'implicit none'
    write(iunit,'(2x,a)') 'private'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(2x,a)') 'public :: evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))
    write(iunit,'(2x,a)') 'contains'
    write(iunit,'(2x,a)') 'subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'(p,amps)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) n
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(tmp,*) this%n_amps
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//'),intent(out) :: amps'
    write(tmp,*) this%max_pp
    write(iunit,'(4x,a)') 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//') :: pp'
    write(tmp,*) this%n_cur
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//') :: val_c'
    write(tmp,*) this%n_vert
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//') :: int_c'
    write(iunit,'(4x,a)') 'amps=(0d0,0d0)'
    write(iunit,'(4x,a)') 'call reset_feynman_numerical_status()'
    write(iunit,'(4x,a)') 'call validate_feynman_momenta(p)'
    write(iunit,'(4x,a)') 'if (.not.feynman_numerical_status_ok()) return'
    write(iunit,'(4x,a)') 'pp=0d0'
    write(iunit,'(4x,a)') 'val_c=(0d0,0d0)'
    write(iunit,'(4x,a)') 'int_c=(0d0,0d0)'
    write(iunit,'(4x,a)') 'call fill_momentum_array(p,pp)'
    write(iunit,'(4x,a)') 'call validate_feynman_momenta(pp)'
    write(iunit,'(4x,a)') 'if (.not.feynman_numerical_status_ok()) return'
    write(iunit,'(4x,a)') 'call compute_external_currents(pp,val_c)'
    write(iunit,'(4x,a)') 'if (.not.feynman_numerical_status_ok()) return'
    do isize=2,n-1
       write(tmp,*) isize
       write(iunit,'(4x,a)') 'call compute_vertices'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
       write(iunit,'(4x,a)') 'if (.not.feynman_numerical_status_ok()) return'
       write(iunit,'(4x,a)') 'call compute_currents'//trim(adjustl(tmp))//'(pp,val_c,int_c)'
       write(iunit,'(4x,a)') 'if (.not.feynman_numerical_status_ok()) return'
    enddo
    write(iunit,'(4x,a)') 'call compute_amps(amps,val_c)'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(2x,a)') 'end subroutine evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))
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

       ! Count the interaction buckets first.  The old implementation used
       ! arrays dimensioned by n_vert for every one of the 225 possible
       ! (type,key) buckets, which could exhaust the stack and reserve orders
       ! of magnitude more memory than the generated source needs.
       icount(0:24,0:8)=0
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          call get_library_vertex_bucket(iv,itype,vkey)
          icount(itype,vkey)=icount(itype,vkey)+1
       enddo
       ! Keep even an empty size bucket allocated.  This gives all later
       ! bucket-emission paths a single, explicit allocation invariant while
       ! retaining zero counts for sizes that have no interactions.
       max_interaction_bucket=max(1,maxval(icount))
       generator_workspace_bytes=int(max_interaction_bucket,kind=8)*25_8*9_8*&
            (5_8*4_8+2_8*8_8)
       if (generator_workspace_bytes.gt.max_amplitude_workspace_bytes) then
          write (*,*) 'Amplitude-library interaction buckets exceed supported workspace:',&
               max_interaction_bucket,generator_workspace_bytes,&
               max_amplitude_workspace_bytes
          stop 1
       endif
       allocate(cur1(max_interaction_bucket,0:24,0:8),&
            cur2(max_interaction_bucket,0:24,0:8),&
            int1(max_interaction_bucket,0:24,0:8),&
            pp1(max_interaction_bucket,0:24,0:8),&
            pp2(max_interaction_bucket,0:24,0:8),&
            coupl(1:2,max_interaction_bucket,0:24,0:8),&
            stat=allocation_status,errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate amplitude-library interaction buckets: ',&
               trim(allocation_message)
          stop 1
       endif
       if (.not.allocated(coupl)) then
          write (*,*) 'Amplitude-library coupling workspace was not allocated'
          stop 1
       endif
       icount=0
       do iv=this%n_vert_start(isize),this%n_vert_end(isize)
          call get_library_vertex_bucket(iv,itype,vkey)
          icount(itype,vkey)=icount(itype,vkey)+1
          cur1(icount(itype,vkey),itype,vkey)=&
               this%interaction_list(iv)%currents(1)
          cur2(icount(itype,vkey),itype,vkey)=&
               this%interaction_list(iv)%currents(2)
          pp1(icount(itype,vkey),itype,vkey)=this%pp_bin_to_i(&
               this%current_list(this%interaction_list(iv)%currents(1))%bin)
          pp2(icount(itype,vkey),itype,vkey)=this%pp_bin_to_i(&
               this%current_list(this%interaction_list(iv)%currents(2))%bin)
          int1(icount(itype,vkey),itype,vkey)=iv
          coupl(1:2,icount(itype,vkey),itype,vkey)=&
               this%interaction_list(iv)%coupl(1:2)
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
                write(iunit,'(6x,a)') 'call LeptonAntileptonToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             else
                write(iunit,'(6x,a)') 'call LeptonAntileptonToVector_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(tmp,*) vkey/3-1
                line='[coupl(2*i-1),coupl(2*i)],'//trim(adjustl(tmp))//','
                write(tmp,*) mod(vkey,3)-1
                line=trim(adjustl(line))//trim(adjustl(tmp))//')'
                write(iunit,'(8x,a)') trim(adjustl(line))
             endif
             line=''
          elseif(itype.eq.22) then
             if (vkey.eq.4) then
                write(iunit,'(6x,a)') 'call AntileptonLeptonToVector(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
                write(iunit,'(8x,a)') '[coupl(2*i-1),coupl(2*i)])'
             else
                write(iunit,'(6x,a)') 'call AntileptonLeptonToVector_weyl(val_c(1,cur1(i)),val_c(1,cur2(i)),int_c(1,int1(i)), &'
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
          if (this%interaction_list(int1(1,itype,vkey))%chirality.ne.0) then
             idim=2
          else
             idim=pm%get_inter_dim(itype)
          endif
          if (idim.lt.1 .or. idim.gt.6) then
             write (*,*) 'Invalid generated interaction dimension:',itype,vkey,idim
             stop 1
          endif
          write(tmp,*) idim
          write(iunit,'(6x,a)') 'call validate_complex_wavefunction(int_c(1:'//&
               trim(adjustl(tmp))//',int1(i)))'
          write(iunit,'(6x,a)') 'if (.not.feynman_numerical_status_ok()) exit'
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
       if (allocated(cur1)) then
          deallocate(cur1,cur2,int1,pp1,pp2,coupl)
       endif

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
          if (this%current_list(ic)%n_vert.lt.1 .or. &
               this%current_list(ic)%n_vert.gt.max_current_vertices) then
             write (*,*) 'Invalid current interaction count while generating library:',&
                  ic,this%current_list(ic)%n_vert
             stop 1
          endif
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

             generator_workspace_bytes=(4_8*int(i,kind=8)+24_8)*&
                  int(icount_type(i,j),kind=8)
             if (generator_workspace_bytes.gt.&
                  max_amplitude_workspace_bytes) then
                write (*,*) 'Amplitude-library current bucket exceeds supported workspace:',&
                     i,j,icount_type(i,j),generator_workspace_bytes
                stop 1
             endif
             allocate(curs(0:i,icount_type(i,j)),pp(icount_type(i,j)),&
                  m(icount_type(i,j)),w(icount_type(i,j)),&
                  stat=allocation_status,errmsg=allocation_message)
             call require_generator_allocation('amplitude-library current bucket')
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
             if (ii.ne.icount_type(i,j)) then
                write (*,*) 'Amplitude-library current bucket count changed:',&
                     i,j,ii,icount_type(i,j)
                stop 1
             endif
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
             if (j.eq.6) then
                idim=6
             elseif (j.eq.5 .or. j.eq.7) then
                idim=1
             elseif (j.ge.8 .and. j.le.11) then
                idim=2
             else
                idim=4
             endif
             write(tmp,*) idim
             write(iunit,'(6x,a)') 'call validate_complex_wavefunction(val_c(1:'//&
                  trim(adjustl(tmp))//',int1(0,i)))'
             write(iunit,'(6x,a)') 'if (.not.feynman_numerical_status_ok()) exit'
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
             write(iunit,'(6x,a)') 'if (.not.feynman_numerical_status_ok()) exit'
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

    write(iunit,'(2x,a)') 'subroutine compute_amps(amps,val_c)'
    write(iunit,'(4x,a)') 'implicit none'
    write(tmp,*) this%n_amps
    write(iunit,'(4x,a)') 'complex(kind=8),dimension('//trim(adjustl(tmp))//'),intent(out) :: amps'
    write(tmp,*) this%n_cur
    write(iunit,'(4x,a)') 'complex(kind=8),dimension(1:6,'//trim(adjustl(tmp))//'),intent(in) :: val_c'

    ! first the 'non-same-flavour' ones
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          if (.not.this%same_flav(iproc)) then
             write(tmp,*) iamp
             line='amps('//trim(adjustl(tmp))//')=ContractFermionCurrents('
             write(tmp,*) this%curr2amp(1,iamp)
             line=trim(adjustl(line))//'val_c(1,'//trim(adjustl(tmp))//'),'
             write(tmp,*) this%current_list(this%curr2amp(1,iamp))%chirality
             line=trim(adjustl(line))//trim(adjustl(tmp))//',val_c(1,'
             write(tmp,*) this%curr2amp(2,iamp)
             line=trim(adjustl(line))//trim(adjustl(tmp))//'),'
             write(tmp,*) this%current_list(this%curr2amp(2,iamp))%chirality
             line=trim(adjustl(line))//trim(adjustl(tmp))//')'
             write(iunit,'(4x,a)') trim(adjustl(line))
          endif
       enddo
    enddo
    ! now the same-flavour ones. They are the "sum" of two non-same-flavour ones.
    do iproc=1,this%nprocs
       do iamp=this%iproc_start(iproc),this%iproc_start(iproc+1)-1
          if (this%same_flav(iproc)) then
             ! same-flavour amps are build from two different-flavour amps
             write(tmp,*) iamp
             line='amps('//trim(adjustl(tmp))//')='
             do idau=1,2
                if (this%same_flavour_sum(iamp,idau).gt.0) then
                   write(tmp,*) this%same_flavour_sum(iamp,idau)
                   if (this%same_flavour_sum_operation(iamp,idau) .eq. 0) then
                      line=trim(adjustl(line))//'+amps('//trim(adjustl(tmp))//')'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 1) then
                      line=trim(adjustl(line))//'-conjg(amps('//trim(adjustl(tmp))//'))'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 2) then
                      line=trim(adjustl(line))//'+conjg(amps('//trim(adjustl(tmp))//'))'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 3) then
                      line=trim(adjustl(line))//'-amps('//trim(adjustl(tmp))//')'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 4) then
                      line=trim(adjustl(line))//'+cmplx(aimag(amps('//trim(adjustl(tmp))// &
                           ')),dble(amps('//trim(adjustl(tmp))//')),kind=8)'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 5) then
                      line=trim(adjustl(line))//'+cmplx(-aimag(amps('//trim(adjustl(tmp))// &
                           ')),dble(amps('//trim(adjustl(tmp))//')),kind=8)'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 6) then
                      line=trim(adjustl(line))//'+cmplx(aimag(amps('//trim(adjustl(tmp))// &
                           ')),-dble(amps('//trim(adjustl(tmp))//')),kind=8)'
                   elseif (this%same_flavour_sum_operation(iamp,idau) .eq. 7) then
                      line=trim(adjustl(line))//'+cmplx(-aimag(amps('//trim(adjustl(tmp))// &
                           ')),-dble(amps('//trim(adjustl(tmp))//')),kind=8)'
                   else
                      write (*,*) 'ERROR: unknown operation in creating library', &
                           this%same_flavour_sum_operation(iamp,idau)
                      stop 1
                   endif
                endif
             enddo
             write(iunit,'(4x,a)') trim(adjustl(line))
          endif
       enddo
    enddo
    write(iunit,'(4x,a)') 'call validate_complex_wavefunction(amps)'
    write(iunit,'(2x,a)') 'end subroutine compute_amps'
    write(tmp,*) igroup
    write(line,*) iint
    write(iunit,'(a)') 'end module amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib'
    close(iunit,iostat=ios,iomsg=io_message)
    call require_generator_io('amplitude-library source')
    deallocate(icount_type)
  contains
    subroutine get_library_vertex_bucket(vertex_index,bucket_type,bucket_key)
      implicit none
      integer,intent(in) :: vertex_index
      integer,intent(out) :: bucket_type,bucket_key
      integer :: first_current,second_current,first_chirality,second_chirality,&
           result_chirality,first_bin,second_bin

      if (vertex_index.lt.1 .or. vertex_index.gt.this%n_vert) then
         write (*,*) 'Out-of-range interaction while creating amplitude library:',&
              vertex_index,this%n_vert
         stop 1
      endif
      bucket_type=this%interaction_list(vertex_index)%type
      if (bucket_type.lt.0 .or. bucket_type.gt.24) then
         write (*,*) 'Unsupported interaction type while creating amplitude library:',&
              vertex_index,bucket_type
         stop 1
      endif
      first_current=this%interaction_list(vertex_index)%currents(1)
      second_current=this%interaction_list(vertex_index)%currents(2)
      if (first_current.lt.1 .or. first_current.gt.this%n_cur .or. &
           second_current.lt.1 .or. second_current.gt.this%n_cur) then
         write (*,*) 'Out-of-range current in generated interaction:',vertex_index,&
              first_current,second_current,this%n_cur
         stop 1
      endif
      if (.not.all(ieee_is_finite(this%interaction_list(vertex_index)%coupl))) then
         write (*,*) 'Non-finite coupling in generated interaction:',vertex_index
         stop 1
      endif
      first_chirality=this%current_list(first_current)%chirality
      second_chirality=this%current_list(second_current)%chirality
      result_chirality=this%interaction_list(vertex_index)%chirality
      if (first_chirality.lt.-1 .or. first_chirality.gt.1 .or. &
           second_chirality.lt.-1 .or. second_chirality.gt.1 .or. &
           result_chirality.lt.-1 .or. result_chirality.gt.1) then
         write (*,*) 'Invalid chirality in generated interaction:',vertex_index,&
              first_chirality,second_chirality,result_chirality
         stop 1
      endif
      first_bin=this%current_list(first_current)%bin
      second_bin=this%current_list(second_current)%bin
      if (first_bin.lt.1 .or. first_bin.gt.size(this%pp_bin_to_i) .or. &
           second_bin.lt.1 .or. second_bin.gt.size(this%pp_bin_to_i)) then
         write (*,*) 'Out-of-range momentum bin in generated interaction:',&
              vertex_index,first_bin,second_bin
         stop 1
      endif
      if (this%pp_bin_to_i(first_bin).lt.1 .or. &
           this%pp_bin_to_i(first_bin).gt.this%max_pp .or. &
           this%pp_bin_to_i(second_bin).lt.1 .or. &
           this%pp_bin_to_i(second_bin).gt.this%max_pp) then
         write (*,*) 'Unmapped momentum bin in generated interaction:',&
              vertex_index,first_bin,second_bin
         stop 1
      endif

      bucket_key=4
      select case(bucket_type)
      case(4:7)
         bucket_key=result_chirality+4
      case(10,11,23,24)
         if (result_chirality.ne.0) then
            bucket_key=result_chirality+4
         elseif (bucket_type.eq.10 .or. bucket_type.eq.11) then
            bucket_key=4*(first_chirality+1)
         else
            bucket_key=4*(second_chirality+1)
         endif
      case(8,9,21,22)
         bucket_key=(first_chirality+1)*3+(second_chirality+1)
      end select
      if (bucket_key.lt.0 .or. bucket_key.gt.8) then
         write (*,*) 'Invalid generated interaction key:',vertex_index,&
              bucket_type,bucket_key
         stop 1
      endif
    end subroutine get_library_vertex_bucket

    subroutine require_generator_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (allocation_status.ne.0) then
         write (*,*) 'Could not allocate ',trim(label),': ',trim(allocation_message)
         stop 1
      endif
    end subroutine require_generator_allocation

    subroutine require_generator_io(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Could not write ',trim(label),': ',ios,trim(io_message)
         stop 1
      endif
    end subroutine require_generator_io
  end subroutine create_library
  subroutine build_helicity_filter(this,samples,include_hel,collapse_equivalent)
    ! Determine removable helicities from several normalized phase-space
    ! samples.  Equivalent helicities must agree at every non-zero sample;
    ! this avoids inferring a symmetry from an accidental one-point equality.
    implicit none
    class(amplitude_QCD),intent(in) :: this
    real(kind=8),intent(in) :: samples(:,:)
    integer,intent(out) :: include_hel(size(samples,1))
    logical,intent(in),optional :: collapse_equivalent
    integer :: ih1,ih2,iproc1,iproc2
    logical :: collapse

    if (this%n_amps.lt.1 .or. &
         this%n_amps.gt.max_amplitude_current_records .or. &
         this%nprocs.lt.1 .or. this%nprocs.gt.this%n_amps) then
       write (*,*) 'Invalid helicity-filter amplitude dimensions:',&
            this%n_amps,this%nprocs
       stop 1
    endif
    if (size(samples,1).ne.this%n_amps .or. size(samples,2).lt.2) then
       write (*,*) 'Invalid helicity-filter sample dimensions',shape(samples),this%n_amps
       stop 1
    endif
    if (.not.allocated(this%iproc_start)) then
       write (*,*) 'Helicity-filter subprocess offsets are missing'
       stop 1
    endif
    if (lbound(this%iproc_start,1).ne.1 .or. &
         ubound(this%iproc_start,1).ne.this%nprocs+1) then
       write (*,*) 'Invalid helicity-filter subprocess-offset bounds:',&
            lbound(this%iproc_start,1),ubound(this%iproc_start,1),this%nprocs
       stop 1
    endif
    if (this%iproc_start(1).ne.1 .or. &
         this%iproc_start(this%nprocs+1).ne.this%n_amps+1 .or. &
         any(this%iproc_start(2:this%nprocs+1).lt.&
         this%iproc_start(1:this%nprocs))) then
       write (*,*) 'Invalid helicity-filter subprocess offsets:',this%iproc_start
       stop 1
    endif
    if (.not.all(ieee_is_finite(samples))) then
       write (*,*) 'Non-finite helicity sample encountered during amplitude optimisation'
       stop 1
    endif
    collapse=.true.
    if (present(collapse_equivalent)) collapse=collapse_equivalent
    include_hel=0
    do ih1=1,this%n_amps
       if (include_hel(ih1).ne.0) cycle
       if (maxval(abs(samples(ih1,:))).le.helicity_zero_tolerance) cycle
       include_hel(ih1)=1
       if (.not.collapse) cycle
       iproc1=helicity_process(ih1)
       do ih2=ih1+1,this%n_amps
          if (include_hel(ih2).ne.0) cycle
          if (maxval(abs(samples(ih2,:))).le.helicity_zero_tolerance) cycle
          iproc2=helicity_process(ih2)
          if (iproc1.ne.iproc2) cycle
          if (.not.helicity_samples_match(ih1,ih2)) cycle
          include_hel(ih2)=-ih1
          include_hel(ih1)=include_hel(ih1)+1
       enddo
    enddo
    ! Ten identically zero warm-up points are not sufficient evidence that
    ! the entire process vanishes.  Keep the unfiltered amplitude in that
    ! exceptional case.
    if (.not.any(include_hel.gt.0)) include_hel=1

  contains

    integer function helicity_process(ihel)
      integer,intent(in) :: ihel
      integer :: iproc
      helicity_process=0
      do iproc=1,this%nprocs
         if (ihel.ge.this%iproc_start(iproc) .and. &
              ihel.lt.this%iproc_start(iproc+1)) then
            helicity_process=iproc
            return
         endif
      enddo
      write (*,*) 'Could not assign helicity to a process',ihel
      stop 1
    end function helicity_process

    logical function helicity_samples_match(first,second)
      integer,intent(in) :: first,second
      integer :: isample
      real(kind=8) :: difference,scale
      real(kind=8),parameter :: equality_tolerance=1d-11
      logical :: compared
      helicity_samples_match=.false.
      compared=.false.
      do isample=1,size(samples,2)
         scale=max(abs(samples(first,isample)),abs(samples(second,isample)))
         if (scale.le.helicity_zero_tolerance) cycle
         compared=.true.
         difference=abs(samples(first,isample)/scale-&
              samples(second,isample)/scale)
         if (difference.gt.equality_tolerance) return
      enddo
      helicity_samples_match=compared
    end function helicity_samples_match

  end subroutine build_helicity_filter

  elemental logical function complex_value_is_finite(value)
    implicit none
    complex(kind=8),intent(in) :: value
    complex_value_is_finite=ieee_is_finite(real(value,kind=8)) .and. &
         ieee_is_finite(aimag(value))
  end function complex_value_is_finite

  elemental logical function complex_amplitude_is_safe(value)
    implicit none
    complex(kind=8),intent(in) :: value
    complex_amplitude_is_safe=.false.
    if (.not.complex_value_is_finite(value)) return
    if (abs(real(value,kind=8)).gt.amplitude_value_limit .or. &
         abs(aimag(value)).gt.amplitude_value_limit) return
    complex_amplitude_is_safe=.true.
  end function complex_amplitude_is_safe

  subroutine filter_helicity(this,n,nhel,include_hel)
    implicit none
    class(amplitude_qcd),intent(inout) :: this
    integer,intent(in) :: n
    integer,intent(inout) :: nhel
    integer,intent(inout),dimension(nhel) :: include_hel
    integer :: nspin,ispin,ic,iv,iamp,max_multiplicity,representative,&
         allocation_status
    integer(kind=8) :: workspace_bytes
    logical,dimension(:),allocatable :: include_current
    integer,dimension(:,:,:),allocatable :: tmp_spin
    integer,dimension(:,:),allocatable :: tmp_perm
    integer,dimension(:),allocatable :: compact_multiplicity
    character(len=256) :: allocation_message

    if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
         nhel.lt.1 .or. nhel.ne.this%n_amps .or. &
         nhel.gt.max_amplitude_current_records .or. this%n_cur.lt.1 .or. &
         this%n_cur.gt.max_amplitude_current_records) then
       write (*,*) 'Invalid helicity-filter dimensions:',n,nhel,this%n_amps,&
            this%n_cur
       stop 1
    endif
    if (.not.allocated(this%spins) .or. .not.allocated(this%perm) .or. &
         .not.allocated(this%same_flavour_sum)) then
       write (*,*) 'Amplitude metadata is incomplete for helicity filtering'
       stop 1
    endif
    if (size(this%spins,1).ne.n .or. size(this%spins,3).lt.nhel .or. &
         size(this%perm,2).lt.nhel .or. &
         size(this%same_flavour_sum,1).lt.nhel .or. &
         size(this%same_flavour_sum,2).ne.2) then
       write (*,*) 'Amplitude metadata has incompatible helicity-filter dimensions'
       stop 1
    endif
    if (any(include_hel.lt.-nhel) .or. any(include_hel.gt.nhel) .or. &
         .not.any(include_hel.gt.0)) then
       write (*,*) 'Invalid helicity-filter map:',include_hel
       stop 1
    endif
    do iamp=1,nhel
       if (include_hel(iamp).gt.0) then
          if (include_hel(iamp).ne.1+count(include_hel.eq.-iamp)) then
             write (*,*) 'Inconsistent helicity-filter multiplicity:',iamp,&
                  include_hel(iamp),1+count(include_hel.eq.-iamp)
             stop 1
          endif
       elseif (include_hel(iamp).lt.0) then
          representative=-include_hel(iamp)
          if (representative.ge.iamp .or. include_hel(representative).le.0) then
             write (*,*) 'Invalid helicity-filter representative:',iamp,&
                  representative
             stop 1
          endif
       endif
       if (this%same_flavour_sum(iamp,1).gt.0) then
          if (any(this%same_flavour_sum(iamp,1:2).lt.1) .or. &
               any(this%same_flavour_sum(iamp,1:2).gt.nhel)) then
             write (*,*) 'Invalid same-flavour map during helicity filtering:',&
                  iamp,this%same_flavour_sum(iamp,1:2)
             stop 1
          endif
       endif
    enddo
    max_multiplicity=maxval(include_hel,mask=include_hel.gt.0)
    workspace_bytes=int(this%n_cur,kind=8)*&
         max(1_8,int(storage_size(.true.),kind=8)/8_8)+&
         4_8*int(nhel,kind=8)+&
         4_8*int(n,kind=8)*int(max_multiplicity,kind=8)*int(nhel,kind=8)+&
         4_8*int(size(this%perm,1),kind=8)*int(nhel,kind=8)
    if (workspace_bytes.lt.0_8 .or. &
         workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Helicity filtering exceeds the supported workspace:',&
            workspace_bytes,max_amplitude_workspace_bytes
       stop 1
    endif
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
    this%evaluation_workspace_ready=.false.
    
    allocate(include_current(this%n_cur),&
         tmp_spin(1:n,1:max_multiplicity,1:nhel),&
         tmp_perm(1:size(this%perm,1),1:nhel),&
         compact_multiplicity(1:nhel),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate helicity-filter workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    include_current=.false.

    if (.not.allocated(this%include_amp)) then
       allocate(this%include_amp(1:this%n_amps),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate helicity inclusion map: ',&
               trim(allocation_message)
          stop 1
       endif
    elseif (size(this%include_amp).lt.this%n_amps) then
       write (*,*) 'Helicity inclusion map is too short'
       stop 1
    endif
    this%include_amp(1:this%n_amps)=.false.
    tmp_spin=0
    tmp_perm=0
    compact_multiplicity=0

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
          tmp_perm(:,nspin)=this%perm(:,iamp)
          ic=1
          do ispin=iamp+1,nhel
             if (-include_hel(ispin).eq.iamp) then
                ic=ic+1
                tmp_spin(1:n,ic,nspin)=this%spins(1:n,1,ispin)
             endif
          enddo
          compact_multiplicity(nspin)=include_hel(iamp)
       endif
    enddo

    call this%filter_dead_trees(n,include_current)
    if (this%n_amps.ne.nspin) then
       write (*,*) 'Helicity filtering retained an inconsistent amplitude count',&
            this%n_amps,nspin
       stop 1
    endif
    deallocate(this%spins)
    allocate(this%spins(1:n,1:max_multiplicity,1:nspin),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate filtered helicities: ',&
            trim(allocation_message)
       stop 1
    endif
    this%spins=tmp_spin(1:n,1:max_multiplicity,1:nspin)
    deallocate(this%perm)
    allocate(this%perm(1:size(tmp_perm,1),1:nspin),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate filtered colour permutations: ',&
            trim(allocation_message)
       stop 1
    endif
    this%perm=tmp_perm(:,1:nspin)
    include_hel=0
    include_hel(1:nspin)=compact_multiplicity(1:nspin)
    nhel=this%n_amps
    this%include_amp(1:this%n_amps)=.true.
    write (99,*) 'Total number of currents, vertices and amplitudes after filtering helicities',this%n_cur,this%n_vert,this%n_amps

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
    class(amplitude_qcd),intent(inout) :: this
    logical,dimension(:),allocatable :: is_needed_cur,is_needed_ver
    integer,dimension(:),allocatable :: where_to_cur,where_to_ver,where_to_amp
    logical,dimension(:),intent(in),optional :: include_current
    integer,intent(in) :: n
    integer :: to_skip,isize,nc,iv,iamp,iproc,i,allocation_status
    integer(kind=8) :: workspace_bytes
    character(len=256) :: allocation_message

    if (n.lt.3 .or. n.gt.max_amplitude_external_particles .or. &
         this%n_cur.lt.1 .or. this%n_cur.gt.max_amplitude_current_records .or. &
         this%n_vert.lt.0 .or. this%n_vert.gt.max_amplitude_interaction_records .or. &
         this%n_amps.lt.1 .or. this%n_amps.gt.max_amplitude_current_records) then
       write (*,*) 'Invalid amplitude dimensions for dead-tree filtering:',&
            n,this%n_cur,this%n_vert,this%n_amps
       stop 1
    endif
    if (.not.allocated(this%current_list) .or. &
         .not.allocated(this%interaction_list) .or. &
         .not.allocated(this%n_cur_start) .or. .not.allocated(this%n_cur_end) .or. &
         .not.allocated(this%n_vert_start) .or. .not.allocated(this%n_vert_end) .or. &
         .not.allocated(this%include_amp)) then
       write (*,*) 'Amplitude metadata is incomplete for dead-tree filtering'
       stop 1
    endif
    if (size(this%current_list).lt.this%n_cur .or. &
         size(this%interaction_list).lt.this%n_vert .or. &
         size(this%n_cur_start).ne.n .or. size(this%n_cur_end).ne.n .or. &
         size(this%n_vert_start).ne.n-2 .or. size(this%n_vert_end).ne.n-2 .or. &
         size(this%include_amp).lt.this%n_amps) then
       write (*,*) 'Amplitude metadata has incompatible dead-tree dimensions'
       stop 1
    endif
    if (present(include_current)) then
       if (size(include_current).lt.this%n_cur) then
          write (*,*) 'Current-inclusion map is too short:',size(include_current),&
               this%n_cur
          stop 1
       endif
    endif
    workspace_bytes=int(this%n_cur+this%n_vert,kind=8)*&
         max(1_8,int(storage_size(.true.),kind=8)/8_8)+&
         4_8*(int(this%n_cur,kind=8)+int(this%n_vert,kind=8)+&
         int(this%n_amps,kind=8)+1_8)
    if (workspace_bytes.lt.0_8 .or. &
         workspace_bytes.gt.max_amplitude_workspace_bytes) then
       write (*,*) 'Dead-tree filtering exceeds the supported workspace:',&
            workspace_bytes,max_amplitude_workspace_bytes
       stop 1
    endif
    allocate(is_needed_cur(this%n_cur),is_needed_ver(this%n_vert),&
         where_to_cur(this%n_cur),where_to_ver(this%n_vert),&
         where_to_amp(0:this%n_amps),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate dead-tree filtering workspace: ',&
            trim(allocation_message)
       stop 1
    endif
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
       do i=1,2
          if (this%curr2amp(i,iamp).ne.0) then
             this%curr2amp(i,where_to_amp(iamp))=where_to_cur(this%curr2amp(i,iamp))
          endif
       enddo
    enddo
    if (allocated(this%three_line_partner_curr2amp)) then
       do iamp=1,this%n_amps
          do i=1,2
             if (this%three_line_partner_curr2amp(i,iamp).ne.0) then
                this%three_line_partner_curr2amp(i,iamp)=&
                     where_to_cur(this%three_line_partner_curr2amp(i,iamp))
             endif
          enddo
       enddo
    endif
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

  subroutine assign_interaction(lhs,rhs)
    ! sets non-custom 'lhs' = 'rhs' for interactions
    use particles
    implicit none
    type(interaction),intent(inout) :: lhs
    type(interaction),intent(in) :: rhs
    integer :: val_size,allocation_status
    character(len=256) :: allocation_message
    lhs%type=rhs%type
    lhs%chirality=rhs%chirality
    lhs%currents(1:2)=rhs%currents(1:2)
    lhs%coupl(1:2)=rhs%coupl(1:2)
    if (allocated(lhs%singlet_mv)) deallocate(lhs%singlet_mv)
    if (allocated(rhs%singlet_mv)) then
       if (lbound(rhs%singlet_mv,1).ne.0 .or. &
            ubound(rhs%singlet_mv,1).lt.0) then
          write (*,*) 'Invalid interaction singlet map during assignment'
          stop 1
       endif
       if (rhs%singlet_mv(0).lt.0 .or. &
            rhs%singlet_mv(0).gt.ubound(rhs%singlet_mv,1)) then
          write (*,*) 'Invalid interaction singlet map during assignment'
          stop 1
       endif
       allocate(lhs%singlet_mv(0:rhs%singlet_mv(0)),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('interaction singlet map',&
            allocation_status,allocation_message)
       lhs%singlet_mv=rhs%singlet_mv(0:rhs%singlet_mv(0))
    endif
    if (allocated(lhs%val_c)) deallocate(lhs%val_c)
    if (allocated(rhs%val_c)) then
       val_size=size(rhs%val_c)
       allocate(lhs%val_c(1:val_size),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('interaction complex values',&
            allocation_status,allocation_message)
       lhs%val_c=rhs%val_c
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       val_size=size(rhs%val_r)
       allocate(lhs%val_r(1:val_size),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('interaction real values',&
            allocation_status,allocation_message)
       lhs%val_r=rhs%val_r
    endif
  end subroutine assign_interaction
  
  subroutine assign_current(lhs,rhs)
    ! sets non-custom 'lhs' = 'rhs' for currents
    use particles
    implicit none
    type(current),intent(inout) :: lhs
    type(current),intent(in) :: rhs
    integer :: isize,val_size,lsize,allocation_status
    character(len=256) :: allocation_message
    lhs%type=rhs%type
    lhs%bin=rhs%bin
    lhs%chirality=rhs%chirality
    if (rhs%bin.lt.0 .or. popcnt(rhs%bin).gt.max_amplitude_external_particles .or. &
         rhs%n_vert.lt.0 .or. rhs%n_vert.gt.max_amplitude_interaction_records) then
       write (*,*) 'Invalid current metadata during assignment',rhs%bin,rhs%n_vert
       stop 1
    endif
    isize=popcnt(lhs%bin)
    lhs%n_vert=rhs%n_vert
    if (allocated(lhs%iproc%bits)) deallocate(lhs%iproc%bits)
    lhs%iproc%n_bits=0
    if (rhs%iproc%n_bits.lt.0 .or. rhs%iproc%n_bits.gt.max_bitset_bits .or. &
         (rhs%iproc%n_bits.gt.0 .and. .not.allocated(rhs%iproc%bits))) then
       write (*,*) 'Invalid current process mask during assignment',rhs%iproc%n_bits
       stop 1
    endif
    if (allocated(rhs%iproc%bits)) then
       call lhs%iproc%init(rhs%iproc%n_bits)
       lhs%iproc%bits=rhs%iproc%bits
    endif
    lhs%mass=rhs%mass
    lhs%width=rhs%width
    lhs%ext_cur=rhs%ext_cur
    if (allocated(lhs%vertices)) deallocate(lhs%vertices)
    if (allocated(lhs%vertex_sign)) deallocate(lhs%vertex_sign)
    if (allocated(rhs%vertices) .and. rhs%n_vert.gt.0) then
       if (.not.allocated(rhs%vertex_sign)) then
          write (*,*) 'Missing current interaction signs during assignment'
          stop 1
       endif
       if (lbound(rhs%vertices,1).gt.1 .or. &
            ubound(rhs%vertices,1).lt.rhs%n_vert .or. &
            lbound(rhs%vertex_sign,1).gt.1 .or. &
            ubound(rhs%vertex_sign,1).lt.rhs%n_vert) then
          write (*,*) 'Invalid current interaction references during assignment'
          stop 1
       endif
       allocate(lhs%vertices(1:lhs%n_vert),&
            lhs%vertex_sign(1:lhs%n_vert),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current interaction references',&
            allocation_status,allocation_message)
       lhs%vertices(1:lhs%n_vert)=rhs%vertices(1:lhs%n_vert)
       lhs%vertex_sign(1:lhs%n_vert)=rhs%vertex_sign(1:lhs%n_vert)
    elseif (rhs%n_vert.ne.0 .or. &
         (allocated(rhs%vertices) .neqv. allocated(rhs%vertex_sign))) then
       write (*,*) 'Inconsistent current interaction references during assignment'
       stop 1
    endif
    if (allocated(lhs%order)) deallocate(lhs%order)
    if (allocated(rhs%order) .and. isize.gt.0) then
       if (lbound(rhs%order,1).gt.1 .or. ubound(rhs%order,1).lt.isize) then
          write (*,*) 'Current order is too short during assignment'
          stop 1
       endif
       allocate(lhs%order(1:isize),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current order',allocation_status,&
            allocation_message)
       lhs%order(1:isize)=rhs%order(1:isize)
    endif
    if (allocated(lhs%spin)) deallocate(lhs%spin)
    if (allocated(rhs%spin) .and. isize.gt.0) then
       if (lbound(rhs%spin,1).gt.1 .or. ubound(rhs%spin,1).lt.isize) then
          write (*,*) 'Current spin list is too short during assignment'
          stop 1
       endif
       allocate(lhs%spin(1:isize),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current spin list',allocation_status,&
            allocation_message)
       lhs%spin(1:isize)=rhs%spin(1:isize)
    endif
    if (allocated(lhs%ext_type)) deallocate(lhs%ext_type)
    if (allocated(rhs%ext_type) .and. isize.gt.0) then
       if (lbound(rhs%ext_type,1).gt.1 .or. &
            ubound(rhs%ext_type,1).lt.isize) then
          write (*,*) 'Current external-type list is too short during assignment'
          stop 1
       endif
       allocate(lhs%ext_type(1:isize),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current external-type list',&
            allocation_status,allocation_message)
       lhs%ext_type(1:isize)=rhs%ext_type(1:isize)
    endif
    if (allocated(lhs%val_c)) deallocate(lhs%val_c)
    if (allocated(rhs%val_c)) then
       val_size=size(rhs%val_c)
       allocate(lhs%val_c(1:val_size),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current complex values',&
            allocation_status,allocation_message)
       lhs%val_c=rhs%val_c
    endif
    if (allocated(lhs%val_r)) deallocate(lhs%val_r)
    if (allocated(rhs%val_r)) then
       val_size=size(rhs%val_r)
       allocate(lhs%val_r(1:val_size),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current real values',&
            allocation_status,allocation_message)
       lhs%val_r=rhs%val_r
    endif
    if (allocated(lhs%fermi_list)) deallocate(lhs%fermi_list)
    if (allocated(rhs%fermi_list)) then
       lsize=size(rhs%fermi_list)
       allocate(lhs%fermi_list(1:lsize),stat=allocation_status,&
            errmsg=allocation_message)
       call require_assignment_allocation('current Fermi list',allocation_status,&
            allocation_message)
       lhs%fermi_list=rhs%fermi_list
    endif
  end subroutine assign_current
  subroutine require_assignment_allocation(label,allocation_status,&
       allocation_message)
    implicit none
    character(len=*),intent(in) :: label
    integer,intent(in) :: allocation_status
    character(len=*),intent(in) :: allocation_message
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate ',trim(label),' during assignment: ',&
            trim(allocation_message)
       stop 1
    endif
  end subroutine require_assignment_allocation
  subroutine finalize_amplitude_QCD(amp)
    type(amplitude_QCD),intent(inout) :: amp
    call reset_amplitude_QCD(amp)
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
