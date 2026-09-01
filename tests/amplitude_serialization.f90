program amplitude_serialization_regression
  use amplitude_QCD_mod
  use particles
  use run_parameters, only: reset_run_parameters
  implicit none
  integer,parameter :: n=4
  type(physics_model) :: model
  type(amplitude_QCD) :: original,restored
  integer :: process(n,1),order(n,1),spin(0:3,n),hel(n),iunit,ipart,index
  integer,allocatable :: replacement_spins(:,:,:)
  real(kind=8) :: p(0:3,n),scale
  complex(kind=8),allocatable :: reference(:)
  character(len=64) :: mode

  call reset_run_parameters()
  mode='roundtrip'
  call get_command_argument(1,mode)
  if (len_trim(mode).eq.0) mode='roundtrip'
  open(unit=99,status='scratch',action='readwrite')
  call model%init_part()
  call model%init_vert()

  process(:,1)=[2,-2,21,21]
  order(:,1)=[2,3,4,1]
  spin=0
  do ipart=1,n
     spin(0,ipart)=model%get_spin(process(ipart,1))
     spin(1:2,ipart)=[-1,1]
  enddo
  hel=spin(1,:)
  call fill_point(p)

  call original%init(1,n,1,process,spin,order,model)
  if (.not.allocated(original%processes)) then
     write(*,*) 'Amplitude lost its process metadata during initialization'
     stop 1
  endif
  call original%evaluate(n,p,hel,.false.,model)
  if (.not.allocated(original%processes)) then
     write(*,*) 'Amplitude evaluation destroyed its process metadata'
     stop 1
  endif
  if (size(original%processes,1).ne.n .or. &
       size(original%processes,2).ne.original%nprocs) then
     write(*,*) 'Amplitude process metadata has inconsistent shape',&
          shape(original%processes),original%nprocs
     stop 1
  endif
  allocate(reference(original%n_amps))
  reference=original%amps

  open(newunit=iunit,status='scratch',action='readwrite',access='stream',&
       form='unformatted')
  select case(trim(mode))
  case('missing_current_ranges')
     deallocate(original%n_cur_start)
  case('short_current_vertices')
     do index=1,original%n_cur
        if (original%current_list(index)%n_vert.gt.0) exit
     enddo
     if (index.gt.original%n_cur) then
        write(*,*) 'Regression setup found no composite current'
        stop 90
     endif
     deallocate(original%current_list(index)%vertices)
     allocate(original%current_list(index)%vertices(1:0))
  case('malformed_singlet_map')
     do index=1,original%n_vert
        if (allocated(original%interaction_list(index)%singlet_mv)) exit
     enddo
     if (index.gt.original%n_vert) then
        write(*,*) 'Regression setup found no interaction singlet map'
        stop 90
     endif
     deallocate(original%interaction_list(index)%singlet_mv)
     allocate(original%interaction_list(index)%singlet_mv(1:1))
     original%interaction_list(index)%singlet_mv=0
  case('invalid_permutation')
     original%perm(1,1)=original%perm(2,1)
  case('invalid_process')
     original%processes(1,1)=-huge(0)-1
  case('bad_offsets')
     original%iproc_start(1)=2
  case('bad_process_mask')
     original%current_list(1)%iproc%n_bits=original%nprocs+1
  case('orphan_same_flavour_operation')
     original%same_flavour_sum(1,1)=-1
     original%same_flavour_sum_operation(1,1)=1
  case('extra_spin_multiplicity')
     allocate(replacement_spins(n,2,original%n_amps))
     replacement_spins=0
     replacement_spins(:,1,:)=original%spins(:,1,:)
     call move_alloc(replacement_spins,original%spins)
  case('roundtrip')
     continue
  case default
     write(*,*) 'Unknown amplitude-serialization regression mode: ',trim(mode)
     stop 90
  end select
  call original%write_init_amps_to_file(n,iunit)
  if (trim(mode).ne.'roundtrip') then
     write(*,*) 'Malformed amplitude state was serialized: ',trim(mode)
     stop 91
  endif
  rewind(iunit)
  ! Reading must safely replace a fully populated object, including numerical
  ! and colour workspaces left by an earlier evaluation.
  call restored%init(1,n,1,process,spin,order,model)
  call restored%evaluate(n,p,hel,.false.,model)
  call restored%read_init_amps_from_file(n,iunit)
  close(iunit)

  if (restored%n_amps.ne.original%n_amps .or. &
       restored%nprocs.ne.original%nprocs .or. &
       any(restored%processes.ne.original%processes) .or. &
       any(restored%spins.ne.original%spins) .or. &
       any(restored%perm.ne.original%perm)) then
     write(*,*) 'Serialized amplitude metadata changed during round trip'
     stop 1
  endif
  call restored%evaluate(n,p,hel,.false.,model)
  scale=max(maxval(abs(reference)),maxval(abs(restored%amps)),tiny(1d0))
  if (maxval(abs(reference-restored%amps)).gt.2d-12*scale) then
     write(*,*) 'Serialized amplitude evaluation changed during round trip',&
          maxval(abs(reference-restored%amps)),scale
     stop 1
  endif

  ! Normal initialisation must likewise be reusable on an already evaluated,
  ! deserialized object.
  call restored%init(1,n,1,process,spin,order,model)
  call restored%evaluate(n,p,hel,.false.,model)
  scale=max(maxval(abs(reference)),maxval(abs(restored%amps)),tiny(1d0))
  if (maxval(abs(reference-restored%amps)).gt.2d-12*scale) then
     write(*,*) 'Reinitialised amplitude evaluation changed',&
          maxval(abs(reference-restored%amps)),scale
     stop 1
  endif

  call check_empty_colour_order()
  call check_two_quark_line_singlet_colour_order()

  close(99)
  write(*,'(a)') 'amplitude serialization regression: PASS'

contains

  subroutine fill_point(momentum)
    real(kind=8),intent(out) :: momentum(0:3,n)
    real(kind=8),parameter :: energy=500d0,costheta=0.37d0,phi=0.61d0
    real(kind=8) :: sintheta,px,py,pz
    sintheta=sqrt(1d0-costheta*costheta)
    px=energy*sintheta*cos(phi)
    py=energy*sintheta*sin(phi)
    pz=energy*costheta
    momentum(:,1)=[energy,0d0,0d0,energy]
    momentum(:,2)=[energy,0d0,0d0,-energy]
    momentum(:,3)=[energy,px,py,pz]
    momentum(:,4)=[energy,-px,-py,-pz]
  end subroutine fill_point

  subroutine check_empty_colour_order()
    type(amplitude_QCD) :: singlet_original,singlet_restored
    integer :: singlet_process(n,1),singlet_order(n,1),singlet_spin(0:3,n)
    integer :: particle,unit,singlet_hel(n)
    complex(kind=8),allocatable :: singlet_reference(:)
    real(kind=8) :: singlet_scale

    singlet_process(:,1)=[11,-11,13,-13]
    singlet_order(:,1)=[1,2,3,4]
    singlet_spin=0
    do particle=1,n
       singlet_spin(0,particle)=model%get_spin(singlet_process(particle,1))
       singlet_spin(1:2,particle)=[-1,1]
    enddo
    call singlet_original%init(1,n,1,singlet_process,singlet_spin,&
         singlet_order,model)
    if (size(singlet_original%perm,1).ne.0) then
       write(*,*) 'All-singlet amplitude did not have an empty colour order'
       stop 1
    endif
    call singlet_original%init_col(n,2)
    if (any(singlet_original%n_col_vals.ne.1) .or. &
         any(abs(singlet_original%diff_col_vals-1d0).gt.0d0)) then
       write(*,*) 'All-singlet amplitude has a nontrivial colour matrix'
       stop 1
    endif
    singlet_hel=[-1,1,-1,1]
    call singlet_original%evaluate(n,p,singlet_hel,.false.,model)
    allocate(singlet_reference(singlet_original%n_amps))
    singlet_reference=singlet_original%amps
    open(newunit=unit,status='scratch',action='readwrite',access='stream',&
         form='unformatted')
    call singlet_original%write_init_amps_to_file(n,unit)
    rewind(unit)
    call singlet_restored%read_init_amps_from_file(n,unit)
    close(unit)
    if (size(singlet_restored%perm,1).ne.0 .or. &
         any(singlet_restored%processes.ne.singlet_process)) then
       write(*,*) 'Empty colour order changed during serialization'
       stop 1
    endif
    call singlet_restored%init_col(n,2)
    call singlet_restored%evaluate(n,p,singlet_hel,.false.,model)
    singlet_scale=max(maxval(abs(singlet_reference)),&
         maxval(abs(singlet_restored%amps)),tiny(1d0))
    if (maxval(abs(singlet_reference-singlet_restored%amps)).gt.&
         2d-12*singlet_scale) then
       write(*,*) 'All-singlet amplitude changed during serialization'
       stop 1
    endif
  end subroutine check_empty_colour_order

  subroutine check_two_quark_line_singlet_colour_order()
    integer,parameter :: n_coloured=4,n_with_singlets=6
    type(amplitude_QCD) :: coloured,singlet
    integer :: coloured_process(n_coloured,1),coloured_order(n_coloured,1)
    integer :: singlet_process(n_with_singlets,1),singlet_order(n_with_singlets,1)
    integer :: coloured_spin(0:3,n_coloured),singlet_spin(0:3,n_with_singlets)
    integer :: particle

    coloured_process(:,1)=[2,-2,1,-1]
    coloured_order(:,1)=[2,1,3,4]
    singlet_process(:,1)=[2,-2,1,-1,22,22]
    singlet_order(:,1)=[2,1,3,5,6,4]
    coloured_spin=0
    do particle=1,n_coloured
       coloured_spin(0,particle)=model%get_spin(coloured_process(particle,1))
       coloured_spin(1:2,particle)=[-1,1]
    enddo
    singlet_spin=0
    do particle=1,n_with_singlets
       singlet_spin(0,particle)=model%get_spin(singlet_process(particle,1))
       singlet_spin(1:2,particle)=[-1,1]
    enddo

    call coloured%init(1,n_coloured,1,coloured_process,coloured_spin,&
         coloured_order,model)
    call singlet%init(1,n_with_singlets,1,singlet_process,singlet_spin,&
         singlet_order,model)
    call coloured%init_col(n_coloured,2)
    call singlet%init_col(n_with_singlets,2)

    if (coloured%nColOrd.ne.singlet%nColOrd .or. &
         any(coloured%n_col_vals.ne.singlet%n_col_vals)) then
       write(*,*) 'Colour singlets changed the two-quark-line colour basis'
       stop 1
    endif
    if (any(shape(coloured%diff_col_vals).ne.shape(singlet%diff_col_vals)) .or. &
         any(shape(coloured%i_col_i).ne.shape(singlet%i_col_i)) .or. &
         any(shape(coloured%row_index).ne.shape(singlet%row_index)) .or. &
         size(coloured%col_index).ne.size(singlet%col_index)) then
       write(*,*) 'Colour singlets changed compressed colour-matrix dimensions'
       stop 1
    endif
    if (any(abs(coloured%diff_col_vals-singlet%diff_col_vals).gt.0d0) .or. &
         any(coloured%i_col_i.ne.singlet%i_col_i) .or. &
         any(coloured%row_index.ne.singlet%row_index) .or. &
         any(coloured%col_index.ne.singlet%col_index)) then
       write(*,*) 'Colour singlets changed the two-quark-line colour matrix'
       stop 1
    endif
  end subroutine check_two_quark_line_singlet_colour_order

end program amplitude_serialization_regression
