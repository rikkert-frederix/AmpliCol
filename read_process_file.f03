module read_process_file
  use handling_processes
  use run_parameters, only: ignore_final_state_width_fix,flavour_scheme,set_flavour_scheme
  implicit none
  integer,parameter :: supported_process_file_version=7
  integer :: sf_nprocs
  integer,dimension(:,:),allocatable :: unique_procs,processes,color_orders,multi_chans
  integer,dimension(:,:),allocatable :: phase_space_permutations
  integer,dimension(:,:,:),allocatable :: multi_chan_permutations
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map,iden_iproc
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: idenCOandMAPfactor
  type :: partner_set_file
     integer :: npairs=0
     integer,dimension(:),allocatable :: map_ids,permutation_ids
  end type partner_set_file
  logical :: process_heft_enabled=.false.
contains
  logical function is_supported_process_file_version(version)
    implicit none
    integer,intent(in) :: version
    is_supported_process_file_version=version.eq.supported_process_file_version
  end function is_supported_process_file_version

  subroutine parse_process_options(options_line,file_flavour_scheme)
    implicit none
    character(len=*),intent(in) :: options_line
    integer,intent(out) :: file_flavour_scheme
    integer :: key_position,value_start,value_length,value_end,ios

    if (index(adjustl(options_line),'# options:').ne.1) then
       write (*,*) 'Missing process options in version-7 process file'
       stop 1
    endif
    key_position=index(options_line,'flavour_scheme=')
    if (key_position.eq.0) then
       write (*,*) 'Process options do not define flavour_scheme'
       stop 1
    endif
    value_start=key_position+len('flavour_scheme=')
    value_length=scan(options_line(value_start:),' ')
    if (value_length.eq.0) then
       value_end=len_trim(options_line)
    else
       value_end=value_start+value_length-2
    endif
    read(options_line(value_start:value_end),*,iostat=ios) file_flavour_scheme
    if (ios.ne.0 .or. file_flavour_scheme.lt.1 .or. file_flavour_scheme.gt.5) then
       write (*,*) 'Invalid flavour_scheme in process options',&
            trim(options_line(value_start:value_end))
       stop 1
    endif
  end subroutine parse_process_options

  subroutine configure_flavour_scheme_from_process_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    character(len=65536) :: buffer
    integer :: iunit,ios,n_external,n_unique,process_file_version
    integer :: file_flavour_scheme,heft_flag

    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not open process file while reading model metadata: ',trim(filename)
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    if (ios.ne.0) then
       write (*,*) 'Could not read process-file header while reading model metadata'
       stop 1
    endif
    read(buffer,*,iostat=ios) n_external,n_unique,process_file_version,heft_flag
    if (ios.ne.0 .or. n_external.lt.3 .or. n_unique.lt.1) then
       write (*,*) 'Malformed process-file header while reading model metadata'
       stop 1
    endif
    if (.not.is_supported_process_file_version(process_file_version)) then
       write (*,*) 'Unsupported process-file version; regenerate processes.txt',&
            process_file_version
       stop 1
    endif
    if (heft_flag.ne.0 .and. heft_flag.ne.1) then
       write (*,*) 'HEFT process flag must be zero or one',heft_flag
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    if (ios.ne.0 .or. index(adjustl(buffer),'# process:').ne.1) then
       write (*,*) 'Missing process provenance in version-7 process file'
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    close(iunit)
    if (ios.ne.0) then
       write (*,*) 'Missing process options in version-7 process file'
       stop 1
    endif
    call parse_process_options(buffer,file_flavour_scheme)
    call set_flavour_scheme(file_flavour_scheme)
    process_heft_enabled=heft_flag.eq.1
    call phys_model%set_heft_enabled(process_heft_enabled)
  end subroutine configure_flavour_scheme_from_process_file

  subroutine configure_model_from_process_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    call configure_flavour_scheme_from_process_file(filename)
  end subroutine configure_model_from_process_file

  subroutine apply_final_state_widths_from_process_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    character(len=65536) :: buffer
    character(len=1) :: record_kind
    integer :: iunit,ios,n_external,n_unique,version,heft_flag,iproc
    integer :: topology_pdg,itopology_species,file_flavour_scheme
    integer,dimension(4),parameter :: topology_species=[6,23,24,25]
    logical,dimension(4) :: mapped_species
    real(kind=8),dimension(4) :: nominal_widths
    integer,dimension(:,:),allocatable :: physical_process

    call configure_model_from_process_file(filename)
    if (ignore_final_state_width_fix) return
    mapped_species=.false.
    do itopology_species=1,size(topology_species)
       nominal_widths(itopology_species)=phys_model%get_width(&
            topology_species(itopology_species))
    enddo
    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not open process file while applying final-state widths: ',trim(filename)
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    if (ios.ne.0) then
       write (*,*) 'Could not read process-file header while applying final-state widths'
       stop 1
    endif
    read(buffer,*,iostat=ios) n_external,n_unique,version,heft_flag
    if (ios.ne.0 .or. n_external.lt.3 .or. n_unique.lt.1) then
       write (*,*) 'Malformed process-file header while applying final-state widths'
       stop 1
    endif
    if (.not.is_supported_process_file_version(version)) then
       write (*,*) 'Unsupported process-file version; regenerate processes.txt',version
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    if (ios.ne.0 .or. index(adjustl(buffer),'# process:').ne.1) then
       write (*,*) 'Missing process provenance in version-7 process file'
       stop 1
    endif
    read(iunit,'(a)',iostat=ios) buffer
    if (ios.ne.0) then
       write (*,*) 'Missing process options in version-7 process file'
       stop 1
    endif
    call parse_process_options(buffer,file_flavour_scheme)
    if (file_flavour_scheme.ne.flavour_scheme) then
       write (*,*) 'Process flavour scheme changed after model initialization',&
            file_flavour_scheme,flavour_scheme
       stop 1
    endif
    do iproc=1,n_unique
       read(iunit,*,iostat=ios)
       if (ios.ne.0) then
          write (*,*) 'Could not read unique subprocess while applying final-state widths',iproc
          stop 1
       endif
    enddo
    allocate(physical_process(n_external,1))
    do
       read(iunit,'(a)',iostat=ios) buffer
       if (ios.ne.0) exit
       if (len_trim(buffer).lt.1) cycle
       record_kind=adjustl(buffer)
       if (record_kind.eq.'N') then
          read(buffer,*,iostat=ios) record_kind,topology_pdg
          if (ios.ne.0) then
             write (*,*) 'Malformed phase-map node while applying final-state widths'
             stop 1
          endif
          do itopology_species=1,size(topology_species)
             if (abs(topology_pdg).eq.topology_species(itopology_species)) then
                mapped_species(itopology_species)=.true.
             endif
          enddo
       elseif (record_kind.eq.'C') then
          read(buffer,*,iostat=ios) record_kind,physical_process(:,1)
          if (ios.ne.0) then
             write (*,*) 'Malformed coefficient row while applying final-state widths'
             stop 1
          endif
          call phys_model%apply_final_state_widths(n_external,1,physical_process)
       endif
    enddo
    close(iunit)
    ! The physics model stores one width per species.  An external occurrence
    ! can therefore only be made stable globally.  Restore the nominal width
    ! for species that also occur as mapped internal resonances; otherwise the
    ! corresponding Breit-Wigner and matrix-element propagator would be lost.
    do itopology_species=1,size(topology_species)
       if (mapped_species(itopology_species)) call phys_model%set_width(&
            topology_species(itopology_species),nominal_widths(itopology_species))
    enddo
    deallocate(physical_process)
  end subroutine apply_final_state_widths_from_process_file

  subroutine apply_final_state_widths_from_loaded_groups()
    implicit none
    integer :: igroup,iproc,iident,imap,itopology,itopology_species
    integer,dimension(4),parameter :: topology_species=[6,23,24,25]
    logical,dimension(4) :: mapped_species
    real(kind=8),dimension(4) :: nominal_widths
    integer,dimension(:,:),allocatable :: physical_processes

    if (ignore_final_state_width_fix) return
    mapped_species=.false.
    do itopology_species=1,size(topology_species)
       nominal_widths(itopology_species)=phys_model%get_width(&
            topology_species(itopology_species))
    enddo
    do imap=1,nphase_maps
       do itopology=1,phase_map_catalogue(imap)%ntopology_nodes
          do itopology_species=1,size(topology_species)
             if (abs(phase_map_catalogue(imap)%topology_pdgs(itopology)).eq.&
                  topology_species(itopology_species)) then
                mapped_species(itopology_species)=.true.
             endif
          enddo
       enddo
    enddo
    do igroup=1,ngroups
       if (.not.allocated(pgl(igroup)%iden_processes)) cycle
       allocate(physical_processes(pgl(igroup)%next,&
            sum(pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc))))
       iident=0
       do iproc=1,pgl(igroup)%nproc
          physical_processes(:,iident+1:iident+pgl(igroup)%iden_iproc(iproc))=&
               pgl(igroup)%iden_processes(:,1:pgl(igroup)%iden_iproc(iproc),iproc)
          iident=iident+pgl(igroup)%iden_iproc(iproc)
       enddo
       call phys_model%apply_final_state_widths(pgl(igroup)%next,iident,&
            physical_processes(:,1:iident))
       deallocate(physical_processes)
    enddo
    do itopology_species=1,size(topology_species)
       if (mapped_species(itopology_species)) call phys_model%set_width(&
            topology_species(itopology_species),nominal_widths(itopology_species))
    enddo
  end subroutine apply_final_state_widths_from_loaded_groups

  subroutine read_next_nonblank(iunit,buffer,ios)
    implicit none
    integer,intent(in) :: iunit
    integer,intent(out) :: ios
    character(len=*),intent(out) :: buffer
    do
       read(iunit,'(a)',iostat=ios) buffer
       if (ios.ne.0 .or. len_trim(buffer).gt.0) return
    enddo
  end subroutine read_next_nonblank

  logical function valid_physical_topology_pdg(pdg)
    implicit none
    integer,intent(in) :: pdg
    valid_physical_topology_pdg=pdg.eq.0 .or. pdg.eq.21 .or.&
         pdg.eq.22 .or. pdg.eq.23 .or. abs(pdg).eq.24 .or.&
         pdg.eq.25 .or. (abs(pdg).ge.1 .and. abs(pdg).le.6) .or.&
         (abs(pdg).ge.11 .and. abs(pdg).le.16)
  end function valid_physical_topology_pdg

  subroutine read_processes_from_file(filename)
    use phase_space_base, only: transform_breit_wigner,&
         transform_massless_pole,transform_massive_power,&
         transform_flat_contact
    implicit none
    character(len=80),intent(in) :: filename
    character(len=65536) :: buffer
    character(len=32) :: marker
    integer :: ios,process_file_version,file_flavour_scheme,heft_flag,ndim
    integer :: iproc,imap,inode,ipermutation,iset,ifamily,iflav
    integer :: map_index,node_pdg,node_kind,node_parameter,node_mask
    integer :: node_left,node_right,permutation_index,partner_set_index
    integer :: nsets,nfamilies,nrows,max_channels,icheck,i,j,overlap
    integer :: final_mask
    real(kind=8) :: idenCOfactor
    type(partner_set_file),dimension(:),allocatable :: partner_sets
    integer,dimension(:),allocatable :: process,order,phase_permutation
    integer,dimension(:),allocatable :: phase_space_orders,ichans
    integer,dimension(:,:),allocatable :: channel_permutations
    integer,dimension(:,:),allocatable :: raw_processes,raw_orders
    real(kind=8),dimension(:),allocatable :: raw_factors

    open(unit=10,file=filename,status='old',action='read',iostat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not open process file: ',trim(filename)
       stop 1
    endif
    read(10,'(a)',iostat=ios) buffer
    if (ios.ne.0) then
       write (*,*) 'Could not read process-file header'
       stop 1
    endif
    read(buffer,*,iostat=ios) next,nproc_unique,process_file_version,heft_flag
    if (ios.ne.0 .or. next.lt.3 .or. nproc_unique.lt.1) then
       write (*,*) 'Malformed process-file header'
       stop 1
    endif
    if (.not.is_supported_process_file_version(process_file_version)) then
       write (*,*) 'Only version-7 process files are supported; regenerate processes.txt'
       stop 1
    endif
    if (heft_flag.ne.0 .and. heft_flag.ne.1) then
       write (*,*) 'HEFT process flag must be zero or one',heft_flag
       stop 1
    endif
    if ((heft_flag.eq.1).neqv.phys_model%heft_enabled) then
       write (*,*) 'Process-file HEFT mode does not match initialized model'
       stop 1
    endif
    read(10,'(a)',iostat=ios) buffer
    if (ios.ne.0 .or. index(adjustl(buffer),'# process:').ne.1) then
       write (*,*) 'Missing process provenance in version-7 process file'
       stop 1
    endif
    read(10,'(a)',iostat=ios) buffer
    if (ios.ne.0) then
       write (*,*) 'Missing process options in version-7 process file'
       stop 1
    endif
    call parse_process_options(buffer,file_flavour_scheme)
    if (file_flavour_scheme.ne.flavour_scheme) then
       write (*,*) 'Process flavour scheme does not match initialized model',&
            file_flavour_scheme,flavour_scheme
       stop 1
    endif

    if (next.eq.3) then
       if (.not.include_pdf) then
          write (*,*) 'One-body cross sections require include_pdf=.true.'
          stop 1
       endif
       ndim=1
    else
       ndim=3*(next-2)-4
       if (include_pdf) ndim=ndim+2
    endif
    allocate(unique_procs(next,nproc_unique))
    do iproc=1,nproc_unique
       read(10,*,iostat=ios) unique_procs(:,iproc)
       if (ios.ne.0) then
          write (*,*) 'Malformed unique-process catalogue',iproc
          stop 1
       endif
    enddo
    allocate(pgl_unique)
    pgl_unique%next=next
    pgl_unique%ndim=ndim
    call check_unique_processes()

    call read_next_nonblank(10,buffer,ios)
    read(buffer,*,iostat=ios) marker,nphase_permutations
    if (ios.ne.0 .or. trim(marker).ne.'PERMUTATIONS' .or.&
         nphase_permutations.lt.1) then
       write (*,*) 'Missing version-7 permutation catalogue'
       stop 1
    endif
    allocate(phase_permutation_catalogue(next,nphase_permutations))
    do ipermutation=1,nphase_permutations
       call read_next_nonblank(10,buffer,ios)
       read(buffer,*,iostat=ios) marker,permutation_index,&
            phase_permutation_catalogue(:,ipermutation)
       if (ios.ne.0 .or. trim(marker).ne.'P' .or.&
            permutation_index.ne.ipermutation) then
          write (*,*) 'Malformed permutation-catalogue entry',ipermutation
          stop 1
       endif
       if (phase_permutation_catalogue(1,ipermutation).ne.1 .or.&
            phase_permutation_catalogue(2,ipermutation).ne.2) then
          write (*,*) 'Phase-map permutations must fix incoming legs'
          stop 1
       endif
       do i=1,next
          if (count(phase_permutation_catalogue(:,ipermutation).eq.i).ne.1) then
             write (*,*) 'Invalid phase-map permutation',ipermutation
             stop 1
          endif
       enddo
    enddo

    call read_next_nonblank(10,buffer,ios)
    read(buffer,*,iostat=ios) marker,nphase_maps
    if (ios.ne.0 .or. trim(marker).ne.'PHASE_MAPS' .or.&
         nphase_maps.lt.1) then
       write (*,*) 'Missing version-7 phase-map catalogue'
       stop 1
    endif
    allocate(phase_map_catalogue(nphase_maps))
    final_mask=ibclr(ibclr(2**next-1,0),1)
    do imap=1,nphase_maps
       call read_next_nonblank(10,buffer,ios)
       allocate(phase_map_catalogue(imap)%order(next))
       read(buffer,*,iostat=ios) marker,map_index,&
            phase_map_catalogue(imap)%ntopology_nodes,&
            phase_map_catalogue(imap)%order
       if (ios.ne.0 .or. trim(marker).ne.'M' .or. map_index.ne.imap .or.&
            phase_map_catalogue(imap)%ntopology_nodes.lt.0) then
          write (*,*) 'Malformed phase-map recipe',imap
          stop 1
       endif
       do i=1,next
          if (count(phase_map_catalogue(imap)%order.eq.i).ne.1) then
             write (*,*) 'Invalid phase-map production order',imap
             stop 1
          endif
       enddo
       if (phase_map_catalogue(imap)%order(1).ne.1) then
          write (*,*) 'Phase-map production order must start at incoming leg 1'
          stop 1
       endif
       i=phase_map_catalogue(imap)%ntopology_nodes
       allocate(phase_map_catalogue(imap)%topology_pdgs(i))
       allocate(phase_map_catalogue(imap)%topology_kinds(i))
       allocate(phase_map_catalogue(imap)%topology_parameters(i))
       allocate(phase_map_catalogue(imap)%topology_masks(i))
       allocate(phase_map_catalogue(imap)%topology_left_masks(i))
       do inode=1,i
          call read_next_nonblank(10,buffer,ios)
          read(buffer,*,iostat=ios) marker,node_pdg,node_kind,node_parameter,&
               node_mask,node_left,node_right
          if (ios.ne.0 .or. trim(marker).ne.'N') then
             write (*,*) 'Malformed phase-map node',imap,inode
             stop 1
          endif
          if (.not.valid_physical_topology_pdg(node_pdg)) then
             write (*,*) 'Auxiliary or unknown field in physical phase map',node_pdg
             stop 1
          endif
          if (node_kind.lt.transform_breit_wigner .or.&
               node_kind.gt.transform_flat_contact) then
             write (*,*) 'Unknown phase-map transform kind',node_kind
             stop 1
          endif
          if ((node_pdg.eq.0) .neqv.&
               (node_kind.eq.transform_flat_contact)) then
             write (*,*) 'Flat contacts must use PDG zero exclusively',node_pdg,node_kind
             stop 1
          endif
          select case (node_kind)
          case (transform_breit_wigner)
             if (abs(node_pdg).ne.node_parameter .or.&
                  all(node_parameter.ne.[6,23,24,25])) then
                write (*,*) 'Invalid Breit-Wigner transform parameter',&
                     node_pdg,node_parameter
                stop 1
             endif
          case (transform_massless_pole)
             if (node_parameter.ne.0 .or.&
                  .not.(node_pdg.eq.21 .or. node_pdg.eq.22 .or.&
                  (abs(node_pdg).ge.11 .and. abs(node_pdg).le.16) .or.&
                  (abs(node_pdg).ge.1 .and.&
                  abs(node_pdg).le.flavour_scheme))) then
                write (*,*) 'Invalid massless-pole transform parameter',&
                     node_pdg,node_parameter
                stop 1
             endif
          case (transform_massive_power)
             if (abs(node_pdg).ne.node_parameter .or.&
                  node_parameter.le.flavour_scheme .or.&
                  node_parameter.gt.6) then
                write (*,*) 'Invalid massive-power transform parameter',&
                     node_pdg,node_parameter
                stop 1
             endif
          case (transform_flat_contact)
             if (node_parameter.ne.0) then
                write (*,*) 'Flat contacts cannot have a transform parameter'
                stop 1
             endif
          end select
          if (node_mask.le.0 .or. iand(node_mask,final_mask).ne.node_mask .or.&
               node_left.le.0 .or. node_right.le.0 .or.&
               iand(node_left,node_right).ne.0 .or.&
               ior(node_left,node_right).ne.node_mask) then
             write (*,*) 'Invalid phase-map node masks',imap,inode
             stop 1
          endif
          do j=1,inode-1
             overlap=iand(node_mask,&
                  phase_map_catalogue(imap)%topology_masks(j))
             if (node_mask.eq.phase_map_catalogue(imap)%topology_masks(j) .or.&
                  (overlap.ne.0 .and. overlap.ne.node_mask .and.&
                  overlap.ne.phase_map_catalogue(imap)%topology_masks(j))) then
                write (*,*) 'Phase-map topology must be distinct and laminar',imap
                stop 1
             endif
          enddo
          phase_map_catalogue(imap)%topology_pdgs(inode)=node_pdg
          phase_map_catalogue(imap)%topology_kinds(inode)=node_kind
          phase_map_catalogue(imap)%topology_parameters(inode)=node_parameter
          phase_map_catalogue(imap)%topology_masks(inode)=node_mask
          phase_map_catalogue(imap)%topology_left_masks(inode)=node_left
       enddo
    enddo

    call read_next_nonblank(10,buffer,ios)
    read(buffer,*,iostat=ios) marker,nsets
    if (ios.ne.0 .or. trim(marker).ne.'PARTNER_SETS' .or. nsets.lt.1) then
       write (*,*) 'Missing version-7 partner-set catalogue'
       stop 1
    endif
    allocate(partner_sets(nsets))
    do iset=1,nsets
       call read_next_nonblank(10,buffer,ios)
       read(buffer,*,iostat=ios) marker,partner_set_index,&
            partner_sets(iset)%npairs
       if (ios.ne.0 .or. trim(marker).ne.'S' .or.&
            partner_set_index.ne.iset .or. partner_sets(iset)%npairs.lt.1) then
          write (*,*) 'Malformed partner-set entry',iset
          stop 1
       endif
       allocate(partner_sets(iset)%map_ids(partner_sets(iset)%npairs))
       allocate(partner_sets(iset)%permutation_ids(partner_sets(iset)%npairs))
       read(buffer,*,iostat=ios) marker,partner_set_index,&
            partner_sets(iset)%npairs,&
            (partner_sets(iset)%map_ids(i),&
             partner_sets(iset)%permutation_ids(i),&
             i=1,partner_sets(iset)%npairs)
       if (ios.ne.0 .or. any(partner_sets(iset)%map_ids.lt.1) .or.&
            any(partner_sets(iset)%map_ids.gt.nphase_maps) .or.&
            any(partner_sets(iset)%permutation_ids.lt.1) .or.&
            any(partner_sets(iset)%permutation_ids.gt.nphase_permutations)) then
          write (*,*) 'Partner-set entry references an unknown catalogue item',iset
          stop 1
       endif
       do i=2,partner_sets(iset)%npairs
          if (partner_sets(iset)%map_ids(i).lt.&
               partner_sets(iset)%map_ids(i-1) .or.&
               (partner_sets(iset)%map_ids(i).eq.&
               partner_sets(iset)%map_ids(i-1) .and.&
               partner_sets(iset)%permutation_ids(i).le.&
               partner_sets(iset)%permutation_ids(i-1))) then
             write (*,*) 'Partner-set entries must be ordered and unique',iset
             stop 1
          endif
       enddo
    enddo

    call read_next_nonblank(10,buffer,ios)
    read(buffer,*,iostat=ios) marker,nfamilies
    if (ios.ne.0 .or. trim(marker).ne.'INTEGRATION_FAMILIES' .or.&
         nfamilies.lt.1 .or. nfamilies.ne.nsets) then
       write (*,*) 'Malformed version-7 integration-family catalogue'
       stop 1
    endif
    ngroups=nfamilies
    allocate(pgl(ngroups))
    allocate(process(next))
    allocate(order(next))
    allocate(phase_permutation(next))
    allocate(phase_space_orders(next))

    do ifamily=1,nfamilies
       call read_next_nonblank(10,buffer,ios)
       read(buffer,*,iostat=ios) marker,icheck,partner_set_index,nrows
       if (ios.ne.0 .or. trim(marker).ne.'F' .or. icheck.ne.ifamily .or.&
            partner_set_index.ne.ifamily .or. nrows.lt.1) then
          write (*,*) 'Malformed integration-family entry',ifamily
          stop 1
       endif
       max_channels=partner_sets(partner_set_index)%npairs
       pgl(ifamily)%nsubmaps=max_channels
       allocate(pgl(ifamily)%phase_maps(max_channels))
       allocate(ichans(0:max_channels))
       allocate(channel_permutations(next,max_channels))
       ichans(0)=max_channels
       do i=1,max_channels
          ichans(i)=i
          map_index=partner_sets(partner_set_index)%map_ids(i)
          permutation_index=&
               partner_sets(partner_set_index)%permutation_ids(i)
          pgl(ifamily)%phase_maps(i)%recipe_id=map_index
          pgl(ifamily)%phase_maps(i)%permutation_id=permutation_index
          allocate(pgl(ifamily)%phase_maps(i)%permutation(next))
          pgl(ifamily)%phase_maps(i)%permutation=&
               phase_permutation_catalogue(:,permutation_index)
          channel_permutations(:,i)=&
               phase_permutation_catalogue(:,permutation_index)
       enddo
       phase_permutation=channel_permutations(:,1)
       phase_space_orders=phase_map_catalogue(&
            partner_sets(partner_set_index)%map_ids(1))%order

       allocate(raw_processes(next,nrows))
       allocate(raw_orders(next,nrows))
       allocate(raw_factors(nrows))
       do iproc=1,nrows
          call read_next_nonblank(10,buffer,ios)
          read(buffer,*,iostat=ios) marker,raw_processes(:,iproc),&
               raw_orders(:,iproc),raw_factors(iproc)
          if (ios.ne.0 .or. trim(marker).ne.'C') then
             write (*,*) 'Malformed coefficient row',ifamily,iproc
             stop 1
          endif
          do i=1,next
             if (count(raw_orders(:,iproc).eq.i).ne.1) then
                write (*,*) 'Invalid colour order in coefficient row',ifamily,iproc
                stop 1
             endif
          enddo
       enddo

       nprocs=0
       sf_nprocs=0
       allocate(iden_iproc(nrows))
       allocate(processes(next,nrows))
       allocate(color_orders(next,nrows))
       allocate(phase_space_permutations(next,nrows))
       phase_space_permutations=0
       allocate(iden_processes(next,nrows,nrows))
       allocate(idenCOandMAPfactor(nrows,nrows))
       allocate(multi_chans(0:max_channels,nrows))
       allocate(multi_chan_permutations(next,max_channels,nrows))
       multi_chan_permutations=0
       do iflav=1,2
          do iproc=1,nrows
             process=raw_processes(:,iproc)
             order=raw_orders(:,iproc)
             idenCOfactor=raw_factors(iproc)
             call add_to_process_list(process,order,phase_permutation,&
                  channel_permutations,idenCOfactor,max_channels,ichans,&
                  iflav.eq.1)
          enddo
          if (.not.decompose_same_flavour_into_two_diff_flavour) exit
       enddo
       if (nprocs.lt.1) then
          write (*,*) 'Integration family has no active coefficient rows',ifamily
          stop 1
       endif

       pgl(ifamily)%next=next
       pgl(ifamily)%nproc=nprocs
       pgl(ifamily)%ndim=ndim
       pgl(ifamily)%multichan%max_channels=max_channels
       if (keep_processes_separate) then
          allocate(pgl(ifamily)%amps(nprocs))
          allocate(pgl(ifamily)%nhel(nprocs))
          allocate(pgl(ifamily)%passed(nprocs))
       else
          allocate(pgl(ifamily)%amps(1))
          allocate(pgl(ifamily)%nhel(1))
          allocate(pgl(ifamily)%passed(1))
       endif
       allocate(pgl(ifamily)%processes(next,nprocs))
       allocate(pgl(ifamily)%color_orders(next,nprocs))
       allocate(pgl(ifamily)%phase_space_permutations(next,nprocs))
       allocate(pgl(ifamily)%phase_space_orders(next))
       allocate(pgl(ifamily)%idenCOandMAPfactor(&
            maxval(iden_iproc(1:nprocs)),nprocs))
       allocate(pgl(ifamily)%iden_iproc(nprocs))
       allocate(pgl(ifamily)%iden_processes(next,&
            maxval(iden_iproc(1:nprocs)),nprocs))
       allocate(pgl(ifamily)%val_procs(&
            maxval(iden_iproc(1:nprocs)),nprocs))
       allocate(pgl(ifamily)%multichan%channels(max_channels,nprocs))
       allocate(pgl(ifamily)%multichan%channel_permutations(&
            next,max_channels,nprocs))
       allocate(pgl(ifamily)%multichan%number_of_channels(nprocs))
       pgl(ifamily)%processes=processes(:,1:nprocs)
       pgl(ifamily)%color_orders=color_orders(:,1:nprocs)
       pgl(ifamily)%phase_space_permutations=&
            phase_space_permutations(:,1:nprocs)
       pgl(ifamily)%phase_space_orders=phase_space_orders
       pgl(ifamily)%idenCOandMAPfactor=&
            idenCOandMAPfactor(1:maxval(iden_iproc(1:nprocs)),1:nprocs)
       pgl(ifamily)%iden_iproc=iden_iproc(1:nprocs)
       pgl(ifamily)%iden_processes=&
            iden_processes(:,1:maxval(iden_iproc(1:nprocs)),1:nprocs)
       pgl(ifamily)%multichan%channels=&
            multi_chans(1:max_channels,1:nprocs)
       pgl(ifamily)%multichan%channel_permutations=&
            multi_chan_permutations(:,:,1:nprocs)
       pgl(ifamily)%multichan%number_of_channels=&
            multi_chans(0,1:nprocs)
       pgl(ifamily)%passed=0
       pgl(ifamily)%ntopology_nodes=0

       deallocate(raw_processes,raw_orders,raw_factors)
       deallocate(iden_iproc,processes,color_orders)
       deallocate(phase_space_permutations,iden_processes)
       deallocate(idenCOandMAPfactor,multi_chans)
       deallocate(multi_chan_permutations,ichans,channel_permutations)
    enddo

    call read_next_nonblank(10,buffer,ios)
    if (ios.ne.0 .or. trim(buffer).ne.'END_PROCESSES') then
       write (*,*) 'Missing END_PROCESSES marker in version-7 process file'
       stop 1
    endif
    close(10)
    deallocate(process,order,phase_permutation,phase_space_orders)
    deallocate(partner_sets)
  end subroutine read_processes_from_file

  subroutine read_processes_from_file_legacy(filename)
    implicit none
    character(len=80) :: filename
    integer :: iproc,igroup,icheck,nproc_in_group,max_channels,iflav,ndim,&
         process_file_version,heft_flag
    integer :: file_flavour_scheme
    integer :: ntopology_nodes,itopology,topology_pdg,nlabels,nleft,label,mask,&
         left_mask,overlap
    real(kind=8) :: idenCOfactor
    integer,dimension(:),allocatable :: process,order,ichans,phase_space_orders,phase_permutation
    integer,dimension(:),allocatable :: topology_labels,left_labels
    integer,dimension(:,:),allocatable :: channel_permutations
    character(len=65536) :: buff
    integer :: i,j,ios
    open(unit=10,file=filename,status='old')
    read(10,'(a)',iostat=ios) buff
    if (ios.ne.0) then
       write (*,*) 'Could not read the process-file header'
       stop 1
    endif
    read(buff,*,iostat=ios) next,nproc_unique,process_file_version,heft_flag
    if (ios.ne.0) then
       write (*,*) 'Malformed process-file header'
       stop 1
    endif
    if (.not.is_supported_process_file_version(process_file_version)) then
       write (*,*) 'Unsupported process-file version; regenerate processes.txt',&
            process_file_version
       stop 1
    endif
    if (next.lt.3 .or. nproc_unique.lt.1) then
       write (*,*) 'Invalid process-file dimensions',next,nproc_unique
       stop 1
    endif
    read(10,'(a)',iostat=ios) buff
    if (ios.ne.0 .or. index(adjustl(buff),'# process:').ne.1) then
       write (*,*) 'Missing process provenance in version-7 process file'
       stop 1
    endif
    read(10,'(a)',iostat=ios) buff
    if (ios.ne.0) then
       write (*,*) 'Missing process options in version-7 process file'
       stop 1
    endif
    call parse_process_options(buff,file_flavour_scheme)
    if (file_flavour_scheme.ne.flavour_scheme) then
       write (*,*) 'Process flavour scheme does not match initialized model',&
            file_flavour_scheme,flavour_scheme
       stop 1
    endif
    if (heft_flag.ne.0 .and. heft_flag.ne.1) then
       write (*,*) 'HEFT process flag must be zero or one',heft_flag
       stop 1
    endif
    process_heft_enabled=heft_flag.eq.1
    call phys_model%set_heft_enabled(process_heft_enabled)
    if (next.eq.3) then
       if (.not.include_pdf) then
          write (*,*) 'One-body cross sections require include_pdf=.true.'
          stop 1
       endif
       ! The partonic invariant mass is fixed by the final-state mass.  Only
       ! the boost rapidity remains to be integrated.
       ndim=1
    else
       ndim=3*(next-2)-4
       if (include_pdf) ndim=ndim+2
    endif
    allocate(unique_procs(1:next,1:nproc_unique))
    do iproc=1,nproc_unique
       read(10,*) unique_procs(1:next,iproc)
    enddo
    allocate(pgl_unique)
    pgl_unique%next=next
    pgl_unique%ndim=ndim
    call check_unique_processes()
    read(10,*)
    read(10,*)
    read (10,*) ngroups
    allocate(pgl(ngroups))

    allocate(process(1:next))
    allocate(order(1:next))
    allocate(phase_permutation(1:next))

    read (10,*) 
    do igroup=1,ngroups
       nprocs=0
       sf_nprocs=0
       allocate(phase_space_orders(1:next))
       read(10,*,iostat=ios) icheck,nproc_in_group,max_channels,phase_space_orders(1:next),ntopology_nodes
       if (ios.ne.0 .or. nproc_in_group.lt.1 .or. max_channels.lt.1 .or.&
            ntopology_nodes.lt.0) then
          write (*,*) 'Malformed version-4 phase-space group header',igroup
          stop 1
       endif
       if (icheck.ne.igroup) then
          write (*,*) 'ERROR in processes file',icheck,igroup
          stop 1
       endif
       pgl(igroup)%ntopology_nodes=ntopology_nodes
       allocate(pgl(igroup)%topology_pdgs(ntopology_nodes))
       allocate(pgl(igroup)%topology_masks(ntopology_nodes))
       allocate(pgl(igroup)%topology_left_masks(ntopology_nodes))
       pgl(igroup)%topology_left_masks=0
       do itopology=1,ntopology_nodes
          read(10,'(a)',iostat=ios) buff
          if (ios.ne.0) then
             write (*,*) 'Could not read topology record',igroup,itopology
             stop 1
          endif
          read(buff,*,iostat=ios) topology_pdg,nlabels
          if (ios.ne.0 .or. nlabels.lt.2 .or. nlabels.gt.next-2) then
             write (*,*) 'Invalid topology descendant count',igroup,itopology
             stop 1
          endif
          if (process_file_version.eq.4) then
             if (.not.(topology_pdg.eq.23 .or. abs(topology_pdg).eq.24 .or.&
                  topology_pdg.eq.25 .or. abs(topology_pdg).eq.6)) then
                write (*,*) 'Unsupported mapped resonance PDG',topology_pdg
                stop 1
             endif
          elseif (.not.((abs(topology_pdg).ge.1 .and. abs(topology_pdg).le.6) .or.&
               (abs(topology_pdg).ge.11 .and. abs(topology_pdg).le.16) .or.&
               topology_pdg.eq.21 .or. topology_pdg.eq.-21 .or.&
               topology_pdg.eq.22 .or. topology_pdg.eq.23 .or.&
               topology_pdg.eq.-23 .or. abs(topology_pdg).eq.24 .or.&
               topology_pdg.eq.25 .or. abs(topology_pdg).eq.26 .or.&
               (topology_pdg.ge.125 .and. topology_pdg.le.127))) then
             write (*,*) 'Unsupported topology-node PDG',topology_pdg
             stop 1
          endif
          allocate(topology_labels(nlabels))
          if (process_file_version.eq.4) then
             read(buff,*,iostat=ios) topology_pdg,nlabels,topology_labels
          else
             read(buff,*,iostat=ios) topology_pdg,nlabels,topology_labels,nleft
             if (ios.eq.0 .and. nleft.ge.1 .and. nleft.lt.nlabels) then
                allocate(left_labels(nleft))
                read(buff,*,iostat=ios) topology_pdg,nlabels,topology_labels,&
                     nleft,left_labels
             else
                ios=1
             endif
          endif
          if (ios.ne.0 .or. any(topology_labels.lt.3) .or.&
               any(topology_labels.gt.next)) then
             write (*,*) 'Invalid topology descendant labels',igroup,itopology
             stop 1
          endif
          mask=0
          do i=1,nlabels
             label=topology_labels(i)
             if (btest(mask,label-1)) then
                write (*,*) 'Repeated topology descendant label',igroup,itopology,label
                stop 1
             endif
             mask=ibset(mask,label-1)
          enddo
          left_mask=0
          if (process_file_version.ge.5) then
             if (any(left_labels.lt.3) .or. any(left_labels.gt.next)) then
                write (*,*) 'Invalid topology left-child labels',igroup,itopology
                stop 1
             endif
             do i=1,nleft
                label=left_labels(i)
                if (.not.btest(mask,label-1) .or. btest(left_mask,label-1)) then
                   write (*,*) 'Invalid topology left-child split',igroup,itopology,label
                   stop 1
                endif
                left_mask=ibset(left_mask,label-1)
             enddo
             if (left_mask.eq.0 .or. left_mask.eq.mask) then
                write (*,*) 'Topology child must be a nonempty proper subset',igroup,itopology
                stop 1
             endif
          endif
          do i=1,itopology-1
             overlap=iand(mask,pgl(igroup)%topology_masks(i))
             if (mask.eq.pgl(igroup)%topology_masks(i) .or.&
                  (overlap.ne.0 .and. overlap.ne.mask .and.&
                  overlap.ne.pgl(igroup)%topology_masks(i))) then
                write (*,*) 'Topology descendants must be distinct and laminar',igroup,itopology
                stop 1
             endif
          enddo
          pgl(igroup)%topology_pdgs(itopology)=topology_pdg
          pgl(igroup)%topology_masks(itopology)=mask
          pgl(igroup)%topology_left_masks(itopology)=left_mask
          deallocate(topology_labels)
          if (allocated(left_labels)) deallocate(left_labels)
       enddo
       allocate(iden_iproc(nproc_in_group))
       allocate(processes(1:next,nproc_in_group))
       allocate(color_orders(1:next,nproc_in_group))
       allocate(phase_space_permutations(1:next,nproc_in_group))
       phase_space_permutations=0
       allocate(iden_processes(1:next,nproc_in_group,nproc_in_group))
       allocate(idenCOandMAPfactor(nproc_in_group,nproc_in_group))
       allocate(multi_chans(0:max_channels,nproc_in_group))
       allocate(multi_chan_permutations(1:next,1:max_channels,nproc_in_group))
       multi_chan_permutations=0
       allocate(ichans(0:max_channels))
       allocate(channel_permutations(1:next,1:max_channels))
       do iflav=1,2
          do iproc=1,nproc_in_group
             read(10,'(a)',iostat=ios) buff
             if (ios.ne.0 .or. len_trim(buff).eq.len(buff)) then
                write (*,*) 'Could not read a complete subprocess row'
                stop 1
             endif
             read(buff,*,iostat=ios) ichans(0)
             if (ios.ne.0 .or. ichans(0).lt.1 .or. ichans(0).gt.max_channels) then
                write (*,*) 'Invalid number of multichannel partners in subprocess row'
                stop 1
             endif
             read(buff,*,iostat=ios) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),&
                  idenCOfactor,phase_permutation(1:next),&
                  channel_permutations(1:next,1:ichans(0))
             if (ios.ne.0) then
                write (*,*) 'Malformed version-7 subprocess row; phase-space maps are required'
                stop 1
             endif
             if (any(ichans(1:ichans(0)).lt.1) .or. any(ichans(1:ichans(0)).gt.ngroups)) then
                write (*,*) 'Multichannel partner outside the phase-space group range'
                stop 1
             endif
             if (phase_permutation(1).ne.1 .or. phase_permutation(2).ne.2) then
                write (*,*) 'Phase-space permutation must keep both incoming legs fixed'
                stop 1
             endif
             if (any(channel_permutations(1:2,1:ichans(0)).ne.&
                  spread([1,2],2,ichans(0)))) then
                write (*,*) 'Multichannel phase-space permutations must keep incoming legs fixed'
                stop 1
             endif
             do i=1,next
                if (count(phase_permutation.eq.i).ne.1) then
                   write (*,*) 'Invalid phase-space permutation',phase_permutation
                   stop 1
                endif
                do j=1,ichans(0)
                   if (count(channel_permutations(:,j).eq.i).ne.1) then
                      write (*,*) 'Invalid multichannel phase-space permutation',&
                           channel_permutations(:,j)
                      stop 1
                   endif
                enddo
             enddo
             if (iflav.eq.1) then
                call add_to_process_list(process,order,phase_permutation,channel_permutations,&
                     idenCOfactor,max_channels,ichans,.true.)
             else
                call add_to_process_list(process,order,phase_permutation,channel_permutations,&
                     idenCOfactor,max_channels,ichans,.false.)
             endif
          enddo
          if (iflav.eq.1) then
             if (.not.decompose_same_flavour_into_two_diff_flavour) exit
             do iproc=1,nproc_in_group
                backspace(10)
             enddo
          endif
       enddo
       pgl(igroup)%next=next
       pgl(igroup)%nproc=nprocs
       pgl(igroup)%ndim=ndim
       pgl(igroup)%multichan%max_channels=max_channels
       if (keep_processes_separate) then
          allocate(pgl(igroup)%amps(1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%nhel(1:pgl(igroup)%nproc))
          allocate(pgl(igroup)%passed(1:pgl(igroup)%nproc))
       else
          allocate(pgl(igroup)%amps(1))
          allocate(pgl(igroup)%nhel(1))
          allocate(pgl(igroup)%passed(1))
       endif
       allocate(pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%phase_space_permutations(1:next,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%phase_space_orders(1:next))
       allocate(pgl(igroup)%idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%val_procs(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%channels(1:max_channels,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%channel_permutations(1:next,1:max_channels,1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%multichan%number_of_channels(1:pgl(igroup)%nproc))
       pgl(igroup)%processes(1:next,1:pgl(igroup)%nproc)=processes(1:next,1:pgl(igroup)%nproc)
       pgl(igroup)%color_orders(1:next,1:pgl(igroup)%nproc)=color_orders(1:next,1:pgl(igroup)%nproc)
       pgl(igroup)%phase_space_permutations(1:next,1:pgl(igroup)%nproc)=&
            phase_space_permutations(1:next,1:pgl(igroup)%nproc)
       pgl(igroup)%phase_space_orders(1:next)=phase_space_orders(1:next)
       pgl(igroup)%idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
            idenCOandMAPfactor(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
       pgl(igroup)%iden_iproc(1:pgl(igroup)%nproc)=iden_iproc(1:pgl(igroup)%nproc)
       pgl(igroup)%iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)=&
            iden_processes(1:next,1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc)
       pgl(igroup)%multichan%channels(1:max_channels,1:pgl(igroup)%nproc)=multi_chans(1:max_channels,1:pgl(igroup)%nproc)
       pgl(igroup)%multichan%channel_permutations(1:next,1:max_channels,1:pgl(igroup)%nproc)=&
            multi_chan_permutations(1:next,1:max_channels,1:pgl(igroup)%nproc)
       pgl(igroup)%multichan%number_of_channels(1:pgl(igroup)%nproc)=multi_chans(0,1:pgl(igroup)%nproc)
       pgl(igroup)%passed=0
       deallocate(iden_iproc)
       deallocate(processes)
       deallocate(color_orders)
       deallocate(phase_space_permutations)
       deallocate(phase_space_orders)
       deallocate(iden_processes)
       deallocate(idenCOandMAPfactor)
       deallocate(multi_chans)
       deallocate(multi_chan_permutations)
       deallocate(ichans)
       deallocate(channel_permutations)
       write (99,*) '****************************************************'
       do iproc=1,pgl(igroup)%nproc
          write(99,*) iproc,':',pgl(igroup)%processes(1:next,iproc),' ; ',&
               pgl(igroup)%color_orders(1:next,iproc),' ; ',pgl(igroup)%iden_iproc(iproc),' ; ',&
               pgl(igroup)%multichan%channels(1:pgl(igroup)%multichan%number_of_channels(iproc),iproc)
       enddo
       write (99,*) '****************************************************'
       read(10,*)
       read(10,*)
       read(10,*)
    enddo
    close(10)
  end subroutine read_processes_from_file_legacy



  subroutine check_unique_processes()
    use phase_space_gen23_mod
    use phase_space_onebody_mod
    use cuts
    implicit none
    integer :: i,j,iproc,ih
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    real(kind=8),dimension(:),allocatable :: mass,width
    real(kind=8),dimension(pgl_unique%ndim) :: x
    real(kind=8),external :: ran2
    type(psv) :: ps
    if (pgl_unique%next.eq.3) then
       allocate(phase_space_onebody :: pgl_unique%phase_space)
    else
       allocate(phase_space_gen23 :: pgl_unique%phase_space)
    endif
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(pgl_unique%amps(1))
    allocate(pgl_unique%hel(next))
    pgl_unique%hel=1
    allocate(mass(next))
    allocate(width(next))
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes(1:next,1:nproc_unique)=unique_procs(1:next,1:nproc_unique)
    do i=1,pgl_unique%next
       mass(i)=phys_model%get_mass(pgl_unique%processes(i,1))
       width(i)=phys_model%get_width(pgl_unique%processes(i,1))
       ! For the general unique-process check, use only massless particles.
       ! A 2->1 map instead needs the physical final mass to resolve its
       ! partonic delta function.
       if (pgl_unique%next.gt.3 .and. mass(i).ne.0d0) mass(i)=0d0
       if (width(i).ne.0d0) width(i)=0d0
    enddo
    call setup_spin(pgl_unique)
    call setup_color_order(pgl_unique)

    do iproc=1,pgl_unique%nproc
       do i=1,2
          pgl_unique%processes(i,iproc)=phys_model%get_antipart(pgl_unique%processes(i,iproc))
       enddo
    enddo

    ! No multi-channel needed to check: simply use the color_orders for the phase-space order
    pgl_unique%phase_space_orders(1:next)=pgl_unique%color_orders(1:next,1)

    call setup_cuts_for_each_particle(pgl_unique,0)
    call pgl_unique%phase_space%init(sqrts,next,mass,pgl_unique%phase_space_orders,&
         pgl_unique%pt_min,pgl_unique%eta_max,pgl_unique%DR_min,pgl_unique%sqrt_s_min,.false.,include_pdf,.true.)
    allocate(ps%x(1:pgl_unique%ndim+pgl_unique%phase_space%ndim_extra))
    allocate(ps%p(0:3,1:pgl_unique%next))
    
    call pgl_unique%amps(1)%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
         pgl_unique%spin,pgl_unique%color_orders,phys_model)
        
    allocate(amp2(nevent,pgl_unique%nproc))
    allocate(pgl_unique%passed(1))

    pgl_unique%passed=0
    do while (pgl_unique%passed(1).lt.nevent)
       do i=1,pgl_unique%ndim+pgl_unique%phase_space%ndim_extra
          ps%x(i)=ran2()
       enddo
       call pgl_unique%phase_space%generate_momenta(ps)
       if (ps%jac.lt.0d0) cycle
       pgl_unique%passed(1)=pgl_unique%passed(1)+1
       call pgl_unique%amps(1)%evaluate(next,ps%p,pgl_unique%hel,read_proc_from_file,phys_model)
       iproc=0
       amp2(pgl_unique%passed(1),:)=0d0
       do ih=1,pgl_unique%amps(1)%n_amps
          do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
          amp2(pgl_unique%passed(1),iproc)=amp2(pgl_unique%passed(1),iproc)+&
               dble(pgl_unique%amps(1)%amps(ih)*dconjg(pgl_unique%amps(1)%amps(ih)))
       enddo
       call find_same_flavour(pgl_unique,nevent,amp2(1,:))
    enddo
    allocate(unique_map(1:pgl_unique%nproc))
    allocate(unique_map_value(1:pgl_unique%nproc))
    call find_unique(pgl_unique,nevent,amp2,unique_map,unique_map_value)

    do iproc=1,pgl_unique%nproc
       write (99,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:pgl_unique%next,iproc),&
            ':',pgl_unique%color_orders(1:pgl_unique%next,iproc)
    enddo
    deallocate(pgl_unique%spin)
    deallocate(pgl_unique%phase_space)
    deallocate(ps%p)
    deallocate(ps%x)
    deallocate(amp2)
  end subroutine check_unique_processes

  subroutine find_unique(pgl,nevent,amp2,unique_map,unique_map_value)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: nevent
    real(kind=8),dimension(nevent,pgl%nproc),intent(in) :: amp2
    real(kind=8),dimension(pgl%nproc),intent(out) :: unique_map_value
    integer,dimension(pgl%nproc),intent(out) :: unique_map
    integer :: i,j,k
    real(kind=8),dimension(nevent) :: ratio
    real(kind=8) :: ave
    real(kind=8),parameter :: tiny=1d-6
    unique_map=-1d0
    do i=1,pgl%nproc
       if (all(amp2(1:nevent,i).eq.0d0)) then
          unique_map(i)=0
          unique_map_value(i)=0d0
          cycle
       endif
       do j=1,i-1
          if (all(amp2(1:nevent,j).eq.0d0)) cycle
          ! A numerical relation found on the deliberately massless probe
          ! point is not a valid runtime reduction if it moves a massive
          ! species to another external leg.
          do k=1,pgl%next
             if (phys_model%get_mass(pgl%processes(k,i)).ne.&
                  phys_model%get_mass(pgl%processes(k,j))) exit
          enddo
          if (k.le.pgl%next) cycle
          ratio(1:nevent)=amp2(1:nevent,i)/amp2(1:nevent,j)
          ave=sum(ratio(1:nevent))/nevent
          if (all(abs(ratio(1:nevent)/ave-1d0).lt.tiny)) then
             unique_map_value(i)=ave
             unique_map(i)=j
             exit
          endif
       enddo
       if (j.eq.i) then
          unique_map(i)=-1
          unique_map_value(i)=1d0
       endif
    enddo
  end subroutine find_unique

  subroutine add_to_process_list(process,order,phase_permutation,channel_permutations,&
       idenCOfactor,max_channels,ichans,skip_same_flavour)
    implicit none
    integer :: max_channels
    integer,dimension(0:max_channels) :: ichans
    integer,dimension(next) :: process,order,process_unique,phase_permutation
    integer,dimension(next,max_channels) :: channel_permutations
    real(kind=8) :: idenCOfactor,idenCOMAPfactor,idenMAPfactor
    integer,dimension(0:6) :: quarks
    logical :: skip_same_flavour,is_same_flavour
    call find_quarks(process,order,quarks)
    call get_unique_process_from_quarks(quarks,process,order,process_unique,idenMAPfactor,is_same_flavour)
    if (decompose_same_flavour_into_two_diff_flavour) then
       if (skip_same_flavour .and. is_same_flavour) return
       if ((.not.skip_same_flavour) .and. (.not.is_same_flavour)) return
    endif
    if (reduce_to_unique_matrix_elements) then
       idenCOMAPfactor=idenMAPfactor*idenCOfactor
    else
       idenCOMAPfactor=idenCOfactor
    endif
    if (idenCOMAPfactor.eq.0d0) return
    call add_to_unique_process_list(process,process_unique,order,phase_permutation,&
         channel_permutations,idenCOMAPfactor,max_channels,ichans)
  end subroutine add_to_process_list


  subroutine get_unique_process_from_quarks(quarks,process,order,process_unique,idenMAPfactor,is_same_flavour)
    implicit none
    integer,dimension(0:6) :: quarks
    integer,dimension(next) :: process,order,process_unique,process_mapped,mapping
    integer :: i,iproc,iq
    real(kind=8) :: idenMAPfactor
    logical :: is_same_flavour
    call map_to_canonical_form(process,process_mapped,mapping)
    do iproc=1,pgl_unique%nproc
       if (pgl_unique%amps(1)%n_qqbar(iproc)*2.ne.quarks(0)) cycle
       if (quarks(0).eq.4) then
          if (all(pgl_unique%processes(1:4,iproc).eq.-abs(quarks(1:4)))) then ! quarks are consistent 
             if (all(process_mapped(5:next).eq.pgl_unique%processes(5:next,iproc))) exit ! and the rest as well
          endif
       elseif (quarks(0).eq.6) then
          if (all(abs(pgl_unique%processes(1:6,iproc)).eq.abs(quarks(1:6)))) then ! quarks are consistent 
             if (all(process_mapped(7:next).eq.pgl_unique%processes(7:next,iproc))) exit ! and the rest as well
          endif
       else
          if (all(process_mapped(1:next).eq.pgl_unique%processes(1:next,iproc))) exit
       endif
    enddo
    if (iproc.gt.pgl_unique%nproc) then
       write (*,*) 'Process not found',quarks
       write (*,*) process
       stop 1
    endif
    process_unique(1:next)=process(1:next)
    idenMAPfactor=unique_map_value(iproc)
    if (unique_map(iproc).gt.0) then
       iq=0
       do i=1,next
          if (phys_model%is_quark(process(order(i))) .or. phys_model%is_antiquark(process(order(i)))) then
             iq=iq+1
             if (quarks(0).eq.4) then
                if (iq.eq.1 .or. iq.eq.4) then
                   process_unique(order(i))=sign(pgl_unique%processes(iq,unique_map(iproc)),quarks(iq))
                elseif (iq.eq.2) then
                   process_unique(order(i))=sign(pgl_unique%processes(3,unique_map(iproc)),quarks(3))
                elseif (iq.eq.3) then
                   process_unique(order(i))=sign(pgl_unique%processes(2,unique_map(iproc)),quarks(2))
                endif
             elseif (quarks(0).eq.6) then
                if (iq.eq.1 .or. iq.eq.6) then
                   process_unique(order(i))=sign(pgl_unique%processes(iq,unique_map(iproc)),quarks(iq))
                elseif (iq.eq.2) then
                   process_unique(order(i))=sign(pgl_unique%processes(4,unique_map(iproc)),quarks(4))
                elseif (iq.eq.3) then
                   process_unique(order(i))=sign(pgl_unique%processes(2,unique_map(iproc)),quarks(2))
                elseif (iq.eq.4) then
                   process_unique(order(i))=sign(pgl_unique%processes(5,unique_map(iproc)),quarks(5))
                elseif (iq.eq.5) then
                   process_unique(order(i))=sign(pgl_unique%processes(3,unique_map(iproc)),quarks(3))
                endif
             elseif (quarks(0).eq.2) then
                process_unique(order(i))=sign(pgl_unique%processes(iq,unique_map(iproc)),quarks(iq))
             endif
          endif
       enddo
    endif
    ! Numerical flavour reduction is constructed on a common massless phase-
    ! space point.  It is valid only when the reduced representative has the
    ! same external mass layout as the physical subprocess.  In particular,
    ! crossing a final-state top into an initial light-quark slot would feed
    ! massive wavefunctions the wrong labelled momentum.  Retain the exact
    ! subprocess in that case; light-flavour reductions remain unchanged.
    do i=1,next
       if (phys_model%get_mass(process_unique(i)).ne.&
            phys_model%get_mass(process(i))) then
          process_unique=process
          idenMAPfactor=1d0
          exit
       endif
    enddo
    if (pgl_unique%amps(1)%same_flav(iproc)) then
       is_same_flavour=.true.
    else
       is_same_flavour=.false.
    endif
  end subroutine get_unique_process_from_quarks

  subroutine add_to_unique_process_list(process,process_unique,order,phase_permutation,&
       channel_permutations,idenCOMAPfactor,max_channels,ichans)
    implicit none
    integer,intent(in) :: max_channels
    integer,dimension(0:max_channels),intent(in) :: ichans
    integer,dimension(next) :: process,process_unique,order,phase_permutation
    integer,dimension(next,max_channels) :: channel_permutations
    real(kind=8) :: idenCOMAPfactor
    integer :: iproc
    call move_colour_singlet_in_order(process,order)
    if (.not.reduce_to_unique_matrix_elements) then
       ! always add the (unaltered) process
       nprocs=nprocs+1
       processes(1:next,nprocs)=process(1:next)
       color_orders(1:next,nprocs)=order(1:next)
       phase_space_permutations(1:next,nprocs)=phase_permutation(1:next)
       iden_iproc(nprocs)=1
       iden_processes(1:next,iden_iproc(nprocs),nprocs)=process(1:next)
       idenCOandMAPfactor(iden_iproc(nprocs),nprocs)=idenCOMAPfactor
       multi_chans(0:ichans(0),nprocs)=ichans(0:ichans(0))
       multi_chan_permutations(1:next,1:ichans(0),nprocs)=&
            channel_permutations(1:next,1:ichans(0))
       return
    endif
    do iproc=1,nprocs
       if (all(process_unique(1:next).eq.processes(1:next,iproc)) &
            .and. all(order(1:next).eq.color_orders(1:next,iproc)) &
            .and. all(phase_permutation(1:next).eq.phase_space_permutations(1:next,iproc))) exit
    enddo

    if (iproc.gt.nprocs) then
       ! new matrix element to generate
       nprocs=nprocs+1
       processes(1:next,nprocs)=process_unique(1:next)
       color_orders(1:next,nprocs)=order(1:next)
       phase_space_permutations(1:next,nprocs)=phase_permutation(1:next)
       iden_iproc(nprocs)=1
       iden_processes(1:next,iden_iproc(nprocs),nprocs)=process(1:next)
       idenCOandMAPfactor(iden_iproc(nprocs),nprocs)=idenCOMAPfactor
       multi_chans(0:ichans(0),nprocs)=ichans(0:ichans(0))
       multi_chan_permutations(1:next,1:ichans(0),nprocs)=&
            channel_permutations(1:next,1:ichans(0))
    else
       ! identical to another matrix element
       iden_iproc(iproc)=iden_iproc(iproc)+1
       iden_processes(1:next,iden_iproc(iproc),iproc)=process(1:next)
       idenCOandMAPfactor(iden_iproc(iproc),iproc)=idenCOMAPfactor
       if (ichans(0).ne.multi_chans(0,iproc)) then
          write (*,*) 'Number of multi-channels not the same among identical processes',&
               ichans(0),multi_chans(0,iproc)
          stop 1
       endif
       if (any(ichans(1:ichans(0)).ne.multi_chans(1:ichans(0),iproc))) then
          write (*,*) 'Multi-channels not the same among identical processes'
          write (*,*) ichans(1:ichans(0))
          write (*,*) multi_chans(1:ichans(0),iproc)
          stop 1
       endif
       if (any(channel_permutations(1:next,1:ichans(0)).ne.&
            multi_chan_permutations(1:next,1:ichans(0),iproc))) then
          write (*,*) 'Phase-space permutations not the same among identical processes'
          stop 1
       endif
    endif
  end subroutine add_to_unique_process_list


  subroutine find_quarks(process,order,quarks)
    implicit none
    integer,dimension(next) :: process,order
    integer,dimension(0:6) :: quarks
    integer :: i
    quarks=0
    do i=1,next
       if (phys_model%is_quark(process(order(i))) .or. phys_model%is_antiquark(process(order(i)))) then
          quarks(0)=quarks(0)+1
          quarks(quarks(0))=process(order(i))
       endif
    enddo
    if (quarks(0).eq.4) then
       quarks(1:4)=[quarks(1),quarks(3),quarks(2),quarks(4)]
    elseif (quarks(0).eq.6) then
       quarks(1:6)=[quarks(1),quarks(3),quarks(5),quarks(2),quarks(4),quarks(6)]
    endif
  end subroutine find_quarks

  subroutine map_to_canonical_form(process,part,mapping)
    ! cross the two initial state particle PDGs, order according to
    ! the PDG value, (and reflip the two initial states again)
    implicit none
    integer,dimension(next) :: process,part,mapping
    part(1:next)=process(1:next)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
    call sort_with_mapping(next,part,mapping)
    ! cross the initial state
    part(1)=phys_model%get_antipart(part(1))
    part(2)=phys_model%get_antipart(part(2))
  end subroutine map_to_canonical_form

  subroutine sort_with_mapping(n,array,mapping)
    !
    ! EXAMPLE:
    !
    ! input:
    ! n=5
    ! array=[4, 1, 8, 2, 3]
    !
    ! output:
    ! array=[1, 2, 3, 4, 8]
    ! mapping=[2, 4, 5, 1, 3]
    !
    implicit none
    integer,intent(in) :: n
    integer,dimension(n),intent(inout) :: array
    integer,dimension(n),intent(out) :: mapping
    integer :: i, j, temp
    ! 'array_map' maps the PDG codes to the 'sort_particle' codes
    ! (see Utilities/process_list.py)
    integer,dimension(-24:25),parameter :: array_map=[83,0,0,0,0,0,0&
         &,0,95,88,93,86,91,84,0,0,0,0,12,11,10,9,8,7,0,1,2,3,4,5,6,0,0&
         &,0,0,85,90,87,92,89,94,0,0,0,0,13,80,81,82,96]
    ! Initialize mapping
    mapping = [(i,i=1,n)]
    ! Sort the array and mapping using a simple bubble sort
    do i=1,n-1
       do j=1,n-i
          if (array_map(array(j)) .gt. array_map(array(j+1))) then
             ! Swap array elements
             temp = array(j)
             array(j) = array(j+1)
             array(j+1) = temp
             ! Swap mapping
             temp = mapping(j)
             mapping(j) = mapping(j+1)
             mapping(j+1) = temp
          endif
       enddo
    enddo
  end subroutine sort_with_mapping


  subroutine move_colour_singlet_in_order(process,order)
    ! Move the colour singlet(s) to go *before* the final anti-quark in the order
    implicit none
    integer,dimension(next),intent(in) :: process
    integer,dimension(next),intent(inout) :: order
    integer :: i,iord,aq,iaq,ipart,is,ic
    integer,dimension(next) :: reordered

    ! A pure-gluon HEFT trace has no antiquark endpoint.  Canonicalise its
    ! colourless Higgs before the cyclic coloured trace, matching the current
    ! recursion while leaving the phase-space density order untouched.
    if (.not.any(abs(process).ge.1 .and. abs(process).le.6)) then
       is=1
       ic=count([(phys_model%is_singlet(process(order(i))),i=1,next)])+1
       do i=1,next
          iord=order(i)
          if (phys_model%is_singlet(process(iord))) then
             reordered(is)=iord
             is=is+1
          else
             reordered(ic)=iord
             ic=ic+1
          endif
       enddo
       order=reordered
       return
    endif
    ! find the final anti-quark
    do i=next,1,-1
       iord=order(i)
       ipart=process(iord)
       if (((iord.le.2 .and. phys_model%is_quark(ipart)).or.(iord.gt.2 .and. phys_model%is_antiquark(ipart)))) then
          aq=i
          iaq=iord
          exit
       endif
    enddo
    ! move the colour_singlet to after the anti-quark
    i=1
    do while (i.le.next)
       iord=order(i)
       ipart=process(iord)
       if(phys_model%is_singlet(ipart)) then
          if (i.gt.aq) then
             order(aq:i)=[order(i),order(aq)]
             aq=aq+1
             i=i+1
          else
             order(i:aq)=[order(i+1:aq),order(i)]
             aq=aq-1
          endif
       else
          i=i+1
       endif
    enddo
  end subroutine move_colour_singlet_in_order


end module read_process_file
