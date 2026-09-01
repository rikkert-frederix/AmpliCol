module read_process_file
  use handling_processes
  use amplitude_QCD_mod, only: max_process_external_particles => max_amplitude_external_particles
  use run_parameters, only: ignore_final_state_width_fix
  use random_number_interface, only: ran2
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: sf_nprocs
  integer,parameter :: max_process_file_unique_records=1000000
  integer,parameter :: max_process_file_groups=100000
  integer,parameter :: max_process_file_group_records=1000000
  integer(kind=8),parameter :: max_process_group_workspace_bytes=536870912_8
  integer(kind=8),parameter :: max_process_reduction_comparisons=100000000_8
  integer,dimension(:,:),allocatable :: unique_procs,processes,color_orders,multi_chans
  integer,dimension(:,:),allocatable :: phase_space_permutations
  integer,dimension(:,:,:),allocatable :: multi_chan_permutations
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map,iden_iproc
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: idenCOandMAPfactor
contains
  subroutine read_process_file_header(iunit,file_next,file_nunique,file_flavour_scheme,file_version)
    implicit none
    integer,intent(in) :: iunit
    integer,intent(out) :: file_next,file_nunique,file_flavour_scheme
    integer,intent(out),optional :: file_version
    character(len=1024) :: line
    integer :: ios,ipos,version

    read(iunit,'(a)',iostat=ios) line
    if (ios.ne.0) then
       write(*,*) 'ERROR: cannot read process-file header'
       stop 1
    endif
    file_flavour_scheme=0
    if (index(adjustl(line),'# AmpliCol process-list v2 nf=').eq.1) then
       ipos=index(line,'nf=')
       read(line(ipos+3:),*,iostat=ios) file_flavour_scheme
       if (ios.ne.0) then
          write(*,*) 'ERROR: malformed process-file flavour scheme header: ',trim(line)
          stop 1
       endif
       if (file_flavour_scheme.lt.1 .or. file_flavour_scheme.gt.5) then
          write(*,*) 'ERROR: invalid process-file flavour scheme header: ',trim(line)
          stop 1
       endif
       read(iunit,'(a)',iostat=ios) line
       if (ios.ne.0) then
          write(*,*) 'ERROR: cannot read numerical process-file header'
          stop 1
       endif
    endif

    version=1
    read(line,*,iostat=ios) file_next,file_nunique,version
    if (ios.ne.0) then
       version=1
       read(line,*,iostat=ios) file_next,file_nunique
    endif
    if (ios.ne.0) then
       write(*,*) 'ERROR: invalid process-file first record: ',trim(line)
       stop 1
    endif
    if (version.lt.1 .or. version.gt.2) then
       write(*,*) 'ERROR: unsupported process-file version: ',version
       stop 1
    endif
    ! The integration phase-space implementations require at least two
    ! final-state particles.  A 2->1 process has zero fixed-energy phase-space
    ! dimensions (one with PDFs) and needs a dedicated resonance mapping that
    ! this generator does not provide.
    if (file_next.lt.4 .or. file_next.gt.max_process_external_particles) then
       write(*,*) 'ERROR: unsupported process-file particle multiplicity:',file_next,&
            max_process_external_particles
       stop 1
    endif
    if (file_nunique.lt.1 .or. file_nunique.gt.max_process_file_unique_records) then
       write(*,*) 'ERROR: invalid or unsupported number of unique subprocesses:',file_nunique
       stop 1
    endif
    if (present(file_version)) file_version=version
  end subroutine read_process_file_header

  subroutine count_process_groups(filename,nfile_groups)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(out) :: nfile_groups
    integer :: file_next,file_nunique,file_flavour_scheme,i,ios,iunit
    character(len=256) :: io_message
    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios,&
         iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: cannot open process file: ',trim(filename),trim(io_message)
       stop 1
    endif
    nfile_groups=0
    call read_process_file_header(iunit,file_next,file_nunique,file_flavour_scheme)
    do i=1,file_nunique
       read(iunit,*,iostat=ios,iomsg=io_message)
       if (ios.ne.0) then
          write (*,*) 'ERROR: truncated unique-process block:',i,trim(io_message)
          stop 1
       endif
    enddo
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: missing separator after unique-process block: ',trim(io_message)
       stop 1
    endif
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: missing process-group preamble: ',trim(io_message)
       stop 1
    endif
    read(iunit,*,iostat=ios,iomsg=io_message) nfile_groups
    if (ios.ne.0) then
       write (*,*) 'ERROR: missing process-group count: ',trim(io_message)
       stop 1
    endif
    if (nfile_groups.lt.1 .or. nfile_groups.gt.max_process_file_groups) then
       write (*,*) 'ERROR: invalid or missing process-group count'
       stop 1
    endif
    close(iunit,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: cannot close process file: ',trim(filename),trim(io_message)
       stop 1
    endif
  end subroutine count_process_groups

  subroutine allocate_process_groups(ntotal_groups)
    implicit none
    integer,intent(in) :: ntotal_groups
    integer :: ios
    character(len=256) :: message
    if (ntotal_groups.lt.1 .or. ntotal_groups.gt.max_process_file_groups) then
       write (*,*) 'ERROR: invalid or unsupported process-group allocation:',ntotal_groups
       stop 1
    endif
    if (allocated(pgl)) deallocate(pgl)
    allocate(pgl(ntotal_groups),stat=ios,errmsg=message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: cannot allocate process groups:',trim(message)
       stop 1
    endif
    ngroups=ntotal_groups
  end subroutine allocate_process_groups

  subroutine clear_process_file_metadata()
    implicit none
    if (allocated(unique_procs)) deallocate(unique_procs)
    if (allocated(pgl_unique)) deallocate(pgl_unique)
    if (allocated(unique_map_value)) deallocate(unique_map_value)
    if (allocated(unique_map)) deallocate(unique_map)
  end subroutine clear_process_file_metadata

  subroutine apply_final_state_widths_from_process_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    character(len=65536) :: buffer
    integer :: iunit,ios,n_external,n_unique,version,iproc,igroup,ngroup_file,file_flavour_scheme,j
    integer :: group_index,nproc_in_group,max_channels,nchannels
    integer,dimension(:),allocatable :: phase_order,channels
    integer,dimension(:,:),allocatable :: physical_process
    character(len=256) :: allocation_message,io_message

    if (ignore_final_state_width_fix) return
    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios,&
         iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not open process file while applying final-state widths: ',&
            trim(filename),trim(io_message)
       stop 1
    endif
    call read_process_file_header(iunit,n_external,n_unique,file_flavour_scheme,version)
    if (n_external.lt.3 .or. n_unique.lt.1) then
       write (*,*) 'Malformed process-file header while applying final-state widths'
       stop 1
    endif
    do iproc=1,n_unique
       read(iunit,*,iostat=ios,iomsg=io_message)
       if (ios.ne.0) then
          write (*,*) 'Could not read unique subprocess while applying final-state widths',iproc
          stop 1
       endif
    enddo
    ! The first block contains canonical amplitude representatives, not the
    ! physical subprocess labelling.  Walk the phase-space-group rows instead;
    ! their first process field is the one later stored in iden_processes and
    ! has physical incoming/final-state positions.
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Missing separator after unique subprocesses while applying widths: ',&
            trim(io_message)
       stop 1
    endif
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Missing phase-space-group preamble while applying widths: ',&
            trim(io_message)
       stop 1
    endif
    read(iunit,*,iostat=ios,iomsg=io_message) ngroup_file
    if (ios.ne.0) then
       write (*,*) 'Could not read phase-space groups while applying final-state widths: ',&
            trim(io_message)
       stop 1
    endif
    if (ngroup_file.lt.1 .or. ngroup_file.gt.max_process_file_groups) then
       write (*,*) 'Could not read phase-space groups while applying final-state widths'
       stop 1
    endif
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Missing separator before phase-space groups while applying final-state widths'
       stop 1
    endif
    allocate(phase_order(n_external),physical_process(n_external,1),stat=ios,&
         errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate final-state-width process workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    do igroup=1,ngroup_file
       group_index=0
       nproc_in_group=0
       max_channels=0
       phase_order=0
       read(iunit,*,iostat=ios,iomsg=io_message) &
            group_index,nproc_in_group,max_channels,phase_order
       if (ios.ne.0) then
          write (*,*) 'Could not read phase-space group while applying final-state widths',&
               igroup,trim(io_message)
          stop 1
       endif
       if (group_index.ne.igroup .or. nproc_in_group.lt.1 .or.&
            nproc_in_group.gt.max_process_file_group_records .or. max_channels.lt.1 .or. &
            max_channels.gt.ngroup_file) then
          write (*,*) 'Malformed phase-space group while applying final-state widths',igroup
          stop 1
       endif
       do j=1,n_external
          if (count(phase_order.eq.j).ne.1) then
             write (*,*) 'Invalid phase-space order while applying final-state widths',igroup
             stop 1
          endif
       enddo
       do iproc=1,nproc_in_group
          buffer=''
          read(iunit,'(a)',iostat=ios,iomsg=io_message) buffer
          if (ios.ne.0) then
             write (*,*) 'Could not read subprocess row while applying final-state widths',&
                  igroup,iproc,trim(io_message)
             stop 1
          endif
          if (len_trim(buffer).eq.len(buffer)) then
             write (*,*) 'Could not read subprocess row while applying final-state widths',igroup,iproc
             stop 1
          endif
          nchannels=0
          read(buffer,*,iostat=ios,iomsg=io_message) nchannels
          if (ios.ne.0) then
             write (*,*) 'Malformed multichannel count while applying final-state widths',&
                  igroup,iproc,trim(io_message)
             stop 1
          endif
          if (nchannels.lt.1 .or. nchannels.gt.max_channels) then
             write (*,*) 'Invalid multichannel count while applying final-state widths',igroup,iproc
             stop 1
          endif
          allocate(channels(nchannels),stat=ios,errmsg=allocation_message)
          if (ios.ne.0) then
             write (*,*) 'Could not allocate final-state-width channel workspace: ',&
                  trim(allocation_message)
             stop 1
          endif
          read(buffer,*,iostat=ios,iomsg=io_message) nchannels,channels,physical_process(:,1)
          deallocate(channels)
          if (ios.ne.0) then
             write (*,*) 'Malformed subprocess process field while applying final-state widths',igroup,iproc
             stop 1
          endif
          call phys_model%apply_final_state_widths(n_external,1,physical_process)
       enddo
       do j=1,3
          read(iunit,*,iostat=ios,iomsg=io_message)
          if (ios.ne.0) then
             write (*,*) 'Missing process-group separator while applying final-state widths',&
                  igroup,j
             stop 1
          endif
       enddo
    enddo
    close(iunit,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not close process file after applying final-state widths: ',&
            trim(io_message)
       stop 1
    endif
    deallocate(phase_order)
    deallocate(physical_process)
  end subroutine apply_final_state_widths_from_process_file

  subroutine apply_final_state_widths_from_loaded_groups()
    implicit none
    integer :: igroup,iproc,iident,ios
    integer(kind=8) :: total_physical
    integer,dimension(:,:),allocatable :: physical_processes
    character(len=256) :: allocation_message

    if (ignore_final_state_width_fix) return
    do igroup=1,ngroups
       if (.not.allocated(pgl(igroup)%iden_processes)) cycle
       if (.not.allocated(pgl(igroup)%iden_iproc)) then
          write (*,*) 'Missing loaded identical-process metadata while applying widths',igroup
          stop 1
       endif
       if (pgl(igroup)%nproc.lt.1 .or. &
            pgl(igroup)%nproc.gt.size(pgl(igroup)%iden_iproc)) then
          write (*,*) 'Invalid loaded identical-process metadata while applying widths',igroup
          stop 1
       endif
       total_physical=0_8
       do iproc=1,pgl(igroup)%nproc
          if (pgl(igroup)%iden_iproc(iproc).lt.1) then
             write (*,*) 'Invalid loaded identical-process count while applying widths',&
                  igroup,iproc,pgl(igroup)%iden_iproc(iproc)
             stop 1
          endif
          if (pgl(igroup)%iden_iproc(iproc).gt.size(pgl(igroup)%iden_processes,2)) then
             write (*,*) 'Loaded identical-process count exceeds its workspace',&
                  igroup,iproc,pgl(igroup)%iden_iproc(iproc)
             stop 1
          endif
          total_physical=total_physical+int(pgl(igroup)%iden_iproc(iproc),kind=8)
          if (total_physical.gt.int(max_process_file_group_records,kind=8)) then
             write (*,*) 'Loaded process group is too large while applying widths',igroup,&
                  total_physical
             stop 1
          endif
       enddo
       allocate(physical_processes(pgl(igroup)%next,int(total_physical)),stat=ios,&
            errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'Could not allocate loaded width-update workspace: ',&
               trim(allocation_message)
          stop 1
       endif
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
  end subroutine apply_final_state_widths_from_loaded_groups

  subroutine read_processes_from_file(filename,igroup_offset,file_flavour_scheme)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(in),optional :: igroup_offset
    integer,intent(out),optional :: file_flavour_scheme
    integer :: iproc,igroup,igroup_global,icheck,nproc_in_group,max_channels,iflav,ndim,&
         nfile_groups,offset,nf,process_file_version,max_iden
    real(kind=8) :: idenCOfactor
    integer,dimension(:),allocatable :: process,order,ichans,phase_space_orders,phase_permutation
    integer,dimension(:,:),allocatable :: channel_permutations
    character(len=65536) :: buff
    integer :: i,j,ios,iunit
    integer(kind=8) :: nrecord8,nexternal8,nchannel8,workspace_bytes
    character(len=256) :: allocation_message,io_message

    offset=0
    if (present(igroup_offset)) offset=igroup_offset
    if (allocated(pgl_unique) .or. allocated(unique_procs) .or. allocated(unique_map) .or.&
         allocated(unique_map_value)) then
       call clear_process_file_metadata()
    endif
    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios,&
         iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: cannot open process file: ',trim(filename),trim(io_message)
       stop 1
    endif
    nfile_groups=0
    call read_process_file_header(iunit,next,nproc_unique,nf,process_file_version)
    if (present(file_flavour_scheme)) file_flavour_scheme=nf
    ndim=3*(next-2)-4
    if (include_pdf) ndim=ndim+2
    allocate(unique_procs(1:next,1:nproc_unique),stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate unique-process metadata: ',trim(allocation_message)
       stop 1
    endif
    do iproc=1,nproc_unique
       read(iunit,*,iostat=ios,iomsg=io_message) unique_procs(1:next,iproc)
       if (ios.ne.0) then
          write (*,*) 'ERROR: malformed unique subprocess:',iproc
          stop 1
       endif
    enddo
    allocate(pgl_unique,stat=ios,errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate unique-process group: ',trim(allocation_message)
       stop 1
    endif
    pgl_unique%next=next
    pgl_unique%ndim=ndim
    call check_unique_processes()
    read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.eq.0) read(iunit,*,iostat=ios,iomsg=io_message)
    if (ios.eq.0) read(iunit,*,iostat=ios,iomsg=io_message) nfile_groups
    if (ios.ne.0) then
       write (*,*) 'ERROR: missing process-group count: ',trim(io_message)
       stop 1
    endif
    if (nfile_groups.lt.1 .or. nfile_groups.gt.max_process_file_groups) then
       write (*,*) 'ERROR: invalid or missing process-group count'
       stop 1
    endif
    if (.not.allocated(pgl)) then
       call allocate_process_groups(nfile_groups)
       offset=0
    endif
    if (offset.lt.0 .or. offset.gt.size(pgl)) then
       write (*,*) 'ERROR: process-file group range is invalid',offset,nfile_groups,size(pgl)
       stop 1
    endif
    if (nfile_groups.gt.size(pgl)-offset) then
       write (*,*) 'ERROR: process-file group range is invalid',offset,nfile_groups,size(pgl)
       stop 1
    endif

    allocate(process(1:next),order(1:next),phase_permutation(1:next),stat=ios,&
         errmsg=allocation_message)
    if (ios.ne.0) then
       write (*,*) 'Could not allocate subprocess-row workspace: ',trim(allocation_message)
       stop 1
    endif

    read (iunit,*,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: missing separator before process groups'
       stop 1
    endif
    do igroup=1,nfile_groups
       igroup_global=offset+igroup
       nprocs=0
       sf_nprocs=0
       icheck=0
       nproc_in_group=0
       max_channels=0
       allocate(phase_space_orders(1:next),stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'Could not allocate phase-space-order workspace: ',&
               trim(allocation_message)
          stop 1
       endif
       read(iunit,*,iostat=ios,iomsg=io_message) &
            icheck,nproc_in_group,max_channels,phase_space_orders(1:next)
       if (ios.ne.0) then
          write (*,*) 'Could not read process-file group header',igroup,&
               trim(io_message)
          stop 1
       endif
       if (icheck.ne.igroup .or. nproc_in_group.lt.1 .or. &
            nproc_in_group.gt.max_process_file_group_records .or. max_channels.lt.1 .or. &
            max_channels.gt.nfile_groups) then
          write (*,*) 'ERROR in process-file group header',igroup,icheck,&
               nproc_in_group,max_channels
          stop 1
       endif
       do i=1,next
          if (count(phase_space_orders.eq.i).ne.1) then
             write (*,*) 'Invalid group phase-space order:',phase_space_orders
             stop 1
          endif
       enddo
       nrecord8=int(nproc_in_group,kind=8)
       nexternal8=int(next,kind=8)
       nchannel8=int(max_channels,kind=8)
       workspace_bytes=12_8*nexternal8*nrecord8+&
            4_8*nrecord8+&
            (4_8*nexternal8+8_8)*nrecord8*nrecord8+&
            4_8*(nchannel8+1_8)*nrecord8+&
            4_8*nexternal8*nchannel8*nrecord8+&
            4_8*(nchannel8+1_8)+4_8*nexternal8*nchannel8
       if (workspace_bytes.gt.max_process_group_workspace_bytes) then
          write (*,*) 'Process group exceeds supported parser workspace:',igroup,&
               workspace_bytes,max_process_group_workspace_bytes
          stop 1
       endif
       allocate(iden_iproc(nproc_in_group),processes(1:next,nproc_in_group),&
            color_orders(1:next,nproc_in_group),&
            phase_space_permutations(1:next,nproc_in_group),&
            iden_processes(1:next,nproc_in_group,nproc_in_group),&
            idenCOandMAPfactor(nproc_in_group,nproc_in_group),&
            multi_chans(0:max_channels,nproc_in_group),&
            multi_chan_permutations(1:next,1:max_channels,nproc_in_group),&
            ichans(0:max_channels),channel_permutations(1:next,1:max_channels),&
            stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'Could not allocate process-group parser workspace: ',&
               trim(allocation_message)
          stop 1
       endif
       phase_space_permutations=0
       multi_chan_permutations=0
       ichans=0
       channel_permutations=0
       do iflav=1,2
          do iproc=1,nproc_in_group
             buff=''
             read(iunit,'(a)',iostat=ios,iomsg=io_message) buff
             if (ios.ne.0) then
                write (*,*) 'Could not read a complete subprocess row: ',trim(io_message)
                stop 1
             endif
             if (len_trim(buff).eq.len(buff)) then
                write (*,*) 'Could not read a complete subprocess row'
                stop 1
             endif
             read(buff,*,iostat=ios) ichans(0)
             if (ios.ne.0) then
                write (*,*) 'Malformed number of multichannel partners in subprocess row'
                stop 1
             endif
             if (ichans(0).lt.1 .or. ichans(0).gt.max_channels) then
                write (*,*) 'Invalid number of multichannel partners in subprocess row'
                stop 1
             endif
             read(buff,*,iostat=ios) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),&
                  idenCOfactor,phase_permutation(1:next),&
                  channel_permutations(1:next,1:ichans(0))
             if (ios.ne.0) then
                if (process_file_version.ge.2) then
                   write (*,*) 'Malformed version-2 subprocess row; phase-space maps are required'
                   stop 1
                endif
                ! Backward compatibility with process files written before
                ! per-density external-leg permutations were introduced.
                read(buff,*,iostat=ios) ichans(0),ichans(1:ichans(0)),process(1:next),order(1:next),idenCOfactor
                if (ios.ne.0) then
                   write (*,*) 'Malformed legacy subprocess row'
                   stop 1
                endif
                phase_permutation=[(i,i=1,next)]
                do i=1,ichans(0)
                   channel_permutations(1:next,i)=[(j,j=1,next)]
                enddo
             endif
             if (.not.ieee_is_finite(idenCOfactor)) then
                write (*,*) 'Non-finite subprocess symmetry or mapping factor'
                stop 1
             endif
             do i=1,next
                if (count(order.eq.i).ne.1) then
                   write (*,*) 'Invalid subprocess colour order:',order
                   stop 1
                endif
             enddo
             if (any(ichans(1:ichans(0)).lt.1) .or. any(ichans(1:ichans(0)).gt.nfile_groups)) then
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
                backspace(iunit,iostat=ios,iomsg=io_message)
                if (ios.ne.0) then
                   write (*,*) 'Could not rewind subprocess row for flavour decomposition: ',&
                        igroup,iproc,trim(io_message)
                   stop 1
                endif
             enddo
          endif
       enddo
       pgl(igroup_global)%next=next
       pgl(igroup_global)%nproc=nprocs
       if (pgl(igroup_global)%nproc.lt.1) then
          write (*,*) 'Process-file group contains no usable subprocesses:',igroup
          stop 1
       endif
       pgl(igroup_global)%ndim=ndim
       pgl(igroup_global)%ndim_extra=0
       pgl(igroup_global)%multichan%max_channels=max_channels
       if (keep_processes_separate) then
          allocate(pgl(igroup_global)%amps(1:pgl(igroup_global)%nproc),&
               pgl(igroup_global)%nhel(1:pgl(igroup_global)%nproc),&
               pgl(igroup_global)%passed(1:pgl(igroup_global)%nproc),&
               stat=ios,errmsg=allocation_message)
       else
          allocate(pgl(igroup_global)%amps(1),pgl(igroup_global)%nhel(1),&
               pgl(igroup_global)%passed(1),stat=ios,errmsg=allocation_message)
       endif
       if (ios.ne.0) then
          write (*,*) 'Could not allocate process-group amplitude workspace: ',&
               trim(allocation_message)
          stop 1
       endif
       max_iden=maxval(iden_iproc(1:pgl(igroup_global)%nproc))
       allocate(pgl(igroup_global)%processes(1:next,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%color_orders(1:next,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%phase_space_permutations(1:next,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%phase_space_orders(1:next),&
            pgl(igroup_global)%idenCOandMAPfactor(1:max_iden,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%iden_iproc(1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%iden_processes(1:next,1:max_iden,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%val_procs(1:max_iden,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%multichan%channels(1:max_channels,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%multichan%channel_permutations(&
            1:next,1:max_channels,1:pgl(igroup_global)%nproc),&
            pgl(igroup_global)%multichan%number_of_channels(1:pgl(igroup_global)%nproc),&
            stat=ios,errmsg=allocation_message)
       if (ios.ne.0) then
          write (*,*) 'Could not allocate loaded process-group metadata: ',&
               trim(allocation_message)
          stop 1
       endif
       pgl(igroup_global)%processes=processes(1:next,1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%color_orders=color_orders(1:next,1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%phase_space_permutations=&
            phase_space_permutations(1:next,1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%phase_space_orders=phase_space_orders(1:next)
       pgl(igroup_global)%idenCOandMAPfactor=&
            idenCOandMAPfactor(1:max_iden,&
            1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%iden_iproc=iden_iproc(1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%iden_processes=&
            iden_processes(1:next,1:max_iden,&
            1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%val_procs=0d0
       pgl(igroup_global)%multichan%channels=0
       pgl(igroup_global)%multichan%channel_permutations=0
       do iproc=1,pgl(igroup_global)%nproc
          pgl(igroup_global)%multichan%channels(1:multi_chans(0,iproc),iproc)=&
               multi_chans(1:multi_chans(0,iproc),iproc)+offset
          pgl(igroup_global)%multichan%channel_permutations(:,1:multi_chans(0,iproc),iproc)=&
               multi_chan_permutations(:,1:multi_chans(0,iproc),iproc)
       enddo
       pgl(igroup_global)%multichan%number_of_channels=&
            multi_chans(0,1:pgl(igroup_global)%nproc)
       pgl(igroup_global)%passed=0
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
       do iproc=1,pgl(igroup_global)%nproc
          write(99,*) iproc,':',pgl(igroup_global)%processes(1:next,iproc),' ; ',&
               pgl(igroup_global)%color_orders(1:next,iproc),' ; ',pgl(igroup_global)%iden_iproc(iproc),' ; ',&
               pgl(igroup_global)%multichan%channels(1:pgl(igroup_global)%multichan%number_of_channels(iproc),iproc)
       enddo
       write (99,*) '****************************************************'
       do j=1,3
          read(iunit,*,iostat=ios,iomsg=io_message)
          if (ios.ne.0) then
             write (*,*) 'ERROR: missing separator after process group',igroup,j
             stop 1
          endif
       enddo
    enddo
    close(iunit,iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'ERROR: cannot close process file: ',trim(filename),trim(io_message)
       stop 1
    endif
  end subroutine read_processes_from_file



  subroutine check_unique_processes()
    use phase_space_gen23_mod
    use cuts
    implicit none
    integer :: i,iproc,ih,evaluation_status,nattempt,isample,allocation_status,&
         coordinate_count
    integer,parameter :: nevent=10,max_probe_attempts=100000
    integer(kind=8) :: workspace_bytes
    real(kind=8),dimension(:,:),allocatable :: amp2
    real(kind=8),dimension(:),allocatable :: mass,width
    type(psv) :: ps
    character(len=256) :: allocation_message

    if (.not.allocated(pgl_unique) .or. .not.allocated(unique_procs) .or. &
         next.lt.4 .or. next.gt.max_process_external_particles .or. &
         nproc_unique.lt.1 .or. nproc_unique.gt.max_process_file_unique_records) then
       write (*,*) 'Invalid unique-process state before numerical reduction:',&
            next,nproc_unique
       stop 1
    endif
    if (size(unique_procs,1).ne.next .or. &
         size(unique_procs,2).ne.nproc_unique) then
       write (*,*) 'Incompatible unique-process array before numerical reduction'
       stop 1
    endif
    workspace_bytes=12_8*int(next,kind=8)*int(nproc_unique,kind=8)+&
         8_8*int(nevent,kind=8)*int(nproc_unique,kind=8)+&
         12_8*int(nproc_unique,kind=8)+16_8*int(next,kind=8)
    if (workspace_bytes.gt.max_process_group_workspace_bytes) then
       write (*,*) 'Unique-process reduction exceeds the supported workspace:',&
            workspace_bytes,max_process_group_workspace_bytes
       stop 1
    endif
    allocate(phase_space_gen23 :: pgl_unique%phase_space,&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process phase space: ',&
            trim(allocation_message)
       stop 1
    endif
    allocate(pgl_unique%processes(next,nproc_unique),&
         pgl_unique%color_orders(next,nproc_unique),&
         pgl_unique%phase_space_orders(next),pgl_unique%amps(1),&
         mass(next),width(next),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process reduction metadata: ',&
            trim(allocation_message)
       stop 1
    endif
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes(1:next,1:nproc_unique)=unique_procs(1:next,1:nproc_unique)
    do i=1,pgl_unique%next
       mass(i)=phys_model%get_mass(pgl_unique%processes(i,1))
       width(i)=phys_model%get_width(pgl_unique%processes(i,1))
       ! For this unique_prcess checks, use only massless particles
       if (mass(i).ne.0d0) mass(i)=0d0
       if (width(i).ne.0d0) width(i)=0d0
    enddo
    call setup_spin(pgl_unique)
    allocate(pgl_unique%hel(pgl_unique%next),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process helicity workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    pgl_unique%hel=pgl_unique%spin(1,1:pgl_unique%next)
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
    coordinate_count=pgl_unique%ndim+pgl_unique%phase_space%ndim_extra
    if (coordinate_count.lt.pgl_unique%ndim .or. coordinate_count.gt.1000000) then
       write (*,*) 'Invalid unique-process phase-space dimension:',&
            pgl_unique%ndim,pgl_unique%phase_space%ndim_extra
       stop 1
    endif
    allocate(ps%x(1:coordinate_count),ps%p(0:3,1:pgl_unique%next),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process point workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    
    call pgl_unique%amps(1)%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
         pgl_unique%spin,pgl_unique%color_orders,phys_model)
        
    allocate(amp2(nevent,pgl_unique%nproc),pgl_unique%passed(1),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process sample workspace: ',&
            trim(allocation_message)
       stop 1
    endif
    amp2=0d0

    pgl_unique%passed=0
    nattempt=0
    do while (pgl_unique%passed(1).lt.nevent)
       nattempt=nattempt+1
       if (nattempt.gt.max_probe_attempts) then
          write (*,*) 'Could not find enough finite process-reduction probe points',&
               pgl_unique%passed(1),nevent
          stop 1
       endif
       do i=1,pgl_unique%ndim+pgl_unique%phase_space%ndim_extra
          ps%x(i)=ran2()
       enddo
       if (.not.all(ieee_is_finite(ps%x))) then
          write (*,*) 'Random-number generator returned a non-finite process-reduction coordinate'
          stop 1
       endif
       if (any(ps%x.lt.0d0) .or. any(ps%x.ge.1d0)) then
          write (*,*) 'Random-number generator returned an invalid process-reduction coordinate'
          stop 1
       endif
       call pgl_unique%phase_space%generate_momenta(ps)
       if (.not.ieee_is_finite(ps%jac)) cycle
       if (ps%jac.le.0d0) cycle
       if (.not.all(ieee_is_finite(ps%p))) cycle
       call pgl_unique%amps(1)%evaluate(next,ps%p,pgl_unique%hel,&
            read_proc_from_file,phys_model,evaluation_status)
       if (evaluation_status.ne.0) cycle
       isample=pgl_unique%passed(1)+1
       iproc=0
       amp2(isample,:)=0d0
       if (use_real_gluons .and. all(pgl_unique%amps(1)%n_qqbar(1:pgl_unique%amps(1)%nprocs).eq.0)) then
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(isample,iproc)=amp2(isample,iproc)+&
                  pgl_unique%amps(1)%amps_r(ih)*pgl_unique%amps(1)%amps_r(ih)
          enddo
       else
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(isample,iproc)=amp2(isample,iproc)+&
                  dble(pgl_unique%amps(1)%amps(ih)*dconjg(pgl_unique%amps(1)%amps(ih)))
          enddo
       endif
       if (.not.all(ieee_is_finite(amp2(isample,:)))) cycle
       if (any(amp2(isample,:).lt.0d0)) cycle
       pgl_unique%passed(1)=isample
       call find_same_flavour(pgl_unique,nevent,amp2(1,:))
    enddo
    allocate(unique_map(1:pgl_unique%nproc),&
         unique_map_value(1:pgl_unique%nproc),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process reduction map: ',&
            trim(allocation_message)
       stop 1
    endif
    call find_unique(pgl_unique,nevent,amp2,unique_map,unique_map_value)

    do iproc=1,pgl_unique%nproc
       write (99,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:pgl_unique%next,iproc),&
            ':',pgl_unique%color_orders(1:pgl_unique%next,iproc)
    enddo
    deallocate(pgl_unique%spin)
    deallocate(pgl_unique%hel)
    call pgl_unique%phase_space%cleanup()
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
    integer :: i,j,k,bucket,bucket_count,allocation_status
    integer,dimension(:),allocatable :: bucket_head,next_candidate
    integer(kind=8) :: hash_value,quantized,comparison_count,workspace_bytes
    real(kind=8) :: ave,scale_i,scale_j
    real(kind=8),parameter :: reduction_tolerance=1d-10
    integer(kind=8),parameter :: hash_modulus=2147483647_8,&
         hash_multiplier=1000003_8
    character(len=256) :: allocation_message
    if (nevent.lt.2 .or. pgl%nproc.lt.1 .or. &
         pgl%nproc.gt.max_process_file_unique_records .or. &
         pgl%next.lt.4 .or. pgl%next.gt.max_process_external_particles .or. &
         .not.allocated(pgl%processes)) then
       write (*,*) 'Invalid process-reduction dimensions:',&
            nevent,pgl%nproc,pgl%next
       stop 1
    endif
    if (size(pgl%processes,1).ne.pgl%next .or. &
         size(pgl%processes,2).ne.pgl%nproc) then
       write (*,*) 'Incompatible process metadata during numerical reduction'
       stop 1
    endif
    if (.not.all(ieee_is_finite(amp2))) then
       write (*,*) 'Non-finite process-reduction matrix-element samples'
       stop 1
    endif
    if (any(amp2.lt.0d0)) then
       write (*,*) 'Invalid process-reduction matrix-element samples'
       stop 1
    endif
    bucket_count=2*pgl%nproc+1
    workspace_bytes=4_8*int(bucket_count+pgl%nproc,kind=8)
    if (workspace_bytes.gt.max_process_group_workspace_bytes) then
       write (*,*) 'Process-reduction lookup exceeds the supported workspace:',&
            workspace_bytes,max_process_group_workspace_bytes
       stop 1
    endif
    allocate(bucket_head(bucket_count),next_candidate(pgl%nproc),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate process-reduction lookup: ',&
            trim(allocation_message)
       stop 1
    endif
    bucket_head=0
    next_candidate=0
    unique_map=-1
    unique_map_value=1d0
    comparison_count=0_8
    do i=1,pgl%nproc
       if (all(amp2(1:nevent,i).eq.0d0)) then
          unique_map(i)=0
          unique_map_value(i)=0d0
          cycle
       endif
       scale_i=maxval(amp2(1:nevent,i))
       if (scale_i.le.0d0) then
          write (*,*) 'Invalid nonzero process-reduction scale:',i,scale_i
          stop 1
       endif
       hash_value=0_8
       do k=1,nevent
          ! Divide first: every sample is in [0,scale_i], whereas multiplying
          ! an otherwise valid near-HUGE sample by 1e9 before normalization
          ! can overflow under trapping arithmetic.
          quantized=nint(1d9*(amp2(k,i)/scale_i),kind=8)
          hash_value=modulo(hash_value*hash_multiplier+&
               modulo(quantized,hash_modulus),hash_modulus)
       enddo
       bucket=1+int(modulo(hash_value,int(bucket_count,kind=8)))
       j=bucket_head(bucket)
       do while (j.gt.0)
          if (comparison_count.ge.max_process_reduction_comparisons) then
             write (*,*) 'Process reduction exceeds the supported comparison budget:',&
                  max_process_reduction_comparisons
             stop 1
          endif
          comparison_count=comparison_count+1_8
          ! A numerical relation found on the deliberately massless probe
          ! point is not a valid runtime reduction if it moves a massive
          ! species to another external leg.
          do k=1,pgl%next
             if (phys_model%get_mass(pgl%processes(k,i)).ne.&
                  phys_model%get_mass(pgl%processes(k,j))) exit
          enddo
          if (k.le.pgl%next) then
             j=next_candidate(j)
             cycle
          endif
          scale_j=maxval(amp2(1:nevent,j))
          if (scale_j.le.0d0) then
             j=next_candidate(j)
             cycle
          endif
          if (scale_j.lt.1d0) then
             if (scale_i.gt.huge(1d0)*scale_j) then
                j=next_candidate(j)
                cycle
             endif
          endif
          ave=scale_i/scale_j
          if (.not.ieee_is_finite(ave)) then
             j=next_candidate(j)
             cycle
          endif
          if (maxval(abs(amp2(1:nevent,i)/scale_i-&
               amp2(1:nevent,j)/scale_j)).le.reduction_tolerance) then
             unique_map_value(i)=ave
             unique_map(i)=j
             exit
          endif
          j=next_candidate(j)
       enddo
       if (unique_map(i).eq.-1) then
          unique_map(i)=-1
          unique_map_value(i)=1d0
          next_candidate(i)=bucket_head(bucket)
          bucket_head(bucket)=i
       endif
    enddo
    deallocate(bucket_head,next_candidate)
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
       if (idenMAPfactor.ne.0d0 .and. idenCOfactor.ne.0d0) then
          if (abs(idenMAPfactor).gt.huge(1d0)/abs(idenCOfactor)) then
             write (*,*) 'Subprocess symmetry or mapping factor overflows',&
                  idenMAPfactor,idenCOfactor
             stop 1
          endif
       endif
       idenCOMAPfactor=idenMAPfactor*idenCOfactor
    else
       idenCOMAPfactor=idenCOfactor
    endif
    if (.not.ieee_is_finite(idenCOMAPfactor)) then
       write (*,*) 'Non-finite combined subprocess symmetry or mapping factor'
       stop 1
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
       if (nprocs.ge.size(processes,2)) then
          write (*,*) 'Process-group workspace exhausted while adding subprocess'
          stop 1
       endif
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
       if (nprocs.ge.size(processes,2)) then
          write (*,*) 'Process-group workspace exhausted while adding unique subprocess'
          stop 1
       endif
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
       if (iden_iproc(iproc).ge.size(iden_processes,2)) then
          write (*,*) 'Identical-subprocess workspace exhausted',iproc,iden_iproc(iproc)
          stop 1
       endif
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
          if (quarks(0).ge.ubound(quarks,1)) then
             write (*,*) 'More than three quark lines in process-file subprocess:',process
             stop 1
          endif
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
    if (any(array.lt.lbound(array_map,1)) .or. any(array.gt.ubound(array_map,1))) then
       write (*,*) 'Unsupported particle identifier while sorting subprocess:',array
       stop 1
    endif
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
    ! Move colour singlets immediately before the final coloured leg.  The
    ! latter closes the recursive current, and therefore must remain last.
    implicit none
    integer,dimension(next),intent(in) :: process
    integer,dimension(next),intent(inout) :: order
    integer,dimension(next) :: reordered
    integer :: i,iord,ncoloured,nentry,closing_leg

    ncoloured=0
    do i=1,next
       iord=order(i)
       if (.not.phys_model%is_singlet(process(iord))) then
          ncoloured=ncoloured+1
          reordered(ncoloured)=iord
       endif
    enddo
    if (ncoloured.eq.0 .or. ncoloured.eq.next) return
    closing_leg=reordered(ncoloured)
    nentry=ncoloured-1
    do i=1,next
       iord=order(i)
       if (phys_model%is_singlet(process(iord))) then
          nentry=nentry+1
          reordered(nentry)=iord
       endif
    enddo
    nentry=nentry+1
    reordered(nentry)=closing_leg
    if (nentry.ne.next) then
       write (*,*) 'Internal colour-order reconstruction failure:',nentry,next
       stop 1
    endif
    order=reordered
  end subroutine move_colour_singlet_in_order


end module read_process_file
