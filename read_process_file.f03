module read_process_file
  use handling_processes
  use run_parameters, only: ignore_final_state_width_fix
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: sf_nprocs
  integer,dimension(:,:),allocatable :: unique_procs,processes,color_orders,multi_chans
  integer,dimension(:,:),allocatable :: phase_space_permutations
  integer,dimension(:,:,:),allocatable :: multi_chan_permutations
  type(phase_space_order_group),allocatable :: pgl_unique
  real(kind=8),dimension(:),allocatable :: unique_map_value
  integer,dimension(:),allocatable :: unique_map,iden_iproc
  integer,dimension(:,:,:),allocatable :: iden_processes
  real(kind=8),dimension(:,:),allocatable :: idenCOandMAPfactor
contains
  subroutine apply_final_state_widths_from_process_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    character(len=65536) :: buffer
    integer :: iunit,ios,n_external,n_unique,version,iproc,igroup,ngroup_file
    integer :: group_index,nproc_in_group,max_channels,nchannels
    integer,dimension(:),allocatable :: phase_order,channels
    integer,dimension(:,:),allocatable :: physical_process

    if (ignore_final_state_width_fix) return
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
    version=1
    read(buffer,*,iostat=ios) n_external,n_unique,version
    if (ios.ne.0) then
       read(buffer,*,iostat=ios) n_external,n_unique
    endif
    if (ios.ne.0 .or. n_external.lt.3 .or. n_unique.lt.1) then
       write (*,*) 'Malformed process-file header while applying final-state widths'
       stop 1
    endif
    do iproc=1,n_unique
       read(iunit,*,iostat=ios)
       if (ios.ne.0) then
          write (*,*) 'Could not read unique subprocess while applying final-state widths',iproc
          stop 1
       endif
    enddo
    ! The first block contains canonical amplitude representatives, not the
    ! physical subprocess labelling.  Walk the phase-space-group rows instead;
    ! their first process field is the one later stored in iden_processes and
    ! has physical incoming/final-state positions.
    read(iunit,*,iostat=ios)
    read(iunit,*,iostat=ios)
    read(iunit,*,iostat=ios) ngroup_file
    if (ios.ne.0 .or. ngroup_file.lt.1) then
       write (*,*) 'Could not read phase-space groups while applying final-state widths'
       stop 1
    endif
    read(iunit,*,iostat=ios)
    allocate(phase_order(n_external))
    allocate(physical_process(n_external,1))
    do igroup=1,ngroup_file
       read(iunit,*,iostat=ios) group_index,nproc_in_group,max_channels,phase_order
       if (ios.ne.0 .or. group_index.ne.igroup .or. nproc_in_group.lt.1 .or.&
            max_channels.lt.1) then
          write (*,*) 'Malformed phase-space group while applying final-state widths',igroup
          stop 1
       endif
       do iproc=1,nproc_in_group
          read(iunit,'(a)',iostat=ios) buffer
          if (ios.ne.0 .or. len_trim(buffer).eq.len(buffer)) then
             write (*,*) 'Could not read subprocess row while applying final-state widths',igroup,iproc
             stop 1
          endif
          read(buffer,*,iostat=ios) nchannels
          if (ios.ne.0 .or. nchannels.lt.1 .or. nchannels.gt.max_channels) then
             write (*,*) 'Invalid multichannel count while applying final-state widths',igroup,iproc
             stop 1
          endif
          allocate(channels(nchannels))
          read(buffer,*,iostat=ios) nchannels,channels,physical_process(:,1)
          deallocate(channels)
          if (ios.ne.0) then
             write (*,*) 'Malformed subprocess process field while applying final-state widths',igroup,iproc
             stop 1
          endif
          call phys_model%apply_final_state_widths(n_external,1,physical_process)
       enddo
       read(iunit,*,iostat=ios)
       read(iunit,*,iostat=ios)
       read(iunit,*,iostat=ios)
       if (ios.ne.0 .and. igroup.ne.ngroup_file) then
          write (*,*) 'Unexpected end of process file while applying final-state widths'
          stop 1
       endif
    enddo
    close(iunit)
    deallocate(phase_order)
    deallocate(physical_process)
  end subroutine apply_final_state_widths_from_process_file

  subroutine apply_final_state_widths_from_loaded_groups()
    implicit none
    integer :: igroup,iproc,iident
    integer,dimension(:,:),allocatable :: physical_processes

    if (ignore_final_state_width_fix) return
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
  end subroutine apply_final_state_widths_from_loaded_groups

  subroutine read_processes_from_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    integer :: iproc,igroup,icheck,nproc_in_group,max_channels,iflav,ndim,process_file_version
    real(kind=8) :: idenCOfactor
    integer,dimension(:),allocatable :: process,order,ichans,phase_space_orders,phase_permutation
    integer,dimension(:,:),allocatable :: channel_permutations
    character(len=65536) :: buff
    integer :: i,j,ios,iunit
    call read_unique_process_catalogue(filename,iunit,process_file_version)
    ndim=3*(next-2)-4
    if (include_pdf) ndim=ndim+2
    allocate(pgl_unique)
    pgl_unique%next=next
    pgl_unique%ndim=ndim
    call check_unique_processes()
    read(iunit,*)
    read(iunit,*)
    read (iunit,*) ngroups
    allocate(pgl(ngroups))

    allocate(process(1:next))
    allocate(order(1:next))
    allocate(phase_permutation(1:next))

    read (iunit,*)
    do igroup=1,ngroups
       nprocs=0
       sf_nprocs=0
       allocate(phase_space_orders(1:next))
       read(iunit,*) icheck,nproc_in_group,max_channels,phase_space_orders(1:next)
       if (icheck.ne.igroup) then
          write (*,*) 'ERROR in processes file',icheck,igroup
          stop 1
       endif
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
             read(iunit,'(a)',iostat=ios) buff
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
             if (.not.decompose_same_flavour_amplitudes) exit
             do iproc=1,nproc_in_group
                backspace(iunit)
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
       read(iunit,*)
       read(iunit,*)
       read(iunit,*)
    enddo
    close(iunit)
  end subroutine read_processes_from_file

  subroutine read_unique_process_catalogue(filename,iunit,version)
    implicit none
    character(len=*),intent(in) :: filename
    integer,intent(out) :: iunit,version
    character(len=65536) :: line
    integer :: ios,iproc,jproc,nquarks,nantiquarks,i,nfields
    integer,dimension(3) :: header

    open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios.ne.0) call process_catalogue_error('could not open '//trim(filename))
    read(iunit,'(a)',iostat=ios) line
    if (ios.ne.0) call process_catalogue_error('could not read the header')
    if (len_trim(line).eq.len(line)) &
         call process_catalogue_error('header exceeds the supported line length')
    call parse_integer_record(line,header,nfields,ios)
    if (ios.ne.0 .or. (nfields.ne.2 .and. nfields.ne.3)) &
         call process_catalogue_error('malformed header')
    next=header(1)
    nproc_unique=header(2)
    version=1
    if (nfields.eq.3) version=header(3)
    if (next.lt.4 .or. nproc_unique.lt.1) &
         call process_catalogue_error('invalid process or record count')
    if (version.lt.1 .or. version.gt.2) &
         call process_catalogue_error('unsupported version',version)

    allocate(unique_procs(next,nproc_unique))
    do iproc=1,nproc_unique
       read(iunit,'(a)',iostat=ios) line
       if (ios.ne.0) call process_catalogue_error(&
            'could not read unique process',iproc)
       if (len_trim(line).eq.len(line)) call process_catalogue_error(&
            'unique process exceeds the supported line length',iproc)
       call parse_integer_record(line,unique_procs(:,iproc),nfields,ios)
       if (ios.ne.0 .or. nfields.ne.next) call process_catalogue_error(&
            'malformed unique process',iproc)
       nquarks=count([(phys_model%is_quark(unique_procs(i,iproc)),i=1,next)])
       nantiquarks=count([(phys_model%is_antiquark(unique_procs(i,iproc)),i=1,next)])
       if (nquarks.ne.nantiquarks .or. nquarks.gt.3) &
            call process_catalogue_error('unbalanced or unsupported quark-line count',iproc)
       if (nquarks.gt.0) then
          if (.not.all([(phys_model%is_quark(unique_procs(i,iproc)),i=1,nquarks)])) &
               call process_catalogue_error('quark endpoints are not first',iproc)
          if (.not.all([(phys_model%is_antiquark(unique_procs(i,iproc)),&
               i=nquarks+1,2*nquarks)])) &
               call process_catalogue_error('antiquark endpoints do not follow quarks',iproc)
       endif
       do jproc=1,iproc-1
          if (all(unique_procs(:,iproc).eq.unique_procs(:,jproc))) &
               call process_catalogue_error('duplicate unique process',iproc)
       enddo
    enddo
  end subroutine read_unique_process_catalogue

  subroutine parse_integer_record(line,values,nvalues,ios)
    implicit none
    character(len=*),intent(in) :: line
    integer,dimension(:),intent(out) :: values
    integer,intent(out) :: nvalues,ios
    integer :: first,last,position
    values=0
    nvalues=0
    ios=0
    last=len_trim(line)
    position=1
    do
       do while (position.le.last .and. is_record_space(line(position:position)))
          position=position+1
       enddo
       if (position.gt.last) return
       first=position
       do while (position.le.last .and. .not.is_record_space(line(position:position)))
          position=position+1
       enddo
       if (nvalues.ge.size(values)) then
          ios=1
          return
       endif
       nvalues=nvalues+1
       read(line(first:position-1),*,iostat=ios) values(nvalues)
       if (ios.ne.0) return
    enddo
  contains
    logical function is_record_space(token_character)
      implicit none
      character(len=1),intent(in) :: token_character
      is_record_space=token_character.eq.' ' .or. &
           iachar(token_character).eq.9 .or. iachar(token_character).eq.13
    end function is_record_space
  end subroutine parse_integer_record

  subroutine process_catalogue_error(message,record)
    implicit none
    character(len=*),intent(in) :: message
    integer,intent(in),optional :: record
    if (present(record)) then
       write (*,'(a,1x,a,1x,i0)') 'Invalid unique-process catalogue:',trim(message),record
    else
       write (*,'(a,1x,a)') 'Invalid unique-process catalogue:',trim(message)
    endif
    stop 1
  end subroutine process_catalogue_error



  subroutine check_unique_processes()
    use phase_space_gen23_mod
    use cuts
    implicit none
    integer :: i,iproc,ih
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    real(kind=8),dimension(:),allocatable :: mass
    real(kind=8),external :: ran2
    type(psv) :: ps
    allocate(phase_space_gen23 :: pgl_unique%phase_space)
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(pgl_unique%amps(1))
    allocate(mass(next))
    pgl_unique%nproc=nproc_unique
    pgl_unique%processes(1:next,1:nproc_unique)=unique_procs(1:next,1:nproc_unique)
    ! Use a common massless phase space for all unique-process probes.  The
    ! physical mass layout is checked separately before a map is accepted.
    mass=0d0
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
       call find_same_flavour(pgl_unique,nevent,&
            amp2(pgl_unique%passed(1),:))
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
    deallocate(mass)
  end subroutine check_unique_processes

  subroutine find_unique(pgl,nevent,amp2,unique_map,unique_map_value)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: nevent
    real(kind=8),dimension(nevent,pgl%nproc),intent(in) :: amp2
    real(kind=8),dimension(pgl%nproc),intent(out) :: unique_map_value
    integer,dimension(pgl%nproc),intent(out) :: unique_map
    integer :: i,j,representative
    real(kind=8),dimension(nevent) :: difference
    real(kind=8) :: denominator,factor,residual,point_residual,best_residual,scale,floor
    real(kind=8),parameter :: relation_tolerance=1d-8
    unique_map=-1
    unique_map_value=1d0
    do i=1,pgl%nproc
       if (maxval(abs(amp2(:,i))).eq.0d0) then
          unique_map(i)=0
          unique_map_value(i)=0d0
          cycle
       endif
       best_residual=huge(1d0)
       do j=1,i-1
          if (maxval(abs(amp2(:,j))).eq.0d0) cycle
          if (.not.compatible_unique_probes(pgl,i,j)) cycle
          denominator=dot_product(amp2(:,j),amp2(:,j))
          if (denominator.le.0d0 .or. .not.ieee_is_finite(denominator)) cycle
          ! Least-squares proportionality avoids divisions by individual
          ! phase-space values (and hence the NaNs produced by the old ratio
          ! test whenever a representative vanished at one sample).
          factor=dot_product(amp2(:,i),amp2(:,j))/denominator
          if (factor.le.0d0 .or. .not.ieee_is_finite(factor)) cycle
          difference=amp2(:,i)-factor*amp2(:,j)
          scale=sqrt(dot_product(amp2(:,i),amp2(:,i)))+&
               abs(factor)*sqrt(denominator)
          if (scale.le.0d0 .or. .not.ieee_is_finite(scale)) cycle
          residual=sqrt(dot_product(difference,difference))/scale
          floor=epsilon(1d0)*max(maxval(abs(amp2(:,i))),&
               abs(factor)*maxval(abs(amp2(:,j))))
          point_residual=maxval(abs(difference)/(&
               abs(amp2(:,i))+abs(factor*amp2(:,j))+floor))
          residual=max(residual,point_residual)
          if (.not.ieee_is_finite(residual)) cycle
          if (residual.lt.relation_tolerance .and. residual.lt.best_residual) then
             best_residual=residual
             unique_map_value(i)=factor
             unique_map(i)=j
          endif
       enddo
       ! Store a direct map to the root representative.  Chained maps made
       ! the reconstructed flavour labels depend on catalogue ordering and
       ! lost one or more proportionality factors.
       if (unique_map(i).gt.0) then
          representative=unique_map(i)
          do while (unique_map(representative).gt.0)
             unique_map_value(i)=unique_map_value(i)*&
                  unique_map_value(representative)
             representative=unique_map(representative)
          enddo
          unique_map(i)=representative
       endif
    enddo
  end subroutine find_unique

  logical function compatible_unique_probes(pgl,left,right)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: left,right
    integer :: i,nlines
    compatible_unique_probes=.false.
    if (pgl%amps(1)%n_qqbar(left).ne.pgl%amps(1)%n_qqbar(right)) return
    nlines=pgl%amps(1)%n_qqbar(left)
    ! The unique-process block varies quark endpoint labels deliberately, but
    ! all other external particles must agree exactly.  This prevents an
    ! accidental numerical relation from changing the requested EW state.
    if (2*nlines.lt.pgl%next) then
       if (any(unique_procs(2*nlines+1:pgl%next,left).ne.&
            unique_procs(2*nlines+1:pgl%next,right))) return
    endif
    do i=1,pgl%next
       if (phys_model%get_mass(pgl%processes(i,left)).ne.&
            phys_model%get_mass(pgl%processes(i,right))) return
       if (phys_model%get_spin(pgl%processes(i,left)).ne.&
            phys_model%get_spin(pgl%processes(i,right))) return
    enddo
    compatible_unique_probes=.true.
  end function compatible_unique_probes

  subroutine add_to_process_list(process,order,phase_permutation,channel_permutations,&
       idenCOfactor,max_channels,ichans,skip_same_flavour)
    implicit none
    integer :: max_channels
    integer,dimension(0:max_channels) :: ichans
    integer,dimension(next) :: process,order,process_unique,phase_permutation
    integer,dimension(next,max_channels) :: channel_permutations
    real(kind=8) :: idenCOfactor,idenCOMAPfactor,idenMAPfactor
    logical :: skip_same_flavour,is_same_flavour
    call get_unique_process(process,order,process_unique,idenMAPfactor,is_same_flavour)
    if (decompose_same_flavour_amplitudes) then
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


  subroutine get_unique_process(process,order,process_unique,idenMAPfactor,is_same_flavour)
    implicit none
    integer,dimension(next),intent(in) :: process,order
    integer,dimension(next),intent(out) :: process_unique
    integer,dimension(next) :: canonical_process,mapping
    integer,dimension(3) :: quark_pdgs,antiquark_pdgs,quark_legs,antiquark_legs
    integer :: i,iproc,nlines,outgoing_pdg
    real(kind=8),intent(out) :: idenMAPfactor
    logical,intent(out) :: is_same_flavour
    call collect_quark_endpoints(process,order,nlines,quark_pdgs,&
         antiquark_pdgs,quark_legs,antiquark_legs)

    ! Build exactly the layout serialized by process_list.py: all outgoing,
    ! quark endpoints first in colour-line order, then antiquark endpoints,
    ! followed by the canonically sorted non-fermion legs.  This replaces the
    ! separate, hand-written index shuffles for two and three quark lines.
    canonical_process=process
    canonical_process(1)=phys_model%get_antipart(canonical_process(1))
    canonical_process(2)=phys_model%get_antipart(canonical_process(2))
    call sort_with_mapping(next,canonical_process,mapping)
    if (nlines.gt.0) then
       canonical_process(1:nlines)=quark_pdgs(1:nlines)
       canonical_process(nlines+1:2*nlines)=antiquark_pdgs(1:nlines)
    endif
    do iproc=1,pgl_unique%nproc
       if (pgl_unique%amps(1)%n_qqbar(iproc).ne.nlines) cycle
       if (all(canonical_process.eq.unique_procs(:,iproc))) exit
    enddo
    if (iproc.gt.pgl_unique%nproc) then
       write (*,*) 'Process not found in unique-process catalogue'
       write (*,*) 'Physical process:',process
       write (*,*) 'Colour order:',order
       write (*,*) 'Canonical process:',canonical_process
       stop 1
    endif
    process_unique=process
    idenMAPfactor=unique_map_value(iproc)
    if (unique_map(iproc).gt.0) then
       do i=1,nlines
          outgoing_pdg=unique_procs(i,unique_map(iproc))
          if (quark_legs(i).le.2) outgoing_pdg=phys_model%get_antipart(outgoing_pdg)
          process_unique(quark_legs(i))=outgoing_pdg
          outgoing_pdg=unique_procs(nlines+i,unique_map(iproc))
          if (antiquark_legs(i).le.2) outgoing_pdg=phys_model%get_antipart(outgoing_pdg)
          process_unique(antiquark_legs(i))=outgoing_pdg
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
    is_same_flavour=pgl_unique%amps(1)%same_flav(iproc)
  end subroutine get_unique_process

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


  subroutine collect_quark_endpoints(process,order,nlines,quark_pdgs,&
       antiquark_pdgs,quark_legs,antiquark_legs)
    implicit none
    integer,dimension(next),intent(in) :: process,order
    integer,intent(out) :: nlines
    integer,dimension(3),intent(out) :: quark_pdgs,antiquark_pdgs,quark_legs,antiquark_legs
    integer :: i,leg,pdg,nquarks,nantiquarks
    quark_pdgs=0
    antiquark_pdgs=0
    quark_legs=0
    antiquark_legs=0
    nquarks=0
    nantiquarks=0
    do i=1,next
       leg=order(i)
       pdg=process(leg)
       if (leg.le.2) pdg=phys_model%get_antipart(pdg)
       if (phys_model%is_quark(pdg)) then
          nquarks=nquarks+1
          if (nquarks.gt.size(quark_pdgs)) then
             write (*,*) 'More than three quark lines are not supported'
             stop 1
          endif
          quark_pdgs(nquarks)=pdg
          quark_legs(nquarks)=leg
       elseif (phys_model%is_antiquark(pdg)) then
          nantiquarks=nantiquarks+1
          if (nantiquarks.gt.size(antiquark_pdgs)) then
             write (*,*) 'More than three quark lines are not supported'
             stop 1
          endif
          antiquark_pdgs(nantiquarks)=pdg
          antiquark_legs(nantiquarks)=leg
       endif
    enddo
    if (nquarks.ne.nantiquarks) then
       write (*,*) 'Unbalanced quark endpoints in process',process
       stop 1
    endif
    nlines=nquarks
  end subroutine collect_quark_endpoints

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
    integer :: i,iord,aq,iaq,ipart
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
