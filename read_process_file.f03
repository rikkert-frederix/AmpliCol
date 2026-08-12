module read_process_file
  use handling_processes
  use run_parameters, only: ignore_final_state_width_fix
  use coupling_orders, only: reset_coupling_order_selection,set_coupling_order_selection,&
       coupling_order_selection,coupling_order_mode_automatic
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
    call read_process_header(buffer,n_external,n_unique,version)
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
    character(len=80) :: filename
    integer :: iproc,igroup,icheck,nproc_in_group,max_channels,iflav,ndim,process_file_version
    real(kind=8) :: idenCOfactor
    integer,dimension(:),allocatable :: process,order,ichans,phase_space_orders,phase_permutation
    integer,dimension(:,:),allocatable :: channel_permutations
    character(len=65536) :: buff
    integer :: i,j,ios
    open(unit=10,file=filename,status='old')
    read(10,'(a)',iostat=ios) buff
    if (ios.ne.0) then
       write (*,*) 'Could not read the process-file header'
       stop 1
    endif
    call read_process_header(buff,next,nproc_unique,process_file_version)
    ndim=3*(next-2)-4
    if (include_pdf) ndim=ndim+2
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
       read(10,*) icheck,nproc_in_group,max_channels,phase_space_orders(1:next)
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
                if (process_file_version.ge.2) then
                   write (*,*) 'Malformed subprocess row; phase-space maps are required for process-file version',&
                        process_file_version
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
       allocate(pgl(igroup)%val_procs_abs(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
       allocate(pgl(igroup)%alias_factors(1:maxval(iden_iproc(1:pgl(igroup)%nproc)),1:pgl(igroup)%nproc))
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
  end subroutine read_processes_from_file


  subroutine read_process_header(buffer,n_external,n_unique,version)
    implicit none
    character(len=*),intent(in) :: buffer
    integer,intent(out) :: n_external,n_unique,version
    integer :: ios,mode,as_min2,as_max2,aew_min2,aew_max2
    logical :: valid_selection

    version=1
    read(buffer,*,iostat=ios) n_external,n_unique,version
    if (ios.ne.0) then
       version=1
       read(buffer,*,iostat=ios) n_external,n_unique
       if (ios.ne.0) then
          write (*,*) 'Malformed process-file header'
          stop 1
       endif
    endif
    if (n_external.lt.3 .or. n_unique.lt.1) then
       write (*,*) 'Invalid process counts in process-file header',n_external,n_unique
       stop 1
    endif
    if (version.lt.1 .or. version.gt.3) then
       write (*,*) 'Unsupported process-file version',version
       stop 1
    endif

    if (version.le.2) then
       ! Historical files contained only the leading-QCD contribution. Their
       ! compatible interpretation is therefore the automatic maximum-aS mode.
       call reset_coupling_order_selection()
       return
    endif

    read(buffer,*,iostat=ios) n_external,n_unique,version,mode,as_min2,as_max2,&
         aew_min2,aew_max2
    if (ios.ne.0) then
       write (*,*) 'Malformed version-3 process-file header; coupling-order metadata is required'
       stop 1
    endif
    call set_coupling_order_selection(mode,as_min2,as_max2,aew_min2,aew_max2,&
         valid_selection)
    if (.not.valid_selection) then
       write (*,*) 'Invalid coupling-order selection in process-file header',&
            mode,as_min2,as_max2,aew_min2,aew_max2
       stop 1
    endif
  end subroutine read_process_header



  subroutine check_unique_processes()
    use phase_space_gen23_mod
    use cuts
    implicit none
    integer :: i,j,iproc,ih,max_as_ew_order
    integer,parameter :: nevent=10
    real(kind=8),dimension(:,:),allocatable :: amp2
    complex(kind=8),dimension(:,:,:),allocatable :: amp_by_order
    real(kind=8),dimension(:),allocatable :: mass,width
    real(kind=8),dimension(pgl_unique%ndim) :: x
    real(kind=8),external :: ran2
    type(psv) :: ps
    allocate(phase_space_gen23 :: pgl_unique%phase_space)
    allocate(pgl_unique%processes(next,nproc_unique))
    allocate(pgl_unique%color_orders(next,nproc_unique))
    allocate(pgl_unique%phase_space_orders(next))
    allocate(pgl_unique%amps(1))
    allocate(mass(next))
    allocate(width(next))
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
    
    max_as_ew_order=compact_unique_max_as_ew_order()
    if (max_as_ew_order.ge.0) then
       call pgl_unique%amps(1)%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
            pgl_unique%spin,pgl_unique%color_orders,phys_model,&
            max_as_ew_order=max_as_ew_order)
    else
       call pgl_unique%amps(1)%init(1,next,pgl_unique%nproc,pgl_unique%processes,&
            pgl_unique%spin,pgl_unique%color_orders,phys_model)
    endif
        
    allocate(amp2(nevent,pgl_unique%nproc))
    allocate(amp_by_order(nevent,pgl_unique%amps(1)%n_amps,&
         pgl_unique%amps(1)%n_sectors))
    allocate(pgl_unique%passed(1))
    allocate(pgl_unique%hel(next))

    pgl_unique%passed=0
    pgl_unique%hel=0
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
       if (use_real_gluons .and. all(pgl_unique%amps(1)%n_qqbar(1:pgl_unique%amps(1)%nprocs).eq.0)) then
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(pgl_unique%passed(1),iproc)=amp2(pgl_unique%passed(1),iproc)+&
                  pgl_unique%amps(1)%amps_r(ih)*pgl_unique%amps(1)%amps_r(ih)
          enddo
       else
          do ih=1,pgl_unique%amps(1)%n_amps
             do while (pgl_unique%amps(1)%iproc_start(iproc+1).eq.ih) ; iproc=iproc+1 ; enddo
             amp2(pgl_unique%passed(1),iproc)=amp2(pgl_unique%passed(1),iproc)+&
                  dble(pgl_unique%amps(1)%amps(ih)*dconjg(pgl_unique%amps(1)%amps(ih)))
          enddo
       endif
       amp_by_order(pgl_unique%passed(1),:,:)=pgl_unique%amps(1)%amps_by_order
       call find_same_flavour(pgl_unique,nevent,amp2(pgl_unique%passed(1),:))
    enddo
    allocate(unique_map(1:pgl_unique%nproc))
    allocate(unique_map_value(1:pgl_unique%nproc))
    call find_unique(pgl_unique,nevent,amp_by_order,unique_map,unique_map_value)

    do iproc=1,pgl_unique%nproc
       write (99,*) unique_map(iproc),unique_map_value(iproc),':',pgl_unique%processes(1:pgl_unique%next,iproc),&
            ':',pgl_unique%color_orders(1:pgl_unique%next,iproc)
    enddo
    deallocate(pgl_unique%spin)
    deallocate(pgl_unique%phase_space)
    deallocate(ps%p)
    deallocate(ps%x)
    deallocate(amp2)
    deallocate(amp_by_order)
  end subroutine check_unique_processes

  integer function compact_unique_max_as_ew_order()
    ! unique_procs is still in the all-outgoing canonical convention.  Return
    ! the EW power of a compact automatic leading-QCD construction, or -1
    ! when unique-process testing needs the complete route history.
    implicit none
    integer :: iproc,flavour,process_ew_order,candidate_order
    integer :: nquarks,nquark_lines,candidate_quark_lines
    integer :: n_wplus,n_wminus,n_charged_vectors,n_higgs
    integer :: n_leptons
    integer,dimension(6) :: flavour_delta
    logical :: compatible_flavour_change

    compact_unique_max_as_ew_order=-1
    if (coupling_order_selection%mode.ne.coupling_order_mode_automatic) return
    candidate_order=-1
    candidate_quark_lines=-1
    do iproc=1,nproc_unique
       nquarks=count(abs(unique_procs(:,iproc)).ge.1 .and. &
            abs(unique_procs(:,iproc)).le.6)
       if (mod(nquarks,2).ne.0) return
       nquark_lines=nquarks/2
       if (nquark_lines.gt.3) return
       if (candidate_quark_lines.lt.0) then
          candidate_quark_lines=nquark_lines
       elseif (nquark_lines.ne.candidate_quark_lines) then
          return
       endif
       if (any(.not.((abs(unique_procs(:,iproc)).ge.1 .and. &
            abs(unique_procs(:,iproc)).le.6) .or. &
            unique_procs(:,iproc).eq.21 .or. unique_procs(:,iproc).eq.22 .or. &
            unique_procs(:,iproc).eq.23 .or. abs(unique_procs(:,iproc)).eq.24 .or. &
            unique_procs(:,iproc).eq.25 .or. &
            (abs(unique_procs(:,iproc)).ge.11 .and. &
            abs(unique_procs(:,iproc)).le.16)))) then
          return
       endif
       if (nquark_lines.eq.0 .and. any(unique_procs(:,iproc).ne.21)) return
       n_wplus=count(unique_procs(:,iproc).eq.24)
       n_wminus=count(unique_procs(:,iproc).eq.-24)
       n_charged_vectors=n_wplus+n_wminus
       n_higgs=count(unique_procs(:,iproc).eq.25)
       n_leptons=count(abs(unique_procs(:,iproc)).ge.11 .and. &
            abs(unique_procs(:,iproc)).le.16)
       if (n_leptons.gt.0 .and. (n_charged_vectors.gt.0 .or. n_higgs.gt.0)) return
       do flavour=11,16
          if (count(unique_procs(:,iproc).eq.flavour).ne.&
               count(unique_procs(:,iproc).eq.-flavour)) return
       enddo
       if (n_higgs.gt.1 .or. &
            (n_charged_vectors.gt.0 .and. n_higgs.gt.0)) return
       process_ew_order=count(unique_procs(:,iproc).eq.22 .or. &
            unique_procs(:,iproc).eq.23 .or. abs(unique_procs(:,iproc)).eq.24 .or. &
            unique_procs(:,iproc).eq.25 .or. &
            (abs(unique_procs(:,iproc)).ge.11 .and. &
            abs(unique_procs(:,iproc)).le.16))
       if (candidate_order.lt.0) candidate_order=process_ew_order
       if (process_ew_order.ne.candidate_order) return
       do flavour=1,6
          flavour_delta(flavour)=count(unique_procs(:,iproc).eq.flavour)-&
               count(unique_procs(:,iproc).eq.-flavour)
       enddo
       if (n_higgs.eq.1 .and. (phys_model%get_mass(6).eq.0d0 .or. &
            count(abs(unique_procs(:,iproc)).eq.6).lt.2)) then
          return
       endif
       if (n_higgs.eq.1 .and. nquark_lines.eq.2) then
          if (count(abs(unique_procs(:,iproc)).eq.6).ne.2) return
       endif
       if (n_charged_vectors.eq.0) then
          if (any(flavour_delta.ne.0)) then
             return
          endif
       else
          compatible_flavour_change=.true.
          do flavour=1,5,2
             if (flavour_delta(flavour+1).ne.-flavour_delta(flavour)) &
                  compatible_flavour_change=.false.
          enddo
          if (sum(max(flavour_delta(1:5:2),0)).gt.n_wplus) &
               compatible_flavour_change=.false.
          if (sum(max(-flavour_delta(1:5:2),0)).gt.n_wminus) &
               compatible_flavour_change=.false.
          if (sum(flavour_delta(1:5:2)).ne.n_wplus-n_wminus) &
               compatible_flavour_change=.false.
          if (.not.compatible_flavour_change) return
       endif
    enddo
    compact_unique_max_as_ew_order=candidate_order
  end function compact_unique_max_as_ew_order

  subroutine find_unique(pgl,nevent,amp_by_order,unique_map,unique_map_value)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,intent(in) :: nevent
    complex(kind=8),dimension(nevent,pgl%amps(1)%n_amps,&
         pgl%amps(1)%n_sectors),intent(in) :: amp_by_order
    real(kind=8),dimension(pgl%nproc),intent(out) :: unique_map_value
    integer,dimension(pgl%nproc),intent(out) :: unique_map
    integer :: i,j,k,first_i,last_i,first_j,last_j,ievent,iamp,isector
    complex(kind=8) :: proportionality
    real(kind=8) :: scale_i,scale_j,comparison_scale,component_scale
    real(kind=8),parameter :: tiny=1d-10
    logical :: compatible
    unique_map=-1d0
    do i=1,pgl%nproc
       first_i=pgl%amps(1)%iproc_start(i)
       last_i=pgl%amps(1)%iproc_start(i+1)-1
       scale_i=maxval(abs(amp_by_order(:,first_i:last_i,:)))
       if (scale_i.le.1d-300) then
          unique_map(i)=0
          unique_map_value(i)=0d0
          cycle
       endif
       do j=1,i-1
          first_j=pgl%amps(1)%iproc_start(j)
          last_j=pgl%amps(1)%iproc_start(j+1)-1
          if (last_i-first_i.ne.last_j-first_j) cycle
          scale_j=maxval(abs(amp_by_order(:,first_j:last_j,:)))
          if (scale_j.le.1d-300) cycle
          ! A numerical relation found on the deliberately massless probe
          ! point is not a valid runtime reduction if it moves a massive
          ! species to another external leg.
          do k=1,pgl%next
             if (phys_model%get_mass(pgl%processes(k,i)).ne.&
                  phys_model%get_mass(pgl%processes(k,j))) exit
          enddo
          if (k.le.pgl%next) cycle
          ! A representative substitution must remain valid after arbitrary
          ! squared-order restrictions and alpha_s running.  Require one
          ! common complex amplitude factor across every probe point,
          ! helicity/colour slot, and coupling sector.  Equality of only the
          ! coherent unit-coupling square is insufficient for this purpose.
          proportionality=(0d0,0d0)
          comparison_scale=0d0
          do ievent=1,nevent
             do iamp=0,last_i-first_i
                do isector=1,pgl%amps(1)%n_sectors
                   if (abs(amp_by_order(ievent,first_j+iamp,isector)).gt.&
                        comparison_scale) then
                      comparison_scale=abs(amp_by_order(ievent,first_j+iamp,isector))
                      proportionality=amp_by_order(ievent,first_i+iamp,isector)/&
                           amp_by_order(ievent,first_j+iamp,isector)
                   endif
                enddo
             enddo
          enddo
          if (comparison_scale.le.1d-300) cycle
          compatible=.true.
          do iamp=0,last_i-first_i
             do isector=1,pgl%amps(1)%n_sectors
                if (pgl%amps(1)%sector_present(first_i+iamp,isector).neqv.&
                     pgl%amps(1)%sector_present(first_j+iamp,isector)) then
                   compatible=.false.
                   exit
                endif
                component_scale=max(&
                     maxval(abs(amp_by_order(:,first_i+iamp,isector))),&
                     abs(proportionality)*&
                     maxval(abs(amp_by_order(:,first_j+iamp,isector))))
                if (component_scale.le.1d-300) cycle
                if (maxval(abs(amp_by_order(:,first_i+iamp,isector)-&
                     proportionality*amp_by_order(:,first_j+iamp,isector))).gt.&
                     tiny*component_scale) then
                   compatible=.false.
                   exit
                endif
             enddo
             if (.not.compatible) exit
          enddo
          if (.not.compatible) cycle
          unique_map_value(i)=abs(proportionality)**2
          unique_map(i)=j
          exit
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
