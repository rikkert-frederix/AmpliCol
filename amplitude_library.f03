module amplitude_library
  use handling_processes
  use read_process_file
  use amplitude_QCD_mod, only: max_library_external => max_amplitude_external_particles
  use pdf_wrap, only: set_ipdgs_for_PDF
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: iostat_end
contains
  subroutine create_amplitude_lib()
    implicit none
    character(len=1024) :: tmp,line,filename
    character(len=256) :: io_message
    integer :: igroup,j,iamp,ios
    real(kind=8),dimension(9) :: stored_model_signature
    if (.not.allocated(pgl_unique) .or. .not.allocated(pgl)) then
       write (*,*) 'Cannot create an amplitude library before process initialisation'
       stop 1
    endif
    if (pgl_unique%next.lt.4 .or. pgl_unique%next.gt.max_library_external) then
       write (*,*) 'Unsupported amplitude-library particle multiplicity:',&
            pgl_unique%next,max_library_external
       stop 1
    endif
    call validate_library_write_state()
    filename='Library/amplib.f03'
    open(unit=14,file=filename,status='replace',action='write',iostat=ios,&
         iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not create amplitude-library dispatcher: ',trim(filename),&
            trim(io_message)
       stop 1
    endif
    call write_dispatcher('module amp_lib')
    call write_dispatcher('use FeynmanRules, only: invalidate_feynman_point')
    do igroup=1,ngroups
       do j=1,size(pgl(igroup)%amps)
          write(tmp,*) igroup
          line=trim(adjustl(tmp))//'_'
          write(tmp,*) j
          line=trim(adjustl(line))//trim(adjustl(tmp))
          call write_dispatcher('use amp'//trim(adjustl(line))//'_lib')
       enddo
    enddo
    call write_dispatcher('implicit none')
    call write_dispatcher('private')
    call write_dispatcher('public :: evaluate_amp')
    call write_dispatcher('contains')
    call write_dispatcher('subroutine evaluate_amp(ichan,iint,p,amps)')
    call write_dispatcher('implicit none')
    call write_dispatcher('integer,intent(in) :: ichan,iint')
    write(tmp,*) maxval(pgl(:)%next)
    call write_dispatcher('real(kind=8),dimension(0:3,'//&
         trim(adjustl(tmp))//'),intent(in) :: p')
    call write_dispatcher('complex(kind=8),intent(out) :: amps(*)')
    call write_dispatcher('logical :: dispatched')
    call write_dispatcher('dispatched=.false.')
    do igroup=1,ngroups
       if (igroup.eq.1) then
          call write_dispatcher('if (ichan.eq.1) then')
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                call write_dispatcher('if (iint.eq.1) then')
                call write_dispatcher('call evaluate_amp1_1(p,amps)')
                call write_dispatcher('dispatched=.true.')
             else
                write(tmp,*) j
                call write_dispatcher('elseif (iint.eq.'//trim(adjustl(tmp))//') then')
                call write_dispatcher('call evaluate_amp1_'//trim(adjustl(tmp))//'(p,amps)')
                call write_dispatcher('dispatched=.true.')
             endif
          enddo
          call write_dispatcher('endif')
       else
          write(tmp,*) igroup
          call write_dispatcher('elseif (ichan.eq.'//trim(adjustl(tmp))//') then')
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                call write_dispatcher('if (iint.eq.1) then')
                call write_dispatcher('call evaluate_amp'//trim(adjustl(tmp))//'_1(p,amps)')
                call write_dispatcher('dispatched=.true.')
             else
                write(line,*) j
                call write_dispatcher('elseif (iint.eq.'//trim(adjustl(line))//') then')
                call write_dispatcher('call evaluate_amp'//trim(adjustl(tmp))//'_'//&
                     trim(adjustl(line))//'(p,amps)')
                call write_dispatcher('dispatched=.true.')
             endif
          enddo
          call write_dispatcher('endif')
       endif
    enddo
    call write_dispatcher('endif')
    call write_dispatcher('if (.not.dispatched) call invalidate_feynman_point()')
    call write_dispatcher('end subroutine evaluate_amp')
    call write_dispatcher('end module amp_lib')
    close(14,iostat=ios,iomsg=io_message)
    call require_library_io('closing amplitude-library dispatcher')
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='replace',&
         action='write',iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not create amplitude-library metadata: ',trim(filename),&
            trim(io_message)
       stop 1
    endif
    write(14,iostat=ios,iomsg=io_message) 4 ! binary amplitude-library format version
    call require_library_io('writing metadata format version')
    stored_model_signature=phys_model%model_signature()
    if (.not.all(ieee_is_finite(stored_model_signature))) then
       write (*,*) 'Cannot store a non-finite amplitude-library model signature'
       stop 1
    endif
    write(14,iostat=ios,iomsg=io_message) stored_model_signature
    call require_library_io('writing model signature')
    write(14,iostat=ios,iomsg=io_message) pgl_unique%next,pgl_unique%nproc
    call require_library_io('writing unique-process dimensions')
    write(14,iostat=ios,iomsg=io_message) unique_map
    call require_library_io('writing unique-process map')
    write(14,iostat=ios,iomsg=io_message) unique_map_value
    call require_library_io('writing unique-process factors')
    write(14,iostat=ios,iomsg=io_message) pgl_unique%processes
    call require_library_io('writing unique processes')
    write(14,iostat=ios,iomsg=io_message) ngroups
    call require_library_io('writing process-group count')
    do igroup=1,ngroups
       ! amplitudes
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%amps)
       call require_library_io('writing amplitude count')
       do iamp=1,size(pgl(igroup)%amps)
          write(14,iostat=ios,iomsg=io_message) pgl(igroup)%amps(iamp)%n_amps
          call require_library_io('writing matrix-element count')
          write(14,iostat=ios,iomsg=io_message) &
               size(pgl(igroup)%amps(iamp)%iproc_start)
          call require_library_io('writing subprocess-offset extent')
          write(14,iostat=ios,iomsg=io_message) pgl(igroup)%amps(iamp)%iproc_start
          call require_library_io('writing subprocess offsets')
          write(14,iostat=ios,iomsg=io_message) &
               size(pgl(igroup)%amps(iamp)%n_sing),pgl(igroup)%amps(iamp)%n_sing
          call require_library_io('writing singlet counts')
          write(14,iostat=ios,iomsg=io_message) &
               size(pgl(igroup)%amps(iamp)%n_qqbar),pgl(igroup)%amps(iamp)%n_qqbar
          call require_library_io('writing quark-line counts')
          write(14,iostat=ios,iomsg=io_message) &
               shape(pgl(igroup)%amps(iamp)%spins),pgl(igroup)%amps(iamp)%spins
          call require_library_io('writing helicity metadata')
          write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%amps(iamp)%amps)
          call require_library_io('writing amplitude-workspace extent')
       enddo
       ! multichannel
       write(14,iostat=ios,iomsg=io_message) &
            pgl(igroup)%multichan%n_unique_channels,&
            pgl(igroup)%multichan%n_unique_channelgroups
       call require_library_io('writing unique-channel counts')
       write(14,iostat=ios,iomsg=io_message) &
            shape(pgl(igroup)%multichan%unique_channelgroup_list),&
            pgl(igroup)%multichan%unique_channelgroup_list
       call require_library_io('writing unique-channel groups')
       write(14,iostat=ios,iomsg=io_message) &
            size(pgl(igroup)%multichan%unique_channel_list),&
            pgl(igroup)%multichan%unique_channel_list
       call require_library_io('writing unique-channel list')
       write(14,iostat=ios,iomsg=io_message) &
            size(pgl(igroup)%multichan%map_proc_to_channelgroup),&
            pgl(igroup)%multichan%map_proc_to_channelgroup
       call require_library_io('writing process-channel map')
       write(14,iostat=ios,iomsg=io_message) &
            size(pgl(igroup)%multichan%number_of_channels),&
            pgl(igroup)%multichan%number_of_channels
       call require_library_io('writing channel counts')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%multichan%channels),&
            pgl(igroup)%multichan%channels
       call require_library_io('writing channel list')
       write(14,iostat=ios,iomsg=io_message) &
            shape(pgl(igroup)%multichan%channel_permutations),&
            pgl(igroup)%multichan%channel_permutations
       call require_library_io('writing channel permutations')
       ! rest
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%processes),&
            pgl(igroup)%processes
       call require_library_io('writing process array')
       write(14,iostat=ios,iomsg=io_message) &
            shape(pgl(igroup)%phase_space_permutations),&
            pgl(igroup)%phase_space_permutations
       call require_library_io('writing phase-space permutations')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%iden_iproc),&
            pgl(igroup)%iden_iproc
       call require_library_io('writing identical-process counts')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%phase_space_orders),&
            pgl(igroup)%phase_space_orders
       call require_library_io('writing phase-space order')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%nhel),&
            pgl(igroup)%nhel
       call require_library_io('writing helicity counts')
       write(14,iostat=ios,iomsg=io_message) pgl(igroup)%nproc
       call require_library_io('writing process count')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%val_procs),&
            pgl(igroup)%val_procs
       call require_library_io('writing process values')
       write(14,iostat=ios,iomsg=io_message) &
            shape(pgl(igroup)%idenCOandMAPfactor),&
            pgl(igroup)%idenCOandMAPfactor
       call require_library_io('writing identical-process factors')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%iden_processes),&
            pgl(igroup)%iden_processes
       call require_library_io('writing identical processes')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%spin),&
            pgl(igroup)%spin
       call require_library_io('writing spin table')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%hel_fac),&
            pgl(igroup)%hel_fac
       call require_library_io('writing helicity factors')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%iden),&
            pgl(igroup)%iden
       call require_library_io('writing identity array')
       write(14,iostat=ios,iomsg=io_message) pgl(igroup)%ipdgs
       call require_library_io('writing PDF flavour mask')
       write(14,iostat=ios,iomsg=io_message) pgl(igroup)%next,pgl(igroup)%ndim
       call require_library_io('writing group dimensions')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%col_fac),&
            pgl(igroup)%col_fac
       call require_library_io('writing colour factors')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%amp2)
       call require_library_io('writing matrix-element-square extent')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%amp2_hel)
       call require_library_io('writing helicity-square extent')
       write(14,iostat=ios,iomsg=io_message) size(pgl(igroup)%passed),&
            pgl(igroup)%passed
       call require_library_io('writing passed counts')
       write(14,iostat=ios,iomsg=io_message) shape(pgl(igroup)%color_orders),&
            pgl(igroup)%color_orders
       call require_library_io('writing colour orders')
    enddo
    close(14,iostat=ios,iomsg=io_message)
    call require_library_io('closing amplitude-library metadata')

  contains

    subroutine write_dispatcher(source_line)
      implicit none
      character(len=*),intent(in) :: source_line
      write(14,'(a)',iostat=ios,iomsg=io_message) trim(source_line)
      call require_library_io('writing amplitude-library dispatcher')
    end subroutine write_dispatcher

    subroutine require_library_io(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Could not finish ',trim(label),': ',ios,trim(io_message)
         stop 1
      endif
    end subroutine require_library_io

    subroutine invalid_library_write_state(label,group_index,amplitude_index)
      implicit none
      character(len=*),intent(in) :: label
      integer,intent(in),optional :: group_index,amplitude_index
      if (present(amplitude_index)) then
         if (.not.present(group_index)) then
            write (*,*) 'Invalid amplitude-library state diagnostic: ',trim(label)
            stop 1
         endif
         write (*,*) 'Invalid amplitude-library state before writing: ',trim(label),&
              group_index,amplitude_index
      elseif (present(group_index)) then
         write (*,*) 'Invalid amplitude-library state before writing: ',trim(label),&
              group_index
      else
         write (*,*) 'Invalid amplitude-library state before writing: ',trim(label)
      endif
      stop 1
    end subroutine invalid_library_write_state

    subroutine validate_library_write_state()
      implicit none
      integer :: group_index,amplitude_index,n_amplitudes,n_processes

      if (ngroups.lt.1 .or. size(pgl).ne.ngroups) &
           call invalid_library_write_state('inconsistent process-group count')
      if (pgl_unique%next.lt.4 .or. &
           pgl_unique%next.gt.max_library_external .or. &
           pgl_unique%nproc.lt.1) &
           call invalid_library_write_state('invalid unique-process dimensions')
      if (.not.allocated(unique_map) .or. .not.allocated(unique_map_value) .or. &
           .not.allocated(pgl_unique%processes)) &
           call invalid_library_write_state('missing unique-process metadata')
      if (size(unique_map).ne.pgl_unique%nproc .or. &
           size(unique_map_value).ne.pgl_unique%nproc .or. &
           size(pgl_unique%processes,1).ne.pgl_unique%next .or. &
           size(pgl_unique%processes,2).ne.pgl_unique%nproc) &
           call invalid_library_write_state('incompatible unique-process metadata')
      if (.not.all(ieee_is_finite(unique_map_value)) .or. &
           any(unique_map.lt.-1) .or. any(unique_map.gt.pgl_unique%nproc)) &
           call invalid_library_write_state('invalid unique-process values')

      do group_index=1,ngroups
         if (pgl(group_index)%next.ne.pgl_unique%next .or. &
              pgl(group_index)%next.lt.4 .or. &
              pgl(group_index)%next.gt.max_library_external .or. &
              pgl(group_index)%nproc.lt.1) &
              call invalid_library_write_state('invalid group dimensions',group_index)
         if (.not.allocated(pgl(group_index)%amps) .or. &
              .not.allocated(pgl(group_index)%processes) .or. &
              .not.allocated(pgl(group_index)%phase_space_permutations) .or. &
              .not.allocated(pgl(group_index)%iden_iproc) .or. &
              .not.allocated(pgl(group_index)%phase_space_orders) .or. &
              .not.allocated(pgl(group_index)%nhel) .or. &
              .not.allocated(pgl(group_index)%val_procs) .or. &
              .not.allocated(pgl(group_index)%idenCOandMAPfactor) .or. &
              .not.allocated(pgl(group_index)%iden_processes) .or. &
              .not.allocated(pgl(group_index)%spin) .or. &
              .not.allocated(pgl(group_index)%hel_fac) .or. &
              .not.allocated(pgl(group_index)%iden) .or. &
              .not.allocated(pgl(group_index)%col_fac) .or. &
              .not.allocated(pgl(group_index)%amp2) .or. &
              .not.allocated(pgl(group_index)%amp2_hel) .or. &
              .not.allocated(pgl(group_index)%passed) .or. &
              .not.allocated(pgl(group_index)%color_orders)) &
              call invalid_library_write_state('missing process-group metadata',&
              group_index)
         n_amplitudes=size(pgl(group_index)%amps)
         n_processes=pgl(group_index)%nproc
         if (n_amplitudes.lt.1 .or. &
              size(pgl(group_index)%processes,1).ne.pgl(group_index)%next .or. &
              size(pgl(group_index)%processes,2).ne.n_processes .or. &
              any(shape(pgl(group_index)%phase_space_permutations).ne.&
              shape(pgl(group_index)%processes)) .or. &
              any(shape(pgl(group_index)%color_orders).ne.&
              shape(pgl(group_index)%processes)) .or. &
              size(pgl(group_index)%iden_iproc).ne.n_processes .or. &
              size(pgl(group_index)%phase_space_orders).ne.pgl(group_index)%next .or. &
              size(pgl(group_index)%nhel).ne.n_amplitudes .or. &
              size(pgl(group_index)%passed).ne.n_amplitudes .or. &
              size(pgl(group_index)%col_fac).ne.n_processes .or. &
              size(pgl(group_index)%iden).ne.n_processes) &
              call invalid_library_write_state('incompatible process-group extents',&
              group_index)
         if (size(pgl(group_index)%val_procs,2).ne.n_processes .or. &
              size(pgl(group_index)%idenCOandMAPfactor,2).ne.n_processes .or. &
              size(pgl(group_index)%iden_processes,1).ne.pgl(group_index)%next .or. &
              size(pgl(group_index)%iden_processes,3).ne.n_processes .or. &
              size(pgl(group_index)%spin,1).ne.4 .or. &
              size(pgl(group_index)%spin,2).ne.pgl(group_index)%next .or. &
              size(pgl(group_index)%hel_fac,1).ne.size(pgl(group_index)%amp2_hel) .or. &
              (size(pgl(group_index)%hel_fac,2).ne.n_amplitudes .and. &
              size(pgl(group_index)%hel_fac,2).ne.n_processes)) &
              call invalid_library_write_state('incompatible group work arrays',&
              group_index)
         if (any(pgl(group_index)%iden_iproc.lt.1) .or. &
              any(pgl(group_index)%nhel.lt.1) .or. &
              any(pgl(group_index)%passed.lt.0) .or. &
              any(pgl(group_index)%col_fac.lt.1) .or. &
              any(pgl(group_index)%iden.lt.1_8) .or. &
              .not.all(ieee_is_finite(pgl(group_index)%val_procs)) .or. &
              .not.all(ieee_is_finite(pgl(group_index)%idenCOandMAPfactor))) &
              call invalid_library_write_state('invalid group values',group_index)
         if (size(pgl(group_index)%val_procs,1).lt.&
              maxval(pgl(group_index)%iden_iproc) .or. &
              size(pgl(group_index)%idenCOandMAPfactor,1).lt.&
              maxval(pgl(group_index)%iden_iproc) .or. &
              size(pgl(group_index)%iden_processes,2).lt.&
              maxval(pgl(group_index)%iden_iproc)) &
              call invalid_library_write_state('short identical-process metadata',&
              group_index)

         if (.not.allocated(pgl(group_index)%multichan%unique_channelgroup_list) .or. &
              .not.allocated(pgl(group_index)%multichan%unique_channel_list) .or. &
              .not.allocated(pgl(group_index)%multichan%map_proc_to_channelgroup) .or. &
              .not.allocated(pgl(group_index)%multichan%number_of_channels) .or. &
              .not.allocated(pgl(group_index)%multichan%channels) .or. &
              .not.allocated(pgl(group_index)%multichan%channel_permutations)) &
              call invalid_library_write_state('missing multichannel metadata',&
              group_index)
         if (pgl(group_index)%multichan%n_unique_channels.lt.1 .or. &
              pgl(group_index)%multichan%n_unique_channelgroups.lt.1 .or. &
              size(pgl(group_index)%multichan%unique_channel_list).ne.&
              pgl(group_index)%multichan%n_unique_channels .or. &
              size(pgl(group_index)%multichan%unique_channelgroup_list,2).ne.&
              pgl(group_index)%multichan%n_unique_channelgroups .or. &
              size(pgl(group_index)%multichan%map_proc_to_channelgroup).ne.&
              n_processes .or. &
              size(pgl(group_index)%multichan%number_of_channels).ne.n_processes .or. &
              size(pgl(group_index)%multichan%channels,2).ne.n_processes .or. &
              size(pgl(group_index)%multichan%channel_permutations,1).ne.&
              pgl(group_index)%next .or. &
              size(pgl(group_index)%multichan%channel_permutations,2).ne.&
              size(pgl(group_index)%multichan%channels,1) .or. &
              size(pgl(group_index)%multichan%channel_permutations,3).ne.&
              n_processes) &
              call invalid_library_write_state('incompatible multichannel metadata',&
              group_index)

         do amplitude_index=1,n_amplitudes
            if (pgl(group_index)%amps(amplitude_index)%n_amps.lt.1) &
                 call invalid_library_write_state('invalid matrix-element count',&
                 group_index,amplitude_index)
            if (.not.allocated(pgl(group_index)%amps(amplitude_index)%iproc_start) .or. &
                 .not.allocated(pgl(group_index)%amps(amplitude_index)%n_sing) .or. &
                 .not.allocated(pgl(group_index)%amps(amplitude_index)%n_qqbar) .or. &
                 .not.allocated(pgl(group_index)%amps(amplitude_index)%spins) .or. &
                 .not.allocated(pgl(group_index)%amps(amplitude_index)%amps)) &
                 call invalid_library_write_state('missing amplitude metadata',&
                 group_index,amplitude_index)
            if (size(pgl(group_index)%amps(amplitude_index)%iproc_start).lt.2 .or. &
                 size(pgl(group_index)%amps(amplitude_index)%n_sing).ne.&
                 size(pgl(group_index)%amps(amplitude_index)%iproc_start)-1 .or. &
                 size(pgl(group_index)%amps(amplitude_index)%n_qqbar).ne.&
                 size(pgl(group_index)%amps(amplitude_index)%iproc_start)-1 .or. &
                 size(pgl(group_index)%amps(amplitude_index)%spins,1).ne.&
                 pgl(group_index)%next .or. &
                 size(pgl(group_index)%amps(amplitude_index)%spins,3).ne.&
                 pgl(group_index)%amps(amplitude_index)%n_amps .or. &
                 size(pgl(group_index)%amps(amplitude_index)%amps).ne.&
                 pgl(group_index)%amps(amplitude_index)%n_amps .or. &
                 pgl(group_index)%nhel(amplitude_index).ne.&
                 pgl(group_index)%amps(amplitude_index)%n_amps) &
                 call invalid_library_write_state('incompatible amplitude metadata',&
                 group_index,amplitude_index)
         enddo
      enddo
    end subroutine validate_library_write_state
  end subroutine create_amplitude_lib

  subroutine read_amplitude_lib(metadata_filename)
    implicit none
    character(len=*),intent(in),optional :: metadata_filename
    character(len=1024) :: filename
    integer :: dim1,dim2,dim3,iamp,igroup,library_version,ios,iproc,label,k,&
         stored_nprocs,self_count,partner,ncoloured,nsinglet
    integer(kind=1) :: trailing_byte
    integer,parameter :: max_library_extent=200000000,max_library_groups=100000
    integer(kind=8),parameter :: max_library_workspace_bytes=2147483648_8
    integer(kind=8) :: library_workspace_bytes
    real(kind=8),dimension(9) :: stored_model_signature,current_model_signature
    real(kind=8),dimension(9) :: comparison_scale,comparison_residual
    character(len=256) :: allocation_message
    logical,allocatable :: seen_channels(:),referenced_unique_channels(:)
    integer,allocatable :: channel_marks(:)
    library_workspace_bytes=0_8
    filename='Library/amplitudes.bin'
    if (present(metadata_filename)) filename=trim(metadata_filename)
    open(unit=14,file=filename,form='unformatted',access='stream',status='old',&
         action='read',iostat=ios)
    if (ios.ne.0) then
       write (*,*) 'Could not open amplitude-library metadata: ',trim(filename)
       stop 1
    endif
    read(14,iostat=ios) library_version
    call require_read('format version')
    if (library_version.ne.4) then
       write (*,*) 'Amplitude library has an incompatible binary format; recreate it'
       stop 1
    endif
    read(14,iostat=ios) stored_model_signature
    call require_read('model signature')
    if (.not.all(ieee_is_finite(stored_model_signature))) then
       write (*,*) 'Amplitude library contains a non-finite model signature'
       stop 1
    endif
    call reserve_library_array(1_8,2048_8,'unique-process descriptor',.false.)
    allocate(pgl_unique,stat=ios,errmsg=allocation_message)
    call require_allocation('unique-process descriptor')
    read(14,iostat=ios) pgl_unique%next,pgl_unique%nproc
    call require_read('unique-process dimensions')
    if (pgl_unique%next.lt.4 .or. pgl_unique%next.gt.max_library_external .or. &
         pgl_unique%nproc.lt.1 .or. pgl_unique%nproc.gt.max_library_extent) then
       write (*,*) 'Invalid unique-process dimensions in amplitude library',&
            pgl_unique%next,pgl_unique%nproc
       stop 1
    endif
    next=pgl_unique%next
    call require_shape2(pgl_unique%next,pgl_unique%nproc,'unique-process array shape')
    call reserve_library_array(int(pgl_unique%nproc,kind=8),4_8,&
         'unique-process map',.true.)
    allocate(unique_map(pgl_unique%nproc),stat=ios,errmsg=allocation_message)
    call require_allocation('unique-process map')
    read(14,iostat=ios) unique_map
    call require_read('unique-process map')
    call reserve_library_array(int(pgl_unique%nproc,kind=8),8_8,&
         'unique-process factors',.true.)
    allocate(unique_map_value(pgl_unique%nproc),stat=ios,errmsg=allocation_message)
    call require_allocation('unique-process factors')
    read(14,iostat=ios) unique_map_value
    call require_read('unique-process factors')
    if (.not.all(ieee_is_finite(unique_map_value))) then
       write (*,*) 'Amplitude library contains non-finite unique-process factors'
       stop 1
    endif
    if (any(unique_map.lt.-1) .or. any(unique_map.gt.pgl_unique%nproc)) then
       write (*,*) 'Amplitude library contains invalid unique-process indices'
       stop 1
    endif
    call reserve_library_array(int(pgl_unique%next,kind=8)*&
         int(pgl_unique%nproc,kind=8),4_8,'unique processes',.true.)
    allocate(pgl_unique%processes(pgl_unique%next,pgl_unique%nproc),&
         stat=ios,errmsg=allocation_message)
    call require_allocation('unique processes')
    read(14,iostat=ios) pgl_unique%processes
    call require_read('unique processes')
    read(14,iostat=ios) ngroups
    call require_read('group count')
    if (ngroups.lt.1 .or. ngroups.gt.max_library_groups) then
       write (*,*) 'Invalid group count in amplitude library',ngroups
       stop 1
    endif
    call reserve_library_array(int(ngroups,kind=8),2048_8,&
         'process-group descriptors',.false.)
    allocate(pgl(ngroups),stat=ios,errmsg=allocation_message)
    call require_allocation('process-group descriptors')
    do igroup=1,ngroups
       ! amplitudes
       read(14,iostat=ios) dim1
       call require_read('amplitude count')
       call require_extent(dim1,'amplitude count')
       call reserve_library_array(int(dim1,kind=8),2048_8,&
            'amplitude descriptors',.false.)
       allocate(pgl(igroup)%amps(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('amplitude descriptors')
       do iamp=1,dim1
          read(14,iostat=ios) pgl(igroup)%amps(iamp)%n_amps
          call require_read('matrix-element count')
          call require_extent(pgl(igroup)%amps(iamp)%n_amps,'matrix-element count')
          read(14,iostat=ios) dim1
          call require_read('subprocess-offset extent')
          call require_extent(dim1,'subprocess-offset extent')
          if (dim1.lt.2 .or. dim1.gt.pgl(igroup)%amps(iamp)%n_amps+1) then
             write (*,*) 'Invalid subprocess-offset extent in amplitude library',&
                  igroup,iamp,dim1
             stop 1
          endif
          stored_nprocs=dim1-1
          call reserve_library_array(int(dim1,kind=8),4_8,&
               'subprocess offsets',.true.)
          allocate(pgl(igroup)%amps(iamp)%iproc_start(dim1),&
               stat=ios,errmsg=allocation_message)
          call require_allocation('subprocess offsets')
          read(14,iostat=ios) pgl(igroup)%amps(iamp)%iproc_start
          call require_read('subprocess offsets')
          if (any(pgl(igroup)%amps(iamp)%iproc_start.lt.1) .or. &
               any(pgl(igroup)%amps(iamp)%iproc_start.gt.&
               pgl(igroup)%amps(iamp)%n_amps+1)) then
             write (*,*) 'Invalid subprocess offsets in amplitude library',igroup,iamp
             stop 1
          endif
          read(14,iostat=ios) dim1
          call require_read('singlet-count extent')
          call require_extent(dim1,'singlet-count extent')
          if (dim1.ne.stored_nprocs) then
             write (*,*) 'Inconsistent singlet-count extent in amplitude library',&
                  igroup,iamp,dim1,stored_nprocs
             stop 1
          endif
          call reserve_library_array(int(dim1,kind=8),4_8,'singlet counts',.true.)
          allocate(pgl(igroup)%amps(iamp)%n_sing(dim1),&
               stat=ios,errmsg=allocation_message)
          call require_allocation('singlet counts')
          read(14,iostat=ios) pgl(igroup)%amps(iamp)%n_sing
          call require_read('singlet counts')
          read(14,iostat=ios) dim1
          call require_read('quark-line-count extent')
          call require_extent(dim1,'quark-line-count extent')
          if (dim1.ne.stored_nprocs) then
             write (*,*) 'Inconsistent quark-line-count extent in amplitude library',&
                  igroup,iamp,dim1,stored_nprocs
             stop 1
          endif
          call reserve_library_array(int(dim1,kind=8),4_8,'quark-line counts',.true.)
          allocate(pgl(igroup)%amps(iamp)%n_qqbar(dim1),&
               stat=ios,errmsg=allocation_message)
          call require_allocation('quark-line counts')
          read(14,iostat=ios) pgl(igroup)%amps(iamp)%n_qqbar
          call require_read('quark-line counts')
          read(14,iostat=ios) dim1,dim2,dim3
          call require_read('spin-array shape')
          call require_shape3(dim1,dim2,dim3,'spin-array shape')
          if (dim1.ne.pgl_unique%next .or. &
               dim3.ne.pgl(igroup)%amps(iamp)%n_amps) then
             write (*,*) 'Inconsistent spin-array shape in amplitude library',&
                  igroup,iamp,dim1,dim2,dim3
             stop 1
          endif
          call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8)*&
               int(dim3,kind=8),4_8,'spin array',.true.)
          allocate(pgl(igroup)%amps(iamp)%spins(dim1,dim2,dim3),&
               stat=ios,errmsg=allocation_message)
          call require_allocation('spin array')
          read(14,iostat=ios) pgl(igroup)%amps(iamp)%spins
          call require_read('spin array')
          read(14,iostat=ios) dim1
          call require_read('amplitude-array extent')
          if (dim1.ne.pgl(igroup)%amps(iamp)%n_amps) then
             write (*,*) 'Inconsistent amplitude-array extent in library',igroup,iamp,dim1,&
                  pgl(igroup)%amps(iamp)%n_amps
             stop 1
          endif
          call reserve_library_array(int(dim1,kind=8),16_8,&
               'amplitude workspace',.false.)
          allocate(pgl(igroup)%amps(iamp)%amps(dim1),&
               stat=ios,errmsg=allocation_message)
          call require_allocation('amplitude workspace')
          pgl(igroup)%amps(iamp)%amps=(0d0,0d0)
       enddo
       ! multichannel
       read(14,iostat=ios) pgl(igroup)%multichan%n_unique_channels,&
            pgl(igroup)%multichan%n_unique_channelgroups
       call require_read('unique-channel counts')
       call require_extent(pgl(igroup)%multichan%n_unique_channels,'unique-channel count')
       call require_extent(pgl(igroup)%multichan%n_unique_channelgroups,&
            'unique-channel-group count')
       read(14,iostat=ios) dim1,dim2
       call require_read('unique-channel-group-list shape')
       call require_shape2(dim1,dim2,'unique-channel-group-list shape')
       if (dim2.ne.pgl(igroup)%multichan%n_unique_channelgroups .or. &
            dim1.gt.ngroups+1) then
          write (*,*) 'Inconsistent unique-channel-group-list shape in amplitude library',&
               igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'unique-channel-group list',.true.)
       allocate(pgl(igroup)%multichan%unique_channelgroup_list(0:dim1-1,dim2),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('unique-channel-group list')
       read(14,iostat=ios) pgl(igroup)%multichan%unique_channelgroup_list
       call require_read('unique-channel-group list')
       read(14,iostat=ios) dim1
       call require_read('unique-channel-list extent')
       call require_extent(dim1,'unique-channel-list extent')
       if (dim1.ne.pgl(igroup)%multichan%n_unique_channels) then
          write (*,*) 'Inconsistent unique-channel-list extent in amplitude library',&
               igroup,dim1,pgl(igroup)%multichan%n_unique_channels
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,&
            'unique-channel list',.true.)
       allocate(pgl(igroup)%multichan%unique_channel_list(dim1),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('unique-channel list')
       read(14,iostat=ios) pgl(igroup)%multichan%unique_channel_list
       call require_read('unique-channel list')
       read(14,iostat=ios) dim1
       call require_read('process-channel-map extent')
       call require_extent(dim1,'process-channel-map extent')
       call reserve_library_array(int(dim1,kind=8),4_8,&
            'process-channel map',.true.)
       allocate(pgl(igroup)%multichan%map_proc_to_channelgroup(dim1),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('process-channel map')
       read(14,iostat=ios) pgl(igroup)%multichan%map_proc_to_channelgroup
       call require_read('process-channel map')
       read(14,iostat=ios) dim1
       call require_read('channel-count-array extent')
       call require_extent(dim1,'channel-count-array extent')
       if (dim1.ne.size(pgl(igroup)%multichan%map_proc_to_channelgroup)) then
          write (*,*) 'Inconsistent channel-count-array extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,&
            'channel-count array',.true.)
       allocate(pgl(igroup)%multichan%number_of_channels(dim1),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('channel-count array')
       read(14,iostat=ios) pgl(igroup)%multichan%number_of_channels
       call require_read('channel-count array')
       read(14,iostat=ios) dim1,dim2
       call require_read('channel-array shape')
       call require_shape2(dim1,dim2,'channel-array shape')
       if (dim1.gt.ngroups .or. &
            dim2.ne.size(pgl(igroup)%multichan%number_of_channels)) then
          write (*,*) 'Inconsistent channel-array shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'channel array',.true.)
       allocate(pgl(igroup)%multichan%channels(dim1,dim2),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('channel array')
       read(14,iostat=ios) pgl(igroup)%multichan%channels
       call require_read('channel array')
       pgl(igroup)%multichan%max_channels=dim1
       read(14,iostat=ios) dim1,dim2,dim3
       call require_read('channel-permutation shape')
       call require_shape3(dim1,dim2,dim3,'channel-permutation shape')
       if (dim1.ne.pgl_unique%next .or. &
            dim2.ne.size(pgl(igroup)%multichan%channels,1) .or. &
            dim3.ne.size(pgl(igroup)%multichan%channels,2)) then
          write (*,*) 'Inconsistent channel-permutation shape in amplitude library',&
               igroup,dim1,dim2,dim3
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8)*&
            int(dim3,kind=8),4_8,'channel permutations',.true.)
       allocate(pgl(igroup)%multichan%channel_permutations(dim1,dim2,dim3),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('channel permutations')
       read(14,iostat=ios) pgl(igroup)%multichan%channel_permutations
       call require_read('channel permutations')
       ! rest
       read(14,iostat=ios) dim1,dim2
       call require_read('process-array shape')
       call require_shape2(dim1,dim2,'process-array shape')
       if (dim1.lt.3 .or. dim1.gt.max_library_external) then
          write (*,*) 'Invalid external multiplicity in amplitude library',dim1
          stop 1
       endif
       if (dim1.ne.pgl_unique%next .or. &
            dim2.ne.size(pgl(igroup)%multichan%number_of_channels)) then
          write (*,*) 'Inconsistent process-array shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'process array',.true.)
       allocate(pgl(igroup)%processes(dim1,dim2),stat=ios,errmsg=allocation_message)
       call require_allocation('process array')
       read(14,iostat=ios) pgl(igroup)%processes
       call require_read('process array')
       read(14,iostat=ios) dim1,dim2
       call require_read('phase-space-permutation shape')
       call require_shape2(dim1,dim2,'phase-space-permutation shape')
       if (dim1.ne.size(pgl(igroup)%processes,1) .or. &
            dim2.ne.size(pgl(igroup)%processes,2)) then
          write (*,*) 'Inconsistent phase-space-permutation shape in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'phase-space permutations',.true.)
       allocate(pgl(igroup)%phase_space_permutations(dim1,dim2),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('phase-space permutations')
       read(14,iostat=ios) pgl(igroup)%phase_space_permutations
       call require_read('phase-space permutations')
       if (any(shape(pgl(igroup)%phase_space_permutations).ne.&
            shape(pgl(igroup)%processes))) then
          write (*,*) 'Inconsistent phase-space-permutation shape in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1
       call require_read('identical-process-count extent')
       call require_extent(dim1,'identical-process-count extent')
       if (dim1.ne.size(pgl(igroup)%processes,2)) then
          write (*,*) 'Inconsistent identical-process-count extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,&
            'identical-process counts',.true.)
       allocate(pgl(igroup)%iden_iproc(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('identical-process counts')
       read(14,iostat=ios) pgl(igroup)%iden_iproc
       call require_read('identical-process counts')
       if (any(pgl(igroup)%iden_iproc.lt.1)) then
          write (*,*) 'Invalid identical-process count in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1
       call require_read('phase-space-order extent')
       call require_extent(dim1,'phase-space-order extent')
       if (dim1.ne.size(pgl(igroup)%processes,1)) then
          write (*,*) 'Inconsistent phase-space-order extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,&
            'phase-space order',.true.)
       allocate(pgl(igroup)%phase_space_orders(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('phase-space order')
       read(14,iostat=ios) pgl(igroup)%phase_space_orders
       call require_read('phase-space order')
       if (dim1.ne.size(pgl(igroup)%processes,1)) then
          write (*,*) 'Inconsistent phase-space-order extent in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1
       call require_read('helicity-count extent')
       call require_extent(dim1,'helicity-count extent')
       if (dim1.ne.size(pgl(igroup)%amps)) then
          write (*,*) 'Inconsistent helicity-count extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,'helicity counts',.true.)
       allocate(pgl(igroup)%nhel(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('helicity counts')
       read(14,iostat=ios) pgl(igroup)%nhel
       call require_read('helicity counts')
       if (any(pgl(igroup)%nhel.lt.1)) then
          write (*,*) 'Invalid helicity count in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) pgl(igroup)%nproc
       call require_read('process count')
       call require_extent(pgl(igroup)%nproc,'process count')
       if (pgl(igroup)%nproc.ne.size(pgl(igroup)%processes,2) .or. &
            pgl(igroup)%nproc.ne.size(pgl(igroup)%iden_iproc)) then
          write (*,*) 'Inconsistent process count in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1,dim2
       call require_read('process-value shape')
       call require_shape2(dim1,dim2,'process-value shape')
       if (dim2.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent process-value shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),8_8,&
            'process values',.true.)
       allocate(pgl(igroup)%val_procs(dim1,dim2),stat=ios,errmsg=allocation_message)
       call require_allocation('process values')
       read(14,iostat=ios) pgl(igroup)%val_procs
       call require_read('process values')
       if (.not.all(ieee_is_finite(pgl(igroup)%val_procs))) then
          write (*,*) 'Non-finite process values in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1,dim2
       call require_read('identical-factor shape')
       call require_shape2(dim1,dim2,'identical-factor shape')
       if (dim2.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent identical-factor shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),8_8,&
            'identical factors',.true.)
       allocate(pgl(igroup)%idenCOandMAPfactor(dim1,dim2),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('identical factors')
       read(14,iostat=ios) pgl(igroup)%idenCOandMAPfactor
       call require_read('identical factors')
       if (.not.all(ieee_is_finite(pgl(igroup)%idenCOandMAPfactor))) then
          write (*,*) 'Non-finite identical-process factors in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1,dim2,dim3
       call require_read('identical-process-array shape')
       call require_shape3(dim1,dim2,dim3,'identical-process-array shape')
       ! pgl%next is read later in format v4; the process array is already
       ! authoritative for the external multiplicity at this point.
       if (dim1.ne.size(pgl(igroup)%processes,1) .or. &
            dim3.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent identical-process-array shape in amplitude library',&
               igroup,dim1,dim2,dim3
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8)*&
            int(dim3,kind=8),4_8,'identical processes',.true.)
       allocate(pgl(igroup)%iden_processes(dim1,dim2,dim3),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('identical processes')
       read(14,iostat=ios) pgl(igroup)%iden_processes
       call require_read('identical processes')
       read(14,iostat=ios) dim1,dim2
       call require_read('spin-table shape')
       call require_shape2(dim1,dim2,'spin-table shape')
       if (dim1.ne.4 .or. dim2.ne.size(pgl(igroup)%processes,1)) then
          write (*,*) 'Invalid spin-table shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'spin table',.true.)
       allocate(pgl(igroup)%spin(0:dim1-1,dim2),stat=ios,errmsg=allocation_message)
       call require_allocation('spin table')
       read(14,iostat=ios) pgl(igroup)%spin
       call require_read('spin table')
       read(14,iostat=ios) dim1,dim2
       call require_read('helicity-factor shape')
       call require_shape2(dim1,dim2,'helicity-factor shape')
       if (dim2.ne.size(pgl(igroup)%amps) .and. &
            dim2.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent helicity-factor shape in amplitude library',igroup,dim1,dim2
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'helicity factors',.true.)
       allocate(pgl(igroup)%hel_fac(dim1,dim2),stat=ios,errmsg=allocation_message)
       call require_allocation('helicity factors')
       read(14,iostat=ios) pgl(igroup)%hel_fac
       call require_read('helicity factors')
       if (any(pgl(igroup)%hel_fac.lt.1)) then
          write (*,*) 'Invalid helicity factor in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1
       call require_read('identity-array extent')
       call require_extent(dim1,'identity-array extent')
       if (dim1.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent identity-array extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),8_8,'identity array',.true.)
       allocate(pgl(igroup)%iden(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('identity array')
       read(14,iostat=ios) pgl(igroup)%iden
       call require_read('identity array')
       read(14,iostat=ios) pgl(igroup)%ipdgs
       call require_read('PDF flavour mask')
       read(14,iostat=ios) pgl(igroup)%next,pgl(igroup)%ndim
       call require_read('group dimensions')
       if (pgl(igroup)%next.ne.size(pgl(igroup)%processes,1) .or. &
            pgl(igroup)%next.lt.4 .or. pgl(igroup)%next.gt.max_library_external .or. &
            pgl(igroup)%ndim.lt.0 .or. &
            (pgl(igroup)%ndim.ne.3*(pgl(igroup)%next-2)-4 .and. &
            pgl(igroup)%ndim.ne.3*(pgl(igroup)%next-2)-2)) then
          write (*,*) 'Invalid group dimensions in amplitude library',igroup,&
               pgl(igroup)%next,pgl(igroup)%ndim
          stop 1
       endif
       call reserve_library_array(int(pgl(igroup)%next,kind=8),4_8,&
            'helicity workspace',.false.)
       allocate(pgl(igroup)%hel(pgl(igroup)%next),stat=ios,errmsg=allocation_message)
       call require_allocation('helicity workspace')
       pgl(igroup)%hel=pgl(igroup)%spin(1,1:pgl(igroup)%next)
       read(14,iostat=ios) dim1
       call require_read('colour-factor extent')
       call require_extent(dim1,'colour-factor extent')
       if (dim1.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent colour-factor extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,'colour factors',.true.)
       allocate(pgl(igroup)%col_fac(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('colour factors')
       read(14,iostat=ios) pgl(igroup)%col_fac
       call require_read('colour factors')
       read(14,iostat=ios) dim1
       call require_read('matrix-element-square extent')
       call require_extent(dim1,'matrix-element-square extent')
       call reserve_library_array(int(dim1,kind=8),8_8,&
            'matrix-element-square workspace',.false.)
       allocate(pgl(igroup)%amp2(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('matrix-element-square workspace')
       pgl(igroup)%amp2=0d0
       read(14,iostat=ios) dim1
       call require_read('helicity-square extent')
       call require_extent(dim1,'helicity-square extent')
       call reserve_library_array(int(dim1,kind=8),8_8,&
            'helicity-square workspace',.false.)
       allocate(pgl(igroup)%amp2_hel(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('helicity-square workspace')
       pgl(igroup)%amp2_hel=0d0
       read(14,iostat=ios) dim1
       call require_read('passed-count extent')
       call require_extent(dim1,'passed-count extent')
       if (dim1.ne.size(pgl(igroup)%amps)) then
          write (*,*) 'Inconsistent passed-count extent in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8),4_8,'passed counts',.true.)
       allocate(pgl(igroup)%passed(dim1),stat=ios,errmsg=allocation_message)
       call require_allocation('passed counts')
       read(14,iostat=ios) pgl(igroup)%passed
       call require_read('passed counts')
       if (any(pgl(igroup)%passed.lt.0)) then
          write (*,*) 'Invalid passed count in amplitude library',igroup
          stop 1
       endif
       read(14,iostat=ios) dim1,dim2
       call require_read('colour-order shape')
       call require_shape2(dim1,dim2,'colour-order shape')
       if (dim1.ne.size(pgl(igroup)%processes,1) .or. &
            dim2.ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent colour-order shape in amplitude library',igroup
          stop 1
       endif
       call reserve_library_array(int(dim1,kind=8)*int(dim2,kind=8),4_8,&
            'colour orders',.true.)
       allocate(pgl(igroup)%color_orders(dim1,dim2),stat=ios,errmsg=allocation_message)
       call require_allocation('colour orders')
       read(14,iostat=ios) pgl(igroup)%color_orders
       call require_read('colour orders')
       if (any(shape(pgl(igroup)%color_orders).ne.shape(pgl(igroup)%processes))) then
          write (*,*) 'Inconsistent colour-order shape in amplitude library',igroup
          stop 1
       endif
       if (any(pgl(igroup)%color_orders.lt.1) .or. &
            any(pgl(igroup)%color_orders.gt.pgl(igroup)%next)) then
          write (*,*) 'Out-of-range colour order in amplitude library',igroup
          stop 1
       endif
       if (pgl(igroup)%next.ne.pgl_unique%next) then
          write (*,*) 'Inconsistent external multiplicity across amplitude-library groups',igroup
          stop 1
       endif
       do label=1,pgl(igroup)%next
          if (count(pgl(igroup)%phase_space_orders.eq.label).ne.1) then
             write (*,*) 'Invalid phase-space order in amplitude library',igroup
             stop 1
          endif
       enddo
       if (size(pgl(igroup)%phase_space_permutations,2).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%color_orders,2).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%nhel).ne.size(pgl(igroup)%amps) .or. &
            size(pgl(igroup)%passed).ne.size(pgl(igroup)%amps) .or. &
            size(pgl(igroup)%col_fac).ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent group metadata extents in amplitude library',igroup
          stop 1
       endif
       if (size(pgl(igroup)%val_procs,2).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%idenCOandMAPfactor,2).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%iden_processes,1).ne.pgl(igroup)%next .or. &
            size(pgl(igroup)%iden_processes,3).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%val_procs,1).lt.maxval(pgl(igroup)%iden_iproc) .or. &
            size(pgl(igroup)%idenCOandMAPfactor,1).lt.maxval(pgl(igroup)%iden_iproc) .or. &
            size(pgl(igroup)%iden_processes,2).lt.maxval(pgl(igroup)%iden_iproc)) then
          write (*,*) 'Inconsistent identical-process metadata in amplitude library',igroup
          stop 1
       endif
       if (size(pgl(igroup)%iden).ne.pgl(igroup)%nproc .or. any(pgl(igroup)%iden.lt.1_8) .or. &
            any(pgl(igroup)%col_fac.lt.1)) then
          write (*,*) 'Invalid process normalisation metadata in amplitude library',igroup
          stop 1
       endif
       if (size(pgl(igroup)%amp2_hel).ne.maxval(pgl(igroup)%nhel) .or. &
            size(pgl(igroup)%hel_fac,1).ne.size(pgl(igroup)%amp2_hel)) then
          write (*,*) 'Inconsistent helicity workspace metadata in amplitude library',igroup
          stop 1
       endif
       if (size(pgl(igroup)%multichan%number_of_channels).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%multichan%channels,2).ne.pgl(igroup)%nproc .or. &
            size(pgl(igroup)%multichan%channel_permutations,1).ne.pgl(igroup)%next .or. &
            size(pgl(igroup)%multichan%channel_permutations,2).ne.&
            pgl(igroup)%multichan%max_channels .or. &
            size(pgl(igroup)%multichan%channel_permutations,3).ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent multichannel metadata in amplitude library',igroup
          stop 1
       endif
       if (size(pgl(igroup)%multichan%unique_channel_list).ne.&
            pgl(igroup)%multichan%n_unique_channels .or. &
            size(pgl(igroup)%multichan%unique_channelgroup_list,1).ne.&
            pgl(igroup)%multichan%max_channels+1 .or. &
            size(pgl(igroup)%multichan%unique_channelgroup_list,2).ne.&
            pgl(igroup)%multichan%n_unique_channelgroups .or. &
            size(pgl(igroup)%multichan%map_proc_to_channelgroup).ne.pgl(igroup)%nproc) then
          write (*,*) 'Inconsistent compressed multichannel metadata in amplitude library',igroup
          stop 1
       endif
       if (any(pgl(igroup)%multichan%unique_channel_list.lt.1) .or. &
            any(pgl(igroup)%multichan%unique_channel_list.gt.ngroups) .or. &
            any(pgl(igroup)%multichan%map_proc_to_channelgroup.lt.1) .or. &
            any(pgl(igroup)%multichan%map_proc_to_channelgroup.gt.&
            pgl(igroup)%multichan%n_unique_channelgroups)) then
          write (*,*) 'Out-of-range compressed multichannel metadata in amplitude library',igroup
          stop 1
       endif
       allocate(seen_channels(ngroups),stat=ios,errmsg=allocation_message)
       call require_allocation('unique-channel validation workspace')
       seen_channels=.false.
       do iamp=1,pgl(igroup)%multichan%n_unique_channels
          label=pgl(igroup)%multichan%unique_channel_list(iamp)
          if (seen_channels(label)) then
             write (*,*) 'Duplicate unique-channel entry in amplitude library',igroup,iamp
             stop 1
          endif
          seen_channels(label)=.true.
       enddo
       deallocate(seen_channels)
       if (any(pgl(igroup)%multichan%number_of_channels.lt.1) .or. &
            any(pgl(igroup)%multichan%number_of_channels.gt.&
            pgl(igroup)%multichan%max_channels)) then
          write (*,*) 'Invalid channel multiplicity in amplitude library',igroup
          stop 1
       endif
       if (any(pgl(igroup)%phase_space_permutations(1,:).ne.1) .or. &
            any(pgl(igroup)%phase_space_permutations(2,:).ne.2) .or. &
            any(pgl(igroup)%multichan%channel_permutations(1,:,:).ne.1) .or. &
            any(pgl(igroup)%multichan%channel_permutations(2,:,:).ne.2)) then
          write (*,*) 'Amplitude library permutes an incoming phase-space leg',igroup
          stop 1
       endif
       do label=1,pgl(igroup)%next
          if (pgl(igroup)%spin(0,label).ne.&
               phys_model%get_spin(pgl(igroup)%processes(label,1))) then
             write (*,*) 'Spin table disagrees with the process model in amplitude library',&
                  igroup,label
             stop 1
          endif
          select case (pgl(igroup)%spin(0,label))
          case (1)
             if (pgl(igroup)%spin(1,label).ne.0) then
                write (*,*) 'Invalid scalar helicity table in amplitude library',igroup,label
                stop 1
             endif
          case (2)
             if (pgl(igroup)%spin(1,label).ne.-1 .or. &
                  pgl(igroup)%spin(2,label).ne.1) then
                write (*,*) 'Invalid massless helicity table in amplitude library',igroup,label
                stop 1
             endif
          case (3)
             if (any(pgl(igroup)%spin(1:3,label).ne.[-1,0,1])) then
                write (*,*) 'Invalid massive-vector helicity table in amplitude library',igroup,label
                stop 1
             endif
          case default
             write (*,*) 'Invalid external-particle spin in amplitude library',igroup,label
             stop 1
          end select
       enddo
       do iproc=1,pgl(igroup)%nproc
          ncoloured=0
          nsinglet=0
          do label=1,pgl(igroup)%next
             if (phys_model%get_spin(pgl(igroup)%processes(label,iproc)).ne.&
                  pgl(igroup)%spin(0,label)) then
                write (*,*) 'Incompatible process spins in amplitude library',igroup,iproc,label
                stop 1
             endif
             if (phys_model%is_singlet(pgl(igroup)%processes(label,iproc))) then
                nsinglet=nsinglet+1
             else
                ncoloured=ncoloured+1
             endif
          enddo
          do label=1,pgl(igroup)%next
             if (count(pgl(igroup)%phase_space_permutations(:,iproc).eq.label).ne.1 .or. &
                  count(pgl(igroup)%color_orders(:,iproc).eq.label).ne.1) then
                write (*,*) 'Invalid process permutation in amplitude library',igroup,iproc
                stop 1
             endif
          enddo
          if (ncoloured.gt.0 .and. nsinglet.gt.0) then
             label=pgl(igroup)%color_orders(pgl(igroup)%next,iproc)
             if (phys_model%is_singlet(pgl(igroup)%processes(label,iproc))) then
                write (*,*) 'Colour order closes on a singlet in amplitude library',igroup,iproc
                stop 1
             endif
          endif
          if (any(pgl(igroup)%multichan%channels(&
               1:pgl(igroup)%multichan%number_of_channels(iproc),iproc).lt.1) .or. &
               any(pgl(igroup)%multichan%channels(&
               1:pgl(igroup)%multichan%number_of_channels(iproc),iproc).gt.ngroups)) then
             write (*,*) 'Out-of-range multichannel group in amplitude library',igroup,iproc
             stop 1
          endif
          do iamp=1,pgl(igroup)%multichan%number_of_channels(iproc)
             do label=1,pgl(igroup)%next
                if (count(pgl(igroup)%multichan%channel_permutations(:,iamp,iproc).eq.label).ne.1) then
                   write (*,*) 'Invalid multichannel permutation in amplitude library',&
                        igroup,iproc,iamp
                   stop 1
                endif
             enddo
          enddo
          iamp=pgl(igroup)%multichan%map_proc_to_channelgroup(iproc)
          if (pgl(igroup)%multichan%unique_channelgroup_list(0,iamp).ne.&
               pgl(igroup)%multichan%number_of_channels(iproc)) then
             write (*,*) 'Inconsistent compressed channel-group size in amplitude library',igroup,iproc
             stop 1
          endif
          do k=1,pgl(igroup)%multichan%number_of_channels(iproc)
             label=pgl(igroup)%multichan%unique_channelgroup_list(k,iamp)
             if (label.lt.1 .or. label.gt.pgl(igroup)%multichan%n_unique_channels) then
                write (*,*) 'Out-of-range compressed channel entry in amplitude library',igroup,iproc
                stop 1
             endif
             if (pgl(igroup)%multichan%unique_channel_list(label).ne.&
                  pgl(igroup)%multichan%channels(k,iproc)) then
                write (*,*) 'Compressed channel group does not reproduce its process list',igroup,iproc
                stop 1
             endif
          enddo
       enddo
       do iamp=1,pgl(igroup)%multichan%n_unique_channelgroups
          dim1=pgl(igroup)%multichan%unique_channelgroup_list(0,iamp)
          if (dim1.lt.1 .or. dim1.gt.pgl(igroup)%multichan%max_channels) then
             write (*,*) 'Invalid compressed channel-group size in amplitude library',igroup,iamp
             stop 1
          endif
          if (any(pgl(igroup)%multichan%unique_channelgroup_list(1:dim1,iamp).lt.1) .or. &
               any(pgl(igroup)%multichan%unique_channelgroup_list(1:dim1,iamp).gt.&
               pgl(igroup)%multichan%n_unique_channels)) then
             write (*,*) 'Invalid compressed channel-group entries in amplitude library',igroup,iamp
             stop 1
          endif
       enddo
       allocate(channel_marks(pgl(igroup)%multichan%n_unique_channels),&
            referenced_unique_channels(pgl(igroup)%multichan%n_unique_channels),&
            stat=ios,errmsg=allocation_message)
       call require_allocation('compressed-channel validation workspace')
       if (.not.allocated(referenced_unique_channels)) then
          write (*,*) 'Compressed-channel reference workspace was not allocated'
          close(14)
          stop 1
       endif
       channel_marks=0
       referenced_unique_channels=.false.
       do iamp=1,pgl(igroup)%multichan%n_unique_channelgroups
          dim1=pgl(igroup)%multichan%unique_channelgroup_list(0,iamp)
          do k=1,dim1
             label=pgl(igroup)%multichan%unique_channelgroup_list(k,iamp)
             if (channel_marks(label).eq.iamp) then
                write (*,*) 'Duplicate channel inside compressed channel group',&
                     igroup,iamp,label
                stop 1
             endif
             channel_marks(label)=iamp
             referenced_unique_channels(label)=.true.
          enddo
       enddo
       if (.not.all(referenced_unique_channels)) then
          write (*,*) 'Amplitude library contains an unused unique-channel entry',igroup
          stop 1
       endif
       deallocate(channel_marks,referenced_unique_channels)
       if (keep_processes_separate) then
          if (size(pgl(igroup)%amps).ne.pgl(igroup)%nproc .or. &
               size(pgl(igroup)%hel_fac,2).ne.pgl(igroup)%nproc .or. &
               size(pgl(igroup)%amp2).ne.1) then
             write (*,*) 'Amplitude library is incompatible with keep_processes_separate=true',igroup
             stop 1
          endif
       else
          if (size(pgl(igroup)%amps).ne.1 .or. size(pgl(igroup)%hel_fac,2).ne.1 .or. &
               size(pgl(igroup)%amp2).ne.pgl(igroup)%nproc) then
             write (*,*) 'Amplitude library is incompatible with keep_processes_separate=false',igroup
             stop 1
          endif
       endif
       do iamp=1,size(pgl(igroup)%amps)
          stored_nprocs=size(pgl(igroup)%amps(iamp)%iproc_start)-1
          if (stored_nprocs.lt.1 .or. &
               size(pgl(igroup)%amps(iamp)%n_sing).ne.stored_nprocs .or. &
               size(pgl(igroup)%amps(iamp)%n_qqbar).ne.stored_nprocs) then
             write (*,*) 'Inconsistent amplitude subprocess metadata in library',igroup,iamp
             stop 1
          endif
          if (pgl(igroup)%amps(iamp)%iproc_start(1).ne.1 .or. &
               pgl(igroup)%amps(iamp)%iproc_start(stored_nprocs+1).ne.&
               pgl(igroup)%amps(iamp)%n_amps+1 .or. &
               any(pgl(igroup)%amps(iamp)%iproc_start(2:stored_nprocs+1).lt.&
               pgl(igroup)%amps(iamp)%iproc_start(1:stored_nprocs))) then
             write (*,*) 'Invalid amplitude subprocess offsets in library',igroup,iamp
             stop 1
          endif
          if (any(pgl(igroup)%amps(iamp)%n_sing.lt.0) .or. &
               any(pgl(igroup)%amps(iamp)%n_sing.gt.pgl(igroup)%next) .or. &
               any(pgl(igroup)%amps(iamp)%n_qqbar.lt.0) .or. &
               any(pgl(igroup)%amps(iamp)%n_qqbar.gt.3)) then
             write (*,*) 'Invalid amplitude particle-count metadata in library',igroup,iamp
             stop 1
          endif
          if (size(pgl(igroup)%amps(iamp)%spins,1).ne.pgl(igroup)%next .or. &
               size(pgl(igroup)%amps(iamp)%spins,3).ne.pgl(igroup)%amps(iamp)%n_amps .or. &
               pgl(igroup)%nhel(iamp).ne.pgl(igroup)%amps(iamp)%n_amps) then
             write (*,*) 'Inconsistent amplitude helicity metadata in library',igroup,iamp
             stop 1
          endif
          do k=1,size(pgl(igroup)%amps(iamp)%spins,2)
             do label=1,pgl(igroup)%next
                select case (pgl(igroup)%spin(0,label))
                case (1)
                   if (any(pgl(igroup)%amps(iamp)%spins(label,k,:).ne.0)) then
                      write (*,*) 'Invalid scalar amplitude helicity in library',&
                           igroup,iamp,label
                      stop 1
                   endif
                case (2)
                   if (any(pgl(igroup)%amps(iamp)%spins(label,k,:).ne.-1 .and. &
                        pgl(igroup)%amps(iamp)%spins(label,k,:).ne.1)) then
                      write (*,*) 'Invalid massless amplitude helicity in library',&
                           igroup,iamp,label
                      stop 1
                   endif
                case (3)
                   if (any(pgl(igroup)%amps(iamp)%spins(label,k,:).lt.-1) .or. &
                        any(pgl(igroup)%amps(iamp)%spins(label,k,:).gt.1)) then
                      write (*,*) 'Invalid vector amplitude helicity in library',&
                           igroup,iamp,label
                      stop 1
                   endif
                end select
             enddo
          enddo
          pgl(igroup)%amps(iamp)%nprocs=stored_nprocs
       enddo
    enddo
    do igroup=1,ngroups
       do iproc=1,pgl(igroup)%nproc
          self_count=0
          do k=1,pgl(igroup)%multichan%number_of_channels(iproc)
             partner=pgl(igroup)%multichan%channels(k,iproc)
             if (pgl(partner)%next.ne.pgl(igroup)%next) then
                write (*,*) 'Incompatible partner multiplicity in amplitude library',igroup,partner
                stop 1
             endif
             if (partner.eq.igroup .and. all(&
                  pgl(igroup)%multichan%channel_permutations(:,k,iproc).eq.&
                  pgl(igroup)%phase_space_permutations(:,iproc))) self_count=self_count+1
          enddo
          if (self_count.ne.1) then
             write (*,*) 'Amplitude library does not contain exactly one self density',&
                  igroup,iproc,self_count
             stop 1
          endif
       enddo
    enddo
    read(14,iostat=ios) trailing_byte
    if (ios.ne.iostat_end) then
       write (*,*) 'Amplitude library contains trailing or unreadable metadata',ios
       stop 1
    endif
    close(14)
    call apply_final_state_widths_from_loaded_groups()
    current_model_signature=phys_model%model_signature()
    if (.not.all(ieee_is_finite(current_model_signature))) then
       write (*,*) 'Current model has a non-finite amplitude-library signature'
       stop 1
    endif
    comparison_scale=max(1d0,abs(current_model_signature),&
         abs(stored_model_signature))
    comparison_residual=abs(current_model_signature/comparison_scale-&
         stored_model_signature/comparison_scale)
    if (any(comparison_residual.gt.1d-13)) then
       write (*,*) 'Amplitude library was created with incompatible model parameters.'
       write (*,*) 'Recreate it with --library=create and the current input card.'
       write (*,*) 'Stored model signature:',stored_model_signature
       write (*,*) 'Current model signature:',current_model_signature
       stop 1
    endif
    ! Phase-space dimensionality and PDF flavour masks are run settings, not
    ! properties of the compiled matrix elements.  Reconstruct them so that a
    ! library can be reused when include_pdf changes.
    pgl_unique%ndim=3*(pgl_unique%next-2)-4
    if (include_pdf) pgl_unique%ndim=pgl_unique%ndim+2
    do igroup=1,ngroups
       pgl(igroup)%ndim=3*(pgl(igroup)%next-2)-4
       if (include_pdf) then
          pgl(igroup)%ndim=pgl(igroup)%ndim+2
          call set_ipdgs_for_PDF(pgl(igroup))
       else
          pgl(igroup)%ipdgs=.false.
       endif
    enddo
  contains
    subroutine require_allocation(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Could not allocate amplitude-library ',trim(label),': ',&
              trim(allocation_message)
         close(14)
         stop 1
      endif
    end subroutine require_allocation

    subroutine reserve_library_array(elements,element_bytes,label,serialized)
      implicit none
      integer(kind=8),intent(in) :: elements,element_bytes
      character(len=*),intent(in) :: label
      logical,intent(in) :: serialized
      integer(kind=8) :: bytes,stream_position,stream_size,remaining

      if (elements.lt.1_8 .or. element_bytes.lt.1_8) then
         write (*,*) 'Invalid amplitude-library allocation for ',trim(label),&
              elements,element_bytes
         close(14)
         stop 1
      endif
      if (elements.gt.huge(bytes)/element_bytes) then
         write (*,*) 'Amplitude-library byte count overflows for ',trim(label),&
              elements,element_bytes
         close(14)
         stop 1
      endif
      bytes=elements*element_bytes
      if (library_workspace_bytes.gt.max_library_workspace_bytes-bytes) then
         write (*,*) 'Amplitude-library metadata exceeds the supported workspace for ',&
              trim(label),library_workspace_bytes,bytes,max_library_workspace_bytes
         close(14)
         stop 1
      endif
      if (serialized) then
         inquire(unit=14,pos=stream_position,size=stream_size,iostat=ios)
         if (ios.ne.0) then
            write (*,*) 'Cannot inspect amplitude-library stream while reading ',trim(label),ios
            close(14)
            stop 1
         endif
         remaining=max(0_8,stream_size-stream_position+1_8)
         if (bytes.gt.remaining) then
            write (*,*) 'Amplitude-library payload is truncated before ',trim(label),&
                 bytes,remaining
            close(14)
            stop 1
         endif
      endif
      library_workspace_bytes=library_workspace_bytes+bytes
    end subroutine reserve_library_array

    subroutine require_read(label)
      implicit none
      character(len=*),intent(in) :: label
      if (ios.ne.0) then
         write (*,*) 'Truncated/corrupt amplitude library while reading ',trim(label),ios
         close(14)
         stop 1
      endif
    end subroutine require_read

    subroutine require_extent(extent,label)
      implicit none
      integer,intent(in) :: extent
      character(len=*),intent(in) :: label
      if (extent.lt.1 .or. extent.gt.max_library_extent) then
         write (*,*) 'Invalid amplitude-library extent for ',trim(label),extent
         close(14)
         stop 1
      endif
    end subroutine require_extent

    subroutine require_shape2(first,second,label)
      implicit none
      integer,intent(in) :: first,second
      character(len=*),intent(in) :: label
      integer(kind=8) :: elements
      call require_extent(first,label)
      call require_extent(second,label)
      elements=int(first,kind=8)*int(second,kind=8)
      if (elements.gt.int(max_library_extent,kind=8)) then
         write (*,*) 'Amplitude-library array is too large for ',trim(label),first,second
         close(14)
         stop 1
      endif
    end subroutine require_shape2

    subroutine require_shape3(first,second,third,label)
      implicit none
      integer,intent(in) :: first,second,third
      character(len=*),intent(in) :: label
      integer(kind=8) :: elements
      call require_extent(first,label)
      call require_extent(second,label)
      call require_extent(third,label)
      elements=int(first,kind=8)*int(second,kind=8)
      if (elements.gt.int(max_library_extent,kind=8)/int(third,kind=8)) then
         write (*,*) 'Amplitude-library array is too large for ',trim(label),first,second,third
         close(14)
         stop 1
      endif
    end subroutine require_shape3
  end subroutine read_amplitude_lib

  subroutine test_lib
    use amp_lib
    use FeynmanRules, only: feynman_numerical_status_ok
    implicit none
    integer :: igroup,iint,i,ios
    integer(kind=1) :: trailing_byte
    real(kind=8),dimension(:,:),allocatable :: p
    real(kind=8) :: comparison_scale
    complex(kind=8),dimension(:),allocatable :: amps_save,amps
    character(len=170) :: line,tmp
    character(len=256) :: allocation_message
    logical :: mismatch
    do igroup=1,ngroups
       do iint=1,size(pgl(igroup)%amps)
          allocate(p(0:3,pgl(igroup)%next),&
               amps(size(pgl(igroup)%amps(iint)%amps)),&
               amps_save(size(pgl(igroup)%amps(iint)%amps)),&
               stat=ios,errmsg=allocation_message)
          if (ios.ne.0) then
             write (*,*) 'Could not allocate amplitude-library test workspace: ',&
                  trim(allocation_message)
             stop 1
          endif
          write(tmp,*) igroup
          write(line,*) iint
          line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
          open(file=line,unit=14,form='unformatted',access='stream',status='old',&
               action='read',iostat=ios)
          if (ios.ne.0) then
             write (*,*) 'Could not open amplitude-library reference data',trim(line)
             stop 1
          endif
          read(14,iostat=ios) p
          if (ios.ne.0) then
             write (*,*) 'Could not read amplitude-library reference momenta',igroup,iint
             stop 1
          endif
          read(14,iostat=ios) amps_save
          if (ios.ne.0) then
             write (*,*) 'Could not read amplitude-library reference amplitudes',igroup,iint
             stop 1
          endif
          read(14,iostat=ios) trailing_byte
          if (ios.ne.iostat_end) then
             write (*,*) 'Amplitude-library reference data contains trailing bytes',igroup,iint,ios
             stop 1
          endif
          close(14)
          if (.not.all(ieee_is_finite(p)) .or. &
               .not.all(ieee_is_finite(real(amps_save,kind=8))) .or. &
               .not.all(ieee_is_finite(aimag(amps_save)))) then
             write (*,*) 'Non-finite amplitude-library reference data',igroup,iint
             stop 1
          endif
          call evaluate_amp(igroup,iint,p,amps)
          mismatch=.not.feynman_numerical_status_ok()
          if (.not.all(ieee_is_finite(real(amps,kind=8))) .or. &
               .not.all(ieee_is_finite(aimag(amps)))) mismatch=.true.
          do i=1,size(amps)
             if (mismatch) exit
             comparison_scale=max(abs(real(amps_save(i),kind=8)),&
                  abs(aimag(amps_save(i))),abs(real(amps(i),kind=8)),&
                  abs(aimag(amps(i))))
             if (comparison_scale.eq.0d0) cycle
             if (abs(amps_save(i)/comparison_scale-amps(i)/comparison_scale).gt.1d-8) &
                  mismatch=.true.
          enddo
          if (mismatch) then
             write (*,*) 'Process library not compatible with saved amplitudes',igroup,iint
             do i=1,size(pgl(igroup)%amps(iint)%amps)
                if (.not.ieee_is_finite(real(amps_save(i),kind=8)) .or. &
                     .not.ieee_is_finite(aimag(amps_save(i))) .or. &
                     .not.ieee_is_finite(real(amps(i),kind=8)) .or. &
                     .not.ieee_is_finite(aimag(amps(i)))) then
                   write (*,*) i,amps_save(i)
                   write (*,*) i,amps(i)
                   cycle
                endif
                comparison_scale=max(abs(real(amps_save(i),kind=8)),&
                     abs(aimag(amps_save(i))),abs(real(amps(i),kind=8)),&
                     abs(aimag(amps(i))))
                if (comparison_scale.eq.0d0) cycle
                if (abs(amps_save(i)/comparison_scale-amps(i)/comparison_scale).gt.1d-8) then
                   write (*,*) i,amps_save(i)
                   write (*,*) i,amps(i)
                endif
             enddo
             stop 1
          endif
          deallocate(p)
          deallocate(amps)
          deallocate(amps_save)
       enddo
    enddo
  end subroutine test_lib
end module amplitude_library
