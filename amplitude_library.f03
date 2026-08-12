module amplitude_library
  use handling_processes
  use read_process_file
  use coupling_orders, only: coupling_order_selection,set_coupling_order_selection,&
       resolve_automatic_coupling_order,coupling_order_mode_automatic,&
       coupling_order_unbounded
  use pdf_wrap, only: set_ipdgs_for_PDF
contains
  subroutine create_amplitude_lib()
    implicit none
    character(len=170) :: tmp,line,filename
    integer :: igroup,j,iamp
    real(kind=8),dimension(9) :: stored_model_signature
    filename='Library/amplib.f03'
    open(unit=14,file=filename,status='unknown')
    write(14,*) 'module amp_lib'
    do igroup=1,ngroups
       do j=1,size(pgl(igroup)%amps)
          write(tmp,*) igroup
          line=trim(adjustl(tmp))//'_'
          write(tmp,*) j
          line=trim(adjustl(line))//trim(adjustl(tmp))
          write(14,*) 'use amp'//trim(adjustl(line))//'_lib'
       enddo
    enddo
    write(14,*) 'implicit none'
    write(14,*) 'private'
    write(14,*) 'public :: evaluate_amp,evaluate_amp_by_order'
    write(14,*) 'contains'
    write(14,*) 'subroutine evaluate_amp(ichan,iint,p,amps)'
    write(14,*) 'implicit none'
    write(14,*) 'integer,intent(in) :: ichan,iint'
    write(tmp,*) maxval(pgl(:)%next)
    write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(14,*) 'complex(kind=8),intent(out) :: amps(*)'
    do igroup=1,ngroups
       if (igroup.eq.1) then
          write(14,*) 'if (ichan.eq.1) then'
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                write(14,*) 'if (iint.eq.1) then'
                write(14,*) 'call evaluate_amp1_1(p,amps)'
             else
                write(tmp,*) j
                write(14,*) 'elseif (iint.eq.'//trim(adjustl(tmp))//') then'
                write(14,*) 'call evaluate_amp1_'//trim(adjustl(tmp))//'(p,amps)'
             endif
          enddo
          write(14,*) 'endif'
       else
          write(tmp,*) igroup
          write(14,*) 'elseif (ichan.eq.'//trim(adjustl(tmp))//') then'
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                write(14,*) 'if (iint.eq.1) then'
                write(14,*) 'call evaluate_amp'//trim(adjustl(tmp))//'_1(p,amps)'
             else
                write(line,*) j
                write(14,*) 'elseif (iint.eq.'//trim(adjustl(line))//') then'
                write(14,*) 'call evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'(p,amps)'
             endif
          enddo
          write(14,*) 'endif'
       endif
    enddo
    write(14,*) 'endif'
    write(14,*) 'end subroutine evaluate_amp'
    write(14,*) 'subroutine evaluate_amp_by_order(ichan,iint,p,amps_by_order)'
    write(14,*) 'implicit none'
    write(14,*) 'integer,intent(in) :: ichan,iint'
    write(tmp,*) maxval(pgl(:)%next)
    write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(14,*) 'complex(kind=8),intent(out) :: amps_by_order(*)'
    do igroup=1,ngroups
       if (igroup.eq.1) then
          write(14,*) 'if (ichan.eq.1) then'
       else
          write(tmp,*) igroup
          write(14,*) 'elseif (ichan.eq.'//trim(adjustl(tmp))//') then'
       endif
       do j=1,size(pgl(igroup)%amps)
          write(tmp,*) igroup
          write(line,*) j
          if (j.eq.1) then
             write(14,*) 'if (iint.eq.1) then'
          else
             write(14,*) 'elseif (iint.eq.'//trim(adjustl(line))//') then'
          endif
          write(14,*) 'call evaluate_amp'//trim(adjustl(tmp))//'_'//&
               trim(adjustl(line))//'_by_order(p,amps_by_order)'
       enddo
       write(14,*) 'endif'
    enddo
    write(14,*) 'endif'
    write(14,*) 'end subroutine evaluate_amp_by_order'
    write(14,*) 'end module amp_lib'
    close(14)
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
    write(14) 5 ! binary amplitude-library format version
    stored_model_signature=phys_model%model_signature()
    write(14) stored_model_signature
    write(14) coupling_order_selection%mode,coupling_order_selection%resolved_as2,&
         coupling_order_selection%as_min2,coupling_order_selection%as_max2,&
         coupling_order_selection%aew_min2,coupling_order_selection%aew_max2
    write(14) pgl_unique%next,pgl_unique%nproc
    write(14) unique_map
    write(14) unique_map_value
    write(14) pgl_unique%processes
    write(14) ngroups
    do igroup=1,ngroups
       ! amplitudes
       write(14) size(pgl(igroup)%amps)
       do iamp=1,size(pgl(igroup)%amps)
          write(14) pgl(igroup)%amps(iamp)%n_amps
          write(14) pgl(igroup)%amps(iamp)%n_sectors
          write(14) shape(pgl(igroup)%amps(iamp)%sector_powers),&
               pgl(igroup)%amps(iamp)%sector_powers
          write(14) shape(pgl(igroup)%amps(iamp)%sector_present),&
               pgl(igroup)%amps(iamp)%sector_present
          write(14) shape(pgl(igroup)%amps(iamp)%sector_sign),&
               pgl(igroup)%amps(iamp)%sector_sign
          write(14) size(pgl(igroup)%amps(iamp)%iproc_start)
          write(14) pgl(igroup)%amps(iamp)%iproc_start
          write(14) size(pgl(igroup)%amps(iamp)%n_sing),pgl(igroup)%amps(iamp)%n_sing
          write(14) size(pgl(igroup)%amps(iamp)%n_qqbar),pgl(igroup)%amps(iamp)%n_qqbar
          write(14) shape(pgl(igroup)%amps(iamp)%spins),pgl(igroup)%amps(iamp)%spins
          write(14) size(pgl(igroup)%amps(iamp)%amps)
       enddo
       ! multichannel
       write(14) pgl(igroup)%multichan%n_unique_channels,pgl(igroup)%multichan%n_unique_channelgroups
       write(14) shape(pgl(igroup)%multichan%unique_channelgroup_list),&
            pgl(igroup)%multichan%unique_channelgroup_list
       write(14) size(pgl(igroup)%multichan%unique_channel_list),pgl(igroup)%multichan%unique_channel_list
       write(14) size(pgl(igroup)%multichan%map_proc_to_channelgroup),&
            pgl(igroup)%multichan%map_proc_to_channelgroup
       write(14) size(pgl(igroup)%multichan%number_of_channels),pgl(igroup)%multichan%number_of_channels
       write(14) shape(pgl(igroup)%multichan%channels),pgl(igroup)%multichan%channels
       write(14) shape(pgl(igroup)%multichan%channel_permutations),&
            pgl(igroup)%multichan%channel_permutations
       ! rest
       write(14) shape(pgl(igroup)%processes),pgl(igroup)%processes
       write(14) shape(pgl(igroup)%phase_space_permutations),pgl(igroup)%phase_space_permutations
       write(14) size(pgl(igroup)%iden_iproc),pgl(igroup)%iden_iproc
       write(14) size(pgl(igroup)%phase_space_orders),pgl(igroup)%phase_space_orders
       write(14) size(pgl(igroup)%nhel),pgl(igroup)%nhel
       write(14) pgl(igroup)%nproc
       write(14) shape(pgl(igroup)%val_procs),pgl(igroup)%val_procs
       write(14) shape(pgl(igroup)%val_procs_abs)
       write(14) shape(pgl(igroup)%alias_factors)
       write(14) shape(pgl(igroup)%idenCOandMAPfactor),pgl(igroup)%idenCOandMAPfactor
       write(14) shape(pgl(igroup)%iden_processes),pgl(igroup)%iden_processes
       write(14) shape(pgl(igroup)%spin),pgl(igroup)%spin
       write(14) shape(pgl(igroup)%hel_fac),pgl(igroup)%hel_fac
       write(14) size(pgl(igroup)%iden),pgl(igroup)%iden
       write(14) pgl(igroup)%ipdgs
       write(14) pgl(igroup)%next,pgl(igroup)%ndim
       write(14) size(pgl(igroup)%col_fac),pgl(igroup)%col_fac
       write(14) size(pgl(igroup)%amp2)
       write(14) size(pgl(igroup)%amp2_abs)
       write(14) size(pgl(igroup)%amp2_hel)
       write(14) size(pgl(igroup)%amp2_hel_abs)
       write(14) size(pgl(igroup)%passed),pgl(igroup)%passed
       write(14) shape(pgl(igroup)%color_orders),pgl(igroup)%color_orders
    enddo
    close(14)
  end subroutine create_amplitude_lib

  subroutine read_amplitude_lib()
    implicit none
    character(len=170) :: filename
    integer :: dim1,dim2,dim3,iamp,igroup,library_version
    integer :: selection_mode,resolved_as2,as_min2,as_max2,aew_min2,aew_max2
    logical :: valid_selection
    real(kind=8),dimension(9) :: stored_model_signature,current_model_signature
    real(kind=8),dimension(9) :: tolerance
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='old')
    read(14) library_version
    if (library_version.ne.5) then
       write (*,*) 'Amplitude library has an incompatible binary format; recreate it'
       stop 1
    endif
    read(14) stored_model_signature
    read(14) selection_mode,resolved_as2,as_min2,as_max2,aew_min2,aew_max2
    call set_coupling_order_selection(selection_mode,as_min2,as_max2,aew_min2,aew_max2,&
         valid_selection)
    if (.not.valid_selection) then
       write (*,*) 'Amplitude library contains an invalid coupling-order selection'
       stop 1
    endif
    if (selection_mode.eq.coupling_order_mode_automatic) then
       call resolve_automatic_coupling_order(resolved_as2,valid_selection)
       if (.not.valid_selection) then
          write (*,*) 'Amplitude library contains an invalid resolved automatic coupling order'
          stop 1
       endif
    elseif (resolved_as2.ne.coupling_order_unbounded) then
       write (*,*) 'Amplitude library contains inconsistent explicit coupling-order metadata'
       stop 1
    endif
    allocate(pgl_unique)
    read(14) pgl_unique%next,pgl_unique%nproc
    next=pgl_unique%next
    allocate(unique_map(pgl_unique%nproc))
    read(14) unique_map
    allocate(unique_map_value(pgl_unique%nproc))
    read(14) unique_map_value
    allocate(pgl_unique%processes(pgl_unique%next,pgl_unique%nproc))
    read(14) pgl_unique%processes
    read(14)ngroups
    allocate(pgl(ngroups))
    do igroup=1,ngroups
       ! amplitudes
       read(14) dim1
       allocate(pgl(igroup)%amps(dim1))
       do iamp=1,dim1
          read(14) pgl(igroup)%amps(iamp)%n_amps
          read(14) pgl(igroup)%amps(iamp)%n_sectors
          read(14) dim1,dim2
          allocate(pgl(igroup)%amps(iamp)%sector_powers(dim1,dim2))
          read(14) pgl(igroup)%amps(iamp)%sector_powers
          read(14) dim1,dim2
          allocate(pgl(igroup)%amps(iamp)%sector_present(dim1,dim2))
          read(14) pgl(igroup)%amps(iamp)%sector_present
          read(14) dim1,dim2
          allocate(pgl(igroup)%amps(iamp)%sector_sign(dim1,dim2))
          read(14) pgl(igroup)%amps(iamp)%sector_sign
          read(14) dim1
          allocate(pgl(igroup)%amps(iamp)%iproc_start(dim1))
          read(14) pgl(igroup)%amps(iamp)%iproc_start
          read(14) dim1
          allocate(pgl(igroup)%amps(iamp)%n_sing(dim1))
          read(14) pgl(igroup)%amps(iamp)%n_sing
          read(14) dim1
          allocate(pgl(igroup)%amps(iamp)%n_qqbar(dim1))
          read(14) pgl(igroup)%amps(iamp)%n_qqbar
          read(14) dim1,dim2,dim3
          allocate(pgl(igroup)%amps(iamp)%spins(dim1,dim2,dim3))
          read(14) pgl(igroup)%amps(iamp)%spins
          read(14) dim1
          allocate(pgl(igroup)%amps(iamp)%amps(dim1))
          allocate(pgl(igroup)%amps(iamp)%amps_by_order(&
               pgl(igroup)%amps(iamp)%n_amps,pgl(igroup)%amps(iamp)%n_sectors))
       enddo
       ! multichannel
       read(14) pgl(igroup)%multichan%n_unique_channels,pgl(igroup)%multichan%n_unique_channelgroups
       read(14) dim1,dim2
       allocate(pgl(igroup)%multichan%unique_channelgroup_list(0:dim1-1,dim2))
       read(14) pgl(igroup)%multichan%unique_channelgroup_list
       read(14) dim1
       allocate(pgl(igroup)%multichan%unique_channel_list(dim1))
       read(14) pgl(igroup)%multichan%unique_channel_list
       read(14) dim1
       allocate(pgl(igroup)%multichan%map_proc_to_channelgroup(dim1))
       read(14)pgl(igroup)%multichan%map_proc_to_channelgroup
       read(14) dim1
       allocate(pgl(igroup)%multichan%number_of_channels(dim1))
       read(14) pgl(igroup)%multichan%number_of_channels
       read(14) dim1,dim2
       allocate(pgl(igroup)%multichan%channels(dim1,dim2))
       read(14) pgl(igroup)%multichan%channels
       pgl(igroup)%multichan%max_channels=dim1
       read(14) dim1,dim2,dim3
       allocate(pgl(igroup)%multichan%channel_permutations(dim1,dim2,dim3))
       read(14) pgl(igroup)%multichan%channel_permutations
       ! rest
       read(14) dim1,dim2
       allocate(pgl(igroup)%processes(dim1,dim2))
       read(14) pgl(igroup)%processes
       read(14) dim1,dim2
       allocate(pgl(igroup)%phase_space_permutations(dim1,dim2))
       read(14) pgl(igroup)%phase_space_permutations
       read(14) dim1
       allocate(pgl(igroup)%iden_iproc(dim1))
       read(14) pgl(igroup)%iden_iproc
       read(14) dim1
       allocate(pgl(igroup)%phase_space_orders(dim1))
       read(14) pgl(igroup)%phase_space_orders
       read(14) dim1
       allocate(pgl(igroup)%nhel(dim1))
       read(14) pgl(igroup)%nhel
       read(14) pgl(igroup)%nproc
       read(14) dim1,dim2
       allocate(pgl(igroup)%val_procs(dim1,dim2))
       read(14) pgl(igroup)%val_procs
       read(14) dim1,dim2
       allocate(pgl(igroup)%val_procs_abs(dim1,dim2))
       read(14) dim1,dim2
       allocate(pgl(igroup)%alias_factors(dim1,dim2))
       read(14) dim1,dim2
       allocate(pgl(igroup)%idenCOandMAPfactor(dim1,dim2))
       read(14) pgl(igroup)%idenCOandMAPfactor
       read(14) dim1,dim2,dim3
       allocate(pgl(igroup)%iden_processes(dim1,dim2,dim3))
       read(14) pgl(igroup)%iden_processes
       read(14) dim1,dim2
       allocate(pgl(igroup)%spin(dim1,dim2))
       read(14) pgl(igroup)%spin
       read(14) dim1,dim2
       allocate(pgl(igroup)%hel_fac(dim1,dim2))
       read(14) pgl(igroup)%hel_fac
       read(14) dim1
       allocate(pgl(igroup)%iden(dim1))
       read(14) pgl(igroup)%iden
       read(14) pgl(igroup)%ipdgs
       read(14) pgl(igroup)%next,pgl(igroup)%ndim
       allocate(pgl(igroup)%hel(pgl(igroup)%next))
       read(14) dim1
       allocate(pgl(igroup)%col_fac(dim1))
       read(14) pgl(igroup)%col_fac
       read(14) dim1
       allocate(pgl(igroup)%amp2(dim1))
       read(14) dim1
       allocate(pgl(igroup)%amp2_abs(dim1))
       read(14) dim1
       allocate(pgl(igroup)%amp2_hel(dim1))
       read(14) dim1
       allocate(pgl(igroup)%amp2_hel_abs(dim1))
       read(14) dim1
       allocate(pgl(igroup)%passed(dim1))
       read(14) pgl(igroup)%passed
       read(14) dim1,dim2
       allocate(pgl(igroup)%color_orders(dim1,dim2))
       read(14) pgl(igroup)%color_orders
    enddo
    close(14)
    call apply_final_state_widths_from_loaded_groups()
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
    current_model_signature=phys_model%model_signature()
    tolerance=1d-13*max(1d0,abs(current_model_signature),&
         abs(stored_model_signature))
    if (any(abs(current_model_signature-stored_model_signature).gt.tolerance)) then
       write (*,*) 'Amplitude library was created with incompatible model parameters.'
       write (*,*) 'Recreate it with --library=create and the current input card.'
       write (*,*) 'Stored model signature:',stored_model_signature
       write (*,*) 'Current model signature:',current_model_signature
       stop 1
    endif
  end subroutine read_amplitude_lib

  subroutine test_lib
    use amp_lib
    implicit none
    integer :: igroup,iint,i
    real(kind=8),dimension(:,:),allocatable :: p
    complex(kind=8),dimension(:),allocatable :: amps_save,amps
    complex(kind=8),dimension(:,:),allocatable :: sectors_save,sectors
    character(len=170) :: line,tmp
    do igroup=1,ngroups
       do iint=1,size(pgl(igroup)%amps)
          allocate(p(0:3,pgl(igroup)%next))
          allocate(amps(size(pgl(igroup)%amps(iint)%amps)))
          allocate(amps_save(size(pgl(igroup)%amps(iint)%amps)))
          allocate(sectors(pgl(igroup)%amps(iint)%n_amps,&
               pgl(igroup)%amps(iint)%n_sectors))
          allocate(sectors_save(pgl(igroup)%amps(iint)%n_amps,&
               pgl(igroup)%amps(iint)%n_sectors))
          write(tmp,*) igroup
          write(line,*) iint
          line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
          open(file=line,unit=14,form='unformatted',access='stream',status='old')
          read(14) p
          read(14) amps_save
          read(14) sectors_save
          close(14)
          call evaluate_amp(igroup,iint,p,amps)
          call evaluate_amp_by_order(igroup,iint,p,sectors)
          if (any(abs(amps_save-amps)/max(1d-30,abs(amps_save),abs(amps)).gt.1d-8)) then
             write (*,*) 'Process library not compatible with saved amplitudes',igroup,iint
             do i=1,size(pgl(igroup)%amps(iint)%amps)
                if (abs(amps_save(i)-amps(i))/(abs(amps_save(i))+abs(amps(i))).gt.1d-8) then
                   write (*,*) i,amps_save(i)
                   write (*,*) i,amps(i)
                endif
             enddo
             stop 1
          endif
          if (any(abs(sectors_save-sectors)/&
               max(1d-30,abs(sectors_save),abs(sectors)).gt.1d-8)) then
             write (*,*) 'Coupling-order sectors in library are incompatible with saved amplitudes',&
                  igroup,iint
             stop 1
          endif
          if (any(abs(amps-sum(sectors,dim=2))/&
               max(1d-30,abs(amps),abs(sum(sectors,dim=2))).gt.1d-8)) then
             write (*,*) 'Compiled library sectors do not reconstruct the complete amplitudes',&
                  igroup,iint
             stop 1
          endif
          deallocate(p)
          deallocate(amps)
          deallocate(amps_save)
          deallocate(sectors)
          deallocate(sectors_save)
       enddo
    enddo
  end subroutine test_lib
end module amplitude_library
