module amplitude_library
  use handling_processes
  use read_process_file
contains
  subroutine create_amplitude_lib(ncalls0,PS_choice)
    implicit none
    integer,intent(in) :: ncalls0,PS_choice
    character(len=170) :: tmp,line,filename
    integer :: igroup,j,iamp
    filename='library/amplib.f03'
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
    write(14,*) 'contains'
    write(14,*) 'subroutine evaluate_amp(ichan,iint,p,amps)'
    write(14,*) 'implicit none'
    write(14,*) 'integer :: ichan,iint'
    write(tmp,*) maxval(pgl(:)%next)
    write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//') :: p'
    write(14,*) 'complex(kind=8),dimension(*) :: amps'
    do igroup=1,ngroups
       if (igroup.eq.1) then
          write(14,*) 'if (ichan.eq.1) then'
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                write(14,*) 'if (iint.eq.1) then'
                write(14,*) 'call evaluate_amp1_1(amps,p)'
             else
                write(tmp,*) j
                write(14,*) 'elseif (iint.eq.'//trim(adjustl(tmp))//') then'
                write(14,*) 'call evaluate_amp1_'//trim(adjustl(tmp))//'(amps,p)'
             endif
          enddo
          write(14,*) 'endif'
       else
          write(tmp,*) igroup
          write(14,*) 'elseif (ichan.eq.'//trim(adjustl(tmp))//') then'
          do j=1,size(pgl(igroup)%amps)
             if (j.eq.1) then
                write(14,*) 'if (iint.eq.1) then'
                write(14,*) 'call evaluate_amp'//trim(adjustl(tmp))//'_1(amps,p)'
             else
                write(line,*) j
                write(14,*) 'elseif (iint.eq.'//trim(adjustl(line))//') then'
                write(14,*) 'call evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'(amps,p)'
             endif
          enddo
          write(14,*) 'endif'
       endif
    enddo
    write(14,*) 'endif'
    write(14,*) 'end subroutine evaluate_amp'
    write(14,*) 'end module amp_lib'
    close(14)
    filename='library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
    write(14) PS_choice
    write(14) ncalls0
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
       ! rest
       write(14) shape(pgl(igroup)%processes),pgl(igroup)%processes
       write(14) size(pgl(igroup)%iden_iproc),pgl(igroup)%iden_iproc
       write(14) size(pgl(igroup)%phase_space_orders),pgl(igroup)%phase_space_orders
       write(14) size(pgl(igroup)%nhel),pgl(igroup)%nhel
       write(14) pgl(igroup)%nproc
       write(14) shape(pgl(igroup)%val_procs),pgl(igroup)%val_procs
       write(14) shape(pgl(igroup)%idenCOandMAPfactor),pgl(igroup)%idenCOandMAPfactor
       write(14) shape(pgl(igroup)%iden_processes),pgl(igroup)%iden_processes
       write(14) shape(pgl(igroup)%spin),pgl(igroup)%spin
       write(14) shape(pgl(igroup)%hel_fac),pgl(igroup)%hel_fac
       write(14) size(pgl(igroup)%iden),pgl(igroup)%iden
       write(14) pgl(igroup)%ipdgs
       write(14) pgl(igroup)%next,pgl(igroup)%ndim
       write(14) size(pgl(igroup)%col_fac),pgl(igroup)%col_fac
       write(14) size(pgl(igroup)%amp2)
       write(14) size(pgl(igroup)%amp2_hel)
       write(14) size(pgl(igroup)%passed),pgl(igroup)%passed
       write(14) shape(pgl(igroup)%color_orders),pgl(igroup)%color_orders
    enddo
    close(14)
  end subroutine create_amplitude_lib

  subroutine read_amplitude_lib(ncalls0,PS_choice)
    implicit none
    integer,intent(out) :: ncalls0,PS_choice
    character(len=170) :: filename
    integer :: dim1,dim2,dim3,iamp,igroup
    filename='library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='old')
    read(14)PS_choice
    read(14)ncalls0
    allocate(pgl_unique)
    read(14) pgl_unique%next,pgl_unique%nproc
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
       ! rest
       read(14) dim1,dim2
       allocate(pgl(igroup)%processes(dim1,dim2))
       read(14) pgl(igroup)%processes
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
       allocate(pgl(igroup)%hel(next))
       read(14) dim1
       allocate(pgl(igroup)%col_fac(dim1))
       read(14) pgl(igroup)%col_fac
       read(14) dim1
       allocate(pgl(igroup)%amp2(dim1))
       read(14) dim1
       allocate(pgl(igroup)%amp2_hel(dim1))
       read(14) dim1
       allocate(pgl(igroup)%passed(dim1))
       read(14) pgl(igroup)%passed
       read(14) dim1,dim2
       allocate(pgl(igroup)%color_orders(dim1,dim2))
       read(14) pgl(igroup)%color_orders
    enddo
    close(14)
  end subroutine read_amplitude_lib

  subroutine test_lib
    use amp_lib
    implicit none
    integer :: igroup,iint
    real(kind=8),dimension(:,:),allocatable :: p
    complex(kind=8),dimension(:),allocatable :: amps_save,amps
    character(len=170) :: line,tmp
    do igroup=1,ngroups
       do iint=1,size(pgl(igroup)%amps)
          allocate(p(0:3,pgl(igroup)%next))
          allocate(amps(size(pgl(igroup)%amps(iint)%amps)))
          allocate(amps_save(size(pgl(igroup)%amps(iint)%amps)))
          write(tmp,*) igroup
          write(line,*) iint
          line='library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
          open(file=line,unit=14,form='unformatted',access='stream',status='old')
          read(14) p
          read(14) amps_save
          close(14)
          call evaluate_amp(igroup,iint,p,amps)
          if (any(abs(amps_save-amps)/(abs(amps_save)+abs(amps)).gt.1d-10)) then
             write (*,*) 'Process library not compatible with saved amplitudes',igroup,iint
             stop 1
          endif
          deallocate(p)
          deallocate(amps)
          deallocate(amps_save)
       enddo
    enddo
  end subroutine test_lib
end module amplitude_library
