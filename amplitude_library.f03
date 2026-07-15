module amplitude_library
  use handling_processes
  use read_process_file
contains
  subroutine create_amplitude_lib_f()
    use handling_events
    use common
    implicit none
    character(len=170) :: tmp,tmp2,line,filename
    integer :: igroup,j,iamp,iproc,ihel,ihel1,iden_iproc
    integer,dimension(2,pgl(1)%next) :: ICOLUP
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
    write(14,*) 'public :: evaluate_amp'
    write(14,*) 'contains'
    write(14,*) 'subroutine evaluate_amp(ichan,iint,p,amps)'
    write(14,*) 'implicit none'
    write(14,*) 'integer,intent(in) :: ichan,iint'
    write(tmp,*) maxval(pgl(:)%next)
    write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
    write(14,*) 'complex(kind=8),intent(out) :: amps(*)'
    do igroup=1,ngroups
      write(tmp,*) igroup
      write(line,*) 'elseif (ichan.eq.'//trim(adjustl(tmp))//') then'
      if(igroup.eq.1) line=line(6:)
      write(14,'(2x,a)') trim(adjustl(line))
      do j=1,size(pgl(igroup)%amps)
         write(tmp2,*) j
         write(line,*) 'elseif (iint.eq.'//trim(adjustl(tmp2))//') then'
         if(j.eq.1) line=line(6:)
         write(14,'(4x,a)') trim(adjustl(line))
         write(14,'(6x,a)') 'call evaluate_amp'//trim(adjustl(tmp))//'_'//trim(adjustl(tmp2))//'(p,amps)'
      enddo
      write(14,'(4x,a)') 'endif'
    enddo
    write(14,'(2x,a)') 'endif'
    write(14,*) 'end subroutine evaluate_amp'
    write(14,*) 'end module amp_lib'
    close(14)
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
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
    filename='Library/amplitudes_light.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
!!$    open(unit=14,file=filename,status='unknown')
    write(14) keep_processes_separate
    write(14) ngroups
    write(14) alphaEW
    do igroup=1,ngroups
       write(14) pgl(igroup)%next
       write(14) size(pgl(igroup)%amp2)
       write(14) pgl(igroup)%amps(1)%n_sing(1)
       write(14) pgl(igroup)%nproc
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac(:,1))
       else
          write(14) pgl(igroup)%nhel(1)
       endif
       write(14) pgl(igroup)%amps(1:size(pgl(igroup)%amps))%n_amps
       if (.not.keep_processes_separate) then
          write(14) size(pgl(igroup)%amps(1)%iproc_start),pgl(igroup)%amps(1)%iproc_start
       endif
       write(14) size(pgl(igroup)%col_fac),dble(pgl(igroup)%col_fac)
       write(14) size(pgl(igroup)%iden),dble(pgl(igroup)%iden)
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac),dble(pgl(igroup)%hel_fac)
       else
          write(14) pgl(igroup)%nhel(1),dble(pgl(igroup)%hel_fac(1:pgl(igroup)%nhel(1),1))
       endif
    enddo
    close(14)

    filename='madspace.json'
    open(unit=14,file=filename,status='unknown')
    write(14,*) '{'
    write(14,*) '  "channels": ['
    do igroup=1,ngroups
       write(14,*) '    {'
       write(14,*) '      "channel": ',igroup,','
       write(14,*) '      "phasespace_order": ['
       write(14,'("         ",*(I0,:,","))') pgl(igroup)%phase_space_orders
       write(14,*) '      ],'
       write(14,*) '      "processes": ['
       do iproc=1,pgl(igroup)%nproc
          write(14,*) '        {'
          write(14,*) '          "process":',iproc,','
          write(14,*) '          "color_order": ['
          write(14,'("             ",*(I0,:,","))') pgl(igroup)%color_orders(1:pgl(igroup)%next,iproc)
          write(14,*) '          ],'
          iproc_picked=iproc
          iproc_iden_picked=1
          call get_col_info(pgl(igroup),ICOLUP)
          write(14,*) '          "color_flows1": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(1,:)
          write(14,*) '          ],'
          write(14,*) '          "color_flows2": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(2,:)
          write(14,*) '          ],'
          write(14,*) '          "multichannels": ['
          write(14,'("             ",*(I0,:,","))') &
               pgl(igroup)%multichan%channels(1:pgl(igroup)%multichan%number_of_channels(iproc),iproc)
          write(14,*) '          ],'
          write(14,*) '          "matrix_elements": ['
          do iden_iproc=1,pgl(igroup)%iden_iproc(iproc)
             write(14,*) '            {'
             write(14,*) '              "label": ', iden_iproc, ','
             write(14,*) '              "factor": ', pgl(igroup)%idenCOandMAPfactor(iden_iproc,iproc), ','
             write(14,*) '              "pdg_ids": ['
             write(14,'("                 ",*(I0,:,","))') pgl(igroup)%iden_processes(1:pgl(igroup)%next,iden_iproc,iproc)
             write(14,*) '              ]'
             if (iden_iproc.eq.pgl(igroup)%iden_iproc(iproc)) then
                write(14,*) '            }'
             else
                write(14,*) '            },'
             endif
          enddo
          write(14,*) '          ],'
          write(14,*) '          "helicities": ['
          if (keep_processes_separate) then
             do ihel=pgl(igroup)%amps(iproc)%iproc_start(1),pgl(igroup)%amps(iproc)%iproc_start(2)-1
                write(14,*) '            ['
                do ihel1=1,pgl(igroup)%hel_fac(ihel,iproc)
                   write(14,*) '              ['
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(iproc)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,iproc)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(iproc)%iproc_start(2)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          else
             do ihel=pgl(igroup)%amps(1)%iproc_start(iproc),pgl(igroup)%amps(1)%iproc_start(iproc+1)-1
                do ihel1=1,pgl(igroup)%hel_fac(ihel,1)
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(1)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,1)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(1)%iproc_start(iproc+1)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          endif
          write(14,*) '          ]'
          if (iproc.eq.pgl(igroup)%nproc) then
             write(14,*) '        }'
          else
             write(14,*) '        },'
          endif
       enddo
       write(14,*) '      ]'
       if (igroup.eq.ngroups) then
          write(14,*) '    }'
       else
          write(14,*) '    },'
       endif
    enddo
    write(14,*) '  ],'
    write(14,'(a)', advance="no") '  "xml_header": "'
    write(14,'(a)', advance="no") '<header>\n'
    write(14,'(a)', advance="no") '   <unique_me>\n'
    write(14,'(I0, 1X, I0, A)', advance="no") pgl_unique%next,pgl_unique%nproc,"\n"
    do iproc=1,pgl_unique%nproc
       write(14,'(I0, 1X, F0.16, *(1X, I0))', advance="no") &
           unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:pgl_unique%next,iproc)
       write(14,'(a)', advance="no") '\n'
    enddo
    write(14,'(a)', advance="no") '   </unique_me>\n'
    write(14,'(a)', advance="no") '   <nevents> ... </nevents>\n'
    write(14,'(a)', advance="no") '   <seed>    ... </seed>\n'
    write(14,'(a)', advance="no") '</header>\n'
    write(14,*) '"'
    write(14,*) '}'
    close(14)
  end subroutine create_amplitude_lib_f

  subroutine create_amplitude_lib_c()
    use handling_events
    use common
    implicit none
    character(len=170) :: tmp,line,filename,headername
    character(len=170) :: tmp2,tmp3,namespace
    integer :: igroup,j,iamp,iproc,ihel,ihel1,iden_iproc
    integer,dimension(2,pgl(1)%next) :: ICOLUP
    integer :: max_amps, k, curr_amps
    filename='Library/amplibc.c'
    headername='Library/amplib.h'
    open(unit=14,file=filename,status='unknown')
   !  write(14,*) 'module amp_lib'
    max_amps=0
    open(unit=15,file=headername,status='unknown')
    write(15,*) '#pragma once'
    write(15,*) '#include "AmpliColTypes.h"'
    write(15,*) '#include "FeynmanRules.h"'
    write(14,*) '#include "amplib.h"'
    do igroup=1,ngroups
       do j=1,size(pgl(igroup)%amps)
          write(tmp,*) igroup
          line=trim(adjustl(tmp))//'_'
          write(tmp,*) j
          line=trim(adjustl(line))//trim(adjustl(tmp))
         !  write(14,*) 'use amp'//trim(adjustl(line))//'_lib'
          write(15,*) 'extern const struct AC_AMP amp'//trim(adjustl(line))//';'
          do k=1,size(pgl(igroup)%amps)
             if (pgl(igroup)%amps(k)%n_amps.gt.max_amps) max_amps=pgl(igroup)%amps(k)%n_amps
          end do
       enddo
    enddo
   !  write(14,*) 'implicit none'
   !  write(14,*) 'private'
   !  write(14,*) 'public :: evaluate_amp'
   !  write(14,*) 'contains'
    write(tmp,*) maxval(pgl(:)%next)
    write(tmp2,*) max_amps
    write(14,*) ''
    write(14,*) 'extern void evaluate_amp(int const* ichan,int const* iint, AC_D_FP const*' // &
     ' p, AC_D_CX *amps) {'
     write(14,'(2x,a)') '// AC_D_FP p_arr['//trim(adjustl(tmp))//'][4];'
      write(14,'(2x,a)') '// AC_D_CX amps_arr['//trim(adjustl(tmp2))//'];'
      ! write(14,'(2x,a)') 'for (int j=0; j<'//trim(adjustl(tmp))//'; j++) for (int i=0; i<4; i++) p_arr[j][i]=p[4*j + i];'
      write(14,'(2x,a)') 'evaluate_amp_impl(*ichan,*iint,(const AC_D_FP (*)[4])p,amps);'
      ! write(14,'(2x,a)') 'for (int i=0; i<'//trim(adjustl(tmp2))//'; i++) amps[i]=amps_arr[i];'
      write(14,'(a)') '}'
      write(14,*) ''
   write(14,*) 'void evaluate_amp_impl(int ichan,int iint, const AC_D_FP p['//trim(adjustl(tmp))//'][4]'// &
    ', AC_D_CX amps['//trim(adjustl(tmp2))//']) {'
    write(15,*) ''
    write(15,*) 'extern void evaluate_amp(int const* ichan,int const* iint, AC_D_FP const*'// &
     ' p, AC_D_CX *amps);'
    write(15,*) ''
   write(15,*) 'void evaluate_amp_impl(int ichan,int iint, const AC_D_FP p['//trim(adjustl(tmp))//'][4]'// &
    ', AC_D_CX amps['//trim(adjustl(tmp2))//']);'
    close(15)
   !  write(14,*) 'implicit none'
   !  write(14,*) 'integer,intent(in) :: ichan,iint'
   !  write(tmp,*) maxval(pgl(:)%next)
   !  write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
   !  write(14,*) 'complex(kind=8),intent(out) :: amps(*)'
    do igroup=1,ngroups
      write(tmp,*) igroup
      write(line,*) 'else if (ichan == '//trim(adjustl(tmp))//')'
      if(igroup.eq.1) line=line(7:)
      write(14,'(2x,a)') trim(adjustl(line))
      write(14,'(2x,a)') '{'
      write(14,'(4x,a)') 'switch (iint) {'
      do j=1,size(pgl(igroup)%amps)
         curr_amps=pgl(igroup)%amps(j)%n_amps
         write(tmp3,*) curr_amps
         write(line,*) j
         write(namespace,*) 'amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'.'
         write(14,'(6x,a)') 'case '//trim(adjustl(line))//': '//trim(adjustl(namespace))//'eval(p,amps); break;'
      enddo
      write(14,'(4x,a)') '}'
      write(14,'(2x,a)') '}'
    enddo
    write(14,*) '} // end of evaluate_amp'
    close(14)
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
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
    filename='Library/amplitudes_light.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
!!$    open(unit=14,file=filename,status='unknown')
    write(14) keep_processes_separate
    write(14) ngroups
    write(14) alphaEW
    do igroup=1,ngroups
       write(14) pgl(igroup)%next
       write(14) size(pgl(igroup)%amp2)
       write(14) pgl(igroup)%amps(1)%n_sing(1)
       write(14) pgl(igroup)%nproc
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac(:,1))
       else
          write(14) pgl(igroup)%nhel(1)
       endif
       write(14) pgl(igroup)%amps(1:size(pgl(igroup)%amps))%n_amps
       if (.not.keep_processes_separate) then
          write(14) size(pgl(igroup)%amps(1)%iproc_start),pgl(igroup)%amps(1)%iproc_start
       endif
       write(14) size(pgl(igroup)%col_fac),dble(pgl(igroup)%col_fac)
       write(14) size(pgl(igroup)%iden),dble(pgl(igroup)%iden)
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac),dble(pgl(igroup)%hel_fac)
       else
          write(14) pgl(igroup)%nhel(1),dble(pgl(igroup)%hel_fac(1:pgl(igroup)%nhel(1),1))
       endif
    enddo
    close(14)

    filename='madspace.json'
    open(unit=14,file=filename,status='unknown')
    write(14,*) '{'
    write(14,*) '  "channels": ['
    do igroup=1,ngroups
       write(14,*) '    {'
       write(14,*) '      "channel": ',igroup,','
       write(14,*) '      "phasespace_order": ['
       write(14,'("         ",*(I0,:,","))') pgl(igroup)%phase_space_orders
       write(14,*) '      ],'
       write(14,*) '      "processes": ['
       do iproc=1,pgl(igroup)%nproc
          write(14,*) '        {'
          write(14,*) '          "process":',iproc,','
          write(14,*) '          "color_order": ['
          write(14,'("             ",*(I0,:,","))') pgl(igroup)%color_orders(1:pgl(igroup)%next,iproc)
          write(14,*) '          ],'
          iproc_picked=iproc
          iproc_iden_picked=1
          call get_col_info(pgl(igroup),ICOLUP)
          write(14,*) '          "color_flows1": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(1,:)
          write(14,*) '          ],'
          write(14,*) '          "color_flows2": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(2,:)
          write(14,*) '          ],'
          write(14,*) '          "multichannels": ['
          write(14,'("             ",*(I0,:,","))') &
               pgl(igroup)%multichan%channels(1:pgl(igroup)%multichan%number_of_channels(iproc),iproc)
          write(14,*) '          ],'
          write(14,*) '          "matrix_elements": ['
          do iden_iproc=1,pgl(igroup)%iden_iproc(iproc)
             write(14,*) '            {'
             write(14,*) '              "label": ', iden_iproc, ','
             write(14,*) '              "factor": ', pgl(igroup)%idenCOandMAPfactor(iden_iproc,iproc), ','
             write(14,*) '              "pdg_ids": ['
             write(14,'("                 ",*(I0,:,","))') pgl(igroup)%iden_processes(1:pgl(igroup)%next,iden_iproc,iproc)
             write(14,*) '              ]'
             if (iden_iproc.eq.pgl(igroup)%iden_iproc(iproc)) then
                write(14,*) '            }'
             else
                write(14,*) '            },'
             endif
          enddo
          write(14,*) '          ],'
          write(14,*) '          "helicities": ['
          if (keep_processes_separate) then
             do ihel=pgl(igroup)%amps(iproc)%iproc_start(1),pgl(igroup)%amps(iproc)%iproc_start(2)-1
                write(14,*) '            ['
                do ihel1=1,pgl(igroup)%hel_fac(ihel,iproc)
                   write(14,*) '              ['
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(iproc)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,iproc)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(iproc)%iproc_start(2)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          else
             do ihel=pgl(igroup)%amps(1)%iproc_start(iproc),pgl(igroup)%amps(1)%iproc_start(iproc+1)-1
                do ihel1=1,pgl(igroup)%hel_fac(ihel,1)
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(1)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,1)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(1)%iproc_start(iproc+1)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          endif
          write(14,*) '          ]'
          if (iproc.eq.pgl(igroup)%nproc) then
             write(14,*) '        }'
          else
             write(14,*) '        },'
          endif
       enddo
       write(14,*) '      ]'
       if (igroup.eq.ngroups) then
          write(14,*) '    }'
       else
          write(14,*) '    },'
       endif
    enddo
    write(14,*) '  ],'
    write(14,'(a)', advance="no") '  "xml_header": "'
    write(14,'(a)', advance="no") '<header>\n'
    write(14,'(a)', advance="no") '   <unique_me>\n'
    write(14,'(I0, 1X, I0, A)', advance="no") pgl_unique%next,pgl_unique%nproc,"\n"
    do iproc=1,pgl_unique%nproc
       write(14,'(I0, 1X, F0.16, *(1X, I0))', advance="no") &
           unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:pgl_unique%next,iproc)
       write(14,'(a)', advance="no") '\n'
    enddo
    write(14,'(a)', advance="no") '   </unique_me>\n'
    write(14,'(a)', advance="no") '   <nevents> ... </nevents>\n'
    write(14,'(a)', advance="no") '   <seed>    ... </seed>\n'
    write(14,'(a)', advance="no") '</header>\n'
    write(14,*) '"'
    write(14,*) '}'
    close(14)
  end subroutine create_amplitude_lib_c

  subroutine create_amplitude_lib_cu()
    use handling_events
    use common
    implicit none
    character(len=170) :: tmp,line,filename,headername
    character(len=170) :: tmp2,tmp3,namespace,groupname
    integer :: igroup,j,iamp,iproc,ihel,ihel1,iden_iproc
    integer,dimension(2,pgl(1)%next) :: ICOLUP
    integer :: max_amps, k
    integer :: max_curr, max_vert, max_pmom
    filename='Library/amplibcu.cu'
    headername='Library/amplib.cuh'
    open(unit=14,file=filename,status='unknown')
   !  write(14,*) 'module amp_lib'
    max_amps=0
    max_curr=0
    max_vert=0
    max_pmom=0
    open(unit=15,file=headername,status='unknown')
    write(15,*) '#pragma once'
    write(15,*) '#include "AmpliColTypes.h"'
    write(15,*) '#include <cuda_runtime.h>'
    write(14,*) '#include "amplib.cuh"'
    write(14,*) '#include <cstdio>'
    write(14,*) ''
    do igroup=1,ngroups
       write(tmp,*) igroup
       groupname=trim(adjustl(tmp))
       filename='Library/amp'//trim(groupname)//'_libcu.cu'
       open(unit=16,file=filename,status='unknown')
       write(16,*) '#include "amplib.cuh"'
       write(16,*) ''
       do j=1,size(pgl(igroup)%amps)
          line=trim(groupname)//'_'
          write(tmp,*) j
          line=trim(adjustl(line))//trim(adjustl(tmp))
          write(16,*) '#include "amp'//trim(adjustl(line))//'_lib.cuh"'
          do k=1,size(pgl(igroup)%amps)
             if (pgl(igroup)%amps(k)%n_amps.gt.max_amps) max_amps=pgl(igroup)%amps(k)%n_amps
             if (pgl(igroup)%amps(k)%n_cur.gt.max_curr) max_curr=pgl(igroup)%amps(k)%n_cur
             if (pgl(igroup)%amps(k)%n_vert.gt.max_vert) max_vert=pgl(igroup)%amps(k)%n_vert
             if (pgl(igroup)%amps(k)%max_pp.gt.max_pmom) max_pmom=pgl(igroup)%amps(k)%max_pp
          end do
       enddo
       write(16,*) ''
       write(16,*) '__device__ void evaluate_amp'//trim(groupname)//'(int iint, const AC_D_FP p[AC_NEXT][4],'// &
        ' AC_D_CX amps[AC_NAMP], AC_D_FP pp[][4], AC_CX val_c[][6], AC_CX int_c[][6]) {'
       write(16,'(2x,a)') 'switch (iint) {'
       do j=1,size(pgl(igroup)%amps)
          write(line,*) j
          write(namespace,*) 'amp'//trim(groupname)//'_'//trim(adjustl(line))
          write(16,'(4x,a)') 'case '//trim(adjustl(line))//': '//trim(adjustl(namespace))//'::evaluate_' &
          //trim(adjustl(namespace))//'(p,amps,pp,val_c,int_c); break;'
       enddo
       write(16,'(2x,a)') '}'
       write(16,*) '} // end of evaluate_amp'//trim(groupname)
       close(16)
    enddo
   !  write(14,*) 'implicit none'
   !  write(14,*) 'private'
   !  write(14,*) 'public :: evaluate_amp'
   !  write(14,*) 'contains'
    write(tmp,*) maxval(pgl(:)%next)
    write(tmp2,*) max_amps
    write(15,*) ''
    write(15,*) '#define AC_NEXT '//trim(adjustl(tmp))
    write(15,*) '#define AC_NAMP '//trim(adjustl(tmp2))//' //max number of amplitudes'
    write(tmp,*) max_curr
    write(tmp2,*) max_vert
    write(tmp3,*) max_pmom
    write(15,*) '#define AC_NCURR '//trim(adjustl(tmp))
    write(15,*) '#define AC_NVERT '//trim(adjustl(tmp2))//' //max number of vertices'
    write(15,*) '#define AC_NPMOM '//trim(adjustl(tmp3))//' //max number of momenta'
    write(15,*) ''
    write(15,*) '#define MOM_STEP ((size_t) 4 * AC_NPMOM * sizeof(AC_D_FP))   // memory size for maximal internal momenta'
    write(15,*) '#define CUR_STEP ((size_t) 6 * AC_NCURR * sizeof(AC_CX))   // memory size for maximal internal currents'
    write(15,*) '#define VER_STEP ((size_t) 6 * AC_NVERT * sizeof(AC_CX))   // memory size for maximal internal vertices'
    write(15,*) '#define STRIDE   (MOM_STEP + CUR_STEP + VER_STEP)          // total memory size per amplitude'
    write(14,*) ''
    do igroup=1,ngroups
       write(tmp,*) igroup
       write(14,*) '__device__ void evaluate_amp'//trim(adjustl(tmp))//'(int iint, const AC_FP p[AC_NEXT][4],'// &
        ' AC_CX amps[AC_NAMP], AC_FP pp[][4], AC_CX val_c[][6], AC_CX int_c[][6]);'
    enddo
    write(14,*) ''
    write(14,*) '__device__ void evaluate_amp_impl(int ichan,int iint, const AC_D_FP p[AC_NEXT][4]'// &
    ', AC_D_CX amps[AC_NAMP], AC_D_FP pp[][4], AC_CX val_c[][6], AC_CX int_c[][6]) {'
    write(15,*) ''
    write(15,*) 'extern "C" void evaluate_amp(int const* ichan,int const* iint, AC_D_FP const*'// &
     ' p, AC_D_CX *amps, int N);'
    write(15,*) ''
    write(15,*) 'void evaluate_amp_device(const int* d_ichan,const int* d_iint, const AC_D_FP* d_p, AC_D_CX* d_amps,'// &
     ' int N, char* d_work, cudaStream_t stream = 0);'
     write(15,*) ''
    write(15,*) '__device__ void evaluate_amp_impl(int ichan,int iint, const AC_D_FP p['// &
     'AC_NEXT][4], AC_D_CX amps[AC_NAMP], AC_D_FP pp[][4], AC_CX val_c[][6], AC_CX int_c[][6]);'
    close(15)
   !  write(14,*) 'implicit none'
   !  write(14,*) 'integer,intent(in) :: ichan,iint'
   !  write(tmp,*) maxval(pgl(:)%next)
   !  write(14,*) 'real(kind=8),dimension(0:3,'//trim(adjustl(tmp))//'),intent(in) :: p'
   !  write(14,*) 'complex(kind=8),intent(out) :: amps(*)'
    do igroup=1,ngroups
      write(tmp,*) igroup
      write(line,*) 'else if (ichan == '//trim(adjustl(tmp))//')'
      if(igroup.eq.1) line=line(7:)
      write(14,'(2x,a)') trim(adjustl(line))
      write(14,'(2x,a)') '{'
      write(14,'(4x,a)') 'evaluate_amp'//trim(adjustl(tmp))//'(iint,p,amps,pp,val_c,int_c);'
      write(14,'(2x,a)') '}'
    enddo
    write(14,*) '} // end of evaluate_amp_impl'
    write(14,*) ''
    write(14,*) '__global__ void evaluate_amp_kernel(const int* d_ichan,const int* d_iint,'// &
     ' const AC_D_FP* d_p, AC_D_CX* d_amps, char* d_work, int N)'
    write(14,*) '{'
    write(14,'(2x,a)') 'for (int i=blockIdx.x*blockDim.x + threadIdx.x; i<N; i+=gridDim.x*blockDim.x)'
    write(14,'(2x,a)') '{'
    write(14,'(4x,a)') 'const AC_D_FP (*p)[4] = (const AC_D_FP (*)[4])&d_p[i*4*AC_NEXT];'
    write(14,'(4x,a)') 'AC_D_CX *amps = &d_amps[i*AC_NAMP];'
    write(14,'(4x,a)') 'char* work = &d_work[i*STRIDE];'
    write(14,'(4x,a)') 'AC_D_FP (*pp)[4] = (AC_D_FP (*)[4])work;'
    write(14,'(4x,a)') 'AC_CX (*val_c)[6] = (AC_CX (*)[6])&work[MOM_STEP];'
    write(14,'(4x,a)') 'AC_CX (*int_c)[6] = (AC_CX (*)[6])&work[MOM_STEP + CUR_STEP];'
    write(14,'(4x,a)') 'evaluate_amp_impl(d_ichan[i],d_iint[i],p,amps,pp,val_c,int_c);'
    write(14,'(2x,a)') '}'
    write(14,*) '}'
    write(14,*) ''
    write(14,*) 'void evaluate_amp_device(const int* d_ichan,const int* d_iint, const AC_D_FP* d_p, AC_D_CX* d_amps,'// &
     ' int N, char* d_work, cudaStream_t stream)'
    write(14,*) '{'
    write(14,'(2x,a)') 'int blockSize = 256;'
    write(14,'(2x,a)') 'int numBlocks = (N + blockSize - 1) / blockSize;'
    write(14,'(2x,a)') 'evaluate_amp_kernel<<<numBlocks, blockSize, 0, stream>>>(d_ichan, d_iint, d_p, d_amps, d_work, N);'
    write(14,'(2x,a)') 'cudaError_t err = cudaGetLastError();'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'
    write(14,'(4x,a)') 'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err));'
    write(14,'(2x,a)') '}'
    write(14,*) '}'
    write(14,*) ''
    write(14,*) 'extern "C" void evaluate_amp(int const* ichan,int const* iint, AC_D_FP const*' // &
     ' p, AC_D_CX *amps, int N) {'
    write(14,'(2x,a)') 'int *d_ichan, *d_iint;'
    write(14,'(2x,a)') 'AC_D_FP *d_p;'
    write(14,'(2x,a)') 'AC_D_CX *d_amps;'
    write(14,'(2x,a)') 'char *d_work;'
    write(14,'(2x,a)') 'cudaError_t err;'
    write(14,'(2x,a)') 'err = cudaMalloc(&d_ichan, N * sizeof(int));'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
     'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); return; }'
    write(14,'(2x,a)') 'err = cudaMalloc(&d_iint, N * sizeof(int));'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
      'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); return; }'
    write(14,'(2x,a)') 'err = cudaMalloc(&d_p, (size_t)N * 4 * AC_NEXT * sizeof(AC_D_FP));'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
      'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); return; }'
    write(14,'(2x,a)') 'err = cudaMalloc(&d_amps, (size_t)N * AC_NAMP * sizeof(AC_D_CX));'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
      'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); return; }'
    write(14,'(2x,a)') 'err = cudaMalloc(&d_work, (size_t)N * STRIDE);'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
      'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); return; }'
    write(14,'(2x,a)') ''
    write(14,'(2x,a)') 'cudaMemset(d_amps, 0, (size_t)N * AC_NAMP * sizeof(AC_D_CX));'
    write(14,'(2x,a)') ''
    write(14,'(2x,a)') 'cudaMemcpy(d_ichan, ichan, N * sizeof(int), cudaMemcpyHostToDevice);'
    write(14,'(2x,a)') 'cudaMemcpy(d_iint, iint, N * sizeof(int), cudaMemcpyHostToDevice);'
    write(14,'(2x,a)') 'cudaMemcpy(d_p, p, (size_t)N * 4 * AC_NEXT * sizeof(AC_D_FP), cudaMemcpyHostToDevice);'
    write(14,'(2x,a)') ''
    write(14,'(2x,a)') 'evaluate_amp_device(d_ichan, d_iint, d_p, d_amps, N, d_work);'
    write(14,'(2x,a)') ''
    write(14,'(2x,a)') 'err = cudaDeviceSynchronize();'
    write(14,'(2x,a)') 'if (err != cudaSuccess) {'//&
      'fprintf(stderr, "CUDA error: %s\\n", cudaGetErrorString(err)); }'
    write(14,'(2x,a)') ''
    write(14,'(2x,a)') 'cudaMemcpy(amps, d_amps, (size_t)N * AC_NAMP * sizeof(AC_D_CX), cudaMemcpyDeviceToHost);'
    write(14,'(2x,a)') 'cudaFree(d_ichan);'
    write(14,'(2x,a)') 'cudaFree(d_iint);'
    write(14,'(2x,a)') 'cudaFree(d_p);'
    write(14,'(2x,a)') 'cudaFree(d_amps);'
    write(14,'(2x,a)') 'cudaFree(d_work);'
    write(14,'(a)') '}'

    close(14)
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
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
    filename='Library/amplitudes_light.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='unknown')
!!$    open(unit=14,file=filename,status='unknown')
    write(14) keep_processes_separate
    write(14) ngroups
    write(14) alphaEW
    do igroup=1,ngroups
       write(14) pgl(igroup)%next
       write(14) size(pgl(igroup)%amp2)
       write(14) pgl(igroup)%amps(1)%n_sing(1)
       write(14) pgl(igroup)%nproc
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac(:,1))
       else
          write(14) pgl(igroup)%nhel(1)
       endif
       write(14) pgl(igroup)%amps(1:size(pgl(igroup)%amps))%n_amps
       if (.not.keep_processes_separate) then
          write(14) size(pgl(igroup)%amps(1)%iproc_start),pgl(igroup)%amps(1)%iproc_start
       endif
       write(14) size(pgl(igroup)%col_fac),dble(pgl(igroup)%col_fac)
       write(14) size(pgl(igroup)%iden),dble(pgl(igroup)%iden)
       if (keep_processes_separate) then
          write(14) size(pgl(igroup)%hel_fac),dble(pgl(igroup)%hel_fac)
       else
          write(14) pgl(igroup)%nhel(1),dble(pgl(igroup)%hel_fac(1:pgl(igroup)%nhel(1),1))
       endif
    enddo
    close(14)

    filename='madspace.json'
    open(unit=14,file=filename,status='unknown')
    write(14,*) '{'
    write(14,*) '  "channels": ['
    do igroup=1,ngroups
       write(14,*) '    {'
       write(14,*) '      "channel": ',igroup,','
       write(14,*) '      "phasespace_order": ['
       write(14,'("         ",*(I0,:,","))') pgl(igroup)%phase_space_orders
       write(14,*) '      ],'
       write(14,*) '      "processes": ['
       do iproc=1,pgl(igroup)%nproc
          write(14,*) '        {'
          write(14,*) '          "process":',iproc,','
          write(14,*) '          "color_order": ['
          write(14,'("             ",*(I0,:,","))') pgl(igroup)%color_orders(1:pgl(igroup)%next,iproc)
          write(14,*) '          ],'
          iproc_picked=iproc
          iproc_iden_picked=1
          call get_col_info(pgl(igroup),ICOLUP)
          write(14,*) '          "color_flows1": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(1,:)
          write(14,*) '          ],'
          write(14,*) '          "color_flows2": ['
          write(14,'("             ",*(I0,:,","))') ICOLUP(2,:)
          write(14,*) '          ],'
          write(14,*) '          "multichannels": ['
          write(14,'("             ",*(I0,:,","))') &
               pgl(igroup)%multichan%channels(1:pgl(igroup)%multichan%number_of_channels(iproc),iproc)
          write(14,*) '          ],'
          write(14,*) '          "matrix_elements": ['
          do iden_iproc=1,pgl(igroup)%iden_iproc(iproc)
             write(14,*) '            {'
             write(14,*) '              "label": ', iden_iproc, ','
             write(14,*) '              "factor": ', pgl(igroup)%idenCOandMAPfactor(iden_iproc,iproc), ','
             write(14,*) '              "pdg_ids": ['
             write(14,'("                 ",*(I0,:,","))') pgl(igroup)%iden_processes(1:pgl(igroup)%next,iden_iproc,iproc)
             write(14,*) '              ]'
             if (iden_iproc.eq.pgl(igroup)%iden_iproc(iproc)) then
                write(14,*) '            }'
             else
                write(14,*) '            },'
             endif
          enddo
          write(14,*) '          ],'
          write(14,*) '          "helicities": ['
          if (keep_processes_separate) then
             do ihel=pgl(igroup)%amps(iproc)%iproc_start(1),pgl(igroup)%amps(iproc)%iproc_start(2)-1
                write(14,*) '            ['
                do ihel1=1,pgl(igroup)%hel_fac(ihel,iproc)
                   write(14,*) '              ['
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(iproc)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,iproc)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(iproc)%iproc_start(2)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          else
             do ihel=pgl(igroup)%amps(1)%iproc_start(iproc),pgl(igroup)%amps(1)%iproc_start(iproc+1)-1
                do ihel1=1,pgl(igroup)%hel_fac(ihel,1)
                   write(14,'("                 ",*(I0,:,","))') pgl(igroup)%amps(1)%spins(1:pgl(igroup)%next,ihel1,ihel)
                   if (ihel1.eq.pgl(igroup)%hel_fac(ihel,1)) then
                      write(14,*) '              ]'
                   else
                      write(14,*) '              ],'
                   endif
                enddo
                if (ihel.eq.pgl(igroup)%amps(1)%iproc_start(iproc+1)-1) then
                   write(14,*) '            ]'
                else
                   write(14,*) '            ],'
                endif
             enddo
          endif
          write(14,*) '          ]'
          if (iproc.eq.pgl(igroup)%nproc) then
             write(14,*) '        }'
          else
             write(14,*) '        },'
          endif
       enddo
       write(14,*) '      ]'
       if (igroup.eq.ngroups) then
          write(14,*) '    }'
       else
          write(14,*) '    },'
       endif
    enddo
    write(14,*) '  ],'
    write(14,'(a)', advance="no") '  "xml_header": "'
    write(14,'(a)', advance="no") '<header>\n'
    write(14,'(a)', advance="no") '   <unique_me>\n'
    write(14,'(I0, 1X, I0, A)', advance="no") pgl_unique%next,pgl_unique%nproc,"\n"
    do iproc=1,pgl_unique%nproc
       write(14,'(I0, 1X, F0.16, *(1X, I0))', advance="no") &
           unique_map(iproc),unique_map_value(iproc),pgl_unique%processes(1:pgl_unique%next,iproc)
       write(14,'(a)', advance="no") '\n'
    enddo
    write(14,'(a)', advance="no") '   </unique_me>\n'
    write(14,'(a)', advance="no") '   <nevents> ... </nevents>\n'
    write(14,'(a)', advance="no") '   <seed>    ... </seed>\n'
    write(14,'(a)', advance="no") '</header>\n'
    write(14,*) '"'
    write(14,*) '}'
    close(14)
  end subroutine create_amplitude_lib_cu

  subroutine read_amplitude_lib()
    implicit none
    character(len=170) :: filename
    integer :: dim1,dim2,dim3,iamp,igroup
    filename='Library/amplitudes.bin'
    open(unit=14,file=filename,form='unformatted',access='stream',status='old')
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
    integer :: igroup,iint,i
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
          line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
          open(file=line,unit=14,form='unformatted',access='stream',status='old')
          read(14) p
          read(14) amps_save
          close(14)
          call evaluate_amp(igroup,iint,p,amps)
          if (any(abs(amps_save-amps)/(abs(amps_save)+abs(amps)).gt.1d-8)) then
             write (*,*) 'Process library not compatible with saved amplitudes',igroup,iint
             do i=1,size(pgl(igroup)%amps(iint)%amps)
                if (abs(amps_save(i)-amps(i))/(abs(amps_save(i))+abs(amps(i))).gt.1d-8) then
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
