program amplicol_library_benchmark
  use amplitude_library, only: read_amplitude_lib
  use amp_lib, only: evaluate_amp
  use common, only: alphaEW,keep_processes_separate
  use handling_processes, only: ngroups,pgl
  implicit none

  integer(kind=8) :: points, i
  integer :: igroup, iint, ih, iproc, argc, j
  real(kind=8) :: t0, t1, t_eval, t_square, t_total
  real(kind=8),dimension(:,:),allocatable :: p
  complex(kind=8),dimension(:),allocatable :: amps, amps_save
  real(kind=8),dimension(:),allocatable :: amp2, amp2_hel
  real(kind=8) :: checksum, norm_factor, selected_value
  character(len=170) :: line,tmp
  character(len=1024) :: arg, momenta_file
  real(kind=8),parameter :: alpha_check=0.118d0
  real(kind=8),parameter :: pi=3.14159265358979323846d0

  points = 100000_8
  igroup = 1
  iint = 1
  momenta_file = ''

  argc = command_argument_count()
  if (argc >= 1) then
     call get_command_argument(1,arg)
     read(arg,*) points
  endif
  if (argc >= 2) then
     call get_command_argument(2,arg)
     read(arg,*) igroup
  endif
  if (argc >= 3) then
     call get_command_argument(3,arg)
     read(arg,*) iint
  endif
  if (argc >= 4) then
     call get_command_argument(4,arg)
     read(arg,'(a)') momenta_file
  endif
  if (points < 1_8) then
     write (*,*) 'points must be positive'
     stop 1
  endif

  call read_amplitude_lib()
  if (igroup < 1 .or. igroup > ngroups) then
     write (*,*) 'invalid group',igroup,'ngroups=',ngroups
     stop 1
  endif
  if (iint < 1 .or. iint > size(pgl(igroup)%amps)) then
     write (*,*) 'invalid integral',iint,'n_integrals=',size(pgl(igroup)%amps)
     stop 1
  endif

  allocate(p(0:3,pgl(igroup)%next))
  allocate(amps(size(pgl(igroup)%amps(iint)%amps)))
  allocate(amps_save(size(pgl(igroup)%amps(iint)%amps)))
  if (keep_processes_separate) then
     allocate(amp2(1))
  else
     allocate(amp2(1:pgl(igroup)%nproc))
  endif
  allocate(amp2_hel(1:maxval(pgl(igroup)%nhel)))

  if (len_trim(momenta_file).eq.0) then
     write(tmp,*) igroup
     write(line,*) iint
     line='Library/amp'//trim(adjustl(tmp))//'_'//trim(adjustl(line))//'_lib.data'
     open(file=line,unit=14,form='unformatted',access='stream',status='old')
     read(14) p
     read(14) amps_save
     close(14)
  else
     open(file=trim(adjustl(momenta_file)),unit=14,status='old')
     do j=1,pgl(igroup)%next
        read(14,*) p(0,j),p(1,j),p(2,j),p(3,j)
     enddo
     close(14)
  endif

  call evaluate_amp(igroup,iint,p,amps)
  if (len_trim(momenta_file).eq.0 .and. &
       any(abs(amps_save-amps).gt.1d-20+1d-8*(abs(amps_save)+abs(amps)))) then
     write (*,*) 'Process library not compatible with saved amplitudes',igroup,iint
     stop 1
  endif
  call square_amplitudes(igroup,iint,amps,amp2,amp2_hel)
  norm_factor = (4*pi*alpha_check)**(pgl(igroup)%next-2-&
       pgl(igroup)%amps(iint)%n_sing(1))&
       *(2d0*4d0*pi*alphaEW)**pgl(igroup)%amps(iint)%n_sing(1)&
       /dble(pgl(igroup)%iden(iint))
  selected_value = norm_factor*amp2(1)

  checksum = 0d0
  call cpu_time(t0)
  do i=1,points
     call evaluate_amp(igroup,iint,p,amps)
     checksum = checksum + dble(amps(1)) + aimag(amps(1))
  enddo
  call cpu_time(t1)
  t_eval = t1 - t0

  call cpu_time(t0)
  do i=1,points
     call square_amplitudes(igroup,iint,amps,amp2,amp2_hel)
     checksum = checksum + sum(amp2)
  enddo
  call cpu_time(t1)
  t_square = t1 - t0
  t_total = t_eval + t_square

  write (*,'(a)') 'AmpliCol direct generated-library benchmark'
  write (*,'(a,i0)') 'points ',points
  write (*,'(a,i0)') 'group ',igroup
  write (*,'(a,i0)') 'integral ',iint
  write (*,'(a,i0)') 'amplitudes ',size(amps)
  write (*,'(a,es24.16)') 'checksum ',checksum
  write (*,'(a,1x,i0,1x,i0,1x,es24.16)') &
       'AMPICOL_SELECTED_FLOW_PROBE_VALUE',igroup,iint,selected_value
  write (*,'(a,1x,i0,1x,*(i0,1x))') &
       'AMPICOL_SELECTED_FLOW_PROBE_PDGS',pgl(igroup)%next,&
       pgl(igroup)%processes(1:pgl(igroup)%next,iint)
  write (*,'(a,1x,i0,1x,*(i0,1x))') &
       'AMPICOL_SELECTED_FLOW_PROBE_COLOR_ORDER',pgl(igroup)%next,&
       pgl(igroup)%color_orders(1:pgl(igroup)%next,iint)
  write (*,'(a,1x,i0)') &
       'AMPICOL_SELECTED_FLOW_PROBE_AMPLITUDES',size(amps)
  write (*,'(a,1x,i0)') &
       'AMPICOL_SELECTED_FLOW_PROBE_COLOR_FACTOR',pgl(igroup)%col_fac(iint)
  write (*,'(a,1x,i0)') &
       'AMPICOL_SELECTED_FLOW_PROBE_IDENTICAL_FACTOR',pgl(igroup)%iden(iint)
  write (*,'(a,1x,i0)') &
       'AMPICOL_SELECTED_FLOW_PROBE_SINGLET_VERTICES',&
       pgl(igroup)%amps(iint)%n_sing(1)
  write (*,'(a,1x,es24.16)') &
       'AMPICOL_SELECTED_FLOW_PROBE_NORMALIZATION',norm_factor
  write (*,'(a)') repeat('-',78)
  write (*,'(a)') 'Timing summary                           seconds    percent  note'
  write (*,'(a)') repeat('-',78)
  call print_row('amplitude evaluation',t_eval,t_total,'direct-library')
  call print_row('squaring amplitudes',t_square,t_total,'direct-library')
  call print_row('total',t_total,t_total,'')
  write (*,'(a)') repeat('-',78)

contains

  subroutine square_amplitudes(ichan,jint,local_amps,local_amp2,local_amp2_hel)
    implicit none
    integer,intent(in) :: ichan,jint
    complex(kind=8),dimension(:),intent(in) :: local_amps
    real(kind=8),dimension(:),intent(inout) :: local_amp2,local_amp2_hel
    integer :: local_ih, local_iproc
    local_iproc = 0
    local_amp2 = 0d0
    if (keep_processes_separate) then
       do local_ih=1,pgl(ichan)%amps(jint)%n_amps
          do while (pgl(ichan)%amps(jint)%iproc_start(local_iproc+1).eq.local_ih)
             local_iproc=local_iproc+1
          enddo
          local_amp2_hel(local_ih)=dble(local_amps(local_ih)*&
               pgl(ichan)%col_fac(jint)*dconjg(local_amps(local_ih)))*&
               pgl(ichan)%hel_fac(local_ih,jint)
          local_amp2(local_iproc)=local_amp2(local_iproc)+local_amp2_hel(local_ih)
       enddo
    else
       do local_ih=1,pgl(ichan)%amps(jint)%n_amps
          do while (pgl(ichan)%amps(jint)%iproc_start(local_iproc+1).eq.local_ih)
             local_iproc=local_iproc+1
          enddo
          local_amp2_hel(local_ih)=dble(local_amps(local_ih)*&
               pgl(ichan)%col_fac(local_iproc)*dconjg(local_amps(local_ih)))*&
               pgl(ichan)%hel_fac(local_ih,jint)
          local_amp2(local_iproc)=local_amp2(local_iproc)+local_amp2_hel(local_ih)
       enddo
    endif
  end subroutine square_amplitudes

  subroutine print_row(label,time,total,note)
    implicit none
    character(len=*),intent(in) :: label,note
    real(kind=8),intent(in) :: time,total
    real(kind=8) :: pct
    character(len=32) :: label_fmt
    if (total.gt.0d0) then
       pct=100d0*time/total
    else
       pct=0d0
    endif
    label_fmt=adjustl(label)
    write(*,'(a32,2x,f14.6,2x,f8.2,a,2x,a)') label_fmt,time,pct,'%',trim(note)
  end subroutine print_row

end program amplicol_library_benchmark
