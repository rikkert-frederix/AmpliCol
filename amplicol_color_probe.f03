program amplicol_color_probe
  use common
  use particles
  use amplitude_QCD_mod
  use read_process_file
  use handling_processes
  implicit none

  type(amplitude_QCD) :: amp
  integer(kind=8) :: points, ipoint
  integer :: igroup, iint, argc, n, iacc_request, col_acc
  integer :: irow, i, ic, icol, ioff
  integer,dimension(:),allocatable :: hel
  integer,dimension(:,:),allocatable :: local_part, local_order, spin_init
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(3) :: matrix2
  real(kind=8) :: t0, t1, t_eval, t_colour, t_total, norm_factor
  complex(kind=8) :: amp_col_c, amp2_c
  character(len=80) :: process_file, momenta_file
  character(len=32) :: color_arg, color_name
  character(len=256) :: arg, env_value
  logical :: print_matrix
  real(kind=8),parameter :: alpha_check=0.118d0
  real(kind=8),parameter :: pi=3.14159265358979323846d0

  process_file = 'processes.txt'
  momenta_file = ''
  color_arg = 'lc'
  color_name = 'lc'
  points = 1_8
  igroup = 1
  iint = 1
  print_matrix = .false.
  call get_environment_variable('AMPICOL_COLOR_PROBE_MATRIX',env_value)
  if (trim(adjustl(env_value)).eq.'1') print_matrix = .true.

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
     read(arg,'(a)') color_arg
  endif
  if (argc >= 5) then
     call get_command_argument(5,arg)
     read(arg,'(a)') process_file
  endif
  if (argc >= 6) then
     call get_command_argument(6,arg)
     read(arg,'(a)') momenta_file
  endif
  if (points < 1_8) then
     write (*,*) 'points must be positive'
     stop 1
  endif

  call parse_color_accuracy()

  open(unit=99,file='amplicol_color_probe.output',status='unknown')
  call phys_model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call phys_model%init_vert()
  call read_processes_from_file(process_file)
  if (igroup < 1 .or. igroup > ngroups) then
     write (*,*) 'invalid group',igroup,'ngroups=',ngroups
     stop 1
  endif
  if (iint < 1 .or. iint > pgl(igroup)%nproc) then
     write (*,*) 'invalid integral',iint,'n_integrals=',pgl(igroup)%nproc
     stop 1
  endif

  n = pgl(igroup)%next
  allocate(p(0:3,n))
  allocate(hel(n))
  allocate(local_part(n,1))
  allocate(local_order(n,1))
  allocate(spin_init(0:3,n))

  call setup_spin(pgl(igroup))
  if (.not.allocated(pgl(igroup)%iden)) allocate(pgl(igroup)%iden(pgl(igroup)%nproc))
  pgl(igroup)%iden(1:pgl(igroup)%nproc) = 1
  call set_final_state_identical_particle_factor(pgl(igroup))
  call set_initial_state_average_factor(pgl(igroup))
  call read_probe_momenta()

  local_part(1:n,1) = pgl(igroup)%processes(1:n,iint)
  local_order(1:n,1) = pgl(igroup)%color_orders(1:n,iint)
  spin_init = 0
  spin_init(0,1:n) = 1
  spin_init(1,1:n) = -9

  call amp%init(2,n,1,local_part,spin_init,local_order,phys_model,&
       pgl(igroup)%lepton_list(1),pgl(igroup)%lepton_list)
  call amp%init_col(n,col_acc)
  if (print_matrix) call print_color_matrix()
  norm_factor = (4*pi*alpha_check)**(n-2-amp%n_sing(1))&
       *(2d0*4d0*pi*alphaEW)**amp%n_sing(1)/dble(pgl(igroup)%iden(iint))

  matrix2 = 0d0
  t_eval = 0d0
  t_colour = 0d0
  call cpu_time(t0)
  do ipoint=1,points
     call loop_helicity(1)
  enddo
  call cpu_time(t1)
  t_total = t1 - t0

  write (*,'(a)') 'AmpliCol colour probe'
  write (*,'(a,i0)') 'points ',points
  write (*,'(a,i0)') 'group ',igroup
  write (*,'(a,i0)') 'integral ',iint
  write (*,'(a,a)') 'color_accuracy ',trim(color_name)
  write (*,'(a,i0)') 'color_orders ',amp%nColOrd
  write (*,'(a,3(1x,es24.16))') 'AMPICOL_COLOR_PROBE_COMPONENTS',&
       norm_factor*matrix2(1)/dble(points),&
       norm_factor*matrix2(2)/dble(points),&
       norm_factor*matrix2(3)/dble(points)
  write (*,'(a,3(1x,es24.16))') 'AMPICOL_COLOR_PROBE_RAW_COMPONENTS',&
       matrix2(1)/dble(points),matrix2(2)/dble(points),matrix2(3)/dble(points)
  write (*,'(a,1x,a,1x,i0,1x,i0,1x,es24.16)') 'AMPICOL_COLOR_PROBE_VALUE',&
       trim(color_name),igroup,iint,norm_factor*matrix2(iacc_request)/dble(points)
  write (*,'(a)') repeat('-',78)
  write (*,'(a)') 'Timing summary                           seconds    percent  note'
  write (*,'(a)') repeat('-',78)
  call print_row('amplitude evaluation',t_eval,t_total,'colour-probe')
  call print_row('colour contraction',t_colour,t_total,'colour-probe')
  call print_row('total',t_total,t_total,'')
  write (*,'(a)') repeat('-',78)

contains

  subroutine parse_color_accuracy()
    implicit none
    character(len=32) :: lowered
    lowered = trim(adjustl(color_arg))
    call lowercase(lowered)
    select case (trim(lowered))
    case ('lc','leading','leading-colour','leading-color','1')
       color_name = 'lc'
       iacc_request = 1
       col_acc = 0
    case ('nlc','next-to-leading-colour','next-to-leading-color','2')
       color_name = 'nlc'
       iacc_request = 2
       col_acc = 1
    case ('full','full-colour','full-color','3')
       color_name = 'full'
       iacc_request = 3
       col_acc = 20
    case default
       write (*,*) 'unknown color accuracy: ',trim(color_arg)
       stop 1
    end select
  end subroutine parse_color_accuracy

  subroutine lowercase(text)
    implicit none
    character(len=*),intent(inout) :: text
    integer :: j, code
    do j=1,len(text)
       code = iachar(text(j:j))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          text(j:j) = achar(code + iachar('a') - iachar('A'))
       endif
    enddo
  end subroutine lowercase

  subroutine read_probe_momenta()
    implicit none
    integer :: j
    character(len=80) :: tmp, default_file
    if (len_trim(momenta_file) == 0) then
       write(tmp,*) igroup
       default_file = 'Utilities/ME_checks/momenta_'//trim(adjustl(tmp))//'_'
       write(tmp,*) iint
       default_file = trim(adjustl(default_file))//trim(adjustl(tmp))//'.txt'
    else
       default_file = momenta_file
    endif
    open(unit=14,file=trim(adjustl(default_file)),status='old')
    do j=1,n
       read(14,*) p(0,j),p(1,j),p(2,j),p(3,j)
    enddo
    close(14)
  end subroutine read_probe_momenta

  recursive subroutine loop_helicity(position)
    implicit none
    integer,intent(in) :: position
    integer :: spin_index
    real(kind=8) :: before, after
    if (position > n) then
       call cpu_time(before)
       call amp%evaluate(n,p,hel,.true.,phys_model)
       call cpu_time(after)
       t_eval = t_eval + after - before
       call cpu_time(before)
       call accumulate_colour()
       call cpu_time(after)
       t_colour = t_colour + after - before
       return
    endif
    do spin_index=1,pgl(igroup)%spin(0,position)
       hel(position) = pgl(igroup)%spin(spin_index,position)
       call loop_helicity(position+1)
    enddo
  end subroutine loop_helicity

  subroutine accumulate_colour()
    implicit none
    integer :: iacc_loop
    ioff = amp%iproc_start(amp%nprocs) - 1
    do iacc_loop=1,3
       if (iacc_loop.eq.2 .and. col_acc.lt.1) cycle
       if (iacc_loop.eq.3 .and. col_acc.lt.2) cycle
       do irow=1,amp%nColOrd
          amp_col_c = (0d0,0d0)
          do i=1,amp%n_col_vals(iacc_loop)
             amp2_c = (0d0,0d0)
             do ic=amp%row_index(irow-1,i,iacc_loop)+1,amp%row_index(irow,i,iacc_loop)
                icol = amp%col_index(amp%i_col_i(i,iacc_loop)+ic)
                amp2_c = amp2_c + amp%amps(ioff+icol)
             enddo
             amp_col_c = amp_col_c + amp2_c*amp%diff_col_vals(i,iacc_loop)
          enddo
          matrix2(iacc_loop) = matrix2(iacc_loop) + &
               dble(amp_col_c*conjg(amp%amps(ioff+irow)))
       enddo
    enddo
  end subroutine accumulate_colour

  subroutine print_color_matrix()
    implicit none
    integer :: iacc_loop, row, val, pos, col, first_pos, last_pos, word_index
    character(len=16),dimension(3) :: labels
    labels(1) = 'lc'
    labels(2) = 'nlc'
    labels(3) = 'full'
    write (*,'(a)') 'AMPICOL_COLOR_MATRIX_BEGIN'
    write (*,'(a,i0)') 'AMPICOL_COLOR_MATRIX_NCOLOR_ORDERS ',amp%nColOrd
    do row=1,amp%nColOrd
       write (*,'(a,i0,1x,*(i0,1x))') 'AMPICOL_COLOR_MATRIX_PERM ',row,&
            (amp%perm(word_index,row),word_index=1,n-amp%n_sing(1))
    enddo
    do iacc_loop=1,3
       if (iacc_loop.eq.2 .and. col_acc.lt.1) cycle
       if (iacc_loop.eq.3 .and. col_acc.lt.2) cycle
       do row=1,amp%nColOrd
          do val=1,amp%n_col_vals(iacc_loop)
             first_pos = amp%row_index(row-1,val,iacc_loop)+1
             last_pos = amp%row_index(row,val,iacc_loop)
             do pos=first_pos,last_pos
                col = amp%col_index(amp%i_col_i(val,iacc_loop)+pos)
                write (*,'(a,1x,a,1x,i0,1x,i0,1x,es24.16)')&
                     'AMPICOL_COLOR_MATRIX_ENTRY',trim(labels(iacc_loop)),&
                     row,col,amp%diff_col_vals(val,iacc_loop)
             enddo
          enddo
       enddo
    enddo
    write (*,'(a)') 'AMPICOL_COLOR_MATRIX_END'
  end subroutine print_color_matrix

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

end program amplicol_color_probe
