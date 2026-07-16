program amplicol_color_library_probe
  use common
  use particles
  use amplitude_QCD_mod
  use amplitude_library, only: read_amplitude_lib
  use amp_lib, only: evaluate_amp
  use handling_processes
  implicit none

  type(amplitude_QCD) :: colour_amp
  integer(kind=8) :: points, ipoint
  integer :: igroup, iint, argc, n, nmax, iacc_request, col_acc
  integer :: n_helicity_combinations
  integer :: irow, i, ic, icol, ioff, max_amp_size
  integer,dimension(:),allocatable :: hel, row_to_integral, lepton_list
  integer,dimension(:,:),allocatable :: local_part, local_order, spin_init, spin_loop
  integer,dimension(:,:),allocatable :: helicity_table, helicity_amp_index
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(3) :: matrix2, matrix2_result
  real(kind=8) :: t0, t1, t_eval, t_colour, t_total, norm_factor
  complex(kind=8),dimension(:),allocatable :: colour_amps
  complex(kind=8),dimension(:,:),allocatable :: order_amps
  complex(kind=8) :: amp_col_c, amp2_c
  character(len=80) :: momenta_file
  character(len=32) :: color_arg, color_name
  character(len=256) :: arg, env_value
  logical :: print_matrix
  real(kind=8),parameter :: alpha_check=0.118d0
  real(kind=8),parameter :: pi=3.14159265358979323846d0

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
     read(arg,'(a)') momenta_file
  endif
  if (points < 1_8) then
     write (*,*) 'points must be positive'
     stop 1
  endif

  call parse_color_accuracy()

  open(unit=99,file='amplicol_color_library_probe.output',status='unknown')
  call phys_model%init_part(173d0,1.491500d0,91.188d0,2.441404d0,&
       80.419002445756163d0,2.0476d0,125d0,0.0063823389999999999d0)
  call phys_model%init_vert()
  call read_amplitude_lib()
  if (igroup < 1 .or. igroup > ngroups) then
     write (*,*) 'invalid group',igroup,'ngroups=',ngroups
     stop 1
  endif
  if (iint < 1 .or. iint > size(pgl(igroup)%amps)) then
     write (*,*) 'invalid integral',iint,'n_integrals=',size(pgl(igroup)%amps)
     stop 1
  endif

  n = pgl(igroup)%next
  nmax = maxval(pgl(:)%next)
  allocate(p(0:3,nmax))
  allocate(hel(n))
  allocate(local_part(n,1))
  allocate(local_order(n,1))
  allocate(spin_init(0:3,n))
  allocate(spin_loop(0:3,n))

  call read_probe_momenta()

  local_part(1:n,1) = pgl(igroup)%processes(1:n,iint)
  local_order(1:n,1) = pgl(igroup)%color_orders(1:n,iint)
  spin_init = 0
  spin_init(0,1:n) = 1
  spin_init(1,1:n) = -9
  call setup_local_spin_loop()
  call build_lepton_list(local_part(1:n,1))

  call colour_amp%init(2,n,1,local_part,spin_init,local_order,phys_model,&
       lepton_list(1),lepton_list)
  call colour_amp%init_col(n,col_acc)
  if (print_matrix) call print_color_matrix()
  call build_row_to_integral()
  call build_helicity_lookup()

  max_amp_size = 0
  do i=1,size(pgl(igroup)%amps)
     max_amp_size = max(max_amp_size,size(pgl(igroup)%amps(i)%amps))
  enddo
  allocate(order_amps(max_amp_size,colour_amp%nColOrd))
  allocate(colour_amps(colour_amp%nColOrd))

  norm_factor = (4*pi*alpha_check)**(n-2-colour_amp%n_sing(1))&
       *(2d0*4d0*pi*alphaEW)**colour_amp%n_sing(1)/dble(pgl(igroup)%iden(iint))

  matrix2 = 0d0
  call cpu_time(t0)
  do ipoint=1,points
     call evaluate_colour_order_amplitudes()
     call contract_all_helicities()
  enddo
  call cpu_time(t1)
  t_total = t1 - t0
  matrix2_result = matrix2

  ! Time the two components in separate outer loops.  Per-point and
  ! per-helicity CPU_TIME calls dominate low-multiplicity processes.
  call cpu_time(t0)
  do ipoint=1,points
     call evaluate_colour_order_amplitudes()
  enddo
  call cpu_time(t1)
  t_eval = t1 - t0

  matrix2 = 0d0
  call cpu_time(t0)
  do ipoint=1,points
     call contract_all_helicities()
  enddo
  call cpu_time(t1)
  t_colour = t1 - t0
  matrix2 = matrix2_result

  write (*,'(a)') 'AmpliCol generated-library colour probe'
  write (*,'(a,i0)') 'points ',points
  write (*,'(a,i0)') 'group ',igroup
  write (*,'(a,i0)') 'integral ',iint
  write (*,'(a,a)') 'color_accuracy ',trim(color_name)
  write (*,'(a,i0)') 'color_orders ',colour_amp%nColOrd
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
  call print_row('amplitude evaluation',t_eval,t_total,'outer-loop-diagnostic')
  call print_row('colour contraction',t_colour,t_total,'outer-loop-diagnostic')
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
       default_file = 'Library/amp'//trim(adjustl(tmp))//'_'
       write(tmp,*) iint
       default_file = trim(adjustl(default_file))//trim(adjustl(tmp))//'_lib.data'
       p = 0d0
       open(unit=14,file=trim(adjustl(default_file)),form='unformatted',&
            access='stream',status='old')
       read(14) p
       close(14)
       return
    else
       default_file = momenta_file
    endif
    p = 0d0
    open(unit=14,file=trim(adjustl(default_file)),status='old')
    do j=1,n
       read(14,*) p(0,j),p(1,j),p(2,j),p(3,j)
    enddo
    close(14)
  end subroutine read_probe_momenta

  subroutine read_legacy_probe_momenta()
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
    p = 0d0
    open(unit=14,file=trim(adjustl(default_file)),status='old')
    do j=1,n
       read(14,*) p(0,j),p(1,j),p(2,j),p(3,j)
    enddo
    close(14)
  end subroutine read_legacy_probe_momenta

  subroutine build_lepton_list(process)
    implicit none
    integer,dimension(n),intent(in) :: process
    integer :: j, nl, nal, nlep
    nlep = 0
    do j=1,n
       if (phys_model%is_lepton(process(j)) .or. &
           phys_model%is_antilepton(process(j))) nlep = nlep + 1
    enddo
    allocate(lepton_list(1+2*nlep))
    lepton_list = 0
    lepton_list(1) = 2*nlep
    nl = 0
    nal = 0
    do j=1,n
       if (phys_model%is_lepton(process(j))) then
          nl = nl + 1
          lepton_list(nl+2) = j
       elseif (phys_model%is_antilepton(process(j))) then
          nal = nal + 1
          lepton_list(nal+2) = -j
       endif
    enddo
  end subroutine build_lepton_list

  subroutine setup_local_spin_loop()
    implicit none
    integer :: j
    spin_loop = 0
    do j=1,n
       spin_loop(0,j) = phys_model%get_spin(local_part(j,1))
       if (spin_loop(0,j).eq.2) then
          spin_loop(1,j) = -1
          spin_loop(2,j) = 1
       elseif (spin_loop(0,j).eq.3) then
          spin_loop(1,j) = -1
          spin_loop(2,j) = 0
          spin_loop(3,j) = 1
       elseif (spin_loop(0,j).eq.1) then
          spin_loop(1,j) = 0
       else
          write (*,*) 'spin state not known',j,local_part(j,1),spin_loop(0,j)
          stop 1
       endif
    enddo
  end subroutine setup_local_spin_loop

  subroutine build_helicity_lookup()
    implicit none
    integer :: combination, position, state, remaining, row, jint, active
    n_helicity_combinations = 1
    do position=1,n
       n_helicity_combinations = n_helicity_combinations*spin_loop(0,position)
    enddo
    allocate(helicity_table(n,n_helicity_combinations))
    allocate(helicity_amp_index(colour_amp%nColOrd,n_helicity_combinations))
    helicity_amp_index = 0
    do combination=1,n_helicity_combinations
       remaining = combination - 1
       do position=1,n
          state = mod(remaining,spin_loop(0,position)) + 1
          remaining = remaining/spin_loop(0,position)
          helicity_table(position,combination) = spin_loop(state,position)
       enddo
       hel(1:n) = helicity_table(1:n,combination)
       do row=1,colour_amp%nColOrd
          jint = row_to_integral(row)
          helicity_amp_index(row,combination) = find_helicity_index(jint)
       enddo
    enddo
    active = 0
    do combination=1,n_helicity_combinations
       if (.not.any(helicity_amp_index(:,combination).gt.0)) cycle
       active = active + 1
       if (active.ne.combination) then
          helicity_table(:,active) = helicity_table(:,combination)
          helicity_amp_index(:,active) = helicity_amp_index(:,combination)
       endif
    enddo
    n_helicity_combinations = active
  end subroutine build_helicity_lookup

  subroutine build_row_to_integral()
    implicit none
    integer :: row, jint
    allocate(row_to_integral(colour_amp%nColOrd))
    row_to_integral = 0
    do row=1,colour_amp%nColOrd
       do jint=1,size(pgl(igroup)%amps)
          if (.not.all(pgl(igroup)%processes(1:n,jint).eq.local_part(1:n,1))) cycle
          if (colour_order_matches(jint,row)) then
             row_to_integral(row) = jint
             exit
          endif
       enddo
       if (row_to_integral(row).eq.0) then
          write (*,*) 'Could not map colour row to generated-library integral',row
          write (*,*) 'row permutation:',colour_amp%perm(1:n-colour_amp%n_sing(1),row)
          stop 1
       endif
    enddo
  end subroutine build_row_to_integral

  logical function colour_order_matches(jint,row)
    implicit none
    integer,intent(in) :: jint,row
    integer :: pos, label, m, nord
    integer,dimension(n) :: candidate
    nord = n - colour_amp%n_sing(1)
    candidate = 0
    m = 0
    do pos=1,n
       label = pgl(igroup)%color_orders(pos,jint)
       if (label.lt.1 .or. label.gt.n) cycle
       if (.not.phys_model%is_singlet(pgl(igroup)%processes(label,jint))) then
          m = m + 1
          candidate(m) = label
       endif
    enddo
    colour_order_matches = .false.
    if (m.ne.nord) return
    if (all(candidate(1:nord).eq.colour_amp%perm(1:nord,row))) then
       colour_order_matches = .true.
       return
    endif
    if (is_pure_gluon_word(candidate,nord)) then
       colour_order_matches = cyclic_order_matches(candidate,colour_amp%perm(1:nord,row),nord)
    endif
  end function colour_order_matches

  logical function is_pure_gluon_word(word,nord)
    implicit none
    integer,intent(in) :: nord
    integer,dimension(n),intent(in) :: word
    integer :: j
    is_pure_gluon_word = .true.
    do j=1,nord
       if (.not.phys_model%is_gluon(local_part(word(j),1))) then
          is_pure_gluon_word = .false.
          return
       endif
    enddo
  end function is_pure_gluon_word

  logical function cyclic_order_matches(candidate,row_perm,nord)
    implicit none
    integer,intent(in) :: nord
    integer,dimension(n),intent(in) :: candidate
    integer,dimension(:),intent(in) :: row_perm
    integer :: shift, j
    cyclic_order_matches = .false.
    do shift=0,nord-1
       cyclic_order_matches = .true.
       do j=1,nord
          if (candidate(mod(j-1+shift,nord)+1).ne.row_perm(j)) then
             cyclic_order_matches = .false.
             exit
          endif
       enddo
       if (cyclic_order_matches) return
    enddo
  end function cyclic_order_matches

  subroutine contract_all_helicities()
    implicit none
    integer :: combination, row, ih
    do combination=1,n_helicity_combinations
       colour_amps = (0d0,0d0)
       do row=1,colour_amp%nColOrd
          ih = helicity_amp_index(row,combination)
          if (ih.gt.0) colour_amps(row) = order_amps(ih,row)
       enddo
       call accumulate_colour()
    enddo
  end subroutine contract_all_helicities

  subroutine evaluate_colour_order_amplitudes()
    implicit none
    integer :: row, jint
    order_amps = (0d0,0d0)
    do row=1,colour_amp%nColOrd
       jint = row_to_integral(row)
       call evaluate_amp(igroup,jint,p,order_amps(:,row))
    enddo
  end subroutine evaluate_colour_order_amplitudes

  integer function find_helicity_index(jint)
    implicit none
    integer,intent(in) :: jint
    integer :: ih, ispin
    find_helicity_index = 0
    do ih=1,pgl(igroup)%amps(jint)%n_amps
       do ispin=1,size(pgl(igroup)%amps(jint)%spins,2)
          if (all(pgl(igroup)%amps(jint)%spins(1:n,ispin,ih).eq.hel(1:n))) then
             find_helicity_index = ih
             return
          endif
       enddo
    enddo
  end function find_helicity_index

  subroutine accumulate_colour()
    implicit none
    integer :: iacc_loop
    ioff = colour_amp%iproc_start(colour_amp%nprocs) - 1
    do iacc_loop=1,3
       if (iacc_loop.eq.2 .and. col_acc.lt.1) cycle
       if (iacc_loop.eq.3 .and. col_acc.lt.2) cycle
       do irow=1,colour_amp%nColOrd
          amp_col_c = (0d0,0d0)
          do i=1,colour_amp%n_col_vals(iacc_loop)
             amp2_c = (0d0,0d0)
             do ic=colour_amp%row_index(irow-1,i,iacc_loop)+1,&
                   colour_amp%row_index(irow,i,iacc_loop)
                icol = colour_amp%col_index(colour_amp%i_col_i(i,iacc_loop)+ic)
                amp2_c = amp2_c + colour_amps(icol)
             enddo
             amp_col_c = amp_col_c + amp2_c*colour_amp%diff_col_vals(i,iacc_loop)
          enddo
          matrix2(iacc_loop) = matrix2(iacc_loop) + &
               dble(amp_col_c*conjg(colour_amps(irow)))
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
    write (*,'(a,i0)') 'AMPICOL_COLOR_MATRIX_NCOLOR_ORDERS ',colour_amp%nColOrd
    do row=1,colour_amp%nColOrd
       write (*,'(a,i0,1x,*(i0,1x))') 'AMPICOL_COLOR_MATRIX_PERM ',row,&
            (colour_amp%perm(word_index,row),word_index=1,n-colour_amp%n_sing(1))
    enddo
    do iacc_loop=1,3
       if (iacc_loop.eq.2 .and. col_acc.lt.1) cycle
       if (iacc_loop.eq.3 .and. col_acc.lt.2) cycle
       do row=1,colour_amp%nColOrd
          do val=1,colour_amp%n_col_vals(iacc_loop)
             first_pos = colour_amp%row_index(row-1,val,iacc_loop)+1
             last_pos = colour_amp%row_index(row,val,iacc_loop)
             do pos=first_pos,last_pos
                col = colour_amp%col_index(colour_amp%i_col_i(val,iacc_loop)+pos)
                write (*,'(a,1x,a,1x,i0,1x,i0,1x,es24.16)')&
                     'AMPICOL_COLOR_MATRIX_ENTRY',trim(labels(iacc_loop)),&
                     row,col,colour_amp%diff_col_vals(val,iacc_loop)
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

end program amplicol_color_library_probe
