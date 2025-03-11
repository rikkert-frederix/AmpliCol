! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module rw_events
  implicit none
  real(kind=8) :: wgt,evt_wgt,weight,amp2,rwgt_NLC,rwgt_full
  integer,dimension(:),allocatable :: helicity,col_order,iPDG
  real(kind=8),dimension(:,:),allocatable :: momenta
end module rw_events
module timings
  implicit none
  real(kind=4) :: tBefore,tAfter,tTot_A=0.,tTot_B=0.,t_amp=0.,t_amp_init=0.,&
       t_mat_LC=0.,t_mat_NLC=0.,t_mat_full=0.,t_all=0.,t_ran=0.
end module timings
module arguments
  implicit none
  integer :: c_o,c_o_t,c_o_i,c_o_j,c_o_k,imode
end module arguments

program matrix_reweight
  use math_functions
  use amplitude_QCD_mod
  use timings
  use particles
  implicit none
  logical,parameter :: use_only_canonical_form=.true.
  integer,parameter :: max_proc=1280
  type(amplitude_QCD),dimension(max_proc) :: amps
  type(physics_model) :: phys_model
  logical :: read_proc_from_file
  integer,parameter :: string_len=150
  integer :: i,j,col_acc,icol,irow,ic,iacc,nColOrd,next,nprocs,iproc,ioff,unique_nproc
  integer,dimension(:),allocatable :: hel,unique_map
  integer,dimension(:,:),allocatable :: spin,o,part,processes,unique_processes
  real(kind=8) :: amp2,amp_col,process_map_value
  real(kind=8),dimension(3) :: matrix2
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(:),allocatable :: unique_map_value
  complex(kind=8) :: amp2_c,amp_col_c
  logical :: done
  character(len=string_len) :: tag,tag_read,add_arg=''
  
  call get_run_arguments()

  call cpu_time(tTot_B)
  
  call phys_model%init_part(173d0,1.491500d0,91.18800d0,2.0d0)

  nprocs=0
  call setup_spin()
  col_acc=20

  do
     call read_event(11,done)
     if (done) exit
     do iproc=1,nprocs
        if (all(part(1:next,1).eq.processes(1:next,iproc)))exit
     enddo
     if (iproc.eq.nprocs+1) then
        call cpu_time(tBefore)
        nprocs=nprocs+1
        processes(1:next,iproc)=part(1:next,1)
        call amps(iproc)%init(2,next,1,part,spin,o,phys_model,read_proc_from_file)
        call amps(iproc)%init_col(next,col_acc)
        call cpu_time(tAfter)
        t_amp_init=t_amp_init+tAfter-tBefore
     endif
     matrix2(1:3)=0d0

     call cpu_time(tBefore)
     
     call amps(iproc)%evaluate(next,p,hel,read_proc_from_file)

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore

     ioff=amps(iproc)%iproc_start(amps(iproc)%nprocs)-1
     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amps(iproc)%n_qqbar(1).eq.0 .and. use_real_gluons) then
           ! same as in the 'else' below, except that all are real variables instead of complex. 
           do irow=1,amps(iproc)%nColOrd
              amp_col=0d0
              do i=1,amps(iproc)%n_col_vals(iacc)
                 amp2=0d0
                 do ic=amps(iproc)%row_index(irow-1,i,iacc)+1,amps(iproc)%row_index(irow,i,iacc)
                    icol=amps(iproc)%col_index(amps(iproc)%i_col_i(i,iacc)+ic)
                    amp2=amp2+amps(iproc)%amps_r(ioff+icol)
                 enddo
                 amp_col=amp_col+amp2*amps(iproc)%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+amp_col*amps(iproc)%amps_r(ioff+irow)
           enddo
        else
           do irow=1,amps(iproc)%nColOrd
              amp_col_c=(0d0,0d0)
              do i=1,amps(iproc)%n_col_vals(iacc)
                 amp2_c=(0d0,0d0)
                 do ic=amps(iproc)%row_index(irow-1,i,iacc)+1,amps(iproc)%row_index(irow,i,iacc)
                    icol=amps(iproc)%col_index(amps(iproc)%i_col_i(i,iacc)+ic)
                    amp2_c=amp2_c+amps(iproc)%amps(ioff+icol)
                 enddo
                 amp_col_c=amp_col_c+amp2_c*amps(iproc)%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps(iproc)%amps(ioff+irow)))
           enddo
        endif
        ! The following line can be removed since 'process_map_value'
        ! will drop out when taking the ratio w.r.t. LC.
!        matrix2(iacc)=matrix2(iacc)*process_map_value
        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo
     call write_event(12)
  enddo
  
  close(11)
  close(12)
  call cpu_time(tTot_a)
  t_all=tTot_a-tTot_b

  write(*,*) 'Time spent in amplitude initialisation',t_Amp_init
  write(*,*) 'Time spent in amplitude evaluation',t_Amp
  write(*,*) 'Time spent in squaring amplitudes (LC)',t_mat_LC
  write(*,*) 'Time spent in squaring amplitudes (NLC)',t_mat_NLC
  write(*,*) 'Time spent in squaring amplitudes (full)',t_mat_full
  write(*,*) 'Time spent in picking random colors',t_ran
  write(*,*) 'Total time:',t_all
contains  


  subroutine get_run_arguments()
    use arguments
    implicit none
    integer :: argc,i,k
    character(len=256) :: argv,filename
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()

    if (argc.eq.1) then
       call get_command_argument(1,argv)
       read(argv,'(a)') filename
       open(unit=11,file=filename,status='old')
       call read_unique_in_file()
       call allocate_process_info()
       open(unit=12,file=trim(adjustl(filename))//'.rwgt',status='unknown')
    elseif (argc.le.8) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'next, *process*, *order*'
       write(*,*) '--------- or: ---------------'
       write(*,*) 'event_file_name_to_reweight'
       stop 1
    else
       read_proc_from_file=.false.
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) then
             read(argv,*) next
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             unique_nproc=0
             call allocate_process_info()
          endif
          do k=0,next-1
             if (i.eq.2+k) then
                read(argv,*) part(k+1,1)
             endif
          enddo
          do k=0,next-1
             if (i.eq.2+next+k) then
                read(argv,*) o(k+1,1)
             endif
          enddo
          if (argc.eq.1+2*next +1 .and. i.eq.argc) then
             ! Special case: we have an additional argument. Use it as a special tag
             read(argv,*) add_arg
          endif
       enddo
       if (next.lt.4) then
          write (*,*) 'Not enough external particles',next
          stop 1
       endif
       call create_run_tag_and_open_files()
    endif
    
  end subroutine get_run_arguments

  subroutine allocate_process_info()
    use rw_events
    implicit none
    if (.not. allocated(hel)) then
       allocate(momenta(0:3,1:next))
       allocate(helicity(1:next))
       allocate(col_order(1:next))
       allocate(iPDG(1:next))
       if (.not.allocated(o)) allocate(o(next,1))
       if (.not.allocated(part)) allocate(part(next,1))
       if (.not.allocated(processes)) allocate(processes(next,max_proc))
       allocate(hel(next))
       allocate(p(0:3,next))
    endif
  end subroutine allocate_process_info
  
  subroutine read_unique_in_file()
    implicit none
    integer :: iproc
    read(11,*) next,unique_nproc
    allocate(unique_map(unique_nproc))
    allocate(unique_map_value(unique_nproc))
    allocate(unique_processes(next,unique_nproc))
    do iproc=1,unique_nproc
       read(11,*) unique_map(iproc),unique_map_value(iproc),unique_processes(1:next,iproc)
    enddo
  end subroutine read_unique_in_file
  
  subroutine setup_spin()
    implicit none
    if (.not. allocated(spin)) allocate(spin(0:3,1:next))
    do i=1,next
       spin(0,i)=1  ! one arbitrary spin state (use '-9')
       spin(1,i)=-9
    enddo
  end subroutine setup_spin
  
  subroutine create_run_tag_and_open_files()
    use arguments
    implicit none
    tag='_'       ! tag of current run
    tag_read='_'  ! same as 'tag', but with previous imode (i.e., defines the file to read the integration grids from)
    call add_to_string(tag,next,.true.)
    call add_to_string(tag_read,next,.true.)
    call add_to_string(tag,2,.true.)
    call add_to_string(tag_read,2,.true.)
    do i=1,next
       call add_to_string(tag,part(i,1),.true.)
       call add_to_string(tag_read,part(i,1),.true.)
    enddo
    do i=1,next-1
       call add_to_string(tag,o(i,1),.true.)
       call add_to_string(tag_read,o(i,1),.true.)
    enddo
    call add_to_string(tag,o(next,1),.false.)
    call add_to_string(tag_read,o(next,1),.false.)
    open(unit=11,file='Outputs'//trim(adjustl(add_arg))//'/events'//trim(adjustl(tag))//'.lhe',status='old')
    open(unit=12,file='Outputs'//trim(adjustl(add_arg))//'/events'//trim(adjustl(tag))//'.lhe.rwgt',status='unknown')
  end subroutine create_run_tag_and_open_files

  subroutine add_to_string(string,inter,add_underscore)
    ! Adds an integer 'inter' to the end of the string 'string' (followed by
    ! an underscore if 'add_underscore=.true.')
    implicit none
    character(len=string_len) :: string
    integer :: inter
    logical :: add_underscore
    character(len=1) :: s1
    character(len=2) :: s2
    character(len=3) :: s3
    if (inter.ge.0 .and. inter.le.9) then
       write(s1,'(i1)') inter
       string=trim(adjustl(string))//trim(adjustl(s1))
       if (add_underscore) string=trim(adjustl(string))//'_'
    elseif(inter.ge.-9 .and. inter.le.99) then
       write(s2,'(i2)') inter
       string=trim(adjustl(string))//trim(adjustl(s2))
       if (add_underscore) string=trim(adjustl(string))//'_'
    elseif(inter.ge.-99 .and. inter.le.999) then
       write(s3,'(i3)') inter
       string=trim(adjustl(string))//trim(adjustl(s3))
       if (add_underscore) string=trim(adjustl(string))//'_'
    else
       write (*,*) 'value too large to add to the run tag',inter
    endif
  end subroutine add_to_string

  subroutine read_event(iunit,done)
    use rw_events
    implicit none
    integer :: i,iunit
    logical :: done
    character :: dummy
    real(kind=8) :: dum
    done=.false.
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) next,evt_wgt!,wgt,amp2,weight
    read (iunit,*,err=99,end=99) helicity(1:next)
    read (iunit,*,err=99,end=99) col_order(1:next)
    do i=1,next
       read (iunit,*,err=99,end=99) iPDG(i),momenta(1:3,i),momenta(0,i)
    enddo
    call map_to_canonical_form()
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event


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
    ! Initialize mapping
    mapping = [(i,i=1,n)]
    ! Sort the array and mapping using a simple bubble sort
    do i=1,n-1
       do j=1,n-i
          if (array(j) .gt. array(j+1)) then
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

  
  subroutine map_to_canonical_form()
    ! cross the two initial state particle PDGs, order according to
    ! the PDG value, (and reflip the two initial states again)
    use rw_events
    implicit none
    integer,dimension(next) :: mapping
    real(kind=8),dimension(0:3,next) :: p_cross
    integer :: i,iproc
    part(1:next,1)=iPDG(1:next)
    if (.not.use_only_canonical_form) then
       ! do no use the mapping to canonical form, but reweight the
       ! events as they are.
       hel(1:next)=helicity(1:next)
       o(1:next,1)=col_order(1:next)
       p(0:3,1:next)=momenta(0:3,1:next)
    else
       ! Map to canonical from to reduce the number of matrix elements
       ! to initialise.
       ! cross the initial state
       part(1,1)=phys_model%get_antipart(part(1,1))
       part(2,1)=phys_model%get_antipart(part(2,1))
       p_cross(0:3,1:2)=-momenta(0:3,1:2)
       p_cross(0:3,3:next)=momenta(0:3,3:next)
       ! determing the mapping
       call sort_with_mapping(next,part(1,1),mapping)
       ! cross the initial state
       part(1,1)=phys_model%get_antipart(part(1,1))
       part(2,1)=phys_model%get_antipart(part(2,1))
       ! apply the mapping to the momenta and helicity.
       do i=1,next
          if (i.le.2) then
             p(0:3,i)=-p_cross(0:3,mapping(i))
          else
             p(0:3,i)=p_cross(0:3,mapping(i))
          endif
          hel(i)=helicity(mapping(i))
          o(i,1)=col_order(mapping(i)) ! this is not correct, but isn't used
       enddo
       ! Convert to 'unique flavour configuration' (if available)
       if (unique_nproc.eq.0) return
       do iproc=1,unique_nproc
          if (all(part(1:next,1).eq.unique_processes(1:next,iproc))) then
             process_map_value=unique_map_value(iproc)
             if (unique_map(iproc).gt.0) then
                part(1:next,1)=unique_processes(1:next,unique_map(iproc))
             endif
             exit
          endif
       enddo
       if (iproc.eq.unique_nproc+1) then
          write (*,*) 'Process not found among unique processes'
          write (*,*) part(1:next,1)
          stop 1
       endif
    endif
  end subroutine map_to_canonical_form

  subroutine write_event(iunit)
    use rw_events
    implicit none
    integer :: i,iunit
    rwgt_NLC=matrix2(2)/matrix2(1)
    rwgt_full=matrix2(3)/matrix2(1)
    write (iunit,*) '<event>'
    write (iunit,*) next,evt_wgt*rwgt_full!,wgt,matrix2,weight
    write (iunit,'(100i3)') helicity(1:next)
    write (iunit,'(100i3)') col_order(1:next)
    write (iunit,*) rwgt_full,rwgt_NLC!,matrix2(1),matrix2(2),matrix2(3)
    write (iunit,*) evt_wgt,evt_wgt*rwgt_NLC,evt_wgt*rwgt_full
    do i=1,next
       write (iunit,*) iPDG(i),momenta(1:3,i),momenta(0,i)
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event



end program matrix_reweight
