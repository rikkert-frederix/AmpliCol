! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module rw_events
  implicit none
  real(kind=8) :: wgt,evt_wgt,weight,amp2,rwgt_NLC,rwgt_full
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
  type(amplitude_QCD) :: amps
  type(physics_model) :: phys_model
  integer,parameter :: string_len=150
  integer :: i,j,col_acc,icol,irow,ic,iacc,nColOrd,next
  integer,dimension(:),allocatable :: hel
  integer,dimension(:,:),allocatable :: spin,o,part
  real(kind=8) :: amp2,amp_col
  real(kind=8),dimension(3) :: matrix2
  real(kind=8),dimension(:,:),allocatable :: p
  complex(kind=8) :: amp2_c,amp_col_c
  logical :: done
  character(len=string_len) :: tag,tag_read,add_arg=''
  
  call get_run_arguments()

  call cpu_time(tTot_B)
  
  call phys_model%init_part(173d0,1.491500d0)
  call setup_spin()

  ! read one event to determine 'next' and allocate the required arrays.
  call read_event(11,done)
  rewind(11)
  

  call cpu_time(tBefore)

  call amps%init(2,next,1,part,spin,o,phys_model)
  col_acc=20
  call amps%init_col(next,col_acc)
  if (amps%n_qqbar(1).eq.2 .and. amps%same_flav(1)) then
     nColOrd=amps%n_amps/2
  else
     nColOrd=amps%n_amps
  endif

  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  do
     call read_event(11,done)
     if (done) exit
     matrix2(1:3)=0d0

     call cpu_time(tBefore)

     call amps%evaluate(next,p,hel)

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore

     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amps%n_qqbar(1).eq.0 .and. use_real_gluons) then
           ! same as in the 'else' below, except that all are real variables instead of complex. 
           do irow=1,nColOrd
              amp_col=0d0
              do i=1,amps%n_col_vals(iacc)
                 amp2=0d0
                 do ic=amps%row_index(irow-1,i,iacc)+1,amps%row_index(irow,i,iacc)
                    icol=amps%col_index(amps%i_col_i(i,iacc)+ic)
                    amp2=amp2+amps%amps_r(icol)
                 enddo
                 amp_col=amp_col+amp2*amps%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+amp_col*amps%amps_r(irow)
           enddo
        else
           do irow=1,nColOrd
              amp_col_c=(0d0,0d0)
              do i=1,amps%n_col_vals(iacc)
                 amp2_c=(0d0,0d0)
                 do ic=amps%row_index(irow-1,i,iacc)+1,amps%row_index(irow,i,iacc)
                    icol=amps%col_index(amps%i_col_i(i,iacc)+ic)
                    amp2_c=amp2_c+amps%amps(icol)
                 enddo
                 amp_col_c=amp_col_c+amp2_c*amps%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps%amps(irow)))
           enddo
        endif
        
        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo

     write (*,*) 'matrix2:', matrix2
     stop 1
     
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
       open(unit=12,file=trim(adjustl(filename))//'.rwgt',status='unknown')
    elseif (argc.le.8) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'next, *process*, *order*'
       write(*,*) '--------- or: ---------------'
       write(*,*) 'event_file_name_to_reweight'
       stop 1
    else
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) then
             read(argv,*) next
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             allocate(part(1:next,1))
             allocate(o(1:next,1))
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
    if (.not. allocated(hel)) then
       allocate(o(next,1))
       allocate(part(next,1))
       allocate(hel(next))
       allocate(p(0:3,next))
    endif
       read (iunit,*,err=99,end=99) hel(1:next)
    read (iunit,*,err=99,end=99) o(1:next,1)
    do i=1,next
       read (iunit,*,err=99,end=99) part(i,1),p(1:3,i),p(0,i)
    enddo
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event


  subroutine write_event(iunit)
    use rw_events
    implicit none
    integer :: i,iunit
    rwgt_NLC=matrix2(2)/matrix2(1)
    rwgt_full=matrix2(3)/matrix2(1)
    write (iunit,*) '<event>'
    write (iunit,*) next,evt_wgt*rwgt_full!,wgt,matrix2,weight
    write (iunit,'(100i3)') hel(1:next)
    write (iunit,'(100i3)') o(1:next,1)
    write (iunit,*) rwgt_full,rwgt_NLC!,matrix2(1),matrix2(2),matrix2(3)
    write (iunit,*) evt_wgt,evt_wgt*rwgt_NLC,evt_wgt*rwgt_full
    do i=1,next
       write (iunit,*) part(i,1),p(1:3,i),p(0,i)
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event



end program matrix_reweight
