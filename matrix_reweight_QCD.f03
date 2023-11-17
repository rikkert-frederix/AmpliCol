! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module common
  use amplitude_QCD_mod
  implicit none
  integer :: next
  type(amplitude_QCD) :: amp_QCD
  real(kind=8),dimension(:,:),allocatable :: p
end module common
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
  integer :: c_o,imode
end module arguments

program matrix_reweight
  use math_functions
  use common
  use timings
  implicit none
  integer :: i,j,col_acc,icol,ihel,hel_picked,irow,ic,iacc
  integer,dimension(:),allocatable :: hel,o,part
  real(kind=8),dimension(3) :: matrix2
  real(kind=8) :: color_wgt,amp,amp2,amp_col
  complex(kind=8) :: amp2_c,amp_col_c
  real(kind=8),dimension(:),allocatable :: mass
  logical :: done
  
  call get_run_arguments()

  call cpu_time(tTot_B)

  allocate(o(next))
  allocate(hel(next))
  allocate(mass(next))
  allocate(p(0:3,next))

  mass(1:next)=0d0
  call create_run_tag_and_open_files()

  call cpu_time(tBefore)

  if (.not.allocated(part)) allocate(part(1:next))
  call read_event(11,done)
  rewind(11)
  call amp_QCD%init(2,next,part,o)
  col_acc=20
  call amp_QCD%init_col2(next,part,o,col_acc)

  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  do 
     call read_event(11,done)
     if (done) exit
     matrix2(1:3)=0d0

     call cpu_time(tBefore)
     ! read helicity from event file
     ihel=hel_picked
     do i=1,next
        if (btest(ihel-1,i-1)) then
           hel(i)=1
        else
           hel(i)=0
        endif
     enddo
     call amp_QCD%evaluate(next,p,ihel)
     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore


     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amp_QCD%n_qqbar.eq.0) then
           do irow=1,factorial(next-1)
              if (use_real_gluons) then
                 amp_col=0d0
              else
                 amp_col_c=(0d0,0d0)
              endif
              do i=1,amp_QCD%n_col_vals(iacc)
                 if (use_real_gluons) then
                    amp2=0d0
                 else
                    amp2_c=(0d0,0d0)
                 endif
                 do ic=amp_QCD%row_index(irow-1,i,iacc)+1,amp_QCD%row_index(irow,i,iacc)
                    icol=amp_QCD%col_index(ic,i,iacc)
                    if (use_real_gluons) then
                       amp2=amp2+amp_QCD%amps_r(icol)
                    else
                       amp2_c=amp2_c+amp_QCD%amps(icol)
                    endif
                 enddo
                 if (use_real_gluons) then
                    amp_col=amp_col+amp2*amp_QCD%diff_col_vals(i,iacc)
                 else
                    amp_col_c=amp_col_c+amp2_c*amp_QCD%diff_col_vals(i,iacc)
                 endif
              enddo
              if (use_real_gluons) then
                 matrix2(iacc)=matrix2(iacc)+amp_col*amp_QCD%amps_r(irow)
              else
                 matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amp_QCD%amps(irow)))
              endif
           enddo
        else
           do irow=1,factorial(next-2)
              amp_col_c=(0d0,0d0)
              do i=1,amp_QCD%n_col_vals(iacc)
                 amp2_c=(0d0,0d0)
                 do ic=amp_QCD%row_index(irow-1,i,iacc)+1,amp_QCD%row_index(irow,i,iacc)
                    icol=amp_QCD%col_index(ic,i,iacc)
                    amp2_c=amp2_c+amp_QCD%amps(icol)
                 enddo
                 amp_col_c=amp_col_c+amp2_c*amp_QCD%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amp_QCD%amps(irow)))
           enddo
        endif
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
    integer :: argc
    character(len=256) :: argv
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    if (argc.ne.3) then
       write (*,*) 'Give number of gluons, imode and color'// &
            ' ordering (number of gluons on first color line):'
       read (*,*) next,imode,c_o
    else
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) read(argv,*) next
          if (i.eq.2) read(argv,*) imode
          if (i.eq.3) read(argv,*) c_o
       enddo
    endif
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif
    if (imode.ne.2) then
       write (*,*) 'Incorrect imode',imode, ' (should be 2)'
       stop
    endif
    if (c_o.lt.0 .or. c_o .gt. next-2) then
       write (*,*) 'inconsistent color-ordering',c_o
       stop
    endif
  end subroutine get_run_arguments

  subroutine create_run_tag_and_open_files()
    use arguments
    implicit none
    character(len=1) :: s1
    character(len=2) :: s2
    character(len=8) :: tag,tag_read
    if (next.le.9) then
       write(s1,'(i1)') next
       tag=trim(adjustl(s1))//'_'
       tag_read=trim(adjustl(s1))//'_'
    else
       write(s2,'(i2)') next
       tag=trim(adjustl(s2))//'_'
       tag_read=trim(adjustl(s2))//'_'
    endif
    write(s1,'(i1)') imode
    tag=trim(adjustl(tag))//trim(adjustl(s1))//'_'
    if (imode.gt.0) write(s1,'(i1)') imode-1
    tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))//'_'
    if (c_o.le.9) then
       write(s1,'(i1)') c_o
       tag=trim(adjustl(tag))//trim(adjustl(s1))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s1))
    else
       write(s2,'(i2)') c_o
       tag=trim(adjustl(tag))//trim(adjustl(s2))
       tag_read=trim(adjustl(tag_read))//trim(adjustl(s2))
    endif
    if (len(trim(tag_read)).lt.8) then
       if (8-len(trim(tag)).eq.1) then
          tag='_'//trim(adjustl(tag))
          tag_read='_'//trim(adjustl(tag_read))
       elseif(8-len(trim(tag)).eq.2) then
          tag='__'//trim(adjustl(tag))
          tag_read='__'//trim(adjustl(tag_read))
       elseif(8-len(trim(tag)).eq.3) then
          tag='___'//trim(adjustl(tag))
          tag_read='___'//trim(adjustl(tag_read))
       endif
    endif
    open(unit=11,file='Outputs/events'//tag//'.lhe',status='old')
    open(unit=12,file='Outputs/events'//tag//'.lhe.rwgt',status='unknown')
  end subroutine create_run_tag_and_open_files

  subroutine read_event(iunit,done)
    use rw_events
    implicit none
    integer :: i,iunit
    logical :: done
    character :: dummy
    real(kind=8) :: dum
    done=.false.
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) dum,hel_picked,evt_wgt,wgt,amp2,weight
    read (iunit,*,err=99,end=99) o(1:next)
    do i=1,next
       read (iunit,*,err=99,end=99) part(i),p(1:3,i),p(0,i)
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
    write (iunit,*) next,evt_wgt,wgt,matrix2,weight
    write (iunit,'(100i3)') o(1:next)
    write (iunit,*) rwgt_full,rwgt_NLC,matrix2(1),matrix2(2),matrix2(3)
    write (iunit,*) evt_wgt,evt_wgt*rwgt_NLC,evt_wgt*rwgt_full
    write (iunit,'(100i3)') hel(1:next)
    do i=1,next
       if (i.le.2) then
          write (iunit,*) part(i),p(1:3,i),p(0,i)
       else
          write (iunit,*) part(i),p(1:3,i),p(0,i)
       endif
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event


end program matrix_reweight
