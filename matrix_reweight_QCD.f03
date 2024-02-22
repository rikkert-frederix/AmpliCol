! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module common
  use amplitude_QCD_mod
  implicit none
  integer :: next
  type(amplitude_QCD) :: amp_QCD
  type(amplitude_QCD),dimension(:),allocatable :: amps
  real(kind=8),dimension(:,:),allocatable :: p
  integer :: string_len=50
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
  integer :: c_o,c_o_t,c_o_i,c_o_j,c_o_k,imode
end module arguments

program matrix_reweight
  use math_functions
  use common
  use timings
  implicit none
  integer :: i,j,col_acc,icol,ihel,hel_picked,irow,ic,iacc
  integer :: icol_mat,irow_mat,ri,ri_end,m,proc_num,iacc_in,k,skip 
  integer,dimension(:),allocatable :: hel,o,part,temp_part
  real(kind=8),dimension(3) :: matrix2
  real(kind=8) :: amp2,amp_col
  complex(kind=8) :: amp2_c,amp_col_c
  real(kind=8),dimension(:),allocatable :: mass
  logical :: done
  
  call get_run_arguments()

  write(*,*) 'itteittenn'
  call cpu_time(tTot_B)

  if (.not.allocated(o)) allocate(o(next))
  allocate(hel(next))
  allocate(mass(next))
  allocate(p(0:3,next))

  mass(1:next)=0d0
  write(*,*) 'fali'
  call create_run_tag_and_open_files()

  write(*,*) 'hahhahoo'
  call cpu_time(tBefore)

  if (.not.allocated(part)) allocate(part(1:next))
  if (.not.allocated(temp_part)) allocate(temp_part(1:next))
  call read_event(11,done)
  rewind(11)

  allocate(amps((next-2)*(next-2))) 
  temp_part=part
  call amps(1)%init(2,next,temp_part,o)
  col_acc=20
  call amps(1)%init_col2(next,o,col_acc)
  if (color_flow) then
        do i=2,2 ! TV: only single external U(1) for now!
         skip = 1
         do k=1,next-2 ! loop through alil gluons
          temp_part = part
          do j=1,next
           if (temp_part(j).eq.21) then
               if (j.le.skip) then
                   cycle
               endif
               temp_part(j) = 22
               skip = j
               exit
           endif           
           enddo
          call amps(i+k-1)%init(2,next,temp_part,o)
          col_acc=20
          call amps(i+k-1)%init_col2(next,o,col_acc)
         enddo
        enddo
  endif


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


     call amps(1)%evaluate(next,p,ihel)
     if (color_flow) then
       do i=2,2
        do k=1,next-2
         call amps(i+k-1)%evaluate(next,p,ihel)
        enddo
       enddo
     endif

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore

     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amps(1)%n_qqbar.eq.0) then
           do irow=1,amp_QCD%nColOrd
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
           ri_end=0
           if (color_flow) ri_end= 1 !next-2 for FC
           do ri=0,ri_end ! loop over no U(1) and one U(1) in the rows
              iacc_in=iacc
           do k=1,ri*(next-2)+1
              proc_num = ri+1+k-1
              if (proc_num.gt.1) iacc_in=min(3,iacc+1)
           do irow=1,amps(proc_num)%nColOrd
              amp_col_c=(0d0,0d0)
              irow_mat = irow
              do i=1,amps(proc_num)%n_col_vals(iacc)
                 amp2_c=(0d0,0d0)
                 do ic=amps(proc_num)%row_index(irow_mat-1,i,iacc)+1,amps(proc_num)%row_index(irow_mat,i,iacc)
                    icol=amps(proc_num)%col_index(ic,i,iacc)
                    icol_mat = icol
                    if (proc_num.eq.1) then
                       amp2_c=amp2_c+amps(proc_num)%amps(icol_mat)
                    else
                       if (icol_mat.eq.irow) then
                          amp2_c=amp2_c+amps(proc_num)%amps(icol_mat)
                       endif
                    endif
                 enddo
                 if (proc_num.gt.1) then
                   amp_col_c=amp_col_c+amp2_c*amps(proc_num)%diff_col_vals(i,iacc)*(-1d0/3d0)
                 else
                   amp_col_c=amp_col_c+amp2_c*amps(proc_num)%diff_col_vals(i,iacc)
                 endif
              enddo
              matrix2(iacc_in)=matrix2(iacc_in)+dble(amp_col_c*conjg(amps(proc_num)%amps(irow)))
           enddo
           enddo
           enddo
        endif

        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo

     !write(*,*) 'matrix LC',matrix2(1)
     !write(*,*) 'matrix NLC',matrix2(2)

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

    if (argc.le.8) then
       write(*,*) 'Inconsistent arguments:'
       write(*,*) '--------- Should be: --------'
       write(*,*) 'next, *process*, *order*'
       stop 2
    else
       do i = 1, argc
          CALL GET_COMMAND_ARGUMENT(i, argv)
          if (i.eq.1) then
             read(argv,*) next
             allocate(part(1:next))
             allocate(o(1:next))
          endif
          do k=0,next
          if (i.eq.2+k) then
             read(argv,*) part(k+1)
          endif
          enddo
          do k=0,next
          if (i.eq.2+next+k) then
               read(argv,*) o(k+1)
          endif
          enddo
       enddo
    endif
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
    endif

    write(*,*) next
    write(*,*) part
    write(*,*) o
  end subroutine get_run_arguments

  subroutine create_run_tag_and_open_files()
    use arguments
    implicit none
    character(len=1) :: s1
    character(len=2) :: s2
    character(len=string_len) :: tag,tag_read
    tag=''       ! tag of current run
    tag_read=''  ! same as 'tag', but with previous imode (i.e., defines the file to read the integration grids from)
    call add_to_string(tag,next,.true.)
    call add_to_string(tag_read,next,.true.)
    call add_to_string(tag,2,.true.)
    call add_to_string(tag_read,2,.true.)
    do i=1,next
       call add_to_string(tag,part(i),.true.)
       call add_to_string(tag_read,part(i),.true.)
    enddo
    do i=1,next
       call add_to_string(tag,o(i),.true.)
       call add_to_string(tag_read,o(i),.true.)
    enddo
    call fill_string(tag,len(trim(tag)))
    call fill_string(tag_read,len(trim(tag)))
    open(unit=11,file='Outputs/events'//trim(tag)//'.lhe',status='old')
    open(unit=12,file='Outputs/events'//trim(tag)//'.lhe.rwgt',status='unknown')
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

  subroutine fill_string(string,size)
    ! Fills the string 'string' with leading underscores until the string has
    ! size 'size'. The declaration of the string must be at least size 'size'.
    implicit none
    character(len=string_len) :: string
    integer :: size,n_to_add
    if (size.gt.len(string)) then
       write (*,*) 'Size greater than string',size,string
       stop 1
    endif
    n_to_add=len(trim(string))+2-len(trim(string))
    do i=1,n_to_add
       string='_'//trim(adjustl(string))
    enddo
  end subroutine fill_string


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
    rwgt_full=matrix2(2)/matrix2(1)
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
