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
  implicit none
  integer :: next
  type(amplitude_QCD),dimension(:),allocatable :: amps
  real(kind=8),dimension(:,:),allocatable :: p
  integer,parameter :: string_len=150
  integer :: i,j,col_acc,icol,irow,ic,iacc
  integer,dimension(:),allocatable :: hel,o,part,part_sf,orig_part,temp_part
  integer,dimension(:,:),allocatable :: spin
  real(kind=8),dimension(3) :: matrix2
  real(kind=8) :: amp2,amp_col
  complex(kind=8) :: amp2_c,amp_col_c
  real(kind=8),dimension(:),allocatable :: mass,width
  logical :: done
  integer :: gi,gi_iperm,gi_prev ! type for 2qq process
  integer,dimension(:),allocatable :: iper_test
  integer :: ic_low,ic_upp
  character(len=string_len) :: tag,tag_read,add_arg=''
  
  call get_run_arguments()

  call cpu_time(tTot_B)

  if (.not.allocated(o)) allocate(o(next))
  allocate(hel(next))
  allocate(mass(next))
  allocate(width(next))
  allocate(p(0:3,next))
  allocate(iper_test(1:next)) ! needed for 2qq

  mass(1:next)=0d0
  width(1:next)=0d0

!  mass(1:2)= 0d0
!  mass(3:4)= 173.d0
!  mass(5)  = 0d0
!  width(1:2)= 0d0
!  width(3:4)= 1.4915d0
!  width(5)  = 0d0


  call create_run_tag_and_open_files()

  call cpu_time(tBefore)

  if (.not.allocated(part)) allocate(part(1:next))
  if (.not.allocated(part_sf)) allocate(part_sf(1:next))
  if (.not.allocated(orig_part)) allocate(orig_part(1:next))
  if (.not.allocated(temp_part)) allocate(temp_part(1:next))
  call read_event(11,done)
  rewind(11)

  allocate(amps((next-2)*(next-2))) 
  orig_part(:)=part(:)

  ! counting of quark flavours in process
  call fill_quark_info()

  call define_symm_2qq(next,part,1)
  call amps(1)%init(2,next,part,spin,mass,width,o)
  col_acc=20
!!$  call amps(1)%init_col(next,orig_part,col_acc)
  call amps(1)%init_col(next,part,col_acc)

  if (amps(1)%n_qqbar.eq.2.and.amps(1)%same_flav) then
     amps(3)%n_qqbar=amps(1)%n_qqbar
     amps(3)%same_flav=amps(1)%same_flav
     part_sf(:)=orig_part(:)
     call define_symm_2qq(next,part_sf,2)
     call amps(3)%init(2,next,part_sf,spin,mass,width,o,amps(1))
  endif


  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  do
     call read_event(11,done)
     if (done) exit
     matrix2(1:3)=0d0

     call cpu_time(tBefore)

     call amps(1)%evaluate(next,p,mass,width,hel,part)
     if (amps(1)%n_qqbar.eq.2 .and. amps(1)%same_flav) then
        call amps(3)%evaluate(next,p,mass,width,hel,part)
        amps(1)%amps(:)=amps(1)%amps(:)+amps(3)%amps(:)
     endif



     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore


!!! ###########################################     
     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amps(1)%n_qqbar.eq.0) then
           do irow=1,amps(1)%n_amps
              if (use_real_gluons) then
                 amp_col=0d0
              else
                 amp_col_c=(0d0,0d0)
              endif
              do i=1,amps(1)%n_col_vals(iacc)
                 if (use_real_gluons) then
                    amp2=0d0
                 else
                    amp2_c=(0d0,0d0)
                 endif
                 do ic=amps(1)%row_index(irow-1,i,iacc)+1,amps(1)%row_index(irow,i,iacc)
                    icol=amps(1)%col_index(amps(1)%i_col_i(i,iacc)+ic)
                    if (use_real_gluons) then
                       amp2=amp2+amps(1)%amps_r(icol)
                    else
                       amp2_c=amp2_c+amps(1)%amps(icol)
                    endif
                 enddo
                 if (use_real_gluons) then
                    amp_col=amp_col+amp2*amps(1)%diff_col_vals(i,iacc)
                 else
                    amp_col_c=amp_col_c+amp2_c*amps(1)%diff_col_vals(i,iacc)
                 endif
              enddo
              if (use_real_gluons) then
                 matrix2(iacc)=matrix2(iacc)+amp_col*amps(1)%amps_r(irow)
              else
                 matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps(1)%amps(irow)))
              endif
           enddo


        elseif (amps(1)%n_qqbar.eq.1 .or. amps(1)%n_qqbar.eq.2) then
           do irow=1,amps(1)%n_amps
              amp_col_c=(0d0,0d0)
              do i=1,amps(1)%n_col_vals(iacc)
                 amp2_c=(0d0,0d0)
                 do ic=amps(1)%row_index(irow-1,i,iacc)+1,amps(1)%row_index(irow,i,iacc)
                    icol=amps(1)%col_index(amps(1)%i_col_i(i,iacc)+ic)
                    amp2_c=amp2_c+amps(1)%amps(icol)
                 enddo
                 amp_col_c=amp_col_c+amp2_c*amps(1)%diff_col_vals(i,iacc)
              enddo
              matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps(1)%amps(irow)))
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

  subroutine fill_quark_info()
    implicit none
    integer,dimension(8) :: flav
    integer :: k

    flav = 0
    k = 1
    amps%n_qqbar= 0
    amps(1)%same_flav=.true.
    do i=1,next
       if (i.le.2) then
          if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
             flav(k) = abs(orig_part(i))
             k= k+1
             if (orig_part(i).lt.0) amps(1)%n_qqbar=amps(1)%n_qqbar+1
          endif
       else
          if (orig_part(i).ne.21 .and. orig_part(i).ne.22) then
             flav(k) = abs(orig_part(i))
             k= k+1
             if (orig_part(i).gt.0) amps(1)%n_qqbar=amps(1)%n_qqbar+1
          endif
       endif
    enddo

    if (any(flav(1:2*amps(1)%n_qqbar).ne.flav(1))) amps(1)%same_flav = .false.
  end subroutine fill_quark_info


  subroutine get_run_arguments()
    use arguments
    implicit none
    integer :: argc,i,k
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
             if (next.le.3) then
                write (*,*) 'Need at least 4 particles (2->2 scattering)',next
                stop 1
             endif
             allocate(part(1:next))
             allocate(o(1:next))
          endif
          do k=0,next-1
             if (i.eq.2+k) then
                read(argv,*) part(k+1)
             endif
          enddo
          do k=0,next-1
             if (i.eq.2+next+k) then
                read(argv,*) o(k+1)
             endif
          enddo
          if (argc.eq.1+2*next +1 .and. i.eq.argc) then
             ! Special case: we have an additional argument. Use it as a special tag
             read(argv,*) add_arg
          endif
       enddo
    endif

    call setup_spin()
    
    if (next.lt.4) then
       write (*,*) 'Not enough external particles',next
       stop 1
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
       call add_to_string(tag,part(i),.true.)
       call add_to_string(tag_read,part(i),.true.)
    enddo
    do i=1,next-1
       call add_to_string(tag,o(i),.true.)
       call add_to_string(tag_read,o(i),.true.)
    enddo
    call add_to_string(tag,o(next),.false.)
    call add_to_string(tag_read,o(next),.false.)
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
    read (iunit,*,err=99,end=99) dum,evt_wgt,wgt,amp2,weight
    read (iunit,*,err=99,end=99) hel(1:next)
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

  subroutine define_symm_2qq(next,part,chan)
    implicit none
    integer :: next,chan
    integer, dimension(next) :: part
    integer :: i,j,sgn
    logical :: first
    if (amps(1)%same_flav) then
       if (chan.eq.2) then
          do i=1,next
             if (abs(part(i)).gt.0.and.abs(part(i)).lt.6) then
                first=.true.
                do j=i+1,next
                   if (i.le.2.and.j.le.2) sgn=-1
                   if (i.le.2.and.j.gt.2) sgn=+1
                   if (i.gt.2.and.j.gt.2) sgn=-1
                   if (part(j).eq.sgn*part(i).and..not.first) then
                      part(i) = sign(abs(part(i))+1,part(i))
                      part(j) = sgn*(part(i))
                      exit
                   endif
                   if (part(j).eq.sgn*part(i).and.first) then
                      first = .false.
                   endif
                enddo
             endif
          enddo
       elseif (chan.eq.1) then
          do i=1,next
             if (abs(part(i)).gt.0.and.abs(part(i)).lt.6) then
                do j=i+1,next
                   if (i.le.2.and.j.le.2) sgn=-1
                   if (i.le.2.and.j.gt.2) sgn=+1
                   if (i.gt.2.and.j.gt.2) sgn=-1
                   if (part(j).eq.sgn*part(i)) then
                      part(i) = sign(abs(part(i))+1,part(i))
                      part(j) = sgn*(part(i))
                      exit
                   endif
                enddo
                exit
             endif
          enddo
       endif
    endif
  end subroutine define_symm_2qq


end program matrix_reweight
