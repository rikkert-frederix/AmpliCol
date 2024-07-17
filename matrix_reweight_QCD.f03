! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module common
  use amplitude_QCD_mod
  implicit none
  integer :: next
  type(amplitude_QCD) :: amp_QCD
  type(amplitude_QCD),dimension(:),allocatable :: amps
  real(kind=8),dimension(:,:),allocatable :: p
  integer,parameter :: string_len=150
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
  integer :: i,j,col_acc,icol,irow,ic,iacc
  integer :: icol_mat,irow_mat,ri,ri_end,m,proc_num,iacc_in,k,skip 
  integer,dimension(:),allocatable :: hel,o,part,part_sf,orig_part,temp_part
  integer,dimension(:,:),allocatable :: spin
  real(kind=8),dimension(3) :: matrix2
  real(kind=8) :: amp2,amp_col
  complex(kind=8) :: amp2_c,amp_col_c
  real(kind=8),dimension(:),allocatable :: mass,width
  logical :: done,first
  integer :: swap_q,swap_aq
  integer :: it,gi,gi_iperm,gi_prev ! type for 2qq process
  integer,dimension(:),allocatable :: iper_test
  integer :: ic_low,ic_upp
  integer,dimension(8) :: flav
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
  it = 1
  call amps(1)%init(2,next,part,spin,mass,width,o,it)
  col_acc=20
  call amps(1)%init_col2(next,orig_part,o,it,col_acc)


  if (amps(1)%n_qqbar.eq.2) then
     amps(2)%n_qqbar=amps(1)%n_qqbar
     amps(2)%same_flav=amps(1)%same_flav
     ! the other colour order with the two anti-quarks interchanged.
     it = 2
     call amps(2)%init(2,next,part,spin,mass,width,o,it)
     col_acc=20
     call amps(2)%init_col2(next,orig_part,o,it,col_acc)
  endif

  if (amps(1)%n_qqbar.eq.2.and.amps(1)%same_flav) then
     amps(3)%n_qqbar=amps(1)%n_qqbar
     amps(3)%same_flav=amps(1)%same_flav
     part_sf(:)=orig_part(:)
     call define_symm_2qq(next,part_sf,2)
     it = 1
     call amps(3)%init(2,next,part_sf,spin,mass,width,o,it)
     col_acc=20
     call amps(3)%init_col2(next,orig_part,o,it,col_acc)

     amps(4)%n_qqbar=amps(1)%n_qqbar
     amps(4)%same_flav=amps(1)%same_flav
     it = 2
     call amps(4)%init(2,next,part_sf,spin,mass,width,o,it)
     col_acc=20
     call amps(4)%init_col2(next,orig_part,o,it,col_acc)
  endif

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
           it = 0! dummy
           call amps(i+k-1)%init(2,next,temp_part,spin,mass,width,o,it)
           col_acc=20
           call amps(i+k-1)%init_col2(next,part,o,it,col_acc)
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

     call amps(1)%evaluate(next,p,mass,width,hel,part)
     if (amps(1)%n_qqbar.eq.2) then
        call amps(2)%evaluate(next,p,mass,width,hel,part)
        if (amps(1)%same_flav) then
           call amps(3)%evaluate(next,p,mass,width,hel,part)
           call amps(4)%evaluate(next,p,mass,width,hel,part)
           amps(1)%amps(:)=amps(1)%amps(:)+(1d0/3d0)*amps(3)%amps(:)
           amps(2)%amps(:)=(1d0/3d0)*amps(2)%amps(:)+amps(4)%amps(:)
        endif
     endif


     if (color_flow) then
        do i=2,2
           do k=1,next-2
              call amps(i+k-1)%evaluate(next,p,mass,width,hel,part)
           enddo
        enddo
     endif

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore


!!! ###########################################     
     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        if (amps(1)%n_qqbar.eq.0) then
           do irow=1,amps(1)%nColOrd
              if (use_real_gluons) then
                 amp_col=0d0
              else
                 amp_col_c=(0d0,0d0)
              endif
              do i=1,amps(1)%n_col_vals(iacc,1)
                 if (use_real_gluons) then
                    amp2=0d0
                 else
                    amp2_c=(0d0,0d0)
                 endif
                 do ic=amps(1)%row_index(irow-1,i,iacc,1)+1,amps(1)%row_index(irow,i,iacc,1)
                    icol=amps(1)%col_index(amps(1)%i_col_i(i,iacc,1)+ic,1)
                    if (use_real_gluons) then
                       amp2=amp2+amps(1)%amps_r(icol)
                    else
                       amp2_c=amp2_c+amps(1)%amps(icol)
                    endif
                 enddo
                 if (use_real_gluons) then
                    amp_col=amp_col+amp2*amps(1)%diff_col_vals(i,iacc,1)
                 else
                    amp_col_c=amp_col_c+amp2_c*amps(1)%diff_col_vals(i,iacc,1)
                 endif
              enddo
              if (use_real_gluons) then
                 matrix2(iacc)=matrix2(iacc)+amp_col*amps(1)%amps_r(irow)
              else
                 matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps(1)%amps(irow)))
              endif
           enddo


        elseif (amps(1)%n_qqbar.eq.1) then
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
                    do i=1,amps(proc_num)%n_col_vals(iacc,1)
                       amp2_c=(0d0,0d0)
                       do ic=amps(proc_num)%row_index(irow_mat-1,i,iacc,1)+1,amps(proc_num)%row_index(irow_mat,i,iacc,1)
                          icol=amps(1)%col_index(amps(1)%i_col_i(i,iacc,1)+ic,1)
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
                          amp_col_c=amp_col_c+amp2_c*amps(proc_num)%diff_col_vals(i,iacc,1)*(-1d0/3d0)
                       else
                          amp_col_c=amp_col_c+amp2_c*amps(proc_num)%diff_col_vals(i,iacc,1)
                       endif
                    enddo
                    matrix2(iacc_in)=matrix2(iacc_in)+dble(amp_col_c*conjg(amps(proc_num)%amps(irow)))
                 enddo
              enddo
           enddo


        elseif (amps(1)%n_qqbar.eq.2) then
           do it=1,2
              do irow=1,amps(it)%nColOrd
                 if (next.ge.5) then ! there is at least one gluon
                    iper_test(1:next)=amps(it)%perm(1:next,irow) !
                    do i=2,next-2
                       if ((abs(part(iper_test(i))).ge.1.and.abs(part(iper_test(i))).le.6)) then
                          gi = i - 2 ! number of gluons in the colour order between the first quark and anti-quark
                          exit
                       endif
                    enddo
                    gi_iperm = gi + 1
                 else
                    gi_iperm=1
                 endif
                 if (next.ge.5) then
                    if (irow.ge.2) then
                       iper_test(1:next)=amps(it)%perm(1:next,irow-1) !
                       do i=2,next-2
                          if ((abs(part(iper_test(i))).ge.1.and.abs(part(iper_test(i))).le.6)) then
                             gi = i - 2 ! number of gluons in the colour order between the first quark and anti-quark
                             exit
                          endif
                       enddo
                       gi_prev = gi + 1
                    else
                       gi_prev = 1
                    endif
                 else
                    gi_prev=1
                 endif

                 amp_col_c=(0d0,0d0)
                 do i=1,amps(it)%n_col_vals(iacc,gi_iperm) 
                    amp2_c=(0d0,0d0)
                    ic_low = amps(it)%row_index(irow-1,i,iacc,gi_prev)+1
                    if(gi_prev.ne.gi_iperm) then
                       ic_low = 1
                    endif
                    ic_upp = amps(it)%row_index(irow,i,iacc,gi_iperm)
                    do ic = ic_low,ic_upp
                       icol=amps(it)%col_index(amps(it)%i_col_i(i,iacc,gi_iperm)+ic,gi_iperm)
                       if (icol.le.amps(it)%nColOrd) then
                          amp2_c=amp2_c+amps(1)%amps(icol)
                       else
                          amp2_c=amp2_c+amps(2)%amps(icol-amps(it)%nColOrd)
                       endif
                 enddo
                    amp_col_c=amp_col_c+amp2_c*amps(it)%diff_col_vals(i,iacc,gi_iperm)
                 enddo
                 matrix2(iacc)=matrix2(iacc)+dble(amp_col_c*conjg(amps(it)%amps(irow)))
              enddo
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
    integer :: argc,nquarks
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
    character(len=1) :: s1
    character(len=2) :: s2
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
