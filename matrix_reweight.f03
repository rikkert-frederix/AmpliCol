
! gfortran -mcmodel=large -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 matrix_reweight.f03

module common
  use amplitude_mod
  implicit none
  integer :: next
  type(amplitude) :: amplitudes_LC,amplitudes_NLC,amplitudes_full
  type(amplitude_cache) :: amplitudes_cache
  real(kind=8),dimension(:,:),allocatable :: p
  type(col_amp) :: color_amp
!!$  type(col_amp),dimension(:),allocatable :: col_amp_list
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
module random_colours
  implicit none
  integer :: n_col
  integer,dimension(:,:),allocatable :: n_cols,col_labels
  integer,dimension(:,:,:),allocatable :: color_labels
  real(kind=8) :: factor_wgt
  real(kind=8),dimension(:),allocatable :: accum_color_probs
end module random_colours

program matrix_reweight
  use common
  use timings
  implicit none
  ! allowed reweight modes
  ! 0 : reweight with matrix elements summed over colors
  ! 1 : reweight with matrix elements limited to the color row used to generate the LC event
  ! 2 : reweight with random color assignment (3**(2*n) random assignments)
  ! 3 : same as reweight_mode=2, but with zero's cycled over
  ! 4 : same as reweight_mode=2, but improved color assigments (e.g., making sure that the number of colors is equal to anti-colors)
  ! 5 : same as reweight_mode=4, but with zero's cycled over)  -- this gives equivalent results to reweight_mode=3
  ! 6 : same as reweight_mode=4, but making unlikely assigments (e.g., all the same color) more likely but with smaller weight
  ! 7 : same as reweight_mode=6, but with zero's cycled over
  ! 8 : pick color assignment compatible with color order used to generate the LC event, and reweight that
  ! 9 : same as reweight_mode=8, but checking compatible colours event-by-event (requires less memory and is faster if n is large).
  !10 : same as reweight_mode=9, but using new code to evaluate amplitudes
  !11 : same as reweight_mode=10, but with new way of assigning random colours
  !12 : same as reweight_mode=11, but with re-using wavefunctions
  !13 : same as reweight_mode=11, but allows for setting number of colours
  !14 : same as reweight_mode=11, but looping over colours
  !15 : same as 0, but smarter summing over the colour matrix (using CSR format)
  !16 : same as 15, but with better caching of interactions
  !17 : same as 16, but with 4-gloun vertex replaced by two 3-vertices
  ! ...
  ! All reweight_modes<=9 can be run for random helicity assignments, or summed over helicities
  integer,parameter :: reweight_mode=17
  logical,parameter :: sum_hel=.false.
  integer,parameter :: n_rwgt=1   ! number of colour configurations to use average over per event (use n_rwgt=1 for reweight_mode<=1 or reweight_mode >=15)
  integer :: i,j,k,col_acc,icol,ih,iperm,jperm,ihel,iperm_ev,hel_picked,n_valid,irow,ic
  integer,dimension(:),allocatable :: hel,list,list_valid_iperm,o
  integer,dimension(:,:),allocatable :: list_orders
  real(kind=8) :: amp2_LC(n_rwgt),amp2_NLC(n_rwgt),amp2_full(n_rwgt),amp2_hel,color_wgt,amp,amp2,amp_col
  real(kind=8),dimension(2) :: ievt_count
  real(kind=8),dimension(:),allocatable :: mass
  real(kind=8),dimension(:,:),allocatable :: icount
  real(kind=8),external :: ran2
  logical :: done

  call get_run_arguments()

  call cpu_time(tTot_B)

  allocate(o(next))
  allocate(hel(next))
  allocate(mass(next))
  allocate(p(0:3,next))

  if (reweight_mode.ge.2 .and. reweight_mode.le.9) allocate(list(0:factorial(next-1)))
  if (reweight_mode.ge.10 .and. reweight_mode.le.14) allocate(list_orders(next,factorial(next-1)))
  if (reweight_mode.ge.2 .and. reweight_mode.le.14) allocate (list_valid_iperm(1:factorial(next-1)))

!!$  allocate(icount(factorial(next-1),factorial(next-1)))
!!$  icount=0d0

  mass(1:next)=0d0
  call create_run_tag_and_open_files()

  if (reweight_mode.eq.8 .or. reweight_mode.eq.14) then
     call read_event(11,done)
     call iperm_encode(next-1,o,iperm_ev)
     rewind(11)
  endif
  if ((reweight_mode.ge.4 .and. reweight_mode.le.8) .or. reweight_mode.eq.14) call init_color_probs(next)

  ievt_count(1:2)=0d0

  call cpu_time(tBefore)
  if (reweight_mode.eq.8 .or. reweight_mode.eq.9) then
     col_acc=0
     call amplitudes_NLC%init(next,col_acc,sum_hel)
  elseif (reweight_mode.ge.0 .and. reweight_mode.le.7) then
     col_acc=0
     call amplitudes_LC%init_onlycol(next,col_acc,sum_hel)
     col_acc=1
     call amplitudes_NLC%init(next,col_acc,sum_hel)
     col_acc=20
     call amplitudes_full%init_onlycol(next,col_acc,sum_hel)
  elseif (reweight_mode.eq.15) then
     col_acc=0
     call amplitudes_LC%init_onlycol_CSR(next,col_acc,sum_hel)
     col_acc=1
     call amplitudes_NLC%init_CSR(next,col_acc,sum_hel)
     col_acc=20
     call amplitudes_full%init_onlycol_CSR(next,col_acc,sum_hel)
  elseif(reweight_mode.eq.16) then
     call amplitudes_cache%setup_imap_cache(next)
     
     col_acc=1
     if (col_acc.le.1) then
        call amplitudes_cache%setup_colmap_cache_NLC(col_acc)
     else
        call amplitudes_cache%setup_colmap_cache(col_acc)
     endif
  elseif(reweight_mode.eq.17) then

     call amplitudes_cache%setup_imap_3vert(next)
     
     col_acc=1
     if (col_acc.le.1) then
        call amplitudes_cache%setup_colmap_cache_NLC(col_acc)
     else
        call amplitudes_cache%setup_colmap_cache(col_acc)
     endif
  endif
  call cpu_time(tAfter)
  t_amp_init=t_amp_init+tAfter-tBefore

  do 
     call read_event(11,done)
     if (done) exit
     amp2_LC(1:n_rwgt)=0d0
     amp2_NLC(1:n_rwgt)=0d0
     amp2_full(1:n_rwgt)=0d0

     if (reweight_mode.eq.1 .or. (reweight_mode.ge.9 .and. reweight_mode.le.13)) call iperm_encode(next-1,o,iperm_ev)

     do k=1,n_rwgt

        if (reweight_mode.ge.2 .and. reweight_mode.le.14) then
           call cpu_time(tBefore)
           call pick_random_color(next)
           call find_permutations(next-1)
           if (reweight_mode.ge.2 .and. reweight_mode.le.9) then
              list(0)=n_valid
              list(1:n_valid)=list_valid_iperm(1:n_valid)
           endif
           call cpu_time(tAfter)
           t_ran=t_ran+tAfter-tBefore
        endif

        call cpu_time(tBefore)
        if (sum_hel) then
           if (reweight_mode.ge.10 .and. reweight_mode.le.14) then
              write (*,*) 'sum_hel=.true. is not a viable'&
                   //' option with reweight_mode == 10, 11, 12, 13'
              stop
           endif
           call amplitudes_NLC%evaluate(p)
        else
           ! read helicity from event file
           ihel=hel_picked
           do i=1,next
              if (btest(ihel,i-1)) then
                 hel(i)=1
              else
                 hel(i)=0
              endif
           enddo
           if ((reweight_mode.ge.0 .and. reweight_mode.le.9) .or. reweight_mode.eq.15) then
              call amplitudes_NLC%evaluate(p,hel)
           elseif (reweight_mode.eq.16) then
              call amplitudes_cache%evaluate_cache(p,hel)
           elseif (reweight_mode.eq.17) then
              call amplitudes_cache%evaluate_3vert(p,hel)
           endif
        endif
        call cpu_time(tAfter)
        t_amp=t_amp+tAfter-tBefore

        call cpu_time(tBefore)
        ! Leading color matrix elements
        if (reweight_mode.le.7) then
           do icol=1,amplitudes_LC%colmap(0,0)
              iperm=amplitudes_LC%colmap(1,icol)
              jperm=amplitudes_LC%colmap(2,icol)
              if (reweight_mode.eq.1) then
                 if (iperm.ne.iperm_ev) cycle
              endif
              do ih=0,amplitudes_NLC%nhel(amplitudes_NLC%isize+1)-1
                 amp2_LC(k)=amp2_LC(k)+amplitudes_NLC%amps(amplitudes_NLC%helmap(iperm,ih),iperm)* &
                      amplitudes_LC%colmap(0,icol)* &
                      amplitudes_NLC%amps(amplitudes_NLC%helmap(jperm,ih),jperm)
              enddo
           enddo
        elseif (reweight_mode.eq.15) then
           ih=0 ! no summing over helicities
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,(next+1)/2
                 amp2=0d0
                 do ic=amplitudes_LC%row_index(irow-1,i)+1,amplitudes_LC%row_index(irow,i)
                    icol=amplitudes_LC%col_index(ic,i)
                    amp2=amp2+amplitudes_NLC%amps(amplitudes_NLC%helmap(icol,ih),icol)
                 enddo
                 if (i.eq.1) then
                    amp_col=amp_col+amp2*amplitudes_LC%col_value(i)
                 else
                    amp_col=amp_col+amp2*amplitudes_LC%col_value(i)*2d0
                 endif
              enddo
              amp2_LC(k)=amp2_LC(k)+amp_col*amplitudes_NLC%amps(amplitudes_NLC%helmap(irow,ih),irow)
           enddo
        elseif (reweight_mode.eq.16 .or. reweight_mode.eq.17) then
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,1
                 amp2=0d0
                 do ic=amplitudes_cache%row_index_LC(irow-1,i)+1,amplitudes_cache%row_index_LC(irow,i)
                    icol=amplitudes_cache%col_index_LC(ic,i)
                    amp2=amp2+amplitudes_cache%amps(icol)
                 enddo
                 amp_col=amp_col+amp2*amplitudes_cache%col_value_LC(i)
              enddo
              amp2_LC(k)=amp2_LC(k)+amp_col*amplitudes_cache%amps(irow)
           enddo
        elseif (reweight_mode.eq.8 .or. reweight_mode.eq.9) then
           do ih=0,amplitudes_NLC%nhel(amplitudes_NLC%isize+1)-1
              amp2_hel=0d0
              do icol=1,list(0)
                 amp2_hel=amp2_hel+amplitudes_NLC%amps(amplitudes_NLC%helmap(list(icol),ih),list(icol))**2
              enddo
              amp2_LC(k)=amp2_LC(k)+amp2_hel
           enddo
        elseif (reweight_mode.ge.10 .and. reweight_mode.le.13 .and. reweight_mode.ne.12) then
           do icol=1,n_valid
              call evaluate_order_v2(next,p,list_orders(1,icol),hel,amp)
              amp2_LC(k)=amp2_LC(k)+amp**2
              amp2_full(k)=amp2_full(k)+amp
           enddo
           amp2_full(k)=amp2_full(k)**2
        elseif (reweight_mode.eq.14) then
!!$           if (k.eq.1) then
!!$              call amplitudes_LC%evaluate(p,hel)
!!$              do icol=1,amplitudes_LC%colmap(0,0)
!!$                 iperm=amplitudes_LC%colmap(1,icol)
!!$                 jperm=amplitudes_LC%colmap(2,icol)
!!$                 do ih=0,amplitudes_LC%nhel(amplitudes_LC%isize+1)-1
!!$                    amp2_LC(k)=amp2_LC(k)+amplitudes_LC%amps(amplitudes_LC%helmap(iperm,ih),iperm)* &
!!$                         amplitudes_LC%colmap(0,icol)* &
!!$                         amplitudes_LC%amps(amplitudes_LC%helmap(jperm,ih),jperm)
!!$                 enddo
!!$              enddo
!!$           endif
           do icol=1,n_valid
              call evaluate_order_v2(next,p,list_orders(1,icol),hel,amp)
              amp2_LC(k)=amp2_LC(k)+amp**2
              amp2_full(k)=amp2_full(k)+amp
           enddo
           amp2_full(k)=amp2_full(k)**2*color_wgt
        elseif (reweight_mode.eq.12) then
           allocate(col_amp_list(n_valid))
           do icol=1,n_valid
              call col_amp_list(icol)%evaluate_order_v3(icol,next,p,list_orders(1,icol),hel,amp)
              amp2_LC(k)=amp2_LC(k)+amp**2
              amp2_full(k)=amp2_full(k)+amp
           enddo
           deallocate(col_amp_list)
           amp2_full(k)=amp2_full(k)**2
        endif
        call cpu_time(tAfter)
        t_mat_LC=t_mat_LC+tAfter-tBefore

        ! Next-to-leading color matrix elements
        if (reweight_mode.eq.0 .or. reweight_mode.eq.1) then
           call cpu_time(tBefore)
           do icol=1,amplitudes_NLC%colmap(0,0)
              iperm=amplitudes_NLC%colmap(1,icol)
              jperm=amplitudes_NLC%colmap(2,icol)
              if (reweight_mode.eq.1) then
                 if (iperm.ne.iperm_ev) cycle
              endif
              do ih=0,amplitudes_NLC%nhel(amplitudes_NLC%isize+1)-1
                 amp2_NLC(k)=amp2_NLC(k)+amplitudes_NLC%amps(amplitudes_NLC%helmap(iperm,ih),iperm)* &
                      amplitudes_NLC%colmap(0,icol)* &
                      amplitudes_NLC%amps(amplitudes_NLC%helmap(jperm,ih),jperm)
              enddo
           enddo
           call cpu_time(tAfter)
           t_mat_NLC=t_mat_NLC+tAfter-tBefore
        elseif (reweight_mode.eq.15) then
           call cpu_time(tBefore)
           ih=0 ! no summing over helicities
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,(next+1)/2
                 amp2=0d0
                 do ic=amplitudes_NLC%row_index(irow-1,i)+1,amplitudes_NLC%row_index(irow,i)
                    icol=amplitudes_NLC%col_index(ic,i)
                    amp2=amp2+amplitudes_NLC%amps(amplitudes_NLC%helmap(icol,ih),icol)
                 enddo
                 if (i.eq.1) then
                    amp_col=amp_col+amp2*amplitudes_NLC%col_value(i)
                 else
                    amp_col=amp_col+amp2*amplitudes_NLC%col_value(i)*2d0
                 endif
              enddo
              amp2_NLC(k)=amp2_NLC(k)+amp_col*amplitudes_NLC%amps(amplitudes_NLC%helmap(irow,ih),irow)
           enddo
           call cpu_time(tAfter)
           t_mat_NLC=t_mat_NLC+tAfter-tBefore
        elseif (reweight_mode.eq.16 .or. reweight_mode.eq.17) then
           call cpu_time(tBefore)
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,2
                 amp2=0d0
                 do ic=amplitudes_cache%row_index_NLC(irow-1,i)+1,amplitudes_cache%row_index_NLC(irow,i)
                    icol=amplitudes_cache%col_index_NLC(ic,i)
                    amp2=amp2+amplitudes_cache%amps(icol)
                 enddo
                 if (i.eq.1) then
                    amp_col=amp_col+amp2*amplitudes_cache%col_value_NLC(i)
                 else
                    amp_col=amp_col+amp2*amplitudes_cache%col_value_NLC(i)*2d0
                 endif
              enddo
              amp2_NLC(k)=amp2_NLC(k)+amp_col*amplitudes_cache%amps(irow)
           enddo
           call cpu_time(tAfter)
           t_mat_NLC=t_mat_NLC+tAfter-tBefore
       endif


        ! Full color matrix elements
        if (reweight_mode.eq.0 .or. reweight_mode.eq.1) then
           call cpu_time(tBefore)
           do icol=1,amplitudes_full%colmap(0,0)
              iperm=amplitudes_full%colmap(1,icol)
              jperm=amplitudes_full%colmap(2,icol)
              if (reweight_mode.eq.1) then
                 if (iperm.ne.iperm_ev) cycle
              endif
              do ih=0,amplitudes_NLC%nhel(amplitudes_NLC%isize+1)-1
                 amp2_full(k)=amp2_full(k)+amplitudes_NLC%amps(amplitudes_NLC%helmap(iperm,ih),iperm)* &
                      amplitudes_full%colmap(0,icol)* &
                      amplitudes_NLC%amps(amplitudes_NLC%helmap(jperm,ih),jperm)
              enddo
           enddo
           call cpu_time(tAfter)
           t_mat_full=t_mat_full+tAfter-tBefore
        elseif (reweight_mode.eq.15) then
           call cpu_time(tBefore)
           ih=0 ! no summing over helicities
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,(next+1)/2
                 amp2=0d0
                 do ic=amplitudes_full%row_index(irow-1,i)+1,amplitudes_full%row_index(irow,i)
                    icol=amplitudes_full%col_index(ic,i)
                    amp2=amp2+amplitudes_NLC%amps(amplitudes_NLC%helmap(icol,ih),icol)
                 enddo
                 if (i.eq.1) then
                    amp_col=amp_col+amp2*amplitudes_full%col_value(i)
                 else
                    amp_col=amp_col+amp2*amplitudes_full%col_value(i)*2d0
                 endif
              enddo
              amp2_full(k)=amp2_full(k)+amp_col*amplitudes_NLC%amps(amplitudes_NLC%helmap(irow,ih),irow)
           enddo
           call cpu_time(tAfter)
           t_mat_full=t_mat_full+tAfter-tBefore
        elseif ((reweight_mode.eq.16 .or. reweight_mode.eq.17) .and. amplitudes_cache%col_acc.ge.2) then
           call cpu_time(tBefore)
           do irow=1,factorial(next-1)
              amp_col=0d0
              do i=1,(next+1)/2
                 amp2=0d0
                 do ic=amplitudes_cache%row_index_full(irow-1,i)+1,amplitudes_cache%row_index_full(irow,i)
                    icol=amplitudes_cache%col_index_full(ic,i)
                    amp2=amp2+amplitudes_cache%amps(icol)
                 enddo
                 if (i.eq.1) then
                    amp_col=amp_col+amp2*amplitudes_cache%col_value_full(i)
                 else
                    amp_col=amp_col+amp2*amplitudes_cache%col_value_full(i)*2d0
                 endif
              enddo
              amp2_full(k)=amp2_full(k)+amp_col*amplitudes_cache%amps(irow)
           enddo
           call cpu_time(tAfter)
           t_mat_full=t_mat_full+tAfter-tBefore
        elseif (reweight_mode.ge.2 .and. reweight_mode.le.9) then
           call cpu_time(tBefore)
           do ih=0,amplitudes_NLC%nhel(amplitudes_NLC%isize+1)-1
              amp2_hel=0d0
              do icol=1,list(0)
                 amp2_hel=amp2_hel+amplitudes_NLC%amps(amplitudes_NLC%helmap(list(icol),ih),list(icol))
              enddo
              amp2_full(k)=amp2_full(k)+amp2_hel**2
           enddo
           amp2_full(k)=amp2_full(k) * color_wgt
           call cpu_time(tAfter)
           t_mat_full=t_mat_full+tAfter-tBefore
        endif

     enddo

!!$     do i=1,factorial(next-1)
!!$        do j=1,factorial(next-1)
!!$           if (any(list(1:list(0)).eq.i .and. any(list(1:list(0)).eq.j))) icount(i,j)=icount(i,j)+color_wgt
!!$        enddo
!!$     enddo

     call write_event(12)
  enddo

  if (reweight_mode.eq.3 .or. reweight_mode.eq.5  .or. reweight_mode.eq.7) then
     write (*,*) 'factor to reweight event weights by:', ievt_count(1)/ievt_count(2)
     write (12,*) 'REWEIGHT_FACTOR', ievt_count(1)/ievt_count(2)
  endif

!!$  do i=1,1!factorial(next-1)
!!$     write (*,*) icount(i,1:factorial(next-1))/ievt_count(2)
!!$  enddo

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
    open(unit=11,file='events'//tag//'.lhe',status='old')
    open(unit=12,file='events'//tag//'.lhe.rwgt',status='unknown')
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
       read (iunit,*,err=99,end=99) dum,p(1:3,i),p(0,i)
    enddo
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event


  subroutine write_event(iunit)
    use rw_events
    implicit none
    integer :: i,iunit
!!$    if (reweight_mode.eq.14) then
!!$       rwgt_NLC=0d0
!!$       rwgt_full=sum(amp2_full(1:n_rwgt))/sum(amp2_LC(1:n_rwgt))
!!$    else
       rwgt_NLC=sum(amp2_NLC(1:n_rwgt)/amp2_LC(1:n_rwgt))/dble(n_rwgt)
       rwgt_full=sum(amp2_full(1:n_rwgt)/amp2_LC(1:n_rwgt))/dble(n_rwgt)
!!$    endif
    write (iunit,*) '<event>'
    write (iunit,*) next,evt_wgt,wgt,amp2,weight
    write (iunit,'(100i3)') o(1:next)
    write (iunit,*) rwgt_full,sum(amp2_LC(1:n_rwgt)),sum(amp2_full(1:n_rwgt))
    write (iunit,*) evt_wgt,evt_wgt*rwgt_NLC,evt_wgt*rwgt_full
    if (sum_hel) then
       write (iunit,'(i3)') 99
    else
       write (iunit,'(100i3)') hel(1:next)
    endif
    do i=1,next
       if (i.le.2) then
          write (iunit,*) ' 21',p(1:3,i),p(0,i)
       else
          write (iunit,*) ' 21',p(1:3,i),p(0,i)
       endif
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event

  subroutine pick_random_color(n)
    ! Randomly selects a colour assignment. For some 'reweight_mode's this
    ! uses information from the 'init_color_probs()' subroutine.
    use random_colours
    implicit none
    integer,intent(in) ::n
    integer :: i,j,nc
    integer,dimension(1:3) :: nrc
    integer,dimension(1:n) :: ida
    real(kind=8) :: random
    integer,save :: col_picked=0
    if (.not. allocated(col_labels)) allocate ( col_labels(1:n,2) )
    if (reweight_mode.eq.8) then
       col_picked=int(ran2()*3**n)+1
       col_labels(1:n,1:2)=color_labels(1:n,1:2,col_picked)
       color_wgt=1d0
    elseif (reweight_mode.eq.14) then
       col_picked=mod(col_picked,(3**(n-1)-1)/2)+1
       col_labels(1:n,1:2)=color_labels(1:n,1:2,col_picked)
       color_wgt=dble((3**(n-1)-1))/dble(3**(n-1))
!!$       color_wgt=1d0
    elseif (reweight_mode.eq.9 .or. reweight_mode.eq.10) then
       color_wgt=1d0
       do
          do i=1,n
             col_labels(i,1)=int(ran2()*3)+1
             col_labels(i,2)=int(ran2()*3)+1
          enddo
          if (is_ok_order(n)) then
             if (check_permutation(iperm_ev,n-1)) return
          endif
       enddo
    elseif (reweight_mode.eq.11 .or. reweight_mode.eq.12) then
       color_wgt=1d0
       do i=1,n
          col_labels(i,1)=int(ran2()*3)+1
       enddo
       call iperm_decode(iperm_ev,n-1,ida)
       ida(n)=n
       do i=1,n-1
          col_labels(ida(i+1),2)=col_labels(ida(i),1)
       enddo
       col_labels(ida(1),2)=col_labels(ida(n),1)
       return
    elseif (reweight_mode.eq.13) then
       color_wgt=1d0
!!$       nrc(1:3)=(/ 4,4,0 /)
!!$       do i=1,n
!!$          nc=int(ran2()*(n+1-i))
!!$          if (nc.lt.nrc(1)) then
!!$             nrc(1)=nrc(1)-1
!!$             col_labels(i,1)=1
!!$          elseif (nc.lt.nrc(1)+nrc(2)) then
!!$             nrc(2)=nrc(2)-1
!!$             col_labels(i,1)=2
!!$          else
!!$             nrc(3)=nrc(3)-1
!!$             col_labels(i,1)=3
!!$          endif
!!$       enddo
       col_labels(1:n,1)= (/ 1,1,1,1,1,1,2,2 /)
       call iperm_decode(iperm_ev,n-1,ida)
       ida(n)=n
       do i=1,n-1
          col_labels(ida(i+1),2)=col_labels(ida(i),1)
       enddo
       col_labels(ida(1),2)=col_labels(ida(n),1)
       return
    else
       do
          if (reweight_mode.eq.2 .or. reweight_mode.eq.3) then
             do i=1,n
                col_labels(i,1)=int(ran2()*3)+1
                col_labels(i,2)=int(ran2()*3)+1
             enddo
             color_wgt=3**(next*2)

          elseif(reweight_mode.ge.4 .and. reweight_mode.le.8) then
             ! pick a random color assignment, weighted by
             ! accum_color_probs
             random=ran2()
             do i=1,n_col
                if (accum_color_probs(i).gt.random) then
                   col_picked=i
                   exit
                endif
             enddo

             do j=1,2
                nrc(1:3)=n_cols(col_picked,1:3)
                do i=1,n
                   nc=int(ran2()*(n+1-i))
                   if (nc.lt.nrc(1)) then
                      nrc(1)=nrc(1)-1
                      col_labels(i,j)=1
                   elseif (nc.lt.nrc(1)+nrc(2)) then
                      nrc(2)=nrc(2)-1
                      col_labels(i,j)=2
                   else
                      nrc(3)=nrc(3)-1
                      col_labels(i,j)=3
                   endif
                enddo
             enddo

             if (reweight_mode.eq.6 .or. reweight_mode.eq.7) then
                color_wgt=dble(n_col)
                if ( n_cols(col_picked,1).eq.n_cols(col_picked,2) .and. &
                     n_cols(col_picked,1).eq.n_cols(col_picked,3) ) then
                   color_wgt=color_wgt*1d0
                elseif ( n_cols(col_picked,1).eq.n_cols(col_picked,2) .or. &
                     n_cols(col_picked,1).eq.n_cols(col_picked,3) .or. &
                     n_cols(col_picked,2).eq.n_cols(col_picked,3) ) then
                   color_wgt=color_wgt*3d0
                else
                   color_wgt=color_wgt*6d0
                endif
                color_wgt=color_wgt*(factorial_dble(next)/(factorial_dble(n_cols(col_picked,1))* &
                     factorial_dble(n_cols(col_picked,2))*factorial_dble(n_cols(col_picked,3))))**2
             else
                color_wgt=factor_wgt
             endif
          endif
          if (reweight_mode.eq.3 .or. reweight_mode.eq.5) then
             ievt_count(2)=ievt_count(2)+1d0
             if (is_ok_order(n)) then
                ievt_count(1)=ievt_count(1)+1d0
                return
             endif
          elseif (reweight_mode.eq.7) then
             ievt_count(2)=ievt_count(2)+color_wgt/factor_wgt
             if (is_ok_order(n)) then
                ievt_count(1)=ievt_count(1)+1d0
                return
             endif
          else
             ievt_count(1)=ievt_count(1)+1d0
             ievt_count(2)=ievt_count(2)+1d0
             return
          endif
       enddo
    endif
  end subroutine pick_random_color

  subroutine init_color_probs(n)
    ! Initialises all possible colour assignments with appropriate weights.
    use random_colours
    implicit none
    integer,intent(in) :: n
    real(kind=8) :: temp,factor
    integer,dimension(1:3) :: itmp
    integer,dimension(1:n) :: ida
    integer :: i,j,jj
    if (reweight_mode.eq.8) then
       if (.not. allocated(color_labels)) allocate(color_labels(1:n,2,3**n))
       if (.not. allocated(col_labels)) allocate(col_labels(1:n,2))
       jj=0
       do i=1,3**(2*n)
          do j=1,2*n
             if (j.le.n) then
                col_labels(j,1)=mod((i-1)/(3**(2*n-j)),3)+1
             else
                col_labels(j-n,2)=mod((i-1)/(3**(2*n-j)),3)+1
             endif
          enddo
          if (is_ok_order(n)) then
             if (check_permutation(iperm_ev,n-1)) then
                jj=jj+1
                color_labels(1:n,1,jj)=col_labels(1:n,1)
                color_labels(1:n,2,jj)=col_labels(1:n,2)
             endif
          endif
       enddo
    elseif(reweight_mode.eq.14) then
       if (.not. allocated(color_labels)) allocate(color_labels(1:n,2,(3**(n-1)-1)/2))
       if (.not. allocated(col_labels)) allocate(col_labels(1:n,2))
       jj=0
       colourloop: do i=1,2*3**(n-2)
          do j=1,n
             if (j.eq.1) then
                col_labels(j,1)=1
             elseif (j.eq.2) then
                if (i.gt.3**(n-2)) then
                   col_labels(j,1)=2
                else
                   col_labels(j,1)=1
                endif
             else
                col_labels(j,1)=mod((i-1)/(3**(n-j)),3)+1
             endif
          enddo
          if (all(col_labels(1:n,1).eq.1)) cycle colourloop
          do j=1,n
             if (col_labels(j,1).eq.2) exit
             if (col_labels(j,1).eq.3) cycle colourloop
          enddo
          call iperm_decode(iperm_ev,n-1,ida)
          ida(n)=n
          do j=1,n-1
             col_labels(ida(j+1),2)=col_labels(ida(j),1)
          enddo
          col_labels(ida(1),2)=col_labels(ida(n),1)
          jj=jj+1
          color_labels(1:n,1,jj)=col_labels(1:n,1)
          color_labels(1:n,2,jj)=col_labels(1:n,2)
       enddo colourloop
    else
       ! Count how many we have so that we can allocated relevant arrays
       n_col=0
       do j=0,n/3
          do jj=j,(n-j)/2
             n_col=n_col+1
          enddo
       enddo
       allocate(n_cols(n_col,1:3))
       allocate(accum_color_probs(0:n_col))
       accum_color_probs(0)=0d0
       n_col=0
       factor_wgt=0d0
       do j=0,n/3
          do jj=j,(n-j)/2
             n_col=n_col+1
             n_cols(n_col,1)=n-j-jj
             n_cols(n_col,2)=jj
             n_cols(n_col,3)=j
             if (reweight_mode.ne.6) then
                ! 'factor' takes into account that we only consider
                ! n_cols(,1)>=n_cols(,2)>=n_cols(,3), which means we have to
                ! increase the relative probability if they are all different
                ! by a factor 6 as compared to them being all identical (and
                ! by a factor 3 if two out of the three are identical)
                if ( n_cols(n_col,1).eq.n_cols(n_col,2) .and. &
                     n_cols(n_col,1).eq.n_cols(n_col,3) ) then
                   factor=1d0
                elseif ( n_cols(n_col,1).eq.n_cols(n_col,2) .or. &
                     n_cols(n_col,1).eq.n_cols(n_col,3) .or. &
                     n_cols(n_col,2).eq.n_cols(n_col,3) ) then
                   factor=3d0
                else
                   factor=6d0
                endif
                ! needs **2, since this applies both to the color and
                ! anti-color indices of the gluons
                factor=factor*(factorial_dble(next)/(factorial_dble(n_cols(n_col,1))* &
                     factorial_dble(n_cols(n_col,2))*factorial_dble(n_cols(n_col,3))))**2
                factor_wgt=factor_wgt+factor
             endif
             if (reweight_mode.le.5) then
                accum_color_probs(n_col)=factor
             elseif (reweight_mode.eq.6 .or. reweight_mode.eq.7) then
                accum_color_probs(n_col)=1d0
             endif
          enddo
       enddo
       ! order them according to their probabilty (largest comes first)
       do j=1,n_col-1
          do jj=j+1,n_col
             if (accum_color_probs(j).lt.accum_color_probs(jj)) then
                temp=accum_color_probs(j)
                accum_color_probs(j)=accum_color_probs(jj)
                accum_color_probs(jj)=temp
                itmp(1:3)=n_cols(j,1:3)
                n_cols(j,1:3)=n_cols(jj,1:3)
                n_cols(jj,1:3)=itmp(1:3)
             endif
          enddo
       enddo
       ! convert probs to accumulated probs
       do j=1,n_col
          accum_color_probs(j)=accum_color_probs(j-1)+accum_color_probs(j)
       enddo
       accum_color_probs(1:n_col)=accum_color_probs(1:n_col)/accum_color_probs(n_col)
    endif
  end subroutine init_color_probs

  logical function is_ok_order(n)
    ! Checks if the colour assignment saved in 'col_labels()' is a reasonable
    ! colour assignment that will result in at least one colour order that
    ! contributes.
    use random_colours
    implicit none
    integer :: n,i,icol,n_color,n_anticolor
    if (reweight_mode.eq.2 .or. reweight_mode.eq.3 .or. reweight_mode.eq.8 .or. reweight_mode.eq.9) then
       ! checks if number of red, green and blue is identical in color,
       ! 'col_label(,1)', and anti-color, 'col_label(,2)', indices
       do icol=1,3 ! red, green and blue
          n_color=count(col_labels(1:n,1).eq.icol)
          n_anticolor=count(col_labels(1:n,2).eq.icol)
          if (n_color.ne.n_anticolor) then
             is_ok_order=.false.
             return
          endif
       enddo
    endif
    is_ok_order=.true.
    do icol=1,3
       if (all(col_labels(1:n,1).ne.icol)) cycle ! deal with special case
       if (all(col_labels(1:n,1).eq.icol)) cycle ! deal with special case
       is_ok_order=.false.
       do i=1,n
          if (col_labels(i,1).eq.icol .and. col_labels(i,2).ne.icol) then
             is_ok_order=.true.
             exit
          endif
       enddo
       if (.not.is_ok_order) return
    enddo
  end function is_ok_order

  logical function check_permutation(iperm,n)
    ! Checks if permutation 'iperm' contributes to colour assignment saved in
    ! 'col_labels()'.
    use random_colours
    implicit none
    integer :: iperm,n
    integer,dimension(n+1) :: ida
    call iperm_decode(iperm,n,ida)
    ida(n+1)=n+1
    check_permutation=.true.
    do i=1,n
       if (col_labels(ida(i),1).ne.col_labels(ida(i+1),2)) then
          check_permutation=.false.
          return
       endif
    enddo
    if (col_labels(ida(n+1),1).ne.col_labels(ida(1),2)) check_permutation=.false.
  end function check_permutation

  subroutine find_permutations(n)
    ! Finds all permutations that contribute to colour assignment saved in
    ! 'col_labels()'. The number of contributing permutations is saved in
    ! 'n_valid' and the permutations themselves in
    ! 'list_valid_iperm(1:n_valid)'.
    use random_colours
    implicit none
    integer :: n,ipos
    integer,dimension(1:n) :: sigma
    logical,dimension(1:n) :: free
    n_valid=0
    free=.true.
    ipos=1
    sigma(1:n)=0
    call check_next_pos(ipos,col_labels(n+1,2),sigma,n,free)
  end subroutine find_permutations

  recursive subroutine check_next_pos(ipos,jj,sigma,n,free)
    ! Works together with the subroutine 'find_permutations()'
    use random_colours
    implicit none
    integer :: jj,n,i,iperm
    integer :: ipos,ipos1
    logical,dimension(1:n) :: free,free1
    logical :: foundone
    integer,dimension(1:n) :: sigma,sigma1
    if (ipos.eq.n+1) then
       if (col_labels(n+1,1).eq.jj) then
          n_valid=n_valid+1
          if (reweight_mode.ge.10 .and. reweight_mode.le.14) then
             list_orders(1:n,n_valid)=sigma(1:n)
             list_orders(n+1,n_valid)=n+1
          else
             call iperm_encode(n,sigma,iperm)
             list_valid_iperm(n_valid)=iperm
          endif
       else
          ipos=ipos-1
       endif
       return
    endif
    foundone=.false.
    do i=1,n
       if (free(i).and.col_labels(i,1).eq.jj) then
          free1(1:n)=free(1:n)
          foundone=.true.
          free1(i)=.false.
          sigma1=sigma
          sigma1(ipos)=i
          ipos1=ipos+1
          call check_next_pos(ipos1,col_labels(i,2),sigma1,n,free1)
       endif
    enddo
    if (.not. foundone) then
       ipos=ipos-1
    endif
  end subroutine check_next_pos
end program matrix_reweight
