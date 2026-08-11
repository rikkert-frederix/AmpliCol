! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module rw_events
  implicit none
  type :: lhef_event
!!$     real(kind=8) :: wgt,evt_wgt,weight,rwgt_NLC,rwgt_full
!!$     integer,dimension(:),allocatable :: helicity,col_order,iPDG
!!$     real(kind=8),dimension(:,:),allocatable :: momenta
     integer :: NUP,IDPRUP
     integer,dimension(:),allocatable :: IDUP,ISTUP
     integer,dimension(:,:),allocatable :: MOTHUP,ICOLUP
     real(kind=8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP
     real(kind=8),dimension(:,:),allocatable :: PUP
     real(kind=8),dimension(:),allocatable :: VTIMUP,SPINUP
     real(kind=8),dimension(3) :: overwgt
     integer,dimension(:),allocatable :: col_order
     real(kind=8),dimension(3) :: matrix2
  end type lhef_event
  type(lhef_event),dimension(:),allocatable :: events
  integer :: nevents
  integer :: IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP
  real(kind=8) :: EBMUP(2),XSECUP,XERRUP,XMAXUP
  character(len=1024) :: generator_string
  logical :: unwgt,keep_comments
end module rw_events
module timings
  implicit none
  real(kind=4) :: tBefore,tAfter,tTot_A=0.,tTot_B=0.,t_amp=0.,t_amp_init=0.,&
       t_mat_LC=0.,t_mat_NLC=0.,t_mat_full=0.,t_all=0.,t_ran=0.
end module timings
module overall
  real(kind=8) :: xsec,xsec_abs,max_wgt
  integer :: nevt
end module overall
module arguments
  implicit none
  integer :: c_o,c_o_t,c_o_i,c_o_j,c_o_k,imode
end module arguments

program amplicol_reweight
  use rw_events
  use math_functions
  use amplitude_QCD_mod
  use timings
  use particles
  use run_parameters
  use overall
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
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(:),allocatable :: unique_map_value
  complex(kind=8) :: amp2_c,amp_col_c
  logical :: done
  character(len=256) :: event_filename,input_filename
  
  call get_run_arguments()
  call read_run_parameters(input_filename)

  call cpu_time(tTot_B)

  call phys_model%init_part()
  open(unit=11,file=trim(event_filename),status='old')
  call read_unique_in_file()
  call allocate_process_info()
  open(unit=12,file=trim(event_filename)//'.rwgt',status='unknown')
  call read_init_and_allocate_events(11)
  do nevt=1,nevents
     call read_event(11,events(nevt))
  enddo
  call apply_final_state_widths_from_events()
  call phys_model%init_vert()

  nprocs=0
  call setup_spin()
  col_acc=20

  xsec=0d0
  xsec_abs=0d0
  max_wgt=0d0
  do nevt=1,nevents
     call map_to_canonical_form(events(nevt))
     do iproc=1,nprocs
        if (all(part(1:next,1).eq.processes(1:next,iproc)))exit
     enddo
     if (iproc.eq.nprocs+1) then
        call cpu_time(tBefore)
        nprocs=nprocs+1
        processes(1:next,iproc)=part(1:next,1)
        call amps(iproc)%init(2,next,1,part,spin,o,phys_model)
        call amps(iproc)%init_col(next,col_acc)
        call cpu_time(tAfter)
        t_amp_init=t_amp_init+tAfter-tBefore
     endif
     events(nevt)%matrix2(1:3)=0d0

     call cpu_time(tBefore)
    
     call amps(iproc)%evaluate(next,p,hel,read_proc_from_file,phys_model)

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
              events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)+amp_col*amps(iproc)%amps_r(ioff+irow)
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
              events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)+dble(amp_col_c*conjg(amps(iproc)%amps(ioff+irow)))
           enddo
        endif
        ! The following line can be removed since 'process_map_value'
        ! will drop out when taking the ratio w.r.t. LC.
!        events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)*process_map_value
        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo
     max_wgt=max(max_wgt,events(nevt)%matrix2(3)/events(nevt)%matrix2(1))
     xsec=xsec+events(nevt)%matrix2(3)/events(nevt)%matrix2(1)*events(nevt)%XWGTUP
     xsec_abs=xsec_abs+abs(events(nevt)%matrix2(3)/events(nevt)%matrix2(1)*events(nevt)%XWGTUP)
  enddo
  xsec=xsec/dble(nevents)
  xsec_abs=xsec_abs/dble(nevents)
  close(11)

  ! write event file
  call write_init(12)
  do nevt=1,nevents
     call write_event(12,events(nevt))
  enddo
  write(12,'(a)') '</LesHouchesEvents>'
  
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
  write(*,*) 'Total FC cross section:',xsec
contains  

  subroutine get_run_arguments()
    use arguments
    implicit none
    integer :: argc,i,k
    character(len=256) :: argv
    logical :: show_help
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()

    show_help=.false.
    unwgt=.false.
    keep_comments=.true.
    event_filename=''
    input_filename='run_card.dat'
    
    if (argc.lt.1) then
       show_help=.true.
    else
       call get_command_argument(1,argv)
       do i=1,argc
          call get_command_argument(i, argv)
          argv = trim(argv)
          if (index(argv,"--help").eq.1 .or. index(argv,"-h").eq.1) then
             show_help = .true.
          elseif (index(argv, "--unwgt").eq.1) then
             unwgt=.true.
          elseif (index(argv, "--remove_comments").eq.1) then
             keep_comments=.false.
          elseif (index(argv, "--input=").eq.1 .or. index(argv, "--card=").eq.1) then
             input_filename=argv(index(argv,"=")+1:)
          elseif (index(argv,"-").eq.1) then
             write (*,*) 'Unknown argument: ',trim(argv)
             stop 1
          else
             if (len_trim(event_filename).ne.0) then
                write (*,*) 'More than one event file was specified'
                stop 1
             endif
             event_filename=argv
          endif
       enddo
    endif
    if (len_trim(event_filename).eq.0) show_help=.true.

    if (show_help) then
       write (*,'(a)') ""
       write (*,'(a)') "Usage: 'amplicol_reweight <event_file> <arguments>', where"//&
              " <event_file> is the leading colour event file to reweight to full colour."
       write (*,'(a)') "The code creates an LHEF, '<event_file>.rwgt' containing the full colour events."
       write (*,'(a)') "Possible arguments are"
       write (*,'(a)') ""
       write (*,'(a)') "  --help,   -h      : Show this message."
       write (*,'(a)') "  --unwgt           : Unweight the reweight events."
       write (*,'(a)') "  --remove_comments : Remove comment lines in the final LHEF."
       write (*,'(a)') "  --input=[X]       : Physics/run input card (default is './run_card.dat')."
       write (*,'(a)') ""
       stop
    endif

  end subroutine get_run_arguments

  subroutine apply_final_state_widths_from_events()
    implicit none
    integer :: ievent,i
    integer,dimension(:,:),allocatable :: event_process
    if (ignore_final_state_width_fix) return
    if (nevents.eq.0) return
    allocate(event_process(events(1)%NUP,nevents))
    do ievent=1,nevents
       if (events(ievent)%NUP.ne.events(1)%NUP) then
          write (*,*) 'Cannot combine LHE events with different particle multiplicities'
          stop 1
       endif
       event_process(:,ievent)=21
       event_process(1:2,ievent)=events(ievent)%IDUP(1:2)
       do i=3,events(ievent)%NUP
          if (events(ievent)%ISTUP(i).eq.1) &
               event_process(i,ievent)=events(ievent)%IDUP(i)
       enddo
    enddo
    call phys_model%apply_final_state_widths(events(1)%NUP,nevents,event_process)
    deallocate(event_process)
  end subroutine apply_final_state_widths_from_events

  subroutine allocate_process_info()
    implicit none
    if (.not. allocated(hel)) then
!!$       allocate(momenta(0:3,1:next))
!!$       allocate(helicity(1:next))
!!$       allocate(col_order(1:next))
!!$       allocate(iPDG(1:next))
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
    character :: dummy
    read(11,*) dummy ! <LesHouchesEvents>-tag
    read(11,*) dummy ! <header>-tag
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

  subroutine read_event(iunit,event)
    implicit none
    type(lhef_event) :: event
    integer :: i,iunit
    character(len=1024) :: string
    real(kind=8) :: dum
    read (iunit,'(a)') string ! <event>-tag
    read(iunit,*) event%NUP,event%IDPRUP,event%XWGTUP,event%SCALUP,event%AQEDUP,event%AQCDUP
    allocate(event%IDUP(event%NUP),event%ISTUP(event%NUP),event%MOTHUP(2,event%NUP)&
         &,event%ICOLUP(2,event%NUP),event%PUP(5,event%NUP),event%VTIMUP(event%NUP)&
         &,event%SPINUP(event%NUP),event%col_order(event%NUP))
    do i=1,event%NUP
       read(iunit,*) event%IDUP(I),event%ISTUP(I),event%MOTHUP(1,I),event%MOTHUP(2,I), &
            event%ICOLUP(1,I),event%ICOLUP(2,I),&
            event%PUP(1,I),event%PUP(2,I),event%PUP(3,I),event%PUP(4,I),event%PUP(5,I),&
            event%VTIMUP(I),event%SPINUP(I)
    enddo
    read(iunit,'(a)') string
    read(string(7:),*) event%col_order(1:event%NUP)
    read(iunit,'(a)') string
    read(string(9:),*) event%overwgt(1:3)

!!$    read (iunit,*,err=99,end=99) next,evt_wgt!,wgt,amp2,weight
!!$    read (iunit,*,err=99,end=99) helicity(1:next)
!!$    read (iunit,*,err=99,end=99) col_order(1:next)
!!$    do i=1,next
!!$       read (iunit,*,err=99,end=99) iPDG(i),momenta(1:3,i),momenta(0,i)
!!$    enddo
    read (iunit,'(a)') string  ! </event>-tag
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
    ! 'array_map' maps the PDG codes to the 'sort_particle' codes
    ! (see Utilities/process_list.py)
    integer,dimension(-24:25),parameter :: array_map=[83,0,0,0,0,0,0&
         &,0,95,88,93,86,91,84,0,0,0,0,12,11,10,9,8,7,0,1,2,3,4,5,6,0,0&
         &,0,0,85,90,87,92,89,94,0,0,0,0,13,80,81,82,96]
    ! Initialize mapping
    mapping = [(i,i=1,n)]
    ! Sort the array and mapping using a simple bubble sort
    do i=1,n-1
       do j=1,n-i
          if (array_map(array(j)) .gt. array_map(array(j+1))) then
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

  
  subroutine map_to_canonical_form(event)
    ! cross the two initial state particle PDGs, order according to
    ! the PDG value, (and reflip the two initial states again)
    implicit none
    type(lhef_event) :: event
    integer,dimension(next) :: mapping
    real(kind=8),dimension(0:3,next) :: p_cross
    integer :: i,iproc,j
    logical :: compatible_mass_layout
    next=event%NUP
    part(1:next,1)=event%IDUP(1:next)
    if (.not.use_only_canonical_form) then
       ! do no use the mapping to canonical form, but reweight the
       ! events as they are.
       hel(1:next)=nint(event%SPINUP(1:next))
       o(1:next,1)=event%col_order(1:next)
       p(1:3,1:next)=event%PUP(1:3,1:next)
       p(0,1:next)=event%PUP(4,1:next)
    else
       ! Map to canonical from to reduce the number of matrix elements
       ! to initialise.
       ! cross the initial state
       part(1,1)=phys_model%get_antipart(part(1,1))
       part(2,1)=phys_model%get_antipart(part(2,1))
       p_cross(1:3,1:2)=-event%PUP(1:3,1:2)
       p_cross(0,1:2)=-event%PUP(4,1:2)
       p_cross(1:3,3:next)=event%PUP(1:3,3:next)
       p_cross(0,3:next)=event%PUP(4,3:next)
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
          hel(i)=nint(event%SPINUP(mapping(i)))
          o(i,1)=event%col_order(mapping(i)) ! this is not correct, but isn't used
       enddo
       ! Convert to 'unique flavour configuration' (if available)
       if (unique_nproc.eq.0) return
       do iproc=1,unique_nproc
          if (all(part(1:next,1).eq.unique_processes(1:next,iproc))) then
             process_map_value=unique_map_value(iproc)
             if (unique_map(iproc).gt.0) then
                compatible_mass_layout=.true.
                do j=1,next
                   if (phys_model%get_mass(unique_processes(j,iproc)).ne.&
                        phys_model%get_mass(unique_processes(j,unique_map(iproc)))) then
                      compatible_mass_layout=.false.
                      exit
                   endif
                enddo
                if (compatible_mass_layout) then
                   part(1:next,1)=unique_processes(1:next,unique_map(iproc))
                else
                   ! Backward compatibility for LHE files written before
                   ! massive-leg layouts were excluded from flavour maps.
                   process_map_value=1d0
                endif
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

  subroutine write_event(iunit,event)
    use overall
    implicit none
    type(lhef_event),intent(in) :: event
    integer :: iunit
    real(kind=8) :: rwgt_NLC,rwgt_full,XWGTUP_new
    real(kind=8),external :: ran2
    rwgt_NLC=event%matrix2(2)/event%matrix2(1)
    rwgt_full=event%matrix2(3)/event%matrix2(1)
    if (unwgt) then
       if (rwgt_full.gt.max_wgt*ran2()) then
          XWGTUP_new=xsec_abs
       else
          return
       endif
    else
       XWGTUP_new=event%XWGTUP*rwgt_full
    endif
    write (iunit,'(a)') '<event>'
    write(iunit,503)event%NUP,event%IDPRUP,XWGTUP_new,event%SCALUP,event%AQEDUP,event%AQCDUP
    do i=1,event%NUP
       write(iunit,504)event%IDUP(I),event%ISTUP(I),event%MOTHUP(1,I),event%MOTHUP(2,I), &
            event%ICOLUP(1,I),event%ICOLUP(2,I),&
            event%PUP(1,I),event%PUP(2,I),event%PUP(3,I),event%PUP(4,I),event%PUP(5,I),&
            event%VTIMUP(I),event%SPINUP(I)
    enddo
    if (keep_comments) then
       write (iunit,506) '#color',event%col_order(1:event%NUP)
       write (iunit,505) '#overwgt',event%overwgt(1:3)
       write (iunit,505) '#color_expansion',event%XWGTUP,event%XWGTUP*rwgt_NLC,event%XWGTUP*rwgt_full
    endif
    write (iunit,'(a)') '</event>'
503 format(1x,i2,1x,i6,4(1x,e14.8))
504 format(1x,i8,1x,i2,4(1x,i4),5(1x,e24.17),2(1x,e10.4))
505 format(a,3(1x,e14.8))
506 format(a,100i3)
  end subroutine write_event

  subroutine write_init(ounit)
    use overall
    implicit none
    integer,intent(in) :: ounit
    real(kind=8) :: rwgt_LCtoFC
    rwgt_LCtoFC=xsec/XSECUP
    write(ounit,'(a)') '<LesHouchesEvents version="3.0">'
    write(ounit,'(a)') '<init>'
    write(ounit,501)IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1)&
         &,PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
    if (unwgt) then
       write(ounit,502)xsec,XERRUP*rwgt_LCtoFC,xsec_abs,-3
    else
       write(ounit,502)xsec,XERRUP*rwgt_LCtoFC,XMAXUP*max_wgt,-4
    endif
    write(ounit,'(a)') trim(generator_string)
    write(ounit,'(a)') '</init>'
501 format(2(1x,i6),2(1x,e14.8),2(1x,i2),2(1x,i8),1x,i2,1x,i3)
502 format(3(1x,e14.8),1x,i6)
  end subroutine write_init

  subroutine read_init_and_allocate_events(iunit)
    implicit none
    integer,intent(in) :: iunit
    character(len=1024) :: string
    integer(kind=8) iseed
    common /to_seed/iseed
    do
       read(iunit,'(a)') string
       if (index(string,"<nevents>").ne.0) then
          read(string(10:),*) nevents
       elseif (index(string,"<seed>").ne.0) then
          read(string(10:),*) iseed
       endif
       if (index(string,"<init>").ne.0) then
          read(iunit,*) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1) &
               &,PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
          read(iunit,*) XSECUP,XERRUP,XMAXUP,LPRUP
          read(iunit,'(a)') generator_string
       endif
       if (index(string,"</init>").ne.0) exit
    enddo
    allocate(events(nevents))
  end subroutine read_init_and_allocate_events
  
end program amplicol_reweight
