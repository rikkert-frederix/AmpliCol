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
     logical :: keep=.true.
  end type lhef_event
  type(lhef_event),dimension(:),allocatable :: events
  integer :: nevents
  integer :: IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP
  real(kind=8) :: EBMUP(2),XSECUP,XERRUP,XMAXUP
  character(len=1024) :: generator_string
  logical :: unwgt,keep_comments
  integer(kind=8) :: input_seed=0_8
  logical :: has_input_seed=.false.
end module rw_events
module timings
  implicit none
  real(kind=4) :: tBefore,tAfter,tTot_A=0.,tTot_B=0.,t_amp=0.,t_amp_init=0.,&
       t_mat_LC=0.,t_mat_NLC=0.,t_mat_full=0.,t_all=0.,t_ran=0.
end module timings
module overall
  real(kind=8) :: xsec,xsec_abs,xsec_error,xsec_sumw2,max_wgt,max_abs_event_weight
  integer :: nevt,nevents_output
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
  use coupling_orders
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
  call resolve_legacy_automatic_selection()

  nprocs=0
  call setup_spin()
  col_acc=20

  xsec=0d0
  xsec_abs=0d0
  xsec_sumw2=0d0
  max_wgt=0d0
  max_abs_event_weight=0d0
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
        call prune_selected_coupling_sectors(amps(iproc))
        call cpu_time(tAfter)
        t_amp_init=t_amp_init+tAfter-tBefore
     endif
     events(nevt)%matrix2(1:3)=0d0

     call cpu_time(tBefore)
    
     call amps(iproc)%evaluate(next,p,hel,read_proc_from_file,phys_model)

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore

     do iacc=1,3 ! LC, NLC and full colour
        call cpu_time(tBefore)
        if (iacc.eq.3 .and. col_acc.lt.2) cycle
        events(nevt)%matrix2(iacc)=selected_colour_matrix2(&
             amps(iproc),iacc,events(nevt)%AQCDUP,events(nevt)%AQEDUP)
        ! The following line can be removed since 'process_map_value'
        ! will drop out when taking the ratio w.r.t. LC.
!        events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)*process_map_value
        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo
     amp2=reweight_ratio(events(nevt)%matrix2(3),events(nevt)%matrix2(1),&
          'full-colour/leading-colour')
     max_wgt=max(max_wgt,abs(amp2))
     max_abs_event_weight=max(max_abs_event_weight,&
          abs(amp2*events(nevt)%XWGTUP))
     xsec=xsec+amp2*events(nevt)%XWGTUP
     xsec_abs=xsec_abs+abs(amp2*events(nevt)%XWGTUP)
     xsec_sumw2=xsec_sumw2+(amp2*events(nevt)%XWGTUP)**2
  enddo
  xsec=xsec/dble(nevents)
  xsec_abs=xsec_abs/dble(nevents)
  if (nevents.gt.1) then
     xsec_error=sqrt(max(0d0,(xsec_sumw2/dble(nevents)-xsec**2)/&
          dble(nevents-1)))
  else
     xsec_error=0d0
  endif
  call prepare_output_selection()
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

  subroutine prune_selected_coupling_sectors(amplitude)
    implicit none
    type(amplitude_QCD),intent(inout) :: amplitude
    logical,dimension(:,:),allocatable :: allowed_sector_pairs
    logical :: valid_pair_mask

    call build_coupling_order_pair_mask(amplitude%sector_powers,&
         allowed_sector_pairs,valid_pair_mask)
    if (.not.valid_pair_mask) then
       write (*,*) 'Coupling-order selection must be resolved before amplitude pruning'
       stop 1
    endif
    call amplitude%prune_coupling_sectors(allowed_sector_pairs)
    deallocate(allowed_sector_pairs)
  end subroutine prune_selected_coupling_sectors

  subroutine resolve_legacy_automatic_selection()
    implicit none
    integer :: iprocess,ileg,ievent,nsing,max_as2
    logical :: valid_selection

    if (coupling_order_selection%mode.ne.coupling_order_mode_automatic) return
    if (coupling_order_selection%resolved_as2.ne.coupling_order_unbounded) return
    ! LHE files predating coupling-order metadata contain the historical
    ! maximum-QCD sample.  Recover that order from their external content.
    max_as2=-1
    do iprocess=1,unique_nproc
       nsing=0
       do ileg=1,next
          if (phys_model%is_singlet(unique_processes(ileg,iprocess))) nsing=nsing+1
       enddo
       max_as2=max(max_as2,2*max(0,next-2-nsing))
    enddo
    ! Some legacy LHE files predate the unique-process table altogether.
    ! Their physical event records are already loaded at this point, so use
    ! those external labels to recover the same historical leading-QCD order.
    ! Taking the maximum keeps heterogeneous legacy samples deterministic.
    if (max_as2.lt.0) then
       do ievent=1,nevents
          nsing=0
          do ileg=1,events(ievent)%NUP
             if (phys_model%is_singlet(events(ievent)%IDUP(ileg))) nsing=nsing+1
          enddo
          max_as2=max(max_as2,2*max(0,events(ievent)%NUP-2-nsing))
       enddo
    endif
    call resolve_automatic_coupling_order(max_as2,valid_selection)
    if (.not.valid_selection) then
       write (*,*) 'Could not infer the automatic coupling order from legacy LHE metadata'
       stop 1
    endif
  end subroutine resolve_legacy_automatic_selection

  real(kind=8) function reweight_ratio(numerator,denominator,label)
    use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
    implicit none
    real(kind=8),intent(in) :: numerator,denominator
    character(len=*),intent(in) :: label
    real(kind=8) :: reference

    if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator)) then
       write (*,*) 'Cannot form ',trim(label),' reweighting ratio from non-finite values',&
            denominator,numerator
       stop 1
    endif
    reference=max(1d-300,abs(numerator),abs(denominator))
    if (abs(denominator).le.1d-14*reference) then
       if (abs(numerator).le.1d-14*reference) then
          reweight_ratio=0d0
          return
       endif
       write (*,*) 'Cannot form ',trim(label),' reweighting ratio: denominator is zero',&
            denominator,numerator
       stop 1
    endif
    reweight_ratio=numerator/denominator
    if (.not.ieee_is_finite(reweight_ratio)) then
       write (*,*) 'Non-finite ',trim(label),' reweighting ratio',reweight_ratio
       stop 1
    endif
  end function reweight_ratio

  subroutine prepare_output_selection()
    implicit none
    integer :: ievent
    real(kind=8) :: ratio
    real(kind=8),external :: ran2
    integer(kind=8) :: iseed
    common /to_seed/iseed

    nevents_output=nevents
    events(:)%keep=.true.
    if (.not.unwgt) return
    nevents_output=0
    do ievent=1,nevents
       ratio=reweight_ratio(events(ievent)%matrix2(3),events(ievent)%matrix2(1),&
            'full-colour/leading-colour')
       events(ievent)%keep=abs(events(ievent)%XWGTUP*ratio).gt.&
            max_abs_event_weight*ran2()
       if (events(ievent)%keep) nevents_output=nevents_output+1
    enddo
    if (.not.has_input_seed) then
       ! RANMAR multiplies the user seed by 31300 on first use.
       input_seed=iseed/31300_8
       has_input_seed=.true.
    endif
  end subroutine prepare_output_selection

  real(kind=8) function selected_colour_matrix2(amplitude,iacc,alpha_s,alpha_ew)
    implicit none
    real(kind=8),parameter :: pi_order=3.1415926535897932384626433832795d0
    type(amplitude_QCD),intent(in) :: amplitude
    integer,intent(in) :: iacc
    real(kind=8),intent(in) :: alpha_s,alpha_ew
    complex(kind=8),dimension(:,:),allocatable :: scaled
    complex(kind=8),dimension(:),allocatable :: combined
    real(kind=8) :: diagonal_i,diagonal_j,gs,gew
    integer :: isector,jsector,as2,aew2

    allocate(scaled(amplitude%n_amps,amplitude%n_sectors))
    allocate(combined(amplitude%n_amps))
    scaled=(0d0,0d0)
    gs=sqrt(4d0*pi_order*alpha_s)
    gew=sqrt(8d0*pi_order*alpha_ew)
    do isector=1,amplitude%n_sectors
       where (amplitude%sector_present(:,isector))
          scaled(:,isector)=amplitude%amps_by_order(:,isector)*&
               gs**amplitude%sector_powers(1,isector)*&
               gew**amplitude%sector_powers(2,isector)
       endwhere
    enddo

    selected_colour_matrix2=0d0
    do isector=1,amplitude%n_sectors
       diagonal_i=colour_quadratic(amplitude,scaled(:,isector),iacc)
       do jsector=isector,amplitude%n_sectors
          as2=amplitude%sector_powers(1,isector)+amplitude%sector_powers(1,jsector)
          aew2=amplitude%sector_powers(2,isector)+amplitude%sector_powers(2,jsector)
          if (.not.coupling_order_selection%allows(as2,aew2)) cycle
          if (jsector.eq.isector) then
             selected_colour_matrix2=selected_colour_matrix2+diagonal_i
          else
             diagonal_j=colour_quadratic(amplitude,scaled(:,jsector),iacc)
             combined=scaled(:,isector)+scaled(:,jsector)
             selected_colour_matrix2=selected_colour_matrix2+&
                  colour_quadratic(amplitude,combined,iacc)-diagonal_i-diagonal_j
          endif
       enddo
    enddo
    deallocate(scaled,combined)
  end function selected_colour_matrix2

  real(kind=8) function colour_quadratic(amplitude,values,iacc)
    implicit none
    type(amplitude_QCD),intent(in) :: amplitude
    complex(kind=8),dimension(:),intent(in) :: values
    integer,intent(in) :: iacc
    complex(kind=8) :: row_value,colour_sum
    integer :: irow,icol,i,ic,ioffset

    colour_quadratic=0d0
    ioffset=amplitude%iproc_start(amplitude%nprocs)-1
    do irow=1,amplitude%nColOrd
       colour_sum=(0d0,0d0)
       do i=1,amplitude%n_col_vals(iacc)
          row_value=(0d0,0d0)
          do ic=amplitude%row_index(irow-1,i,iacc)+1,amplitude%row_index(irow,i,iacc)
             icol=amplitude%col_index(amplitude%i_col_i(i,iacc)+ic)
             row_value=row_value+values(ioffset+icol)
          enddo
          colour_sum=colour_sum+row_value*amplitude%diff_col_vals(i,iacc)
       enddo
       colour_quadratic=colour_quadratic+dble(colour_sum*conjg(values(ioffset+irow)))
    enddo
  end function colour_quadratic

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
    integer :: iproc,mode,resolved_as2,as_min2,as_max2,aew_min2,aew_max2,ios,&
         content_start,content_end
    character :: dummy
    character(len=1024) :: string
    logical :: valid_selection
    read(11,*) dummy ! <LesHouchesEvents>-tag
    read(11,*) dummy ! <header>-tag
    read(11,*) next,unique_nproc
    allocate(unique_map(unique_nproc))
    allocate(unique_map_value(unique_nproc))
    allocate(unique_processes(next,unique_nproc))
    do iproc=1,unique_nproc
       read(11,*) unique_map(iproc),unique_map_value(iproc),unique_processes(1:next,iproc)
    enddo
    read(11,'(a)',iostat=ios) string
    if (ios.eq.0 .and. index(adjustl(string),'<coupling_orders>').eq.1) then
       content_start=index(string,'<coupling_orders>')+len('<coupling_orders>')
       content_end=index(string,'</coupling_orders>')
       if (content_end.le.content_start) then
          write (*,*) 'Malformed coupling-order metadata in LHE header'
          stop 1
       endif
       read(string(content_start:content_end-1),*,iostat=ios) mode,resolved_as2,&
            as_min2,as_max2,aew_min2,aew_max2
       if (ios.ne.0) then
          write (*,*) 'Malformed coupling-order metadata in LHE header'
          stop 1
       endif
       call set_coupling_order_selection(mode,as_min2,as_max2,aew_min2,aew_max2,&
            valid_selection)
       if (.not.valid_selection) then
          write (*,*) 'Invalid coupling-order metadata in LHE header'
          stop 1
       endif
       if (mode.eq.coupling_order_mode_automatic) then
          call resolve_automatic_coupling_order(resolved_as2,valid_selection)
          if (.not.valid_selection) then
             write (*,*) 'Invalid resolved automatic coupling order in LHE header'
             stop 1
          endif
       elseif (resolved_as2.ne.coupling_order_unbounded) then
          write (*,*) 'Explicit coupling-order metadata must not contain an automatic order'
          stop 1
       endif
    else
       if (ios.eq.0) backspace(11)
       call reset_coupling_order_selection()
    endif
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
    integer :: i,iunit,ios
    character(len=1024) :: string
    character(len=1024) :: line
    logical :: found_colour
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
    found_colour=.false.
    event%overwgt=0d0
    do
       read(iunit,'(a)',iostat=ios) string
       if (ios.ne.0) then
          write (*,*) 'Unexpected end of LHE event while reading metadata'
          stop 1
       endif
       line=adjustl(string)
       if (index(line,'</event>').eq.1) exit
       if (index(line,'#color_expansion').eq.1) cycle
       if (index(line,'#color').eq.1) then
          read(line(len('#color')+1:),*,iostat=ios) event%col_order(1:event%NUP)
          if (ios.ne.0) then
             write (*,*) 'Malformed #color metadata in LHE event'
             stop 1
          endif
          found_colour=.true.
       elseif (index(line,'#overwgt').eq.1) then
          read(line(len('#overwgt')+1:),*,iostat=ios) event%overwgt(1:3)
          if (ios.ne.0) then
             write (*,*) 'Malformed #overwgt metadata in LHE event'
             stop 1
          endif
       endif
    enddo
    if (.not.found_colour) then
       write (*,*) 'LHE event has no #color metadata and cannot be colour-reweighted'
       stop 1
    endif

!!$    read (iunit,*,err=99,end=99) next,evt_wgt!,wgt,amp2,weight
!!$    read (iunit,*,err=99,end=99) helicity(1:next)
!!$    read (iunit,*,err=99,end=99) col_order(1:next)
!!$    do i=1,next
!!$       read (iunit,*,err=99,end=99) iPDG(i),momenta(1:3,i),momenta(0,i)
!!$    enddo
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
    if (.not.event%keep) return
    rwgt_NLC=reweight_ratio(event%matrix2(2),event%matrix2(1),&
         'next-to-leading-colour/leading-colour')
    rwgt_full=reweight_ratio(event%matrix2(3),event%matrix2(1),&
         'full-colour/leading-colour')
    if (unwgt) then
       XWGTUP_new=sign(xsec_abs,event%XWGTUP*rwgt_full)
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
    integer :: iproc,output_idwtup
    if (unwgt) then
       output_idwtup=-3
    else
       output_idwtup=-4
    endif
    write(ounit,'(a)') '<LesHouchesEvents version="3.0">'
    write(ounit,'(a)') '<header>'
    write(ounit,*) next,unique_nproc
    do iproc=1,unique_nproc
       write(ounit,*) unique_map(iproc),unique_map_value(iproc),unique_processes(1:next,iproc)
    enddo
    write(ounit,'(a,6(1x,i0),1x,a)') '<coupling_orders>',&
         coupling_order_selection%mode,coupling_order_selection%resolved_as2,&
         coupling_order_selection%as_min2,coupling_order_selection%as_max2,&
         coupling_order_selection%aew_min2,coupling_order_selection%aew_max2,&
         '</coupling_orders>'
    write(ounit,'(a,1x,i12,1x,a)') '<nevents>',nevents_output,'</nevents>'
    if (has_input_seed) write(ounit,'(a,1x,i12,1x,a)') '<seed>   ',input_seed,'</seed>'
    write(ounit,'(a)') '</header>'
    write(ounit,'(a)') '<init>'
    write(ounit,501)IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1)&
         &,PDFGUP(2),PDFSUP(1),PDFSUP(2),output_idwtup,NPRUP
    if (unwgt) then
       write(ounit,502)xsec,xsec_error,xsec_abs,LPRUP
    else
       write(ounit,502)xsec,xsec_error,XMAXUP*max_wgt,LPRUP
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
    nevents=0
    has_input_seed=.false.
    input_seed=0_8
    iseed=0_8
    do
       read(iunit,'(a)') string
       if (index(string,"<nevents>").ne.0) then
          read(string(10:),*) nevents
       elseif (index(string,"<seed>").ne.0) then
          read(string(10:),*) iseed
          input_seed=iseed
          has_input_seed=.true.
       endif
       if (index(string,"<init>").ne.0) then
          read(iunit,*) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1) &
               &,PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
          read(iunit,*) XSECUP,XERRUP,XMAXUP,LPRUP
          read(iunit,'(a)') generator_string
       endif
       if (index(string,"</init>").ne.0) exit
    enddo
    if (nevents.lt.1) then
       write (*,*) 'The input LHE file contains no events to reweight'
       stop 1
    endif
    allocate(events(nevents))
  end subroutine read_init_and_allocate_events
  
end program amplicol_reweight
