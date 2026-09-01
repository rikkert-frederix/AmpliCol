! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 amplitude_real.f03 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module rw_events
  implicit none
  type :: lhef_event
!!$     real(kind=8) :: wgt,evt_wgt,weight,rwgt_NLC,rwgt_full
!!$     integer,dimension(:),allocatable :: helicity,col_order,iPDG
!!$     real(kind=8),dimension(:,:),allocatable :: momenta
     integer :: NUP=0,IDPRUP=0
     integer,dimension(:),allocatable :: IDUP,ISTUP
     integer,dimension(:,:),allocatable :: MOTHUP,ICOLUP
     real(kind=8) :: XWGTUP=0d0,SCALUP=0d0,AQEDUP=0d0,AQCDUP=0d0
     real(kind=8),dimension(:,:),allocatable :: PUP
     real(kind=8),dimension(:),allocatable :: VTIMUP,SPINUP
     real(kind=8),dimension(3) :: overwgt=0d0
     integer,dimension(:),allocatable :: col_order
     real(kind=8),dimension(3) :: matrix2=0d0
  end type lhef_event
  type(lhef_event),dimension(:),allocatable :: events
  integer :: nevents=0
  integer :: IDBMUP(2)=0,PDFGUP(2)=0,PDFSUP(2)=0,IDWTUP=0,NPRUP=0,LPRUP=0
  real(kind=8) :: EBMUP(2)=0d0,XSECUP=0d0,XERRUP=0d0,XMAXUP=0d0
  character(len=1024) :: generator_string=''
  logical :: unwgt=.false.,keep_comments=.true.
end module rw_events
module timings
  implicit none
  real(kind=4) :: tBefore,tAfter,tTot_A=0.,tTot_B=0.,t_amp=0.,t_amp_init=0.,&
       t_mat_LC=0.,t_mat_NLC=0.,t_mat_full=0.,t_all=0.,t_ran=0.
end module timings
module overall
  real(kind=8) :: xsec,xsec_abs,max_wgt,max_abs_event_contribution
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
  use random_number_interface, only: ran2
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: iostat_eor,iostat_end
  implicit none
  logical,parameter :: use_only_canonical_form=.true.
  integer,parameter :: max_proc=1280,max_external_particles=max_amplitude_external_particles
  integer(kind=8),parameter :: max_reweight_workspace_bytes=2_8*1024_8*1024_8*1024_8
  real(kind=8),parameter :: reweight_value_limit=0.25d0*sqrt(huge(1d0))
  type(amplitude_QCD),dimension(max_proc) :: amps
  type(physics_model) :: phys_model
  logical :: read_proc_from_file=.false.
  integer :: i,col_acc,icol,irow,ic,iacc,next,nprocs,iproc,ioff,unique_nproc
  integer :: evaluation_status,ios
  integer,dimension(:),allocatable :: hel,unique_map
  integer,dimension(:,:),allocatable :: spin,o,part,processes,unique_processes
  real(kind=8) :: amp2,amp_col,amplitude_scale,rwgt_full,event_contribution,updated_sum
  real(kind=8) :: matrix_term,roundoff_limit
  real(kind=8),dimension(3) :: contraction_abs
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(:),allocatable :: unique_map_value
  complex(kind=8) :: amp2_c,amp_col_c
  character(len=256) :: event_filename,input_filename,io_message
  
  call get_run_arguments()
  io_message=''
  open(unit=99,status='scratch',action='write',form='formatted',iostat=ios,&
       iomsg=io_message)
  if (ios.ne.0) then
     write (*,*) 'Could not open the reweighter diagnostic stream: ',trim(io_message)
     stop 1
  endif
  call read_run_parameters(input_filename)

  call cpu_time(tTot_B)

  call phys_model%init_part()
  io_message=''
  open(unit=11,file=trim(event_filename),status='old',action='read',&
       form='formatted',iostat=ios,iomsg=io_message)
  if (ios.ne.0) then
     write (*,*) 'Could not open input event file: ',trim(event_filename),&
          trim(io_message)
     stop 1
  endif
  call read_unique_in_file()
  call allocate_process_info()
  call read_init_and_allocate_events(11)
  do nevt=1,nevents
     call read_event(11,events(nevt))
  enddo
  call validate_event_file_end(11)
  call apply_final_state_widths_from_events()
  call phys_model%init_vert()

  nprocs=0
  call setup_spin()
  col_acc=20

  xsec=0d0
  xsec_abs=0d0
  max_wgt=0d0
  max_abs_event_contribution=0d0
  do nevt=1,nevents
     call map_to_canonical_form(events(nevt))
     do iproc=1,nprocs
        if (all(part(1:next,1).eq.processes(1:next,iproc)))exit
     enddo
     if (iproc.eq.nprocs+1) then
        if (nprocs.ge.max_proc) then
           write (*,*) 'Too many distinct subprocesses in event file; limit is',max_proc
           stop 1
        endif
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
    
     call amps(iproc)%evaluate(next,p,hel,read_proc_from_file,phys_model,evaluation_status)
     if (evaluation_status.ne.0) then
        write (*,*) 'Amplitude evaluation failed while reweighting event',nevt,evaluation_status
        stop 1
     endif

     call cpu_time(tAfter)
     t_amp=t_amp+tAfter-tBefore

     if (amps(iproc)%nprocs.lt.1 .or. amps(iproc)%nColOrd.lt.1) then
        write (*,*) 'Invalid amplitude dimensions while reweighting event',nevt,&
             amps(iproc)%nprocs,amps(iproc)%nColOrd
        stop 1
     endif
     if (.not.allocated(amps(iproc)%iproc_start)) then
        write (*,*) 'Missing subprocess offsets while reweighting event',nevt
        stop 1
     endif
     if (size(amps(iproc)%iproc_start).lt.amps(iproc)%nprocs) then
        write (*,*) 'Truncated subprocess offsets while reweighting event',nevt
        stop 1
     endif
     ioff=amps(iproc)%iproc_start(amps(iproc)%nprocs)-1
     if (ioff.lt.0) then
        write (*,*) 'Invalid amplitude offset while reweighting event',nevt,ioff
        stop 1
     endif
     ! Only ratios of the three colour contractions are used.  Normalise the
     ! freshly evaluated amplitude once before the quadratic contraction so
     ! extremely small or large, but otherwise valid, amplitudes cannot
     ! underflow/overflow while leaving those ratios unchanged.
     if (use_real_gluons .and. amps(iproc)%n_qqbar(1).eq.0) then
        if (.not.allocated(amps(iproc)%amps_r)) then
           write (*,*) 'Missing real-amplitude workspace while reweighting event',nevt
           stop 1
        endif
        if (size(amps(iproc)%amps_r).lt.ioff+amps(iproc)%nColOrd) then
           write (*,*) 'Truncated real-amplitude workspace while reweighting event',&
                nevt,size(amps(iproc)%amps_r),ioff,amps(iproc)%nColOrd
           stop 1
        endif
        if (.not.all(ieee_is_finite(amps(iproc)%amps_r))) then
           write (*,*) 'Non-finite real amplitude while reweighting event',nevt
           stop 1
        endif
        amplitude_scale=maxval(abs(amps(iproc)%amps_r))
        if (amplitude_scale.le.0d0) then
           write (*,*) 'Zero real amplitude while reweighting event',nevt
           stop 1
        endif
        amps(iproc)%amps_r=amps(iproc)%amps_r/amplitude_scale
     else
        if (.not.allocated(amps(iproc)%amps)) then
           write (*,*) 'Missing complex-amplitude workspace while reweighting event',nevt
           stop 1
        endif
        if (size(amps(iproc)%amps).lt.ioff+amps(iproc)%nColOrd) then
           write (*,*) 'Truncated complex-amplitude workspace while reweighting event',&
                nevt,size(amps(iproc)%amps),ioff,amps(iproc)%nColOrd
           stop 1
        endif
        if (.not.all(ieee_is_finite(real(amps(iproc)%amps,kind=8))) .or. &
             .not.all(ieee_is_finite(aimag(amps(iproc)%amps)))) then
           write (*,*) 'Non-finite complex amplitude while reweighting event',nevt
           stop 1
        endif
        amplitude_scale=max(maxval(abs(real(amps(iproc)%amps,kind=8))),&
             maxval(abs(aimag(amps(iproc)%amps))))
        if (amplitude_scale.le.0d0) then
           write (*,*) 'Zero complex amplitude while reweighting event',nevt
           stop 1
        endif
        amps(iproc)%amps=amps(iproc)%amps/amplitude_scale
     endif
     contraction_abs=0d0
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
              matrix_term=amp_col*amps(iproc)%amps_r(ioff+irow)
              events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)+matrix_term
              contraction_abs(iacc)=contraction_abs(iacc)+abs(matrix_term)
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
              matrix_term=dble(amp_col_c*conjg(amps(iproc)%amps(ioff+irow)))
              events(nevt)%matrix2(iacc)=events(nevt)%matrix2(iacc)+matrix_term
              contraction_abs(iacc)=contraction_abs(iacc)+abs(matrix_term)
           enddo
        endif
        call cpu_time(tAfter)
        if (iacc.eq.1) t_mat_LC=t_mat_LC+tAfter-tBefore
        if (iacc.eq.2) t_mat_NLC=t_mat_NLC+tAfter-tBefore
        if (iacc.eq.3) t_mat_full=t_mat_full+tAfter-tBefore
     enddo
     if (.not.all(ieee_is_finite(events(nevt)%matrix2)) .or. &
          .not.all(ieee_is_finite(contraction_abs))) then
        write (*,*) 'Non-finite colour-summed matrix element for event',nevt
        stop 1
     endif
     do iacc=1,3,2
        roundoff_limit=512d0*epsilon(1d0)*contraction_abs(iacc)
        if (events(nevt)%matrix2(iacc).lt.0d0 .and. &
             events(nevt)%matrix2(iacc).ge.-roundoff_limit) &
             events(nevt)%matrix2(iacc)=0d0
     enddo
     if (events(nevt)%matrix2(1).le.0d0 .or. events(nevt)%matrix2(3).lt.0d0) then
        write (*,*) 'Invalid leading/full-colour matrix elements for event',nevt,&
             events(nevt)%matrix2(1),events(nevt)%matrix2(3)
        stop 1
     endif
     if (.not.safe_real_ratio(events(nevt)%matrix2(3),events(nevt)%matrix2(1),rwgt_full)) then
        write (*,*) 'Unsafe full-colour reweighting ratio for event',nevt
        stop 1
     endif
     if (.not.safe_real_product(rwgt_full,events(nevt)%XWGTUP,event_contribution)) then
        write (*,*) 'Unsafe reweighted event contribution for event',nevt
        stop 1
     endif
     if (.not.safe_real_sum(xsec,event_contribution,updated_sum)) then
        write (*,*) 'Full-colour cross-section accumulator overflow'
        stop 1
     endif
     xsec=updated_sum
     if (.not.safe_real_sum(xsec_abs,abs(event_contribution),updated_sum)) then
        write (*,*) 'Absolute cross-section accumulator overflow'
        stop 1
     endif
     xsec_abs=updated_sum
     max_wgt=max(max_wgt,rwgt_full)
     max_abs_event_contribution=max(max_abs_event_contribution,&
          abs(event_contribution))
  enddo
  xsec=xsec/dble(nevents)
  xsec_abs=xsec_abs/dble(nevents)
  io_message=''
  close(11,iostat=ios,iomsg=io_message)
  call require_reweight_io(ios,'closing the input event file',io_message)

  ! write event file
  io_message=''
  open(unit=12,file=trim(event_filename)//'.rwgt',status='replace',action='write',&
       form='formatted',iostat=ios,iomsg=io_message)
  if (ios.ne.0) then
     write (*,*) 'Could not open reweighted output file: ',&
          trim(event_filename)//'.rwgt',trim(io_message)
     stop 1
  endif
  call write_init(12)
  do nevt=1,nevents
     call write_event(12,events(nevt))
  enddo
  write(12,'(a)',iostat=ios,iomsg=io_message) '</LesHouchesEvents>'
  call require_reweight_io(ios,'writing the Les Houches closing tag',io_message)
  
  io_message=''
  close(12,iostat=ios,iomsg=io_message)
  call require_reweight_io(ios,'closing the output event file',io_message)
  io_message=''
  close(99,iostat=ios,iomsg=io_message)
  call require_reweight_io(ios,'closing the diagnostic stream',io_message)
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

  logical function safe_real_ratio(numerator,denominator,value)
    implicit none
    real(kind=8),intent(in) :: numerator,denominator
    real(kind=8),intent(out) :: value
    value=0d0
    safe_real_ratio=.false.
    if (.not.ieee_is_finite(numerator) .or. .not.ieee_is_finite(denominator)) return
    if (denominator.eq.0d0) return
    if (abs(denominator).lt.1d0) then
       if (abs(numerator).gt.huge(1d0)*abs(denominator)) return
    endif
    value=numerator/denominator
    if (.not.ieee_is_finite(value)) then
       value=0d0
       return
    endif
    safe_real_ratio=.true.
  end function safe_real_ratio

  logical function safe_real_product(first,second,value)
    implicit none
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    real(kind=8) :: abs_first,abs_second
    value=0d0
    safe_real_product=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (first.eq.0d0 .or. second.eq.0d0) then
       safe_real_product=.true.
       return
    endif
    abs_first=abs(first)
    abs_second=abs(second)
    if (abs_second.gt.1d0) then
       if (abs_first.gt.huge(1d0)/abs_second) return
    elseif (abs_first.gt.1d0) then
       if (abs_second.gt.huge(1d0)/abs_first) return
    endif
    value=first*second
    if (.not.ieee_is_finite(value)) then
       value=0d0
       return
    endif
    safe_real_product=.true.
  end function safe_real_product

  logical function safe_real_sum(first,second,value)
    implicit none
    real(kind=8),intent(in) :: first,second
    real(kind=8),intent(out) :: value
    value=0d0
    safe_real_sum=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second)) return
    if (second.gt.0d0) then
       if (first.gt.huge(1d0)-second) return
    elseif (second.lt.0d0) then
       if (first.lt.-huge(1d0)-second) return
    endif
    value=first+second
    if (.not.ieee_is_finite(value)) then
       value=0d0
       return
    endif
    safe_real_sum=.true.
  end function safe_real_sum

  subroutine get_run_arguments()
    use arguments
    implicit none
    integer :: argc,i,arg_length,arg_status,equal_position
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
       do i=1,argc
          arg_length=-1
          arg_status=-1
          call get_command_argument(i,length=arg_length,status=arg_status)
          if (arg_status.ne.0) then
             write (*,*) 'Could not query command-line argument',i
             stop 1
          endif
          if (arg_length.lt.1 .or. arg_length.gt.len(argv)) then
             write (*,*) 'Invalid or overlong command-line argument',i,arg_length
             stop 1
          endif
          argv=''
          arg_status=-1
          call get_command_argument(i,argv,status=arg_status)
          if (arg_status.ne.0) then
             write (*,*) 'Could not read command-line argument',i
             stop 1
          endif
          argv = trim(argv)
          if (argv.eq.'--help' .or. argv.eq.'-h') then
             show_help = .true.
          elseif (argv.eq.'--unwgt') then
             unwgt=.true.
          elseif (argv.eq.'--remove_comments') then
             keep_comments=.false.
          elseif (index(argv, "--input=").eq.1 .or. index(argv, "--card=").eq.1) then
             equal_position=index(argv,'=')
             if (equal_position.eq.len_trim(argv)) then
                write (*,*) 'The input-card filename cannot be empty'
                stop 1
             endif
             if (len_trim(argv)-equal_position.gt.len(input_filename)) then
                write (*,*) 'The input-card filename is too long'
                stop 1
             endif
             input_filename=argv(equal_position+1:len_trim(argv))
          elseif (index(argv,"-").eq.1) then
             write (*,*) 'Unknown argument: ',trim(argv)
             stop 1
          else
             if (len_trim(event_filename).ne.0) then
                write (*,*) 'More than one event file was specified'
                stop 1
             endif
             if (len_trim(argv).gt.len(event_filename)) then
                write (*,*) 'The event filename is too long'
                stop 1
             endif
             event_filename=argv(1:len_trim(argv))
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
    integer :: ievent,i,allocation_status
    integer,dimension(:,:),allocatable :: event_process
    character(len=256) :: allocation_message
    if (ignore_final_state_width_fix) return
    if (nevents.eq.0) return
    if (int(events(1)%NUP,kind=8)*int(nevents,kind=8).gt.&
         max_reweight_workspace_bytes/4_8) then
       write (*,*) 'Final-state width workspace exceeds the supported limit'
       stop 1
    endif
    allocate(event_process(1:events(1)%NUP,1:nevents),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate final-state width workspace: ',&
            trim(allocation_message)
       stop 1
    endif
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
    integer :: allocation_status
    character(len=256) :: allocation_message
    if (.not. allocated(hel)) then
!!$       allocate(momenta(0:3,1:next))
!!$       allocate(helicity(1:next))
!!$       allocate(col_order(1:next))
!!$       allocate(iPDG(1:next))
       allocate(o(1:next,1),part(1:next,1),processes(1:next,1:max_proc),&
            hel(1:next),p(0:3,1:next),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate reweighting process workspace: ',&
               trim(allocation_message)
          stop 1
       endif
       o=0
       part=0
       processes=0
       hel=0
       p=0d0
    endif
  end subroutine allocate_process_info
  
  subroutine read_unique_in_file()
    implicit none
    integer :: iproc,ios,label,mapped_process,allocation_status
    integer(kind=8) :: workspace_bytes,bytes_per_process,iseed
    logical :: nevents_seen,seed_seen,header_closed
    character(len=1024) :: line,stripped
    character(len=64) :: trailing_token
    character(len=256) :: allocation_message
    common /to_seed/iseed

    call read_required_reweight_line(11,line,'Les Houches opening tag')
    stripped=trim(adjustl(line))
    if (stripped.ne.'<LesHouchesEvents version="3.0">') then
       write (*,*) 'Invalid Les Houches opening tag: ',trim(line)
       stop 1
    endif
    call read_required_reweight_line(11,line,'Les Houches header opening tag')
    if (trim(adjustl(line)).ne.'<header>') then
       write (*,*) 'Invalid Les Houches header opening tag: ',trim(line)
       stop 1
    endif
    call read_required_reweight_line(11,line,'AmpliCol process dimensions')
    next=huge(0)
    unique_nproc=huge(0)
    trailing_token=''
    read(line,*,iostat=ios) next,unique_nproc,trailing_token
    call validate_reweight_list_parse(ios,trailing_token,line,&
         'AmpliCol process dimensions')
    if (next.lt.4 .or. next.gt.max_external_particles .or. &
         unique_nproc.lt.1 .or. unique_nproc.gt.max_proc) then
       write (*,*) 'Invalid AmpliCol event-file process header',next,unique_nproc
       stop 1
    endif
    bytes_per_process=12_8+4_8*int(next,kind=8)
    if (int(unique_nproc,kind=8).gt.&
         max_reweight_workspace_bytes/bytes_per_process) then
       write (*,*) 'Unique-process header exceeds the reweighter workspace limit',&
            next,unique_nproc
       stop 1
    endif
    workspace_bytes=int(unique_nproc,kind=8)*bytes_per_process
    allocate(unique_map(1:unique_nproc),unique_map_value(1:unique_nproc),&
         unique_processes(1:next,1:unique_nproc),stat=allocation_status,&
         errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate unique-process metadata: ',&
            trim(allocation_message),workspace_bytes
       stop 1
    endif
    do iproc=1,unique_nproc
       call read_required_reweight_line(11,line,'unique subprocess entry')
       unique_map(iproc)=huge(0)
       unique_map_value(iproc)=ieee_value(0d0,ieee_quiet_nan)
       unique_processes(:,iproc)=huge(0)
       trailing_token=''
       read(line,*,iostat=ios) unique_map(iproc),unique_map_value(iproc),&
            unique_processes(1:next,iproc),trailing_token
       call validate_reweight_list_parse(ios,trailing_token,line,&
            'unique subprocess entry')
       if (.not.ieee_is_finite(unique_map_value(iproc))) then
          write (*,*) 'Non-finite unique-process mapping factor',iproc
          stop 1
       endif
       if (unique_map(iproc).lt.-1 .or. unique_map(iproc).gt.unique_nproc .or. &
            unique_map_value(iproc).lt.0d0 .or. &
            unique_map_value(iproc).gt.reweight_value_limit) then
          write (*,*) 'Invalid unique-process mapping entry',iproc,&
               unique_map(iproc),unique_map_value(iproc)
          stop 1
       endif
       do label=1,next
          if (.not.supported_reweight_particle(unique_processes(label,iproc))) then
             write (*,*) 'Unsupported particle in unique-process entry',&
                  iproc,label,unique_processes(label,iproc)
             stop 1
          endif
       enddo
       do mapped_process=1,iproc-1
          if (all(unique_processes(:,iproc).eq.&
               unique_processes(:,mapped_process))) then
             write (*,*) 'Duplicate subprocess in unique-process header',&
                  mapped_process,iproc
             stop 1
          endif
       enddo
    enddo
    do iproc=1,unique_nproc
       mapped_process=unique_map(iproc)
       if (mapped_process.eq.-1) then
          if (unique_map_value(iproc).ne.1d0) then
             write (*,*) 'Canonical unique process has an invalid factor',iproc
             stop 1
          endif
       elseif (mapped_process.eq.0) then
          if (unique_map_value(iproc).ne.0d0) then
             write (*,*) 'Zero unique process has a nonzero factor',iproc
             stop 1
          endif
       else
          if (mapped_process.ge.iproc) then
             write (*,*) 'Unique-process mapping is forward or cyclic',iproc,mapped_process
             stop 1
          endif
          if (unique_map(mapped_process).ne.-1 .or. &
               unique_map_value(iproc).le.0d0) then
             write (*,*) 'Unique-process mapping target is not canonical',&
                  iproc,mapped_process
             stop 1
          endif
          do label=1,next
             if (phys_model%get_mass(unique_processes(label,iproc)).ne.&
                  phys_model%get_mass(unique_processes(label,mapped_process))) then
                write (*,*) 'Unique-process mapping changes the mass layout',&
                     iproc,mapped_process,label
                stop 1
             endif
          enddo
       endif
    enddo

    nevents=-1
    iseed=-1_8
    nevents_seen=.false.
    seed_seen=.false.
    header_closed=.false.
    do
       call read_required_reweight_line(11,line,'AmpliCol header metadata')
       stripped=trim(adjustl(line))
       if (stripped.eq.'</header>') then
          header_closed=.true.
          exit
       elseif (index(stripped,'<nevents>').eq.1) then
          if (nevents_seen) then
             write (*,*) 'Duplicate event-count record in AmpliCol header'
             stop 1
          endif
          call parse_tagged_default_integer(stripped,'nevents',nevents)
          nevents_seen=.true.
       elseif (index(stripped,'<seed>').eq.1) then
          if (seed_seen) then
             write (*,*) 'Duplicate random-seed record in AmpliCol header'
             stop 1
          endif
          call parse_tagged_int64(stripped,'seed',iseed)
          seed_seen=.true.
       elseif (len_trim(stripped).eq.0 .or. stripped(1:1).eq.'#') then
          cycle
       else
          write (*,*) 'Unexpected record in AmpliCol header: ',trim(line)
          stop 1
       endif
    enddo
    if (.not.header_closed .or. .not.nevents_seen .or. .not.seed_seen .or. &
         nevents.lt.1 .or. iseed.lt.0_8 .or. iseed.gt.904866561_8) then
       write (*,*) 'Incomplete or invalid AmpliCol header metadata',&
            nevents,iseed
       stop 1
    endif
  end subroutine read_unique_in_file
  
  subroutine setup_spin()
    implicit none
    integer :: allocation_status
    character(len=256) :: allocation_message
    if (.not.allocated(spin)) then
       allocate(spin(0:3,1:next),stat=allocation_status,&
            errmsg=allocation_message)
       if (allocation_status.ne.0) then
          write (*,*) 'Could not allocate reweighting spin workspace: ',&
               trim(allocation_message)
          stop 1
       endif
    endif
    spin=0
    do i=1,next
       spin(0,i)=1  ! one arbitrary spin state (use '-9')
       spin(1,i)=-9
    enddo
  end subroutine setup_spin

  subroutine read_event(iunit,event)
    implicit none
    type(lhef_event),intent(out) :: event
    integer,intent(in) :: iunit
    integer :: i,ios,allocation_status,spin_label
    character(len=1024) :: string,stripped
    character(len=64) :: trailing_token
    character(len=256) :: allocation_message
    real(kind=8) :: model_mass,comparison_scale

    call read_required_reweight_line(iunit,string,'event opening tag')
    if (trim(adjustl(string)).ne.'<event>') then
       write (*,*) 'Invalid event opening tag: ',trim(string)
       stop 1
    endif
    call read_required_reweight_line(iunit,string,'event header')
    event%NUP=huge(0)
    event%IDPRUP=huge(0)
    event%XWGTUP=ieee_value(0d0,ieee_quiet_nan)
    event%SCALUP=ieee_value(0d0,ieee_quiet_nan)
    event%AQEDUP=ieee_value(0d0,ieee_quiet_nan)
    event%AQCDUP=ieee_value(0d0,ieee_quiet_nan)
    trailing_token=''
    read(string,*,iostat=ios) event%NUP,event%IDPRUP,event%XWGTUP,&
         event%SCALUP,event%AQEDUP,event%AQCDUP,trailing_token
    call validate_reweight_list_parse(ios,trailing_token,string,'event header')
    if (event%NUP.ne.next .or. event%IDPRUP.ne.1) then
       write (*,*) 'Invalid event header or particle multiplicity',event%NUP,next
       stop 1
    endif
    if (.not.ieee_is_finite(event%XWGTUP) .or. &
         .not.ieee_is_finite(event%SCALUP) .or. &
         .not.ieee_is_finite(event%AQEDUP) .or. &
         .not.ieee_is_finite(event%AQCDUP)) then
       write (*,*) 'Non-finite event-header value'
       stop 1
    endif
    if (abs(event%XWGTUP).gt.reweight_value_limit .or. &
         event%SCALUP.lt.0d0 .or. event%SCALUP.gt.reweight_value_limit .or. &
         event%AQEDUP.lt.0d0 .or. event%AQEDUP.gt.reweight_value_limit .or. &
         event%AQCDUP.lt.0d0 .or. event%AQCDUP.gt.reweight_value_limit) then
       write (*,*) 'Unsafe event-header value'
       stop 1
    endif
    allocate(event%IDUP(1:event%NUP),event%ISTUP(1:event%NUP),&
         event%MOTHUP(1:2,1:event%NUP),event%ICOLUP(1:2,1:event%NUP),&
         event%PUP(1:5,1:event%NUP),event%VTIMUP(1:event%NUP),&
         event%SPINUP(1:event%NUP),event%col_order(1:event%NUP),&
         stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate event payload: ',trim(allocation_message)
       stop 1
    endif
    do i=1,event%NUP
       event%IDUP(i)=huge(0)
       event%ISTUP(i)=huge(0)
       event%MOTHUP(:,i)=-1
       event%ICOLUP(:,i)=-1
       event%PUP(:,i)=ieee_value(0d0,ieee_quiet_nan)
       event%VTIMUP(i)=ieee_value(0d0,ieee_quiet_nan)
       event%SPINUP(i)=ieee_value(0d0,ieee_quiet_nan)
       call read_required_reweight_line(iunit,string,'event particle row')
       trailing_token=''
       read(string,*,iostat=ios) event%IDUP(I),event%ISTUP(I),&
            event%MOTHUP(1,I),event%MOTHUP(2,I), &
            event%ICOLUP(1,I),event%ICOLUP(2,I),&
            event%PUP(1,I),event%PUP(2,I),event%PUP(3,I),event%PUP(4,I),event%PUP(5,I),&
            event%VTIMUP(I),event%SPINUP(I),trailing_token
       call validate_reweight_list_parse(ios,trailing_token,string,'event particle row')
       if (.not.supported_reweight_particle(event%IDUP(i))) then
          write (*,*) 'Unsupported particle in event row',i,event%IDUP(i)
          stop 1
       endif
       if (i.le.2) then
          if (event%ISTUP(i).ne.-1 .or. any(event%MOTHUP(:,i).ne.0) .or. &
               .not.((abs(event%IDUP(i)).ge.1 .and. abs(event%IDUP(i)).le.6) .or. &
               event%IDUP(i).eq.21 .or. event%IDUP(i).eq.22)) then
             write (*,*) 'Invalid incoming particle row',i,event%IDUP(i),&
                  event%ISTUP(i),event%MOTHUP(:,i)
             stop 1
          endif
       else
          if (event%ISTUP(i).ne.1 .or. any(event%MOTHUP(:,i).ne.[1,2])) then
             write (*,*) 'Invalid final-state particle row',i,event%ISTUP(i),&
                  event%MOTHUP(:,i)
             stop 1
          endif
       endif
       if (any(event%ICOLUP(:,i).lt.0)) then
          write (*,*) 'Negative colour label in event row',i,event%ICOLUP(:,i)
          stop 1
       endif
       if (.not.all(ieee_is_finite(event%PUP(:,i))) .or. &
            .not.ieee_is_finite(event%VTIMUP(i)) .or. &
            .not.ieee_is_finite(event%SPINUP(i))) then
          write (*,*) 'Non-finite particle data in event row',i
          stop 1
       endif
       if (any(abs(event%PUP(:,i)).gt.reweight_value_limit) .or. &
            event%PUP(4,i).lt.0d0 .or. event%PUP(5,i).lt.0d0 .or. &
            event%VTIMUP(i).lt.0d0 .or. &
            event%VTIMUP(i).gt.reweight_value_limit .or. &
            abs(event%SPINUP(i)).gt.1d0) then
          write (*,*) 'Unsafe particle data in event row',i
          stop 1
       endif
       spin_label=nint(event%SPINUP(i))
       if (event%SPINUP(i).ne.dble(spin_label)) then
          write (*,*) 'Non-integral helicity in event row',i,event%SPINUP(i)
          stop 1
       endif
       model_mass=phys_model%get_mass(event%IDUP(i))
       if (.not.ieee_is_finite(model_mass) .or. model_mass.lt.0d0) then
          write (*,*) 'Invalid model mass for event particle',i,event%IDUP(i)
          stop 1
       endif
       comparison_scale=max(1d0,abs(model_mass),abs(event%PUP(5,i)))
       if (abs(event%PUP(5,i)-model_mass).gt.1d-10*comparison_scale) then
          write (*,*) 'Event particle mass disagrees with the physics model',&
               i,event%PUP(5,i),model_mass
          stop 1
       endif
    enddo
    call validate_reweight_colour_flow(event)
    call validate_reweight_event_kinematics(event)

    call read_required_reweight_line(iunit,string,&
         'AmpliCol colour-order event record')
    stripped=adjustl(string)
    if (index(stripped,'#color').ne.1) then
       write (*,*) 'Invalid AmpliCol colour-order event record'
       stop 1
    endif
    event%col_order=huge(0)
    trailing_token=''
    read(stripped(7:),*,iostat=ios) event%col_order(1:event%NUP),trailing_token
    call validate_reweight_list_parse(ios,trailing_token,string,&
         'AmpliCol colour-order event record')
    if (any(event%col_order.lt.1) .or. any(event%col_order.gt.event%NUP)) then
       write (*,*) 'AmpliCol colour-order record contains an invalid label'
       stop 1
    endif
    do i=1,event%NUP
       if (count(event%col_order.eq.i).ne.1) then
          write (*,*) 'AmpliCol colour-order record is not a permutation'
          stop 1
       endif
    enddo
    call read_required_reweight_line(iunit,string,'AmpliCol overweight event record')
    stripped=adjustl(string)
    if (index(stripped,'#overwgt').ne.1) then
       write (*,*) 'Invalid AmpliCol overweight event record'
       stop 1
    endif
    event%overwgt=ieee_value(0d0,ieee_quiet_nan)
    trailing_token=''
    read(stripped(9:),*,iostat=ios) event%overwgt(1:3),trailing_token
    call validate_reweight_list_parse(ios,trailing_token,string,&
         'AmpliCol overweight event record')
    if (.not.all(ieee_is_finite(event%overwgt))) then
       write (*,*) 'Invalid AmpliCol overweight event record'
       stop 1
    endif
    if (any(event%overwgt.lt.0d0) .or. &
         any(event%overwgt.gt.reweight_value_limit)) then
       write (*,*) 'Invalid AmpliCol overweight event record'
       stop 1
    endif
    if (.not.reweight_values_agree(event%overwgt(1),abs(event%XWGTUP),1d-12)) then
       write (*,*) 'Event weight and AmpliCol overweight record disagree',&
            event%XWGTUP,event%overwgt(1)
       stop 1
    endif
    call read_required_reweight_line(iunit,string,'event closing tag')
    if (trim(adjustl(string)).ne.'</event>') then
       write (*,*) 'Invalid event closing tag: ',trim(string)
       stop 1
    endif
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
    if (any(array.lt.lbound(array_map,1)) .or. any(array.gt.ubound(array_map,1))) then
       write (*,*) 'Unsupported PDG code in canonical event mapping',array
       stop 1
    endif
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
    integer :: i,iproc
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
       do iproc=1,unique_nproc
          if (all(part(1:next,1).eq.unique_processes(1:next,iproc))) then
             if (unique_map(iproc).eq.0) then
                write (*,*) 'Event belongs to a subprocess declared identically zero',&
                     iproc
                stop 1
             endif
             ! The integration-time flavour map was inferred from
             ! proportional colour-summed matrix-element samples.  That is
             ! sufficient for the integration weight, but it does not prove
             ! proportionality of each colour-ordered amplitude.  Mapping
             ! here could therefore change the full/LC ratio.  Keep the exact
             ! canonical event flavour; use the header only to authenticate
             ! that it belongs to the generated process set.
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
    integer,intent(in) :: iunit
    integer :: i,io_status
    character(len=256) :: message
    real(kind=8) :: rwgt_NLC,rwgt_full,XWGTUP_new,random_value,threshold
    real(kind=8) :: weighted_nlc,weighted_full
    if (.not.safe_real_ratio(event%matrix2(2),event%matrix2(1),rwgt_NLC) .or. &
         .not.safe_real_ratio(event%matrix2(3),event%matrix2(1),rwgt_full)) then
       write (*,*) 'Cannot write event with an unsafe colour-expansion ratio'
       stop 1
    endif
    if (.not.safe_real_product(event%XWGTUP,rwgt_NLC,weighted_nlc) .or. &
         .not.safe_real_product(event%XWGTUP,rwgt_full,weighted_full)) then
       write (*,*) 'Unsafe colour-expansion event weight'
       stop 1
    endif
    if (unwgt) then
       random_value=ran2()
       if (.not.ieee_is_finite(random_value)) then
          write (*,*) 'Random-number generator returned a non-finite unweighting value'
          stop 1
       endif
       if (random_value.lt.0d0 .or. random_value.ge.1d0) then
          write (*,*) 'Random-number generator returned an invalid unweighting value',random_value
          stop 1
       endif
       if (.not.safe_real_product(max_abs_event_contribution,random_value,threshold)) then
          write (*,*) 'Unsafe unweighting threshold'
          stop 1
       endif
       if (abs(weighted_full).gt.threshold .and. weighted_full.ne.0d0) then
          XWGTUP_new=sign(xsec_abs,weighted_full)
       else
          return
       endif
    else
       XWGTUP_new=weighted_full
    endif
    message=''
    write (iunit,'(a)',iostat=io_status,iomsg=message) '<event>'
    call require_reweight_io(io_status,'writing an event opening tag',message)
    write(iunit,503,iostat=io_status,iomsg=message) event%NUP,event%IDPRUP,&
         XWGTUP_new,event%SCALUP,event%AQEDUP,event%AQCDUP
    call require_reweight_io(io_status,'writing an event header',message)
    do i=1,event%NUP
       write(iunit,504,iostat=io_status,iomsg=message) event%IDUP(I),&
            event%ISTUP(I),event%MOTHUP(1,I),event%MOTHUP(2,I), &
            event%ICOLUP(1,I),event%ICOLUP(2,I),&
            event%PUP(1,I),event%PUP(2,I),event%PUP(3,I),event%PUP(4,I),event%PUP(5,I),&
            event%VTIMUP(I),event%SPINUP(I)
       call require_reweight_io(io_status,'writing an event particle row',message)
    enddo
    if (keep_comments) then
       write (iunit,506,iostat=io_status,iomsg=message) &
            '#color',event%col_order(1:event%NUP)
       call require_reweight_io(io_status,'writing an event colour order',message)
       write (iunit,505,iostat=io_status,iomsg=message) '#overwgt',event%overwgt(1:3)
       call require_reweight_io(io_status,'writing an event overweight record',message)
       write (iunit,505,iostat=io_status,iomsg=message) &
            '#color_expansion',event%XWGTUP,weighted_nlc,weighted_full
       call require_reweight_io(io_status,'writing a colour-expansion record',message)
    endif
    write (iunit,'(a)',iostat=io_status,iomsg=message) '</event>'
    call require_reweight_io(io_status,'writing an event closing tag',message)
503 format(1x,i4,1x,i10,4(1x,es24.16e3))
504 format(1x,i10,1x,i4,4(1x,i10),5(1x,es26.17e3),2(1x,es16.8e3))
505 format(a,3(1x,es24.16e3))
506 format(a,100i3)
  end subroutine write_event

  subroutine write_init(ounit)
    use overall
    implicit none
    integer,intent(in) :: ounit
    integer :: io_status,output_weight_strategy
    real(kind=8) :: rwgt_LCtoFC,reweighted_error,reweighted_maximum
    character(len=256) :: message
    if (.not.safe_real_ratio(xsec,XSECUP,rwgt_LCtoFC)) then
       write (*,*) 'Cannot construct output init block from a zero/unsafe input cross section'
       stop 1
    endif
    if (.not.safe_real_product(XERRUP,abs(rwgt_LCtoFC),reweighted_error)) then
       write (*,*) 'Unsafe reweighted cross-section uncertainty'
       stop 1
    endif
    if (unwgt) then
       output_weight_strategy=-3
    else
       output_weight_strategy=-4
    endif
    message=''
    write(ounit,'(a)',iostat=io_status,iomsg=message) &
         '<LesHouchesEvents version="3.0">'
    call require_reweight_io(io_status,'writing the Les Houches opening tag',message)
    write(ounit,'(a)',iostat=io_status,iomsg=message) '<init>'
    call require_reweight_io(io_status,'writing the init opening tag',message)
    write(ounit,501,iostat=io_status,iomsg=message) IDBMUP(1),IDBMUP(2),&
         EBMUP(1),EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),&
         output_weight_strategy,NPRUP
    call require_reweight_io(io_status,'writing the beam-init row',message)
    if (unwgt) then
       write(ounit,502,iostat=io_status,iomsg=message) &
            xsec,reweighted_error,xsec_abs,LPRUP
    else
       if (.not.safe_real_product(XMAXUP,max_wgt,reweighted_maximum)) then
          write (*,*) 'Unsafe reweighted maximum event weight'
          stop 1
       endif
       write(ounit,502,iostat=io_status,iomsg=message) &
            xsec,reweighted_error,reweighted_maximum,LPRUP
    endif
    call require_reweight_io(io_status,'writing the subprocess-init row',message)
    write(ounit,'(a)',iostat=io_status,iomsg=message) trim(generator_string)
    call require_reweight_io(io_status,'writing the generator record',message)
    write(ounit,'(a)',iostat=io_status,iomsg=message) '</init>'
    call require_reweight_io(io_status,'writing the init closing tag',message)
501 format(2(1x,i10),2(1x,es24.16e3),2(1x,i4),2(1x,i10),1x,i4,1x,i6)
502 format(3(1x,es24.16e3),1x,i10)
  end subroutine write_init

  subroutine read_init_and_allocate_events(iunit)
    implicit none
    integer,intent(in) :: iunit
    character(len=1024) :: string,stripped
    character(len=64) :: trailing_token
    character(len=256) :: allocation_message
    type(lhef_event) :: event_template
    integer(kind=8) :: descriptor_bytes,per_event_bytes,event_workspace_bytes
    integer :: ios,allocation_status,generator_close,opening_close
    logical :: generator_seen,init_closed

    call read_required_reweight_line(iunit,string,'Les Houches init opening tag')
    if (trim(adjustl(string)).ne.'<init>') then
       write (*,*) 'Expected Les Houches init block but found: ',trim(string)
       stop 1
    endif
    IDBMUP=huge(0)
    EBMUP=ieee_value(0d0,ieee_quiet_nan)
    PDFGUP=huge(0)
    PDFSUP=huge(0)
    IDWTUP=huge(0)
    NPRUP=huge(0)
    call read_required_reweight_line(iunit,string,'Les Houches beam-init row')
    trailing_token=''
    read(string,*,iostat=ios) IDBMUP,EBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,&
         trailing_token
    call validate_reweight_list_parse(ios,trailing_token,string,&
         'Les Houches beam-init row')
    if (any(IDBMUP.eq.0) .or. any(IDBMUP.eq.huge(0)) .or. &
         any(PDFGUP.eq.huge(0)) .or. any(PDFSUP.eq.huge(0)) .or. &
         (IDWTUP.ne.-3 .and. IDWTUP.ne.-4) .or. NPRUP.ne.1) then
       write (*,*) 'Invalid discrete values in Les Houches beam-init row'
       stop 1
    endif
    if (.not.all(ieee_is_finite(EBMUP))) then
       write (*,*) 'Non-finite beam energy in Les Houches init block'
       stop 1
    endif
    if (any(EBMUP.le.0d0) .or. any(EBMUP.gt.reweight_value_limit)) then
       write (*,*) 'Invalid beam energy in Les Houches init block',EBMUP
       stop 1
    endif

    XSECUP=ieee_value(0d0,ieee_quiet_nan)
    XERRUP=ieee_value(0d0,ieee_quiet_nan)
    XMAXUP=ieee_value(0d0,ieee_quiet_nan)
    LPRUP=huge(0)
    call read_required_reweight_line(iunit,string,'Les Houches subprocess-init row')
    trailing_token=''
    read(string,*,iostat=ios) XSECUP,XERRUP,XMAXUP,LPRUP,trailing_token
    call validate_reweight_list_parse(ios,trailing_token,string,&
         'Les Houches subprocess-init row')
    if (LPRUP.ne.1 .or. .not.ieee_is_finite(XSECUP) .or. &
         .not.ieee_is_finite(XERRUP) .or. .not.ieee_is_finite(XMAXUP)) then
       write (*,*) 'Invalid values in Les Houches subprocess-init row'
       stop 1
    endif
    if (XSECUP.eq.0d0 .or. abs(XSECUP).gt.reweight_value_limit .or. &
         XERRUP.lt.0d0 .or. XERRUP.gt.reweight_value_limit .or. &
         XMAXUP.lt.0d0 .or. XMAXUP.gt.reweight_value_limit) then
       write (*,*) 'Unsafe values in Les Houches subprocess-init row',&
            XSECUP,XERRUP,XMAXUP
       stop 1
    endif

    generator_string=''
    generator_seen=.false.
    init_closed=.false.
    do
       call read_required_reweight_line(iunit,string,'Les Houches init metadata')
       stripped=trim(adjustl(string))
       if (stripped.eq.'</init>') then
          init_closed=.true.
          exit
       elseif (index(stripped,'<generator').eq.1) then
          if (generator_seen) then
             write (*,*) 'More than one generator record in Les Houches init block'
             stop 1
          endif
          generator_close=index(stripped,'</generator>')
          opening_close=index(stripped,'>')
          if (len_trim(stripped).lt.11 .or. &
               (stripped(11:11).ne.' ' .and. stripped(11:11).ne.'>') .or. &
               opening_close.lt.11 .or. generator_close.le.opening_close .or. &
               generator_close+len('</generator>')-1.ne.len_trim(stripped)) then
             write (*,*) 'Malformed generator record in Les Houches init block: ',&
                  trim(string)
             stop 1
          endif
          generator_string=trim(stripped)
          generator_seen=.true.
       elseif (len_trim(stripped).eq.0 .or. stripped(1:1).eq.'#') then
          cycle
       else
          write (*,*) 'Unexpected Les Houches init metadata: ',trim(string)
          stop 1
       endif
    enddo
    if (.not.generator_seen .or. .not.init_closed) then
       write (*,*) 'Incomplete Les Houches init block'
       stop 1
    endif
    if (nevents.lt.1) then
       write (*,*) 'Invalid or missing positive event count',nevents
       stop 1
    endif
    descriptor_bytes=int((storage_size(event_template)+7)/8,kind=8)
    per_event_bytes=descriptor_bytes+512_8+84_8*int(next,kind=8)
    if (per_event_bytes.le.0_8 .or. int(nevents,kind=8).gt.&
         max_reweight_workspace_bytes/per_event_bytes) then
       write (*,*) 'Event file exceeds the reweighter workspace limit',&
            nevents,per_event_bytes,max_reweight_workspace_bytes
       stop 1
    endif
    event_workspace_bytes=int(nevents,kind=8)*per_event_bytes
    allocate(events(1:nevents),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status.ne.0) then
       write (*,*) 'Could not allocate storage for events: ',nevents,&
            event_workspace_bytes,trim(allocation_message)
       stop 1
    endif
    events%NUP=0
    events%IDPRUP=0
    events%XWGTUP=0d0
    events%SCALUP=0d0
    events%AQEDUP=0d0
    events%AQCDUP=0d0
  end subroutine read_init_and_allocate_events

  subroutine read_required_reweight_line(iunit,line,record_name)
    implicit none
    integer,intent(in) :: iunit
    character(len=*),intent(out) :: line
    character(len=*),intent(in) :: record_name
    logical :: reached_end

    call read_optional_reweight_line(iunit,line,record_name,reached_end)
    if (reached_end) then
       write (*,*) 'Unexpected end of file while reading ',trim(record_name)
       stop 1
    endif
  end subroutine read_required_reweight_line

  subroutine read_optional_reweight_line(iunit,line,record_name,reached_end)
    implicit none
    integer,intent(in) :: iunit
    character(len=*),intent(out) :: line
    character(len=*),intent(in) :: record_name
    logical,intent(out) :: reached_end
    integer :: io_status,character_count,extra_count
    character(len=1) :: extra_character
    character(len=256) :: message

    line=''
    message=''
    reached_end=.false.
    read(iunit,'(a)',advance='no',iostat=io_status,iomsg=message,&
         size=character_count) line
    if (io_status.eq.iostat_eor) return
    if (io_status.eq.iostat_end) then
       reached_end=.true.
       return
    endif
    if (io_status.ne.0) then
       write (*,*) 'Could not read ',trim(record_name),': ',io_status,trim(message)
       stop 1
    endif
    extra_character=''
    message=''
    read(iunit,'(a)',advance='no',iostat=io_status,iomsg=message,&
         size=extra_count) extra_character
    if (io_status.eq.iostat_eor .and. extra_count.eq.0) return
    if (io_status.eq.0 .or. (io_status.eq.iostat_eor .and. extra_count.gt.0)) then
       write (*,*) trim(record_name),' exceeds the supported line length',len(line)
       stop 1
    endif
    write (*,*) 'Could not finish reading ',trim(record_name),': ',&
         io_status,trim(message)
    stop 1
  end subroutine read_optional_reweight_line

  subroutine validate_reweight_list_parse(io_status,trailing_token,line,record_name)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: trailing_token,line,record_name

    if (io_status.gt.0) then
       write (*,*) 'Could not parse ',trim(record_name),': ',trim(line)
       stop 1
    endif
    if (io_status.eq.0) then
       if (len_trim(trailing_token).eq.0 .or. trailing_token(1:1).ne.'#') then
          write (*,*) trim(record_name),' contains trailing data: ',trim(line)
          stop 1
       endif
    endif
  end subroutine validate_reweight_list_parse

  subroutine parse_tagged_default_integer(line,tag,value)
    implicit none
    character(len=*),intent(in) :: line,tag
    integer,intent(out) :: value
    character(len=128) :: opening_tag,closing_tag,payload,trailing_token
    integer :: io_status,payload_end

    opening_tag='<'//trim(tag)//'>'
    closing_tag='</'//trim(tag)//'>'
    payload=''
    value=huge(0)
    if (index(line,trim(opening_tag)).ne.1 .or. &
         len_trim(line).lt.len_trim(opening_tag)+len_trim(closing_tag)+1) then
       write (*,*) 'Malformed tagged integer record: ',trim(line)
       stop 1
    endif
    payload_end=len_trim(line)-len_trim(closing_tag)
    if (line(payload_end+1:len_trim(line)).ne.trim(closing_tag)) then
       write (*,*) 'Malformed tagged integer record: ',trim(line)
       stop 1
    endif
    payload=line(len_trim(opening_tag)+1:payload_end)
    trailing_token=''
    read(payload,*,iostat=io_status) value,trailing_token
    if (io_status.ge.0 .or. value.eq.huge(0)) then
       write (*,*) 'Invalid tagged integer value: ',trim(line)
       stop 1
    endif
  end subroutine parse_tagged_default_integer

  subroutine parse_tagged_int64(line,tag,value)
    implicit none
    character(len=*),intent(in) :: line,tag
    integer(kind=8),intent(out) :: value
    character(len=128) :: opening_tag,closing_tag,payload,trailing_token
    integer :: io_status,payload_end

    opening_tag='<'//trim(tag)//'>'
    closing_tag='</'//trim(tag)//'>'
    payload=''
    value=huge(0_8)
    if (index(line,trim(opening_tag)).ne.1 .or. &
         len_trim(line).lt.len_trim(opening_tag)+len_trim(closing_tag)+1) then
       write (*,*) 'Malformed tagged 64-bit integer record: ',trim(line)
       stop 1
    endif
    payload_end=len_trim(line)-len_trim(closing_tag)
    if (line(payload_end+1:len_trim(line)).ne.trim(closing_tag)) then
       write (*,*) 'Malformed tagged 64-bit integer record: ',trim(line)
       stop 1
    endif
    payload=line(len_trim(opening_tag)+1:payload_end)
    trailing_token=''
    read(payload,*,iostat=io_status) value,trailing_token
    if (io_status.ge.0 .or. value.eq.huge(0_8)) then
       write (*,*) 'Invalid tagged 64-bit integer value: ',trim(line)
       stop 1
    endif
  end subroutine parse_tagged_int64

  logical function supported_reweight_particle(ipdg)
    implicit none
    integer,intent(in) :: ipdg
    supported_reweight_particle=(ipdg.ge.1 .and. ipdg.le.6) .or. &
         (ipdg.le.-1 .and. ipdg.ge.-6) .or. &
         (ipdg.ge.11 .and. ipdg.le.16) .or. &
         (ipdg.le.-11 .and. ipdg.ge.-16) .or. &
         ipdg.eq.21 .or. ipdg.eq.22 .or. ipdg.eq.23 .or. &
         ipdg.eq.24 .or. ipdg.eq.-24 .or. ipdg.eq.25
  end function supported_reweight_particle

  logical function reweight_values_agree(first,second,tolerance)
    implicit none
    real(kind=8),intent(in) :: first,second,tolerance
    real(kind=8) :: comparison_scale
    reweight_values_agree=.false.
    if (.not.ieee_is_finite(first) .or. .not.ieee_is_finite(second) .or. &
         .not.ieee_is_finite(tolerance)) return
    if (tolerance.lt.0d0 .or. tolerance.gt.1d0) return
    comparison_scale=max(1d0,abs(first),abs(second))
    reweight_values_agree=abs(first-second).le.tolerance*comparison_scale
  end function reweight_values_agree

  subroutine validate_reweight_colour_flow(event)
    implicit none
    type(lhef_event),intent(in) :: event
    integer :: i,label

    if (.not.allocated(event%IDUP) .or. .not.allocated(event%ICOLUP)) then
       write (*,*) 'Cannot validate an event without particle colour data'
       stop 1
    endif
    if (size(event%IDUP).ne.event%NUP .or. size(event%ICOLUP,1).ne.2 .or. &
         size(event%ICOLUP,2).ne.event%NUP) then
       write (*,*) 'Cannot validate malformed event colour storage'
       stop 1
    endif
    do i=1,event%NUP
       if (event%IDUP(i).ge.1 .and. event%IDUP(i).le.6) then
          if (event%ICOLUP(1,i).le.0 .or. event%ICOLUP(2,i).ne.0) then
             write (*,*) 'Invalid quark colour assignment in event row',i,&
                  event%IDUP(i),event%ICOLUP(:,i)
             stop 1
          endif
       elseif (event%IDUP(i).le.-1 .and. event%IDUP(i).ge.-6) then
          if (event%ICOLUP(1,i).ne.0 .or. event%ICOLUP(2,i).le.0) then
             write (*,*) 'Invalid antiquark colour assignment in event row',i,&
                  event%IDUP(i),event%ICOLUP(:,i)
             stop 1
          endif
       elseif (event%IDUP(i).eq.21) then
          if (any(event%ICOLUP(:,i).le.0) .or. &
               event%ICOLUP(1,i).eq.event%ICOLUP(2,i)) then
             write (*,*) 'Invalid gluon colour assignment in event row',i,&
                  event%ICOLUP(:,i)
             stop 1
          endif
       elseif (any(event%ICOLUP(:,i).ne.0)) then
          write (*,*) 'Colourless particle has a nonzero colour label',i,&
               event%IDUP(i),event%ICOLUP(:,i)
          stop 1
       endif
    enddo
    do i=1,event%NUP
       do label=1,2
          if (event%ICOLUP(label,i).eq.0) cycle
          if (count(event%ICOLUP.eq.event%ICOLUP(label,i)).ne.2) then
             write (*,*) 'Unpaired or multiply used event colour label',&
                  event%ICOLUP(label,i)
             stop 1
          endif
       enddo
    enddo
  end subroutine validate_reweight_colour_flow

  subroutine validate_reweight_event_kinematics(event)
    implicit none
    type(lhef_event),intent(in) :: event
    integer :: i,component
    real(kind=8) :: momentum_square,energy_square,mass_square,shell_value,&
         comparison_scale
    real(kind=8),dimension(4) :: incoming_momentum,outgoing_momentum

    if (event%NUP.lt.4) then
       write (*,*) 'Cannot validate malformed event kinematics'
       stop 1
    endif
    if (.not.allocated(event%PUP)) then
       write (*,*) 'Cannot validate event kinematics without momenta'
       stop 1
    endif
    if (size(event%PUP,1).ne.5 .or. size(event%PUP,2).ne.event%NUP) then
       write (*,*) 'Cannot validate malformed event momentum storage'
       stop 1
    endif
    do i=1,event%NUP
       momentum_square=sum(event%PUP(1:3,i)*event%PUP(1:3,i))
       energy_square=event%PUP(4,i)*event%PUP(4,i)
       mass_square=event%PUP(5,i)*event%PUP(5,i)
       shell_value=energy_square-momentum_square
       comparison_scale=max(1d0,energy_square,momentum_square,mass_square)
       if (abs(shell_value-mass_square).gt.1d-8*comparison_scale) then
          write (*,*) 'Event momentum is not on its declared mass shell',&
               i,shell_value,mass_square
          stop 1
       endif
    enddo
    incoming_momentum(1:3)=event%PUP(1:3,1)+event%PUP(1:3,2)
    incoming_momentum(4)=event%PUP(4,1)+event%PUP(4,2)
    outgoing_momentum(1:3)=sum(event%PUP(1:3,3:event%NUP),dim=2)
    outgoing_momentum(4)=sum(event%PUP(4,3:event%NUP))
    do component=1,4
       comparison_scale=max(1d0,abs(incoming_momentum(component)),&
            abs(outgoing_momentum(component)),&
            sum(abs(event%PUP(component,1:event%NUP))))
       if (abs(incoming_momentum(component)-outgoing_momentum(component)).gt.&
            1d-8*comparison_scale) then
          write (*,*) 'Event violates four-momentum conservation',component,&
               incoming_momentum(component),outgoing_momentum(component)
          stop 1
       endif
    enddo
  end subroutine validate_reweight_event_kinematics

  subroutine validate_event_file_end(iunit)
    implicit none
    integer,intent(in) :: iunit
    logical :: closing_tag_seen,reached_end
    character(len=1024) :: line,stripped

    closing_tag_seen=.false.
    do
       call read_optional_reweight_line(iunit,line,'end of event file',reached_end)
       if (reached_end) exit
       stripped=trim(adjustl(line))
       if (.not.closing_tag_seen) then
          if (len_trim(stripped).eq.0) cycle
          if (stripped.ne.'</LesHouchesEvents>') then
             write (*,*) 'Unexpected record after the declared events: ',trim(line)
             stop 1
          endif
          closing_tag_seen=.true.
       elseif (len_trim(stripped).ne.0) then
          write (*,*) 'Trailing data after the Les Houches closing tag: ',trim(line)
          stop 1
       endif
    enddo
    if (.not.closing_tag_seen) then
       write (*,*) 'Missing Les Houches closing tag after the declared events'
       stop 1
    endif
  end subroutine validate_event_file_end

  subroutine require_reweight_io(io_status,operation,message)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: operation,message
    if (io_status.ne.0) then
       write (*,*) 'Reweighted-event I/O failure while ',trim(operation),': ',&
            io_status,trim(message)
       stop 1
    endif
  end subroutine require_reweight_io
  
end program amplicol_reweight
