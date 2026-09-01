module handling_events
  use common
  use handling_processes
  use amplitude_QCD_mod, only: max_event_external_particles => max_amplitude_external_particles
  use random_number_interface, only: ran2
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  integer :: iproc_picked=0,iproc_iden_picked=0
  integer,dimension(2) :: hel_picked=0
  real(kind=8) :: evt_sign=0d0
  logical :: process_selection_ready=.false.,helicity_selection_ready=.false.
  integer(kind=8) :: event_selection_counter=0_8,active_event_selection_epoch=0_8
  real(kind=8),parameter :: event_value_limit=0.25d0*sqrt(huge(1d0))
  integer,parameter :: max_event_process_records=1280
contains
  subroutine write_event(iunit,pgl,wgt)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,iunit,ios
    real(kind=8),intent(in) :: wgt
    integer :: NUP,IDPRUP
    integer,dimension(max_event_external_particles) :: IDUP,ISTUP
    integer,dimension(2,max_event_external_particles) :: MOTHUP,ICOLUP
    real(kind=8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP
    real(kind=8),dimension(5,max_event_external_particles) :: PUP
    real(kind=8),dimension(max_event_external_particles) :: VTIMUP,SPINUP
    logical :: unit_open
    character(len=20) :: unit_form,unit_action
    character(len=256) :: io_message
    unit_form=''
    unit_action=''
    io_message=''
    inquire(unit=iunit,opened=unit_open,form=unit_form,action=unit_action,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect the requested event-output unit:',&
            iunit,ios,trim(io_message)
       stop 1
    endif
    if (.not.unit_open .or. trim(unit_form).ne.'FORMATTED' .or. &
         trim(unit_action).eq.'READ') then
       write (*,*) 'Cannot write an event to the requested unit:',iunit,ios,&
            trim(unit_form),trim(unit_action),trim(io_message)
       stop 1
    endif
    if (.not.process_selection_ready .or. .not.helicity_selection_ready) then
       write (*,*) 'Cannot write an event without a fresh process/helicity selection'
       stop 1
    endif
    if (active_event_selection_epoch.le.0_8 .or. &
         pgl%event_selection_epoch.ne.active_event_selection_epoch) then
       write (*,*) 'Cannot write an event from a different or stale process group'
       stop 1
    endif
    if (pgl%next.lt.3 .or. pgl%next.gt.max_event_external_particles) then
       write (*,*) 'Cannot write an event with invalid phase-space state'
       stop 1
    endif
    if (pgl%nproc.lt.1 .or. pgl%nproc.gt.max_process_records) then
       write (*,*) 'Cannot write an event without subprocesses'
       stop 1
    endif
    if (.not.allocated(pgl%iden_iproc)) then
       write (*,*) 'Cannot write an event without identical-process counts'
       stop 1
    endif
    if (lbound(pgl%iden_iproc,1).ne.1 .or. &
         size(pgl%iden_iproc).ne.pgl%nproc) then
       write (*,*) 'Cannot write an event with incomplete identical-process counts'
       stop 1
    endif
    if (.not.allocated(pgl%iden_processes)) then
       write (*,*) 'Cannot write an event without process identities'
       stop 1
    endif
    if (lbound(pgl%iden_processes,1).ne.1 .or. &
         lbound(pgl%iden_processes,2).ne.1 .or. &
         lbound(pgl%iden_processes,3).ne.1 .or. &
         size(pgl%iden_processes,1).ne.pgl%next .or. &
         size(pgl%iden_processes,3).ne.pgl%nproc) then
       write (*,*) 'Cannot write an event with incomplete process identities'
       stop 1
    endif
    if (.not.allocated(pgl%color_orders)) then
       write (*,*) 'Cannot write an event without colour orders'
       stop 1
    endif
    if (lbound(pgl%color_orders,1).ne.1 .or. &
         lbound(pgl%color_orders,2).ne.1 .or. &
         size(pgl%color_orders,1).ne.pgl%next .or. &
         size(pgl%color_orders,2).ne.pgl%nproc) then
       write (*,*) 'Cannot write an event with incomplete colour orders'
       stop 1
    endif
    if (.not.allocated(pgl%ps)) then
       write (*,*) 'Cannot write an event without phase-space storage'
       stop 1
    endif
    if (size(pgl%ps).lt.1) then
       write (*,*) 'Cannot write an event with empty phase-space storage'
       stop 1
    endif
    if (.not.allocated(pgl%ps(1)%p)) then
       write (*,*) 'Cannot write an event without momenta'
       stop 1
    endif
    if (lbound(pgl%ps(1)%p,1).ne.0 .or. ubound(pgl%ps(1)%p,1).lt.3 .or. &
         size(pgl%ps(1)%p,1).ne.4 .or. &
         lbound(pgl%ps(1)%p,2).ne.1 .or. &
         size(pgl%ps(1)%p,2).ne.pgl%next) then
       write (*,*) 'Cannot write an event with invalid momentum dimensions'
       stop 1
    endif
    if (.not.all(ieee_is_finite(pgl%ps(1)%p(0:3,1:pgl%next))) .or. &
         .not.ieee_is_finite(wgt) .or. .not.ieee_is_finite(evt_sign)) then
       write (*,*) 'Cannot write an event with non-finite kinematics or weight'
       stop 1
    endif
    if (wgt.lt.0d0 .or. wgt.gt.event_value_limit .or. abs(evt_sign).ne.1d0 .or. &
         any(abs(pgl%ps(1)%p(0:3,1:pgl%next)).gt.event_value_limit)) then
       write (*,*) 'Cannot write an event with invalid magnitude/sign',wgt,evt_sign
       stop 1
    endif
    if (iproc_picked.lt.1 .or. iproc_picked.gt.pgl%nproc) then
       write (*,*) 'Cannot write an event with invalid process index',iproc_picked
       stop 1
    endif
    if (pgl%iden_iproc(iproc_picked).lt.1 .or. &
         pgl%iden_iproc(iproc_picked).gt.size(pgl%iden_processes,2)) then
       write (*,*) 'Cannot write an event with an invalid identical-process count',&
            pgl%iden_iproc(iproc_picked)
       stop 1
    endif
    if (iproc_iden_picked.lt.1 .or. &
         iproc_iden_picked.gt.pgl%iden_iproc(iproc_picked)) then
       write (*,*) 'Cannot write an event with invalid identical-process index',iproc_iden_picked
       stop 1
    endif
    if (hel_picked(1).lt.1 .or. hel_picked(2).lt.1) then
       write (*,*) 'Cannot write an event with invalid helicity selection',hel_picked
       stop 1
    endif
    if (keep_processes_separate) then
       if (.not.allocated(pgl%amps)) then
          write (*,*) 'Cannot write an event without process amplitudes'
          stop 1
       endif
       if (iproc_picked.gt.size(pgl%amps)) then
          write (*,*) 'Cannot write an event without its process amplitude',iproc_picked
          stop 1
       endif
       if (.not.allocated(pgl%amps(iproc_picked)%spins)) then
          write (*,*) 'Cannot write an event without process helicities',iproc_picked
          stop 1
       endif
       if (lbound(pgl%amps(iproc_picked)%spins,1).ne.1 .or. &
            lbound(pgl%amps(iproc_picked)%spins,2).ne.1 .or. &
            lbound(pgl%amps(iproc_picked)%spins,3).ne.1 .or. &
            size(pgl%amps(iproc_picked)%spins,1).ne.pgl%next .or. &
            hel_picked(1).gt.size(pgl%amps(iproc_picked)%spins,2) .or. &
            hel_picked(2).gt.size(pgl%amps(iproc_picked)%spins,3)) then
          write (*,*) 'Cannot write an event with out-of-range helicity indices',hel_picked
          stop 1
       endif
    else
       if (.not.allocated(pgl%amps)) then
          write (*,*) 'Cannot write an event without amplitudes'
          stop 1
       endif
       if (size(pgl%amps).lt.1) then
          write (*,*) 'Cannot write an event without an amplitude'
          stop 1
       endif
       if (.not.allocated(pgl%amps(1)%spins)) then
          write (*,*) 'Cannot write an event without helicities'
          stop 1
       endif
       if (lbound(pgl%amps(1)%spins,1).ne.1 .or. &
            lbound(pgl%amps(1)%spins,2).ne.1 .or. &
            lbound(pgl%amps(1)%spins,3).ne.1 .or. &
            size(pgl%amps(1)%spins,1).ne.pgl%next .or. &
            hel_picked(1).gt.size(pgl%amps(1)%spins,2) .or. &
            hel_picked(2).gt.size(pgl%amps(1)%spins,3)) then
          write (*,*) 'Cannot write an event with out-of-range helicity indices',hel_picked
          stop 1
       endif
    endif
    ! determine LHEF info
    NUP=pgl%next
    IDPRUP=1
    XWGTUP=sign(wgt,evt_sign)
    SCALUP=scale_shower
    AQEDUP=alphaEW
    AQCDUP=alphaS
    if (.not.ieee_is_finite(SCALUP) .or. .not.ieee_is_finite(AQEDUP) .or. &
         .not.ieee_is_finite(AQCDUP)) then
       write (*,*) 'Cannot write an event with invalid scales or couplings',SCALUP,AQEDUP,AQCDUP
       stop 1
    endif
    if (SCALUP.lt.0d0 .or. AQEDUP.lt.0d0 .or. AQCDUP.lt.0d0) then
       write (*,*) 'Cannot write an event with invalid scales or couplings',SCALUP,AQEDUP,AQCDUP
       stop 1
    endif
    if (SCALUP.gt.event_value_limit .or. AQEDUP.gt.event_value_limit .or. &
         AQCDUP.gt.event_value_limit) then
       write (*,*) 'Cannot write an event with unsafe scales or couplings',SCALUP,AQEDUP,AQCDUP
       stop 1
    endif
    do i=1,pgl%next
       IDUP(i)=pgl%iden_processes(i,iproc_iden_picked,iproc_picked)
       if (.not.supported_event_particle(IDUP(i))) then
          write (*,*) 'Cannot write an event with an unsupported particle identity',&
               i,IDUP(i)
          stop 1
       endif
       if(i.le.2) then
          ISTUP(i)=-1
          MOTHUP(1:2,i)=0
       else
          ISTUP(i)=+1
          MOTHUP(1:2,i)=[1,2]
       endif
       PUP(1:3,i)=pgl%ps(1)%p(1:3,i)
       PUP(4,i)=pgl%ps(1)%p(0,i)
       PUP(5,i)=phys_model%get_mass(IDUP(i))
       if (.not.ieee_is_finite(PUP(5,i))) then
          write (*,*) 'Cannot write an event with an invalid particle mass',IDUP(i),PUP(5,i)
          stop 1
       endif
       if (PUP(5,i).lt.0d0) then
          write (*,*) 'Cannot write an event with an invalid particle mass',IDUP(i),PUP(5,i)
          stop 1
       endif
       if (PUP(4,i).lt.0d0 .or. PUP(5,i).gt.event_value_limit) then
          write (*,*) 'Cannot write an event with unsafe particle kinematics',IDUP(i),PUP(:,i)
          stop 1
       endif
       VTIMUP(i)=0d0
       if (keep_processes_separate) then
          SPINUP(i)=dble(pgl%amps(iproc_picked)%spins(i,hel_picked(1),hel_picked(2)))
       else
          SPINUP(i)=dble(pgl%amps(1)%spins(i,hel_picked(1),hel_picked(2)))
       endif
    enddo
    call get_col_info(pgl,ICOLUP)
    ! write event to file
    write (iunit,'(a)',iostat=ios,iomsg=io_message) '<event>'
    call require_event_io(ios,'writing the event opening tag',io_message)
    write(iunit,503,iostat=ios,iomsg=io_message)NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    call require_event_io(ios,'writing the event header',io_message)
    do i=1,NUP
       write(iunit,504,iostat=ios,iomsg=io_message)IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I), &
            ICOLUP(1,I),ICOLUP(2,I),&
            PUP(1,I),PUP(2,I),PUP(3,I),PUP(4,I),PUP(5,I),&
            VTIMUP(I),SPINUP(I)
       call require_event_io(ios,'writing an event particle',io_message)
    enddo
    write (iunit,506,iostat=ios,iomsg=io_message) &
         '#color',pgl%color_orders(1:pgl%next,iproc_picked)
    call require_event_io(ios,'writing the event colour order',io_message)
    write (iunit,'(a)',iostat=ios,iomsg=io_message) '</event>'
    call require_event_io(ios,'writing the event closing tag',io_message)
    process_selection_ready=.false.
    helicity_selection_ready=.false.
    active_event_selection_epoch=0_8
    pgl%event_selection_epoch=0_8
503 format(1x,i4,1x,i10,4(1x,es24.16e3))
504 format(1x,i10,1x,i4,4(1x,i10),5(1x,es26.17e3),2(1x,es16.8e3))
506 format(a,100i3)
  end subroutine write_event

  subroutine event_update_wgt(iunit,ounit,wgt,first_event)
    implicit none
    character(len=1024) :: string
    integer,intent(in) :: iunit,ounit
    real(kind=8),dimension(3),intent(in) :: wgt
    logical,intent(in) :: first_event
    integer :: ios,nbody,record_kind
    logical :: init_seen,input_open,output_open
    integer :: NUP,IDPRUP
    integer :: LPRUP,IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP
    real(kind=8) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,XSECUP,XERRUP,XMAXUP
    real(kind=8) :: EBMUP(2)
    character(len=20) :: input_form,input_action,output_form,output_action
    character(len=64) :: trailing_token
    character(len=256) :: io_message
    input_form=''
    input_action=''
    output_form=''
    output_action=''
    io_message=''
    inquire(unit=iunit,opened=input_open,form=input_form,action=input_action,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect the provisional-event input unit:',&
            iunit,ios,trim(io_message)
       stop 1
    endif
    if (.not.input_open .or. trim(input_form).ne.'FORMATTED' .or. &
         trim(input_action).eq.'WRITE') then
       write (*,*) 'Cannot read provisional events from the requested unit:',&
            iunit,ios,trim(input_form),trim(input_action),trim(io_message)
       stop 1
    endif
    output_form=''
    output_action=''
    io_message=''
    inquire(unit=ounit,opened=output_open,form=output_form,action=output_action,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect the final-event output unit:',&
            ounit,ios,trim(io_message)
       stop 1
    endif
    if (.not.output_open .or. trim(output_form).ne.'FORMATTED' .or. &
         trim(output_action).eq.'READ') then
       write (*,*) 'Cannot write final events to the requested unit:',&
            ounit,ios,trim(output_form),trim(output_action),trim(io_message)
       stop 1
    endif
    if (.not.all(ieee_is_finite(wgt))) then
       write (*,*) 'Cannot update an event with invalid final weights',wgt
       stop 1
    endif
    if (any(wgt.lt.0d0)) then
       write (*,*) 'Cannot update an event with invalid final weights',wgt
       stop 1
    endif
    if (wgt(1).eq.0d0 .and. any(wgt(2:3).ne.0d0)) then
       write (*,*) 'Cannot discard an event with nonzero auxiliary weights',wgt
       stop 1
    endif
    if (first_event) then
       if (.not.all(ieee_is_finite(simple_integrator%res)) .or. &
            .not.all(ieee_is_finite(simple_integrator%unc))) then
          write (*,*) 'Cannot write non-finite integration results to the event header'
          stop 1
       endif
       if (simple_integrator%res(1).lt.0d0 .or. any(simple_integrator%unc.lt.0d0)) then
          write (*,*) 'Cannot write invalid integration results to the event header',&
               simple_integrator%res,simple_integrator%unc
          stop 1
       endif
       if (any(abs(simple_integrator%res).gt.event_value_limit) .or. &
            any(simple_integrator%unc.gt.event_value_limit)) then
          write (*,*) 'Cannot write unsafe integration results to the event header'
          stop 1
       endif
       init_seen=.false.
       do
          read(iunit,'(a)',iostat=ios) string
          if (ios.ne.0) then
             write (*,*) 'Unexpected end or read error before the first event',ios
             stop 1
          endif
          if (trim(adjustl(string)).eq.'<init>') then
             if (init_seen) then
                write (*,*) 'Found more than one init block before the first event'
                stop 1
             endif
             init_seen=.true.
             ! Update the cross section
             write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
             call require_event_io(ios,'copying the init opening tag',io_message)
             read(iunit,'(a)',iostat=ios) string
             if (ios.ne.0) then
                write (*,*) 'Could not read the LHE beam-init record',ios
                stop 1
             endif
             IDBMUP=huge(0)
             EBMUP=ieee_value(0d0,ieee_quiet_nan)
             PDFGUP=huge(0)
             PDFSUP=huge(0)
             IDWTUP=huge(0)
             NPRUP=huge(0)
             trailing_token=''
             read(string,*,iostat=ios) IDBMUP,EBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,&
                  trailing_token
             call validate_lhe_trailing_parse(ios,trailing_token,string,&
                  'beam-init record')
             if (any(IDBMUP.eq.0) .or. any(IDBMUP.eq.huge(0)) .or. &
                  any(PDFGUP.eq.huge(0)) .or. any(PDFSUP.eq.huge(0)) .or. &
                  IDWTUP.eq.0 .or. IDWTUP.lt.-4 .or. IDWTUP.gt.4 .or. &
                  NPRUP.ne.1) then
                write (*,*) 'Invalid discrete values in the LHE beam-init record:',&
                     trim(string)
                stop 1
             endif
             if (.not.all(ieee_is_finite(EBMUP))) then
                write (*,*) 'Invalid beam energies in the LHE beam-init record:',&
                     trim(string)
                stop 1
             endif
             if (any(EBMUP.le.0d0) .or. any(EBMUP.gt.event_value_limit)) then
                write (*,*) 'Invalid beam energies in the LHE beam-init record:',&
                     trim(string)
                stop 1
             endif
             write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
             call require_event_io(ios,'copying the LHE beam-init record',io_message)
             read(iunit,'(a)',iostat=ios) string
             if (ios.ne.0) then
                write (*,*) 'Could not read the LHE subprocess-init record',ios
                stop 1
             endif
             XSECUP=ieee_value(0d0,ieee_quiet_nan)
             XERRUP=ieee_value(0d0,ieee_quiet_nan)
             XMAXUP=ieee_value(0d0,ieee_quiet_nan)
             LPRUP=huge(0)
             trailing_token=''
             read(string,*,iostat=ios) XSECUP,XERRUP,XMAXUP,LPRUP,trailing_token
             call validate_lhe_trailing_parse(ios,trailing_token,string,&
                  'subprocess-init record')
             if (LPRUP.eq.huge(0) .or. .not.ieee_is_finite(XSECUP) .or. &
                  .not.ieee_is_finite(XERRUP) .or. .not.ieee_is_finite(XMAXUP)) then
                write (*,*) 'Invalid values in the LHE subprocess-init record:',&
                     trim(string)
                stop 1
             endif
             if (abs(XSECUP).gt.event_value_limit .or. XERRUP.lt.0d0 .or. &
                  XERRUP.gt.event_value_limit .or. XMAXUP.lt.0d0 .or. &
                  XMAXUP.gt.event_value_limit) then
                write (*,*) 'Unsafe values in the LHE subprocess-init record:',&
                     trim(string)
                stop 1
             endif
             XSECUP=simple_integrator%res(2)
             XERRUP=simple_integrator%unc(2)
             XMAXUP=simple_integrator%res(1)
             LPRUP=1
             write(ounit,502,iostat=ios,iomsg=io_message)XSECUP,XERRUP,XMAXUP,LPRUP
             call require_event_io(ios,'writing the LHE subprocess-init record',&
                  io_message)
             ! LHEF 3 permits metadata such as the generator declaration
             ! between the subprocess records and </init>.  Preserve those
             ! records verbatim; the writer itself emits such a declaration.
             do
                read(iunit,'(a)',iostat=ios) string
                if (ios.ne.0) then
                   write (*,*) 'Could not find the LHE init closing tag',ios
                   stop 1
                endif
                if (trim(adjustl(string)).eq.'</init>') exit
                if (trim(adjustl(string)).eq.'<event>' .or. &
                     trim(adjustl(string)).eq.'<init>') then
                   write (*,*) 'Malformed LHE init block before: ',trim(string)
                   stop 1
                endif
                write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
                call require_event_io(ios,'copying LHE init metadata',io_message)
             enddo
          endif
          if (trim(adjustl(string)).eq.'<event>') then
             if (.not.init_seen) then
                write (*,*) 'Found an event before the LHE init block'
                stop 1
             endif
             if (wgt(1).ne.0d0) then
                write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
                call require_event_io(ios,'copying the first event opening tag',&
                     io_message)
             endif
             exit
          endif
          if (trim(adjustl(string)).eq.'</LesHouchesEvents>') then
             write (*,*) 'Reached the end of the LHE file before the first event'
             stop 1
          endif
          write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
          call require_event_io(ios,'copying the LHE preamble',io_message)
       enddo
    else
       read(iunit,'(a)',iostat=ios) string
       if (ios.ne.0) then
          write (*,*) 'Could not read the next LHE event start',ios
          stop 1
       endif
       if (trim(adjustl(string)).ne.'<event>') then
          write (*,*) 'Expected an LHE event start but found: ',trim(string)
          stop 1
       endif
       if (wgt(1).ne.0d0) then
          write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
          call require_event_io(ios,'copying an event opening tag',io_message)
       endif
    endif
    read(iunit,'(a)',iostat=ios) string
    if (ios.ne.0) then
       write (*,*) 'Could not read the LHE event header',ios
       stop 1
    endif
    NUP=huge(0)
    IDPRUP=huge(0)
    XWGTUP=ieee_value(0d0,ieee_quiet_nan)
    SCALUP=ieee_value(0d0,ieee_quiet_nan)
    AQEDUP=ieee_value(0d0,ieee_quiet_nan)
    AQCDUP=ieee_value(0d0,ieee_quiet_nan)
    trailing_token=''
    read(string,*,iostat=ios) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,&
         trailing_token
    if (ios.gt.0) then
       write (*,*) 'Could not parse the LHE event header',ios,trim(string)
       stop 1
    endif
    if (ios.eq.0) then
       if (len_trim(trailing_token).eq.0 .or. trailing_token(1:1).ne.'#') then
          write (*,*) 'LHE event header contains trailing data: ',trim(string)
          stop 1
       endif
    endif
    if (NUP.lt.1 .or. NUP.gt.max_event_external_particles) then
       write (*,*) 'Invalid particle count in the LHE event header',NUP
       stop 1
    endif
    if (IDPRUP.ne.1 .or. .not.ieee_is_finite(XWGTUP) .or. &
         .not.ieee_is_finite(SCALUP) .or. .not.ieee_is_finite(AQEDUP) .or. &
         .not.ieee_is_finite(AQCDUP)) then
       write (*,*) 'Invalid numerical values in the LHE event header',XWGTUP,SCALUP,AQEDUP,AQCDUP
       stop 1
    endif
    if (SCALUP.lt.0d0 .or. AQEDUP.lt.0d0 .or. AQCDUP.lt.0d0 .or. &
         abs(XWGTUP).gt.event_value_limit .or. SCALUP.gt.event_value_limit .or. &
         AQEDUP.gt.event_value_limit .or. AQCDUP.gt.event_value_limit) then
       write (*,*) 'Invalid numerical values in the LHE event header',XWGTUP,SCALUP,AQEDUP,AQCDUP
       stop 1
    endif
    if (wgt(1).ne.0d0 .and. XWGTUP.eq.0d0) then
       write (*,*) 'Cannot preserve the sign of an event whose provisional weight is zero'
       stop 1
    endif
    ! The provisional event weight carries the subprocess/PDF sign.  The
    ! integrator returns positive absolute-envelope weights, so replacing the
    ! value outright would silently turn negative-weight events positive.
    XWGTUP=sign(wgt(1),XWGTUP)
    if (wgt(1).ne.0d0) then
       write(ounit,503,iostat=ios,iomsg=io_message) &
            NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
       call require_event_io(ios,'writing an updated event header',io_message)
    endif
    nbody=0
    do
       read(iunit,'(a)',iostat=ios) string
       if (ios.ne.0) then
          write (*,*) 'Unexpected end or read error inside an LHE event',ios
          stop 1
       endif
       if (trim(adjustl(string)).eq.'</event>') exit
       if (trim(adjustl(string)).eq.'<event>' .or. &
            trim(adjustl(string)).eq.'<init>' .or. &
            trim(adjustl(string)).eq.'</init>' .or. &
            trim(adjustl(string)).eq.'</LesHouchesEvents>' .or. &
            index(trim(adjustl(string)),'<LesHouchesEvents ').eq.1) then
          write (*,*) 'Malformed or unterminated LHE event near: ',trim(string)
          stop 1
       endif
       record_kind=classify_lhe_event_record(string,NUP)
       if (record_kind.lt.0) then
          write (*,*) 'Malformed LHE particle record: ',trim(string)
          stop 1
       elseif (record_kind.gt.0) then
          nbody=nbody+1
          if (nbody.gt.NUP) then
             write (*,*) 'LHE event contains more particle records than declared',nbody,NUP
             stop 1
          endif
       endif
       if (wgt(1).ne.0d0) then
          write(ounit,'(a)',iostat=ios,iomsg=io_message) trim(string)
          call require_event_io(ios,'copying an event record',io_message)
       endif
    enddo
    if (nbody.ne.NUP) then
       write (*,*) 'LHE event particle-record count differs from NUP',nbody,NUP
       stop 1
    endif
    if (wgt(1).ne.0d0) then
       write (ounit,505,iostat=ios,iomsg=io_message) '#overwgt',wgt(1:3)
       call require_event_io(ios,'writing final event weights',io_message)
       write (ounit,'(a)',iostat=ios,iomsg=io_message) '</event>'
       call require_event_io(ios,'writing an event closing tag',io_message)
    endif
502 format(3(1x,es24.16e3),1x,i10)
503 format(1x,i4,1x,i10,4(1x,es24.16e3))
505 format(a,3(1x,es24.16e3))
  end subroutine event_update_wgt

  integer function classify_lhe_event_record(line,nup) result(record_kind)
    ! Classify a line inside an LHE <event> block.  Comments, blank lines and
    ! optional XML metadata return zero; a validated particle record returns
    ! one; any other plain record is malformed and returns minus one.
    implicit none
    character(len=*),intent(in) :: line
    integer,intent(in) :: nup
    character(len=len(line)) :: stripped
    integer :: ios,idup,istup,mothup(2),icolup(2)
    real(kind=8) :: pup(5),vtimup,spinup
    character(len=64) :: trailing_token

    record_kind=-1
    if (nup.lt.1 .or. nup.gt.max_event_external_particles) return
    stripped=adjustl(line)
    if (len_trim(stripped).eq.0) then
       record_kind=0
       return
    endif
    if (stripped(1:1).eq.'#' .or. stripped(1:1).eq.'<') then
       record_kind=0
       return
    endif
    idup=huge(0)
    istup=huge(0)
    mothup=-1
    icolup=-1
    pup=ieee_value(0d0,ieee_quiet_nan)
    vtimup=ieee_value(0d0,ieee_quiet_nan)
    spinup=ieee_value(0d0,ieee_quiet_nan)
    trailing_token=''
    read(line,*,iostat=ios) idup,istup,mothup,icolup,pup,vtimup,spinup,&
         trailing_token
    if (ios.gt.0) return
    if (ios.eq.0) then
       if (len_trim(trailing_token).eq.0) return
       if (trailing_token(1:1).ne.'#') return
    endif
    if (idup.eq.0 .or. idup.eq.huge(0) .or. istup.eq.0 .or. &
         istup.lt.-9 .or. istup.gt.9) return
    if (any(mothup.lt.0) .or. any(mothup.gt.nup) .or. any(icolup.lt.0)) return
    if (.not.all(ieee_is_finite(pup)) .or. .not.ieee_is_finite(vtimup) .or. &
         .not.ieee_is_finite(spinup)) return
    if (any(abs(pup).gt.event_value_limit) .or. &
         abs(vtimup).gt.event_value_limit .or. &
         abs(spinup).gt.event_value_limit) return
    if (pup(4).lt.0d0 .or. pup(5).lt.0d0) return
    record_kind=1
  end function classify_lhe_event_record
  
  subroutine unwgt_process(pgl,iint)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer,intent(in) :: iint
    integer :: i,iproc,iproc_start,iproc_end
    real(kind=8) :: random,accum,target,contribution
    logical :: selected
    process_selection_ready=.false.
    helicity_selection_ready=.false.
    iproc_picked=0
    iproc_iden_picked=0
    hel_picked=0
    evt_sign=0d0
    active_event_selection_epoch=0_8
    pgl%event_selection_epoch=0_8
    if (pgl%next.lt.3 .or. pgl%next.gt.max_event_external_particles .or. &
         pgl%nproc.lt.1 .or. pgl%nproc.gt.max_event_process_records) then
       write (*,*) 'Cannot unweight a process group with invalid dimensions',&
            pgl%next,pgl%nproc
       stop 1
    endif
    if (.not.allocated(pgl%iden_iproc)) then
       write (*,*) 'Cannot unweight without identical-process counts'
       stop 1
    endif
    if (lbound(pgl%iden_iproc,1).ne.1 .or. &
         size(pgl%iden_iproc).ne.pgl%nproc) then
       write (*,*) 'Cannot unweight incomplete identical-process counts'
       stop 1
    endif
    if (.not.allocated(pgl%val_procs)) then
       write (*,*) 'Cannot unweight without subprocess contributions'
       stop 1
    endif
    if (lbound(pgl%val_procs,1).ne.1 .or. &
         lbound(pgl%val_procs,2).ne.1 .or. &
         size(pgl%val_procs,1).lt.1 .or. &
         size(pgl%val_procs,2).ne.pgl%nproc) then
       write (*,*) 'Cannot unweight incomplete subprocess contributions'
       stop 1
    endif
    if (.not.allocated(pgl%iden_processes)) then
       write (*,*) 'Cannot unweight without identical-process identities'
       stop 1
    endif
    if (lbound(pgl%iden_processes,1).ne.1 .or. &
         lbound(pgl%iden_processes,2).ne.1 .or. &
         lbound(pgl%iden_processes,3).ne.1 .or. &
         size(pgl%iden_processes,1).ne.pgl%next .or. &
         size(pgl%iden_processes,2).ne.size(pgl%val_procs,1) .or. &
         size(pgl%iden_processes,3).ne.pgl%nproc) then
       write (*,*) 'Cannot unweight incomplete identical-process identities'
       stop 1
    endif
    do iproc=1,pgl%nproc
       if (pgl%iden_iproc(iproc).lt.1 .or. &
            pgl%iden_iproc(iproc).gt.size(pgl%val_procs,1)) then
          write (*,*) 'Cannot unweight an invalid identical-process count',iproc,pgl%iden_iproc(iproc)
          stop 1
       endif
    enddo
    if (keep_processes_separate) then
       if (iint.lt.1 .or. iint.gt.pgl%nproc) then
          write (*,*) 'Cannot unweight an invalid process index',iint,pgl%nproc
          stop 1
       endif
       iproc_start=iint
       iproc_end=iint
    else
       iproc_start=1
       iproc_end=pgl%nproc
    endif
    target=0d0
    do iproc=iproc_start,iproc_end
       if (.not.all(ieee_is_finite(pgl%val_procs(1:pgl%iden_iproc(iproc),iproc)))) then
          write (*,*) 'Cannot unweight non-finite subprocess contributions',iproc
          stop 1
       endif
       do i=1,pgl%iden_iproc(iproc)
          contribution=abs(pgl%val_procs(i,iproc))
          if (contribution.gt.huge(1d0)-target) then
             write (*,*) 'Cannot unweight subprocess contributions whose total overflows',iproc
             stop 1
          endif
          target=target+contribution
       enddo
    enddo
    if (.not.ieee_is_finite(target)) then
       write (*,*) 'Cannot unweight a point with a non-finite total absolute weight',target
       stop 1
    endif
    if (target.le.0d0) then
       write (*,*) 'Cannot unweight a point with non-positive total absolute weight',target
       stop 1
    endif
    random=ran2()
    if (.not.ieee_is_finite(random)) then
       write (*,*) 'Cannot unweight with a non-finite random number'
       stop 1
    endif
    if (random.lt.0d0 .or. random.ge.1d0) then
       write (*,*) 'Cannot unweight with a random number outside [0,1)',random
       stop 1
    endif
    random=random*target
    accum=0d0
    selected=.false.
    do iproc=iproc_start,iproc_end
       do i=1,pgl%iden_iproc(iproc)
          accum=accum+abs(pgl%val_procs(i,iproc))
          if (accum.gt.random .or. accum.ge.target) then
             iproc_picked=iproc
             iproc_iden_picked=i
             selected=.true.
             exit
          endif
       enddo
       if (selected) exit
    enddo
    if (.not.selected) then
       write (*,*) 'Could not select a subprocess during unweighting',random,accum,target
       stop 1
    endif
    if (pgl%val_procs(iproc_iden_picked,iproc_picked).lt.0d0) then
       evt_sign=-1d0
    else
       evt_sign=+1d0
    endif
    if (keep_processes_separate .and. iproc_picked.ne.iint) then
       write (*,*) 'Could not unweight process correctly (keep_processes_separate=true)'
       stop 1
    endif
    if (event_selection_counter.eq.huge(event_selection_counter)) then
       write (*,*) 'Event-selection sequence counter exhausted'
       stop 1
    endif
    event_selection_counter=event_selection_counter+1_8
    active_event_selection_epoch=event_selection_counter
    pgl%event_selection_epoch=event_selection_counter
    process_selection_ready=.true.
  end subroutine unwgt_process

  subroutine unwgt_helicity(pgl)
    implicit none
    type(phase_space_order_group),intent(inout) :: pgl
    integer :: i,istart,iend,multiplicity,amplitude_index,offset_index,&
         multiplicity_column,expected_amplitudes,expected_offsets
    integer(kind=8) :: helicity_copy
    real(kind=8) :: random,total,accum
    logical :: selected
    helicity_selection_ready=.false.
    hel_picked=0
    if (.not.process_selection_ready) then
       write (*,*) 'Cannot unweight helicity without a fresh subprocess selection'
       stop 1
    endif
    if (active_event_selection_epoch.le.0_8 .or. &
         pgl%event_selection_epoch.ne.active_event_selection_epoch) then
       write (*,*) 'Cannot unweight helicity for a different or stale process group'
       stop 1
    endif
    if (pgl%next.lt.3 .or. pgl%next.gt.max_event_external_particles .or. &
         pgl%nproc.lt.1 .or. pgl%nproc.gt.max_event_process_records .or. &
         iproc_picked.lt.1 .or. iproc_picked.gt.pgl%nproc) then
       write (*,*) 'Cannot unweight helicity for an invalid subprocess',iproc_picked
       stop 1
    endif
    if (.not.allocated(pgl%amps)) then
       write (*,*) 'Cannot unweight helicity without amplitudes'
       stop 1
    endif
    if (.not.allocated(pgl%amp2_hel)) then
       write (*,*) 'Cannot unweight helicity without helicity contributions'
       stop 1
    endif
    if (.not.allocated(pgl%hel_fac)) then
       write (*,*) 'Cannot unweight helicity without multiplicities'
       stop 1
    endif
    if (keep_processes_separate) then
       expected_amplitudes=pgl%nproc
       expected_offsets=2
       amplitude_index=iproc_picked
       offset_index=1
       multiplicity_column=iproc_picked
    else
       expected_amplitudes=1
       expected_offsets=pgl%nproc+1
       amplitude_index=1
       offset_index=iproc_picked
       multiplicity_column=1
    endif
    if (lbound(pgl%amps,1).ne.1 .or. &
         size(pgl%amps).ne.expected_amplitudes) then
       write (*,*) 'Cannot unweight helicity with incompatible amplitude storage',&
            size(pgl%amps),expected_amplitudes
       stop 1
    endif
    if (lbound(pgl%amp2_hel,1).ne.1 .or. size(pgl%amp2_hel).lt.1 .or. &
         lbound(pgl%hel_fac,1).ne.1 .or. lbound(pgl%hel_fac,2).ne.1 .or. &
         size(pgl%hel_fac,1).ne.size(pgl%amp2_hel) .or. &
         size(pgl%hel_fac,2).ne.expected_amplitudes) then
       write (*,*) 'Cannot unweight helicity with incompatible contribution storage'
       stop 1
    endif
    if (.not.allocated(pgl%amps(amplitude_index)%iproc_start)) then
       write (*,*) 'Cannot unweight helicity without process offsets',amplitude_index
       stop 1
    endif
    if (lbound(pgl%amps(amplitude_index)%iproc_start,1).ne.1 .or. &
         size(pgl%amps(amplitude_index)%iproc_start).ne.expected_offsets .or. &
         pgl%amps(amplitude_index)%nprocs.ne.expected_offsets-1 .or. &
         pgl%amps(amplitude_index)%n_amps.lt.1 .or. &
         pgl%amps(amplitude_index)%n_amps.eq.huge(0) .or. &
         pgl%amps(amplitude_index)%n_amps.gt.size(pgl%amp2_hel)) then
       write (*,*) 'Cannot unweight helicity with incomplete process offsets',&
            amplitude_index
       stop 1
    endif
    if (pgl%amps(amplitude_index)%iproc_start(1).ne.1 .or. &
         pgl%amps(amplitude_index)%iproc_start(expected_offsets).ne.&
         pgl%amps(amplitude_index)%n_amps+1 .or. &
         any(pgl%amps(amplitude_index)%iproc_start(2:expected_offsets).lt.&
         pgl%amps(amplitude_index)%iproc_start(1:expected_offsets-1))) then
       write (*,*) 'Cannot unweight helicity with invalid process offsets',amplitude_index
       stop 1
    endif
    if (.not.allocated(pgl%amps(amplitude_index)%spins)) then
       write (*,*) 'Cannot unweight helicity without spin-copy storage',amplitude_index
       stop 1
    endif
    if (lbound(pgl%amps(amplitude_index)%spins,1).ne.1 .or. &
         lbound(pgl%amps(amplitude_index)%spins,2).ne.1 .or. &
         lbound(pgl%amps(amplitude_index)%spins,3).ne.1 .or. &
         size(pgl%amps(amplitude_index)%spins,1).ne.pgl%next .or. &
         size(pgl%amps(amplitude_index)%spins,2).lt.1 .or. &
         size(pgl%amps(amplitude_index)%spins,3).ne.&
         pgl%amps(amplitude_index)%n_amps) then
       write (*,*) 'Cannot unweight helicity with invalid spin-copy storage',amplitude_index
       stop 1
    endif
    istart=pgl%amps(amplitude_index)%iproc_start(offset_index)
    iend=pgl%amps(amplitude_index)%iproc_start(offset_index+1)-1
    if (istart.lt.1 .or. iend.lt.istart .or. iend.gt.size(pgl%amp2_hel)) then
       write (*,*) 'Cannot unweight an invalid helicity range',istart,iend,size(pgl%amp2_hel)
       stop 1
    endif
    if (.not.all(ieee_is_finite(pgl%amp2_hel(istart:iend)))) then
       write (*,*) 'Cannot unweight invalid helicity contributions',istart,iend
       stop 1
    endif
    if (any(pgl%amp2_hel(istart:iend).lt.0d0)) then
       write (*,*) 'Cannot unweight negative helicity contributions',istart,iend
       stop 1
    endif
    total=0d0
    do i=istart,iend
       if (pgl%amp2_hel(i).gt.huge(1d0)-total) then
          write (*,*) 'Cannot unweight helicity contributions whose total overflows',istart,iend
          stop 1
       endif
       total=total+pgl%amp2_hel(i)
    enddo
    if (.not.ieee_is_finite(total)) then
       write (*,*) 'Cannot unweight a non-finite helicity total',total
       stop 1
    endif
    if (total.le.0d0) then
       write (*,*) 'Cannot unweight a non-positive helicity total',total
       stop 1
    endif
    random=ran2()
    if (.not.ieee_is_finite(random)) then
       write (*,*) 'Cannot unweight helicity with a non-finite random number'
       stop 1
    endif
    if (random.lt.0d0 .or. random.ge.1d0) then
       write (*,*) 'Cannot unweight helicity with a random number outside [0,1)',random
       stop 1
    endif
    random=random*total
    accum=0d0
    selected=.false.
    do i=istart,iend
       accum=accum+pgl%amp2_hel(i)
       if (accum.gt.random .or. accum.ge.total) then
          selected=.true.
          exit
       endif
    enddo
    if (.not.selected) then
       write (*,*) 'Could not select a helicity during unweighting',random,accum,total
       stop 1
    endif
    hel_picked(2)=i
    if (hel_picked(2).lt.istart .or. hel_picked(2).gt.iend) then
       write (*,*) 'Could not unweight helicity',hel_picked,iproc_picked,istart,iend
       stop 1
    endif
    multiplicity=pgl%hel_fac(hel_picked(2),multiplicity_column)
    if (multiplicity.lt.1 .or. &
         multiplicity.gt.size(pgl%amps(amplitude_index)%spins,2)) then
       write (*,*) 'Cannot unweight an invalid helicity multiplicity',&
            multiplicity,size(pgl%amps(amplitude_index)%spins,2)
       stop 1
    endif
    if (multiplicity.gt.1) then
       random=ran2()
       if (.not.ieee_is_finite(random)) then
          write (*,*) 'Cannot select a helicity copy with a non-finite random number'
          stop 1
       endif
       if (random.lt.0d0 .or. random.ge.1d0) then
          write (*,*) 'Cannot select a helicity copy with a random number outside [0,1)',random
          stop 1
       endif
       helicity_copy=1_8+int(random*dble(multiplicity),kind=8)
       hel_picked(1)=int(min(int(multiplicity,kind=8),helicity_copy))
    else
       hel_picked(1)=1
    endif
    helicity_selection_ready=.true.
  end subroutine unwgt_helicity

  subroutine write_unique_in_file(iunit,pgl_unique,unique_map,unique_map_value,nevents)
    implicit none
    type(phase_space_order_group),allocatable,intent(in) :: pgl_unique
    real(kind=8),dimension(:),intent(in) :: unique_map_value
    integer,dimension(:),intent(in) :: unique_map
    integer,intent(in) :: iunit,nevents
    integer :: iproc,ios,label,mapped_process
    integer :: IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP
    real(kind=8) :: EBMUP(2),XSECUP,XERRUP,XMAXUP
    integer(kind=8) iseed
    logical :: unit_open
    character(len=20) :: unit_form,unit_action
    character(len=256) :: io_message
    common /to_seed/iseed
    unit_form=''
    unit_action=''
    io_message=''
    inquire(unit=iunit,opened=unit_open,form=unit_form,action=unit_action,&
         iostat=ios,iomsg=io_message)
    if (ios.ne.0) then
       write (*,*) 'Could not inspect the requested event-header unit:',&
            iunit,ios,trim(io_message)
       stop 1
    endif
    if (.not.unit_open .or. trim(unit_form).ne.'FORMATTED' .or. &
         trim(unit_action).eq.'READ') then
       write (*,*) 'Cannot write an event header to the requested unit:',&
            iunit,ios,trim(unit_form),trim(unit_action),trim(io_message)
       stop 1
    endif
    if (.not.allocated(pgl_unique)) then
       write (*,*) 'Cannot write an event header without unique processes'
       stop 1
    endif
    if (pgl_unique%next.lt.4 .or. pgl_unique%next.gt.max_event_external_particles .or. &
         pgl_unique%nproc.lt.1 .or. &
         pgl_unique%nproc.gt.max_event_process_records) then
       write (*,*) 'Cannot write an event header with invalid process dimensions',&
            pgl_unique%next,pgl_unique%nproc
       stop 1
    endif
    if (.not.allocated(pgl_unique%processes)) then
       write (*,*) 'Cannot write an event header without unique process identities'
       stop 1
    endif
    if (lbound(pgl_unique%processes,1).ne.1 .or. &
         lbound(pgl_unique%processes,2).ne.1 .or. &
         size(pgl_unique%processes,1).ne.pgl_unique%next .or. &
         size(pgl_unique%processes,2).ne.pgl_unique%nproc) then
       write (*,*) 'Cannot write an event header with incomplete unique process identities'
       stop 1
    endif
    if (lbound(unique_map,1).ne.1 .or. lbound(unique_map_value,1).ne.1 .or. &
         size(unique_map).ne.pgl_unique%nproc .or. &
         size(unique_map_value).ne.pgl_unique%nproc) then
       write (*,*) 'Cannot write an event header with inconsistent unique-process maps'
       stop 1
    endif
    if (.not.all(ieee_is_finite(unique_map_value))) then
       write (*,*) 'Cannot write an event header with invalid unique-process maps'
       stop 1
    endif
    if (any(unique_map.lt.-1) .or. any(unique_map.gt.pgl_unique%nproc) .or. &
         any(unique_map_value.lt.0d0) .or. &
         any(unique_map_value.gt.event_value_limit)) then
       write (*,*) 'Cannot write an event header with invalid unique-process maps'
       stop 1
    endif
    do iproc=1,pgl_unique%nproc
       do label=1,pgl_unique%next
          if (.not.supported_event_particle(pgl_unique%processes(label,iproc))) then
             write (*,*) 'Cannot write an event header with an unsupported particle identity',&
                  iproc,label,pgl_unique%processes(label,iproc)
             stop 1
          endif
       enddo
       mapped_process=unique_map(iproc)
       if (mapped_process.eq.-1) then
          if (unique_map_value(iproc).ne.1d0) then
             write (*,*) 'Canonical unique process has an invalid map factor',iproc,&
                  unique_map_value(iproc)
             stop 1
          endif
       elseif (mapped_process.eq.0) then
          if (unique_map_value(iproc).ne.0d0) then
             write (*,*) 'Zero unique process has a nonzero map factor',iproc,&
                  unique_map_value(iproc)
             stop 1
          endif
       else
          if (mapped_process.ge.iproc .or. unique_map(mapped_process).ne.-1 .or. &
               unique_map_value(iproc).le.0d0) then
             write (*,*) 'Unique-process map is cyclic or non-canonical',iproc,mapped_process
             stop 1
          endif
          do label=1,pgl_unique%next
             if (phys_model%get_mass(pgl_unique%processes(label,iproc)).ne.&
                  phys_model%get_mass(pgl_unique%processes(label,mapped_process))) then
                write (*,*) 'Unique-process map changes the external mass layout',&
                     iproc,mapped_process,label
                stop 1
             endif
          enddo
       endif
    enddo
    if (.not.ieee_is_finite(sqrts)) then
       write (*,*) 'Cannot write an event header with invalid collision energy',sqrts
       stop 1
    endif
    if (sqrts.le.0d0 .or. sqrts.gt.event_value_limit) then
       write (*,*) 'Cannot write an event header with invalid collision energy',sqrts
       stop 1
    endif
    if (nevents.lt.1) then
       write (*,*) 'Cannot write an event header with a non-positive event count',nevents
       stop 1
    endif
    if (iseed.lt.0_8 .or. iseed.gt.904866561_8) then
       write (*,*) 'Cannot write an event header with an invalid random seed',iseed
       stop 1
    endif
    IDBMUP(1:2)=2212     ! two protons
    EBMUP(1:2)=sqrts/2d0 ! half of collision energy
    PDFGUP(1:2)=-1
    PDFSUP(1:2)=pdf_lhaid
    IDWTUP=-3
    NPRUP=1
    XSECUP=0d0
    XERRUP=0d0
    XMAXUP=0d0
    LPRUP=1
    write(iunit,'(a)',iostat=ios,iomsg=io_message) &
         '<LesHouchesEvents version="3.0">'
    call require_event_io(ios,'writing the LHE opening tag',io_message)
    write(iunit,'(a)',iostat=ios,iomsg=io_message) '<header>'
    call require_event_io(ios,'writing the LHE header opening tag',io_message)
    write(iunit,*,iostat=ios,iomsg=io_message) pgl_unique%next,pgl_unique%nproc
    call require_event_io(ios,'writing unique-process dimensions',io_message)
    do iproc=1,pgl_unique%nproc
       write(iunit,*,iostat=ios,iomsg=io_message) &
            unique_map(iproc),unique_map_value(iproc),&
            pgl_unique%processes(1:pgl_unique%next,iproc)
       call require_event_io(ios,'writing a unique process',io_message)
    enddo
    write(iunit,'(a,1x,i12,1x,a)',iostat=ios,iomsg=io_message) &
         '<nevents>',nevents,'</nevents>'
    call require_event_io(ios,'writing the requested event count',io_message)
    write(iunit,'(a,1x,i12,1x,a)',iostat=ios,iomsg=io_message) &
         '<seed>   ',iseed,  '</seed>'
    call require_event_io(ios,'writing the event seed',io_message)
    write(iunit,'(a)',iostat=ios,iomsg=io_message) '</header>'
    call require_event_io(ios,'writing the LHE header closing tag',io_message)
    write(iunit,'(a)',iostat=ios,iomsg=io_message) '<init>'
    call require_event_io(ios,'writing the LHE init opening tag',io_message)
    write(iunit,501,iostat=ios,iomsg=io_message) &
         IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),PDFGUP(1)&
         &,PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
    call require_event_io(ios,'writing the LHE beam-init record',io_message)
    write(iunit,502,iostat=ios,iomsg=io_message)XSECUP,XERRUP,XMAXUP,LPRUP
    call require_event_io(ios,'writing the LHE subprocess-init record',io_message)
    write(iunit,'(a)',iostat=ios,iomsg=io_message) &
         "<generator name='AmpliCol' version='1.0'>please cite arXiv:2601.19483</generator>"
    call require_event_io(ios,'writing the LHE generator declaration',io_message)
    write(iunit,'(a)',iostat=ios,iomsg=io_message) '</init>'
    call require_event_io(ios,'writing the LHE init closing tag',io_message)
501 format(2(1x,i10),2(1x,es24.16e3),2(1x,i4),2(1x,i10),1x,i4,1x,i6)
502 format(3(1x,es24.16e3),1x,i10)
  end subroutine write_unique_in_file


  subroutine get_col_info(pgl,ICOLUP)
    implicit none
    type(phase_space_order_group),intent(in) :: pgl
    integer,dimension(:,:),intent(out) :: ICOLUP
    integer,dimension(max_event_external_particles) :: ipdg,ord
    integer :: i,label,ncoloured,nsteps
    ICOLUP=0
    if (pgl%next.lt.3 .or. pgl%next.gt.max_event_external_particles .or. &
         lbound(ICOLUP,1).ne.1 .or. size(ICOLUP,1).ne.2 .or. &
         lbound(ICOLUP,2).ne.1 .or. size(ICOLUP,2).lt.pgl%next) then
       write (*,*) 'Cannot construct event colours with invalid dimensions',&
            pgl%next,shape(ICOLUP)
       stop 1
    endif
    ! Take anti-particles for the initial state
    do i=1,2
       ipdg(i)=phys_model%get_antipart(pgl%iden_processes(i,iproc_iden_picked,iproc_picked))
    enddo
    ipdg(3:pgl%next)=pgl%iden_processes(3:pgl%next,iproc_iden_picked,iproc_picked)
    ord(1:pgl%next)=pgl%color_orders(1:pgl%next,iproc_picked)
    if (any(ord.lt.1) .or. any(ord.gt.pgl%next)) then
       write (*,*) 'Cannot construct event colours from an out-of-range order',ord
       stop 1
    endif
    do i=1,pgl%next
       if (count(ord.eq.i).ne.1) then
          write (*,*) 'Cannot construct event colours from a non-permutation',ord
          stop 1
       endif
    enddo
    ncoloured=0
    do i=1,pgl%next
       if (phys_model%is_quark(ipdg(i)) .or. phys_model%is_antiquark(ipdg(i)) .or. &
            phys_model%is_gluon(ipdg(i))) ncoloured=ncoloured+1
    enddo
    if (ncoloured.eq.0) return
    i=1
    label=501
    nsteps=0
    do
       nsteps=nsteps+1
       if (nsteps.gt.2*pgl%next+ncoloured) then
          write (*,*) 'Cannot close the event colour flow',ipdg,ord
          stop 1
       endif
       if (phys_model%is_quark(iPDG(ord(i))) .or. phys_model%is_gluon(iPDG(ord(i)))) then
          if (ICOLUP(1,ord(i)).ne.0) exit ! already filled. We have made the full loop and are done.
          ! found a colour line.
          ICOLUP(1,ord(i))=label
          ! find the corresponding anti-colour line. This should be
          ! the next particle in the colour ordering (as long as
          ! that's not a colour singlet.
          do
             nsteps=nsteps+1
             if (nsteps.gt.2*pgl%next+ncoloured) then
                write (*,*) 'Cannot find the matching event anti-colour',ipdg,ord
                stop 1
             endif
             i=mod(i,pgl%next)+1
             if (phys_model%is_antiquark(iPDG(ord(i))) .or. phys_model%is_gluon(iPDG(ord(i)))) then
                ICOLUP(2,ord(i))=label
                label=label+1
                exit
             elseif (phys_model%is_quark(iPDG(ord(i)))) then
                write (*,*) 'Cannot have a quark after another one'
                stop 1
             endif
          enddo
       else
          i=mod(i,pgl%next)+1
       endif
    enddo
    ! Flip initial states
    ICOLUP(1:2,1:2)=ICOLUP(2:1:-1,1:2)
  end subroutine get_col_info

  logical function supported_event_particle(ipdg)
    implicit none
    integer,intent(in) :: ipdg
    supported_event_particle=(ipdg.ge.1 .and. ipdg.le.6) .or. &
         (ipdg.le.-1 .and. ipdg.ge.-6) .or. &
         (ipdg.ge.11 .and. ipdg.le.16) .or. &
         (ipdg.le.-11 .and. ipdg.ge.-16) .or. &
         ipdg.eq.21 .or. ipdg.eq.22 .or. ipdg.eq.23 .or. &
         ipdg.eq.24 .or. ipdg.eq.-24 .or. ipdg.eq.25
  end function supported_event_particle

  subroutine validate_lhe_trailing_parse(io_status,trailing_token,line,record_name)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: trailing_token,line,record_name

    if (io_status.gt.0) then
       write (*,*) 'Could not parse the LHE ',trim(record_name),': ',trim(line)
       stop 1
    endif
    if (io_status.eq.0) then
       if (len_trim(trailing_token).eq.0 .or. trailing_token(1:1).ne.'#') then
          write (*,*) 'LHE ',trim(record_name),' contains trailing data: ',trim(line)
          stop 1
       endif
    endif
  end subroutine validate_lhe_trailing_parse

  subroutine require_event_io(io_status,operation,io_message)
    implicit none
    integer,intent(in) :: io_status
    character(len=*),intent(in) :: operation,io_message

    if (io_status.ne.0) then
       write (*,*) 'Event-file I/O failure while ',trim(operation),': ',&
            io_status,trim(io_message)
       stop 1
    endif
  end subroutine require_event_io
end module handling_events
