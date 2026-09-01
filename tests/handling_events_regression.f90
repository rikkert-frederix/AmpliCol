program handling_events_regression
  use common, only: simple_integrator
  use handling_processes, only: phase_space_order_group
  use handling_events, only: event_update_wgt,classify_lhe_event_record,&
       unwgt_process,unwgt_helicity
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer :: iunit,ounit,ios,nup,idprup,lprup
  integer(kind=8) :: iseed
  real(kind=8) :: xwgtup,scalup,aqedup,aqcdup,xsecup,xerrup,xmaxup
  real(kind=8) :: overwgt(3)
  character(len=1024) :: line
  character(len=64) :: mode
  logical :: found_init,found_event,found_overwgt,found_generator
  type(phase_space_order_group) :: selection_group,other_group
  common /to_seed/ iseed

  mode='success'
  if (command_argument_count().ge.1) call get_command_argument(1,mode)
  if (trim(mode).ne.'success') then
     select case (trim(mode))
     case ('nan-weight')
        call open_empty_units(iunit,ounit)
        call event_update_wgt(iunit,ounit,&
             [ieee_value(0d0,ieee_quiet_nan),0d0,0d0],.false.)
     case ('nan-results')
        call open_empty_units(iunit,ounit)
        simple_integrator%res=[ieee_value(0d0,ieee_quiet_nan),1d0]
        simple_integrator%unc=[0d0,0d0]
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.true.)
     case ('header-missing')
        call open_single_event(iunit,ounit,'3 1 -1 91 0.0073')
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.false.)
     case ('header-extra')
        call open_single_event(iunit,ounit,'3 1 -1 91 0.0073 0.118 7')
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.false.)
     case ('header-nan')
        call open_single_event(iunit,ounit,'3 1 NaN 91 0.0073 0.118')
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.false.)
     case ('header-huge')
        call open_single_event(iunit,ounit,'3 1 1e300 91 0.0073 0.118')
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.false.)
     case ('header-process')
        call open_single_event(iunit,ounit,'3 2 1 91 0.0073 0.118')
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.false.)
     case ('init-nproc')
        call open_bad_init(iunit,ounit,&
             '2212 2212 7000 7000 -1 -1 244800 244800 -3 2',&
             '1 0 1 1')
        simple_integrator%res=[1d0,1d0]
        simple_integrator%unc=[0d0,0d0]
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.true.)
     case ('init-nan-beam')
        call open_bad_init(iunit,ounit,&
             '2212 2212 NaN 7000 -1 -1 244800 244800 -3 1',&
             '1 0 1 1')
        simple_integrator%res=[1d0,1d0]
        simple_integrator%unc=[0d0,0d0]
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.true.)
     case ('init-extra')
        call open_bad_init(iunit,ounit,&
             '2212 2212 7000 7000 -1 -1 244800 244800 -3 1 extra',&
             '1 0 1 1')
        simple_integrator%res=[1d0,1d0]
        simple_integrator%unc=[0d0,0d0]
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.true.)
     case ('init-bad-process')
        call open_bad_init(iunit,ounit,&
             '2212 2212 7000 7000 -1 -1 244800 244800 -3 1',&
             '1 NaN 1 1')
        simple_integrator%res=[1d0,1d0]
        simple_integrator%unc=[0d0,0d0]
        call event_update_wgt(iunit,ounit,[1d0,0d0,0d0],.true.)
     case ('bad-multiplicity')
        call setup_selection_group(selection_group,2,1)
        call initialize_test_random_number()
        call unwgt_process(selection_group,1)
        call unwgt_helicity(selection_group)
     case ('shifted-process-bounds')
        call setup_selection_group(selection_group,1,1)
        deallocate(selection_group%iden_iproc)
        allocate(selection_group%iden_iproc(0:0))
        selection_group%iden_iproc=1
        call initialize_test_random_number()
        call unwgt_process(selection_group,1)
     case ('stale-process-group')
        call setup_selection_group(selection_group,1,1)
        call setup_selection_group(other_group,1,1)
        call initialize_test_random_number()
        call unwgt_process(selection_group,1)
        call unwgt_helicity(other_group)
     case default
        error stop 'unknown handling-events regression mode'
     end select
     error stop 'invalid handling-events input was unexpectedly accepted'
  endif

  if (classify_lhe_event_record('#color 1 2 3',3).ne.0) &
       error stop 'an LHE comment was misclassified as a particle'
  if (classify_lhe_event_record('<rwgt>',3).ne.0) &
       error stop 'LHE XML metadata were misclassified as a particle'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 200 91 0 0',3).ne.1) &
       error stop 'a valid LHE particle was rejected'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 200 91 0',3).ne.-1) &
       error stop 'an incomplete LHE particle was accepted'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 200 91 0 0 7',3).ne.-1) &
       error stop 'an LHE particle with trailing data was accepted'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 200 91 0 0 # comment',3).ne.1) &
       error stop 'an LHE particle with a trailing comment was rejected'
  if (classify_lhe_event_record(' malformed body text',3).ne.-1) &
       error stop 'malformed LHE body text was accepted'
  if (classify_lhe_event_record(' 23 1 1 9 0 0 0 0 0 200 91 0 0',3).ne.-1) &
       error stop 'an out-of-range LHE mother index was accepted'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 NaN 91 0 0',3).ne.-1) &
       error stop 'a non-finite LHE particle was accepted'
  if (classify_lhe_event_record(' 23 1 1 2 0 0 0 0 0 1e300 91 0 0',3).ne.-1) &
       error stop 'an unsafe finite LHE particle was accepted'

  simple_integrator%res=[3d0,-1.25d0]
  simple_integrator%unc=[0.2d0,0.125d0]

  open(newunit=iunit,status='scratch',action='readwrite',form='formatted')
  write(iunit,'(a)') '<LesHouchesEvents version="3.0">'
  write(iunit,'(a)') '<header>'
  write(iunit,'(a)') '</header>'
  write(iunit,'(a)') '<init>'
  write(iunit,*) 2212,2212,7000d0,7000d0,-1,-1,244800,244800,-3,1
  write(iunit,*) 99d0,9d0,999d0,1
  write(iunit,'(a)') "<generator name='AmpliCol' version='1.0'>test</generator>"
  write(iunit,'(a)') '</init>'
  write(iunit,'(a)') '<event>'
  write(iunit,*) 3,1,-1d0,91d0,7.3d-3,0.118d0
  write(iunit,'(a)') ' 1 -1 0 0 0 0 0 0 100 100 0 0 -1'
  write(iunit,'(a)') ' -1 -1 0 0 0 0 0 0 -100 100 0 0 1'
  write(iunit,'(a)') ' 23 1 1 2 0 0 0 0 0 200 91 0 0'
  write(iunit,'(a)') '#color 1 2 3'
  write(iunit,'(a)') '</event>'
  rewind(iunit)

  open(newunit=ounit,status='scratch',action='readwrite',form='formatted')
  call event_update_wgt(iunit,ounit,[2d0,2.5d0,0.25d0],.true.)
  rewind(ounit)

  found_init=.false.
  found_event=.false.
  found_overwgt=.false.
  found_generator=.false.
  do
     read(ounit,'(a)',iostat=ios) line
     if (ios.lt.0) exit
     if (ios.ne.0) error stop 'read failure in event-output regression'
     if (index(adjustl(line),'<init>').eq.1) then
        read(ounit,'(a)',iostat=ios) line
        if (ios.ne.0) error stop 'missing beam init record'
        read(ounit,*,iostat=ios) xsecup,xerrup,xmaxup,lprup
        if (ios.ne.0) error stop 'invalid rewritten subprocess init record'
        if (abs(xsecup+1.25d0).gt.1d-14 .or. abs(xerrup-0.125d0).gt.1d-14 .or. &
             abs(xmaxup-3d0).gt.1d-14 .or. lprup.ne.1) &
             error stop 'incorrect rewritten subprocess init values'
        found_init=.true.
     elseif (index(adjustl(line),'<event>').eq.1) then
        read(ounit,*,iostat=ios) nup,idprup,xwgtup,scalup,aqedup,aqcdup
        if (ios.ne.0) error stop 'invalid rewritten event header'
        if (nup.ne.3 .or. idprup.ne.1 .or. abs(xwgtup+2d0).gt.1d-14) &
             error stop 'negative provisional event sign was not preserved'
        found_event=.true.
     elseif (index(adjustl(line),'#overwgt').eq.1) then
        read(line(index(line,'#overwgt')+9:),*,iostat=ios) overwgt
        if (ios.ne.0) error stop 'invalid overweight record'
        if (maxval(abs(overwgt-[2d0,2.5d0,0.25d0])).gt.1d-14) &
             error stop 'incorrect overweight record'
        found_overwgt=.true.
     elseif (index(adjustl(line),'<generator ').eq.1) then
        found_generator=.true.
     endif
  enddo

  if (.not.found_init) error stop 'rewritten init block was not found'
  if (.not.found_event) error stop 'rewritten event was not found'
  if (.not.found_overwgt) error stop 'rewritten overweight record was not found'
  if (.not.found_generator) error stop 'generator metadata was not preserved'
  close(iunit)
  close(ounit)
  write(*,'(a)') 'Handling-events regression: PASS'

contains

  subroutine open_empty_units(input_unit,output_unit)
    integer,intent(out) :: input_unit,output_unit
    open(newunit=input_unit,status='scratch',action='readwrite',form='formatted')
    rewind(input_unit)
    open(newunit=output_unit,status='scratch',action='readwrite',form='formatted')
  end subroutine open_empty_units

  subroutine open_single_event(input_unit,output_unit,header)
    integer,intent(out) :: input_unit,output_unit
    character(len=*),intent(in) :: header
    open(newunit=input_unit,status='scratch',action='readwrite',form='formatted')
    write(input_unit,'(a)') '<event>'
    write(input_unit,'(a)') trim(header)
    rewind(input_unit)
    open(newunit=output_unit,status='scratch',action='readwrite',form='formatted')
  end subroutine open_single_event

  subroutine open_bad_init(input_unit,output_unit,beam_record,process_record)
    integer,intent(out) :: input_unit,output_unit
    character(len=*),intent(in) :: beam_record,process_record
    open(newunit=input_unit,status='scratch',action='readwrite',form='formatted')
    write(input_unit,'(a)') '<LesHouchesEvents version="3.0">'
    write(input_unit,'(a)') '<init>'
    write(input_unit,'(a)') trim(beam_record)
    write(input_unit,'(a)') trim(process_record)
    write(input_unit,'(a)') '</init>'
    write(input_unit,'(a)') '<event>'
    write(input_unit,'(a)') '3 1 1 91 0.0073 0.118'
    rewind(input_unit)
    open(newunit=output_unit,status='scratch',action='readwrite',form='formatted')
  end subroutine open_bad_init

  subroutine setup_selection_group(group,multiplicity,spin_copies)
    type(phase_space_order_group),intent(out) :: group
    integer,intent(in) :: multiplicity,spin_copies
    group%next=4
    group%nproc=1
    allocate(group%iden_iproc(1),group%val_procs(1,1),&
         group%iden_processes(4,1,1))
    group%iden_iproc=1
    group%val_procs=1d0
    group%iden_processes(:,1,1)=[1,-1,23,23]
    allocate(group%amps(1),group%amp2_hel(1),group%hel_fac(1,1))
    group%amps(1)%nprocs=1
    group%amps(1)%n_amps=1
    allocate(group%amps(1)%iproc_start(2),&
         group%amps(1)%spins(4,spin_copies,1))
    group%amps(1)%iproc_start=[1,2]
    group%amps(1)%spins=0
    group%amp2_hel=1d0
    group%hel_fac=multiplicity
  end subroutine setup_selection_group

  subroutine initialize_test_random_number()
    iseed=1_8
    open(unit=99,status='scratch',action='write')
  end subroutine initialize_test_random_number
end program handling_events_regression
