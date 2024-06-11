! without fastjet:
      ! gfortran -fcheck=all -o plot_events HwU.f analysis.f plot_events.f90

program plot_events
  implicit none
  integer ifile,nPSpoints
  logical :: done
  character*140 filename
  character*50 weights_info(10)
  double precision dummy
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: evt_wgt_LC,evt_wgt_NLC,evt_wgt_full,rwgt_factor
  integer :: next
  character(len=256) :: argv

  CALL GET_COMMAND_ARGUMENT(1, argv)
  read(argv,'(a)') filename
  write(*,*) filename
  ifile=11
  open(unit=ifile,file=filename,status='OLD')

  weights_info(1)="central value               "
  rwgt_factor=1d0

  call set_error_estimation(0)
  call analysis_begin(1,weights_info)

  nPSpoints=0
  do
     call read_event(ifile,done)
     if (done) exit
     nPSpoints=nPSpoints+1
     call plot_event()
  enddo
  close (ifile)
  call finalize_histograms(nint(nPSpoints/rwgt_factor))
  call analysis_end(dummy)

contains

  subroutine plot_event()
    implicit none
    call analysis_fill(next,p,evt_wgt_LC,evt_wgt_NLC,evt_wgt_full)
  end subroutine plot_event

  subroutine read_event(iunit,done)
    implicit none
    integer :: i,iunit
    logical :: done
    character(len=100) :: dummy
    character(len=15) :: dummy2
    real(kind=8) :: dum,evt_wgt,wgt,amp2,weight
    done=.false.
    read (iunit,*,err=99,end=99) dummy
    if (index(dummy,'REWEIGHT_FACTOR').ne.0) then
       backspace(iunit)
       read(iunit,*) dummy2,rwgt_factor
       done=.true.
       return
    endif
    read (iunit,*,err=99,end=99) next,evt_wgt,wgt,amp2,weight
    read (iunit,*,err=99,end=99) dum
    read (iunit,*) dum,dum,dum
    read (iunit,*) evt_wgt_LC,evt_wgt_NLC,evt_wgt_full
    read (iunit,*) dum ! helicity
    if (.not.allocated(p)) allocate(p(0:3,next))
    do i=1,next
       read (iunit,*,err=99,end=99) dum,p(1:3,i),p(0,i)
    enddo
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event


end program plot_events

