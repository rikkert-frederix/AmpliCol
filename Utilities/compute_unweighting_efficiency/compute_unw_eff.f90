! without fastjet:
      ! gfortran -o compute_unw_eff compute_unw_eff.f90 


program plot_events
  implicit none
  integer ifile,nPSpoints
  logical :: done
  character*140 :: filename
  character(:),allocatable :: outfile
  character*50 weights_info(10)
  double precision dummy
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8) :: evt_wgt_LC,evt_wgt_NLC,evt_wgt_full,rwgt_factor,unw_eff,max_wgt
  integer :: next
  integer :: argc
  character(len=456) :: argv
  integer :: i

  CALL GET_COMMAND_ARGUMENT(1, argv)
  read(argv,*) filename
  ifile=11
  write(*,*) filename
  open(unit=ifile,file=filename,status='OLD')
  allocate(character(len(trim(filename))) :: outfile)  ! Note the correct form
  outfile=filename(7:)
  do i=1,len(outfile)
     if (outfile(i:i).eq.'l') exit
  enddo
  outfile=adjustl(outfile(1:i-2))
  write(*,*) 'out_unwgt_'//outfile//'.txt'
  outfile='out_unwgt_'//outfile//'.txt'
  open(unit=20,file=outfile)
  max_wgt=0d0
  do
     call read_event(ifile,done)
     if (done) exit
     if (evt_wgt_full.lt.0d0) then
        write (*,*) 'Found negative weight. Stopping...'
        stop
     endif
     max_wgt=max(max_wgt,evt_wgt_full)
  enddo

  rewind(ifile)

  unw_eff=0d0
  nPSpoints=0
  do
     call read_event(ifile,done)
     if (done) exit
     nPSpoints=nPSpoints+1
     unw_eff=unw_eff+evt_wgt_full/max_wgt
  enddo

  write(20,*) 'unweighting efficiency is',unw_eff/dble(nPSpoints)
  !write (*,*) 'unweighting efficiency is',unw_eff/dble(nPSpoints)
  close (ifile)
  close(20)

contains
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
