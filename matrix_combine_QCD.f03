! gfortran -ffast-math -O3 -o matrix_reweight random.f color_algebra.f95 math_functions.f03 feynmanrules.f03 amplitude_QCD.f03 matrix_reweight.f03

module comb_events
  implicit none
  real(kind=8) :: wgt,evt_wgt,evt_scale,evt_aEW,evt_aS,weight,amp2,rwgt_NLC,rwgt_full
  real(kind=8) :: evt_NLC,evt_full
  real(kind=8) :: xsec, xerr, xsecabs, xsecerr
  integer :: integ
  integer,dimension(:),allocatable :: helicity,col_order,iPDG,istat
  real(kind=8),dimension(:,:),allocatable :: momenta
  integer,dimension(:,:),allocatable :: col_info,moth
  real(kind=8),dimension(:),allocatable :: imass,ivit
end module comb_events
module timings
  implicit none
  real(kind=4) :: tBefore,tAfter,tTot_A=0.,tTot_B=0.,t_amp=0.,t_amp_init=0.,&
       t_mat_LC=0.,t_mat_NLC=0.,t_mat_full=0.,t_all=0.,t_ran=0.
end module timings
module arguments
  implicit none
  integer :: c_o,c_o_t,c_o_i,c_o_j,c_o_k,imode
end module arguments

program matrix_combine
  use math_functions
  use amplitude_QCD_mod
  use timings
  use particles
  use comb_events
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
  real(kind=8),dimension(3) :: matrix2
  real(kind=8),dimension(:,:),allocatable :: p
  real(kind=8),dimension(:),allocatable :: unique_map_value,process_map_value
  complex(kind=8) :: amp2_c,amp_col_c
  logical :: done
  character(len=string_len) :: tag,tag_read,add_arg=''
  real(kind=8) :: abs_wgt,ra,ran,xsec_evt,sum_wgt,sum_unwgt
  logical :: accept_event
  real(kind=8) :: n_evt,n_tot,xsec_ave
  real(kind=8),parameter :: tot_evt_wanted=100
  integer :: nfiles

  call get_run_arguments()

  call cpu_time(tTot_B)

  n_evt=0
  do j=1,nfiles
     do
     call read_event(10+j,done)
     if (done) exit
     call write_event(100)
     n_evt=n_evt+1
     enddo
  enddo

  close(100)
  do j=1,nfiles
     close(10+j)
  enddo
  call cpu_time(tTot_a)
  t_all=tTot_a-tTot_b

  write(*,*) 'Total number of events in combined file:',n_evt
  write(*,*) 'Total time:',t_all
contains  


  subroutine get_run_arguments()
    use arguments
    implicit none
    integer :: argc,i,k
    character(len=256) :: argv,filename
    character(len=256),dimension(:),allocatable :: files
    ! integration steps:
    ! imode=0  (Setting up grids)
    ! imode=-1 (same as imode=0, but starting from existing grids)
    ! imode=1  (computing bounding envelope)
    ! imode=2  (event generation)
    argc = COMMAND_ARGUMENT_COUNT()
    
    do i = 1, argc
       CALL GET_COMMAND_ARGUMENT(i, argv)
       if (i.eq.1) read(argv,*) nfiles
       if (.not.allocated(files)) allocate(files(nfiles))
       do j=1,nfiles
          if (i-1.eq.j) read(argv,*) files(j)
       enddo
    enddo

    do j=1,nfiles
       open(unit=10+j,file=files(j),status='old')
    enddo

    open(unit=100,file=trim(adjustl(files(1)))//'.comb',status='unknown')

  end subroutine get_run_arguments

  subroutine read_event(iunit,done)
    use comb_events
    implicit none
    integer :: i,iunit
    logical :: done
    character :: dummy
    real(kind=8) :: dum
    done=.false.
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) next,evt_full !,wgt,amp2,weight
    if (.not.allocated(helicity)) allocate(helicity(next))
    if (.not.allocated(col_order)) allocate(col_order(next))
    if (.not.allocated(IPDG)) allocate(IPDG(next))
    if (.not.allocated(momenta)) allocate(momenta(1:3,next))
    read (iunit,*,err=99,end=99) helicity(1:next)
    read (iunit,*,err=99,end=99) col_order(1:next)
    read (iunit,*,err=99,end=99) rwgt_full,rwgt_NLC!,matrix2(1),matrix2(2),matrix2(3)
    read (iunit,*,err=99,end=99) evt_wgt,evt_NLC,evt_full
    do i=1,next
       read (iunit,*,err=99,end=99) iPDG(i),momenta(1:3,i),momenta(0,i)
    enddo
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event

  subroutine read_event_init(iunit,done)
    use comb_events
    implicit none
    integer :: i,iunit
    logical :: done
    character :: dummy
    real(kind=8) :: dum
    done=.false.
    !call read_unique_in_file(iunit)
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) xsec, xerr, xsecabs, integ
    read (iunit,*,err=99,end=99) dummy
    read (iunit,*,err=99,end=99) dummy
    return
99  done=.true.
  end subroutine read_event_init

  subroutine read_unique_in_file(iunit)
    implicit none
    integer :: iproc,iunit
    read(iunit,*) next,unique_nproc
    if (.not.allocated(unique_map)) then
    allocate(unique_map(unique_nproc))
    allocate(unique_map_value(unique_nproc))
    allocate(unique_processes(next,unique_nproc))
    endif
    do iproc=1,unique_nproc
       read(iunit,*) unique_map(iproc),unique_map_value(iproc),unique_processes(1:next,iproc)
    enddo
  end subroutine read_unique_in_file

  subroutine write_event(iunit)
    use comb_events
    implicit none
    integer :: i,iunit
    write (iunit,*) '<event>'
    write (iunit,*) next,evt_full !,wgt,matrix2,weight
    write (iunit,'(100i3)') helicity(1:next)
    write (iunit,'(100i3)') col_order(1:next)
    write (iunit,*) rwgt_full,rwgt_NLC !,matrix2(1),matrix2(2),matrix2(3)
    write (iunit,*) evt_wgt,evt_NLC,evt_full
    do i=1,next
       write (iunit,*) iPDG(i),momenta(1:3,i),momenta(0,i)
    enddo
    write (iunit,*) '</event>'
  end subroutine write_event
  subroutine write_inits(iunit)
     use comb_events
     implicit none
     integer :: iunit
     write(iunit,'(a)') '<LesHouchesEvents version="3.0">'
     write(iunit,'(a)') '<init>'
     write(iunit,'(a)') '2212 2212 6.500000e+03 6.500000e+03 0 0 315200 315200 -4 1'
     write(iunit,'(es13.7,es14.7,es14.7,i2)') xsec, xerr, xsecabs, integ
     write(iunit,'(a)') "<generator name='BG_colour' version='1.0'>please cite 2409.12128 </generator>"
     write(iunit,'(a)') '</init>'
  end subroutine write_inits
  subroutine write_end(iunit)
    implicit none
    integer :: iunit
    write(iunit,*) '</LesHouchesEvents>'
  end subroutine write_end





end program matrix_combine
